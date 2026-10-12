"""
Mapa de favorabilidade para fosfato com Maxent (ee.Classifier.amnhMaxent) sobre
um mosaico final ASTER — presença + fundo, saída contínua 0–1 (degradê do pixel
mais parecido com as ocorrências ao menos parecido).

Entrada
  ASSET_MOSAICO  mosaico final (ee.Image) — PREENCHER antes de rodar. Usa as
                 bandas de BANDAS que existirem nele; se tiver a banda
                 'saturado', os pixels saturados saem (nas cenas zero-cloud eram
                 nuvem).
  ASSET_FOSFATO  ocorrências/depósitos de fosfato (pontos SGB, erro de posição
                 250–1.000 m).

Features (FEATURES_MODO)
  'bandas'          só as bandas B01–B14 do mosaico
  'indices'         só os índices de features_process/calculo_spectrals_index.js
  'bandas+indices'  os dois (padrão)
  Índices (sem duplicatas — os do JS que são a mesma razão entram uma vez):
    vegetação   NDVI, SAVI
    ferro       ferric_iron_2_1 (B02/B01), ferric_ox_4_3 (B04/B3N),
                gossan_4_2 (B04/B02), ferrous_sil_5_4 (B05/B04)
    argilas     sericite_illite_smectite ((B05+B07)/B06),
                alunite_kaolinite_pyrophyllite ((B04+B06)/B05),
                muscovite_7_6, kaolinite_7_5
    carbonato   dolomite ((B06+B08)/B07), carbonate_8_6 (B08/B06),
                carbonate_13_14 (B13/B14, TIR)
    sílica      quartz_14_12 (B14/B12, TIR)
    composições sultan_r (B04/B07), sultan_g (B04/B01), sultan_b (B02/B3N·B04/B3N)
    textura     GLCM da B3N (janela 3×3): contrast, var, ent, hom (= idm do GEE)
    fosforito   phosphorite_rbd = (ε11 + ε13) / (2·ε12) — profundidade relativa da
                feição de emissividade do apatita em ~9 µm (B12; ombros B11 e B13),
                Guha et al. 2019 (Remote Sensing 11:1003). Emissividade pelo método
                de normalização (ε_max = 0,96) sobre a radiância TOA do TIR (DN → L
                pelos UCC do ASTER L1T) — sem correção atmosférica.
  Um índice só entra se o mosaico tiver as bandas dele (sem SWIR → sem os de
  SWIR). Razões não dependem da escala × 1000; o SAVI usa reflectância (/1000).
  Os de TIR são razões de DN, não de emissividade — valem como contraste
  relativo.

Amostras
  presença (1)   até AMOSTRAS_POR_PONTO pixels sorteados num buffer de BUFFER_M
                 em volta de cada ponto (o erro de posição é maior que o pixel).
  fundo (0)      N_FUNDO pixels aleatórios da grid, distribuídos pelas células.
                 Não são "ausências" — Maxent contrasta presença × paisagem.

Etapas (menu)
  1) amostras → ASSET_AMOSTRAS (task toAsset). A amostragem é pesada demais
     para rodar interativa ("Too many concurrent aggregations"); validação e
     mapa leem a tabela salva.
  2) validação e 3) mapa, abaixo — só depois que a task 1 terminar.

Validação (opção 2): bloco espacial, deixa UMA CÉLULA DA GRID de fora por vez.
  Treina sem as presenças e sem o fundo da célula; prediz na célula deixada de
  fora e calcula a AUC presença × fundo dessa célula e a probabilidade mediana
  nos pontos deixados de fora. Pontos da mesma célula são vizinhos (vários
  depósitos agrupados) — deixá-los juntos de fora evita AUC inflada.
  Resultado → models/maxent_fosfato_validacao_<mosaico>.csv

Mapa (opção 3): treina com todas as amostras e exporta a banda 'probability'
  (float 0–1) para ASSET_SAIDA, na região da grid, a ESCALA m.

Conta: a ativa (precisa escrever em projects/mapbiomas-arida).

Uso: python maxent_fosfato.py
"""

import logging
import sys
from pathlib import Path

import ee
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from configure_account_projects_ee import get_current_account   # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
ASSET_MOSAICO = ''      # ← mosaico final (ee.Image). Ex.: 'projects/mapbiomas-arida/...'
ASSET_FOSFATO = 'projects/mapbiomas-arida/mine/pontos_fosfato'
GRID_ASSET    = 'projects/mapbiomas-arida/mine/grid_irece_5x5_solape100'
PASTA_SAIDA   = 'projects/mapbiomas-arida/mine'
MODELS_DIR    = Path(__file__).resolve().parent / 'models'

BANDAS = ['B01', 'B02', 'B3N',
          'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
          'B10', 'B11', 'B12', 'B13', 'B14']
BUFFER_M           = 500     # metade do erro máximo de posição dos pontos
AMOSTRAS_POR_PONTO = 20
N_FUNDO            = 10000
ESCALA             = 30      # m — amostragem e export
SEED               = 42
FEATURES_MODO      = 'bandas+indices'   # 'bandas', 'indices' ou 'bandas+indices'
ESCALA_TOA         = 1000               # valor no mosaico = reflectância × 1000

# Maxent: poucas ocorrências → regularização maior e só funções simples
MAXENT = dict(
    autoFeature=False,
    linear=True, quadratic=True,
    product=False, threshold=False, hinge=False,
    betaMultiplier=2.0,
    addSamplesToBackground=True,
    seed=SEED,
)


# ═══════════════════════════════════════════════════════════════════════════
# DADOS
# ═══════════════════════════════════════════════════════════════════════════
def ler_mosaico():
    if not ASSET_MOSAICO:
        raise SystemExit('Preencha ASSET_MOSAICO no início do script com o mosaico final.')
    img = ee.Image(ASSET_MOSAICO)
    presentes = img.bandNames().getInfo()
    bandas = [b for b in BANDAS if b in presentes]
    faltando = [b for b in BANDAS if b not in presentes]
    if not bandas:
        raise SystemExit(f'Nenhuma banda de {BANDAS} em {ASSET_MOSAICO} (tem: {presentes}).')
    if faltando:
        log.warning(f'Bandas ausentes no mosaico (ficam de fora): {faltando}')
    if 'saturado' in presentes:
        img = img.updateMask(img.select('saturado').eq(0))
    feats, nomes = calcular_features(img.select(bandas).toFloat(), bandas)
    log.info(f'Mosaico {ASSET_MOSAICO}: {len(bandas)} bandas; {len(nomes)} features ({FEATURES_MODO}): {nomes}')
    return feats, nomes


# nome → (bandas necessárias, função img → ee.Image)
def _razao(a, b):
    return lambda i: i.select(a).divide(i.select(b))

INDICES = {
    'NDVI': (['B3N', 'B02'], lambda i: i.normalizedDifference(['B3N', 'B02'])),
    'SAVI': (['B3N', 'B02'], lambda i: i.expression(
        '1.5 * (n - r) / (n + r + 0.5)',
        {'n': i.select('B3N').divide(ESCALA_TOA), 'r': i.select('B02').divide(ESCALA_TOA)})),
    'ferric_iron_2_1': (['B02', 'B01'], _razao('B02', 'B01')),
    'ferric_ox_4_3':   (['B04', 'B3N'], _razao('B04', 'B3N')),
    'gossan_4_2':      (['B04', 'B02'], _razao('B04', 'B02')),
    'ferrous_sil_5_4': (['B05', 'B04'], _razao('B05', 'B04')),
    'sericite_illite_smectite': (['B05', 'B06', 'B07'],
                                 lambda i: i.select('B05').add(i.select('B07')).divide(i.select('B06'))),
    'alunite_kaolinite_pyrophyllite': (['B04', 'B05', 'B06'],
                                       lambda i: i.select('B04').add(i.select('B06')).divide(i.select('B05'))),
    'muscovite_7_6':   (['B07', 'B06'], _razao('B07', 'B06')),
    'kaolinite_7_5':   (['B07', 'B05'], _razao('B07', 'B05')),
    'dolomite':        (['B06', 'B07', 'B08'],
                        lambda i: i.select('B06').add(i.select('B08')).divide(i.select('B07'))),
    'carbonate_8_6':   (['B08', 'B06'], _razao('B08', 'B06')),
    'carbonate_13_14': (['B13', 'B14'], _razao('B13', 'B14')),
    'quartz_14_12':    (['B14', 'B12'], _razao('B14', 'B12')),
    'sultan_r':        (['B04', 'B07'], _razao('B04', 'B07')),
    'sultan_g':        (['B04', 'B01'], _razao('B04', 'B01')),
    'sultan_b':        (['B02', 'B3N', 'B04'],
                        lambda i: i.select('B02').divide(i.select('B3N'))
                                   .multiply(i.select('B04').divide(i.select('B3N')))),
}
# GLCM da B3N: nome no GEE → nome da feature ('idm' é a homogeneidade; o GEE não tem 'hom')
TEXTURAS = {'contrast': 'contrast', 'var': 'var', 'ent': 'ent', 'idm': 'hom'}

# TIR: DN → radiância (W m⁻² sr⁻¹ µm⁻¹) = (DN − 1) × UCC; λ central (µm)
TIR_UCC   = {'B10': 6.822e-3, 'B11': 6.780e-3, 'B12': 6.590e-3, 'B13': 5.693e-3, 'B14': 5.225e-3}
TIR_LAMBDA = {'B10': 8.291, 'B11': 8.634, 'B12': 9.075, 'B13': 10.657, 'B14': 11.318}
EPS_MAX   = 0.96
C1, C2    = 1.191042e8, 1.4387752e4   # Planck: W µm⁴ m⁻² sr⁻¹ ; µm K


def emissividade_tir(img):
    """Método de normalização de emissividade: T de cada banda com ε = EPS_MAX,
    T_max por pixel, ε_i = L_i / B(λ_i, T_max)."""
    L = {b: img.select(b).subtract(1).max(1).multiply(u) for b, u in TIR_UCC.items()}
    T = {b: L[b].expression('c2 / (lam * log(eps * c1 / (pow(lam, 5) * L) + 1))',
                            {'c1': C1, 'c2': C2, 'lam': TIR_LAMBDA[b], 'eps': EPS_MAX, 'L': L[b]})
         for b in TIR_UCC}
    tmax = ee.Image.cat(list(T.values())).reduce(ee.Reducer.max())
    return {b: L[b].divide(tmax.expression('c1 / (pow(lam, 5) * (exp(c2 / (lam * t)) - 1))',
                                           {'c1': C1, 'c2': C2, 'lam': TIR_LAMBDA[b], 't': tmax}))
            for b in TIR_UCC}


def calcular_features(img, bandas: list[str]):
    """Bandas e/ou índices conforme FEATURES_MODO. Devolve (imagem float, nomes)."""
    partes, nomes = [], []
    if 'bandas' in FEATURES_MODO:
        partes.append(img.select(bandas))
        nomes += bandas
    if 'indices' in FEATURES_MODO:
        for nome, (precisa, f) in INDICES.items():
            if all(b in bandas for b in precisa):
                partes.append(f(img).toFloat().rename(nome))
                nomes.append(nome)
        if all(b in bandas for b in TIR_UCC):
            eps = emissividade_tir(img)
            partes.append(eps['B11'].add(eps['B13']).divide(eps['B12'].multiply(2))
                          .toFloat().rename('phosphorite_rbd'))
            nomes.append('phosphorite_rbd')
        if 'B3N' in bandas:
            # B3N (refl. × 1000) 0–350 → 0–255 em byte para a GLCM
            b3n = img.select('B3N').unitScale(0, 0.35 * ESCALA_TOA).multiply(255).clamp(0, 255).toByte()
            glcm = b3n.glcmTexture(size=1)
            partes.append(glcm.select([f'B3N_{t}' for t in TEXTURAS],
                                      [f'glcm_{n}' for n in TEXTURAS.values()]).toFloat())
            nomes += [f'glcm_{n}' for n in TEXTURAS.values()]
    feats = ee.Image.cat(partes).toFloat()
    # pixel só entra com todas as features válidas (razão com divisor 0 → mascarada)
    return feats.updateMask(feats.mask().reduce(ee.Reducer.min())), nomes


def ler_grid() -> dict[int, dict]:
    info = ee.FeatureCollection(GRID_ASSET).getInfo()
    return {int(f['properties']['id']): {'nome': f['properties']['nome'],
                                         'geom': ee.Geometry(f['geometry'], None, False)}
            for f in info['features']}


def celula_dos_pontos(pontos, grid) -> dict[int, int]:
    """CODIGO_OBJ → id da célula (a primeira que contém o ponto; a grid tem solape)."""
    d = ee.Dictionary({str(cid): pontos.filterBounds(c['geom']).aggregate_array('CODIGO_OBJ')
                       for cid, c in grid.items()}).getInfo()
    cel = {}
    for cid in sorted(d, key=int):
        for cod in d[cid]:
            cel.setdefault(int(cod), int(cid))
    return cel


def amostras_presenca(feats, pontos, cel_ponto):
    cel = ee.Dictionary({str(k): v for k, v in cel_ponto.items()})

    def por_ponto(f):
        cod = f.get('CODIGO_OBJ')
        props = {'classe': 1, 'CODIGO_OBJ': cod, 'celula': cel.get(ee.Number(cod).format('%d'))}
        return (feats.sample(region=f.geometry().buffer(BUFFER_M), scale=ESCALA,
                             numPixels=AMOSTRAS_POR_PONTO * 5, seed=SEED, geometries=True)
                     .limit(AMOSTRAS_POR_PONTO)
                     .map(lambda s: s.set(props)))

    pontos = pontos.filter(ee.Filter.inList('CODIGO_OBJ', list(cel_ponto)))
    return ee.FeatureCollection(pontos.map(por_ponto)).flatten()


def amostras_fundo(feats, grid):
    n = N_FUNDO // len(grid)
    partes = [feats.sample(region=c['geom'], scale=ESCALA, numPixels=n * 3,
                           seed=SEED + cid, geometries=True)
                   .limit(n)
                   .map(lambda s, cid=cid: s.set({'classe': 0, 'celula': cid}))
              for cid, c in grid.items()]
    return ee.FeatureCollection(partes).flatten()


def treinar(amostras, bandas):
    return (ee.Classifier.amnhMaxent(**MAXENT)
              .setOutputMode('PROBABILITY')
              .train(features=amostras, classProperty='classe', inputProperties=bandas))


# ═══════════════════════════════════════════════════════════════════════════
# VALIDAÇÃO — UMA CÉLULA DE FORA POR VEZ
# ═══════════════════════════════════════════════════════════════════════════
def auc(pres: np.ndarray, fundo: np.ndarray) -> float:
    """P(prob. presença > prob. fundo), empates contam meio (Mann-Whitney)."""
    if not len(pres) or not len(fundo):
        return float('nan')
    maior = (pres[:, None] > fundo[None, :]).mean()
    igual = (pres[:, None] == fundo[None, :]).mean()
    return float(maior + 0.5 * igual)


def validar(pres, fundo, bandas, grid, cel_ponto):
    linhas = []
    celulas = sorted(set(cel_ponto.values()))
    log.info(f'Validação: {len(celulas)} células com pontos, uma de fora por vez')
    todas_p, todas_f = [], []
    for cid in celulas:
        fora = ee.Filter.eq('celula', cid)
        clf = treinar(pres.filter(fora.Not()).merge(fundo.filter(fora.Not())), bandas)
        r = ee.Dictionary({
            'p': pres.filter(fora).classify(clf).aggregate_array('probability'),
            'f': fundo.filter(fora).classify(clf).aggregate_array('probability'),
        }).getInfo()
        p, f = np.array(r['p']), np.array(r['f'])
        todas_p.append(p)
        todas_f.append(f)
        n_pts = sum(1 for c in cel_ponto.values() if c == cid)
        linhas.append({'celula': grid[cid]['nome'], 'pontos': n_pts, 'amostras_presenca': len(p),
                       'amostras_fundo': len(f), 'auc': auc(p, f),
                       'prob_mediana_presenca': float(np.median(p)) if len(p) else np.nan,
                       'prob_mediana_fundo': float(np.median(f)) if len(f) else np.nan})
        log.info(f'  {grid[cid]["nome"]}: {n_pts} pontos  AUC={linhas[-1]["auc"]:.3f}  '
                 f'mediana presença={linhas[-1]["prob_mediana_presenca"]:.3f}  '
                 f'fundo={linhas[-1]["prob_mediana_fundo"]:.3f}')
    df = pd.DataFrame(linhas)
    geral = auc(np.concatenate(todas_p), np.concatenate(todas_f))
    log.info(f'AUC geral (todas as células deixadas de fora): {geral:.3f}')
    MODELS_DIR.mkdir(parents=True, exist_ok=True)
    arq = MODELS_DIR / f'maxent_fosfato_validacao_{nome_mosaico()}.csv'
    df.to_csv(arq, index=False)
    log.info(f'Validação → {arq}')
    print(df.round(3).to_string(index=False))


# ═══════════════════════════════════════════════════════════════════════════
# MAPA
# ═══════════════════════════════════════════════════════════════════════════
def nome_mosaico() -> str:
    return ASSET_MOSAICO.rstrip('/').rsplit('/', 1)[-1]


def asset_amostras() -> str:
    modo = FEATURES_MODO.replace('+', '_')
    return f'{PASTA_SAIDA}/maxent_fosfato_amostras_{nome_mosaico()}_{modo}'


def exportar_amostras(feats, pontos, grid, cel_ponto):
    asset = asset_amostras()
    amostras = amostras_presenca(feats, pontos, cel_ponto).merge(amostras_fundo(feats, grid))
    nome = asset.rsplit('/', 1)[-1]
    amostras = amostras.map(lambda f: f.set('features_modo', FEATURES_MODO))
    task = ee.batch.Export.table.toAsset(collection=amostras, description=nome, assetId=asset)
    task.start()
    log.info(f'Amostras: task {task.id} → {asset}  (rode a validação/mapa quando terminar)')


def ler_amostras():
    asset = asset_amostras()
    try:
        ee.data.getAsset(asset)
    except ee.EEException:
        raise SystemExit(f'{asset} não existe — rode a opção 1 e espere a task terminar.')
    fc = ee.FeatureCollection(asset)
    n = fc.aggregate_histogram('classe').getInfo()
    log.info(f'Amostras {asset}: presença {n.get("1", 0)}, fundo {n.get("0", 0)}')
    return fc.filter(ee.Filter.eq('classe', 1)), fc.filter(ee.Filter.eq('classe', 0))


def exportar(feats, pres, fundo, bandas, grid):
    clf = treinar(pres.merge(fundo), bandas)
    regiao = ee.Geometry.MultiPolygon([c['geom'].coordinates() for c in grid.values()]).dissolve()
    mapa = (feats.classify(clf).select('probability').toFloat()
                 .set({'mosaico': ASSET_MOSAICO, 'pontos': ASSET_FOSFATO, 'bandas': ','.join(bandas),
                       'buffer_m': BUFFER_M, 'amostras_por_ponto': AMOSTRAS_POR_PONTO,
                       'n_fundo': N_FUNDO, 'maxent': str(MAXENT)}))
    nome = f'maxent_fosfato_{nome_mosaico()}'
    asset = f'{PASTA_SAIDA}/{nome}'
    task = ee.batch.Export.image.toAsset(image=mapa, description=nome, assetId=asset, region=regiao,
                                         scale=ESCALA, crs='EPSG:4326', maxPixels=1e13,
                                         pyramidingPolicy={'.default': 'mean'})
    task.start()
    log.info(f'Mapa: task {task.id} → {asset}')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ee.Initialize(project=get_current_account())
    feats, bandas = ler_mosaico()
    grid = ler_grid()
    pontos = ee.FeatureCollection(ASSET_FOSFATO)
    cel_ponto = celula_dos_pontos(pontos, grid)
    log.info(f'{len(cel_ponto)} pontos de fosfato dentro da grid, em {len(set(cel_ponto.values()))} células')

    print(f'\nO que fazer?\n  1) amostras (presença + fundo) → {asset_amostras()}'
          '\n  2) validação (uma célula de fora por vez)\n  3) exportar o mapa\n  4) validação + mapa')
    acao = input('Opção [1]: ').strip() or '1'
    if acao == '1':
        return exportar_amostras(feats, pontos, grid, cel_ponto)
    pres, fundo = ler_amostras()
    if acao in ('2', '4'):
        validar(pres, fundo, bandas, grid, cel_ponto)
    if acao in ('3', '4'):
        exportar(feats, pres, fundo, bandas, grid)


if __name__ == '__main__':
    main()
