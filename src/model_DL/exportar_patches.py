"""
Exporta patches 64 × 64 px (30 m) do mosaico ASTER v1 para treinar a CNN de
favorabilidade de fosfato (plano em plano_cnn.md). Sem tasks: cada patch vem
direto do GEE por ee.data.computePixels (array NumPy) e é gravado em .npz.

Amostras (centro do patch)
  ocorrencia_treino  ocorrências de fosfato com papel = treino (bacia):
                     o ponto + K_JITTER centros num raio de RAIO_JITTER_M,
                     mantidos só onde RF e GTB têm p ≥ P_POS (área magenta).
                     y = 1, peso = 1.
  ocorrencia_teste   ocorrências com papel = validacao (leste/norte): o ponto +
                     K_JITTER centros, SEM filtro de p (avaliação não pode
                     depender do que RF/GTB já acertam). y = 1, peso = 1.
  ambos              pseudo-positivo: RF e GTB com p ≥ P_POS. y = média(p_RF,
                     p_GTB), peso = 1 (rótulo suave: a referência não é 100%).
  negativo           RF e GTB com p ≤ P_NEG, a mais de DIST_NEG_OCORR_M das
                     ocorrências. y = média(p_RF, p_GTB) ≈ 0, peso = 1.
  Só RF / só GTB e a faixa P_NEG < p < P_POS ficam de fora.

Divisão (por célula da grid 5×5 sem solape, divisao.json)
  teste    células com ocorrências de validação — usadas uma vez, no fim
  fold k   demais células em N_FOLDS grupos (validação cruzada por grupo); as
           células com ocorrências de treino vão para folds diferentes
  Centros a menos de MARGEM_BORDA_M de uma borda com célula de outro grupo
  são descartados (o patch de 64 px não pode ter pixel nos dois grupos).

Canais (baixados uma vez, gravados num arquivo por modelo)
  modelo A  15 índices espectrais (mesmas fórmulas de modeling/maxent_fosfato.py)
            + swir_estimado → 16 canais
  modelo B1 log da reflectância B01–B09 (reflectância = valor / 1000) e
            emissividade B10–B14 (normalização, ε_max 0,96) + swir_estimado
            → 15 canais
  origem    não entra na rede: vai em 'origem_patch' (64 × 64) nos dois
            arquivos, para conferir bordas de cena
  VNIR (15 m) agregado pela média para 30 m; grade alinhada ao mosaico v1.
  Pixel sem dado = NaN.

Saídas (PASTA_SAIDA = /srv/almacen/db_images/patches_cnn ou $PATCHES_DIR)
  modelo_A/patches_<LxCy>.npz   X (N, 16, 64, 64) float16
  modelo_B1/patches_<LxCy>.npz  X (N, 15, 64, 64) float16
     os dois com: canais, origem_patch (N, 64, 64), y, peso, origem_amostra,
     id, conjunto, fold, celula, lon, lat, p_rf, p_gtb, papel, codigo_obj, jitter
     (mesma ordem de patches nos dois arquivos, mesmo id)
  indice_patches.csv   todas as amostras e metadados (sem X)
  divisao.json         células de teste, folds e parâmetros

Atenção — vazamento: os mapas de favorabilidade_ml usados aqui foram
treinados com todos os pontos (propriedade papeis = 'None'), inclusive os de
validação que vão para o teste. Para um teste independente, refazer os mapas
com PAPEIS = ['treino'] numa coleção separada e apontar ASSET_FAV para ela.

Uso: python exportar_patches.py   (1 = só sortear e gravar índice/divisão;
                                   2 = sortear e baixar os patches)
"""

import json
import logging
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import ee
import numpy as np
import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
sys.path.append(str(RAIZ))
sys.path.append(str(RAIZ / 'modeling'))
import maxent_fosfato as mx   # noqa: E402  (fórmulas dos índices e emissividade)

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
PROJETO_EE    = 'mapbiomas-brazil'
MINE          = 'projects/mapbiomas-arida/mine/'
ASSET_MOSAICO = MINE + 'aster-mosaico-final'          # mosaico v1
ASSET_FAV     = MINE + 'favorabilidade_ml'            # RF/GTB (probability × 10000)
MOSAICO_FAV   = 'aster-mosaico-final'
ALVO          = 'fosfato'
ASSET_PONTOS  = MINE + 'pontos_referencia'
ASSET_GRID    = MINE + 'grid_irece_5x5'               # sem solape
# servidor: /srv/almacen/db_images/patches_cnn; outra pasta com a variável de
# ambiente PATCHES_DIR (ex.: PATCHES_DIR=~/dados/patches python exportar_patches.py)
PASTA_SAIDA   = Path(os.environ.get('PATCHES_DIR', '/srv/almacen/db_images/patches_cnn')).expanduser()

PATCH         = 64          # px
RES_15M       = 0.00013474729261792824                # grau, pixel do mosaico v1
RES           = 2 * RES_15M                           # 30 m, alinhado aos blocos 2×2
X0, Y0        = -44.4606182671107, -9.416390538076989 # canto NW da grade do mosaico

P_POS            = 0.5      # positivo: RF e GTB ≥ P_POS
P_NEG            = 0.05     # negativo: RF e GTB ≤ P_NEG
DIST_NEG_OCORR_M = 2000     # negativos longe das ocorrências
N_POS_CELULA     = 120      # pseudo-positivos sorteados por célula (antes da rarefação)
N_NEG_CELULA     = 150      # negativos sorteados por célula
DIST_MIN_POS_M   = 1000     # distância mínima entre pseudo-positivos
DIST_MIN_NEG_M   = 1500     # distância mínima entre negativos
K_JITTER         = 8        # centros extras por ocorrência
RAIO_JITTER_M    = 500      # erro de posição dos pontos SGB
MARGEM_BORDA_M   = 1000     # meio patch (64 px × 30 m / 2)
N_FOLDS          = 5
SEED             = 42
N_THREADS        = 8
TENTATIVAS       = 4
NODATA           = -9999.0

CANAIS_A = ['ferric_iron_2_1', 'ferric_ox_4_3', 'gossan_4_2', 'ferrous_sil_5_4', 'SAVI',
            'sericite_illite_smectite', 'kaolinite_7_5', 'muscovite_7_6',
            'carbonate_8_6', 'dolomite', 'carbonate_13_14',
            'quartz_14_12', 'phosphorite_rbd', 'sultan_r', 'sultan_g']
REFLETIVAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09']
TERMAIS    = ['B10', 'B11', 'B12', 'B13', 'B14']
CANAIS_B1  = [f'log_{b}' for b in REFLETIVAS] + [f'emis_{b}' for b in TERMAIS]
EXTRAS     = ['swir_estimado', 'origem']
CANAIS     = CANAIS_A + CANAIS_B1 + EXTRAS          # o que é baixado do GEE
MODELOS    = {'A':  CANAIS_A + ['swir_estimado'],    # o que vai em cada arquivo
              'B1': CANAIS_B1 + ['swir_estimado']}


# ═══════════════════════════════════════════════════════════════════════════
# IMAGENS NO GEE
# ═══════════════════════════════════════════════════════════════════════════
def imagem_canais():
    """Imagem com todos os CANAIS, na projeção de 15 m do mosaico; a média 2×2
    para 30 m sai do reduceResolution quando o computePixels pede a grade de 30 m."""
    col = ee.ImageCollection(ASSET_MOSAICO)
    proj15 = col.first().select('B01').projection()
    v1 = col.mosaic().setDefaultProjection(proj15)
    bandas = v1.select(REFLETIVAS + TERMAIS).toFloat()

    indices = [mx.INDICES[n][1](bandas).toFloat().rename(n) for n in CANAIS_A if n in mx.INDICES]
    eps = mx.emissividade_tir(bandas)
    fosforito = eps['B11'].add(eps['B13']).divide(eps['B12'].multiply(2)).toFloat().rename('phosphorite_rbd')
    ind = ee.Image.cat(indices + [fosforito]).select(CANAIS_A)

    logs = ee.Image.cat([bandas.select(b).divide(mx.ESCALA_TOA).max(1e-4).log().rename(f'log_{b}')
                         for b in REFLETIVAS])
    emis = ee.Image.cat([eps[b].toFloat().rename(f'emis_{b}') for b in TERMAIS])

    continuas = ind.addBands(logs).addBands(emis).reduceResolution(ee.Reducer.mean(), maxPixels=16)
    categ = v1.select(EXTRAS).toFloat().reduceResolution(ee.Reducer.mode(), maxPixels=16).rename(EXTRAS)
    return continuas.addBands(categ).select(CANAIS).unmask(NODATA)


def imagem_p():
    """p_rf e p_gtb (0–1) do alvo no mosaico escolhido."""
    def p(modelo):
        return (ee.ImageCollection(ASSET_FAV)
                  .filter(ee.Filter.eq('modelo', modelo))
                  .filter(ee.Filter.eq('alvo', ALVO))
                  .filter(ee.Filter.eq('mosaico', MOSAICO_FAV))
                  .mosaic().divide(10000))
    return ee.Image.cat([p('RF').rename('p_rf'), p('GTB').rename('p_gtb')])


# ═══════════════════════════════════════════════════════════════════════════
# GEOMETRIA LOCAL
# ═══════════════════════════════════════════════════════════════════════════
def metros(lon1, lat1, lon2, lat2):
    dx = (lon2 - lon1) * 111320 * math.cos(math.radians((lat1 + lat2) / 2))
    dy = (lat2 - lat1) * 110540
    return math.hypot(dx, dy)


def rarefazer(df, dist_m):
    """Mantém o primeiro de cada grupo a menos de dist_m (ordem embaralhada)."""
    df = df.sample(frac=1, random_state=SEED)
    ficam = []
    for r in df.itertuples():
        if all(metros(r.lon, r.lat, o.lon, o.lat) >= dist_m for o in ficam):
            ficam.append(r)
    return pd.DataFrame(ficam).drop(columns='Index', errors='ignore')


def celula_de(lon, lat, grid):
    for nome, c in grid.items():
        if c['oeste'] <= lon < c['leste'] and c['sul'] < lat <= c['norte']:
            return nome
    return None


def perto_de_outro_grupo(lon, lat, nome, grid, grupo):
    """True se o centro está a menos de MARGEM_BORDA_M de uma célula vizinha de outro grupo."""
    c = grid[nome]
    vizinhos = {(-1, 0): c['norte'], (1, 0): c['sul'], (0, -1): c['oeste'], (0, 1): c['leste']}
    for (dl, dc), borda in vizinhos.items():
        viz = next((n for n, v in grid.items()
                    if v['linha'] == c['linha'] + dl and v['coluna'] == c['coluna'] + dc), None)
        if viz is None or grupo[viz] == grupo[nome]:
            continue
        d = metros(lon, lat, lon, borda) if dl else metros(lon, lat, borda, lat)
        if d < MARGEM_BORDA_M:
            return True
    return False


# ═══════════════════════════════════════════════════════════════════════════
# SORTEIO
# ═══════════════════════════════════════════════════════════════════════════
def ler_grid():
    info = ee.FeatureCollection(ASSET_GRID).getInfo()['features']
    return {f['properties']['nome']: {**f['properties'], 'geom': ee.Geometry(f['geometry'])} for f in info}


def ler_ocorrencias():
    fc = ee.FeatureCollection(ASSET_PONTOS).filter(ee.Filter.eq('alvo', ALVO))
    return pd.DataFrame([{'lon': f['geometry']['coordinates'][0], 'lat': f['geometry']['coordinates'][1],
                          'papel': f['properties'].get('papel', ''),
                          'codigo_obj': f['properties'].get('CODIGO_OBJ', -1)}
                         for f in fc.getInfo()['features']])


def amostrar_p(df, pimg):
    """p_rf e p_gtb nos centros (getInfo em lotes)."""
    saidas = []
    for i in range(0, len(df), 2000):
        lote = df.iloc[i:i + 2000]
        fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon, r.lat]), {'k': int(k)})
                                   for k, r in lote.iterrows()])
        res = pimg.reduceRegions(collection=fc, reducer=ee.Reducer.first(), scale=30).getInfo()['features']
        saidas += [{'k': f['properties']['k'], 'p_rf': f['properties'].get('p_rf'),
                    'p_gtb': f['properties'].get('p_gtb')} for f in res]
    p = pd.DataFrame(saidas).set_index('k')
    return df.join(p)


def centros_ocorrencias(oc, rng):
    linhas = []
    for r in oc.itertuples():
        linhas.append({'lon': r.lon, 'lat': r.lat, 'papel': r.papel, 'codigo_obj': r.codigo_obj, 'jitter': 0})
        for j in range(1, K_JITTER + 1):
            ang, rr = rng.uniform(0, 2 * math.pi), RAIO_JITTER_M * math.sqrt(rng.uniform())
            linhas.append({'lon': r.lon + rr * math.cos(ang) / (111320 * math.cos(math.radians(r.lat))),
                           'lat': r.lat + rr * math.sin(ang) / 110540,
                           'papel': r.papel, 'codigo_obj': r.codigo_obj, 'jitter': j})
    return pd.DataFrame(linhas)


def sortear_pseudo(grid, pimg, oc):
    """Positivos (ambos ≥ P_POS) e negativos (ambos ≤ P_NEG, longe das ocorrências) por célula."""
    pts = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r.lon, r.lat])) for r in oc.itertuples()])
    longe = pts.distance(searchRadius=DIST_NEG_OCORR_M + 100).gt(DIST_NEG_OCORR_M).unmask(1)
    ambos_pos = pimg.select('p_rf').gte(P_POS).And(pimg.select('p_gtb').gte(P_POS))
    ambos_neg = pimg.select('p_rf').lte(P_NEG).And(pimg.select('p_gtb').lte(P_NEG)).And(longe)
    classe = (ee.Image(0).where(ambos_pos, 3).where(ambos_neg, 1)
                .updateMask(ambos_pos.Or(ambos_neg)).rename('classe').toInt())
    img = classe.addBands(pimg)
    linhas = []
    for nome, c in sorted(grid.items()):
        res = img.stratifiedSample(numPoints=0, classBand='classe', region=c['geom'], scale=30,
                                   classValues=[1, 3], classPoints=[N_NEG_CELULA, N_POS_CELULA],
                                   seed=SEED, geometries=True, tileScale=4).getInfo()['features']
        for f in res:
            pr = f['properties']
            linhas.append({'lon': f['geometry']['coordinates'][0], 'lat': f['geometry']['coordinates'][1],
                           'origem_amostra': 'ambos' if pr['classe'] == 3 else 'negativo',
                           'p_rf': pr.get('p_rf'), 'p_gtb': pr.get('p_gtb')})
        log.info(f'  {nome}: {sum(1 for f in res if f["properties"]["classe"] == 3)} positivos, '
                 f'{sum(1 for f in res if f["properties"]["classe"] == 1)} negativos sorteados')
    return pd.DataFrame(linhas)


def dividir(grid, oc):
    """Células de teste (ocorrências de validação) e folds das demais."""
    rng = np.random.default_rng(SEED)
    cel_teste = sorted({celula_de(r.lon, r.lat, grid) for r in oc.itertuples() if r.papel == 'validacao'} - {None})
    cel_treino_oc = sorted({celula_de(r.lon, r.lat, grid) for r in oc.itertuples() if r.papel == 'treino'}
                           - {None} - set(cel_teste))
    resto = [n for n in sorted(grid) if n not in cel_teste and n not in cel_treino_oc]
    rng.shuffle(resto)
    ordem = cel_treino_oc + resto          # células com ocorrência de treino caem em folds diferentes
    fold = {n: i % N_FOLDS for i, n in enumerate(ordem)}
    grupo = {n: 'teste' if n in cel_teste else f'fold{fold[n]}' for n in grid}
    return cel_teste, fold, grupo


def montar_amostras(grid, pimg):
    rng = np.random.default_rng(SEED)
    oc = ler_ocorrencias()
    log.info(f'{len(oc)} ocorrências de {ALVO} ({oc.papel.value_counts().to_dict()})')
    cel_teste, fold, grupo = dividir(grid, oc)
    log.info(f'Teste: {cel_teste}')

    # ocorrências: ponto + jitter; treino só em magenta, teste sem filtro
    co = amostrar_p(centros_ocorrencias(oc, rng), pimg)
    magenta = (co.p_rf >= P_POS) & (co.p_gtb >= P_POS)
    tr = co[(co.papel == 'treino') & magenta].assign(origem_amostra='ocorrencia_treino')
    te = co[co.papel != 'treino'].assign(origem_amostra='ocorrencia_teste')
    sem = sorted(set(oc[oc.papel == 'treino'].codigo_obj) - set(tr.codigo_obj))
    if sem:
        log.warning(f'Ocorrências de treino sem nenhum centro em área magenta (ficam fora): {sem}')

    log.info('Sorteando pseudo-positivos e negativos por célula…')
    ps = sortear_pseudo(grid, pimg, oc)
    pos = rarefazer(ps[ps.origem_amostra == 'ambos'], DIST_MIN_POS_M)
    neg = rarefazer(ps[ps.origem_amostra == 'negativo'], DIST_MIN_NEG_M)

    df = pd.concat([tr, te, pos, neg], ignore_index=True)
    df['y'] = np.where(df.origem_amostra.str.startswith('ocorrencia'), 1.0,
                       df[['p_rf', 'p_gtb']].mean(axis=1))
    df['peso'] = 1.0
    df['celula'] = [celula_de(r.lon, r.lat, grid) for r in df.itertuples()]
    df = df[df.celula.notna()]
    borda = [perto_de_outro_grupo(r.lon, r.lat, r.celula, grid, grupo) for r in df.itertuples()]
    log.info(f'{sum(borda)} centros descartados por estar a < {MARGEM_BORDA_M} m de célula de outro grupo')
    df = df[~np.array(borda)].reset_index(drop=True)
    df['conjunto'] = np.where(df.celula.isin(cel_teste), 'teste', 'treino')
    df['fold'] = [-1 if c in cel_teste else fold[c] for c in df.celula]
    for col, padrao in [('papel', ''), ('codigo_obj', -1), ('jitter', -1)]:
        df[col] = df[col].fillna(padrao) if col in df else padrao
    df.insert(0, 'id', [f'{c}_{i:05d}' for i, c in enumerate(df.celula)])

    divisao = {'teste': cel_teste, 'folds': {n: int(f) for n, f in fold.items() if n not in cel_teste},
               'grupo': grupo, 'parametros': {k: v for k, v in globals().items()
                                               if k.isupper() and isinstance(v, (int, float, str))}}
    return df, divisao


# ═══════════════════════════════════════════════════════════════════════════
# DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════
def grade_patch(lon, lat):
    col = math.floor((lon - X0) / RES)
    lin = math.floor((Y0 - lat) / RES)
    return {'dimensions': {'width': PATCH, 'height': PATCH},
            'affineTransform': {'scaleX': RES, 'shearX': 0, 'translateX': X0 + (col - PATCH // 2) * RES,
                                'shearY': 0, 'scaleY': -RES, 'translateY': Y0 - (lin - PATCH // 2) * RES},
            'crsCode': 'EPSG:4326'}


def baixar_patch(img, lon, lat):
    for t in range(TENTATIVAS):
        try:
            arr = ee.data.computePixels({'expression': img, 'fileFormat': 'NUMPY_NDARRAY',
                                         'grid': grade_patch(lon, lat)})
            x = np.stack([arr[c] for c in CANAIS]).astype(np.float32)   # (C, H, W)
            x[x == NODATA] = np.nan
            return x
        except ee.EEException as e:
            if t == TENTATIVAS - 1:
                raise
            time.sleep(2 ** t + np.random.rand())
            log.debug(f'nova tentativa ({e})')


def baixar_celula(img, df_cel, cel):
    """Baixa os patches da célula uma vez e grava um .npz por modelo."""
    X = np.full((len(df_cel), len(CANAIS), PATCH, PATCH), np.nan, dtype=np.float16)
    ok = np.zeros(len(df_cel), bool)
    with ThreadPoolExecutor(N_THREADS) as ex:
        futuros = {ex.submit(baixar_patch, img, r.lon, r.lat): i for i, r in enumerate(df_cel.itertuples())}
        for f in as_completed(futuros):
            i = futuros[f]
            try:
                X[i] = f.result().astype(np.float16)
                ok[i] = True
            except Exception as e:   # patch perdido não para a célula
                log.warning(f'  {df_cel.iloc[i].id}: falhou ({e})')
    X, d = X[ok], df_cel[ok].reset_index(drop=True)
    # colunas de texto como str (np.load sem allow_pickle)
    meta = {c: (d[c].astype(str).to_numpy() if d[c].dtype == object else d[c].to_numpy()) for c in d.columns}
    origem = X[:, CANAIS.index('origem')]
    for modelo, canais in MODELOS.items():
        idx = [CANAIS.index(c) for c in canais]
        np.savez_compressed(arquivo_modelo(modelo, cel), X=X[:, idx], canais=np.array(canais),
                            origem_patch=origem, **meta)
    return int(ok.sum())


def arquivo_modelo(modelo, cel):
    return PASTA_SAIDA / f'modelo_{modelo}' / f'patches_{cel}.npz'


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ee.Initialize(project=PROJETO_EE)
    PASTA_SAIDA.mkdir(parents=True, exist_ok=True)
    papeis = ee.ImageCollection(ASSET_FAV).filter(ee.Filter.eq('alvo', ALVO)).aggregate_array('papeis').distinct().getInfo()
    if 'None' in papeis:
        log.warning(f'{ASSET_FAV} foi treinado com todos os pontos (papeis = None): o teste não é '
                    'independente. Ver docstring.')

    print('\n  1) só sortear (índice + divisão, sem baixar)\n  2) sortear e baixar os patches')
    acao = input('Opção [1]: ').strip() or '1'

    grid = ler_grid()
    pimg = imagem_p()
    df, divisao = montar_amostras(grid, pimg)
    df.to_csv(PASTA_SAIDA / 'indice_patches.csv', index=False)
    (PASTA_SAIDA / 'divisao.json').write_text(json.dumps(divisao, indent=2, ensure_ascii=False))
    print('\n' + pd.crosstab(df.origem_amostra, df.conjunto, margins=True).to_string())
    for modelo, canais in MODELOS.items():
        print(f'modelo {modelo}: {len(df)} patches × {len(canais)} canais × {PATCH} × {PATCH} '
              f'≈ {len(df) * len(canais) * PATCH * PATCH * 2 / 1e9:.2f} GB em float16')
    if acao != '2' or input('Baixar agora? [s/N]: ').strip().lower() != 's':
        return

    img = imagem_canais()
    for modelo in MODELOS:
        (PASTA_SAIDA / f'modelo_{modelo}').mkdir(parents=True, exist_ok=True)
    for cel, g in df.groupby('celula'):
        if all(arquivo_modelo(m, cel).exists() for m in MODELOS):
            log.info(f'{cel}: já existe — pulado (apague para refazer)')
            continue
        t0 = time.time()
        n = baixar_celula(img, g, cel)
        log.info(f'{cel}: {n}/{len(g)} patches → modelo_A/ e modelo_B1/ ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
