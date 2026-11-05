"""
VERSÃO 2 (swir_rf_v3_seco) — mesmas regras de estimar_swir_rf.py, mas as entradas
são os mosaicos da estação seca: ASTER_0cloud_v2 (build_mosaicos_zero_cloud_v2.py)
e ASTER_v5 (build_mosaicos_semestrais_v2.py), e a saída vai para
projects/mapbiomas-arida/mine/ASTER_MOSAIC_v2. O grupo dos semestrais passa a se
chamar 'semestral_v5'. Treino (ASSET_PONTOS) e coeficientes de normalização por
cena são os mesmos (as cenas secas são um subconjunto das cenas do v1).

Grupo 'sazonal_seco' (build_mosaico_sazonal_seco.py): entrada é o mosaico seco
quality do período B (ASTER_sazonal_B_seco_quality_<LxCy>, 2009–2020) e o RF
treina nos pontos do mosaico seco do período A (pontos_sazonal_seco_pareados.py,
etapa 1 → PONTOS['sazonal_seco']). O mosaico não tem banda 'data' (cada pixel
vem do qualityMosaic de várias cenas), então não há reta por cena; as opções são
'nenhuma', 'global' ou 'celula' — reta OLS B → A por CÉLULA nos pontos pareados
(points_sazonal_seco_pareados.csv), célula com < MIN_PONTOS_RETA pares usa a
global. A normalização desse grupo é perguntada na execução.

Grupo 'sazonal_pareado': mesma entrada do sazonal_seco, mas o RF treina nos
PARES (join por pid de points_sazonal_seco_B com points_sazonal_seco_A):
features VNIR + TIR do mosaico B, alvo o SWIR do mosaico A no mesmo ponto.
O modelo aprende direto a diferença radiométrica B × A, então não normaliza
('nenhuma'). Ponto do A sem dado no B fica de fora. Saída com o nome do grupo
(ASTER_MOSAIC_sazonal_pareado_2009_2020_<LxCy>), para comparar com o sazonal_seco.

Fonte dos grupos sazonais: pontos_sazonal_seco_pareados.FONTE ('hm_razao' por
padrão = mosaicos harmonizados contra o Landsat). Mosaico de entrada, pontos,
hiperparâmetros (models/sazonal_seco<SUFIXO>/) e nome da saída levam o sufixo da
fonte: ASTER_MOSAIC_sazonal_seco_hm_razao_2009_2020_<LxCy>. Com FONTE = 'original'
o sufixo é vazio e tudo fica como antes.

Texto original (v1):
SWIR estimado (B04–B09) por Random Forest no GEE, aplicado por célula
da grid 5 × 5 aos mosaicos de 2009–2020 (sem SWIR), em dois GRUPOS:

  zero_cloud    mosaico das cenas zero-cloud 2009–2020 da célula
                (build_mosaicos_zero_cloud.py → PASTA_ZERO_CLOUD)
  semestral_v5  mosaico dos semestrais v4 da célula (cenas com nuvem filtradas
                por CC < 77 e quality), 2009–2020 compostos por qualityMosaic
                — o pixel de maior 'quality' entre os semestres (PASTA_V4)

Treino
  Pontos de ASSET_PONTOS (pontos_zero_cloud_pareados.py, etapa 1): VNIR + SWIR
  + TIR do mosaico zero-cloud 2000–2008. Um smileRandomForest (REGRESSION) por banda
  SWIR, features VNIR + TIR (B01, B02, B3N, B10–B14), alvo a banda SWIR.
  TREINO = 'global' → todos os pontos; 'celula' → só os pontos da célula (se
  tiver ≥ MIN_PONTOS_CELULA). Célula sem pontos suficientes usa os pontos das
  células VIZINHAS (anel 1, e se faltar o anel 2…) — mais perto e mais
  parecidas que o bloco todo.
  Todo treino é limitado a N_MAX_TREINO pontos (amostra aleatória): o RF com
  os 34 mil pontos e minLeafPopulation 2 passava do limite do GEE
  ("Image.classify: Computed value is too large").

Normalização radiométrica 2009–2020 → 2000–2008 (train_swir_models.ipynb, seção 10)
  Antes do RF, cada banda VNIR/TIR do mosaico 2009–2020 é levada à escala de
  2000–2008: valor_norm = offset + ganho × valor. Ganho e offset por banda vêm de
  uma reta de mínimos quadrados (OLS) nos pontos pareados (CSV_PAREADOS),
  ajustada POR CENA (data AAAAMMDD da cena 2009–2020 que deu o pixel). Cena com
  menos de MIN_PONTOS_RETA pares usa a reta global. No notebook (validação uma
  célula de fora) o R² do SWIR estimado no teste temporal foi de −0,44 para
  +0,49, mantendo 74–85% do desvio padrão das bandas.
  NORMALIZACAO por grupo: 'cena' precisa da banda 'data' (zero-cloud tem);
  os semestrais v4 não têm → 'nenhuma' (ou 'global', que achata a imagem).
  O RF recebe as bandas normalizadas; a saída guarda VNIR/TIR MEDIDOS (sem
  normalização) e o SWIR estimado.

Saída (ImageCollection DESTINO, conta ativa ou CONTA)
  ASTER_MOSAIC_<grupo>_2009_2020_<LxCy>, ex.: ASTER_MOSAIC_zero_cloud_2009_2020_L3C3
  Bandas: B01, B02, B3N (medidas), B04–B09 (ESTIMADAS), B10–B14 (medidas) —
  INT32, reflectância × 1000, TIR em DN; banda 'normalizacao' (UINT8: 0 nenhuma,
  1 reta global, 2 reta da cena). Pixels saturados (VNIR DN ≥ 254; nas cenas
  zero-cloud eram nuvem) ficam mascarados.
  Propriedades: grupo, periodo, celula_id, celula_nome, origem, swir
  ('estimado_rf'), treino, n_pontos_treino, rf_params, SCALE_FACTORS…

Hiperparâmetros: os do RandomizedSearchCV do train_swir_models.ipynb (seção 5,
RF ≡ smileRandomForest), lidos de models/zero_cloud/best_params.json
(chave rf.<banda>.gee) se existir; senão RF_PARAMS abaixo (o resultado da
última execução do notebook — iguais para as 6 bandas).

Sem normalização, o teste temporal (seção 8) deu R² negativo: o TIR em DN
carrega a temperatura da data da cena. Ver modeling/README.md, seções 4.2–4.3.

Uso: python estimar_swir_rf.py
"""

import csv
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import ee
import numpy as np
import pandas as pd

PASTA_BUILD = Path(__file__).resolve().parents[1] / 'building_dataset'
sys.path.append(str(PASTA_BUILD))
import build_mosaicos_semestrais_v2 as b   # noqa: E402
import pontos_sazonal_seco_pareados as ps   # noqa: E402

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CAATINGA         = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA'
DESTINO          = 'projects/mapbiomas-arida/mine/ASTER_MOSAIC_v2'
PASTA_ZERO_CLOUD = 'projects/mapbiomas-arida/mine/ASTER_0cloud_v2'      # ASTER_zero_cloud_2009_2020_<LxCy>
PASTA_V4         = 'projects/mapbiomas-arida/mine/ASTER_v5'             # ASTER_QualityMosaic_<ano>_semestre_2_<LxCy>_v5
ASSET_PONTOS     = 'projects/mapbiomas-arida/mine/points_zeroCC_00_08'
PASTA_SAZONAL    = ps.COL_FONTE          # ASTER_sazonal_B_seco_quality_<LxCy> (original ou harmonizado)
# pontos de treino (VNIR + SWIR + TIR de 2000–2008) de cada grupo
# PONTOS_PAREADO não é um asset: é o join A × B por pid (pontos_treino())
PONTOS_PAREADO   = f'{ps.ASSET_A}+{ps.ASSET_B.rsplit("/", 1)[-1]}'
PONTOS = {'zero_cloud': ASSET_PONTOS, 'semestral_v5': ASSET_PONTOS, 'sazonal_seco': ps.ASSET_A,
          'sazonal_pareado': PONTOS_PAREADO}
# Conta do GEE: None = usa a conta ativa (~/.config/earthengine/current_user.txt)
# e o projeto dela; um nome (ex.: 'superconta') troca as credenciais antes.
# A conta precisa ter permissão de escrita em DESTINO.
CONTA            = None
# Projeto do GEE onde as tasks rodam (cota de processamento). None = projeto da conta.
# 'mapbiomas-brazil' com caatinga04: fila mais rápida; o destino continua em DESTINO.
# A lista de tasks ativas é por projeto: tasks enviadas em outro projeto não aparecem
# aqui e seriam enviadas de novo — cancele-as antes de trocar.
PROJETO_EE       = 'mapbiomas-brazil'

PERIODO   = '2009_2020'
ANOS      = (2009, 2020)
GRUPOS    = ['zero_cloud', 'semestral_v5', 'sazonal_seco', 'sazonal_pareado']
PREFIXO   = 'ASTER_MOSAIC_'
VERSAO    = 'swir_rf_v3_seco'

VNIR     = ['B01', 'B02', 'B3N']
SWIR     = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']
TIR      = ['B10', 'B11', 'B12', 'B13', 'B14']
FEATURES = VNIR + TIR
SAIDA    = VNIR + SWIR + TIR

TREINO            = 'global'   # 'global' ou 'celula'
MIN_PONTOS_CELULA = 500
N_MAX_TREINO      = 6000       # teto de pontos por RF (2000 por célula funcionou; 34000 estourou)
MAX_ANEL          = 3          # até quantos anéis de vizinhas procurar
# train_swir_models.ipynb (RF, R² teste 0,94–0,96): mesmos para as 6 bandas
RF_PARAMS = {'numberOfTrees': 108, 'variablesPerSplit': 5, 'minLeafPopulation': 2,
             'bagFraction': 0.983, 'maxNodes': None, 'seed': 42}
ARQ_PARAMS = Path(__file__).resolve().parent / 'models' / 'zero_cloud' / 'best_params.json'
# saídas do train_swir_models_sazonal.ipynb da mesma fonte (models/sazonal_seco_hm_razao, …)
DIR_SAZONAL = ARQ_PARAMS.parents[1] / f'sazonal_seco{ps.SUFIXO}'
# hiperparâmetros por grupo (train_swir_models*.ipynb); arquivo ausente → RF_PARAMS
ARQ_PARAMS_GRUPO = {
    'zero_cloud':      ARQ_PARAMS,
    'semestral_v5':    ARQ_PARAMS,
    'sazonal_seco':    DIR_SAZONAL / 'best_params.json',
    'sazonal_pareado': DIR_SAZONAL / 'best_params_pareado.json',
}

# normalização 2009–2020 → 2000–2008 (seção 10 do notebook)
CSV_PAREADOS    = (Path(__file__).resolve().parents[1] / 'Dados' / 'ASTER_Samples_0cloud'
                   / 'pareados' / 'points_zeroCC_pareados.csv')
ARQ_COEFS       = Path(__file__).resolve().parent / 'models' / 'zero_cloud' / 'normalizacao_ols_cena_coef.csv'
# 'cena' (precisa da banda 'data'), 'celula', 'global' ou 'nenhuma'; sazonal_seco é perguntado
NORMALIZACAO    = {'zero_cloud': 'cena', 'semestral_v5': 'nenhuma', 'sazonal_seco': 'nenhuma',
                   'sazonal_pareado': 'nenhuma'}
MIN_PONTOS_RETA = 200          # pares mínimos para uma cena/célula ter reta própria
CODIGO_NORM     = {'nenhuma': 0, 'global': 1, 'cena': 2, 'celula': 3}
# pontos pareados de cada grupo: sufixos das colunas (referência 2000–2008, alvo
# 2009–2020) e a coluna que separa as retas próprias (cena ou célula)
PAREADOS = {
    'zero_cloud':   {'csv': CSV_PAREADOS, 'ref': '_00_08', 'alvo': '_09_20',
                     'chave': 'data_09_20', 'arq': ARQ_COEFS},
    'sazonal_seco': {'csv': ps.SAIDA / 'points_sazonal_seco_pareados.csv', 'ref': '_A', 'alvo': '_B',
                     'chave': 'celula_id',
                     'arq': DIR_SAZONAL / 'normalizacao_ols_celula_coef.csv'},
}
PAREADOS['semestral_v5'] = PAREADOS['zero_cloud']


# ═══════════════════════════════════════════════════════════════════════════
# MOSAICOS DE ENTRADA (por célula)
# ═══════════════════════════════════════════════════════════════════════════
def mosaico_zero_cloud(cel: dict):
    img = ee.Image(f'{PASTA_ZERO_CLOUD}/ASTER_zero_cloud_{PERIODO}_{cel["nome"]}')
    return img, ee.String(img.get('cenas'))


def mosaico_semestral_v5(cel: dict):
    col = (ee.ImageCollection(PASTA_V4)
             .filter(ee.Filter.eq('celula_nome', cel['nome']))
             .filter(ee.Filter.rangeContains('year', ANOS[0], ANOS[1])))
    img = col.select(VNIR + TIR + ['quality', 'saturado']).qualityMosaic('quality')
    return img, col.aggregate_array('system:index').join(',')


def nome_saida(grupo: str, celula: str) -> str:
    """Grupos sazonais levam o sufixo da fonte (ex.: _hm_razao); os outros não."""
    suf = ps.SUFIXO if grupo.startswith('sazonal') else ''
    return f'{PREFIXO}{grupo}{suf}_{PERIODO}_{celula}'


def nome_sazonal(cel: dict) -> str:
    return ps.s.nome_mosaico('B', cel['nome'])


def mosaico_sazonal_seco(cel: dict):
    """Mosaico seco quality 2009–2020; saturado já saiu no build (banda 0 para estimar())."""
    img = ee.Image(f'{PASTA_SAZONAL}/{nome_sazonal(cel)}')
    img = img.select(VNIR + TIR).addBands(ee.Image.constant(0).rename('saturado'))
    return img, ee.String(nome_sazonal(cel))


ENTRADAS = {'zero_cloud': mosaico_zero_cloud, 'semestral_v5': mosaico_semestral_v5,
            'sazonal_seco': mosaico_sazonal_seco, 'sazonal_pareado': mosaico_sazonal_seco}


def existe_entrada(grupo: str, cel: dict) -> ee.Number:
    if grupo == 'zero_cloud':
        return (ee.ImageCollection(PASTA_ZERO_CLOUD)
                  .filter(ee.Filter.eq('system:index', f'ASTER_zero_cloud_{PERIODO}_{cel["nome"]}'))
                  .size())
    if grupo in ('sazonal_seco', 'sazonal_pareado'):
        return (ee.ImageCollection(PASTA_SAZONAL)
                  .filter(ee.Filter.eq('system:index', nome_sazonal(cel))).size())
    return (ee.ImageCollection(PASTA_V4)
              .filter(ee.Filter.eq('celula_nome', cel['nome']))
              .filter(ee.Filter.rangeContains('year', ANOS[0], ANOS[1])).size())


# ═══════════════════════════════════════════════════════════════════════════
# MODELO
# ═══════════════════════════════════════════════════════════════════════════
def pontos_treino(asset: str) -> ee.FeatureCollection:
    """Pontos de treino de PONTOS[grupo]; PONTOS_PAREADO = features do B + SWIR do A (join por pid)."""
    if asset != PONTOS_PAREADO:
        return ee.FeatureCollection(asset)
    a = ee.FeatureCollection(ps.ASSET_A).select(['pid'] + SWIR)
    bb = ee.FeatureCollection(ps.ASSET_B).select(ps.PROPS_PONTO + FEATURES)
    pares = ee.Join.inner().apply(bb, a, ee.Filter.equals(leftField='pid', rightField='pid'))
    return ee.FeatureCollection(pares.map(
        lambda f: ee.Feature(f.get('primary')).copyProperties(ee.Feature(f.get('secondary')), SWIR)))


def params_por_banda(grupo: str) -> dict[str, dict]:
    arq = ARQ_PARAMS_GRUPO[grupo]
    if arq.exists():
        salvos = json.loads(arq.read_text()).get('rf', {})
        if all(s in salvos for s in SWIR):
            log.info(f'RF ({grupo}): params por banda de {arq}')
            return {s: {**salvos[s]['gee'], 'seed': RF_PARAMS['seed']} for s in SWIR}
    log.info(f'RF ({grupo}): {arq.name} não encontrado — params padrão {RF_PARAMS}')
    return {s: RF_PARAMS for s in SWIR}


def coeficientes(grupo: str) -> dict:
    """OLS por banda (2009–2020 → 2000–2008) nos pontos pareados do grupo: global e
    por cena (zero-cloud) ou por célula (sazonal_seco).

    {'global': {banda: (offset, ganho)}, 'grupos': {cena|celula: {banda: (offset, ganho)}}},
    nas unidades do mosaico (reflectância × 1000, TIR em DN). Grava o CSV de PAREADOS[grupo].
    """
    cfg = PAREADOS[grupo]
    if not cfg['csv'].exists():
        raise SystemExit(f'{cfg["csv"]} não existe — baixe os pontos pareados '
                         '(baixar_pontos_pareados.py ou pontos_sazonal_seco_pareados.py → baixar)')
    ref, alvo, chave = cfg['ref'], cfg['alvo'], cfg['chave']
    df = pd.read_csv(cfg['csv']).dropna(subset=[f'{x}{p}' for x in FEATURES for p in (ref, alvo)])

    def reta(d):
        return {x: tuple(np.polyfit(d[f'{x}{alvo}'], d[f'{x}{ref}'], 1)[::-1]) for x in FEATURES}

    coefs = {'global': reta(df), 'grupos': {}}
    for k, d in df.groupby(chave):
        if len(d) >= MIN_PONTOS_RETA:
            coefs['grupos'][int(k)] = reta(d)
    linhas = [{chave: k, 'n_pontos': int((df[chave] == k).sum()) if k != 'global' else len(df),
               **{f'{x}_{c}': v[x][i] for x in FEATURES for i, c in enumerate(('offset', 'ganho'))}}
              for k, v in [('global', coefs['global']), *coefs['grupos'].items()]]
    cfg['arq'].parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(linhas).to_csv(cfg['arq'], index=False)
    log.info(f'Normalização OLS ({grupo}): {len(coefs["grupos"])} de {df[chave].nunique()} '
             f'{chave} com reta própria (≥ {MIN_PONTOS_RETA} pares); demais usam a global '
             f'→ {cfg["arq"].name}')
    return coefs


def normalizar(img, coefs: dict, modo: str, celula_id: int | None = None):
    """(features VNIR+TIR normalizadas, banda 'normalizacao')."""
    feats = img.select(FEATURES).toFloat()
    if modo == 'nenhuma':
        return feats, ee.Image.constant(0).toUint8().rename('normalizacao')
    glob = coefs['global']
    if modo == 'global':
        partes = [feats.select(x).multiply(glob[x][1]).add(glob[x][0]).rename(x) for x in FEATURES]
        return ee.Image.cat(partes), ee.Image.constant(1).toUint8().rename('normalizacao')
    if modo == 'celula':
        # reta da célula (constantes); célula sem reta própria → global
        r = coefs['grupos'].get(celula_id, glob)
        partes = [feats.select(x).multiply(r[x][1]).add(r[x][0]).rename(x) for x in FEATURES]
        codigo = CODIGO_NORM['celula'] if celula_id in coefs['grupos'] else CODIGO_NORM['global']
        return ee.Image.cat(partes), ee.Image.constant(codigo).toUint8().rename('normalizacao')
    # por cena: ganho/offset escolhidos pela banda 'data' (remap; cena sem reta → global)
    datas = sorted(coefs['grupos'])
    data = img.select('data').toInt32()
    partes = []
    for x in FEATURES:
        ganho = data.remap(datas, [coefs['grupos'][d][x][1] for d in datas], glob[x][1]).toFloat()
        offset = data.remap(datas, [coefs['grupos'][d][x][0] for d in datas], glob[x][0]).toFloat()
        partes.append(feats.select(x).multiply(ganho).add(offset).rename(x))
    flag = data.remap(datas, [2] * len(datas), 1).toUint8().rename('normalizacao')
    return ee.Image.cat(partes), flag


def vizinhas(c: int, grid: dict, pts_cel: dict) -> list[int]:
    """Células do menor anel em volta de c que somam ≥ MIN_PONTOS_CELULA pontos."""
    lc = {cid: tuple(map(int, cel['nome'][1:].split('C'))) for cid, cel in grid.items()}
    l0, c0 = lc[c]
    for anel in range(1, MAX_ANEL + 1):
        ids = [cid for cid, (l, cc) in lc.items()
               if cid != c and max(abs(l - l0), abs(cc - c0)) <= anel and pts_cel.get(cid, 0) > 0]
        if sum(pts_cel.get(i, 0) for i in ids) >= MIN_PONTOS_CELULA:
            return ids
    return []


def limitar(pontos: ee.FeatureCollection) -> ee.FeatureCollection:
    return pontos.randomColumn('_r', 42).sort('_r').limit(N_MAX_TREINO)


def treinar(pontos: ee.FeatureCollection, params: dict[str, dict]) -> dict:
    pontos = limitar(pontos)
    return {s: (ee.Classifier.smileRandomForest(**params[s])
                  .setOutputMode('REGRESSION')
                  .train(features=pontos, classProperty=s, inputProperties=FEATURES))
            for s in SWIR}


def estimar(img, modelos: dict, coefs: dict, modo: str, celula_id: int | None = None):
    """Saída: VNIR + SWIR estimado + TIR (INT32) + 'normalizacao'; saturado mascarado.

    O RF recebe as features normalizadas; VNIR/TIR da saída são os medidos.
    """
    img = img.updateMask(img.select('saturado').eq(0))
    medidas = img.select(FEATURES).toFloat()
    mascara = medidas.mask().reduce(ee.Reducer.min())
    feats, flag = normalizar(img, coefs, modo, celula_id)
    feats = feats.updateMask(mascara)
    swir = [feats.classify(modelos[s], s).round().toInt32().rename(s) for s in SWIR]
    medidas = medidas.updateMask(mascara).round().toInt32()
    return ee.Image.cat([medidas.select(VNIR), *swir, medidas.select(TIR),
                         flag.updateMask(mascara)]).select(SAIDA + ['normalizacao'])


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ativa = b.ARQ_CONTA_ATUAL.read_text().strip()
    conta = CONTA or ativa
    if conta != ativa:
        b.switch_user(conta)
    projeto = PROJETO_EE or b.get_project_from_account(conta)
    ee.Initialize(project=projeto)
    print(f'Conta GEE: {conta} (projeto {projeto})')
    grid = b.ler_grid()
    existentes = b.listar_destino(DESTINO, criar=True)
    ativas = {t['description'] for t in ee.data.getTaskList()
              if t['state'] in ('READY', 'RUNNING') and t.get('description', '').startswith(PREFIXO)}

    print(f'\n══════ SWIR estimado (RF) — mosaicos {PERIODO} por célula → {DESTINO} ══════')
    print(f'Fonte dos grupos sazonais: {ps.FONTE} ({PASTA_SAZONAL})')
    print('Normalização 2009–2020 → 2000–2008: ' + ', '.join(f'{g} = {m}' for g, m in NORMALIZACAO.items()))
    grupos = b.perguntar_opcao('Grupos:', [
        ('sazonal_seco (mosaico seco quality 2009–2020)', ['sazonal_seco']),
        ('sazonal_pareado (mesmo mosaico; RF nos pares B → SWIR do A)', ['sazonal_pareado']),
        ('sazonal_seco e sazonal_pareado (comparar)', ['sazonal_seco', 'sazonal_pareado']),
        ('zero_cloud e semestral_v5', ['zero_cloud', 'semestral_v5']),
        ('só zero_cloud', ['zero_cloud']),
        ('só semestral_v5', ['semestral_v5']),
        ('todos', GRUPOS),
    ])
    if 'sazonal_seco' in grupos:
        NORMALIZACAO['sazonal_seco'] = b.perguntar_opcao('Normalização B → A do sazonal_seco:', [
            ('nenhuma', 'nenhuma'),
            ('reta por célula (pontos pareados; sem reta própria → global)', 'celula'),
            ('reta global', 'global'),
        ])
    celulas = b.perguntar_celulas(grid)
    treino = b.perguntar_opcao('Treino do RF:', [
        ('global — todos os pontos 2000–2008 do grupo', 'global'),
        (f'por célula — pontos da célula (≥ {MIN_PONTOS_CELULA}; senão das vizinhas)', 'celula'),
    ], padrao=1 if TREINO == 'global' else 2)
    dry_run = b.perguntar_opcao('Modo:', [
        ('só conferir entradas e pontos (dry-run, não exporta)', True),
        ('exportar os mosaicos', False),
    ])

    itens = [(g, c) for g in grupos for c in celulas]
    nomes = {(g, c): nome_saida(g, grid[c]['nome']) for g, c in itens}
    pendentes = [i for i in itens if dry_run or (nomes[i] not in existentes and nomes[i] not in ativas)]
    print(f'{len(itens)} itens; {len(itens) - len(pendentes)} já salvos ou em processamento; '
          f'{len(pendentes)} a {"conferir" if dry_run else "exportar"}.')
    if not pendentes or not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Nada a fazer.')

    assets_pts = sorted({PONTOS[g] for g in grupos})
    contagens = ee.Dictionary({
        'pontos': ee.List([pontos_treino(a).aggregate_histogram('celula_id') for a in assets_pts]),
        'entradas': ee.List([existe_entrada(g, grid[c]) for g, c in pendentes]),
    }).getInfo()
    pts_asset = {a: {int(k): int(v) for k, v in h.items()}
                 for a, h in zip(assets_pts, contagens['pontos'])}
    for a, pc in pts_asset.items():
        log.info(f'{a}: {sum(pc.values())} pontos em {len(pc)} células')

    params_grupo = {g: params_por_banda(g) for g in grupos}
    # retas só dos grupos que normalizam (o CSV pareado do grupo precisa existir)
    coefs_grupo = {g: coeficientes(g) if NORMALIZACAO[g] != 'nenhuma' else {} for g in grupos}
    # um modelo global por (asset de pontos, grupo): grupos podem ter params diferentes
    modelo_global = {(PONTOS[g], g): treinar(pontos_treino(PONTOS[g]), params_grupo[g]) for g in grupos}
    registro = []
    for (g, c), n_ent in zip(pendentes, contagens['entradas']):
        cel = {'id': c, **grid[c]}
        nome = nomes[(g, c)]
        asset_pts = PONTOS[g]
        pontos, pts_cel = pontos_treino(asset_pts), pts_asset[asset_pts]
        n_pts = pts_cel.get(c, 0)
        por_cel = treino == 'celula' and n_pts >= MIN_PONTOS_CELULA
        viz = vizinhas(c, grid, pts_cel) if treino == 'celula' and not por_cel else []
        if por_cel:
            modo_treino, ids_treino = 'celula', [c]
        elif viz:
            modo_treino, ids_treino = 'vizinhas:' + '+'.join(grid[i]['nome'] for i in viz), viz
        else:
            modo_treino, ids_treino = 'global', list(pts_cel)
        n_treino = min(N_MAX_TREINO, sum(pts_cel.get(i, 0) for i in ids_treino))
        modo = NORMALIZACAO[g]
        linha = {'grupo': g, 'celula_id': c, 'celula_nome': cel['nome'], 'entradas': n_ent,
                 'normalizacao': modo,
                 'treino': modo_treino, 'n_pontos_treino': n_treino}
        if n_ent == 0:
            log.info(f'{nome}: sem mosaico de entrada — pulado')
            registro.append({**linha, 'status': 'sem_entrada'})
            continue
        if dry_run:
            log.info(f'{nome}: {n_ent} imagem(ns) de entrada, treino {linha["treino"]} '
                     f'({linha["n_pontos_treino"]} pontos), normalização {modo}')
            registro.append({**linha, 'status': 'dry_run'})
            continue

        params = params_grupo[g]
        modelos = (modelo_global[(asset_pts, g)] if modo_treino == 'global'
                   else treinar(pontos.filter(ee.Filter.inList('celula_id', ids_treino)), params))
        entrada, origem = ENTRADAS[g](cel)
        coefs = coefs_grupo[g]
        saida = (estimar(entrada, modelos, coefs, modo, c).clip(cel['geom'])
                 .set({'grupo': g, 'periodo': PERIODO,
                       'fonte': ps.FONTE if g.startswith('sazonal') else '', 'celula_id': c, 'celula_nome': cel['nome'],
                       'origem': origem, 'swir': 'estimado_rf', 'treino': linha['treino'],
                       'n_pontos_treino': linha['n_pontos_treino'], 'pontos_treino': asset_pts,
                       'rf_params': str(params[SWIR[0]] if len({str(p) for p in params.values()}) == 1
                                         else params),
                       'bandas_estimadas': ','.join(SWIR), 'grid_asset': b.GRID_ASSET,
                       'SCALE_FACTORS': {**dict.fromkeys(VNIR + SWIR, b.ESCALA_TOA),
                                         **dict.fromkeys(TIR, 1)},
                       'normalizacao': modo,
                       'normalizacao_n_retas_proprias': len(coefs.get('grupos', {})),
                       'normalizacao_coefs': PAREADOS[g]['arq'].name if modo != 'nenhuma' else '',
                       'DATA_TYPE': 'INT32', 'versao': VERSAO}))
        task = ee.batch.Export.image.toAsset(
            image=ee.Image(saida), description=nome, assetId=f'{DESTINO}/{nome}',
            region=cel['geom'], crs=b.CRS_EXPORT, crsTransform=b.CRS_TRANSFORM,
            maxPixels=1e13, pyramidingPolicy={'.default': 'mean', 'normalizacao': 'mode'})
        try:
            task.start()
        except ee.EEException as e:
            log.error(f'{nome}: falhou ao iniciar: {e}')
            registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}'})
            continue
        log.info(f'{nome}: task {task.id} (treino {linha["treino"]}, normalização {modo})')
        registro.append({**linha, 'status': 'exportando', 'task_id': task.id})

    arq = Path(__file__).resolve().parent / f'registro_swir_rf_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['grupo', 'celula_id', 'celula_nome', 'entradas', 'normalizacao', 'treino',
                                          'n_pontos_treino', 'status', 'task_id'], extrasaction='ignore')
        w.writeheader()
        w.writerows(registro)
    log.info(f'Registro → {arq}')


if __name__ == '__main__':
    main()
