"""
SWIR estimado (B04–B09) por Random Forest no GEE, aplicado por célula
da grid 5 × 5 aos mosaicos de 2009–2020 (sem SWIR), em dois GRUPOS:

  zero_cloud    mosaico das cenas zero-cloud 2009–2020 da célula
                (build_mosaicos_zero_cloud.py → PASTA_ZERO_CLOUD)
  semestral_v4  mosaico dos semestrais v4 da célula (cenas com nuvem filtradas
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
  menos de MIN_PONTOS_CENA pares usa a reta global. No notebook (validação uma
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
import build_mosaicos_semestrais as b   # noqa: E402

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CAATINGA         = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA'
DESTINO          = f'{CAATINGA}/ASTER_MOSAIC'
PASTA_ZERO_CLOUD = f'{CAATINGA}/ASTER_0cloud'      # ASTER_zero_cloud_2009_2020_<LxCy>
PASTA_V4         = f'{CAATINGA}/ASTER_v4'          # ASTER_QualityMosaic_<ano>_semestre_<s>_<LxCy>_v4
ASSET_PONTOS     = 'projects/mapbiomas-arida/mine/points_zeroCC_00_08'
# Conta do GEE: None = usa a conta ativa (~/.config/earthengine/current_user.txt)
# e o projeto dela; um nome (ex.: 'superconta') troca as credenciais antes.
# A conta precisa ter permissão de escrita em DESTINO.
CONTA            = None

PERIODO   = '2009_2020'
ANOS      = (2009, 2020)
GRUPOS    = ['zero_cloud', 'semestral_v4']
PREFIXO   = 'ASTER_MOSAIC_'
VERSAO    = 'swir_rf_v2_norm'

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

# normalização 2009–2020 → 2000–2008 (seção 10 do notebook)
CSV_PAREADOS    = (Path(__file__).resolve().parents[1] / 'Dados' / 'ASTER_Samples_0cloud'
                   / 'pareados' / 'points_zeroCC_pareados.csv')
ARQ_COEFS       = Path(__file__).resolve().parent / 'models' / 'zero_cloud' / 'normalizacao_ols_cena_coef.csv'
NORMALIZACAO    = {'zero_cloud': 'cena', 'semestral_v4': 'nenhuma'}   # 'cena', 'global' ou 'nenhuma'
MIN_PONTOS_CENA = 200
CODIGO_NORM     = {'nenhuma': 0, 'global': 1, 'cena': 2}


# ═══════════════════════════════════════════════════════════════════════════
# MOSAICOS DE ENTRADA (por célula)
# ═══════════════════════════════════════════════════════════════════════════
def mosaico_zero_cloud(cel: dict):
    img = ee.Image(f'{PASTA_ZERO_CLOUD}/ASTER_zero_cloud_{PERIODO}_{cel["nome"]}')
    return img, ee.String(img.get('cenas'))


def mosaico_semestral_v4(cel: dict):
    col = (ee.ImageCollection(PASTA_V4)
             .filter(ee.Filter.eq('celula_nome', cel['nome']))
             .filter(ee.Filter.rangeContains('year', ANOS[0], ANOS[1])))
    img = col.select(VNIR + TIR + ['quality', 'saturado']).qualityMosaic('quality')
    return img, col.aggregate_array('system:index').join(',')


ENTRADAS = {'zero_cloud': mosaico_zero_cloud, 'semestral_v4': mosaico_semestral_v4}


def existe_entrada(grupo: str, cel: dict) -> ee.Number:
    if grupo == 'zero_cloud':
        return (ee.ImageCollection(PASTA_ZERO_CLOUD)
                  .filter(ee.Filter.eq('system:index', f'ASTER_zero_cloud_{PERIODO}_{cel["nome"]}'))
                  .size())
    return (ee.ImageCollection(PASTA_V4)
              .filter(ee.Filter.eq('celula_nome', cel['nome']))
              .filter(ee.Filter.rangeContains('year', ANOS[0], ANOS[1])).size())


# ═══════════════════════════════════════════════════════════════════════════
# MODELO
# ═══════════════════════════════════════════════════════════════════════════
def params_por_banda() -> dict[str, dict]:
    if ARQ_PARAMS.exists():
        salvos = json.loads(ARQ_PARAMS.read_text()).get('rf', {})
        if all(s in salvos for s in SWIR):
            log.info(f'RF: params por banda de {ARQ_PARAMS}')
            return {s: {**salvos[s]['gee'], 'seed': RF_PARAMS['seed']} for s in SWIR}
    log.info(f'RF: params padrão {RF_PARAMS}')
    return {s: RF_PARAMS for s in SWIR}


def coeficientes() -> dict:
    """OLS por banda (2009–2020 → 2000–2008) nos pontos pareados: global e por cena.

    {'global': {banda: (offset, ganho)}, 'cenas': {data: {banda: (offset, ganho)}}},
    nas unidades do mosaico (reflectância × 1000, TIR em DN). Grava ARQ_COEFS.
    """
    if not CSV_PAREADOS.exists():
        raise SystemExit(f'{CSV_PAREADOS} não existe — rode building_dataset/baixar_pontos_pareados.py')
    df = pd.read_csv(CSV_PAREADOS).dropna(subset=[f'{x}_{p}' for x in FEATURES for p in ('00_08', '09_20')])

    def reta(d):
        return {x: tuple(np.polyfit(d[f'{x}_09_20'], d[f'{x}_00_08'], 1)[::-1]) for x in FEATURES}

    coefs = {'global': reta(df), 'cenas': {}}
    for data, d in df.groupby('data_09_20'):
        if len(d) >= MIN_PONTOS_CENA:
            coefs['cenas'][int(data)] = reta(d)
    linhas = [{'data': k, 'n_pontos': int((df['data_09_20'] == k).sum()) if k != 'global' else len(df),
               **{f'{x}_{c}': v[x][i] for x in FEATURES for i, c in enumerate(('offset', 'ganho'))}}
              for k, v in [('global', coefs['global']), *coefs['cenas'].items()]]
    ARQ_COEFS.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(linhas).to_csv(ARQ_COEFS, index=False)
    n_cenas = df['data_09_20'].nunique()
    log.info(f'Normalização OLS: {len(coefs["cenas"])} de {n_cenas} cenas com reta própria '
             f'(≥ {MIN_PONTOS_CENA} pares); demais usam a global → {ARQ_COEFS.name}')
    return coefs


def normalizar(img, coefs: dict, modo: str):
    """(features VNIR+TIR normalizadas, banda 'normalizacao')."""
    feats = img.select(FEATURES).toFloat()
    if modo == 'nenhuma':
        return feats, ee.Image.constant(0).toUint8().rename('normalizacao')
    glob = coefs['global']
    if modo == 'global':
        partes = [feats.select(x).multiply(glob[x][1]).add(glob[x][0]).rename(x) for x in FEATURES]
        return ee.Image.cat(partes), ee.Image.constant(1).toUint8().rename('normalizacao')
    # por cena: ganho/offset escolhidos pela banda 'data' (remap; cena sem reta → global)
    datas = sorted(coefs['cenas'])
    data = img.select('data').toInt32()
    partes = []
    for x in FEATURES:
        ganho = data.remap(datas, [coefs['cenas'][d][x][1] for d in datas], glob[x][1]).toFloat()
        offset = data.remap(datas, [coefs['cenas'][d][x][0] for d in datas], glob[x][0]).toFloat()
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


def estimar(img, modelos: dict, coefs: dict, modo: str):
    """Saída: VNIR + SWIR estimado + TIR (INT32) + 'normalizacao'; saturado mascarado.

    O RF recebe as features normalizadas; VNIR/TIR da saída são os medidos.
    """
    img = img.updateMask(img.select('saturado').eq(0))
    medidas = img.select(FEATURES).toFloat()
    mascara = medidas.mask().reduce(ee.Reducer.min())
    feats, flag = normalizar(img, coefs, modo)
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
        b.usar_conta(conta)
    else:
        ee.Initialize(project=b.get_project_from_account(conta))
    print(f'Conta GEE: {conta} (projeto {b.get_project_from_account(conta)})')
    grid = b.ler_grid()
    existentes = b.listar_destino(DESTINO, criar=True)
    ativas = {t['description'] for t in ee.data.getTaskList()
              if t['state'] in ('READY', 'RUNNING') and t.get('description', '').startswith(PREFIXO)}

    print(f'\n══════ SWIR estimado (RF) — mosaicos {PERIODO} por célula → {DESTINO} ══════')
    print('Normalização 2009–2020 → 2000–2008: ' + ', '.join(f'{g} = {m}' for g, m in NORMALIZACAO.items()))
    grupos = b.perguntar_opcao('Grupos:', [
        ('os dois (zero_cloud e semestral_v4)', GRUPOS),
        ('só zero_cloud', ['zero_cloud']),
        ('só semestral_v4', ['semestral_v4']),
    ])
    celulas = b.perguntar_celulas(grid)
    treino = b.perguntar_opcao('Treino do RF:', [
        ('global — todos os pontos 2000–2008', 'global'),
        (f'por célula — pontos da célula (≥ {MIN_PONTOS_CELULA}; senão das vizinhas)', 'celula'),
    ], padrao=1 if TREINO == 'global' else 2)
    dry_run = b.perguntar_opcao('Modo:', [
        ('só conferir entradas e pontos (dry-run, não exporta)', True),
        ('exportar os mosaicos', False),
    ])

    itens = [(g, c) for g in grupos for c in celulas]
    nomes = {(g, c): f'{PREFIXO}{g}_{PERIODO}_{grid[c]["nome"]}' for g, c in itens}
    pendentes = [i for i in itens if dry_run or (nomes[i] not in existentes and nomes[i] not in ativas)]
    print(f'{len(itens)} itens; {len(itens) - len(pendentes)} já salvos ou em processamento; '
          f'{len(pendentes)} a {"conferir" if dry_run else "exportar"}.')
    if not pendentes or not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Nada a fazer.')

    pontos = ee.FeatureCollection(ASSET_PONTOS)
    contagens = ee.Dictionary({
        'pontos': pontos.aggregate_histogram('celula_id'),
        'entradas': ee.List([existe_entrada(g, grid[c]) for g, c in pendentes]),
    }).getInfo()
    pts_cel = {int(k): int(v) for k, v in contagens['pontos'].items()}
    log.info(f'{ASSET_PONTOS}: {sum(pts_cel.values())} pontos em {len(pts_cel)} células')

    params = params_por_banda()
    coefs = coeficientes()
    modelo_global = treinar(pontos, params)
    registro = []
    for (g, c), n_ent in zip(pendentes, contagens['entradas']):
        cel = {'id': c, **grid[c]}
        nome = nomes[(g, c)]
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

        modelos = (modelo_global if modo_treino == 'global'
                   else treinar(pontos.filter(ee.Filter.inList('celula_id', ids_treino)), params))
        entrada, origem = ENTRADAS[g](cel)
        saida = (estimar(entrada, modelos, coefs, modo).clip(cel['geom'])
                 .set({'grupo': g, 'periodo': PERIODO, 'celula_id': c, 'celula_nome': cel['nome'],
                       'origem': origem, 'swir': 'estimado_rf', 'treino': linha['treino'],
                       'n_pontos_treino': linha['n_pontos_treino'], 'pontos_treino': ASSET_PONTOS,
                       'rf_params': str(params[SWIR[0]] if len({str(p) for p in params.values()}) == 1
                                         else params),
                       'bandas_estimadas': ','.join(SWIR), 'grid_asset': b.GRID_ASSET,
                       'SCALE_FACTORS': {**dict.fromkeys(VNIR + SWIR, b.ESCALA_TOA),
                                         **dict.fromkeys(TIR, 1)},
                       'normalizacao': modo,
                       'normalizacao_n_cenas_reta_propria': len(coefs['cenas']) if modo == 'cena' else 0,
                       'normalizacao_coefs': ARQ_COEFS.name,
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
