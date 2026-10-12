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
  tiver ≥ MIN_PONTOS_CELULA; senão cai no global).

Saída (ImageCollection DESTINO, conta superconta)
  ASTER_MOSAIC_<grupo>_2009_2020_<LxCy>, ex.: ASTER_MOSAIC_zero_cloud_2009_2020_L3C3
  Bandas: B01, B02, B3N (medidas), B04–B09 (ESTIMADAS), B10–B14 (medidas) —
  INT32, reflectância × 1000, TIR em DN. Pixels saturados (VNIR DN ≥ 254; nas
  cenas zero-cloud eram nuvem) ficam mascarados.
  Propriedades: grupo, periodo, celula_id, celula_nome, origem, swir
  ('estimado_rf'), treino, n_pontos_treino, rf_params, SCALE_FACTORS…

Hiperparâmetros: os do RandomizedSearchCV do train_swir_models.ipynb (seção 5,
RF ≡ smileRandomForest), lidos de models/zero_cloud/best_params.json
(chave rf.<banda>.gee) se existir; senão RF_PARAMS abaixo (o resultado da
última execução do notebook — iguais para as 6 bandas).

ATENÇÃO: no teste temporal (train_swir_models.ipynb, seção 8) os modelos
treinados em 2000–2008 tiveram R² negativo ao receber VNIR + TIR de 2009–2020 —
o TIR em DN carrega a temperatura da data da cena. Ver modeling/README.md.

Uso: python estimar_swir_rf.py
"""

import csv
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import ee

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
CONTA            = 'superconta'                    # escreve em mapbiomas-brazil

PERIODO   = '2009_2020'
ANOS      = (2009, 2020)
GRUPOS    = ['zero_cloud', 'semestral_v4']
PREFIXO   = 'ASTER_MOSAIC_'
VERSAO    = 'swir_rf_v1'

VNIR     = ['B01', 'B02', 'B3N']
SWIR     = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']
TIR      = ['B10', 'B11', 'B12', 'B13', 'B14']
FEATURES = VNIR + TIR
SAIDA    = VNIR + SWIR + TIR

TREINO            = 'global'   # 'global' ou 'celula'
MIN_PONTOS_CELULA = 500
# train_swir_models.ipynb (RF, R² teste 0,94–0,96): mesmos para as 6 bandas
RF_PARAMS = {'numberOfTrees': 108, 'variablesPerSplit': 5, 'minLeafPopulation': 2,
             'bagFraction': 0.983, 'maxNodes': None, 'seed': 42}
ARQ_PARAMS = Path(__file__).resolve().parent / 'models' / 'zero_cloud' / 'best_params.json'


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


def treinar(pontos: ee.FeatureCollection, params: dict[str, dict]) -> dict:
    return {s: (ee.Classifier.smileRandomForest(**params[s])
                  .setOutputMode('REGRESSION')
                  .train(features=pontos, classProperty=s, inputProperties=FEATURES))
            for s in SWIR}


def estimar(img, modelos: dict):
    """Mosaico de saída: VNIR + SWIR estimado + TIR (INT32), saturado mascarado."""
    img = img.updateMask(img.select('saturado').eq(0))
    feats = img.select(FEATURES).toFloat()
    feats = feats.updateMask(feats.mask().reduce(ee.Reducer.min()))
    swir = [feats.classify(modelos[s], s).round().toInt32().rename(s) for s in SWIR]
    medidas = feats.round().toInt32()
    return ee.Image.cat([medidas.select(VNIR), *swir, medidas.select(TIR)]).select(SAIDA)


# ═══════════════════════════════════════════════════════════════════════════
def main():
    if b.ARQ_CONTA_ATUAL.read_text().strip() != CONTA:
        b.usar_conta(CONTA)
    else:
        ee.Initialize(project=b.get_project_from_account(CONTA))
    grid = b.ler_grid()
    existentes = b.listar_destino(DESTINO, criar=True)
    ativas = {t['description'] for t in ee.data.getTaskList()
              if t['state'] in ('READY', 'RUNNING') and t.get('description', '').startswith(PREFIXO)}

    print(f'\n══════ SWIR estimado (RF) — mosaicos {PERIODO} por célula → {DESTINO} ══════')
    grupos = b.perguntar_opcao('Grupos:', [
        ('os dois (zero_cloud e semestral_v4)', GRUPOS),
        ('só zero_cloud', ['zero_cloud']),
        ('só semestral_v4', ['semestral_v4']),
    ])
    celulas = b.perguntar_celulas(grid)
    treino = b.perguntar_opcao('Treino do RF:', [
        ('global — todos os pontos 2000–2008', 'global'),
        (f'por célula — pontos da célula (≥ {MIN_PONTOS_CELULA}; senão global)', 'celula'),
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
    modelo_global = treinar(pontos, params)
    registro = []
    for (g, c), n_ent in zip(pendentes, contagens['entradas']):
        cel = {'id': c, **grid[c]}
        nome = nomes[(g, c)]
        n_pts = pts_cel.get(c, 0)
        por_cel = treino == 'celula' and n_pts >= MIN_PONTOS_CELULA
        linha = {'grupo': g, 'celula_id': c, 'celula_nome': cel['nome'], 'entradas': n_ent,
                 'treino': 'celula' if por_cel else 'global',
                 'n_pontos_treino': n_pts if por_cel else sum(pts_cel.values())}
        if n_ent == 0:
            log.info(f'{nome}: sem mosaico de entrada — pulado')
            registro.append({**linha, 'status': 'sem_entrada'})
            continue
        if dry_run:
            log.info(f'{nome}: {n_ent} imagem(ns) de entrada, treino {linha["treino"]} '
                     f'({linha["n_pontos_treino"]} pontos)')
            registro.append({**linha, 'status': 'dry_run'})
            continue

        modelos = (treinar(pontos.filter(ee.Filter.eq('celula_id', c)), params) if por_cel
                   else modelo_global)
        entrada, origem = ENTRADAS[g](cel)
        saida = (estimar(entrada, modelos).clip(cel['geom'])
                 .set({'grupo': g, 'periodo': PERIODO, 'celula_id': c, 'celula_nome': cel['nome'],
                       'origem': origem, 'swir': 'estimado_rf', 'treino': linha['treino'],
                       'n_pontos_treino': linha['n_pontos_treino'], 'pontos_treino': ASSET_PONTOS,
                       'rf_params': str(params[SWIR[0]] if len({str(p) for p in params.values()}) == 1
                                         else params),
                       'bandas_estimadas': ','.join(SWIR), 'grid_asset': b.GRID_ASSET,
                       'SCALE_FACTORS': {**dict.fromkeys(VNIR + SWIR, b.ESCALA_TOA),
                                         **dict.fromkeys(TIR, 1)},
                       'DATA_TYPE': 'INT32', 'versao': VERSAO}))
        task = ee.batch.Export.image.toAsset(
            image=ee.Image(saida), description=nome, assetId=f'{DESTINO}/{nome}',
            region=cel['geom'], crs=b.CRS_EXPORT, crsTransform=b.CRS_TRANSFORM,
            maxPixels=1e13, pyramidingPolicy={'.default': 'mean'})
        try:
            task.start()
        except ee.EEException as e:
            log.error(f'{nome}: falhou ao iniciar: {e}')
            registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}'})
            continue
        log.info(f'{nome}: task {task.id} (treino {linha["treino"]})')
        registro.append({**linha, 'status': 'exportando', 'task_id': task.id})

    arq = Path(__file__).resolve().parent / f'registro_swir_rf_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['grupo', 'celula_id', 'celula_nome', 'entradas', 'treino',
                                          'n_pontos_treino', 'status', 'task_id'], extrasaction='ignore')
        w.writeheader()
        w.writerows(registro)
    log.info(f'Registro → {arq}')


if __name__ == '__main__':
    main()
