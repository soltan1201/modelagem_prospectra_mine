"""
Estimativa de bandas SWIR (B04–B09) via GEE usando smileGradientTreeBoost.

Objetivo: estimar bandas SWIR a partir de VNIR + TIR para mosaicos onde
o detector SWIR do ASTER está inoperante (pós-2008) ou com cobertura ruim.

Fluxo:
  1. Carrega FeatureCollection do asset ASSET_POINTS (pontos amostrados)
  2. Aplica filtro de qualidade
  3. Para cada banda SWIR (B04–B09):
       a. Split 80/20 train/val com randomColumn
       b. Grid search manual sobre PARAM_GRID → escolhe params por RMSE no val
       c. Treina modelo final no conjunto completo
  4. Aplica os 6 modelos ao mosaico final
  5. Exporta imagem com bandas SWIR estimadas para asset GEE

Notas:
  - O grid search chama .getInfo() por combinação → pode levar 15–30 min
  - Os melhores params são salvos em models/gee_best_params.json
  - Para re-executar sem grid search: defina USE_SAVED_PARAMS = True
"""

import ee
import json
import logging
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
ASSET_POINTS  = 'projects/mapbiomas-arida/mine/points'
ASSET_MOSAIC  = 'projects/mapbiomas-arida/MOSAICO_FINAL_ASTER_IRECE'
ASSET_OUTPUT  = 'projects/mapbiomas-arida/mine/swir_estimated'

MODELS_DIR    = Path('models')   # onde salvar/ler gee_best_params.json
QUALITY_MIN   = 5000             # limiar INT16 (≥ 5000 ≈ quality float ≥ 0.5)
VAL_FRAC      = 0.2              # fração de validação
SEED          = 42

# True → lê best_params de models/gee_best_params.json (pula grid search)
USE_SAVED_PARAMS = False

FEATURES = ['B01', 'B02', 'B3N', 'B10', 'B11', 'B12', 'B13', 'B14']
TARGETS  = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']

# ── Grid de hiperparâmetros (smileGradientTreeBoost) ───────────────────────
# numberOfTrees ↔ n_estimators
# shrinkage     ↔ learning_rate
# samplingRate  ↔ subsample
# maxNodes      ↔ aprox. 2^max_depth  (None = sem limite)
PARAM_GRID = [
    {'numberOfTrees': 100, 'shrinkage': 0.10, 'samplingRate': 0.7, 'maxNodes': 64},
    {'numberOfTrees': 200, 'shrinkage': 0.10, 'samplingRate': 0.8, 'maxNodes': 128},
    {'numberOfTrees': 300, 'shrinkage': 0.05, 'samplingRate': 0.8, 'maxNodes': 128},
    {'numberOfTrees': 300, 'shrinkage': 0.10, 'samplingRate': 0.8, 'maxNodes': 128},
    {'numberOfTrees': 500, 'shrinkage': 0.05, 'samplingRate': 0.8, 'maxNodes': 256},
    {'numberOfTrees': 500, 'shrinkage': 0.10, 'samplingRate': 0.9, 'maxNodes': 256},
]

# ═══════════════════════════════════════════════════════════════════════════
# INICIALIZAÇÃO
# ═══════════════════════════════════════════════════════════════════════════
projAccount = 'mapbiomas-caatinga-cloud02'
try:
    ee.Initialize(project=projAccount)
    log.info('Earth Engine inicializado.')
except Exception as e:
    log.error(f'Erro de inicialização: {e}')
    raise

MODELS_DIR.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════
# FUNÇÕES
# ═══════════════════════════════════════════════════════════════════════════

def carregar_pontos() -> ee.FeatureCollection:
    """Concatena todos os assets de pontos e aplica filtro de qualidade."""
    fc_list = ee.data.listAssets({'parent': ASSET_POINTS})['assets']
    log.info(f'Assets de pontos encontrados: {len(fc_list)}')

    fcs = [ee.FeatureCollection(a['id']) for a in fc_list]
    fc  = ee.FeatureCollection(fcs).flatten()
    fc  = fc.filter(ee.Filter.gte('quality', QUALITY_MIN))

    n = fc.size().getInfo()
    log.info(f'Pontos após quality >= {QUALITY_MIN}: {n:,}')
    return fc


def split_train_val(
    fc: ee.FeatureCollection,
    val_frac: float = VAL_FRAC,
) -> tuple[ee.FeatureCollection, ee.FeatureCollection]:
    fc  = fc.randomColumn('_rand', seed=SEED)
    return (
        fc.filter(ee.Filter.gt('_rand', val_frac)),   # treino
        fc.filter(ee.Filter.lte('_rand', val_frac)),  # validação
    )


def treinar_gbt(
    fc_tr: ee.FeatureCollection,
    target: str,
    params: dict,
) -> ee.Classifier:
    return (
        ee.Classifier.smileGradientTreeBoost(
            numberOfTrees=params['numberOfTrees'],
            shrinkage=params['shrinkage'],
            samplingRate=params['samplingRate'],
            maxNodes=params.get('maxNodes'),
            seed=SEED,
        )
        .setOutputMode('REGRESSION')
        .train(
            features=fc_tr,
            classProperty=target,
            inputProperties=FEATURES,
        )
    )


def calcular_rmse_val(
    fc_val: ee.FeatureCollection,
    clf: ee.Classifier,
    target: str,
) -> float:
    """Aplica o classificador ao conjunto de validação e retorna RMSE."""
    fc_pred = fc_val.classify(clf, '_pred')
    fc_err  = fc_pred.map(
        lambda f: f.set(
            '_err_sq',
            f.getNumber('_pred').subtract(f.getNumber(target)).pow(2),
        )
    )
    mse = (
        fc_err
        .reduceColumns(ee.Reducer.mean(), ['_err_sq'])
        .get('mean')
        .getInfo()
    )
    return float(mse) ** 0.5


def grid_search(
    fc_tr: ee.FeatureCollection,
    fc_val: ee.FeatureCollection,
    target: str,
) -> tuple[dict, float]:
    """Testa PARAM_GRID e retorna os melhores params + RMSE de validação."""
    melhor_params = PARAM_GRID[0]
    melhor_rmse   = float('inf')

    for i, params in enumerate(PARAM_GRID, 1):
        try:
            clf  = treinar_gbt(fc_tr, target, params)
            rmse = calcular_rmse_val(fc_val, clf, target)
            log.info(f'    [{i}/{len(PARAM_GRID)}] {params}  →  RMSE={rmse:.4f}')
            if rmse < melhor_rmse:
                melhor_rmse, melhor_params = rmse, params
        except Exception as exc:
            log.warning(f'    [{i}/{len(PARAM_GRID)}] Falhou: {exc}')
        time.sleep(0.3)

    return melhor_params, melhor_rmse


# ═══════════════════════════════════════════════════════════════════════════
# PIPELINE PRINCIPAL
# ═══════════════════════════════════════════════════════════════════════════
log.info('Carregando pontos...')
fc_all        = carregar_pontos()
fc_tr, fc_val = split_train_val(fc_all)

n_tr  = fc_tr.size().getInfo()
n_val = fc_val.size().getInfo()
log.info(f'Treino: {n_tr:,}  |  Validação: {n_val:,}')

# Carregar params salvos ou executar grid search
gee_best_params_path = MODELS_DIR / 'gee_best_params.json'

if USE_SAVED_PARAMS and gee_best_params_path.exists():
    with open(gee_best_params_path, encoding='utf-8') as f:
        saved = json.load(f)
    melhores_params = {b: saved[b]['params'] for b in TARGETS}
    log.info('Params carregados de gee_best_params.json (grid search ignorado).')
else:
    melhores_params = {}
    gee_metrics     = {}

    for banda in TARGETS:
        log.info(f'\n─── Grid search: {banda} ───')
        params, rmse_val = grid_search(fc_tr, fc_val, banda)
        melhores_params[banda] = params
        gee_metrics[banda]     = {'val_rmse': rmse_val, 'params': params}
        log.info(f'  Melhor: {params}  RMSE_val={rmse_val:.4f}')

    # Persistir resultados do grid search
    with open(gee_best_params_path, 'w', encoding='utf-8') as f:
        json.dump(gee_metrics, f, indent=2)
    log.info(f'\nParams salvos → {gee_best_params_path}')

# Treinar modelos finais no conjunto completo
log.info('\nTreinando modelos finais (conjunto completo)...')
classifiers: dict[str, ee.Classifier] = {}
for banda in TARGETS:
    log.info(f'  {banda}')
    classifiers[banda] = treinar_gbt(fc_all, banda, melhores_params[banda])

# ═══════════════════════════════════════════════════════════════════════════
# APLICAR AO MOSAICO E EXPORTAR
# ═══════════════════════════════════════════════════════════════════════════
log.info('\nAplicando modelos ao mosaico...')
mosaico      = ee.Image(ASSET_MOSAIC)
img_features = mosaico.select(FEATURES)
area         = mosaico.geometry()

bandas_estimadas = [
    img_features.classify(classifiers[banda], f'{banda}_est')
    for banda in TARGETS
]

img_estimada = (
    ee.Image.cat(bandas_estimadas)
    .int16()
    .set('description', 'SWIR estimado via GBT (VNIR+TIR)',
         'modelos', str(melhores_params))
)

log.info(f'Exportando → {ASSET_OUTPUT}')
task = ee.batch.Export.image.toAsset(
    image=img_estimada,
    description='swir_estimated_gbt',
    assetId=ASSET_OUTPUT,
    region=area,
    scale=30,
    maxPixels=1e13,
    pyramidingPolicy={'.default': 'mean'},
)
task.start()
log.info('Task iniciada.')
log.info('Acompanhe em: https://code.earthengine.google.com/tasks')
