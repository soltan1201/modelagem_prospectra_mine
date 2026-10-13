---
name: estimate-swir-gee
description: Use when the user wants to estimate ASTER SWIR (B04-B09) in GEE — the current per-cell smileRandomForest scripts (estimar_swir_rf/_v2, with per-scene OLS normalization) or the legacy smileGradientTreeBoost path (grid search, apply to the final mosaic, export swir_estimated). Triggers on "estimate_swir_gee", "smileGradientTreeBoost", "grid search GEE", "swir_estimated", "aplicar modelo ao mosaico final", "estimar_swir_rf", "smileRandomForest", "normalização OLS por cena", "ASTER_MOSAIC_v2".
---

# Estimate SWIR (GEE-side)

Frente GEE-side da etapa 4 do pipeline (ver
[src/building_dataset/readme_notas.md](../../../src/building_dataset/readme_notas.md)):
treina classificadores `ee.Classifier.smileGradientTreeBoost` (regressão)
por banda SWIR alvo, diretamente com os pontos amostrados já como asset,
depois aplica os 6 modelos ao mosaico final e exporta.

Ver [SPEC.md](SPEC.md) para o contrato formal.

## Caminho atual: Random Forest por célula (`estimar_swir_rf*.py`)

O GBT abaixo é legado. O SWIR estimado em uso vem de
`src/modeling/estimar_swir_rf.py` (v1) e `estimar_swir_rf_v2.py`
(VERSAO `swir_rf_v3_seco`):
- treino: `projects/mapbiomas-arida/mine/points_zeroCC_00_08` (pontos do
  zero-cloud 2000–2008); um `smileRandomForest` (REGRESSION) por banda
  B04–B09, features B01, B02, B3N, B10–B14; `TREINO='global'` ou
  `'celula'` (com vizinhas); limite `N_MAX_TREINO=6000`
  ("Computed value is too large");
- hiperparâmetros: `models/zero_cloud/best_params.json` ou `RF_PARAMS`;
- **normalização OLS por cena** (2009–2020 → 2000–2008) com a banda `data`;
  coeficientes em `models/zero_cloud/normalizacao_ols_cena_coef.csv`;
  cena com < `MIN_PONTOS_CENA` usa reta global. Sem normalização o teste
  temporal deu R² negativo (TIR em DN carrega temperatura da data);
- v2 entradas: `ASTER_0cloud_v2` (grupo `zero_cloud`, norm. `cena`) e
  `ASTER_v5` (grupo `semestral_v5`, norm. `nenhuma`); saída
  `projects/mapbiomas-arida/mine/ASTER_MOSAIC_v2/ASTER_MOSAIC_<grupo>_2009_2020_<LxCy>`
  com banda `normalizacao` (0/1/2);
- `rf.estimar(...)` é reutilizado por `build_mosaico_unico.py` (etapa 2).

## Arquivo de referência

[src/modeling/estimate_swir_gee.py](../../../src/modeling/estimate_swir_gee.py)
(mirror em `.ipynb`). Roda localmente via `python3`, mas todo o cômputo
pesado acontece no lado do servidor GEE.

## Pré-condição

Pontos já exportados como asset em
`projects/mapbiomas-arida/mine/points` (produzidos pela skill
`sample-mosaico-dataset` com `EXPORTAR_PARA` incluindo `'asset'`) —
`carregar_pontos()` lista e concatena **todos** os assets sob esse
diretório via `ee.data.listAssets`.

## Exportar pontos para o Drive (só no notebook)

O `.ipynb` (não o `.py`) tem uma célula logo após `carregar_pontos()` que
exporta `fc_all` (o conjunto mesclado e já filtrado por `quality`, o mesmo
usado no grid search e no treino final) como um único CSV para
`Drive/ASTER_SWIR_GEE_Points/`, controlada pela flag
`EXPORTAR_PONTOS_DRIVE` (default `True`). É útil para inspecionar/arquivar
exatamente o conjunto de pontos que alimentou o modelo GEE-side desta
rodada — diferente das exportações por semestre da skill
`sample-mosaico-dataset`, que geram um CSV por mosaico, antes da fusão e
do filtro de qualidade final. Se o usuário pedir a mesma exportação no
`.py`, replique o mesmo padrão de `Export.table.toDrive` logo após
`carregar_pontos()`.

## Fluxo e pontos de atenção

1. `carregar_pontos()` filtra `quality >= QUALITY_MIN` com
   `QUALITY_MIN = 5000` — **em unidade INT16 bruta** (≈ quality física
   0.5). Não confundir com o `QUALITY_MIN = 0.5` (física) de
   `train_swir_models.py` — são o mesmo limiar físico, escalas diferentes.
2. Split 80/20 treino/validação via `randomColumn('_rand', seed=SEED)`.
3. **Grid search manual** sobre `PARAM_GRID` (6 combinações de
   `numberOfTrees`/`shrinkage`/`samplingRate`/`maxNodes`), uma por banda
   (6 bandas × 6 combinações = 36 treinos + 36 `.getInfo()` de RMSE) —
   **leva de 15 a 30 minutos**. Avise o usuário antes de rodar do zero.
4. **`USE_SAVED_PARAMS = True`** pula o grid search e lê
   `models/gee_best_params.json` (gerado pela rodada anterior) — ofereça
   essa opção para reexecuções.
5. Treina o modelo final de cada banda no conjunto **completo** (não só
   treino), depois aplica aos 6 alvos sobre `ASSET_MOSAIC`
   (`projects/mapbiomas-arida/MOSAICO_FINAL_ASTER_IRECE`).
6. Exporta como asset INT16 em `ASSET_OUTPUT`
   (`projects/mapbiomas-arida/mine/swir_estimated`) via
   `Export.image.toAsset` com `pyramidingPolicy={'.default': 'mean'}` —
   dispara uma task real no GEE (`task.start()`), consome cota.

## O que este skill NÃO faz

Não usa os `.joblib` treinados localmente (skill `train-swir-models`) —
os dois caminhos de modelagem (scikit-learn local vs. smileGradientTreeBoost
GEE-side) são independentes e comparáveis via `metrics.json` vs.
`gee_best_params.json`, mas não compartilham artefatos.
