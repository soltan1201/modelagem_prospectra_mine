# Spec — estimate-swir-gee

## Objetivo
Rodar/editar `src/modeling/estimate_swir_gee.py`, que treina
`smileGradientTreeBoost` no GEE para estimar SWIR (B04-B09), faz grid
search de hiperparâmetros, aplica ao mosaico final e exporta o resultado
como asset.

## Pré-condições
- `ee.Initialize(project='mapbiomas-caatinga-cloud02')` deve suceder
  (mesma conta/projeto usado em `sample_quality_mosaic.py`).
- Asset(s) de pontos existem sob
  `projects/mapbiomas-arida/mine/points` (produzidos pela skill
  `sample-mosaico-dataset`).
- Asset `projects/mapbiomas-arida/MOSAICO_FINAL_ASTER_IRECE` existe
  (produzido pela skill `merge-mosaico-final`) e contém as bandas
  `FEATURES` (`B01,B02,B3N,B10-B14`).
- Diretório `models/` gravável (para `gee_best_params.json`).

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| `USE_SAVED_PARAMS` (True/False) | usuário — pular ou não o grid search | sim, se já houve rodada anterior |
| `QUALITY_MIN` (default 5000, INT16) | usuário (raro mudar) | não |
| `PARAM_GRID` | usuário (raro mudar sem motivo) | não |

## Saídas
- `models/gee_best_params.json` (se grid search rodou): RMSE de validação
  + params por banda.
- Task de export iniciada (`task.start()`) → asset
  `projects/mapbiomas-arida/mine/swir_estimated` (INT16, 6 bandas
  `{banda}_est`).
- Log no console com progresso do grid search por banda/combinação.

## Efeitos colaterais
- **Grid search é lento (15-30 min)** e faz muitas chamadas `.getInfo()`
  síncronas — avisar o usuário antes de iniciar sem `USE_SAVED_PARAMS`.
- **Export final consome cota de tasks do GEE** — é uma ação real e
  visível no projeto `mapbiomas-arida`, não reversível localmente (embora
  o asset possa ser deletado depois no Code Editor/`earthengine`).

## Modos de falha e como reagir
- **`ee.data.listAssets` retorna lista vazia**: os pontos ainda não foram
  exportados como asset — direcionar para a skill `sample-mosaico-dataset`
  com `EXPORTAR_PARA` incluindo `'asset'`, não inventar dados.
- **RMSE de validação muito alto ou grid search falhando em todas as
  combinações**: reportar os erros capturados (o loop já faz
  `try/except` por combinação e loga `Falhou: {exc}`) em vez de assumir
  sucesso silencioso.
- **Confusão de escala do `QUALITY_MIN`**: se o usuário pedir para "usar o
  mesmo QUALITY_MIN do outro script", esclarecer que 5000 (aqui) e 0.5
  (`train_swir_models.py`) já são equivalentes — não os iguale
  literalmente.

## Não-objetivos
- Não treina nem lê os modelos `.joblib` locais — é um pipeline de
  modelagem paralelo e independente ao da skill `train-swir-models`.
- Não decide sozinho quando reexecutar o grid search vs. usar
  `USE_SAVED_PARAMS` — perguntar a intenção do usuário.
