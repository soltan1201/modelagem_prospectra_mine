# Spec — train-swir-models

## Objetivo
Rodar/editar `src/modeling/train_swir_models.py`, que treina PLS, XGBoost
e GPR por banda alvo (B04-B09) a partir de CSVs de amostras locais e
persiste modelos + métricas em `models/`.

## Pré-condições
- Diretório `data/samples/` existe e contém ao menos um `*.csv` (exportado
  pela skill `sample-mosaico-dataset` com `EXPORTAR_PARA` incluindo
  `'drive'`, depois baixado manualmente para o repo).
- Pacotes instalados: `numpy`, `pandas`, `joblib`, `scikit-learn`, `scipy`
  (confirmados disponíveis) e **`xgboost`** (confirmado ausente — instalar
  antes de rodar).
- CSVs contêm as colunas `FEATURES + TARGETS` (`B01,B02,B3N,B10-B14,
  B04-B09`) e, se filtro de qualidade for aplicado, a coluna `quality`.

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| CSVs em `data/samples/*.csv` | usuário (download manual do Drive) | sim |
| `QUALITY_MIN` (default 0.5, física) | usuário (raro mudar) | não |
| `GPR_NMAX`, `N_ITER_XGB`, `N_CV` | usuário (custo vs. precisão) | não |

## Saídas
- `models/scaler.joblib`
- `models/pls_{B04..B09}.joblib`, `models/xgb_{B04..B09}.joblib`,
  `models/gpr_{B04..B09}.joblib` (18 arquivos de modelo)
- `models/best_params.json` (hiperparâmetros PLS/XGB por banda)
- `models/metrics.json` (R², RMSE, MAE por modelo × banda)
- Log no console com tabela-resumo final de R² por banda/modelo

## Efeitos colaterais
- Sobrescreve qualquer artefato existente em `models/` com o mesmo nome
  sem backup automático.
- Uso intenso de CPU durante `RandomizedSearchCV` (XGBoost) — pode saturar
  todos os núcleos (`N_JOBS = -1`).

## Modos de falha e como reagir
- **`ModuleNotFoundError: xgboost`**: instalar (`pip install xgboost`)
  antes de rodar — não tentar contornar removendo o bloco XGBoost do
  script sem que o usuário peça isso explicitamente.
- **`FileNotFoundError` em `data/samples`**: confirmar com o usuário que os
  CSVs já foram baixados do Drive (etapa anterior do pipeline); não gerar
  dados sintéticos para "fazer o script rodar".
- **Poucos registros após filtro de qualidade**: se `len(df)` cair para
  perto de zero após `quality >= QUALITY_MIN`, avisar — pode indicar
  `QUALITY_MIN` alto demais para os dados amostrados ou mosaico de baixa
  qualidade geral.

## Não-objetivos
- Não decide o `QUALITY_MIN` ideal sozinho — é uma escolha de trade-off
  entre volume de dados e pureza dos exemplos, do usuário.
- Não aplica os modelos treinados a nenhuma imagem/mosaico (isso é feito
  fora deste script, tipicamente em notebook separado usando os
  `.joblib`).
