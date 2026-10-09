---
name: train-swir-models
description: Use when the user wants to train or debug the local PLS/XGBoost/GPR regression models that estimate ASTER SWIR bands (B04-B09) from VNIR+TIR, inspect metrics.json/best_params.json, or troubleshoot missing xgboost. Triggers on "train_swir_models", "treinar modelos SWIR", "PLS XGBoost GPR", "estimar bandas SWIR local".
---

# Train SWIR Models (local)

Frente local da etapa 4 do pipeline (ver
[src/building_dataset/readme_notas.md](../../../src/building_dataset/readme_notas.md)):
treina, por banda alvo, três famílias de modelo (PLS, XGBoost, GPR) para
estimar as bandas SWIR (B04-B09) do ASTER a partir de VNIR (B01, B02, B3N) +
TIR (B10-B14).

Ver [SPEC.md](SPEC.md) para o contrato formal.

## Arquivo de referência

[src/modeling/train_swir_models.py](../../../src/modeling/train_swir_models.py)
(mirror idêntico em `train_swir_models.ipynb`). Roda localmente:
`python3 src/modeling/train_swir_models.py`.

## Pré-condição crítica: `xgboost` não está instalado

Neste ambiente, `python3 -c "import xgboost"` falha com
`ModuleNotFoundError` (confirmado — `ee`, `scikit-learn`, `pandas`, `numpy`,
`joblib` estão presentes, `xgboost` não). **Antes de rodar o script**,
avise o usuário e rode `pip install xgboost` (ou peça confirmação para
instalar). Se o script for executado sem isso, ele treina PLS com sucesso,
salva `scaler.joblib` e os `pls_*.joblib`, e só quebra no bloco XGBoost —
não deixe isso acontecer silenciosamente sem avisar antes.

## Fluxo de dados e escala

- Entrada: todos os `*.csv` em `data/samples/` (não versionado no repo —
  vêm do export Drive da skill `sample-mosaico-dataset`). Se o diretório
  não existir ou estiver vazio, o script levanta `FileNotFoundError` de
  propósito — não crie dados sintéticos para contornar.
- Conversão INT16 → unidade física via dicionário `SCALE`: reflectância
  (B01,B02,B3N,B04-B09) `× 1e-4`; TIR (B10-B14) `× 0.1` (kelvin); `quality`
  também `× 1e-4`.
- Filtro `quality >= QUALITY_MIN` com `QUALITY_MIN = 0.5` — **já em unidade
  física** (pós-conversão). Isso é diferente do `QUALITY_MIN = 5000` em
  `estimate_swir_gee.py`, que filtra em INT16 bruto — mesmo limiar físico
  (~0.5), escalas diferentes. Não copie o valor literal entre os dois
  scripts.

## Notas de custo/tempo

- **GPR**: subamostrado a `GPR_NMAX = 2000` pontos (custo O(n³)); kernel
  ARD-RBF + White noise, `n_restarts_optimizer=2`.
- **XGBoost**: `RandomizedSearchCV` com `N_ITER_XGB = 60` × `N_CV = 5`
  folds × 6 bandas = até 1800 fits — pode levar minutos a dezenas de
  minutos dependendo do hardware; `n_jobs=-1` no search (não no estimador,
  para evitar conflito de paralelismo aninhado).
- **PLS**: varre `n_components` de 1 até `min(len(FEATURES), 12)` com CV
  manual (5 folds) — mais rápido, é o baseline.

## Saídas (sobrescrevem `models/` a cada rodada)

`models/scaler.joblib`, `models/pls_{banda}.joblib`,
`models/xgb_{banda}.joblib`, `models/gpr_{banda}.joblib`,
`models/best_params.json`, `models/metrics.json`. Se já existir uma rodada
anterior que o usuário queira preservar (ex.: para comparar métricas),
avise antes de rodar — o script não versiona/renomeia rodadas antigas.

## O que este skill NÃO faz

Não aplica os modelos treinados a nenhum mosaico GEE — isso é
responsabilidade da skill `estimate-swir-gee`, que treina seus próprios
classificadores GEE-side (`smileGradientTreeBoost`) e é independente
destes `.joblib` locais.
