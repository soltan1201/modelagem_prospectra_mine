"""
Treinamento local — estimativa de bandas SWIR (B04–B09) do ASTER
a partir de VNIR (B01, B02, B3N) + TIR (B10–B14).

Modelos treinados por banda alvo (6 bandas × 3 modelos = 18 modelos):
  PLS  — Partial Least Squares  (n_components otimizado por CV)
  XGB  — XGBoost                (RandomizedSearchCV + KFold)
  GPR  — Gaussian Process       (subamostrado, kernel ARD-RBF)

Entrada : data/samples/*.csv   (CSVs baixados do Google Drive)
Saídas  : models/
            scaler.joblib
            pls_{banda}.joblib / xgb_{banda}.joblib / gpr_{banda}.joblib
            best_params.json   → melhores hiperparâmetros por modelo e banda
            metrics.json       → R², RMSE, MAE por modelo e banda
"""

import os
import json
import glob
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import joblib

from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (
    RBF, WhiteKernel, ConstantKernel as C,
)
from sklearn.model_selection import (
    KFold, RandomizedSearchCV, train_test_split,
)
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from scipy.stats import uniform, randint

import xgboost as xgb

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CSV_DIR    = 'data/samples'   # pasta com CSVs baixados do Google Drive
MODELS_DIR = Path('models')   # saída: modelos + artefatos

N_JOBS      = -1    # -1 = todos os núcleos da CPU
N_CV        = 5     # folds para KFold
N_ITER_XGB  = 60    # iterações RandomizedSearchCV XGBoost
GPR_NMAX    = 2000  # máx. amostras para GPR (custo O(n³))
QUALITY_MIN = 0.5   # limiar float pós-conversão (quality INT16 → float ÷ 10 000)
SEED        = 42

FEATURES = ['B01', 'B02', 'B3N', 'B10', 'B11', 'B12', 'B13', 'B14']
TARGETS  = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']

# Conversão INT16 → unidade física
# VNIR/SWIR: reflectância × 10 000  |  TIR: temperatura (K) × 10
SCALE: dict = {
    **dict.fromkeys(
        ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09'],
        1e-4,
    ),
    **dict.fromkeys(['B10', 'B11', 'B12', 'B13', 'B14'], 0.1),
    'quality': 1e-4,
}

# ── Espaço de hiperparâmetros XGBoost ──────────────────────────────────────
XGB_PARAM_DIST = {
    'n_estimators':     randint(100, 700),
    'max_depth':        randint(3, 10),
    'learning_rate':    uniform(0.01, 0.29),
    'subsample':        uniform(0.5, 0.5),
    'colsample_bytree': uniform(0.5, 0.5),
    'min_child_weight': randint(1, 10),
    'gamma':            uniform(0.0, 0.5),
    'reg_alpha':        uniform(0.0, 1.0),
    'reg_lambda':       uniform(0.5, 1.5),
}

# ═══════════════════════════════════════════════════════════════════════════
# UTILITÁRIOS
# ═══════════════════════════════════════════════════════════════════════════

MODELS_DIR.mkdir(parents=True, exist_ok=True)


def _rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.sqrt(mean_squared_error(y_true, y_pred)))


def avaliar(y_true, y_pred, modelo: str, banda: str) -> dict:
    m = {
        'r2':   float(r2_score(y_true, y_pred)),
        'rmse': _rmse(y_true, y_pred),
        'mae':  float(mean_absolute_error(y_true, y_pred)),
    }
    log.info(
        f'    [{modelo}] {banda}  '
        f'R²={m["r2"]:.4f}  RMSE={m["rmse"]:.6f}  MAE={m["mae"]:.6f}'
    )
    return m


# ═══════════════════════════════════════════════════════════════════════════
# CARGA E PREPARAÇÃO DOS DADOS
# ═══════════════════════════════════════════════════════════════════════════

def carregar_dados() -> tuple[np.ndarray, np.ndarray]:
    arquivos = sorted(glob.glob(os.path.join(CSV_DIR, '*.csv')))
    if not arquivos:
        raise FileNotFoundError(f'Nenhum CSV encontrado em {CSV_DIR!r}')
    log.info(f'Carregando {len(arquivos)} CSV(s)...')

    df = pd.concat([pd.read_csv(f) for f in arquivos], ignore_index=True)
    log.info(f'Registros brutos: {len(df):,}')

    # Converter INT16 → unidade física
    for col, fator in SCALE.items():
        if col in df.columns:
            df[col] = df[col] * fator

    # Filtro de qualidade
    if 'quality' in df.columns:
        n_antes = len(df)
        df = df[df['quality'] >= QUALITY_MIN].copy()
        log.info(f'quality >= {QUALITY_MIN}: {n_antes:,} → {len(df):,} registros')

    df = df[FEATURES + TARGETS].dropna()
    log.info(f'Registros após limpeza: {len(df):,}')

    return (
        df[FEATURES].values.astype(np.float32),
        df[TARGETS].values.astype(np.float32),
    )


X, Y = carregar_dados()

X_tr, X_te, Y_tr, Y_te = train_test_split(
    X, Y, test_size=0.2, random_state=SEED
)
log.info(f'Treino: {X_tr.shape[0]:,}  |  Teste: {X_te.shape[0]:,}')

scaler   = StandardScaler()
X_tr_s   = scaler.fit_transform(X_tr)
X_te_s   = scaler.transform(X_te)
joblib.dump(scaler, MODELS_DIR / 'scaler.joblib')
log.info('Scaler salvo → models/scaler.joblib')

kf = KFold(n_splits=N_CV, shuffle=True, random_state=SEED)

all_metrics: dict = {}
all_params:  dict = {}


# ═══════════════════════════════════════════════════════════════════════════
# 1 — PLS (Partial Least Squares)
# ═══════════════════════════════════════════════════════════════════════════
log.info('\n──────────────────── PLS ────────────────────')
pls_metrics, pls_params = {}, {}
MAX_COMP = min(len(FEATURES), 12)

for j, banda in enumerate(TARGETS):
    log.info(f'  Banda {banda}')
    y_tr, y_te = Y_tr[:, j], Y_te[:, j]

    # Encontra n_components ótimo por CV-RMSE
    best_nc, best_cv_rmse = 1, float('inf')
    for nc in range(1, MAX_COMP + 1):
        fold_rmses = []
        for tr_idx, val_idx in kf.split(X_tr_s):
            m = PLSRegression(n_components=nc)
            m.fit(X_tr_s[tr_idx], y_tr[tr_idx])
            fold_rmses.append(_rmse(y_tr[val_idx], m.predict(X_tr_s[val_idx]).ravel()))
        media = float(np.mean(fold_rmses))
        if media < best_cv_rmse:
            best_cv_rmse, best_nc = media, nc

    log.info(f'    n_components={best_nc}  CV-RMSE={best_cv_rmse:.6f}')
    modelo = PLSRegression(n_components=best_nc)
    modelo.fit(X_tr_s, y_tr)

    pls_metrics[banda] = avaliar(y_te, modelo.predict(X_te_s).ravel(), 'PLS', banda)
    pls_params[banda]  = {'n_components': best_nc, 'cv_rmse': best_cv_rmse}
    joblib.dump(modelo, MODELS_DIR / f'pls_{banda}.joblib')

all_metrics['pls'] = pls_metrics
all_params['pls']  = pls_params


# ═══════════════════════════════════════════════════════════════════════════
# 2 — XGBoost (RandomizedSearchCV)
# ═══════════════════════════════════════════════════════════════════════════
log.info('\n──────────────────── XGBoost ────────────────────')
xgb_metrics, xgb_params = {}, {}

for j, banda in enumerate(TARGETS):
    log.info(f'  Banda {banda}')
    y_tr, y_te = Y_tr[:, j], Y_te[:, j]

    # n_jobs=1 no estimador evita conflito com n_jobs=-1 no SearchCV
    base = xgb.XGBRegressor(
        tree_method='hist',
        random_state=SEED,
        n_jobs=1,
        verbosity=0,
    )
    search = RandomizedSearchCV(
        estimator=base,
        param_distributions=XGB_PARAM_DIST,
        n_iter=N_ITER_XGB,
        cv=kf,
        scoring='neg_root_mean_squared_error',
        random_state=SEED,
        n_jobs=N_JOBS,
        refit=True,
        verbose=0,
    )
    search.fit(X_tr_s, y_tr)

    melhor = search.best_estimator_
    xgb_metrics[banda] = avaliar(y_te, melhor.predict(X_te_s), 'XGB', banda)
    xgb_params[banda]  = {
        **search.best_params_,
        'cv_rmse': float(-search.best_score_),
    }
    log.info(f'    Melhores params: {search.best_params_}')
    joblib.dump(melhor, MODELS_DIR / f'xgb_{banda}.joblib')

all_metrics['xgb'] = xgb_metrics
all_params['xgb']  = xgb_params


# ═══════════════════════════════════════════════════════════════════════════
# 3 — GPR (Gaussian Process Regression)
# ═══════════════════════════════════════════════════════════════════════════
log.info('\n──────────────────── GPR ────────────────────')
gpr_metrics = {}

# Subamostrar para viabilidade (GPR = O(n³); 2 000 pontos ≈ aceitável)
np.random.seed(SEED)
idx_gpr = np.random.choice(len(X_tr_s), min(GPR_NMAX, len(X_tr_s)), replace=False)
X_gpr_s = X_tr_s[idx_gpr]
log.info(f'  Subamostrado: {len(X_gpr_s):,} pontos de treino')

# Kernel ARD-RBF: comprimento de escala independente por feature + ruído branco
kernel = (
    C(1.0, (1e-3, 1e3))
    * RBF(
        length_scale=np.ones(len(FEATURES)),
        length_scale_bounds=(1e-2, 1e2),
    )
    + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-5, 1.0))
)

for j, banda in enumerate(TARGETS):
    log.info(f'  Banda {banda}')
    y_gpr = Y_tr[idx_gpr, j]
    modelo = GaussianProcessRegressor(
        kernel=kernel,
        alpha=0.0,
        normalize_y=True,
        n_restarts_optimizer=2,
        random_state=SEED,
    )
    modelo.fit(X_gpr_s, y_gpr)
    gpr_metrics[banda] = avaliar(Y_te[:, j], modelo.predict(X_te_s), 'GPR', banda)
    joblib.dump(modelo, MODELS_DIR / f'gpr_{banda}.joblib')

all_metrics['gpr'] = gpr_metrics


# ═══════════════════════════════════════════════════════════════════════════
# SALVAR ARTEFATOS
# ═══════════════════════════════════════════════════════════════════════════
with open(MODELS_DIR / 'best_params.json', 'w', encoding='utf-8') as f:
    json.dump(all_params, f, indent=2, default=str)

with open(MODELS_DIR / 'metrics.json', 'w', encoding='utf-8') as f:
    json.dump(all_metrics, f, indent=2)

# ── Tabela resumo final ────────────────────────────────────────────────────
log.info('\n══════════════════ RESUMO R² ══════════════════')
log.info(f"{'Banda':>6} | {'PLS':>7} | {'XGB':>7} | {'GPR':>7}")
log.info('-' * 38)
for banda in TARGETS:
    r2s = {m: all_metrics[m][banda]['r2'] for m in ('pls', 'xgb', 'gpr')}
    log.info(
        f'{banda:>6} | {r2s["pls"]:7.4f} | {r2s["xgb"]:7.4f} | {r2s["gpr"]:7.4f}'
    )

log.info(f'\nModelos    → {MODELS_DIR}/')
log.info(f'Parâmetros → {MODELS_DIR / "best_params.json"}')
log.info(f'Métricas   → {MODELS_DIR / "metrics.json"}')
