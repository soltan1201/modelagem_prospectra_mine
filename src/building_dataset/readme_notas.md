# Notas — Pipeline de construção do dataset ASTER

## Objetivo geral

Construir um dataset de amostras espectrais multitemporais do sensor **ASTER** (Terra/NASA) sobre a região de **Irecê (BA)**, usando mosaicos semestrais de qualidade para extrair pontos representativos de solo exposto / litologia / vegetação em contexto semiárido.

O pipeline está dividido em três etapas sequenciais:

---

## Etapas do pipeline

### 1. `select_black_list_save_mosaicSem.js` — Mosaico semestral por qualidade

**Entrada:** Coleção `ASTER/AST_L1T_003` (radiância L1T)  
**Saída:** Asset INT16 em `projects/mapbiomas-arida/mosaic_aster/ASTER_QualityMosaic_<ano>_semestre_<N>_INT16`

Processa as imagens ASTER de um semestre específico:

- Converte **radiância → reflectância TOA** (bandas VNIR/SWIR)
- Remove imagens problemáticas via **blacklist** (IDs explícitos)
- Aplica **máscara de nuvens e sombras** com critérios espectrais e térmicos
- Calcula **score de qualidade pixel a pixel** (`addQualityASTER_v5_shadow_fix`) com penalização de sombras, nuvens e brilho excessivo
- Compõe o mosaico com **`qualityMosaic('quality')`** — cada pixel vem da imagem com maior score naquele local
- Exporta em **INT16** (fator escala = 10 000 para reflectância; 10 000 para quality)

> Parâmetros ajustáveis no topo: `ano`, `mes`, `blacklist`, limiares de nuvem/sombra.

---

### 2. `merge_mosaic_semestral.js` — Mosaico anual/multitemporal

**Entrada:** `ImageCollection` em `projects/mapbiomas-arida/mosaic_aster`  
**Saída:** Asset `projects/mapbiomas-arida/MOSAICO_FINAL_ASTER_IRECE`

Combina todos os mosaicos semestrais disponíveis em um único mosaico:

- Carrega a coleção de assets INT16 (cada um já tem a banda `quality`)
- Aplica **`qualityMosaic('quality')`** diretamente — seleciona o melhor pixel entre todos os semestres
- Visualiza: RGB Natural, SWIR Falsa-cor, banda quality
- Exporta o mosaico final para asset

> Não recalcula quality — usa a banda já salva em cada semestre.

---

### 3. `sample_quality_mosaic.py` — Amostragem para treinamento/análise

**Entrada:** `ImageCollection` em `projects/mapbiomas-arida/mosaic_aster`  
**Saída:** CSVs no Google Drive em `ASTER_Samples/<system:index>.csv`

Para cada mosaico semestral:

1. **Harmoniza todas as bandas para 30 m** antes de amostrar:
   - VNIR (15 m → 30 m): `reduceResolution(mean)` — média de 4 pixels (2×2)
   - SWIR (30 m): sem alteração
   - TIR (90 m → 30 m): `resample('bilinear')` — interpolação bilinear
   - `quality` (15 m → 30 m): `reduceResolution(mean)` — preservada no CSV
2. Coleta **30 000 pontos aleatórios** sobre todos os pixels válidos (sem filtro prévio de qualidade)
3. Exporta como **CSV** nomeado com o `system:index` da imagem

> O CSV inclui todas as bandas espectrais + `quality`. O filtro por qualidade é feito em pós-processamento (pandas/R) sobre a coluna `quality`.

**Colunas do CSV:** `B01, B02, B3N, B04, B05, B06, B07, B08, B09, B10, B11, B12, B13, B14, quality`

**Parâmetros principais:**

| Parâmetro | Valor | Descrição |
|---|---|---|
| `N_PONTOS` | 30 000 | Pontos amostrados por imagem |
| `ESCALA_AMOST` | 30 m | Resolução de amostragem |
| `TILE_SCALE` | 4 | Particionamento interno GEE |
| `SEED` | 42 | Semente aleatória (reprodutibilidade) |
| `INCLUIR_COORDS` | False | `True` adiciona lat/lon no CSV |

---

## Área de estudo

**Bloco:** Irecê (BA)  
**Coordenadas (bounding box):**

```
[-44.460483, -14.790793]  (SO)
[-39.203525, -14.790793]  (SE)
[-39.203525, -9.416525]   (NE)
[-44.460483, -9.416525]   (NO)
```

---

## Bandas ASTER

| Nome | Descrição                         | Min | Max  | Resolução | Comprimento de onda |
|------|-----------------------------------|-----|------|-----------|---------------------|
| B01  | VNIR_Band1 (visible green/yellow) | 1   | 255  | 15 m      | 0.520–0.600 µm      |
| B02  | VNIR_Band2 (visible red)          | 1   | 255  | 15 m      | 0.630–0.690 µm      |
| B3N  | VNIR_Band3N (near infrared)       | 1   | 255  | 15 m      | 0.780–0.860 µm      |
| B04  | SWIR_Band4 (short-wave infrared)  | 1   | 255  | 30 m      | 1.600–1.700 µm      |
| B05  | SWIR_Band5 (short-wave infrared)  | 1   | 255  | 30 m      | 2.145–2.185 µm      |
| B06  | SWIR_Band6 (short-wave infrared)  | 1   | 255  | 30 m      | 2.185–2.225 µm      |
| B07  | SWIR_Band7 (short-wave infrared)  | 1   | 255  | 30 m      | 2.235–2.285 µm      |
| B08  | SWIR_Band8 (short-wave infrared)  | 1   | 255  | 30 m      | 2.295–2.365 µm      |
| B09  | SWIR_Band9 (short-wave infrared)  | 1   | 255  | 30 m      | 2.360–2.430 µm      |
| B10  | TIR_Band10 (thermal infrared)     | 1   | 4095 | 90 m      | 8.125–8.475 µm      |
| B11  | TIR_Band11 (thermal infrared)     | 1   | 4095 | 90 m      | 8.475–8.825 µm      |
| B12  | TIR_Band12 (thermal infrared)     | 1   | 4095 | 90 m      | 8.925–9.275 µm      |
| B13  | TIR_Band13 (thermal infrared)     | 1   | 4095 | 90 m      | 10.250–10.950 µm    |
| B14  | TIR_Band14 (thermal infrared)     | 1   | 4095 | 90 m      | 10.950–11.650 µm    |

> Assets exportados em **INT16** com fator de escala **10 000** (reflectância) e **10** (TIR temperatura).  
> Banda `quality` exportada com fator **10 000** (range float 0–2 → INT16 0–20 000).

---

## Assets GEE

| Descrição | Caminho |
|---|---|
| Mosaicos semestrais | `projects/mapbiomas-arida/mosaic_aster/` |
| Mosaico final multitemporal | `projects/mapbiomas-arida/MOSAICO_FINAL_ASTER_IRECE` |

---

## Modelagem — Estimativa de bandas SWIR

Tarefa de **regressão espectral**: estimar B04–B09 (SWIR) a partir de B01, B02, B3N (VNIR) + B10–B14 (TIR).  
Scripts em `src/modeling/`.

### Grupos de modelos (do mais simples ao mais complexo)

#### 1. Regressão Linear e Variantes

| Método | Característica |
|---|---|
| OLS (mínimos quadrados) | Baseline. Rápido, interpretável, assume relação linear entre bandas |
| Ridge / Lasso | Regulariza quando há multicolinearidade entre B01–B3N |
| Elastic Net | Combinação Ridge+Lasso, bom ponto de partida regularizado |

#### 2. Partial Least Squares (PLS)

O mais usado em **espectroscopia e sensoriamento remoto** para predição de bandas. Reduz dimensionalidade e faz regressão em um único passo — lida bem com preditores correlacionados. `n_components` otimizado por cross-validation.

#### 3. Regressão por Árvores (Ensemble)

| Método | Característica |
|---|---|
| Random Forest | Robusto a outliers, captura não-linearidades, indica importância das bandas |
| Gradient Boosting (XGBoost / LightGBM) | Melhor acurácia, mais sensível a hiperparâmetros |

Com 30 000 pontos × 15 imagens, esses modelos têm dados suficientes para generalizar bem.

#### 4. Redes Neurais

| Arquitetura | Uso |
|---|---|
| MLP (Multi-Layer Perceptron) | Simples, captura relações não-lineares complexas entre bandas |
| Autoencoder | Aprende representação latente do espectro e reconstrói bandas faltantes |

#### 5. Regressão por Processos Gaussianos (GPR)

Fornece **intervalo de confiança** na predição — útil para identificar onde a estimativa é menos confiável. Subamostrado por custo O(n³).

### Abordagem implementada

```
OLS baseline → PLS → XGBoost (RandomizedSearchCV + KFold) → GPR (subamostrado)
```

| Script | Ambiente | Descrição |
|---|---|---|
| `src/modeling/train_swir_models.py` | Local (notebook) | Treina PLS + XGBoost + GPR; salva modelos e métricas em `models/` |
| `src/modeling/estimate_swir_gee.py` | GEE Python | Grid search de smileGradientTreeBoost nos pontos do asset; aplica ao mosaico final |

**Saídas de `train_swir_models.py`:**

| Arquivo | Conteúdo |
|---|---|
| `models/scaler.joblib` | StandardScaler ajustado às features |
| `models/pls_{banda}.joblib` | Modelo PLS por banda SWIR |
| `models/xgb_{banda}.joblib` | Modelo XGBoost por banda SWIR |
| `models/gpr_{banda}.joblib` | Modelo GPR por banda SWIR |
| `models/best_params.json` | Melhores hiperparâmetros (PLS: n\_components; XGB: todos os params) |
| `models/metrics.json` | R², RMSE, MAE por modelo e banda |
| `models/gee_best_params.json` | Melhores params do smileGradientTreeBoost (gerado pelo script GEE) |

## Notas

<!-- Adicione observações, decisões e resultados abaixo -->

- [ ]  
