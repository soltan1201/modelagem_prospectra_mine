# Notas — Pipeline de construção do dataset ASTER

## Objetivo geral

Construir um dataset de amostras espectrais multitemporais do sensor **ASTER** (Terra/NASA) sobre a região de **Irecê (BA)**, usando mosaicos semestrais de qualidade para extrair pontos representativos de solo exposto / litologia / vegetação em contexto semiárido.

O pipeline está dividido em três etapas sequenciais (fluxo original, em JS; o fluxo
atual é o sazonal seco, descrito na seção seguinte):

---

## Fluxo atual — mosaico sazonal seco (2026-10)

**Por quê:** os mosaicos v1 (prioridade por origem) e único (prioridade à seca)
ainda mostram emendas entre cenas ("colcha de retalhos"). O `qualityMosaic` pega
uma cena por pixel e os mosaicos eram montados semestre a semestre, com poucas
cenas candidatas por célula. No fluxo sazonal, todas as cenas secas do período
inteiro entram juntas antes de escolher o pixel.

**Teste que definiu o fluxo:** `mosaico_sazonal_homogeneo.js` (Code Editor) compara,
por célula, seco × chuvoso e três composições — `quality` (qualityMosaic),
`mediana` e `mediana_norm` (cada cena ajustada à mediana da estação) — com um
índice de costura (gradiente nas bordas de cena ÷ gradiente fora). Melhor
resultado visual: **A | seco | quality**. A mediana mistura datas no mesmo pixel e
perde contraste. O JS tem as mesmas listas de cenas do Python (`BLACKLIST_SECO`,
`ZERO_CLOUD_SECO`), geradas dos `.txt`; regenerar se os `.txt` mudarem.

**Períodos** (separados porque só A tem SWIR medido):

| Período | Datas | Bandas | Blacklist (seca) | Zero-cloud (seca) |
|---|---|---|---|---|
| A | 2000-03-01 a 2008-04-01 (só cenas com B09) | VNIR + SWIR + TIR | 263 cenas | 34 cenas |
| B | 2009-01-01 a 2021-01-01 | VNIR + TIR | 782 cenas | 229 cenas |

**Regras por cena** (iguais ao JS): estação seca = jun–out de todos os anos do
período; `CLOUDCOVER < 60`; DN → radiância (GAIN_COEFFICIENT) → TOA 0–1; pixel com
VNIR DN ≥ 254 mascarado; máscara de nuvem/sombra só com VNIR + TIR; quality v5
(SAVI alto soma pontos). Cenas de `lista_imagens_Notaceites.txt` sempre saem.
Cenas de `list_zero_cloud.txt` não passam pela máscara de nuvem e ganham
`BONUS_ZERO_CLOUD = 3.0` no quality (o v5 vai até 2,1: sempre vencem); cena nas duas
listas fica fora.

**Etapas** (o usuário roda; todas interativas, com dry-run):

| # | Script | O que faz | Saída |
|---|---|---|---|
| 1 | `build_mosaico_sazonal_seco.py` | qualityMosaic seco de A e B por célula da grid 5 × 5 | `ASTER_sazonal/ASTER_sazonal_<A\|B>_seco_quality_<LxCy>` |
| 2 | `pontos_sazonal_seco_pareados.py` → etapa 1 | 2000 pontos/célula a 30 m onde A tem as 14 bandas | `points_sazonal_seco_A` |
| 3 | idem → etapa 2 | mesmos pontos amostrados no B (`sampleRegions`) | `points_sazonal_seco_B` |
| 4 | idem → baixar | CSV por tabela e junção pelo `pid` (sufixos `_A` / `_B`) | `src/Dados/ASTER_Samples_sazonal_seco/points_sazonal_seco_pareados.csv` |
| 5 | `../modeling/estimar_swir_rf_v2.py`, grupo `sazonal_seco` | RF (VNIR+TIR → SWIR) treinado nos pontos A, aplicado ao mosaico B | `ASTER_MOSAIC_v2/ASTER_MOSAIC_sazonal_seco_2009_2020_<LxCy>` |
| 6 | `build_mosaico_final_v2.py`, fonte `sazonal` | origem 1 = A (SWIR medido), origem 2 = B estimado só onde A não tem pixel | `aster-mosaico-final_sazonal/aster-mosaico-final-sazonal_<LxCy>` |

**Saída dos mosaicos (etapa 1):** reflectância × 1000 (B01–B09), TIR em DN, banda
`n_obs` (nº de cenas válidas no pixel), INT32, 15 m em EPSG:4326 com o
`CRS_TRANSFORM` ancorado no canto NW (células encaixam pixel a pixel), pirâmide
`mean`. Propriedades: `periodo`, `anos`, `meses`, `num_images`, `num_zero_cloud`,
`bonus_zero_cloud`, `blacklist`, `swir`…

**Normalização B → A (etapa 5):** o mosaico sazonal não tem banda `data` (cada
pixel vem de uma cena diferente), então não há reta por cena. Opções, perguntadas
na execução: `nenhuma` (padrão), `celula` (reta OLS por banda e célula nos pares;
célula com < 200 pares usa a global; coeficientes em
`src/modeling/models/sazonal_seco/normalizacao_ols_celula_coef.csv`; banda
`normalizacao` = 3) ou `global`. Ainda não avaliado: como A e B agora são da mesma
estação, espera-se um efeito da data da cena menor que no zero-cloud.

**Conta e destino:** etapas 1–4 rodam com `caatinga05` inicializada no projeto
`mapbiomas-brazil` (`CONTA` / `PROJETO` em `build_mosaico_sazonal_seco.py`);
assets em `projects/mapbiomas-arida/mine/` (`DESTINO`). Etapas 5–6 usam a conta
ativa / o rodízio de contas dos próprios scripts.

**Pendências:** exportar A e B de todas as células e conferir a cobertura do A;
índice de costura do mosaico produzido contra o v1 e o único; avaliar o SWIR
estimado do B nos pares (com e sem normalização por célula); testar o mosaico
final sazonal na modelagem de favorabilidade. Relatório:
`docs/relatorio_mosaico_aster/` (seção "Mosaico sazonal seco").

---

## Etapas do pipeline original (JS, legado)

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
| Mosaicos sazonais secos A/B (por célula) | `projects/mapbiomas-arida/mine/ASTER_sazonal/` |
| Pontos pareados sazonais | `projects/mapbiomas-arida/mine/points_sazonal_seco_A`, `..._B` |
| Mosaico B com SWIR estimado | `projects/mapbiomas-arida/mine/ASTER_MOSAIC_v2/ASTER_MOSAIC_sazonal_seco_2009_2020_<LxCy>` |
| Mosaico final sazonal | `projects/mapbiomas-arida/mine/aster-mosaico-final_sazonal/` |

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
