# Prospectra — modelagem mineral com ASTER (bloco Irecê, BA)

Pipeline em Google Earth Engine (Python `earthengine-api` + scripts JS do Code
Editor) para montar um mosaico ASTER multiespectral livre de nuvem sobre o bloco
de Irecê, estimar o SWIR que falta depois de 2008 e gerar mapas de
favorabilidade para fosfato (Maxent, Random Forest, Gradient Tree Boosting),
com tratamento explícito da incerteza.

- **Área:** retângulo -44.4605 / -14.7908 / -39.2035 / -9.4165, processado por
  célula da grid 5 × 5 (`grid_irece_5x5_solape100`, células `L1C1`…`L5C5`).
- **Sensor:** ASTER L1T (2000–2020). O SWIR (B04–B09) parou em 2008-S2; de
  2009 em diante só há VNIR + TIR.
- **Convenção dos dados:** VNIR e SWIR em reflectância TOA × 1000, TIR em DN,
  INT32, 15 m, EPSG:4326. `saturado = 1` marca VNIR com DN ≥ 254.
- **Assets:** `projects/mapbiomas-arida/mine/` (todas as contas exceto a
  superconta escrevem lá).

Lista de tarefas para o Kanri (a fazer / em processo / feitas):
[docs/kanban_kanri.md](docs/kanban_kanri.md).

## Fluxograma das etapas

Cores: verde = concluído, amarelo = em processo, cinza = a fazer, branco
tracejado = alternativa/legado.

```mermaid
flowchart TD
    classDef feito fill:#c8e6c9,stroke:#2e7d32,color:#1b1b1b
    classDef andamento fill:#fff3c4,stroke:#b58900,color:#1b1b1b
    classDef afazer fill:#eceff1,stroke:#607d8b,color:#1b1b1b
    classDef alternativa fill:#ffffff,stroke:#90a4ae,stroke-dasharray:4 3,color:#1b1b1b

    subgraph E0["0 · Base"]
        L1T[ASTER L1T 2000–2020]:::feito
        GRID[Grid 5×5 + área de estudo]:::feito
        LISTAS[Blacklist + lista zero-cloud<br/>visualizar_zero_cloud.js]:::feito
        SAT[Piloto saturação VNIR × Landsat<br/>correção pior que o cortado em 6/8<br/>saturados = nuvem → marcar e excluir]:::feito
    end

    subgraph E1["1 · Mosaicos v1 (ano todo)"]
        V4[Semestrais v4<br/>ASTER_v4]:::feito
        ZC[Zero-cloud 00–08 / 09–20<br/>ASTER_0cloud]:::feito
    end

    subgraph E2["2 · Modelo de SWIR"]
        PAR[Pontos pareados zero-cloud<br/>00–08 × 09–20]:::feito
        NB[Notebook PLS × XGBoost × RF<br/>RF escolhido + norm. OLS por cena]:::feito
        TST[Testes pendentes: época, mesma estação,<br/>sem TIR, validação espacial]:::afazer
        RF[estimar_swir_rf.py<br/>ASTER_MOSAIC 2009–2020]:::feito
    end

    subgraph E3["3 · Mosaico final"]
        FIN[build_mosaico_final.py<br/>aster-mosaico-final 25/25<br/>origem 1–4 + swir_estimado]:::feito
        V2[Cadeia v2 seco jul–out<br/>ASTER_v5 · 0cloud_v2 · MOSAIC_v2 · final_v2<br/>scripts prontos, não exportados]:::alternativa
        UB1[Mosaico único etapa 1<br/>com_swir 25/25 · sem_swir 25/25]:::feito
        UB2[Mosaico único etapa 2<br/>sem_swir_estimado 21/25]:::andamento
        DEC[Mosaico único etapa 3<br/>com_swir + sem_swir, SWIR re-estimado<br/>retas entre os mosaicos + RF<br/>ASTER_mosaico_unico_final 25/25]:::feito
    end

    subgraph E4["4 · Dados complementares (Acervo SGB/CPRM)"]
        ACV[Catálogo + download + .grd → GeoTIFF]:::feito
        ING[acervo_cprm no GEE<br/>gama eU, mag, litologia,<br/>estruturas, ocorrências, P2O5]:::feito
        STK[Pilha de covariáveis por célula<br/>aerogeofísica + geologia + MDE]:::afazer
    end

    subgraph E5["5 · Favorabilidade e incerteza"]
        HIP[Hipóteses de depósito H1–H3]:::afazer
        PTS[pontos_referencia<br/>fosfato + 4 alvos · fonte · papel]:::feito
        AMS[Amostras nos 2 mosaicos<br/>maxent_amostras_*]:::feito
        AVA[avaliar_maxent.ipynb<br/>Maxent × RF × GTB, célula de fora]:::feito
        SHP[SHAP fora da amostra<br/>21 de 36 variáveis]:::feito
        MAX[exportar_favorabilidade_ml.py<br/>RF/GTB → favorabilidade_ml<br/>fosfato v1 25/25 cada]:::feito
        DL[CNN · model_DL/<br/>patches 64 px · modelos A índices e B1 log-bandas<br/>scripts prontos, não rodados]:::andamento
        ENS[Ensemble Maxent<br/>P média, σ, P10/P50/P90, entropia]:::afazer
        RSK[Mapa de confiança +<br/>curva área-favorabilidade]:::afazer
        CMP[Campo em rodadas<br/>alvos P≈0,5 e σ alto]:::afazer
        ALV[Relatório de alvos]:::afazer
    end

    L1T --> V4 & ZC
    LISTAS --> V4 & ZC & UB1
    GRID --> V4 & ZC & UB1
    SAT -.-> UB1
    ZC --> PAR --> NB --> RF
    NB -.-> TST -.-> RF
    V4 --> RF
    ZC --> FIN
    V4 --> FIN
    RF --> FIN
    L1T -.-> V2 -.-> DEC
    L1T --> UB1 --> UB2 --> DEC
    UB1 -- com_swir --> DEC
    RF -. rf.estimar .-> UB2
    ING -- ocorrências --> PTS
    PTS --> AMS
    FIN --> AMS
    DEC --> AMS
    AMS --> AVA --> SHP --> MAX
    MAX -- pseudo-rótulos --> DL
    ACV --> ING --> STK --> MAX
    HIP --> MAX --> ENS --> RSK --> CMP --> ALV
    CMP -- presenças / ausências --> ENS
```

O ciclo da etapa 5 (ensemble → confiança → campo → reassimilação) segue o plano
de [docs/analise_fluxo_prospeccao_e_plano_irece.md](docs/analise_fluxo_prospeccao_e_plano_irece.md)
(fases 0–5).

## Estado (2026-10-22)

| Coleção em `projects/mapbiomas-arida/mine/` | Script | Situação |
|---|---|---|
| `ASTER_v4` | `build_mosaicos_semestrais.py` | 548 semestrais por célula |
| `ASTER_0cloud` | `build_mosaicos_zero_cloud.py` | 42 (células sem cenas puladas) |
| `ASTER_MOSAIC` | `modeling/estimar_swir_rf.py` | 50 (zero_cloud + semestral_v4, 2009–2020) |
| `aster-mosaico-final` | `build_mosaico_final.py` | 25/25 |
| `ASTER_v5`, `ASTER_0cloud_v2`, `ASTER_MOSAIC_v2`, `aster-mosaico-final_v2` | scripts `_v2` | não existem (nunca exportados) |
| `ASTER_mosaico_unico_base` | `build_mosaico_unico.py` etapa 1 | com_swir 25/25; sem_swir 25/25 |
| `ASTER_mosaico_unico` | `build_mosaico_unico.py` etapa 2 | 21/25 (faltam L3C2–L3C5; a etapa 3 não depende dela) |
| `ASTER_mosaico_unico_final` | `build_mosaico_unico.py` etapa 3 | 25/25 |
| `acervo_cprm/` | `dowloadQGIS/exportar_acervo_gee.py` | 10 assets |
| `pontos_referencia` | `modeling/preparar_pontos_referencia.py` | 1.103 pontos: fosfato 47 (antigo+cprm / cprm; treino / validação), carbonato 327, mvt 252, ouro 208, manganês 269 |
| `maxent_amostras_<mosaico>_bandas_indices` | `modeling/maxent_fosfato.py` etapa 1 | 2 (único e v1), 31.860 amostras cada |
| `favorabilidade_ml` | `modeling/exportar_favorabilidade_ml.py` | fosfato no mosaico v1: RF 25/25 e GTB 25/25 (21 variáveis do SHAP); visualizador `show_mosaics/visualizar_favorabilidade_ml.js` |
| `pontos_fosfato`, `maxent_fosfato_amostras_teste_zero_cloud_2000_2008` | — | antigos; substituídos pelos de cima |

### Avaliação dos modelos (`src/modeling/avaliar_maxent.ipynb`)

AUC por ponto, validação deixando uma célula da grid de fora (fundo aleatório):

| alvo | único Maxent | único RF | único GTB | v1 Maxent | v1 RF | v1 GTB |
|---|---|---|---|---|---|---|
| fosfato | 0,58 | 0,61 | 0,63 | 0,70 | 0,83 | 0,84 |
| carbonato | 0,60 | 0,66 | 0,67 | 0,64 | 0,64 | 0,67 |
| mvt | 0,68 | 0,71 | 0,72 | 0,61 | 0,68 | 0,67 |
| ouro | 0,34 | 0,46 | 0,50 | 0,56 | 0,63 | 0,65 |
| manganês | 0,64 | 0,70 | 0,71 | 0,64 | 0,62 | 0,66 |

- RF/GTB ganham do Maxent na maioria; dentro da amostra os dois passam de 0,98
  (decoram o lugar) — o número que vale é o de célula de fora.
- p(outros alvos) como covariável não melhora o fosfato.
- SHAP fora da amostra (RF + GTB, 20 combinações): 21 de 36 variáveis sem
  pares |r| > 0,9 (`FEATURES_DL`); topo: `quartz_14_12`, `B10`,
  `ferric_iron_2_1`, `glcm_contrast`, `B09`, `sultan_g`, `phosphorite_rbd`.
- Mapas RF/GTB do fosfato (v1): concordam nas áreas de p alto em volta dos
  pontos e mostram manchas novas, sem pontos, no sul e no centro da área. As
  áreas com p ≥ 0,5 nos dois modelos são as candidatas a pseudo-amostras do
  deep learning (descontando as manchas com borda de cena).

### Relatórios (LaTeX + PDF em `docs/`)

| Pasta | Tema | Figuras |
|---|---|---|
| [relatorio_mosaico_aster/](docs/relatorio_mosaico_aster/) | sensor, disponibilidade, nuvens, saturação do VNIR (teste com Landsat), SWIR estimado, harmonização, mosaicos v1 e único, índices | `docs/figuras_contagem/` |
| [relatorio_avaliacao_modelos/](docs/relatorio_avaliacao_modelos/) | Maxent × RF × GTB, validação espacial, SHAP, mapas RF/GTB | `relatorio_avaliacao_modelos/image/` |
| [relatorio_acervo_cprm/](docs/relatorio_acervo_cprm/) | aquisição do Acervo SIG SGB/CPRM até o GEE (4 páginas) | `docs/figuras_acervo/` |

Compilar com `latexmk -pdf <arquivo>.tex` dentro da pasta.

### CNN de favorabilidade (`src/model_DL/`, plano em `plano_cnn.md`)

| Script | Faz |
|---|---|
| `exportar_patches.py` | sorteia centros (ocorrências, pseudo-positivos onde RF e GTB têm p ≥ 0,5, negativos com p ≤ 0,05), divide por célula (teste = células das ocorrências de validação; 5 folds nas demais) e baixa patches 64 × 64 px do mosaico v1 por `computePixels` (sem task) em `.npz`, um arquivo por modelo |
| `treinar_cnn.py` | PyTorch + timm: CNN pequena ou RegNetY-016 / EfficientNetV2-S / ResNet-50 / CoAtNet-0; BCE com rótulo suave; validação cruzada por fold; `auc_ponto` contra o GTB (0,84); `--teste` no fim |
| `inferir_cnn.py` | janela deslizante por célula (passo 3 = 90 m ou 1 = 30 m), ensemble dos folds, GeoTIFF de p por célula |
| `cnn_colab.ipynb` | os três scripts no Colab (importa os mesmos `.py`): Drive, GEE, exportação, treino, comparação das rodadas, teste, inferência e mapa |

Modelo **A** = 15 índices espectrais + `swir_estimado` (16 canais); modelo
**B1** = log da reflectância B01–B09 + emissividade B10–B14 + `swir_estimado`
(15 canais). Dados no servidor em `/srv/almacen/db_images/patches_cnn` e
rodadas em `/srv/almacen/db_images/runs_cnn` (variáveis `PATCHES_DIR`,
`RUNS_DIR`).

## Estrutura

```
src/
├── building_dataset/    mosaicos (semestrais, zero-cloud, final, único), listas de cenas,
│                        pontos pareados, amostragem, análises de histograma e saturação
├── modeling/            SWIR (notebook de treino, estimar_swir_rf*.py); favorabilidade:
│                        preparar_pontos_referencia.py, maxent_fosfato.py (amostras,
│                        Maxent por alvo/mosaico), maxent_alvos.js (prévia),
│                        avaliar_maxent.ipynb (Maxent × RF × GTB + SHAP),
│                        exportar_favorabilidade_ml.py (RF/GTB no GEE);
│                        README.md = diário dos testes de modelagem
├── model_DL/            CNN: exportar_patches.py, treinar_cnn.py, inferir_cnn.py, plano_cnn.md,
│                        cnn_colab.ipynb (mesmos scripts no Colab)
├── dowloadQGIS/         Acervo SIG SGB/CPRM: catálogo, download, .grd → tif, envio ao GEE
│                        (ACERVO_CPRM_IRECE.md = documento central)
├── features_process/    biblioteca de índices espectrais (JS); exportar_mosaicos_indices
│                        (.js / .py / _colab.ipynb) = bandas + índices para o Drive
├── show_mosaics/        visualizadores do Code Editor (mosaicos, comparação v1 × único,
│                        favorabilidade RF/GTB, mapa de localização)
├── tools/               contas, tasks, cópia/remoção de assets, cota
├── revision_blibiografic/  triagem Scopus
├── Dados/               saídas locais (CSV, figuras) — maioria fora do git
├── gee_tools.py         troca de conta, tasks
└── configure_account_projects_ee.py   conta → projeto GEE
docs/                    plano de prospecção sob incerteza, kanban, pontos de referência,
                         relatórios LaTeX (mosaico, avaliação, acervo), pastas de figuras
                         + textos de referência
resumo_metodologia.md    resumo da metodologia
```

Os scripts de export são interativos (perguntam contas, células, conferir ou
exportar), pulam assets já salvos e tasks ativas, e gravam
`registro_<etapa>_<data>.csv` na pasta do script.

## Skills do Claude Code (`.claude/skills/`)

| Skill | Cobre |
|---|---|
| `build-mosaico-semestral` | semestrais v4/v5, zero-cloud, blacklist, quality score |
| `build-mosaico-unico` | mosaico único com/sem SWIR, etapas 1–3 |
| `merge-mosaico-final` | mosaico final por prioridade (v1/v2) |
| `estimate-swir-gee` | SWIR por RF por célula (e o GBT legado) |
| `train-swir-models` | treino local PLS/XGBoost/GPR (script antigo) |
| `sample-mosaico-dataset` | amostragem de pontos a 30 m |
| `add-indice-espectral` | novos índices espectrais |
| `acervo-cprm` | dados SGB/CPRM até o GEE |
| `modelagem-favorabilidade` | pontos de referência, amostras, Maxent × RF × GTB, SHAP, `favorabilidade_ml`, relatório |
| `cnn-favorabilidade` | CNN em `src/model_DL/`: patches, pseudo-rótulos, divisão por célula, treino, inferência, Colab |
| `gee-contas-assets` | rodízio de contas, tasks, cópia e remoção de assets |
| `triagem-bibliografica` | filtro de artigos Scopus |

Cada skill tem `SKILL.md` (como usar) e `SPEC.md` (contrato).

## Requisitos

- Python 3 com `earthengine-api`, `pandas`, `numpy`, `scikit-learn`, `joblib`,
  `scipy`, `matplotlib` (o `xgboost` só para o notebook de treino; o `shap` só
  para a seção 9 do `avaliar_maxent.ipynb`).
- CNN (`src/model_DL/`): `torch` (com a CUDA da GPU), `timm`; `rasterio` ou GDAL
  para gravar os GeoTIFFs da inferência.
- Credenciais GEE por conta em `~/.config/earthengine/<conta>/credentials`
  (ver `src/configure_account_projects_ee.py`). A superconta
  (`mapbiomas-brazil`) não lê nem escreve `mapbiomas-arida`.
- Envio do acervo ao GEE: bucket `gs://cengine17-acervo-cprm`
  (`src/dowloadQGIS/GUIA_GCLOUD_BUCKET_GEE.md`).

## Licença

Ver [LICENSE](LICENSE).
