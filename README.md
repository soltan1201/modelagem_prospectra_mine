# Prospectra — modelagem mineral com ASTER (bloco Irecê, BA)

Pipeline em Google Earth Engine (Python `earthengine-api` + scripts JS do Code
Editor) para montar um mosaico ASTER multiespectral livre de nuvem sobre o bloco
de Irecê, estimar o SWIR que falta depois de 2008 e gerar mapas de
favorabilidade para fosfato (Maxent), com tratamento explícito da incerteza.

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
        SAT[Piloto saturação VNIR × Landsat<br/>saturados = nuvem]:::feito
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
        UB1[Mosaico único etapa 1<br/>com_swir 25/25 · sem_swir 21/25]:::andamento
        UB2[Mosaico único etapa 2<br/>sem_swir_estimado 21/25]:::andamento
        DEC[Mosaico único etapa 3<br/>com_swir + sem_swir, SWIR re-estimado<br/>retas OLS entre os mosaicos + RF<br/>ASTER_mosaico_unico_final]:::andamento
    end

    subgraph E4["4 · Dados complementares (Acervo SGB/CPRM)"]
        ACV[Catálogo + download + .grd → GeoTIFF]:::feito
        ING[acervo_cprm no GEE<br/>gama eU, mag, litologia,<br/>estruturas, ocorrências, P2O5]:::feito
        STK[Pilha de covariáveis por célula<br/>aerogeofísica + geologia + MDE]:::afazer
    end

    subgraph E5["5 · Favorabilidade e incerteza"]
        HIP[Hipóteses de depósito H1–H3]:::afazer
        MAX[maxent_fosfato.py<br/>amostras → validação por célula → mapa]:::andamento
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
    FIN --> MAX
    DEC --> MAX
    ACV --> ING --> STK --> MAX
    ING -- ocorrências fosfato --> MAX
    HIP --> MAX --> ENS --> RSK --> CMP --> ALV
    CMP -- presenças / ausências --> ENS
```

O ciclo da etapa 5 (ensemble → confiança → campo → reassimilação) segue o plano
de [docs/analise_fluxo_prospeccao_e_plano_irece.md](docs/analise_fluxo_prospeccao_e_plano_irece.md)
(fases 0–5).

## Estado (2026-10-17)

| Coleção em `projects/mapbiomas-arida/mine/` | Script | Situação |
|---|---|---|
| `ASTER_v4` | `build_mosaicos_semestrais.py` | 548 semestrais por célula |
| `ASTER_0cloud` | `build_mosaicos_zero_cloud.py` | 42 (células sem cenas puladas) |
| `ASTER_MOSAIC` | `modeling/estimar_swir_rf.py` | 50 (zero_cloud + semestral_v4, 2009–2020) |
| `aster-mosaico-final` | `build_mosaico_final.py` | 25/25 |
| `ASTER_v5`, `ASTER_0cloud_v2`, `ASTER_MOSAIC_v2`, `aster-mosaico-final_v2` | scripts `_v2` | não existem (nunca exportados) |
| `ASTER_mosaico_unico_base` | `build_mosaico_unico.py` etapa 1 | com_swir 25/25; sem_swir 21/25 (L3C2 rodando, L3C3–C5 na fila em solkanCengine) |
| `ASTER_mosaico_unico` | `build_mosaico_unico.py` etapa 2 | 21/25 (faltam L3C2–L3C5) |
| `ASTER_mosaico_unico_final` | `build_mosaico_unico.py` etapa 3 | código pronto e testado (só leitura); não exportado |
| `acervo_cprm/` | `dowloadQGIS/exportar_acervo_gee.py` | 10 assets |
| `pontos_fosfato`, `maxent_fosfato_amostras_teste_zero_cloud_2000_2008` | `modeling/maxent_fosfato.py` | teste sem presenças; refazer com mosaico final |

## Estrutura

```
src/
├── building_dataset/    mosaicos (semestrais, zero-cloud, final, único), listas de cenas,
│                        pontos pareados, amostragem, análises de histograma e saturação
├── modeling/            SWIR (notebook de treino, estimar_swir_rf*.py), Maxent fosfato,
│                        README.md = diário dos testes de modelagem
├── dowloadQGIS/         Acervo SIG SGB/CPRM: catálogo, download, .grd → tif, envio ao GEE
│                        (ACERVO_CPRM_IRECE.md = documento central)
├── features_process/    biblioteca de índices espectrais (JS)
├── show_mosaics/        visualizadores do Code Editor
├── tools/               contas, tasks, cópia/remoção de assets, cota
├── revision_blibiografic/  triagem Scopus
├── Dados/               saídas locais (CSV, figuras) — maioria fora do git
├── gee_tools.py         troca de conta, tasks
└── configure_account_projects_ee.py   conta → projeto GEE
docs/                    plano de prospecção sob incerteza + textos de referência
resumo_metodologia.md    resumo da metodologia
```

Os scripts de export são interativos (perguntam contas, células, conferir ou
exportar), pulam assets já salvos e tasks ativas, e gravam
`registro_<etapa>_<data>.csv` na pasta do script.

## Skills do Claude Code (`.claude/skills/`)

| Skill | Cobre |
|---|---|
| `build-mosaico-semestral` | semestrais v4/v5, zero-cloud, blacklist, quality score |
| `build-mosaico-unico` | mosaico único com/sem SWIR, etapas 1 e 2 |
| `merge-mosaico-final` | mosaico final por prioridade (v1/v2) |
| `estimate-swir-gee` | SWIR por RF por célula (e o GBT legado) |
| `train-swir-models` | treino local PLS/XGBoost/GPR (script antigo) |
| `sample-mosaico-dataset` | amostragem de pontos a 30 m |
| `add-indice-espectral` | novos índices espectrais |
| `acervo-cprm` | dados SGB/CPRM até o GEE |
| `gee-contas-assets` | rodízio de contas, tasks, cópia e remoção de assets |
| `triagem-bibliografica` | filtro de artigos Scopus |

Cada skill tem `SKILL.md` (como usar) e `SPEC.md` (contrato).

## Requisitos

- Python 3 com `earthengine-api`, `pandas`, `numpy`, `scikit-learn`, `joblib`
  (o `xgboost` só para o notebook de treino).
- Credenciais GEE por conta em `~/.config/earthengine/<conta>/credentials`
  (ver `src/configure_account_projects_ee.py`). A superconta
  (`mapbiomas-brazil`) não lê nem escreve `mapbiomas-arida`.
- Envio do acervo ao GEE: bucket `gs://cengine17-acervo-cprm`
  (`src/dowloadQGIS/GUIA_GCLOUD_BUCKET_GEE.md`).

## Licença

Ver [LICENSE](LICENSE).
