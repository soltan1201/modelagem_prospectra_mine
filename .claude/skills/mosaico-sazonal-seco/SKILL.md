---
name: mosaico-sazonal-seco
description: Use when the user works on the ASTER seasonal dry mosaic for Irecê (periods A 2000-03–2008-03 with SWIR and B 2009–2020 without) — building it, copying it brazil → arida, viewing it, harmonizing scene seams against MODIS MCD43A4 (scene map, ASTER × MODIS pairs to Drive, per-scene lines, apply), paired points A × B, the sazonal notebook, or the step checklist. Triggers on "mosaico sazonal", "sazonal seco", "ASTER_sazonal", "período A/B", "build_mosaico_sazonal_seco", "harmonizar", "harmonização MODIS", "MCD43A4", "mapa de cenas", "cena_t", "confere", "pares_modis", "coef_modis", "ASTER_sazonal_harmonizado", "pontos_sazonal_seco_pareados", "train_swir_models_sazonal", "CHECKLIST_sazonal_seco".
---

# Mosaico sazonal seco ASTER (A × B) e harmonização

Ver [SPEC.md](SPEC.md). Roteiro de comandos para o usuário:
[README_guia_sazonal_seco.md](../../../src/building_dataset/README_guia_sazonal_seco.md).
Andamento (o usuário marca): [CHECKLIST_sazonal_seco.md](../../../src/building_dataset/CHECKLIST_sazonal_seco.md)
— ao receber resultado de uma etapa, atualize o checklist.

O usuário roda tudo (scripts, notebook, tasks). Claude escreve código e lê os logs colados.

## Períodos e mosaico (`src/building_dataset/build_mosaico_sazonal_seco.py`)

- A: 2000-03 a 2008-03, VNIR + SWIR + TIR (só cenas com B09). B: 2009–2020, VNIR + TIR.
- Estação seca jun–out de todos os anos; `CC_MAX = 60`; blacklist
  `lista_imagens_Notaceites.txt`; zero-cloud `list_zero_cloud.txt` (sem máscara de
  nuvem, `BONUS_ZERO_CLOUD = 3.0`); quality v5; `qualityMosaic` por célula 5 × 5.
- Saída: TOA × 1000, TIR em DN, banda `n_obs`, INT32, `CRS_TRANSFORM` da grid.
  Nome `ASTER_sazonal_<A|B>_seco_quality_<LxCy>`. Não tem banda de data.
- `colecao_seca(periodo, area, blacklist, zero)` devolve `(base, filtrada, cc0, col)`;
  `montar_seco` usa ela. Reaproveitada pela harmonização — não mudar filtros/ordem
  sem refazer os mapas de cenas.
- Conta superconta, projeto mapbiomas-brazil → export temporário em
  `.../CAATINGA/ASTER_sazonal`; cópia para `projects/mapbiomas-arida/mine/ASTER_sazonal`
  com `src/tools/transferir_assets_brazil_para_arida.py` (caatinga04).
- Visualizador: `src/show_mosaics/visualizar_mosaico_sazonal_seco.js`,
  `FONTE = 'arida' | 'brazil' | 'harmonizado'`.
  Comparação lado a lado: `src/show_mosaics/visualizar_harmonizado_sazonal.js`
  (mapas ligados, faixas da direita = faixas da esquerda × reta mediana; camadas
  bordas, ganho aplicado, cena por `cena_t`, `harmonizado = 2`; clique com data).

Problema conhecido: emendas por cena (colcha de retalhos, faixa ciano na borda,
pontos brancos). Causa: bônus de cena do quality v5 (zero-cloud, cloud_bonus) domina
o termo por pixel → cenas inteiras vencem. Saídas discutidas: remontar (pixel-driven
quality, medoide, janela ago–out) ou harmonizar o que já existe (escolhida: S1).

## Harmonização S1 (`src/building_dataset/harmonizar_sazonal_modis.py`)

caatinga04 no projeto `PROJETO = 'mapbiomas-brazil'`. Menu:

1. **mapa de cenas** (task por célula) → `mine/ASTER_sazonal_cenas`: `cena_t`
   (segundos desde 2000-01-01, `T0_MS`) e `confere` (1 = reconstrução igual ao salvo).
2. **pares → Drive** (task `pares_modis_<P>_<ano>_<mês>` por mês; A 40, B 60) →
   `Drive/ASTER_sazonal_harmonizacao/<P>_<ano>_<mês>.csv`. Por cena: só TOA (sem
   máscara de nuvem), ASTER agregado 15/30 → 90 → 500 m na grade MCD43A4, pixels
   homogêneos (`CV_MAX`, `FRAC_MIN`) e QA = 0, `stratifiedSample` em 10 faixas iguais
   do vermelho MODIS 0–10000 × 200 px.
3. **retas (local)**: lê `src/Dados/ASTER_sazonal_harmonizacao/ASTER_sazonal_harmonizacao/*.csv`, tira
   nuvem (`DIF_NUVEM`) e sombra (`DIF_SOMBRA`), OLS em 2 passadas (corte 3 × RMSE) →
   `coef_modis_<P>.csv`. Refazível sem GEE.
4. **aplicar** → `mine/ASTER_sazonal_harmonizado` (mesmos nomes): ganho/offset por
   `remap` de `cena_t` (× 1e6, inteiros); reta inválida, cena sem pares ou
   `confere = 0` → mediana do período; bandas extras `cena_t`, `harmonizado` (1/2).
5. **status**: fração de `confere` por célula (amostra 2000 px) e resumo das retas.

Pares de bandas: B01↔MODIS 4, B02↔1, B3N↔2, B04↔6 (ganho + offset); B05–B09 um ganho
só, sem offset, ajustado em B05↔7 (preserva razões dos índices minerais). TIR não é
corrigido. Saída em escala NBAR (superfície × 1000) — A e B precisam ir juntos.

Lições (não repetir):
- tabela de todas as cenas numa task → out of memory; `toAsset` com
  `Feature(None)` → "Unable to export features with null geometry" (Drive CSV aceita);
- `getInfo` por cena com a máscara de nuvem do build → > 1 h por mês: o `focalMax`
  da sombra (até 2,3 km a 15 m) é o gargalo, e `stratifiedSample` calcula a região
  inteira de qualquer jeito. Por isso: sem máscara + batch toDrive + ajuste local.

Estado (2026-10-30): etapa 1 concluída (todas as células existentes); etapa 2 nova
(batch Drive) ainda não rodada; falta B L2C1 no mosaico original.

## 2ª passada relativa (`src/building_dataset/harmonizar_sazonal_relativo.py`)

Para células que ainda destoam depois da S1 (fenologia da data: a reta contra o MODIS
do mesmo dia não corrige cena verde). Nas células de referência (padrão L1C2, L1C3,
L2C2, L2C3) ajusta harmonizado = a × MCD43A4 mediano seco + b; em cada cena das
células alvo (padrão L1C1, L2C1) ajusta z = a·m + b contra o harmonizado → reta por
cena, aplicada só nas alvo. Etapas: 1 amostras → `Drive/ASTER_sazonal_harmonizacao_rel/rel_<P>.csv`,
2 retas locais → `coef_rel_<P>.csv`, 3 aplicar → `mine/ASTER_sazonal_harm_rel`
(banda extra `ajuste_rel`). Visualizador junta as coleções com `USAR_REL`.

## Depois da harmonização (fonte = harmonizado Landsat, decidido 2026-11-01)

Harmonização adotada: `harmonizar_sazonal_landsat.py` (README_harmonizacao_landsat.md) →
`mine/ASTER_sazonal_hm_razao` e `_hm_residuo`, mesmos nomes do original. A fonte é
`FONTE` em `pontos_sazonal_seco_pareados.py` (padrão `'hm_razao'`; `'original'` volta ao
TOA); o sufixo (`_hm_razao`) vai em tudo que vem depois, para não misturar:

- `pontos_sazonal_seco_pareados.py` (caatinga04): `points_sazonal_seco<SUF>_A/_B`, pasta
  `src/Dados/ASTER_Samples_sazonal_seco<SUF>/` (CSVs `points_sazonal_seco_A/_B/_pareados.csv`).
- `src/modeling/train_swir_models_sazonal.ipynb`: `SUFIXO_FONTE` na célula 2 →
  `models/sazonal_seco<SUF>/` (`best_params.json`, `best_params_pareado.json`).
- `estimar_swir_rf_v2.py`, grupos `sazonal_seco`/`sazonal_pareado` (skill
  `estimate-swir-gee`): lê `ps.COL_FONTE` e `DIR_SAZONAL`; saída
  `ASTER_MOSAIC_sazonal_seco<SUF>_2009_2020_<LxCy>`.
- `build_mosaico_final_v2.py`, fonte `sazonal_hm`: A por cima; B onde falta A (origem 2)
  ou onde o A parece nuvem (origem 3: B01 do A − B01 do B > `DIF_NUVEM_B01` E B13 do B −
  B13 do A > `DIF_NUVEM_B13`). Só funciona porque A e B estão na mesma escala (superfície).
  Saída `mine/aster-mosaico-final_sazonal<SUF>/aster-mosaico-final-sazonal<SUF>_<LxCy>`.
