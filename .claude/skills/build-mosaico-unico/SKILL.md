---
name: build-mosaico-unico
description: Use when the user wants to build, run-plan, or adjust the single-script ASTER mosaic per grid cell (zero-cloud by minimum NDVI on top, other scenes only in gaps), split into com_swir / sem_swir assets, with the SWIR estimated for sem_swir in a second stage. Triggers on "mosaico único", "build_mosaico_unico", "com_swir", "sem_swir", "ASTER_mosaico_unico", "zero cloud + preenchimento", "estimar SWIR do mosaico sem swir".
---

# Build Mosaico Único (ASTER, com/sem SWIR)

Script: [src/building_dataset/build_mosaico_unico.py](../../../src/building_dataset/build_mosaico_unico.py).
Ver [SPEC.md](SPEC.md) para o contrato.

Importa `build_mosaicos_semestrais_v2 as b` (contas, máscaras, quality v6),
`build_mosaicos_zero_cloud_v2 as z` (`ndvi_cena`, lista zero-cloud) e
`modeling/estimar_swir_rf as rf` (treino, `estimar`). Mudanças de regra
nesses módulos afetam este script.

## Regras (não mudar sem pedido explícito)

Etapa 1 — dois assets por célula em
`projects/mapbiomas-arida/mine/ASTER_mosaico_unico_base`:

| Asset | Cenas | Bandas ópticas |
|---|---|---|
| `ASTER_base_com_swir_<LxCy>` | com `GAIN_COEFFICIENT_B04` (até 2008-S1) | VNIR + SWIR + TIR |
| `ASTER_base_sem_swir_<LxCy>` | sem (2008-S2 em diante) | VNIR + TIR |

Os dois com a mesma composição:
- camada 1: cenas zero-cloud, ordenadas por `ndvi_cena` (a mais seca por
  cima), TOA intacta exceto saturados (VNIR DN ≥ 254) → lacuna; quality
  `5000 + 500 × (1 − SAVI norm.)`;
- camada 2: demais cenas (CC < 77, diurnas, fora da blacklist), máscara de
  nuvem/sombra/saturação, `qualityMosaic` com quality v6; só nas lacunas;
- `ImageCollection([camada2, camada1]).mosaic()`;
- INT32: ópticas + `saturado`, `quality`, `data` (AAAAMMDD), `camada` (1/2).

Separar por SWIR evita misturar bandas de cenas diferentes no mesmo pixel
(o `mosaic()` compõe banda a banda). Filtro de SWIR é lazy
(`filtro_swir(com_swir)`): objetos `ee` não podem ser criados no import.

Etapa 2 — só o `sem_swir`: `rf.estimar(..., 'cena')` → B04–B09 estimados,
normalização OLS por cena via banda `data`. Saída
`projects/mapbiomas-arida/mine/ASTER_mosaico_unico/ASTER_sem_swir_estimado_<LxCy>`
(B01–B14, quality, data, camada, normalizacao). O `com_swir` não passa
pela etapa 2.

## Fluxo de uso (o usuário roda)

`python build_mosaico_unico.py` pergunta: etapa (1/2), contas, tipos
(ambos/com/sem), meses (jul–out padrão ou ano todo), células, conferência
ou exportação. Pula assets salvos e tasks ativas. Log
`registro_mosaico_unico_etapa<N>_<data>.csv`. Etapa 2 só depois da 1
terminar para a célula.

## Teste só leitura (referência, L3C3, jul–out)

com_swir: 3 zero-cloud + 14 preenchimento, cobertura 98,9 %, camada 1
17,5 %, NDVI 0,295. sem_swir: 22 + 56, 99,6 %, camada 1 91,5 %, NDVI 0,232.
Poucas zero-cloud secas antes de 2008 → com_swir depende da camada 2.

## Em aberto

Combinar com_swir e sem_swir_estimado num produto final (ex.: com_swir
onde NDVI parecido). Não implementar sem decisão do usuário.
Visualizador e `maxent_fosfato.py` (`ASSET_MOSAICO`) ainda não apontam
para estes assets.
