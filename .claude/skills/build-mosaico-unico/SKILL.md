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

Etapa 3 — mosaico único juntando os dois bases (decidido em 2026-10-17):
`projects/mapbiomas-arida/mine/ASTER_mosaico_unico_final/ASTER_mosaico_unico_<LxCy>`.
- SWIR do sem_swir estimado de novo, como no notebook: reta OLS por banda
  VNIR/TIR (2009–2020 → 2000–2008) + RF da etapa 2. Retas `mosaicos`
  (padrão): pares pixel a pixel sem_swir → com_swir na célula, |ΔNDVI| ≤
  `LIMIAR_NDVI` (0,05), por cena (`data`, ≥ `MIN_PIX_CENA` = 200 pares a
  `ESCALA_AJUSTE` = 90 m), senão reta da célula (vai em `coefs['global']`,
  aplicada por `rf.normalizar(..., 'cena')`). Opção `pareados` = retas da
  etapa 2.
- Retas calculadas por `getInfo` ANTES do export (`retas_mosaicos`) e
  gravadas em `retas_mosaico_unico_<data>.csv`; dentro da imagem, reduções
  aninhadas dão "Too many concurrent aggregations".
- Vence no pixel: camada 1 sobre camada 2; empate → menor NDVI. com_swir
  vencedor = SWIR medido.
- Bandas: B01–B14, quality, data, camada, `swir_estimado` (0/1),
  `normalizacao` (0 medido, 1 célula, 2 cena).
- Reta `METODO_RETA` 'rma' (padrão, ganho σ_com/σ_sem) ou 'ols'; pares
  `PARES_RETA` 'todos' (padrão) ou 'zero_cloud' (camada 1 nos dois —
  célula sem pares → `sem_retas`). Teste L2C2/L1C1: RMA leva ganhos de
  0,53–0,71 a 0,83–1,12 mas o viés do SWIR quase não muda (L2C2 ≈ −9 %,
  L1C1 +3/+11 %); zero_cloud+RMA em L2C2 −4/−7 % com só 6 mil pares. O
  viés restante é do RF (sinal muda por célula), não da reta.
- Teste só leitura L2C2 (recorte de 40 km, 300 m): viés do SWIR estimado
  × medido −7/−10/−9 % (B04/B06/B08) com retas `mosaicos` contra
  −23/−32/−33 % com `pareados`; 7 de 8 cenas com reta própria; ganhos
  OLS da célula 0,61–0,76.

## Superconta

Menu "Contas": rodízio arida (→ `mine/`) ou superconta (→
`CAATINGA_BR/<coleção>`, mapbiomas-brazil). A superconta não lê nem
escreve arida: serve para a etapa 1; o script liga `all_users_can_read`
na coleção nova no brazil. Bases são procuradas nos dois projetos
(`listar_bases`, arida primeiro); dedup também olha a mesma coleção no
outro projeto. Etapa 3 deve rodar no rodízio arida (lê as duas).

Projeto EE: `usar_conta` do script (não `b.usar_conta`) inicializa TODAS as
contas em `PROJETO_EE = 'mapbiomas-brazil'` (cota de tasks do projeto);
destinos não mudam. `CONTAS_PROJETO_EE`: caatinga01–05, solkan1201,
superconta (testado 2026-10-17: inicializam e leem/escrevem arida).
solkanCengine (sem serviceusage no projeto) e solkanGeodatin (não lê
arida por ele) são tiradas do rodízio com aviso. Passos
atuais no topo do docstring do script.

## Fluxo de uso (o usuário roda)

`python build_mosaico_unico.py` pergunta: etapa (1/2/3), retas (etapa 3), contas, tipos
(ambos/com/sem), meses (jul–out padrão ou ano todo), células, conferência
ou exportação. Pula assets salvos e tasks ativas. Log
`registro_mosaico_unico_etapa<N>_<data>.csv`. Etapa 2 só depois da 1
terminar para a célula.

## Diagnóstico: etapa 2 pula tudo com `sem_base`

Log `INFO: ASTER_sem_swir_estimado_<LxCy>: sem mosaico base sem_swir (rode a
etapa 1) — pulado` e `ASTER_mosaico_unico` vazio = falta
`ASTER_base_sem_swir_<LxCy>` em `ASTER_mosaico_unico_base`. Não é bug: a
etapa 2 só lê o base sem_swir. Conferir nesta ordem:
1. coluna `tipo` do `registro_mosaico_unico_etapa1_*.csv` — se só há
   `com_swir` (nem `sem_cenas`), a etapa 1 rodou com "só com_swir";
2. listar o destino (`b.listar_destino(DESTINO_BASE, criar=False)`);
3. tasks ainda `RUNNING` nas contas (base ainda não salvo).
Correção: etapa 1 com "os dois" ou "só sem_swir" (com_swir salvos são
pulados), esperar as tasks, depois etapa 2.

## Estado (2026-10-17)

`ASTER_mosaico_unico_base`: 25/25 `com_swir`; 21/25 `sem_swir` (faltam
L3C2–L3C5, tasks em solkanCengine: L3C2 `RUNNING`, L3C3–C5 `READY`).
solkanCengine ainda tem `READY` para L4C2–L5C5 sem_swir, que já existem —
duplicatas que vão falhar; sugerir ao usuário cancelar (não cancelar sem
pedido).
`ASTER_mosaico_unico`: 21/25 `sem_swir_estimado` (L3C1 feito na rodada
de 2026-10-16 23:50, conta caatinga05; L3C2–C5 `sem_base`). Pendente:
etapa 2 para L3C2–L3C5 após os bases terminarem.

Conferir estado sem trocar conta: `listAssets` nos dois destinos com
credencial de caatinga05 (ver skill `gee-contas-assets`).

## Teste só leitura (referência, L3C3, jul–out)

com_swir: 3 zero-cloud + 14 preenchimento, cobertura 98,9 %, camada 1
17,5 %, NDVI 0,295. sem_swir: 22 + 56, 99,6 %, camada 1 91,5 %, NDVI 0,232.
Poucas zero-cloud secas antes de 2008 → com_swir depende da camada 2.

## Em aberto

Etapa 3 ainda não exportada. Viés residual de ~−9 % (ganho OLS < 1
achata); RMA ou retas só com camada 1 são testes possíveis — não mudar
sem pedido. `maxent_fosfato.py` (`ASSET_MOSAICO`) ainda não aponta para estes assets.

Visualizador de decisão: `src/show_mosaics/visualizar_mosaico_unico.js`
(camada/ano/mês por pixel, ΔNDVI, viés % do SWIR estimado, prévia de 4
regras de combinação, estatísticas por célula, dispersão).

Achados de 2026-10-17 (600 m, |ΔNDVI| ≤ 0,05):
- os assets `com_swir` têm `meses = 7-10`; os `sem_swir` e os estimados,
  `1-12` (rodada com ano todo). Fração jul–out no sem_swir < 50 % em
  L3C5, L4C5, L5C4 e L5C5;
- com_swir: camada 1 ≈ 0 em 11 células; cobertura < 60 % nas colunas
  C3–C5;
- SWIR estimado sistematicamente abaixo do medido: viés de −15 a −29 % em
  B06/B08 na maioria das células; r(B08) 0,07–0,87. Causa não verificada.
