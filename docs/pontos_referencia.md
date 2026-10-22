# Pontos de referência disponíveis — Irecê

Levantado em 2026-10-19 direto dos assets no GEE (`projects/mapbiomas-arida/mine/`)
e dos CSVs em `src/Dados/maxent/`. "Na área" = dentro do retângulo da área de
estudo (`na_area = 1`); distâncias a partir de Irecê (-41,86, -11,30).
Códigos R1–R7 seguem `src/dowloadQGIS/ACERVO_CPRM_IRECE.md`, seção 7.

## 0. Asset unificado `pontos_referencia`

`src/modeling/preparar_pontos_referencia.py` junta tudo num asset só
(`projects/mapbiomas-arida/mine/pontos_referencia`, e CSV em
`src/Dados/referencias/pontos_referencia.csv`), com `alvo`, `fonte`, `papel`,
`subst` padronizada e `CODIGO_OBJ` único. Duplicatas do mesmo alvo a ≤ 500 m
viram um ponto só.

| alvo | pontos | origem |
|---|---|---|
| fosfato | 47 (29 `antigo+cprm`, 18 `cprm`; 13 treino na bacia, 34 validação fora) | `pontos_fosfato` + `ocorrencias_fosfato` |
| carbonato | 327 | calcário, dolomito, mármore, magnesita, calcita, Ca |
| mvt | 252 | Ba, barita, F, fluorita, Pb, Zn |
| ouro | 208 | Au, ouro |
| manganes | 269 | Mn, manganês |

Os 30 pontos antigos estão todos na base CPRM (a ≤ 66 m).
`maxent_fosfato.py` roda um Maxent por alvo e exporta `probability` e
`complemento` (1 − p) como covariáveis.

## 1. Fosfato (alvo do modelo)

| Asset / arquivo | Pontos | Na área | Conteúdo | Papel |
|---|---|---|---|---|
| `acervo_cprm/ocorrencias_fosfato` | 50 | 39 | recmin/ocormin com P, sem duplicatas; campo `grupo` | R1 + R2 |
| ↳ `grupo = bacia_irece` | 14 | 14 | 1 mina (Irecê/Faz. Juazeiro, 5 km), 12 não explotados a 7–13 km, Junco (63 km), Lagedinho (33 km) | **R1 presença (treino)** |
| ↳ `grupo = leste_outros` | 26 | 21 | ~200–245 km a leste | **R2 validação independente** |
| ↳ `grupo = norte` | 10 | 4 | Sobradinho, ~225–247 km | **R2 validação independente** |
| `pontos_fosfato` | 30 | 29 | conjunto antigo (23 não explotados, 7 garimpos); é o que `maxent_fosfato.py` e `baixar_presenca_maxent.py` ainda leem | presença (versão anterior) |
| `src/Dados/maxent/maxent_fosfato_presenca_mosaico_prioridade.csv` | 320 pixels | — | pixels do mosaico v1 em volta de 16 ocorrências (`CODIGO_OBJ`), origem 2–3 | presença extraída |

Status em `ocorrencias_fosfato`: 29 não explotados, 13 garimpos, 4 minas, 4 sem status.

## 2. Geoquímica de rocha (P₂O₅)

`acervo_cprm/geoquimica_p2o5`: 380 amostras, 334 na área (Calcários da Bahia
2021: 333; Agrominerais Irecê–Jaguarari 2016: 47).

| Classe | Na área | Até 100 km | Papel |
|---|---|---|---|
| P₂O₅ ≥ 5 % | 5 | 4 | **R3** presença complementar ou validação (falta decidir) |
| 1 ≤ P₂O₅ < 5 % | 2 | — | anômalo, sem papel definido |
| P₂O₅ < 1 % | 327 | 90 | **R4** ausência / fundo informado |

## 3. Outras ocorrências minerais

`acervo_cprm/recursos_minerais`: 7.992 pontos, 4.180 na área, sem duplicatas.
Papel principal: **R5**, fundo com o mesmo viés de amostragem do fosfato.
Também são referências para outros alvos (ver `ACERVO_CPRM_IRECE.md` §11):

| Grupo | Na área | Observação |
|---|---|---|
| Calcário / dolomito / mármore / magnesita / calcita | 320 (75 até 150 km) | Próximo alvo sugerido; mesmo hospedeiro do fosfato |
| Ba / F / Pb / Zn (estilo MVT) | 337 | Contando grafias diferentes (Ba, Bário…) |
| Ouro (Au + Ouro) | 333 | Jacobina / Chapada Diamantina |
| Manganês (Mn + Manganês) | 303 | — |
| Talco / vermiculita | 122 / 74 | Longe de Irecê |
| Diamante, quartzo/cristal | 254 / 534 | Não recomendados com ASTER |

Status (na área): 1.656 não explotados, ~1.350 garimpos, ~390 minas, 676 sem status.
O campo `subst` tem grafias diferentes para a mesma substância; padronizar antes de usar.

## 4. Fundo / amostras do Maxent

| Asset / arquivo | Pontos | Conteúdo | Papel |
|---|---|---|---|
| `maxent_fosfato_amostras_teste_zero_cloud_2000_2008` | 4.925 | todos `classe = 0`, com bandas do zero-cloud 2000–2008 por célula | fundo (teste) |
| `src/Dados/maxent/celulas/…/celula_<id>.csv` | — | o mesmo, separado por célula (17 células) | fundo (teste) |

## 5. Pontos do modelo de SWIR (não são referência geológica)

| Asset | Pontos | Conteúdo |
|---|---|---|
| `points_zeroCC_00_08` | 34.000 | bandas VNIR+SWIR+TIR, cenas zero-cloud 2000–2008 |
| `points_zeroCC_09_20` | 30.631 | VNIR+TIR, cenas zero-cloud 2009–2020 |

Servem para treinar o RF que estima o SWIR; não entram no Maxent.

## 6. Referências que não são pontos

- **R6 afloramentos descritos** (~19 mil na área): só em shapefile local
  (`/mnt/Almacen/shpIRECE/*/afloramento.shp`), **não estão no GEE**.
- **R7 unidades mapeadas:** `acervo_cprm/lito_jacobina_250k` (653 polígonos) e
  `acervo_cprm/lito_bahia` (4.516); filtrar `NOME_UNIDA` (Salitre, Irecê, Una…).
- Apoio: `acervo_cprm/dolinas_morro_chapeu` (89 polígonos),
  `estrut_jacobina_250k` (1.979 linhas), `estrut_bahia` (3.603 linhas).

## Pendências

- `baixar_presenca_maxent.py` ainda lê `pontos_fosfato` (30).
- Decidir se as 5 amostras R3 (P₂O₅ ≥ 5 %) entram no treino ou só na validação.
- Subir os afloramentos (R6) para o GEE se forem usados para validar os índices ASTER.
