---
name: acervo-cprm
description: Use when the user works with SGB/CPRM "Acervo SIG" data for the Irecê study area — listing plugin layers in QGIS, downloading selected packages, converting Geosoft .grd aerogeophysics (eU, U/Th, U/K), classifying references vs covariates, uploading to GCS and ingesting into GEE assets, or visualizing them. Triggers on "acervo CPRM", "Acervo SIG", "SGB", "1042", "geosoft grd", "exportar_acervo_gee", "bucket", "ocorrências de fosfato".
---

# Acervo SIG — CPRM/SGB (Irecê)

Documento central: [src/dowloadQGIS/ACERVO_CPRM_IRECE.md](../../../src/dowloadQGIS/ACERVO_CPRM_IRECE.md)
(seções 1–12: achados, temas, referências R1–R7 × covariáveis C1–C9,
decisões, export GEE, outros alvos, passo a passo gcloud). Atualize-o
quando mudar seleção, regra ou destino.

Área: retângulo -44.4605 / -14.7908 / -39.2035 / -9.4165.

## Scripts (`src/dowloadQGIS/`)

| Script | Faz |
|---|---|
| `listar_QGIS.py` | cataloga camadas do plugin (roda no QGIS ou fora); CSVs em `src/Dados/acervo_cprm/` |
| `inspecionar_pacotes.py` | lê zips remotos por HTTP Range (`ArquivoHTTP`) sem baixar |
| `baixar_acervo_irece.py` | baixa seleção (regras + `IDS_FIXOS`, filtro `RE_TEMA`) para `/mnt/Almacen/shpIRECE/`, `manifesto.csv` |
| `geosoft_grd_para_tif.py` | Geosoft .grd (zlib por bloco) → GeoTIFF; razões eU, U/Th, U/K |
| `exportar_acervo_gee.py` | prepara, envia ao bucket e ingere em `projects/mapbiomas-arida/mine/acervo_cprm` |
| `ver_acervo_cprm.js` | visualizador GEE |

Visualizador combinado: `src/show_mosaics/visualizar_mosaico_final_acervo.js`.

## Regras

- Rasters de 500 m ficam em 500 m; extensão = área baixada.
- Ocorrências: campos de substância por prefixo (`CAMPOS_SUBST`, inclui
  `RECURSO`, `DESCRI`); `CAMADAS_NAO_RECMIN` exclui camadas que não são
  ocorrência (Morro do Chapéu).
- Bucket `gs://cengine17-acervo-cprm`, projeto de cobrança
  `project-fadb9cba-032b-4fb9-882` (teste gratuito de solkan.cengine17).
  `cliente_storage` força `with_quota_project(PROJETO_STORAGE)` — sem isso,
  403 "billing account ... absent". Nunca sugerir "Fazer upgrade" (perde o
  crédito gratuito). Credencial ADC nunca vai para o git.
- Guia leigo: `GUIA_GCLOUD_BUCKET_GEE.md`.
- Usuário roda os scripts de envio ("sem rodar, eu mesmo rodo").

## Pontos de referência para modelos

`ocorrencias_fosfato` e `recursos_minerais` alimentam
`projects/mapbiomas-arida/mine/pontos_referencia` (fosfato + carbonato,
mvt, ouro, manganês) via `src/modeling/preparar_pontos_referencia.py` —
ver skill `modelagem-favorabilidade` e `docs/pontos_referencia.md`.

## Próximos alvos (só se pedido)

Seção 11 do md: calcário/dolomito primeiro, depois Pb–Zn–F–Ba;
geoquímica multielementar.
