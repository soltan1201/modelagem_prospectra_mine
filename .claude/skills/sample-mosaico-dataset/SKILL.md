---
name: sample-mosaico-dataset
description: Use when the user wants to run or adjust the point-sampling script that harmonizes ASTER mosaic bands to 30m and exports random-point CSVs/assets for model training, change N_PONTOS/scale/seed, or debug ee.Initialize / export tasks. Triggers on "amostragem", "sample_quality_mosaic", "harmonizar bandas 30m", "exportar pontos GEE", "amostras para treino".
---

# Sample Mosaico Dataset

Etapa 3 do pipeline (ver
[src/building_dataset/readme_notas.md](../../../src/building_dataset/readme_notas.md)):
para cada mosaico semestral em `projects/mapbiomas-arida/mosaic_aster`,
harmoniza todas as bandas para 30 m e extrai pontos aleatórios exportados
como CSV (Drive) e/ou `FeatureCollection` (asset).

Ver [SPEC.md](SPEC.md) para o contrato formal.

## Arquivo de referência

[src/building_dataset/sample_quality_mosaic.py](../../../src/building_dataset/sample_quality_mosaic.py)
— **este script roda localmente** (não é Code Editor JS), via
`python3 src/building_dataset/sample_quality_mosaic.py`, desde que a conta
GEE (`mapbiomas-caatinga-cloud02`) esteja autenticada
(`earthengine authenticate` / ADC já configurado — confirmado disponível
neste ambiente: pacote `ee` instalado em
`~/.local/lib/python3.14/site-packages/ee`).

## Regras de harmonização (não altere sem motivo explícito)

- **VNIR (B01, B02, B3N) 15m → 30m**: `reduceResolution(mean, maxPixels=4)`
  — agrega exatamente 2×2 pixels.
- **SWIR (B04-B09) 30m**: sem alteração, é a resolução nativa/referência
  (`proj_30m = img.select('B04').projection()`).
- **TIR (B10-B14) 90m → 30m**: `resample('bilinear')`, não `reduceResolution`
  — é upsampling, não agregação.
- **`quality` 15m → 30m**: mesma regra do VNIR (`reduceResolution` com
  média) — importante para o filtro de qualidade em pós-processamento.
- A máscara de qualidade deve ser aplicada **antes** da harmonização (o
  comentário no código é explícito sobre isso) para que a média do VNIR já
  exclua pixels inválidos.

## Parâmetros ajustáveis (topo do arquivo)

`N_PONTOS` (30000), `ESCALA_AMOST` (30m), `TILE_SCALE` (4), `SEED` (42),
`INCLUIR_COORDS`, `EXPORTAR_PARA` (`'drive' | 'asset' | 'both'`). Mudar
`SEED` quebra a reprodutibilidade dos experimentos existentes — confirme
com o usuário antes.

## Fluxo de trabalho

1. Confirmar `EXPORTAR_PARA`: `'asset'` é necessário se o próximo passo for
   o skill `estimate-swir-gee` (que lê pontos de
   `projects/mapbiomas-arida/mine/points`); `'drive'` é necessário se o
   próximo passo for `train-swir-models` (que lê CSVs locais).
2. Rodar o script; ele itera por todas as imagens do asset
   `mosaic_aster` e dispara uma `Export.table` por imagem (CSV nomeado com
   `system:index`).
3. As tasks são assíncronas — o script apenas as inicia (`task.start()`) e
   imprime um resumo; monitoramento real é em
   https://code.earthengine.google.com/tasks (há um bloco comentado no
   final do arquivo para polling automático, desativado por padrão).
4. Se destino for Drive, os CSVs saem na pasta `ASTER_Samples/` — para
   treinar localmente (skill `train-swir-models`), baixe-os para
   `data/samples/` no repo.

## O que este skill NÃO faz

Não decide destino de export sozinho quando ambíguo — pergunte se é para
treino local (CSV/Drive), GEE-side (asset/pontos), ou ambos.
