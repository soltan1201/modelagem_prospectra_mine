---
name: cnn-favorabilidade
description: Use when the user works on the deep-learning (CNN) phosphate favorability model in src/model_DL — exporting patches from GEE with computePixels, pseudo-labels from the RF/GTB maps, train/validation/test split by grid cell, training runs (small CNN or timm backbones), downloading channel tiles, sliding-window inference, uploading the maps to GEE, the Colab notebook, or writing the deep-learning report. Triggers on "CNN", "deep learning", "patches", "exportar_patches", "treinar_cnn", "inferir_cnn", "baixar_tiles", "subir_predicao_gee", "rodadas", "folds", "cnn_colab", "pseudo-rótulo", "modelo A/B1", "backbone", "janela deslizante", "relatório DL".
---

# CNN de favorabilidade (fosfato)

Plano e decisões: [src/model_DL/plano_cnn.md](../../../src/model_DL/plano_cnn.md)
(fonte única — atualize-o quando uma decisão mudar). Comandos, nomes das
rodadas e justificativa dos folds:
[src/model_DL/rodadas_treino.md](../../../src/model_DL/rodadas_treino.md).
Contrato: [SPEC.md](SPEC.md). Origem dos pseudo-rótulos e métricas de
referência: skill `modelagem-favorabilidade`.

**O usuário roda tudo** (servidor `/srv/almacen/projects/modelagem_prospectra_mine`,
Arch Linux, conda `py_gee`/`py_DL`, GPU local ou Colab). Claude escreve e
confere a sintaxe; não roda treino, exportação, download nem inferência.

## Scripts (`src/model_DL/`)

| Script | Papel |
|---|---|
| `exportar_patches.py` | sorteio dos centros, divisão por célula, patches 64 × 64 px por `ee.data.computePixels` (sem task) → `.npz` por célula e por modelo; destino servidor ou Drive (rclone) |
| `exportar_patches_drive.py` | mesmo sorteio (importa `exportar_patches`), patches por `neighborhoodToArray` + `Export.table.toDrive` TFRecord (tasks, sem servidor/rclone) → pasta `patches_cnn_gee` no Drive da conta do EE |
| `tfrecord_para_npz.py` | no Colab: TFRecord do Drive → `.npz` por célula e modelo + `indice_patches.csv`/`divisao.json` (célula 6b) |
| `exportar_tiles_drive.py` | canais A/B1 por célula (com borda) por `Export.image.toDrive` → pasta `tiles_cnn_gee`; uma task por célula |
| `cortar_tiles_drive.py` | no Colab: corta os GeoTIFF por célula nos pedaços do `baixar_tiles` em `TILES_DIR` (célula 5 do `cnn_predict_colab.ipynb`) |
| `treinar_cnn.py` | treino com validação cruzada em 5 folds (`--nome` = pasta da rodada); `--teste` no fim |
| `baixar_tiles.py` | canais do modelo (A ou B1) da área toda em pedaços GeoTIFF float32 (6 × 6 por célula, borda 32 px); `--destino drive` sobe cada pedaço por rclone |
| `inferir_cnn.py` | janela deslizante; `--tiles padrao` lê os pedaços locais da rodada (A ou B1 pelo `config.json`) e grava p × 10000 UInt16 por célula; sem `--tiles`, baixa do GEE |
| `subir_predicao_gee.py` | bucket + ingestão na coleção `cnn_favorabilidade` |
| `cnn_colab.ipynb` | só treino no Colab, **importando** os `.py`: TFRecord → npz (`tfrecord_para_npz`), rodadas do `rodadas_treino.md` por menu, comparação, teste. Exportações rodam como `.py` na máquina local; predição no `cnn_predict_colab.ipynb` |
| `cnn_predict_colab.ipynb` | mapa no Colab: corta os tiles (`cortar_tiles_drive`), `inferir_cnn --tiles` (passo 3/1), mapa com os pontos de fosfato e percentil do p, `subir_predicao_gee` opcional |
| `copiar_tifs_colab.ipynb` | Colab, uma célula: copia `.tif` de pasta compartilhada por outra conta para pasta própria (liberar a cota da origem) |
| `rodadas_treino.md` | convenção `<modelo>_<backbone curto>_<patch>` (+ `_teste`), comandos, tmux, refazer rodada, VRAM |
| `guia_gdrive.md` | subir `patches_cnn/` e `runs_cnn/` do servidor para o Drive com rclone |

Os scripts compartilham funções: `baixar_tiles` e `inferir_cnn` usam
`exportar_patches.imagem_canais` (canais idênticos ao treino, todos float32) e
`treinar_cnn.criar_modelo`/`CANAL_BIN`; `exportar_patches` usa
`modeling/maxent_fosfato.INDICES` e `emissividade_tir`. Mudar uma fórmula num
lugar muda em todos.

## Regras do desenho (não mudar sem pedido)

- **Rótulos:** positivos só onde RF **e** GTB têm p ≥ 0,5 (magenta), inclusive
  os centros das ocorrências de treino; só RF / só GTB fora. Negativos: p ≤
  0,05 nos dois, > 2 km das ocorrências; negativos difíceis 0,05 < p ≤ 0,30.
  `y` = rótulo suave (média dos p; 1 nas ocorrências), `peso` por amostra.
- **Avaliação sem circularidade:** ocorrências de teste/validação sem filtro de
  p e fundo aleatório (só avaliação) — avaliar só onde RF/GTB já acertam infla
  o resultado.
- **Divisão por célula da grid sem solape**, nunca por patch: teste = células
  das ocorrências de validação; 5 folds nas demais. Justificativa dos folds (para
  o relatório): seção "Por que validação cruzada" do `rodadas_treino.md`.
- **Modelos com kernel 2D:** A = 15 índices + `swir_estimado` (16 canais);
  B1 = log da reflectância B01–B09 + emissividade B10–B14 + `swir_estimado`
  (15). Sem GLCM, sem bandas escolhidas por baixa correlação (perderia as
  diferenças entre bandas vizinhas, que são o sinal mineral).
- **Saída escalar por patch**; mapa por janela deslizante, sem decoder/U-Net.
  Passo ímpar (3 para comparar rodadas, 1 no mapa final).
- **Patch:** 64 px com backbones pré-treinados; 32 px só na CNN pequena.
- **Métricas, nunca acurácia:** decide `ap_ponto` (precisão média, ocorrências
  × fundo, fora do fold) com `ap_ganho`, `boyce`, `top10` e `auc_ponto` (GTB =
  0,84). Só os folds 0 e 1 têm ocorrências reais; nos outros a parada usa o
  `ap` dos pseudo-rótulos. `--teste` uma vez, com épocas fixas (média das
  melhores do CV).
- **Estabilidade:** AMP bf16 (fp16 com GradScaler se não houver bf16), clip
  de gradiente 1,0, lotes não finitos pulados; divergiu → refazer com
  `--lr 3e-4 --amp 0`.

## Vazamento conhecido

`favorabilidade_ml` (fonte dos pseudo-rótulos) foi treinado com todos os
pontos (`papeis = None`), inclusive os de teste. Teste independente exige
mapas RF/GTB com `PAPEIS = ['treino']` numa coleção separada e `ASSET_FAV`
apontando para ela. Citar no relatório como limitação.

## Relatório de deep learning (a escrever)

Mesmo padrão dos outros relatórios em `docs/` (LaTeX, figuras em `image/`,
dpi 350–450). Fontes: `plano_cnn.md` (desenho e decisões), `rodadas_treino.md`
(folds, rodadas, tabela comparativa preenchida), `metricas_folds.csv` /
`oof_validacao.csv` / `metricas_teste.json` de cada rodada, mapas da coleção
`cnn_favorabilidade` comparados a RF/GTB. Incluir a justificativa dos folds,
a escolha das métricas (por que não acurácia) e o vazamento.

## Caminhos

Servidor: patches `/srv/almacen/db_images/patches_cnn/modelo_{A,B1}/`,
rodadas `/srv/almacen/db_images/runs_cnn/<nome>/`, pedaços
`/srv/almacen/db_images/tiles_cnn/canais_{A,B1}/` (`PATCHES_DIR`, `RUNS_DIR`,
`TILES_DIR`). Pedaços: 900 por modelo (25 × 36), ~0,75 GB por célula
(≈ 18 GB por modelo). Colab: `MyDrive/prospectra/patches_cnn` e `.../runs_cnn`
(célula 2); servidor → Drive por `rclone copy`, ver `guia_gdrive.md`.
Sem servidor: patches e pedaços baixados em qualquer máquina direto para o
Drive (destino drive), treino e mapa no Colab (`USAR_TILES`, `TILES_DIR`) —
seção 0 do `guia_gdrive.md`.
