---
name: cnn-favorabilidade
description: Use when the user works on the deep-learning (CNN) phosphate favorability model in src/model_DL — exporting patches from GEE with computePixels, pseudo-labels from the RF/GTB maps, train/validation/test split by grid cell, training (small CNN or timm backbones), sliding-window inference, or the Colab notebook. Triggers on "CNN", "deep learning", "patches", "exportar_patches", "treinar_cnn", "inferir_cnn", "cnn_colab", "pseudo-rótulo", "modelo A/B1", "backbone", "janela deslizante".
---

# CNN de favorabilidade (fosfato)

Plano e decisões: [src/model_DL/plano_cnn.md](../../../src/model_DL/plano_cnn.md)
(fonte única — atualize-o quando uma decisão mudar). Contrato: [SPEC.md](SPEC.md).
Origem dos pseudo-rótulos e métricas de referência: skill `modelagem-favorabilidade`.

**O usuário roda tudo** (servidor `/srv/almacen/projects/modelagem_prospectra_mine`,
GPU local ou Colab). Claude escreve e confere a sintaxe; não roda treino,
exportação nem inferência.

## Scripts (`src/model_DL/`)

| Script | Papel |
|---|---|
| `exportar_patches.py` | sorteio dos centros, divisão por célula, patches 64 × 64 px por `ee.data.computePixels` (sem task) → `.npz` por célula e por modelo |
| `treinar_cnn.py` | treino com validação cruzada por fold; `--teste` no fim |
| `inferir_cnn.py` | janela deslizante por célula → GeoTIFF de p |
| `cnn_colab.ipynb` | os três no Colab, **importando** os `.py` (não copia código) |
| `guia_gdrive.md` | subir `patches_cnn/` e `runs_cnn/` do servidor para o Drive com rclone |

Os três compartilham funções: `inferir_cnn` usa `exportar_patches.imagem_canais`
(canais idênticos ao treino) e `treinar_cnn.criar_modelo`/`CANAL_BIN`;
`exportar_patches` usa `modeling/maxent_fosfato.INDICES` e `emissividade_tir`.
Mudar uma fórmula num lugar muda nos três.

## Regras do desenho (não mudar sem pedido)

- **Rótulos:** positivos só onde RF **e** GTB têm p ≥ 0,5 (magenta), inclusive
  os centros das ocorrências de treino; só RF / só GTB fora. Negativos: p ≤
  0,05 nos dois, > 2 km das ocorrências. `y` = rótulo suave (média dos p; 1 nas
  ocorrências), `peso` por amostra.
- **Ocorrências de teste (papel = validacao) sem filtro de p** — avaliar só
  onde RF/GTB já acertam infla o resultado.
- **Divisão por célula da grid sem solape**, nunca por patch: teste = células
  das ocorrências de validação; 5 folds nas demais (células da bacia em folds
  diferentes); centros a < 1 km de célula de outro grupo descartados.
- **Modelos com kernel 2D:** A = 15 índices + `swir_estimado` (16 canais);
  B1 = log da reflectância B01–B09 + emissividade B10–B14 + `swir_estimado`
  (15). Sem GLCM, sem bandas escolhidas por baixa correlação (perderia as
  diferenças entre bandas vizinhas, que são o sinal mineral). `origem` não é
  canal (vai em `origem_patch`).
- **Saída escalar por patch** (classificação); mapa por janela deslizante,
  sem decoder/U-Net. Passo ímpar (1, 3, 5).
- **Patch:** 64 px com backbones pré-treinados (reduzem 32×); 32 px só na CNN
  pequena.
- **Métrica que decide:** `auc_ponto` (mediana dos jitters por ocorrência ×
  negativos, célula de fora) contra o GTB = 0,84 (mosaico v1). Só 2 folds têm
  ocorrências reais. `--teste` uma vez, com épocas fixas (média das melhores
  do CV) — o teste nunca entra na parada antecipada.

## Vazamento conhecido

`favorabilidade_ml` (fonte dos pseudo-rótulos) foi treinado com todos os
pontos (`papeis = None`), inclusive os de teste; o export avisa. Teste
independente exige mapas RF/GTB com `PAPEIS = ['treino']` numa coleção
separada e `ASSET_FAV` apontando para ela.

## Caminhos

Servidor: patches `/srv/almacen/db_images/patches_cnn/modelo_{A,B1}/`,
rodadas `/srv/almacen/db_images/runs_cnn/`; trocar com `PATCHES_DIR` e
`RUNS_DIR`. Colab: `MyDrive/prospectra/patches_cnn` e `.../runs_cnn` (célula
2), cópia local dos patches antes do treino (`COPIAR_LOCAL`). Servidor → Drive:
`rclone copy` (remote `gdrive:`), passo a passo em `guia_gdrive.md`; com os
patches no Drive, o notebook pula a exportação (células 6–7).

Tamanho: N patches por `.npz` = amostras da célula (~140 a ~400); total
~4.500–5.000 por modelo, ~0,6–0,7 GB em float16. Número exato: opção 1 do
`exportar_patches.py` (tabela origem × conjunto e `indice_patches.csv`).
