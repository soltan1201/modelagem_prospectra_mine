# Spec — cnn-favorabilidade

## Objetivo
Manter o fluxo da CNN de favorabilidade de fosfato: patches do mosaico ASTER
v1 com pseudo-rótulos dos mapas RF/GTB → treino com validação por célula →
teste único → mapa de p por janela deslizante.

## Pré-condições
- Assets: `aster-mosaico-final`, `favorabilidade_ml` (RF e GTB, fosfato, v1),
  `pontos_referencia`, `grid_irece_5x5`.
- Conta EE que lê `projects/mapbiomas-arida/mine` (inicializa com
  `mapbiomas-brazil`).
- Árvore `src/` completa no servidor/Drive (`model_DL/` importa
  `../modeling/maxent_fosfato.py`).
- Python: `earthengine-api`, `numpy`, `pandas`, `torch` (CUDA), `timm`,
  `scikit-learn`; `rasterio` ou GDAL para os GeoTIFFs.

## Entradas
| Campo | Onde |
|---|---|
| quantidades e limiares do sorteio | topo de `exportar_patches.py` (`N_POS_CELULA`, `P_POS`, `P_NEG`, …) |
| modelo, backbone, patch, épocas, lr, folds | argumentos de `treinar_cnn.py` / formulário do notebook |
| rodada, passo, células | argumentos de `inferir_cnn.py` |

## Saídas
- `patches_cnn/modelo_{A,B1}/patches_<LxCy>.npz`, `indice_patches.csv`,
  `divisao.json`.
- `runs_cnn/<modelo>_<backbone>_<patch>px_<data>/`: `config.json`,
  `fold<k>.pt`, `historico_fold<k>.csv`, `metricas_folds.csv`,
  `oof_validacao.csv`; com `--teste`: `final.pt`, `teste.csv`,
  `metricas_teste.json`.
- `<rodada>/mapas_passo<k>/p_<LxCy>.tif` + `config_inferencia.json`.
- Cópia no Drive (`gdrive:prospectra/patches_cnn`, `.../runs_cnn`) via rclone
  (`guia_gdrive.md`), feita pelo usuário.

## Efeitos colaterais
Muitas chamadas síncronas ao GEE (`computePixels`, 8 threads, com novas
tentativas); arquivos grandes no disco. Nenhuma task de export.

## Modos de falha
- `computePixels` > 48 MB por pedido: reduzir `BLOCO` na inferência.
- Formato do retorno do `computePixels` diferente (versão do
  `earthengine-api`): ajustar a leitura em `baixar_patch` / `baixar_bloco`.
- Ocorrência de treino sem centro em magenta: fica fora (aviso no log).
- Fold sem ocorrência real: a parada usa `auc_amostra` (concordância com
  RF/GTB), não `auc_ponto`.

## Não-objetivos
- Não roda nada; não decide quantidades, limiares nem o modelo final sem o
  usuário.
- Não refaz os mapas RF/GTB (skill `modelagem-favorabilidade`).
