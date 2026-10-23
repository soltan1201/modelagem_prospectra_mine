# Spec — cnn-favorabilidade

## Objetivo
Manter o fluxo da CNN de favorabilidade de fosfato: patches do mosaico ASTER
v1 com pseudo-rótulos dos mapas RF/GTB → treino com validação cruzada por
célula (5 folds) → teste único → pedaços locais dos canais → mapa de p por
janela deslizante → coleção no GEE → relatório.

## Pré-condições
- Assets: `aster-mosaico-final`, `favorabilidade_ml` (RF e GTB, fosfato, v1),
  `pontos_referencia`, `grid_irece_5x5`.
- Conta EE que lê/escreve `projects/mapbiomas-arida/mine` (inicializa com
  `mapbiomas-brazil`); bucket `gs://cengine17-acervo-cprm` para a subida.
- Árvore `src/` completa no servidor/Drive (`model_DL/` importa
  `../modeling/maxent_fosfato.py`).
- Python: `earthengine-api`, `numpy`, `pandas`, `torch` (CUDA), `timm`,
  `scikit-learn`, `rasterio`, `google-cloud-storage`.

## Entradas
| Campo | Onde |
|---|---|
| quantidades e limiares do sorteio | topo de `exportar_patches.py` (`N_POS_CELULA`, `P_POS`, `P_NEG`, `P_DIF_MAX`, …) |
| modelo, backbone, patch, épocas, lr, folds, amp, clip, nome | argumentos de `treinar_cnn.py` (comandos em `rodadas_treino.md`) |
| modelo, n, borda, células, threads | argumentos de `baixar_tiles.py` |
| rodada, passo, células, tiles, sub, lote | argumentos de `inferir_cnn.py` |

## Saídas
- `patches_cnn/modelo_{A,B1}/patches_<LxCy>.npz`, `indice_patches.csv`,
  `divisao.json`.
- `runs_cnn/<nome>/` (`<modelo>_<backbone curto>_<patch>`, `_teste`):
  `config.json`, `fold<k>.pt`, `historico_fold<k>.csv`, `metricas_folds.csv`,
  `oof_validacao.csv`; com `--teste`: `final.pt`, `teste.csv`,
  `metricas_teste.json`.
- `tiles_cnn/canais_{A,B1}/<LxCy>/<LxCy>_r<i>c<j>.tif`, `indice_tiles.csv`.
- `<rodada>/mapas_passo<k>/p_<LxCy>.tif` (UInt16, p × 10000, 0 = sem dado) e
  `pedacos/`.
- `projects/mapbiomas-arida/mine/cnn_favorabilidade/<rodada>_p<k>_<LxCy>`;
  `<rodada>/registro_subida_gee_<data>.csv`.

## Efeitos colaterais
Muitas chamadas síncronas ao GEE (`computePixels`, com novas tentativas);
dezenas de GB no disco (pedaços); uploads ao bucket e tasks de ingestão.
Nenhuma task de export.

## Modos de falha
- `computePixels` > 48 MB por pedido: aumentar `--n` (pedaços) ou reduzir
  `BLOCO` (inferência pelo GEE). Canal em float64 dobra o pedido —
  `imagem_canais()` termina em `.toFloat()`.
- `Object of type int64 is not JSON serializable`: valores numpy na grade do
  `computePixels`; converter para `int` do Python.
- 503 do GEE: a biblioteca repete sozinha; só `ERRO:` no log conta — rodar de
  novo baixa o que falta.
- Treino divergindo (NaN): refazer a rodada com `--lr 3e-4 --amp 0`.
- Ocorrência de treino sem centro em magenta: fica fora (aviso no log).
- Fold sem ocorrência real: a parada usa o `ap` dos pseudo-rótulos.

## Não-objetivos
- Não roda nada; não decide quantidades, limiares nem o modelo final sem o
  usuário.
- Não refaz os mapas RF/GTB (skill `modelagem-favorabilidade`).
