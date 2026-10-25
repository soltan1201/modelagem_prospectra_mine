# Plano — CNN para favorabilidade de fosfato (bloco de Irecê)

Decisões tomadas na discussão de 2026-10-22. Contexto e números:
`docs/relatorio_avaliacao_modelos/` (avaliação Maxent × RF × GTB, SHAP, mapas RF/GTB).

## Objetivo

Prever p(fosfato) por pixel a partir de um patch em volta do pixel, usando
contexto espacial que o RF/GTB (pixel a pixel) não veem. A rede só vale se
passar o GTB na mesma validação: AUC por ponto com uma célula da grid de fora
(GTB: 0,84 no mosaico v1, 0,63 no único).

## Rótulos

- **Positivos:** áreas com p ≥ 0,5 no RF **e** no GTB (`favorabilidade_ml`),
  sem manchas com borda de cena, rarefeitas; alvo y = média(p_RF, p_GTB)
  (rótulo suave), peso = 1.
- **Negativos:** áreas com p baixo nos dois modelos (fáceis ≤ 0,05 e difíceis
  em (0,05; 0,30]).
- **Ocorrências reais** (`pontos_referencia`, fosfato): as de teste/validação
  só avaliam, sem filtro de p; as de treino entram no treino só onde RF e GTB
  ≥ 0,5 (magenta), com y = 1.
- Cada patch guarda a célula da grid (validação célula de fora).

## Dois modelos, ambos com kernels 2D

Mesmo backbone, mesmos patches, mesma validação; muda só a entrada.

### Modelo A — índices espectrais (~15 canais + `swir_estimado`)

| Grupo | Índices |
|---|---|
| Ferro | `ferric_iron_2_1`, `ferric_ox_4_3`, `gossan_4_2`, `ferrous_sil_5_4` |
| Vegetação | `SAVI` |
| Al-OH | `sericite_illite_smectite`, `kaolinite_7_5`, `muscovite_7_6` |
| Carbonato | `carbonate_8_6`, `dolomite`, `carbonate_13_14` |
| Sílica / fosfato | `quartz_14_12`, `phosphorite_rbd` |
| Composições | `sultan_r`, `sultan_g` |

Fora por repetirem outra razão ou |r| > 0,9: `abrams_b`, `abrams_g`, `NDVI`,
`alunite_kaolinite_pyrophyllite`, `sultan_b`.
Normalização: corte nos percentis 1–99 (razões têm cauda longa) e padronização
por canal.

### Modelo B1 — bandas em log (14 canais + `swir_estimado`)

- B01, B02, B3N, B04–B09: log da reflectância TOA.
- B10–B14: emissividade normalizada (ε_max 0,96) em vez de DN — tira parte do
  efeito da temperatura da data.
- Com log, razão vira diferença (log B08 − log B06): a 1ª convolução 1×1
  pode reaprender os índices e achar combinações novas.
- Variação a testar: log centrado por pixel (cada banda − média do log do
  pixel), que remove o brilho geral e reduz a assinatura de cena.

Por que não "bandas pouco correlacionadas": com |r| ≤ 0,9 sobrariam ~4 bandas
(uma SWIR, uma TIR) e se perderiam as diferenças entre bandas vizinhas, que são
o sinal mineral.

### Fora dos dois

Texturas GLCM (a CNN aprende textura da B3N), NDVI, índices repetidos.

## Arquitetura

Backbone 2D pré-treinado (`timm`), 1ª camada adaptada ao número de canais
(`in_chans`). Ordem de teste: RegNetY-016 (atenção de canal SE embutida),
EfficientNetV2-S, ResNet-50/ResNet-RS-50, CoAtNet-0 (só se contexto ajudar).
Congelar as primeiras camadas, treinar a cabeça, descongelar aos poucos.

## Decisões de execução (2026-10-22)

- **Ambiente:** duas versões, Colab e GPU local, com o mesmo código:
  módulos `.py` em `src/model_DL/` e um notebook que só chama esses módulos,
  com `AMBIENTE = 'colab' | 'local'` definindo caminhos (Drive × disco).
- **Mosaico:** só o v1 (`aster-mosaico-final`), onde estão os pseudo-rótulos
  e onde o RF/GTB foram melhores.
- **Patch:** 64 px (~1,9 km) com backbones pré-treinados; 32 px (~1 km) só na
  CNN pequena treinada do zero. Motivo: os backbones reduzem a imagem 32×; com
  32 px a última camada vira 1×1 e os pesos do ImageNet (treinados a 224 px)
  ajudam pouco.
- **Saída:** um valor p ∈ [0, 1] por patch (classificação). Sem decoder e sem
  U-Net.
- **Mapa:** janela deslizante — o patch centrado em cada pixel dá o p daquele
  pixel. Estimativa (RegNetY-016 a 64 px, ~0,26 GFLOP por patch; área ≈ 376
  milhões de pixels de 30 m): passo de 1 pixel ≈ 10¹⁷ FLOP (~1–3 h numa GPU
  com fp16); passo de 3 pixels (90 m) ≈ 9× menos. Começar com passo 3 por
  célula da grid; passo 1 na versão final.

## Ainda em aberto

- Pesos pré-treinados em satélite (SSL4EO, ResNet-50) além do ImageNet.

## Exportação dos patches (desenho, 2026-10-22)

Sem tasks: `ee.data.computePixels` (ou `getDownloadURL` em formato NPY) devolve
o recorte direto como array NumPy; os patches são gravados em `.npz` locais.

### Amostras e rótulos

Como está em `exportar_patches.py` (atualizado em 2026-10-23; o desenho
inicial dava peso 0,5 a só RF / só GTB, que acabaram fora do treino).

| Origem | Como é sorteada | `y` (alvo) | `peso` | `uso_treino` |
|---|---|---|---|---|
| `ocorrencia_treino` | ponto de fosfato com `papel` = treino + 8 centros deslocados até 500 m | 1 | 1,0 | 1 só se em magenta (RF e GTB ≥ 0,5); senão 0 |
| `ocorrencia_teste` | ponto com `papel` = validação + 8 centros, sem filtro de p | 1 | 1,0 | 0 (só avaliação) |
| `ambos` (magenta) | RF e GTB com p ≥ 0,5; 120 por célula, rarefeitos a 1 km | média(p_RF, p_GTB) | 1,0 | 1 |
| `negativo` | RF e GTB com p ≤ 0,05, > 2 km das ocorrências; 150 por célula | média(p_RF, p_GTB) ≈ 0 | 1,0 | 1 |
| `negativo_dificil` | RF e GTB com p em (0,05; 0,30], longe das ocorrências; 100 por célula | média(p_RF, p_GTB) | 1,0 | 1 |
| `fundo` | pixels aleatórios sem filtro de p; 80 por célula | média(p_RF, p_GTB) | 1,0 | 0 (só avaliação) |

Negativos (fáceis, difíceis) e fundo rarefeitos a 1,5 km. Fora do treino: só
RF / só GTB e a faixa 0,30 < p < 0,5 (zona incerta).

A incerteza do rótulo vai só no `y` contínuo (rótulo suave); `peso` = 1 em
todas as amostras (o campo continua nos `.npz` para pesos por origem no
futuro). O treino usa BCE com alvo suave, ponderada por `peso` × fator de
classe: `treinar_cnn.fatores_classe` iguala o peso total de positivos
(y ≥ 0,5) e negativos no treino (a₊ = 0,5/π₊, a₋ = 0,5/(1 − π₊)).

### Canais por patch (64 × 64 px, 30 m, grade do mosaico v1)

Um arquivo serve aos dois modelos: 15 índices (modelo A) + 14 bandas em log /
emissividade (modelo B1) + `swir_estimado` + `origem` (para conferir borda de
cena). Gravados em float16.

### Arquivos

Um `.npz` por célula da grid: `X (N, C, 64, 64)`, `y`, `peso`, `origem_amostra`,
`celula`, `lon`, `lat`, `p_rf`, `p_gtb`, `id`; mais `indice_patches.csv` com
todos os metadados e a lista de canais. Retomável: pula ids já gravados.

## Divisão treino / validação / teste (confirmada 2026-10-22)

Por célula da grid (sem solape), nunca por patch: `divisao.json`.

| Conjunto | Células | Uso |
|---|---|---|
| teste | as que têm ocorrências de validação (leste/norte) | uma vez, no fim; ocorrências reais sem filtro de p |
| fold 0–4 | as demais, em 5 grupos (as da bacia em folds diferentes) | validação cruzada por grupo: parada, hiperparâmetros, A × B1 |

Positivos de treino só em magenta (RF e GTB ≥ 0,5), inclusive os centros das
ocorrências de treino; só RF / só GTB ficam fora. Centros a < 1 km de célula de
outro grupo são descartados. Script: `exportar_patches.py`.

Vazamento conhecido: `favorabilidade_ml` foi treinado com todos os pontos;
os mapas não serão refeitos agora — o teste fica marcado como não independente
até existir uma versão com `PAPEIS = ['treino']`.
