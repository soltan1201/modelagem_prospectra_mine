# Rodadas de treino da CNN — comandos e nomes das saídas

Todas as combinações planejadas, com um **nome fixo por rodada** (`--nome`).
O nome vira a pasta da rodada em `RUNS_DIR` e é o que se passa ao
`inferir_cnn.py` depois. Se a pasta já tiver `.pt`, o treino para (não
sobrescreve): use outro nome ou apague a pasta.

- Patches: `/srv/almacen/db_images/patches_cnn/modelo_{A,B1}/` (`PATCHES_DIR`)
- Rodadas: `/srv/almacen/db_images/runs_cnn/<nome>/` (`RUNS_DIR`)
- Modelo **A** = 15 índices espectrais + `swir_estimado` (16 canais)
- Modelo **B1** = log da reflectância B01–B09 + emissividade B10–B14 + `swir_estimado` (15 canais)
- Rodar de dentro de `src/model_DL/`, de preferência num `tmux`.

```bash
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
```

## Convenção de nomes

`<modelo>_<backbone curto>_<patch>` — e `_teste` para a rodada final.

| Backbone (`--backbone`) | Nome curto |
|---|---|
| `pequena` (CNN do zero, conv + SE) | `pequena` |
| `regnety_016` | `regnety` |
| `tf_efficientnetv2_s` | `effv2s` |
| `resnet50` | `resnet50` |
| `coatnet_0_rw_224` | `coatnet0` |

## Arquivos de cada rodada (`runs_cnn/<nome>/`)

| Arquivo | O que é |
|---|---|
| `fold0.pt` … `fold4.pt` | modelo de cada fold (pesos + normalização + canais + argumentos) — usados juntos (ensemble) pelo `inferir_cnn.py` |
| `historico_fold<k>.csv` | métricas por época de cada fold (época da melhor `ap_ponto` → escolher `--epocas` do teste) |
| `metricas_folds.csv` | melhor época de cada fold |
| `oof_validacao.csv` | predição de cada patch fora do seu fold |
| `config.json` | argumentos, canais, número de patches, dispositivo |
| `final.pt` | só nas rodadas `_teste`: modelo treinado em todos os folds |
| `teste.csv`, `metricas_teste.json` | só nas rodadas `_teste`: predições e métricas no conjunto de teste |
| `mapas_passo<k>/p_<LxCy>.tif` | criados depois pelo `inferir_cnn.py` |

## 1. Linha de base — CNN pequena, 32 px

| Nome | Comando |
|---|---|
| `A_pequena_32` | `python treinar_cnn.py --modelo A --backbone pequena --patch 32 --nome A_pequena_32` |
| `B1_pequena_32` | `python treinar_cnn.py --modelo B1 --backbone pequena --patch 32 --nome B1_pequena_32` |

## 2. Backbones pré-treinados (ImageNet), 64 px

| Nome | Comando |
|---|---|
| `A_regnety_64` | `python treinar_cnn.py --modelo A --backbone regnety_016 --patch 64 --nome A_regnety_64` |
| `B1_regnety_64` | `python treinar_cnn.py --modelo B1 --backbone regnety_016 --patch 64 --nome B1_regnety_64` |
| `A_effv2s_64` | `python treinar_cnn.py --modelo A --backbone tf_efficientnetv2_s --patch 64 --nome A_effv2s_64` |
| `B1_effv2s_64` | `python treinar_cnn.py --modelo B1 --backbone tf_efficientnetv2_s --patch 64 --nome B1_effv2s_64` |
| `A_resnet50_64` | `python treinar_cnn.py --modelo A --backbone resnet50 --patch 64 --nome A_resnet50_64` |
| `B1_resnet50_64` | `python treinar_cnn.py --modelo B1 --backbone resnet50 --patch 64 --nome B1_resnet50_64` |
| `A_coatnet0_64` | `python treinar_cnn.py --modelo A --backbone coatnet_0_rw_224 --patch 64 --nome A_coatnet0_64` |
| `B1_coatnet0_64` | `python treinar_cnn.py --modelo B1 --backbone coatnet_0_rw_224 --patch 64 --nome B1_coatnet0_64` |

O CoAtNet é o que mais pede dado: rode só se os outros mostrarem que o contexto
ajuda. Se faltar memória na GPU, acrescente `--lote 32`.

Para rodar as 10 em sequência (num `tmux`):

```bash
mkdir -p /srv/almacen/db_images/runs_cnn
for combo in "A pequena 32 pequena" "B1 pequena 32 pequena" \
             "A regnety_016 64 regnety" "B1 regnety_016 64 regnety" \
             "A tf_efficientnetv2_s 64 effv2s" "B1 tf_efficientnetv2_s 64 effv2s" \
             "A resnet50 64 resnet50" "B1 resnet50 64 resnet50" \
             "A coatnet_0_rw_224 64 coatnet0" "B1 coatnet_0_rw_224 64 coatnet0"; do
  set -- $combo
  python treinar_cnn.py --modelo $1 --backbone $2 --patch $3 --nome ${1}_${4}_${3} \
      2>&1 | tee /srv/almacen/db_images/runs_cnn/log_${1}_${4}_${3}.txt
done
```

## 3. Comparar as rodadas

Olhar a linha **FORA DO FOLD (todos)** de cada rodada (ou a célula 9 do
`cnn_colab.ipynb`): `ap_ponto`, `ap_ganho`, `boyce`, `top10` e `auc_ponto`
(referência GTB = 0,84). Escolher o melhor modelo (A ou B1) e backbone.

| Nome | ap_ponto | ap_ganho | boyce | top10 | auc_ponto | melhor época média |
|---|---|---|---|---|---|---|
| `A_pequena_32` | | | | | | |
| `B1_pequena_32` | | | | | | |
| `A_regnety_64` | | | | | | |
| `B1_regnety_64` | | | | | | |
| `A_effv2s_64` | | | | | | |
| `B1_effv2s_64` | | | | | | |
| `A_resnet50_64` | | | | | | |
| `B1_resnet50_64` | | | | | | |
| `A_coatnet0_64` | | | | | | |
| `B1_coatnet0_64` | | | | | | |

## 4. Teste final (uma vez, só com o escolhido)

Épocas fixas = média das melhores épocas do CV (`historico_fold*.csv` /
célula 9). Exemplo com o RegNetY no modelo A:

| Nome | Comando |
|---|---|
| `A_regnety_64_teste` | `python treinar_cnn.py --modelo A --backbone regnety_016 --patch 64 --teste --epocas 20 --nome A_regnety_64_teste` |

Troque modelo/backbone/épocas pelo escolhido, mantendo o sufixo `_teste`.
Atenção: enquanto os mapas RF/GTB usados nos pseudo-rótulos forem os treinados
com todos os pontos, o teste não é independente (ver `plano_cnn.md`).

## 5. Mapa

```bash
# primeiro: passo 3 (90 m) em duas células
python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --passo 3 --celulas L2C3,L3C3
# depois: todas as células, passo 1 (30 m)
python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --passo 1
```

Com uma rodada sem `_teste` (sem `final.pt`), a inferência usa a média dos
5 `fold*.pt`. Os mapas ficam em `<rodada>/mapas_passo<k>/p_<LxCy>.tif`.
