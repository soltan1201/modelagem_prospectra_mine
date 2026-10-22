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

### Refazer uma rodada (travou, divergiu ou mudou algo)

O `--nome` não sobrescreve: se a pasta já tiver `.pt`, o treino para. Para
refazer, apague a pasta da rodada antes (apaga os modelos e históricos dela —
não dá para desfazer):

```bash
rm -r /srv/almacen/db_images/runs_cnn/A_regnety_64
python treinar_cnn.py --modelo A --backbone regnety_016 --patch 64 --nome A_regnety_64
```

Troque `A_regnety_64` pelo nome da rodada. Alternativa sem apagar: rodar com
outro nome (ex.: `--nome A_regnety_64_v2`).

Se o modelo divergir (predições NaN, aviso "modelo divergiu"), refaça com lr
menor e float32 puro:

```bash
rm -r /srv/almacen/db_images/runs_cnn/A_regnety_64
python treinar_cnn.py --modelo A --backbone regnety_016 --patch 64 --nome A_regnety_64 --lr 3e-4 --amp 0
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

## 5. Mapa (pedaços locais → predição → GEE)

```bash
df -h /srv/almacen                       # espaço livre (≈ 10–24 GB por modelo nos pedaços)
du -sh /srv/almacen/db_images/*           # quanto cada pasta ocupa

# 1) canais em pedaços 6×6 por célula — uma rodada por modelo (só o escolhido)
python baixar_tiles.py --modelo A          # → tiles_cnn/canais_A/<LxCy>/<LxCy>_r<i>c<j>.tif
python baixar_tiles.py --modelo B1         # → tiles_cnn/canais_B1/...

# 2) predição pedaço a pedaço (sem GEE); --sub 2 se faltar memória, --lote menor na GPU
python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --tiles padrao --passo 3 --celulas L2C3,L3C3
python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --tiles padrao --passo 1
#    → <rodada>/mapas_passo1/pedacos/<LxCy>/p_<LxCy>_r<i>c<j>.tif e <rodada>/mapas_passo1/p_<LxCy>.tif

# 3) subir ao GEE (bucket + ingestão) → projects/mapbiomas-arida/mine/cnn_favorabilidade
python subir_predicao_gee.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --passo 1
```

Mapas em UInt16 (p × 10000, 0 = sem dado), como o `favorabilidade_ml`. A
rodada deve ter sido treinada com o mesmo modelo de canais dos pedaços (A ou B1).
Sem `--tiles`, o `inferir_cnn.py` ainda baixa do GEE na hora (bom para testar
poucas células).

### Baixar os pedaços (`baixar_tiles.py`)

Baixe só o modelo de canais da rodada escolhida (A **ou** B1). Rode num `tmux`
(são 25 células × 36 pedaços = 900 pedidos ao GEE) e guarde o log:

```bash
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
df -h /srv/almacen                                    # espaço livre antes
mkdir -p /srv/almacen/db_images/tiles_cnn

# teste rápido: 1 célula (36 pedaços) — confere tempo, tamanho e % com dado
python baixar_tiles.py --modelo A --celulas L2C3

# área toda (pula os pedaços já gravados)
python baixar_tiles.py --modelo A 2>&1 | tee /srv/almacen/db_images/tiles_cnn/log_baixar_A.txt
python baixar_tiles.py --modelo B1 2>&1 | tee /srv/almacen/db_images/tiles_cnn/log_baixar_B1.txt
```

Com `tmux` (continua rodando se a conexão SSH cair):

```bash
tmux new -s tiles                         # abre a sessão "tiles"

# dentro do tmux:
conda activate py_gee
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
python baixar_tiles.py --modelo A 2>&1 | tee /srv/almacen/db_images/tiles_cnn/log_baixar_A.txt

# sair sem parar o processo: Ctrl+b, solta, d
```

```bash
tmux ls                                   # sessões abertas
tmux attach -t tiles                      # voltar à sessão
tail -f /srv/almacen/db_images/tiles_cnn/log_baixar_A.txt   # acompanhar sem entrar (Ctrl+c sai do tail)
tmux kill-session -t tiles                # fechar a sessão quando terminar (para o que estiver rodando)
```

Os dois modelos em sequência numa só sessão, sem precisar entrar:

```bash
tmux new -d -s tiles "bash -lc 'conda activate py_gee && \
  cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL && \
  python baixar_tiles.py --modelo A 2>&1 | tee /srv/almacen/db_images/tiles_cnn/log_baixar_A.txt && \
  python baixar_tiles.py --modelo B1 2>&1 | tee /srv/almacen/db_images/tiles_cnn/log_baixar_B1.txt; exec bash'"
```

O `exec bash` no fim deixa a sessão aberta depois de terminar, para ver a
saída com `tmux attach -t tiles`. Se o `conda activate` falhar dentro do
`bash -lc`, troque por `source ~/miniconda3/bin/activate py_gee` (ajuste o
caminho do conda).

Opções:

| Opção | Padrão | Para quê |
|---|---|---|
| `--modelo A\|B1` | — (obrigatório) | canais do modelo (16 ou 15) |
| `--celulas L2C3,L3C3` | todas | só algumas células |
| `--n 6` | 6 | pedaços por lado em cada célula; aumente (ex.: `--n 8`) se aparecer "Pedaço passa de 48 MB" |
| `--borda 32` | 32 (meio patch de 64) | px extras de cada lado; com patch 32 basta `--borda 16` |
| `--threads 4` | 4 | pedidos simultâneos ao GEE; reduza para 2 se vierem muitos erros de cota |
| `--dir /outro/lugar` | `/srv/almacen/db_images/tiles_cnn` | pasta de saída (ou `export TILES_DIR=...`) |

Conferir o que já baixou e o espaço usado:

```bash
ls /srv/almacen/db_images/tiles_cnn/canais_A/*/*.tif | wc -l     # 900 = completo (25 × 36)
du -sh /srv/almacen/db_images/tiles_cnn/canais_A
grep ERRO /srv/almacen/db_images/tiles_cnn/log_baixar_A.txt      # pedaços que falharam
```

Pedaço com erro não para o resto: rode o mesmo comando de novo e ele baixa só
os que faltam. Para refazer uma célula inteira (apaga os pedaços dela — não dá
para desfazer):

```bash
rm -r /srv/almacen/db_images/tiles_cnn/canais_A/L2C3
python baixar_tiles.py --modelo A --celulas L2C3
```

Se mudar `--n` ou `--borda`, apague a pasta `canais_<modelo>` inteira antes:
os nomes dos pedaços (`r<i>c<j>`) mudam de significado e misturar dá mapa errado.
