# Rodadas de treino da CNN — comandos e nomes das saídas

Todas as combinações planejadas, com um **nome fixo por rodada** (`--nome`).
O nome vira a pasta da rodada em `RUNS_DIR` e é o que se passa ao
`inferir_cnn.py` depois. Se a pasta já tiver `.pt`, o treino para (não
sobrescreve): use outro nome ou apague a pasta.

## Por que validação cruzada em 5 folds (e não uma divisão treino/validação)

Os folds não são do PyTorch: é validação cruzada, uma escolha do desenho (no
TensorFlow/Keras seria o mesmo laço `for k in range(5): fit(treino ≠ k,
validation_data = fold k)`). Os tutoriais usam uma divisão única
(`validation_split`) porque têm muitos dados e amostras independentes; aqui
nenhuma das duas condições vale.

1. **Poucas ocorrências reais.** Só 2 dos 5 folds têm ocorrências de fosfato.
   Numa divisão única a validação poderia cair em células sem nenhuma, e a
   métrica mediria só a concordância com o RF/GTB, não a capacidade de achar
   fosfato.
2. **Autocorrelação espacial.** Patches vizinhos são quase iguais: dividir por
   patch ao acaso põe vizinhos no treino e na validação e infla a métrica. A
   divisão é por célula da grid 5 × 5, sem solape (centros a < 1 km de célula de
   outro grupo descartados); com tão poucas células, uma divisão só dependeria
   muito de quais células caíram na validação.
3. **Estimativa mais estável.** Cada célula fora do teste é validada uma vez,
   por um modelo que nunca a viu; o `oof_validacao.csv` junta essas predições e
   dá a linha **FORA DO FOLD (todos)** usada para comparar backbones e A × B1.
4. **Épocas do teste.** A média das melhores épocas dos folds fixa o
   `--epocas` da rodada `_teste`; o conjunto de teste (células das ocorrências
   de validação) nunca entra na parada antecipada.
5. **Ensemble no mapa.** Sem `final.pt`, o `inferir_cnn.py` faz a média dos 5
   `fold<k>.pt`, o que deixa o mapa menos ruidoso.

Custo: cada rodada treina 5 modelos (≈ 5× o tempo de uma divisão única). Para
triagem rápida dá para treinar só os folds com ocorrências (`--folds 0,1`); a
escolha final usa os 5.

## Caminhos e modelos

- Patches: `/srv/almacen/db_images/patches_cnn/modelo_{A,B1}/` (`PATCHES_DIR`)
- Rodadas: `/srv/almacen/db_images/runs_cnn/<nome>/` (`RUNS_DIR`)
- Modelo **A** = 15 índices espectrais + `swir_estimado` (16 canais)
- Modelo **B1** = log da reflectância B01–B09 + emissividade B10–B14 + `swir_estimado` (15 canais)
- Rodar de dentro de `src/model_DL/`, de preferência num `tmux`.
- Dois servidores (.250 com 8 GB e .251 com 12 GB): cópia dos dados e divisão
  das rodadas em [`guia_dois_servidores.md`](guia_dois_servidores.md).

```bash
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
```

### Patches pelo Google Drive (sem servidor)

Alternativa ao `exportar_patches.py` no servidor: a base de treino/validação/teste
sai do GEE direto para o Drive, com o mesmo sorteio e a mesma divisão por célula
(importados do `exportar_patches.py`). Fluxo completo em
[`guia_gdrive.md`](guia_gdrive.md) (seção 0).

```bash
# máquina local: tasks do GEE (Export.table.toDrive) → Drive patches_cnn_gee/ (TFRecord)
python exportar_patches_drive.py --so-sortear                          # só mostra as quantidades
python exportar_patches_drive.py --projeto ee-solkancengine17
python exportar_patches_drive.py --projeto ee-solkancengine17 --celulas L2C3,L3C3   # teste
```

Retomável: tasks com a mesma descrição já pendentes, rodando ou concluídas são
puladas. Depois, no Colab (célula 6b do `cnn_colab.ipynb`), converter para o
formato do treino (`modelo_A/` e `modelo_B1/` `patches_<LxCy>.npz`,
`indice_patches.csv`, `divisao.json`, ~1,3 GB):

```bash
!python tfrecord_para_npz.py --origem /content/drive/MyDrive/patches_cnn_gee \
    --destino /content/drive/MyDrive/prospectra/patches_cnn
```

Confere o número de patches de cada célula contra o `indice_patches.csv`.
Para treinar no servidor, trazer com
`rclone copy gdrive:prospectra/patches_cnn /srv/almacen/db_images/patches_cnn --progress`.
Com o mesmo `SEED` e os mesmos assets, o sorteio sai igual; se algum asset
mudou, as rodadas antigas deixam de ser comparáveis.

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
| `mapas_passo<k>_preenchido/` | bordas preenchidas pelo `preencher_bordas.py --preencher` (+ `relatorio_bordas.csv`) |
| `registro_subida_gee_p<k>[_preenchido]_<data>.csv` | assets e tasks do `subir_predicao_gee.py` |

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

# 3) bordas: diagnóstico e preenchimento com as células vizinhas → mapas_passo<k>_preenchido/
python preencher_bordas.py --preencher --interpolar 10

# 4) conferir no navegador (opcional)
streamlit run app_visualizar.py --server.port 8501 --server.headless true

# 5) subir ao GEE → projects/mapbiomas-arida/mine/cnn_favorabilidade
python juntar_mapas_gee.py            # 1 GeoTIFF por rodada e passo em runs_cnn/upload_gee/
#    upload à mão no Code Editor, depois:
python propriedades_mapas_gee.py
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

### Canais pelo Google Drive (caminho usado)

**Caminho usado: Opção 1.** Os canais foram exportados por tasks do GEE para
o Drive e cortados no Colab; o `baixar_tiles.py` no servidor (acima) e a Opção 2
ficam como alternativas. Fluxo completo e cotas em
[`guia_gdrive.md`](guia_gdrive.md) (seção 0).

**Opção 1 (usada): tasks do GEE (`Export.image.toDrive`), sem servidor e sem rclone.**
Roda em qualquer máquina; uma task por célula (~1 GB cada, o GEE pode dividir
em mais de um arquivo), na pasta `tiles_cnn_gee/` do Drive da conta EE:

```bash
python exportar_tiles_drive.py --modelo A --projeto ee-solkancengine17
python exportar_tiles_drive.py --modelo B1 --projeto ee-solkancengine17
python exportar_tiles_drive.py --modelo A --projeto ee-solkancengine17 --celulas L2C3   # teste
```

Retomável: tasks com a mesma descrição (pendentes, rodando ou concluídas) são
puladas. Acompanhar no Task Manager do GEE ou com `tools/revisarTasksXconta.py`.

Depois, cortar cada célula nos 6 × 6 pedaços com borda que o `inferir_cnn.py
--tiles` lê (no Colab, célula 5 do `cnn_predict_colab.ipynb` ou 10c do
`cnn_colab.ipynb`):

```bash
!python cortar_tiles_drive.py --modelo A --origem /content/drive/MyDrive/tiles_cnn_gee \
    --destino /content/drive/MyDrive/prospectra/tiles_cnn
```

Use o mesmo `--borda` (e `--n`) nos dois scripts. Depois de cortar, os `.tif`
de `tiles_cnn_gee/` podem ser apagados (~18 GB por modelo).

**Opção 2: `baixar_tiles.py` direto para o Drive (rclone).** Cada pedaço sobe
para `gdrive:prospectra/tiles_cnn/canais_<modelo>` e é apagado do disco:

```bash
python baixar_tiles.py --modelo A --destino drive        # outro remote: --remoto ou TILES_REMOTO
```

**Trazer do Drive para o servidor** (para rodar `inferir_cnn.py` no servidor):

```bash
rclone copy gdrive:prospectra/tiles_cnn/canais_A /srv/almacen/db_images/tiles_cnn/canais_A --progress
ls /srv/almacen/db_images/tiles_cnn/canais_A/*/*.tif | wc -l     # 900
```

São ~18 GB por modelo: passa da cota gratuita do Drive (15 GB). Para os pedaços
ficarem na cota de outra conta (pasta compartilhada → pasta própria), use o
`copiar_tifs_colab.ipynb`.

### Predição por rodada (`inferir_cnn.py --tiles padrao`)

Cada rodada lê os pedaços do **seu** modelo de canais sozinha (`--tiles padrao`
olha o `config.json` da rodada): rodadas `A_*` usam `tiles_cnn/canais_A/`,
rodadas `B1_*` usam `tiles_cnn/canais_B1/`. Os pedaços do modelo têm de estar
completos (900) antes.

Sem `final.pt` (rodadas de CV), a predição é a média dos 5 `fold<k>.pt`
(ensemble); nas rodadas `_teste`, usa o `final.pt`.
`--folds auto` usa 1 modelo só (5× mais rápido): o `final.pt` se existir, senão o
1º fold sem ocorrência na validação (`metricas_folds.csv` sem `n_pontos`), que
treinou com todas as ocorrências — o mesmo fold em todas as rodadas, porque a
divisão é a mesma. Não escolher o fold pelo `ap_ponto`: cada fold é medido em
células diferentes, só 2 têm ocorrência e o vencedor treinou sem metade delas.
`--folds 2` força um fold. Não misture modelos na mesma
`mapas_passo<k>/`: se já houver pedaços feitos com outros, o script para —
apague a pasta ou use `--saida`.

Ordem sugerida: `--passo 3` numa célula (teste rápido, confere o mapa) →
`--passo 3` na área toda (comparar as rodadas) → `--passo 1` só na escolhida.

| Rodada | Pedaços | Teste (1 célula, passo 3) | Área toda (passo 3) |
|---|---|---|---|
| `A_pequena_32` | `canais_A` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_pequena_32 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_pequena_32 --tiles padrao --passo 3 --folds auto` |
| `B1_pequena_32` | `canais_B1` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_pequena_32 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_pequena_32 --tiles padrao --passo 3 --folds auto` |
| `A_regnety_64` | `canais_A` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64 --tiles padrao --passo 3 --folds auto` |
| `B1_regnety_64` | `canais_B1` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_regnety_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_regnety_64 --tiles padrao --passo 3 --folds auto` |
| `A_effv2s_64` | `canais_A` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_effv2s_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_effv2s_64 --tiles padrao --passo 3 --folds auto` |
| `B1_effv2s_64` | `canais_B1` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_effv2s_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_effv2s_64 --tiles padrao --passo 3 --folds auto` |
| `A_resnet50_64` | `canais_A` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_resnet50_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_resnet50_64 --tiles padrao --passo 3 --folds auto` |
| `B1_resnet50_64` | `canais_B1` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_resnet50_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_resnet50_64 --tiles padrao --passo 3 --folds auto` |
| `A_coatnet0_64` | `canais_A` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_coatnet0_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_coatnet0_64 --tiles padrao --passo 3 --folds auto` |
| `B1_coatnet0_64` | `canais_B1` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_coatnet0_64 --tiles padrao --passo 3 --folds auto --celulas L2C3` | `python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/B1_coatnet0_64 --tiles padrao --passo 3 --folds auto` |

Saída: `runs_cnn/<rodada>/mapas_passo3/p_<LxCy>.tif` (e `pedacos/` por célula).
Se faltar memória: `--sub 2` (divide cada pedaço em 2 × 2); na GPU, `--lote 128`.

Todas as rodadas treinadas, em sequência, num `tmux` (pula as que não têm `.pt`):

```bash
tmux new -s predict
conda activate py_DL
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
for nome in A_pequena_32 B1_pequena_32 A_regnety_64 B1_regnety_64 A_effv2s_64 B1_effv2s_64 \
            A_resnet50_64 B1_resnet50_64 A_coatnet0_64 B1_coatnet0_64; do
  run=/srv/almacen/db_images/runs_cnn/$nome
  ls $run/*.pt >/dev/null 2>&1 || { echo "pulando $nome (sem .pt)"; continue; }
  python inferir_cnn.py --run $run --tiles padrao --passo 3 --folds auto 2>&1 | tee $run/log_inferir_passo3.txt
done
# Ctrl+b, d para sair sem parar
```

Mapa final (passo 1) — só a rodada escolhida, de preferência a `_teste`:

```bash
python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --tiles padrao --passo 1 \
    2>&1 | tee /srv/almacen/db_images/runs_cnn/A_regnety_64_teste/log_inferir_passo1.txt
```

Conferir quantas células saíram: `ls /srv/almacen/db_images/runs_cnn/<rodada>/mapas_passo3/p_*.tif | wc -l` (25).

### Bordas das células (`preencher_bordas.py`)

Faixas sem p na beira das células: o diagnóstico diz se é buraco do mosaico
(os canais também não têm dado) ou falha da predição (os canais têm dado e o p
não — aí o certo é apagar os pedaços afetados e rodar `inferir_cnn.py` de
novo). O preenchimento grava `mapas_passo<k>_preenchido/` sem tocar nos
originais.

```bash
python preencher_bordas.py                                   # diagnóstico, todas as rodadas com .pt
python preencher_bordas.py --preencher --interpolar 10       # preenche todas (vizinhas + até 10 px interpolados)
python preencher_bordas.py --mapas /srv/almacen/db_images/runs_cnn/B1_coatnet0_64/mapas_passo3 \
    --tiles /srv/almacen/db_images/tiles_cnn/canais_B1      # uma pasta
```

Resumo em `runs_cnn/relatorio_bordas_rodadas.csv`. A interpolação é só
estética (p inventado, não previsto); `--interpolar 0` deixa só as vizinhas.

### Visualizar no servidor (`app_visualizar.py`)

Mapa de p de uma rodada ao lado de um RGB, um canal, outra rodada ou a
diferença entre rodadas, com zoom sincronizado e as ocorrências por cima. Lê os
arquivos locais (sem GEE). Instalar uma vez: `pip install streamlit plotly` no
`py_DL`.

```bash
tmux new -s app
conda activate py_DL
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
streamlit run app_visualizar.py --server.port 8501 --server.headless true
```

No seu computador: `ssh -p 2222 -N -L 8501:localhost:8501 superusuario@192.168.25.250`
e abrir http://localhost:8501.

### Subir ao GEE

Destino: `projects/mapbiomas-arida/mine/cnn_favorabilidade`, banda p × 10000
(0 = sem dado), como o `favorabilidade_ml`. Nos dois caminhos, quando existe
`mapas_passo<k>_preenchido/` sobe só ela (`--com-originais` sobe as duas);
cada imagem leva `rodada`, `passo`, `preenchido` (0/1), `modelo`, `canais`,
`patch`.

**Caminho 1 — à mão pelo Code Editor (sem bucket; o usado agora, porque o
Storage da solkanCengine17 passou da cota).** Uma imagem por rodada e passo
(as 25 células juntas):

```bash
python juntar_mapas_gee.py                 # todas as rodadas → runs_cnn/upload_gee/<rodada>_p<k>[_preenchido].tif
python juntar_mapas_gee.py --run /srv/almacen/db_images/runs_cnn/B1_regnety_64 --passo 3
```

Gera também `upload_gee/propriedades_upload.csv` (asset e propriedades de cada
arquivo). Tamanho esperado: ~30–60 MB por arquivo no passo 3, ~250–500 MB no
passo 1. Baixe os `.tif` para o seu computador e, no Code Editor, Assets → NEW
→ GeoTIFF:

- Asset ID: `projects/mapbiomas-arida/mine/cnn_favorabilidade/<nome do arquivo>`
  (coluna `asset` do CSV)
- nodata 0, pyramiding MEAN

Depois dos uploads, grave as propriedades a partir do CSV (pula os que ainda
não subiram; pode rodar de novo a cada leva):

```bash
python propriedades_mapas_gee.py --conferir     # só mostra
python propriedades_mapas_gee.py
```

A conta EE ativa do servidor precisa escrever em `mapbiomas-arida/mine`; o
projeto EE é o dela (`configure_account_projects_ee.py`; solkanCengine →
`ee-solkancengine17`) ou `--projeto`. O upload chama a banda de `b1`:

```javascript
var p = ee.ImageCollection('projects/mapbiomas-arida/mine/cnn_favorabilidade')
  .filter(ee.Filter.eq('rodada', 'B1_regnety_64'))
  .filter(ee.Filter.eq('passo', 3))
  .filter(ee.Filter.eq('preenchido', 1))
  .first()
  .select(['b1'], ['probability'])
  .divide(10000);
```

**Caminho 2 — bucket + ingestão (`subir_predicao_gee.py`).** Uma imagem por
célula (25 por rodada e passo, `celula` = `LxCy`; no GEE, `.mosaic()`).
Sem argumentos sobe todas as rodadas e passos de `RUNS_DIR`; pula as imagens que
já existem (`--sobrescrever` reenvia).

```bash
python subir_predicao_gee.py                                    # tudo
python subir_predicao_gee.py --passo 3                          # passo 3 de todas
python subir_predicao_gee.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --passo 1
```

Arquivos em `gs://cengine17-acervo-cprm/cnn_favorabilidade/<rodada>/<pasta>/`,
enviados pelo gcloud ADC do servidor (service account
`mapbiomas-energia@mapbiomas-agua`), que precisa de `roles/storage.objectAdmin`
no bucket:

```bash
# na máquina local, conta cengine17
gcloud storage buckets add-iam-policy-binding gs://cengine17-acervo-cprm \
    --member=serviceAccount:mapbiomas-energia@mapbiomas-agua.iam.gserviceaccount.com \
    --role=roles/storage.objectAdmin
```

A ingestão lê o `gs://` com a conta EE ativa: se não for a dona do bucket, dê
a ela `roles/storage.objectViewer`. As tasks aparecem no Task Manager do GEE.
