# Rodar a CNN em dois servidores

| Servidor | IP | GPU | Papel |
|---|---|---|---|
| serverGeoI9 | `192.168.25.250` | RTX 3070, 8 GB | origem: dados, GEE (`py_gee`), subida dos mapas |
| geoserverI5 | `192.168.25.251` | RTX 3060, 12 GB | predição (e treino) de parte das rodadas |

Os comandos supõem o mesmo usuário (`superusuario`) e o mesmo caminho
`/srv/almacen` nos dois; ajuste se for diferente. Os dois são Arch Linux.

**SSH na porta 2222:** o `sshd` do .251 escuta na porta 2222, não na 22. Sem a
porta, a conexão falha com `connect to host 192.168.25.251 port 22: Connection
refused`. O passo 2 cria o apelido `geoserverI5` com a porta; os passos
seguintes usam esse apelido.

**Regra:** cada rodada roda num servidor só. Cada um grava os mapas no próprio
disco, então a mesma rodada nos dois é trabalho em dobro (e arquivos
diferentes para a mesma coisa).

## 1. Preparar o .251

```bash
ssh -p 2222 superusuario@192.168.25.251
nvidia-smi                       # GPU de 12 GB e driver funcionando
df -h /srv/almacen               # ~60 GB livres (pedaços A+B1 ≈ 36 GB + patches + rodadas)
sudo pacman -S --needed rsync tmux
mkdir -p /srv/almacen/db_images /srv/almacen/projects
```

## 2. Acesso por chave do .250 para o .251

O .250 usa a chave `~/.ssh/keyi9` (privada, só no .250) e o .251 só aceita
login por chave (`PasswordAuthentication no`). Por isso o `ssh-copy-id` não
funciona: dá `Permission denied (publickey)`. A chave pública é colada à mão.

**No .250**, chave no lugar e com a permissão certa:

```bash
mkdir -p ~/.ssh && chmod 700 ~/.ssh
mv keyi9 keyi9.pub ~/.ssh/ 2>/dev/null   # se ainda estiverem fora de ~/.ssh
chmod 600 ~/.ssh/keyi9
cat ~/.ssh/keyi9.pub                     # copiar a linha inteira (Ctrl+Shift+c)
```

**No .251** (console ou outra máquina que já entra por chave), colar a linha no
fim de `~/.ssh/authorized_keys`, sem apagar as chaves que já estão lá:

```bash
mkdir -p ~/.ssh && chmod 700 ~/.ssh
nano ~/.ssh/authorized_keys              # Ctrl+End, linha nova, colar, Ctrl+o, Enter, Ctrl+x
chmod 600 ~/.ssh/authorized_keys
ls -l ~/.ssh/authorized_keys             # dono tem que ser superusuario
```

Sem `sudo` no `nano`. Se o arquivo ficar com dono `root`, o `sshd` ignora o
arquivo: `sudo chown superusuario: ~/.ssh/authorized_keys`.

**No .250**, apelido com porta e chave em `~/.ssh/config`:

```bash
cat >> ~/.ssh/config <<'EOF2'
Host geoserverI5 192.168.25.251
    HostName 192.168.25.251
    Port 2222
    User superusuario
    IdentityFile ~/.ssh/keyi9
    IdentitiesOnly yes
EOF2
chmod 600 ~/.ssh/config
ssh geoserverI5 hostname                 # tem que responder sem pedir senha
```

Com o IP também na linha `Host`, `ssh superusuario@192.168.25.251` usa a porta
2222 e a `keyi9` sem opções extras. `IdentitiesOnly yes` impede o `ssh` de
oferecer antes as chaves do agente (`ssh-add -L`), que o .251 recusa.

Sem o `~/.ssh/config`, porta e chave vão em cada comando:
`ssh -i ~/.ssh/keyi9 -p 2222`, `scp -i ~/.ssh/keyi9 -P 2222` (P maiúsculo) e
`rsync -e 'ssh -i ~/.ssh/keyi9 -p 2222'`. Se o login falhar, o motivo aparece
no .251 em `sudo journalctl -u sshd -n 20`.

## 3. Copiar código, imagens e modelos (no .250, dentro de um `tmux`)

```bash
tmux new -s copia

# código (o projeto inteiro, com as mudanças que ainda não estão no git)
rsync -avh --progress /srv/almacen/projects/modelagem_prospectra_mine \
    geoserverI5:/srv/almacen/projects/

# patches e pedaços dos canais
rsync -avh --progress /srv/almacen/db_images/patches_cnn /srv/almacen/db_images/tiles_cnn \
    geoserverI5:/srv/almacen/db_images/

# modelos treinados (.pt, históricos, config), sem os mapas que o .250 está gerando
rsync -avh --progress --exclude 'mapas_passo*' /srv/almacen/db_images/runs_cnn \
    geoserverI5:/srv/almacen/db_images/
```

~40 GB numa rede de 1 Gb/s levam uns 10 a 15 minutos. Dá para copiar com o
loop rodando no .250. Se a cópia for interrompida, rode o mesmo comando: o
`rsync` continua de onde parou.

Depois de mudar um script no .250, copie só o código de novo (primeiro
`rsync`).

## 4. Ambiente Python igual no .251

No .250:

```bash
conda env export -n py_DL --no-builds > /tmp/py_DL.yml
scp /tmp/py_DL.yml geoserverI5:/tmp/
```

No .251:

```bash
conda env create -f /tmp/py_DL.yml
conda activate py_DL
python -c "import torch, timm, rasterio; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
```

Tem que mostrar `True` e o nome da GPU. Se o `conda env create` falhar em algum
pacote, crie o ambiente à mão com `torch` (CUDA), `timm`, `rasterio`, `numpy`,
`pandas` e `scikit-learn`. A predição com `--tiles` não precisa do
`earthengine-api`.

## 5. Dividir as rodadas

| Servidor | Rodadas |
|---|---|
| .250 (8 GB) | `A_pequena_32`, `B1_pequena_32`, `A_regnety_64`, `A_effv2s_64`, `A_resnet50_64`, `A_coatnet0_64` |
| .251 (12 GB) | `B1_regnety_64`, `B1_effv2s_64`, `B1_resnet50_64`, `B1_coatnet0_64` |

**No .250:** se o loop com a lista completa ainda estiver rodando, pare com
`Ctrl+c` depois que a rodada em curso terminar e reinicie só com as rodadas
dele. O que já foi gravado é pulado.

```bash
for nome in A_regnety_64 A_effv2s_64 A_resnet50_64 A_coatnet0_64; do
  run=/srv/almacen/db_images/runs_cnn/$nome
  ls $run/*.pt >/dev/null 2>&1 || { echo "pulando $nome (sem .pt)"; continue; }
  python inferir_cnn.py --run $run --tiles padrao --passo 3 2>&1 | tee -a $run/log_inferir_passo3.txt
done
```

**No .251:**

```bash
tmux new -s predict
conda activate py_DL
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
for nome in B1_regnety_64 B1_effv2s_64 B1_resnet50_64 B1_coatnet0_64; do
  run=/srv/almacen/db_images/runs_cnn/$nome
  ls $run/*.pt >/dev/null 2>&1 || { echo "pulando $nome (sem .pt)"; continue; }
  python inferir_cnn.py --run $run --tiles padrao --passo 3 --lote 512 2>&1 | tee $run/log_inferir_passo3.txt
done
# Ctrl+b, d para sair sem parar
```

Com 12 GB dá para usar `--lote 512`; se aparecer "CUDA out of memory", volte
para 256. VRAM e temperatura: seção "Acompanhar a VRAM" do `README.md`.

### Treinar no .251

Rodadas ainda não treinadas também podem ir para o .251: com o código e os
patches copiados (passo 3), treine lá com `--nome` (comandos em
`rodadas_treino.md`) e traga a pasta da rodada de volta (passo 6). Uma rodada
treinada no .251 não pode ser treinada também no .250.

## 6. Trazer os resultados de volta (no .250, quando o .251 terminar)

```bash
rsync -avh --progress geoserverI5:/srv/almacen/db_images/runs_cnn/ \
    /srv/almacen/db_images/runs_cnn/
```

Sem `--delete`: só acrescenta e atualiza, nada do .250 é apagado. A subida ao
GEE (`subir_predicao_gee.py`) fica no .250, que tem o `py_gee` com as
credenciais.

Conferir se cada rodada tem as 25 células:

```bash
for r in /srv/almacen/db_images/runs_cnn/*/mapas_passo3; do
  echo "$r: $(ls $r/p_L*C*.tif 2>/dev/null | wc -l)"
done
```
