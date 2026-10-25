# Guia — subir os patches (e as rodadas) do servidor para o Google Drive

Para treinar no Colab com os patches exportados no servidor
(`exportar_patches.py`). Caminho recomendado: **rclone** — copia pastas
inteiras, retoma se a conexão cair e funciona num servidor sem navegador.

Destino usado pelo `cnn_colab.ipynb` (célula 2):
`/content/drive/MyDrive/prospectra/patches_cnn` → no Drive, `prospectra/patches_cnn`.

---

## 0. Fluxo completo sem os servidores (máquina local + Drive + Colab)

Quando as GPUs dos servidores estão ocupadas: os downloads do GEE rodam em
qualquer máquina (só CPU e rede) e vão direto para o Drive; treino e mapa
rodam no Colab lendo do Drive.

| Etapa | Onde | Comando | Resultado no Drive |
|---|---|---|---|
| 1. patches do zero | máquina local | `python exportar_patches_drive.py --projeto ee-solkancengine17` (tasks do GEE, sem rclone) | `patches_cnn_gee/` (TFRecord) |
| 1b. converter | Colab | célula 6b (`tfrecord_para_npz.py`) | `prospectra/patches_cnn/` (~1,3 GB) |
| 2. tiles modelo A | máquina local + Colab | `python exportar_tiles_drive.py --modelo A --projeto ee-solkancengine17`, depois célula 10c (`cortar_tiles_drive.py`) | `prospectra/tiles_cnn/canais_A/` (~18 GB) |
| 3. tiles modelo B1 | máquina local + Colab | `python exportar_tiles_drive.py --modelo B1 --projeto ee-solkancengine17`, depois célula 10c (`cortar_tiles_drive.py`) | `prospectra/tiles_cnn/canais_B1/` (~18 GB) |
| 4. treino e teste | Colab | `cnn_colab.ipynb` células 2–5, pule 6–7, depois 8–10 | `prospectra/runs_cnn/<rodada>/` |
| 5. mapa | Colab | célula 11 com `USAR_TILES` marcado | `runs_cnn/<rodada>/mapas_passo<k>/` |
| 6. subir ao GEE | máquina local | `rclone copy` da pasta de mapas + `subir_predicao_gee.py --run <pasta local>` | coleção `cnn_favorabilidade` |

- O remote do rclone tem de ser a conta Google que o Colab monta. Os scripts
  usam `gdrive:`; se o remote tem outro nome, crie um apelido uma vez:
  `rclone config create gdrive alias remote=<nome_do_remote>:`
- Projeto do GEE (cota dos pedidos): `export PROJETO_EE=ee-solkancengine17`
  antes das etapas 1–3 (vale para `exportar_patches.py`, `baixar_tiles.py` e
  `inferir_cnn.py`; `--projeto` no `baixar_tiles.py` também serve). A conta
  precisa ler os assets de `projects/mapbiomas-arida/mine/`.
- A pasta local é só de passagem (`~/.cache/prospectra/...`): cada arquivo é
  apagado depois de subir. Rode cada etapa em `tmux`; se cair, rode o mesmo
  comando de novo e ele continua de onde parou.
- Etapas 1–3 são independentes: podem rodar ao mesmo tempo em terminais
  separados (cuidado com a cota de pedidos do GEE: se aparecer erro de cota,
  baixe `--threads`).
- Espaço: ~38 GB no total — exige conta com mais que os 15 GB gratuitos.
- A opção 3 da etapa 1 apaga os patches, o índice e a divisão antigos no Drive
  (pede confirmação digitando `apagar`). Com o mesmo `SEED` e os mesmos assets
  (pontos, grid, `favorabilidade_ml`), o sorteio e a divisão saem iguais aos
  anteriores; se algum asset mudou, as rodadas antigas de `runs_cnn/` deixam de
  ser comparáveis com as novas.

---

## 1. Instalar o rclone no servidor

```bash
sudo pacman -S rclone          # servidor é Arch Linux
rclone version
```

## 2. Conectar ao Drive (uma vez só)

O servidor não tem navegador: a autorização é feita na sua máquina.

**No servidor:**

```bash
rclone config
#  n) New remote          → nome: gdrive
#  Storage                → drive
#  client_id / secret     → Enter (vazio)
#  scope                  → 1 (acesso total)
#  service_account_file   → Enter
#  Edit advanced config?  → n
#  Use web browser to automatically authenticate? → n   (servidor sem navegador)
#  ele mostra um comando:  rclone authorize "drive" "..."
```

**Na sua máquina** (com navegador e rclone instalado):

```bash
rclone authorize "drive" "..."     # cole o comando mostrado pelo servidor
# abre o navegador → autorize com a conta Google → ele imprime um token
```

**De volta ao servidor:** cole o token, responda `n` para *Shared Drive* e
confirme. Teste:

```bash
rclone lsd gdrive:
```

## 3. Subir os patches

```bash
rclone copy /srv/almacen/db_images/patches_cnn gdrive:prospectra/patches_cnn \
    --transfers 8 --checkers 8 --progress
```

- `copy` só envia o que falta ou mudou: depois de exportar mais células, rode
  de novo e sobem só as novas.
- Conferir se chegou tudo igual:

  ```bash
  rclone check /srv/almacen/db_images/patches_cnn gdrive:prospectra/patches_cnn
  ```

- Ou exporte direto para o Drive: `python exportar_patches.py`, opção 2 e
  destino 2 (Google Drive). Cada `.npz` sobe por rclone para
  `gdrive:prospectra/patches_cnn` (outro destino: variável `PATCHES_REMOTO`) e
  sai do disco; `indice_patches.csv` e `divisao.json` ficam no servidor e são
  copiados. Rodar de novo retoma de onde parou.

- Rodadas de treino feitas no servidor (para comparar no Colab):

  ```bash
  rclone copy /srv/almacen/db_images/runs_cnn gdrive:prospectra/runs_cnn --progress
  ```

- Pedaços de canais para a predição (`baixar_tiles.py`): em vez de baixar no
  servidor e copiar, mande direto — cada pedaço sobe e é apagado do disco:

  ```bash
  python baixar_tiles.py --modelo A --destino drive     # → gdrive:prospectra/tiles_cnn/canais_A
  ```

  São ~18 GB por modelo: passa da cota gratuita (15 GB).

- Para não depender da sessão SSH: rode dentro de `tmux` (ou com `nohup … &`).

Tamanho esperado dos patches: ~0,6–0,7 GB por modelo (A e B1 ≈ 1,3 GB) — cabe na cota
gratuita do Drive (15 GB).

## 4. No Colab

- Célula 2 do `cnn_colab.ipynb`: `PATCHES_DIR` já aponta para
  `/content/drive/MyDrive/prospectra/patches_cnn` (ajuste se usou outro destino).
- Com os patches no Drive, **pule as células 6–7** (exportação).
- A célula 8 copia os patches para o disco local do Colab antes do treino
  (`COPIAR_LOCAL`), porque ler do Drive a cada época é lento.

---

## Alternativas sem rclone

**Navegador (manual):**

```bash
# no servidor: junta tudo num arquivo (os .npz já são comprimidos)
tar -cf patches_cnn.tar -C /srv/almacen/db_images patches_cnn
```

Baixe o `.tar` para a sua máquina, envie pelo drive.google.com e, no Colab:

```python
!tar -xf /content/drive/MyDrive/prospectra/patches_cnn.tar -C /content/
```

(e use `PATCHES_DIR = '/content/patches_cnn'`). Mais trabalhoso a cada
atualização.

**Exportar direto no Colab:** as células 6–7 do `cnn_colab.ipynb` baixam os
patches do GEE já para o Drive, sem passar pelo servidor — mais lento, mas não
exige nada no servidor.

## Problemas comuns

| Sintoma | Causa / solução |
|---|---|
| `rclone lsd gdrive:` pede login de novo | token expirou ou foi revogado: `rclone config reconnect gdrive:` |
| `403 rateLimitExceeded` | muitas transferências ao mesmo tempo: `--transfers 4 --tpslimit 5` |
| pasta aparece em "Compartilhados comigo" e não em "Meu Drive" | o destino era um Shared Drive; refaça o remote respondendo `n` para *Shared Drive* |
| Colab não acha os arquivos | caminho da célula 2 diferente do destino do `rclone copy` |
