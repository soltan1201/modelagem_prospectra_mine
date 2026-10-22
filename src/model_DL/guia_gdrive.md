# Guia — subir os patches (e as rodadas) do servidor para o Google Drive

Para treinar no Colab com os patches exportados no servidor
(`exportar_patches.py`). Caminho recomendado: **rclone** — copia pastas
inteiras, retoma se a conexão cair e funciona num servidor sem navegador.

Destino usado pelo `cnn_colab.ipynb` (célula 2):
`/content/drive/MyDrive/prospectra/patches_cnn` → no Drive, `prospectra/patches_cnn`.

---

## 1. Instalar o rclone no servidor

```bash
curl https://rclone.org/install.sh | sudo bash      # ou: sudo apt install rclone
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

- Rodadas de treino feitas no servidor (para comparar no Colab):

  ```bash
  rclone copy /srv/almacen/db_images/runs_cnn gdrive:prospectra/runs_cnn --progress
  ```

- Para não depender da sessão SSH: rode dentro de `tmux` (ou com `nohup … &`).

Tamanho esperado: ~0,6–0,7 GB por modelo (A e B1 ≈ 1,3 GB) — cabe na cota
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
