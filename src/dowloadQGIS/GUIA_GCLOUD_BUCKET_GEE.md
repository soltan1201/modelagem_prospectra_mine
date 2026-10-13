# Guia passo a passo: enviar os dados do Acervo CPRM para o Google Earth Engine

Este guia é para quem nunca configurou o Google Cloud. Siga os passos **na ordem**, um de cada vez. Cada passo diz o que digitar, o que vai aparecer na tela e como saber se deu certo.

Conta usada em todo o guia: **`solkan.cengine17@gmail.com`**

Essa conta usa **dois projetos** do Google Cloud, cada um com uma função:

| Projeto | Para que serve | Faturamento |
|---|---|---|
| **Projeto do teste gratuito** (`project-fadb9cba-032b-4fb9-882`) | guarda o **bucket** (a pasta na nuvem) | **Sim:** usa os créditos do teste gratuito (R$ 1.779, até 11/01/2027) |
| **`ee-solkancengine17`** | é o projeto do **Earth Engine** | **Não, e deve continuar assim** |

Por que separar:
- O bucket só pode ser criado num projeto com faturamento. O do teste gratuito já tem.
- O Earth Engine fica num projeto **sem** faturamento. Assim ele nunca tem como cobrar nada: o uso não comercial do GEE é gratuito.
- O único gasto é o bucket: uns 150 MB por poucos dias, **centavos**, descontados dos créditos. **Nada é cobrado no cartão.**

> **Regra de ouro:** no console do Google Cloud, **nunca clique em "Fazer upgrade"**. Enquanto não fizer isso, o teste gratuito continua e só os créditos são usados.

---

## Antes de começar: o que vamos fazer e por quê

Os dados baixados do SGB/CPRM (mapas geológicos, ocorrências de fosfato, grids de geofísica) estão no seu computador, em `/mnt/Almacen/shpIRECE/`. Queremos que eles apareçam no **Google Earth Engine (GEE)** para usar nos modelos.

O GEE não pega arquivos direto do seu computador. O caminho é este:

```
 seu computador  ──►  "bucket" no Google Cloud  ──►  Google Earth Engine
 (arquivos)           (uma pasta na nuvem)            (assets prontos para usar)
```

- **Bucket** é uma pasta na nuvem do Google (Google Cloud Storage). Os arquivos ficam lá só de passagem, e você pode apagá-los no final.
- **gcloud** é o programa de terminal que conversa com o Google Cloud. Ele já está instalado neste computador.
- O **script** `exportar_acervo_gee.py` faz o envio: prepara os arquivos, coloca no bucket e manda o GEE importar.

Para isso, você precisa: (1) entrar com a sua conta no gcloud, (2) dar ao Python permissão para usar essa conta, (3) criar o bucket e (4) rodar o script.

**Tempo estimado:** 20 a 30 minutos.

### Como digitar os comandos

- Abra um **terminal** (a janela preta onde se digitam comandos).
- Copie **uma linha de cada vez**, cole no terminal e aperte **Enter**.
- As linhas que começam com `#` são explicações. **Não precisa digitá-las.**
- Se estiver usando o Claude Code, coloque um `!` antes do comando (exemplo: `! gcloud init`).

---

## Passo 1 — Entrar com a sua conta no gcloud

Digite:

```bash
gcloud init
```

O programa vai fazer algumas perguntas. Responda assim:

| O que aparece na tela | O que responder |
|---|---|
| `Pick configuration to use:` com opções numeradas | Digite o número da opção **`Create a new configuration`** e aperte Enter |
| `Enter configuration name.` | Digite **`cengine17`** e aperte Enter |
| `Choose the account you want to use…` | Escolha o número de **`Log in with a new account`** |
| Abre uma página no navegador | Entre com **`solkan.cengine17@gmail.com`** e clique em **Permitir** |
| `Pick cloud project to use:` | Escolha o projeto do **teste gratuito**: **`project-fadb9cba-032b-4fb9-882`**. Se não aparecer na lista, escolha a opção de digitar o nome e digite `project-fadb9cba-032b-4fb9-882` |
| `Do you want to configure a default Compute Region and Zone?` | Digite **`n`** e aperte Enter |

> **Cuidado:** na primeira pergunta, **não** escolha `Re-initialize this configuration [default]`. Essa opção apaga a configuração da conta da Geodatin, que já existe neste computador.

> **Já tinha feito este passo escolhendo `ee-solkancengine17`?** Não precisa refazer tudo. Digite só:
> ```bash
> gcloud config set project project-fadb9cba-032b-4fb9-882
> ```

**Deu certo se:** ao digitar o comando abaixo, aparecer uma estrela `*` ao lado de `solkan.cengine17@gmail.com`.

```bash
gcloud auth list
```

---

## Passo 2 — Dar ao Python permissão para usar a sua conta

O passo 1 deu acesso ao **terminal**. O script é feito em **Python**, que precisa de uma autorização própria. Digite:

```bash
gcloud auth application-default login
```

Vai abrir o navegador de novo. Entre com **`solkan.cengine17@gmail.com`** e clique em **Permitir**.

Depois digite:

```bash
gcloud auth application-default set-quota-project project-fadb9cba-032b-4fb9-882
```

O código do projeto do teste gratuito é **`project-fadb9cba-032b-4fb9-882`**. Se quiser conferir, `gcloud projects list` mostra esse código na coluna `PROJECT_ID`.

**Deu certo se:** aparecer uma mensagem dizendo que as credenciais foram salvas (`Credentials saved to file…`).

> **Importante:** isso cria um arquivo que funciona como uma **senha** da sua conta (`~/.config/gcloud/application_default_credentials.json`). Não envie esse arquivo para ninguém e não coloque no GitHub.

---

## Passo 3 — Conferir o faturamento dos dois projetos

**3.1 — O projeto do teste gratuito precisa ter faturamento.** Digite:

```bash
gcloud billing projects describe project-fadb9cba-032b-4fb9-882
```

Deve aparecer **`billingEnabled: true`**. É o faturamento do teste gratuito, pago com os créditos.

**3.2 — O projeto do Earth Engine deve ficar SEM faturamento.** Digite:

```bash
gcloud billing projects describe ee-solkancengine17
```

Deve aparecer **`billingEnabled: false`**. **Está certo assim; não ligue o faturamento nele.**

**3.3 — Ativar o serviço de armazenamento no projeto do teste:**

```bash
gcloud services enable storage.googleapis.com --project=project-fadb9cba-032b-4fb9-882
```

**Deu certo se:** o comando terminar sem mensagem de erro (às vezes demora um pouco).

**3.4 — (Recomendado) Criar um alerta de gastos.** No console (<https://console.cloud.google.com>), na tela inicial do teste gratuito:
1. Clique em **"Configurar alertas de orçamento"**.
2. Crie um orçamento de **R$ 10**.

Se algum dia houver um gasto acima do esperado, chega um e-mail.

**3.5 — (Opcional) Confirmar que o Earth Engine está como uso não comercial.** Abra <https://console.cloud.google.com/earth-engine/configuration?project=ee-solkancengine17>. Deve aparecer **Noncommercial** (não comercial).

---

## Passo 4 — Criar o bucket (a pasta na nuvem)

Digite (é **um comando só**, em várias linhas; copie tudo junto):

```bash
gcloud storage buckets create gs://cengine17-acervo-cprm \
    --project=project-fadb9cba-032b-4fb9-882 \
    --location=US-CENTRAL1 \
    --uniform-bucket-level-access
```

- `cengine17-acervo-cprm` é o nome do bucket. Ele precisa ser **único no mundo inteiro**, só com letras minúsculas, números e hífen.
- Se aparecer erro dizendo que o nome **já existe** (`409`), troque por outro, por exemplo `cengine17-acervo-cprm-irece`. Use esse novo nome no resto do guia.

**Teste se o bucket funciona.** Digite os três comandos, um de cada vez:

```bash
echo teste | gcloud storage cp - gs://cengine17-acervo-cprm/acervo_cprm/teste.txt
gcloud storage ls gs://cengine17-acervo-cprm/acervo_cprm/
gcloud storage rm gs://cengine17-acervo-cprm/acervo_cprm/teste.txt
```

**Deu certo se:** o segundo comando mostrar `gs://cengine17-acervo-cprm/acervo_cprm/teste.txt`.

---

## Passo 5 — Dizer ao script qual é o bucket e o projeto dele

1. Abra o arquivo `src/dowloadQGIS/exportar_acervo_gee.py` no editor (VS Code, por exemplo).
2. Perto do começo, procure estas duas linhas:
   ```python
   BUCKET     = ''
   PROJETO_STORAGE = ''
   ```
3. Troque por (**já deixei preenchido assim no script**; só confira):
   ```python
   BUCKET     = 'gs://cengine17-acervo-cprm/acervo_cprm'
   PROJETO_STORAGE = 'project-fadb9cba-032b-4fb9-882'
   ```
4. Salve o arquivo.

- O `BUCKET` diz **onde** colocar os arquivos. Se ficar vazio, o script pergunta ao rodar.
- O `PROJETO_STORAGE` diz **qual projeto** é dono do bucket: o do teste gratuito. **Não deixe vazio:** vazio, o script usaria o `ee-solkancengine17`, que não tem faturamento.
- O Earth Engine continua usando o `ee-solkancengine17`; isso não muda.

---

## Passo 6 — Preparar os arquivos (ainda não envia nada)

No terminal, entre na pasta do script e rode:

```bash
cd /home/superuser/Dados/proj/prosprectra/modelagem_prospectra_mine/src/dowloadQGIS
python exportar_acervo_gee.py
```

Quando perguntar `Etapa [1]:`, apenas aperte **Enter** (a etapa 1 é só preparar).

**Deu certo se:** aparecer a pasta `/mnt/Almacen/shpIRECE/gee_upload/` com estes arquivos:

- `1042_gama.tif` e `1042_mag.tif` (geofísica)
- `ocorrencias_fosfato.zip`, `recursos_minerais.zip`, `geoquimica_p2o5.zip`
- `lito_jacobina_250k.zip`, `estrut_jacobina_250k.zip`, `lito_bahia.zip`, `estrut_bahia.zip`, `dolinas_morro_chapeu.zip`

**Recomendado:** abra alguns desses arquivos no **QGIS** para ver se estão no lugar certo do mapa antes de enviar. Os `.zip` abrem arrastando para o QGIS.

---

## Passo 7 — Enviar para o Google Earth Engine

Rode o script de novo:

```bash
python exportar_acervo_gee.py
```

Responda assim:

| Pergunta | Resposta |
|---|---|
| `Etapa [1]:` | **`2`** |
| `Bucket do Cloud Storage…` | Aperte **Enter** (se fez o passo 5) ou digite `gs://cengine17-acervo-cprm/acervo_cprm` |
| `Conta [4]:` | **`8`** (é a `solkanCengine`, ou seja, `solkan.cengine17@gmail.com`) |
| `Se o asset já existir, sobrescrever? (s/n) [n]:` | Aperte **Enter** (não) |

O que vai aparecer, e é normal:

- `gee.PREVIOUS_USER … gee.NEW_USER: solkanCengine`: o script trocou para a conta certa.
- `credenciais do EE sem acesso ao Storage … tentando gcloud ADC`: **é esperado.** O script passa a usar a autorização do passo 2.
- Uma linha por arquivo, terminando em `(task …)`: cada arquivo foi enviado e o GEE começou a importar.

No fim aparece o nome de um arquivo `registro_exportar_gee_<data>.csv`, com a lista do que foi enviado.

---

## Passo 8 — Esperar o GEE terminar e conferir

A importação no GEE demora de alguns minutos a algumas horas. Para acompanhar:

- **Pelo navegador:** abra <https://code.earthengine.google.com>, entre com `solkan.cengine17@gmail.com` e clique na aba **Tasks** (à direita). Cada arquivo aparece como uma tarefa. Espere todas ficarem **azuis/concluídas** (`COMPLETED`).
- **Pelo terminal:**
  ```bash
  cd ../tools
  python revisarTasksXconta.py
  ```
  Responda: única conta = `Y`, número da conta = **`9`** (solkanCengine), sem filtro, mostrar = `2` (todas).

Quando terminar, os dados estarão no GEE em:

```
projects/mapbiomas-arida/mine/acervo_cprm/
```

---

## Passo 9 — Limpar o bucket (opcional, mas recomendado)

Depois que **todas** as tarefas estiverem concluídas, os arquivos no bucket não são mais necessários. Apague-os para não ter custo:

```bash
gcloud storage rm -r gs://cengine17-acervo-cprm/acervo_cprm
```

O bucket vazio pode ficar; ele não custa nada e serve para a próxima vez.

---

## Se algo der errado

| Mensagem que aparece | O que significa | O que fazer |
|---|---|---|
| `Reauthentication required` ou `RefreshError` | a autorização expirou | Refaça os passos 1 e 2 |
| `billing account … disabled in state absent` ao criar o bucket | o bucket foi criado no projeto errado (`ee-solkancengine17`, sem faturamento) | Crie no projeto do teste: `--project=project-fadb9cba-032b-4fb9-882` (passo 4) |
| `billing account … disabled in state absent` ao **enviar** (passo 7) | a requisição foi cobrada do `ee-solkancengine17` (sem faturamento) | Confira o `PROJETO_STORAGE` (passo 5). O script já força esse projeto nas credenciais; por garantia, rode `gcloud auth application-default set-quota-project project-fadb9cba-032b-4fb9-882` |
| O console avisa que o teste gratuito acabou | acabaram os créditos ou passou de 11/01/2027 | Não faça upgrade; use o Code Editor do GEE para subir os arquivos à mão |
| `409 … already exists` ao criar o bucket | alguém no mundo já usa esse nome | Escolha outro nome no passo 4 e use-o no passo 5 |
| `403 … storage.objects.create` ao enviar | o Python está usando outra conta | Refaça o passo 2 com `solkan.cengine17@gmail.com` |
| Tarefa no GEE com `FAILED … cannot read gs://` | o GEE não conseguiu ler o bucket | Confira se escolheu a conta **8** no passo 7 |
| `Asset already exists` | esse dado já foi enviado antes | Tudo bem, ele foi pulado. Para enviar de novo, responda `s` em "sobrescrever" |
| `Permission denied` na pasta `mapbiomas-arida/mine` | a conta não pode gravar nessa pasta | Peça acesso de escrita a quem administra o projeto `mapbiomas-arida` |
| `Nada em … gee_upload` | a etapa 1 não foi feita | Faça o passo 6 antes do passo 7 |

---

## Anexo A — Se o teste gratuito não estiver disponível

Dá para usar um bucket de **outra conta que já tenha faturamento** (por exemplo, a da Geodatin, projeto `geo-data-s`) e dar permissão para a `solkan.cengine17` usá-lo.

1. Volte para a configuração da Geodatin e entre com ela:
   ```bash
   gcloud config configurations activate default
   gcloud auth login solkan@geodatin.com
   ```
2. Crie o bucket nesse projeto:
   ```bash
   gcloud storage buckets create gs://geodatin-acervo-cprm --project=geo-data-s --location=US-CENTRAL1 --uniform-bucket-level-access
   ```
3. Dê permissão à `solkan.cengine17` para gravar e ler nesse bucket:
   ```bash
   gcloud storage buckets add-iam-policy-binding gs://geodatin-acervo-cprm \
       --member=user:solkan.cengine17@gmail.com --role=roles/storage.objectAdmin
   ```
4. Volte para a configuração da `cengine17`:
   ```bash
   gcloud config configurations activate cengine17
   ```
5. Continue do **passo 5**, usando `gs://geodatin-acervo-cprm/acervo_cprm` como bucket.

## Anexo B — Voltar a usar a conta da Geodatin no gcloud

```bash
gcloud config configurations activate default
gcloud auth application-default login                     # entrar com solkan@geodatin.com
gcloud auth application-default set-quota-project geo-data-s
```

## Anexo C — Se o Earth Engine reclamar da conta `solkanCengine`

O computador já tem a autorização do Earth Engine dessa conta guardada em `~/.config/earthengine/solkanCengine/`. Se algum dia aparecer erro de credencial do Earth Engine:

```bash
earthengine authenticate --force          # entrar com solkan.cengine17@gmail.com
cp ~/.config/earthengine/credentials ~/.config/earthengine/solkanCengine/credentials
```
