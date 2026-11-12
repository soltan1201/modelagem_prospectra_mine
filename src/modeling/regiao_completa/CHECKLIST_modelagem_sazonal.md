# Checklist — modelagem de favorabilidade no mosaico sazonal harmonizado

Mosaico de entrada: `projects/mapbiomas-arida/mine/aster-mosaico-final_sazonal_hm_razao`
(ImageCollection, uma imagem por célula, 15 m, mesma grade do v1). Reflectância de
**superfície** × 1000 (o v1 é TOA × 1000), TIR em DN, banda `origem` (1 = A, 2 = B onde
falta A, 3 = B no lugar do A com nuvem) e `swir_estimado`.

Nome curto do mosaico nos scripts: `aster-mosaico-final_sazonal_hm_razao` (no notebook,
`'sazonal'`). Já é a opção padrão dos menus (`models_ML_fosfato.MOSAICO_PADRAO`).

O usuário roda tudo; marque `[x]` e anote data e números em **Notas**.

---

## Parte 1 — Maxent, RF e GTB (região completa, grid 5 × 5)

Fluxo: amostras (2) → CSV (3) → exploração `avaliar_maxent.ipynb` (4) → seleção de
variáveis e do fundo `selecao_variaveis_RC.ipynb` (4b) → decisão do vazamento (5) →
mapas RF/GTB por célula (6). Tudo em `src/modeling/regiao_completa/`; validação sempre
com uma célula da grid 5 × 5 de fora; export uma task por célula.

### 0. Mosaico pronto

- [x] Coleção com as 25 células (tasks `aster-mosaico-final-sazonal_hm_razao_<LxCy>` concluídas em 2026-11-02; Code Editor:
      `ee.ImageCollection('projects/mapbiomas-arida/mine/aster-mosaico-final_sazonal_hm_razao').size()`)
- [ ] Conferir a banda `origem` = 3 nas nuvens do A (`apresentacao_seco_x_harmonizado.js`,
      composição "Origem do pixel")
- [ ] Ciente: `swir_estimado` desse mosaico está 0 na origem 3 (bug corrigido no build);
      nas amostras é só coluna de análise, e o DL recalcula pela `origem`

### 1. Pontos de referência

- [x] `pontos_referencia` (1.103 pontos) — não muda com o mosaico

### 2. Amostras (`models_ML_fosfato.py`, opção 1)

```bash
cd src/modeling/regiao_completa
python models_ML_fosfato.py      # Mosaico: aster-mosaico-final_sazonal_hm_razao · O que fazer: 1 (amostras)
```

- [x] Task concluída (2026-11-03 02:40 UTC) → `mine/maxent_amostras_aster-mosaico-final_sazonal_hm_razao_bandas_indices`
      (31.849 linhas: fosfato 740, carbonato 6.540, mvt 5.040, ouro 4.160, manganês 5.369, fundo 10.000;
      v1 31.860, só manganês −11)
      (presença 20 px por ponto, fundo 10 mil px, 36 features + `origem`/`swir_estimado`)
- [x] As 14 bandas, 17 índices, 4 GLCM, `origem` e `swir_estimado` na tabela

### 3. Baixar o CSV (`baixar_amostras_maxent.py`)

```bash
python baixar_amostras_maxent.py projects/mapbiomas-arida/mine/maxent_amostras_aster-mosaico-final_sazonal_hm_razao_bandas_indices
```

- [x] `src/Dados/maxent/maxent_amostras_aster-mosaico-final_sazonal_hm_razao_bandas_indices.csv`
      (2026-11-03 08:32, da task de 02:40 UTC: 31.849 linhas, 0 duplicadas, NaN só em
      `fonte`/`papel` do fundo)

### 4. Avaliação (`avaliar_maxent.ipynb`)

Célula [2]: `MOSAICO = 'sazonal'`, `MOSAICOS_USAR = ['v1', 'sazonal']` (as seções que rodam
todos os mosaicos comparam v1 × sazonal; o único saiu da lista).

- [ ] Seções 2–6 no sazonal: pontos, separação por feature, Maxent célula de fora
- [ ] [23] / [24b]: AUC por ponto, Boyce e top 10 % — Maxent, RF e GTB, v1 × sazonal
      (referência v1, fosfato: RF 0,83, GTB 0,84)
- [ ] [25] / [26]: correlação entre alvos e p(outros alvos) — só para conferir se mudou
- [ ] [30]–[34] SHAP no sazonal → ranking novo
      (`models/ranking_features_shap_<data>.csv`)
- [ ] Decidir `FEATURES`: manter as 21 do v1 ou trocar pelo ranking do sazonal. Se
      trocar, mudar **junto** `exportar_favorabilidade_ml.py` (`FEATURES`) e
      `maxent_alvos.js` (`FEATURES_DL`). Para comparar v1 × sazonal no mesmo pé, manter
      as 21 na primeira rodada.
- [ ] Anotar a métrica de referência do sazonal (AUC por ponto do GTB): é a barra que a
      CNN precisa passar

### 4b. Seleção de variáveis e do fundo (`selecao_variaveis_RC.ipynb`)

Mesmo desenho do `selecao_variaveis_AR.ipynb`, com folds = células da grid 5 × 5.
Célula [2]: `MOSAICO = 'aster-mosaico-final_sazonal_hm_razao'` (sufixo `shm`),
`FUNDO_MODO = 'cluster'`, `PAPEIS = None`.

- [x] [2]–[10e]: folds por célula, fundo por cluster → `models/cluster_fundo_RC_shm.json`
      (K = 10, regra `n_distantes`, N = 2; treino nos clusters 0, 2, 5, 7, 8, 9 = 3.834 de
      10.000 px). Atenção: o c8 tem só 2 px em 2 células e o c2 53 px em 9 — os 2 mais
      distantes quase não têm pixels; o fundo vem da completação até 20 % (c0, c5, c7, c9)
- [x] [12] ablação (ver Notas) — falta preencher a tabela da [12a]
- [x] [14]–[15] ranking estável na base `sem_swir` → `SELECAO` com 11 variáveis (as 11
      do `sem_swir` têm freq_topo = 1; o `N_MIN_SEL = 15` não tem de onde completar)
- [x] [17] avaliação nos pontos → `models/avaliacao_pontos_RC_shm_selecao_cluster.csv`
- [x] [19] validação R2 (treina na bacia, testa fora): `selecao` × `fixas_v1`
- [x] [21] → `models/variaveis_RC_shm.json` (recomendado = `selecao`, 2026-11-03 22:50)
- [ ] Opcional: rodar de novo com `FUNDO_MODO = 'aleatorio'` para comparar (CSVs e
      figuras levam o modo no nome; o JSON é sobrescrito — rodar a [21] por último no
      modo escolhido)

### 5. Decisão: vazamento no teste da CNN

Os mapas RF/GTB viram pseudo-rótulos da CNN. Com `PAPEIS = None` eles usam também os
pontos de validação (R2), e o teste da CNN deixa de ser independente.

- [ ] Escolher: (a) `PAPEIS = None`, como no v1, e citar a limitação; ou (b) mapas com
      `PAPEIS = ['treino']` numa coleção separada e `ASSET_FAV` do `exportar_patches.py`
      apontando para ela. O (b) exige um ajuste no `exportar_favorabilidade_ml.py`
      (`ASSET_SAIDA`) — pedir antes de rodar.

### 6. Mapas RF e GTB (`exportar_favorabilidade_ml.py`)

```bash
python exportar_favorabilidade_ml.py
#   Modelos: os dois · Mosaicos: aster-mosaico-final_sazonal_hm_razao
#   Conjunto: selecao (recomendado do variaveis_RC_shm.json; 'fixas' = 21 do v1) · Células: vazio (todas)
```

- [ ] Plano impresso: conjunto, variáveis, fundo `cluster [...]`, exemplo de nome
      `RF_fosfato_aster-mosaico-final_sazonal_hm_razao_selecao_fcl_L1C1`
- [ ] Log sem o aviso "fundo elegível no GEE ≠ notebook" (se aparecer, o asset de
      amostras mudou depois do notebook)

- [ ] Contas: células repartidas entre caatinga02 (13) e caatinga05 (12), cada uma no
      projeto dela (`CONTAS` no `models_ML_fosfato.py`); o registro tem a coluna `conta`
- [ ] 50 tasks (2 modelos × 25 células) → `mine/favorabilidade_ml/<RF|GTB>_fosfato_aster-mosaico-final_sazonal_hm_razao_selecao_fcl_<LxCy>`
- [ ] Registro `registro_favorabilidade_ml_<data>.csv` com n_presenca e n_fundo
- [ ] Contagem no GEE: 25 RF + 25 GTB com `mosaico = aster-mosaico-final_sazonal_hm_razao` e
      `conjunto = selecao`
- [ ] Se também exportar `fixas` no mesmo mosaico: `exportar_patches.py` e
      `visualizar_favorabilidade_ml.js` filtram só por modelo/alvo/mosaico — pedir o filtro
      por `conjunto` antes de usar os mapas

### 7. Conferir os mapas (`js/visualizar_favorabilidade_ml.js`)

- [ ] Mosaico do seletor = sazonal (padrão)
- [ ] Áreas magenta (RF **e** GTB ≥ 0,5): seguem a geologia ou limites de cena/célula?
      No v1 havia manchas de borda reta; o sazonal deve ter menos
- [ ] Curva área × limiar e histogramas: anotar a % da área com p ≥ 0,5 nos dois
      (define quantos pseudo-positivos a CNN vai ter)
- [ ] Prints para o relatório (`docs/relatorio_avaliacao_modelos/image/`)

### 7b. Área reduzida Irecê — ASTER + geofísica + geologia (`*_AR`)

Agora com o harmonizado (`amostras_AR.MOSAICO`); nomes com `irece_shm` (os do v1, `irece`,
ficam intactos). Contas: caatinga02 / caatinga05.

```bash
python amostras_AR.py                       # task (caatinga02) → mine/amostras_AR_irece_shm
python amostras_AR.py --baixar              # → src/Dados/maxent/amostras_AR_irece_shm.csv + _grupos.json
```

- [x] Amostras exportadas e baixadas (task 2026-11-03 00:06 UTC; `amostras_AR_irece_shm.csv`
      18,6 MB + `_grupos.json`)
- [x] `selecao_variaveis_AR.ipynb` inteiro (`AREA = 'irece_shm'`): cluster do fundo [10e],
      ablação, ranking estável, conjunto `selecao` → `models/variaveis_AR_irece_shm.json`
      e `models/cluster_fundo_AR_irece_shm.json`
- [x] Comparar a `selecao` nova com a do v1 (ver Notas)
- [x] `python exportar_favorabilidade_AR.py` (conjunto `selecao`, os dois modelos): RF na
      caatinga02, GTB na caatinga05 → `favorabilidade_AR/<RF|GTB>_fosfato_AR_irece_shm_selecao_fcl`
- [x] Tasks concluídas: `RF_...`, `GTB_..._shm_selecao_fcl` e `cluster_fundo_AR_irece_shm`
- [x] `favorabilidade_AR` com 4 imagens (2 do v1 + 2 `_shm`)
- [ ] Conferir no visualizador

### 8. Opcional

- [x] Maxent no GEE (`models_ML_fosfato.py` opção 3) → `mine/maxent_fosfato_aster-mosaico-final_sazonal_hm_razao`
      (IMAGE, 2026-11-03); validação da opção 2 sem CSV em `models/` (não rodada ou não salvou)
- [ ] Outros alvos (`ALVOS` no `exportar_favorabilidade_ml.py`)
- [ ] Relatório `docs/relatorio_avaliacao_modelos/`: seção v1 × sazonal

---

## Parte 2 — CNN (`src/model_DL/regiao_completa/`; `.md` em `src/model_DL/`)

Entradas já apontadas para o mosaico novo (`exportar_patches.py`: `MOSAICO`,
`ASSET_MOSAICO`, `MOSAICO_FAV`); `baixar_tiles`, `exportar_tiles_drive`,
`exportar_patches_drive` e `inferir_cnn` herdam de lá. Detalhes em
`plano_cnn.md` (seção "Modelo no mosaico sazonal") e `rodadas_treino.md`.

### 9. Guardar os dados do v1 (antes de exportar qualquer coisa)

Os scripts gravam nas mesmas pastas; sem isto os patches do sazonal sobrescrevem os do v1.

```bash
cd /srv/almacen/db_images
mv patches_cnn patches_cnn_v1
mv tiles_cnn   tiles_cnn_v1
mv runs_cnn    runs_cnn_v1
# Drive (se usou o caminho do Colab)
rclone moveto gdrive:prospectra/patches_cnn gdrive:prospectra/patches_cnn_v1
rclone moveto gdrive:prospectra/tiles_cnn   gdrive:prospectra/tiles_cnn_v1
rclone moveto gdrive:prospectra/runs_cnn    gdrive:prospectra/runs_cnn_v1
```

- [ ] Servidor `.250` (e `.251`, se tiver cópia — `guia_dois_servidores.md`).
      2026-11-05: **não feito** — `patches_cnn/` já tem os patches do harmonizado (os do v1
      foram sobrescritos; refazer com `MOSAICO = 'aster-mosaico-final'` se precisar);
      `runs_cnn/` tem v1 e `_shm` lado a lado (nomes não colidem); `tiles_cnn/` ainda é o
      v1 → mover antes do passo 14
- [ ] Drive (inclui `patches_cnn_gee/` e `tiles_cnn_gee/` do GEE, se existirem)
- [ ] Para olhar o v1 depois: `RUNS_DIR=.../runs_cnn_v1 TILES_DIR=.../tiles_cnn_v1`

### 10. Sorteio dos patches (`exportar_patches.py`, opção 1)

Pseudo-positivos = RF **e** GTB do conjunto `selecao` com p ≥ `P_POS = 0.94` nos dois (o
magenta do visualizador). Negativos = `NEG_MODO = 'cluster'`: pixels nos clusters distantes do
`cluster_fundo_RC_shm.json` (os mesmos do fundo do RF/GTB), a > 2 km das ocorrências e com
RF e GTB < 0,30. Fundo aleatório continua só para avaliação.

- [x] Só sortear: 10.435 centros — ambos 2.518, negativo (cluster) 3.335, negativo_dificil
      2.332, fundo 1.917, ocorrencia_treino 117 (13 × 9), ocorrencia_teste 216 (24 × 9);
      treino 9.023, teste 1.412
- [x] `P_POS = 0.94` mantido (área RF e GTB ≥ 0,94 = 1,8 % da região)
- [x] `indice_patches.csv` e `divisao.json` gravados (`MOSAICO` harmonizado, `NEG_MODO = cluster`,
      `CONJUNTO_FAV = selecao`); células de teste L1C4, L2C5, L3C5

### 11. Patches

Servidor: `python building_db/exportar_patches.py` opção 2. Sem servidor: `exportar_patches_drive.py`
e `tfrecord_para_npz.py` no Colab (`rodadas_treino.md`, "Patches pelo Google Drive").

- [x] `modelo_A/patches_<LxCy>.npz` e `modelo_B1/...` em `/srv/almacen/db_images/patches_cnn` (.250)
- [ ] Conferir um patch: canal `swir_estimado` = 1 onde `origem_patch` ≥ 2

### 12. Rodadas (nomes com `_shm`)

Nome: `<modelo>_<backbone curto>_<patch>_shm` (`_shm_teste`, `_shm_1fold`). O sufixo vai
para os assets no GEE e define a propriedade `mosaico` dos mapas.

- [x] As 10 rodadas pelo menu de lote do `training/treinar_cnn.py` (A e B1 × 5 backbones),
      2026-11-05 07:42–09:41 no .250 (≈ 2 h, RTX 3070) →
      `runs_cnn/resumo_lote_20261105_0941.csv`; todas `ok`, 5 folds cada
- [ ] Tabela de comparação (`rodadas_treino.md`, seção 3, tabela do sazonal) — números nas Notas
- [ ] Escolher a rodada para o teste final (candidatas: `A_coatnet0_64_shm`, `B1_effv2s_64_shm`)
- [ ] Conferir os `historico_fold<k>.csv` das candidatas (melhor época por fold → `--epocas`
      do teste)

### 13. Teste final (uma vez)

- [ ] `<escolhido>_shm_teste` com `--epocas` = média das melhores épocas do CV
- [ ] `metricas_teste.json`; comparar com o teste do v1 (mesmas células de teste)

### 14. Mapa (`src/model_DL/regiao_completa/`, no .250)

Pedaços de canais do harmonizado em pasta própria: `tiles_cnn_shm/` (rodadas `_shm`);
`tiles_cnn/` fica com o v1. `baixar_tiles` e `inferir_cnn` pulam o que já existe, por isso
não misturar. Drive: `tiles_cnn_gee_shm/` e `prospectra/tiles_cnn_shm/`.

```bash
cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL/regiao_completa
# 1. canais da área toda (GEE → GeoTIFF; conda py_gee) — só o conjunto da rodada escolhida
python building_db/baixar_tiles.py --modelo A          # → tiles_cnn_shm/canais_A/<LxCy>/
python building_db/baixar_tiles.py --modelo B1         # se a escolhida for B1
# 2. prévia rápida (90 m) numa célula, depois a área toda a 30 m (conda py_DL, GPU)
#    sem argumentos abre menu: python predict/inferir_cnn.py
python predict/inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/<escolhido>_shm --tiles padrao --passo 3 --celulas L2C3
python predict/inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/<escolhido>_shm --tiles padrao --passo 1
# 3. bordas entre pedaços/células
python predict/preencher_bordas.py --mapas /srv/almacen/db_images/runs_cnn/<escolhido>_shm/mapas_passo1 --preencher --interpolar 10
# 4. GEE (mosaico = harmonizado pelo _shm do nome)
python predict/subir_predicao_gee.py --run /srv/almacen/db_images/runs_cnn/<escolhido>_shm --passo 1
#    ou: juntar_mapas_gee.py → upload à mão → propriedades_mapas_gee.py
```

- [ ] Rodada escolhida: com `<escolhido>_shm_teste/final.pt` o `inferir_cnn` usa o modelo final;
      sem ele, a média dos 5 `fold<k>.pt` da rodada de CV
- [ ] `baixar_tiles` do conjunto da escolhida (A ou B1); log diz `tiles_cnn_shm`
- [ ] `inferir_cnn` imprime `Pedaços de canais: .../tiles_cnn_shm/canais_<A|B1>`
- [ ] Prévia `--passo 3` na L2C3 conferida (`predict/app_visualizar.py`, fundo = harmonizado)
- [ ] `--passo 1` nas 25 células → `mapas_passo1/p_<LxCy>_<rodada>.tif`
- [ ] `preencher_bordas` → `mapas_passo1_preenchido/`
- [ ] Subida → coleção `cnn_favorabilidade` com `mosaico = aster-mosaico-final_sazonal_hm_razao`
- [ ] `model_DL/regiao_completa/js/visualizar_cnn_favorabilidade.js`: hoje o fundo é o mosaico v1; trocar para o
      harmonizado se for comparar

### 15. Relatório de deep learning

- [ ] v1 × sazonal: métricas CV e teste, mapas, limitação do vazamento (passo 5)

---

## Notas

- 2026-11-02: scripts apontados para o mosaico sazonal (`models_ML_fosfato.MOSAICOS`,
  `exportar_favorabilidade_ml`, `baixar_amostras_maxent`, `avaliar_maxent.ipynb`,
  `maxent_alvos.js`, `visualizar_favorabilidade_ml.js`, `exportar_patches.py`;
  `_shm` em `juntar_mapas_gee`, `subir_predicao_gee`, `propriedades_mapas_gee`).
- 2026-11-03: área reduzida (7b) rodada no harmonizado. `selecao` nova com 12 variáveis
  (v1: 13): dist_dobra_km, lito_nova_america, MAG_nT, eU_ppm, CT, B09, eTh_ppm, U_K,
  mag_thg, phosphorite_rbd, ferric_ox_4_3, mag_res. Ficam do v1: dist_dobra_km,
  lito_nova_america, MAG_nT, eU_ppm, CT, eTh_ppm, mag_thg; saem gossan_4_2, B10,
  quartz_14_12, dist_falha_km, sultan_g, ferric_iron_2_1; entram B09, U_K,
  phosphorite_rbd, ferric_ox_4_3, mag_res. Ablação, 11 ocorrências, AUC por ponto
  harmonizado × v1: RF 0,986 × 0,990, GTB 0,983 × 0,985; ap_ganho 68,9 / 68,6 (v1 133 / 105);
  Boyce 0,46 / 0,52 (v1 0,59 / 0,71). Só ASTER: RF 0,914, GTB 0,899.
- 2026-11-03: passos 2–7 (área toda) ainda sem saída no GEE: não existe
  `maxent_amostras_aster-mosaico-final_sazonal_hm_razao_bandas_indices` e
  `favorabilidade_ml` só tem os mapas do v1. A limpeza do `src/tools/deletfilesFOLDERs.py`
  (grupos A + B, 967 assets) ainda não foi rodada.
- 2026-11-03: opção 1 rodada; a opção 3 (Maxent no GEE) também gerou a task
  `maxent_fosfato_aster-mosaico-final_sazonal_hm_razao` (SUCCEEDED 03:10 UTC).
- 2026-11-03: `selecao_variaveis_RC.ipynb` rodado (fosfato, 37 ocorrências, fundo cluster).
  `selecao` (11): B02, B10, B14, SAVI, ferric_iron_2_1, B3N, phosphorite_rbd,
  quartz_14_12, glcm_contrast, glcm_hom, glcm_ent — nenhuma variável de SWIR. Ablação,
  célula de fora, AUC por ponto RF / GTB: `vnir_tir` 0,735 / 0,750; `bandas` 0,689 / 0,648;
  `indices` 0,616 / 0,577; `bandas_indices` 0,690 / 0,635; `todas` 0,701 / 0,692;
  `sem_swir` 0,756 / 0,764; `fixas_v1` 0,683 / 0,682; `selecao` 0,758 / 0,764
  (Boyce 0,99 / 0,93, top10 0,32 / 0,35, ap_ganho 2,35 / 3,03). Variáveis de SWIR pioram
  no harmonizado. Validação R2 (treino células 8 e 13, teste 4, 10, 15): `selecao` RF 0,607 /
  GTB 0,630 (Boyce 0,47 / 0,66); `fixas_v1` 0,545 / 0,539. Pontos: percentil mediano
  RF / GTB nas ocorrências de validação 90 / 92, nas de treino 65 / 61; célula 8 é a pior
  (56 / 52); 12 (RF) e 13 (GTB) das 37 no top 10 %. Referência v1 (fosfato, avaliar_maxent):
  RF 0,83, GTB 0,84, mas com fundo aleatório; não é comparável direto com o fundo cluster.
- 2026-11-05: CNN no harmonizado, 10 rodadas `_shm` (patches `P_POS = 0.94`, negativos por
  cluster). Fora do fold, 13 ocorrências de treino × fundo aleatório (auc_ponto / ap_ganho /
  Boyce / top10):
  `A_coatnet0_64` 0,968 / 16,6 / 0,59 / 0,92 · `B1_effv2s_64` 0,958 / 14,1 / 0,61 / 0,85 ·
  `A_pequena_32` 0,947 / 9,2 / 0,54 / 1,00 · `B1_regnety_64` 0,930 / 7,9 / 0,71 / 0,85 ·
  `B1_coatnet0_64` 0,915 / 7,6 / 0,59 / 0,69 · `A_regnety_64` 0,910 / 6,5 / 0,69 / 0,54 ·
  `B1_pequena_32` 0,903 / 5,4 / 0,66 / 0,54 · `A_resnet50_64` 0,898 / 4,8 / 0,66 / 0,31 ·
  `B1_resnet50_64` 0,884 / 4,3 / 0,69 / 0,38 · `A_effv2s_64` 0,877 / 6,5 / 0,74 / 0,69.
  Todas acima da referência GTB (0,84), mas a comparação é otimista: os pseudo-rótulos vêm
  de mapas RF/GTB treinados com todos os pontos (passo 5) e só 13 ocorrências entram na
  validação. Boyce baixo (0,54–0,74) em todas: a ordenação no topo é boa, a curva
  intermediária não.
