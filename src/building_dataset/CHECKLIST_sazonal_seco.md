# Checklist — mosaico sazonal seco, harmonização e SWIR

Marque `[x]` ao concluir. Anote data, conta e observações na linha (ex.: task id, células
que falharam). Detalhes de cada comando: `README_guia_sazonal_seco.md`.

Legenda: `[ ]` a fazer · `[~]` rodando / em andamento · `[x]` feito · `[-]` pulado

## 1. Mosaicos originais (`build_mosaico_sazonal_seco.py`, superconta)

- [x] Período A exportado — 25/25 células
- [x] Período B exportado — 25/25 (L2C1 na 2ª task `5CKEYC7KPTP2ZPNZ36WMSWQL`)
  - [x] Rodar de novo: `só B`, célula `6` (2026-10-29 11:24)
  - [x] Task concluída — asset existe (registro 2026-10-31 20:27)

## 2. Cópia brazil → arida (`transferir_assets_brazil_para_arida.py`)

- [x] 1ª cópia — 49 imagens (2026-10-29, `registro_transferencia_copia_20261029_0748.csv`)
- [x] Copiar B L2C1 depois do item 1 (etapa `2`, copia só o que falta)
- [x] `conferir` mostra 50 na origem e 50 no destino (2026-11-01)
- [ ] (opcional) Apagar originais da superconta — só no fim, pede `APAGAR`

## 3. Inspeção visual (`visualizar_mosaico_sazonal_seco.js`, `FONTE = 'arida'`)

- [x] Emendas entre cenas identificadas (colcha de retalhos, faixa ciano, pontos brancos)

## 3b. Harmonização com MODIS (`harmonizar_sazonal_modis.py`, caatinga04 / projeto mapbiomas-brazil)

Piloto: etapa 1 rodada em todas as células (não só no piloto)

- [x] Etapa 1 — mapas de cenas (`ASTER_sazonal_cenas`) — 2026-10-30
  - [x] Tasks concluídas — todas as `cenas_A_*` e `cenas_B_*` salvas
  - [x] B L2C1: mapa de cenas — feito 2026-11-01 (era: falta (mosaico chegou na arida depois da etapa 1;
        2026-11-01: `ASTER_sazonal_cenas` 49/50). `harmonizar_sazonal_modis.py` etapa 1, só B, célula `6`)
- [x] Etapa 2 — pares ASTER × MODIS → Drive (`pares_modis_<P>_<ano>_<mês>`)
  - [x] Tentativas anteriores: tabela por período (out of memory), por ano (null geometry),
        getInfo por cena (> 1 h no 1º mês). Gargalo: máscara de nuvem a 15 m.
  - [x] Interromper (Ctrl+C) a etapa 2 antiga se ainda estiver rodando
  - [x] A — 38 CSVs (faltam 2000_06 e 2000_07: conferir se as tasks falharam ou se
        não há cenas nesses meses) · 139 cenas
  - [x] B — 60 CSVs · 744 cenas
  - [x] CSVs baixados (2026-10-31, 366 482 pares, ~415 por cena; só as faixas 0–3 do
        vermelho têm pixels) para `src/Dados/ASTER_sazonal_harmonizacao/ASTER_sazonal_harmonizacao/`
- [x] Etapa 3 — retas (local) → `coef_modis_A.csv` (139 cenas, 54 842 pares) / `_B.csv`
      (744 cenas, 311 640 pares) — 2026-10-31
- [ ] Etapa 5 — status
  - [ ] `confere` ≥ 0,99 em todas as células do piloto
  - [ ] Cenas com reta válida — A: __/__ · B: __/__ (prévia: A 124/139, B 658/744)
  - [ ] Ganhos medianos plausíveis e R² mediano ≥ 0,7 (prévia: ganho 1,12–1,28,
        offset −0,01 a −0,05, R² 0,86–0,95 — esperado de TOA → superfície)
- [ ] Etapa 4 — aplicar no piloto (`ASTER_sazonal_harmonizado`)
- [ ] Comparar no visualizador `src/show_mosaics/visualizar_harmonizado_sazonal.js`
      (tela dividida original × harmonizado)
  - [ ] Degrau nas emendas diminuiu (B3N, NDVI)
  - [ ] Faixa ciano da borda — resolvida? (se não: S3, máscara da borda)
  - [ ] Banda `harmonizado` = 2 (mediana) em pouca área
- [ ] 2ª passada relativa (`harmonizar_sazonal_relativo.py`): células que ainda destoam
      (2026-10-31: coluna 1, linhas 1 e 2) puxadas para as vizinhas boas (colunas 2 e 3),
      com o MCD43A4 mediano seco como ponte
  - [ ] Etapa 1 — task `rel_<P>` → `Drive/ASTER_sazonal_harmonizacao_rel/rel_<P>.csv`
  - [ ] CSVs em `src/Dados/ASTER_sazonal_harmonizacao/relativo/`
  - [ ] Etapa 2 — `coef_rel_<P>.csv` (cenas alvo válidas: A __/__ · B __/__)
  - [ ] Etapa 3 — `ASTER_sazonal_harm_rel` (tasks `rel_<P>_<LxCy>`)
  - [ ] Visualizador com `USAR_REL = true`: coluna 1 parecida com as colunas 2 e 3?
  - 2026-10-31: A L1C1/L2C1 salvos; só 7/24 cenas alvo com reta válida (B3N R² ≈ 0
    nas cenas de 2000) — maior parte caiu na mediana. Substituída pela 3c.

## 3c. Harmonização contra o Landsat (`README_harmonizacao_landsat.md`)

- [x] `exportar_landsat_sazonal.py` — `Landsat_sazonal_<P>_seco_mediana_<LxCy>` (A 25/25, B pronto) — 2026-10-31
- [x] `harmonizar_sazonal_landsat.py` etapa 1 — amostras A 25/25, B 24/25 (L2C1 sem mapa de cenas) — 2026-10-31
  - [x] B L2C1 (mosaico e Landsat existem; parou no mapa de cenas, item 3b etapa 1):
    - [x] Mapa de cenas: `harmonizar_sazonal_modis.py` etapa 1, só B, célula `6` → task `cenas_B_L2C1`
    - [x] Amostras: `harmonizar_sazonal_landsat.py` etapa 1, só B, célula `6`
    - [x] Curvas: etapa 2, só B (refaz `curvas_B.json` com 25 células) — 2026-11-01;
          TIR ref B: DN13 = 1057,9 + 22,56·st (antes 1061,2 + 22,53·st, ≈ 2 DN de diferença)
    - [x] Aplicar: etapa 3, só B, célula `6`, `razao` e depois `residuo` — 2026-11-01
    - [-] (opcional, pulado: mudança ≈ 2 DN no B13, bem abaixo de `DIF_NUVEM_B13` = 60) Referência B mudou um pouco: para tudo igual, apagar as 24 B de
          `hm_razao`/`hm_residuo` e refazer a etapa 3 só B
- [x] Etapa 2 — `curvas_A.json` / `curvas_B.json` — 2026-10-31
  - A: 260 pedaços com curva própria, 22 pequenos (curva da célula); TIR DN13 = 1071,2 + 21,92·st
  - B: 609 pedaços com curva própria, 169 pequenos; TIR DN13 = 1061,2 + 22,53·st
- [x] Etapa 4 — relatório; teste TIR × SWIR — 2026-10-31
  - Ganho p50 A: B01 0,69 · B02 0,80 · B3N 1,04 · B04 0,94 · B13 1,00 (padrão TOA → superfície)
  - Ganho p50 B: B01 0,69 · B02 0,85 · B3N 1,09 · B13 1,02 (faixa p10–p90 mais larga que no A)
  - ρ desloc. B13 × ganho_B04 (A) = +0,58 ✔ · × ganho_B3N (B) = +0,49 · × d_evi2 (B) = +0,41
  - Pedaços sem st (TIR fica original): A 11 · B 83 — candidato a TIR previsto pelo SWIR/VNIR
  - d_evi2 fora do padrão: L4C5 e L5C5 (≈ −0,11 nos dois períodos), A L1C1 (+0,01)
- [x] Etapa 3 — `ASTER_sazonal_hm_razao` (A 25, B 25) — 2026-11-01
- [x] Etapa 3 — `ASTER_sazonal_hm_residuo` (A 25, B 25) — 2026-11-01
- [ ] Visualizador (seletor "Harmonização" no painel, Landsat à esquerda):
  - [x] Emendas sumiram no VNIR / EVI2 (A, 2026-11-01: "ficou perfeito")
  - [ ] Al-OH e carbonato contínuos nas emendas; anomalias conhecidas mantidas
  - [ ] `hm = 2` em pouca área

- [x] **Decisão (2026-11-01):** usar o harmonizado Landsat. Pontos e RF refeitos no A
      corrigido. Fonte padrão `hm_razao` (`FONTE` em `pontos_sazonal_seco_pareados.py`);
      trocar para `hm_residuo` se o Al-OH/carbonato mostrar que ele é melhor
- [x] Scripts 4–7 apontados para o harmonizado (sufixo `_hm_razao` em pontos, pasta local,
      modelos, SWIR estimado e mosaico final)

## 4. Pontos pareados (`pontos_sazonal_seco_pareados.py`, caatinga04, `FONTE = 'hm_razao'`)

- [x] Etapa 1 — `points_sazonal_seco_hm_razao_A` (A tem as 25 células) — 50 000 pontos em 25 células
- [x] Etapa 2 — `points_sazonal_seco_hm_razao_B` (B L2C1 harmonizado pronto desde 2026-11-01)
- [x] Status — pontos por célula conferidos
- [x] Baixar (2026-11-01 17:18) — `src/Dados/ASTER_Samples_sazonal_seco_hm_razao/` (`points_sazonal_seco_A.csv`,
      `_B.csv`, `_pareados.csv`)

## 5. Notebook (`modeling/train_swir_models_sazonal.ipynb`, `SUFIXO_FONTE = '_hm_razao'`)

- [x] Seções 2–7 — modelos no A, `models/sazonal_seco_hm_razao/best_params.json` (2026-11-01;
      R² no A: RF 0,88 · XGB 0,80–0,87 · PLS 0,79–0,84)
- [x] Seção 8 — teste temporal (`r2_B`): RF 0,30–0,39 · XGB 0,32–0,40 · PLS 0,25–0,34.
      Teto (VNIR/TIR do B × A nos mesmos pontos): R² 0,09–0,30 — paisagem mudou entre os períodos
- [ ] Seção 9 — modelo pareado, `best_params_pareado.json` (não rodada; sem ela os grupos
      `sazonal_pareado` do passo 6 falham)
- [ ] Seção 10 — normalização (global × célula) (não rodada; passo 6 usou célula sem comparar)
- [~] **Decisão:** ☐ `sazonal_pareado` ☒ `sazonal_seco` + norm. célula (provisório, 2026-11-01) ☐ os dois

## 6. SWIR estimado (`modeling/estimar_swir_rf_v2.py`)

- [x] Dry-run (confere "Fonte dos grupos sazonais: hm_razao")
- [~] Exportar — 25 células `ASTER_MOSAIC_sazonal_seco_hm_razao_2009_2020_<LxCy>` em `ASTER_MOSAIC_v2`
      — 2026-11-01 19:37, sazonal_seco, treino global, normalização por célula (25/25 com reta própria),
      `registro_swir_rf_20261101_1937.csv`; tasks enviadas no projeto da conta;
      as que faltavam foram canceladas → reenviar com `PROJETO_EE = 'mapbiomas-brazil'` (caatinga04)
- [ ] Conferir no mapa

## 7. Mosaico final (`build_mosaico_final_v2.py`, fonte `sazonal_hm`)

A por cima (SWIR medido); B onde falta A (origem 2) ou onde o A parece nuvem (origem 3):
A mais claro no B01 (> `DIF_NUVEM_B01` = 50) **e** mais frio no B13 (> `DIF_NUVEM_B13` = 60 DN).

- [x] Dry-run — 25 células com camadas 1 e 2
- [x] Exportar — `aster-mosaico-final_sazonal_hm_razao/aster-mosaico-final-sazonal_hm_razao_<LxCy>`
      (2026-11-02). Bug: `swir_estimado` saiu 1 só na origem 2 (usava `SWIR_ESTIMADO = {2, 4}`
      do v1); a origem 3 também é SWIR estimado. Corrigido no script para as próximas
      exportações; o DL recalcula a banda pela `origem` (`exportar_patches.imagem_canais`),
      então não precisa reexportar
- [ ] Conferir: banda `origem` = 3 cobre as nuvens do A? solo claro / minas ficaram com o A?
      Ajustar `DIF_NUVEM_*` se trocou demais ou de menos
- [x] Apontar os próximos passos (índices, favorabilidade) para a coleção nova —
      continua em `src/modeling/regiao_completa/CHECKLIST_modelagem_sazonal.md`

## Notas

- 2026-10-29: B L2C1 — task `KT3QISZAL42WM6F56FTCB5XI` iniciada, asset não existe.
- 2026-10-30: etapa 1 da harmonização — todos os mapas de cenas A e B salvos.
- 2026-11-01: B L2C1 — mosaico existe (brazil e arida, 2ª task `5CKEYC7KPTP2ZPNZ36WMSWQL`) e Landsat
  existe; falta o mapa de cenas, por isso as etapas Landsat pularam a célula.
