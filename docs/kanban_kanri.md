# Quadro Kanri — Prospectra Irecê

Atualizado em 2026-10-17. Cada item é um cartão: o título em negrito vai no
cartão, o texto depois do travessão vai na descrição. Etiquetas sugeridas entre
colchetes: [mosaico], [swir], [acervo], [maxent], [incerteza], [infra].

---

## Tarefas a fazer

### Urgente (destrava o que está em processo)

- **Cancelar tasks duplicadas em solkanCengine** [infra] — `ASTER_base_sem_swir`
  L4C2–L4C5 e L5C1–L5C5 estão `READY`, mas esses assets já existem; vão falhar e
  ocupam a fila na frente/ao lado de L3C3–C5. Manter só L3C2–L3C5.
- **Mosaico único etapa 2 para L3C2–L3C5** [mosaico] — rodar
  `build_mosaico_unico.py` etapa 2 só depois que os bases sem_swir terminarem.
- **Conferir os 25 `sem_swir_estimado`** [mosaico] — cobertura, banda
  `normalizacao`, emendas entre células.

### Decisões

- **Exportar etapa 3 do mosaico único** [mosaico] — `build_mosaico_unico.py`
  etapa 3, retas `mosaicos`; depois das bases L3C2–L3C5. Conferir
  `retas_mosaico_unico_*.csv` e o viés no `visualizar_mosaico_unico.js`.
- **Viés residual do SWIR (~−9 %)** [swir] — testar RMA ou retas só com
  camada 1 na etapa 3.
- **Cadeia v2 seco (jul–out)** [mosaico] — exportar `ASTER_v5`,
  `ASTER_0cloud_v2`, `ASTER_MOSAIC_v2`, `aster-mosaico-final_v2` ou abandonar em
  favor do mosaico único.
- **Malha do modelo** [acervo] — como juntar covariáveis de ~500 m (gama 1042)
  com ASTER 30/90 m.
- **Extensão do modelo** [acervo] — área toda ou só a coberta pelo 1042.
- **Amostras R3 (P₂O₅ ≥ 5%)** [maxent] — treino ou só validação.
- **Fundo do Maxent** [maxent] — aleatório, restrito aos carbonatos ou
  target-group (R5).

### Modelo de SWIR

- **Época como feature** [swir] — mês/dia do ano (seno, cosseno) da cena.
- **Normalização + mesma estação no treino** [swir] — retas por cena só com pares
  de estação próxima.
- **Sem TIR / com razões** [swir] — medir quanto o TIR ajuda após a normalização.
- **Validação espacial do RF** [swir] — uma célula de fora; o R² 0,94–0,96 do
  split por ponto está inflado.
- **Banda `data` nos semestrais v4** [swir] — sem ela o grupo `semestral_v4`
  não é normalizado.
- **Refazer tabela 4.3 com a busca completa** [swir].

### Fase 0 — fechar a base

- **Apontar visualizador e `maxent_fosfato.py` para o mosaico escolhido**
  [maxent] — `ASSET_MOSAICO` está vazio.
- **`stack_variaveis.py`** [acervo] — pilha por célula: eU, eTh, K, eU/eTh,
  eU/K, campo magnético, geologia (Fm. Salitre), SRTM/Copernicus DEM.
- **Revisar os pontos de fosfato** [maxent] — depósito × ocorrência, erro de
  posição, unidade geológica.

### Fase 1 — hipóteses

- **`docs/hipoteses_fosfato.md`** [incerteza] — H1 fosforito em carbonato (B08/B06,
  B13/B14, eU alto); H2 supergênico/cárstico; H3 controle estrutural. Variáveis
  e restrições de cada uma; `HIPOTESES` no Maxent.

### Fase 2 — ensemble

- **Maxent: amostras com o mosaico final** [maxent] — `FEATURES_MODO =
  'bandas+indices'`; o teste com zero-cloud 2000–2008 ficou sem presenças.
- **Maxent: validação por célula de fora** [maxent] — AUC por célula.
- **Maxent: mapa de probabilidade** [maxent].
- **Modo ensemble no Maxent** [incerteza] — N = 100–500 realizações: bootstrap e
  jitter das presenças, fundo, variáveis, `betaMultiplier`; saída
  `fosfato_ensemble` (p_media, p_dp, p_p10, p_p90, entropia, hipotese).
- **K sementes no `estimar_swir_rf.py`** [incerteza] — K versões do SWIR e dos
  índices entram no ensemble.

### Fase 3 — confiança e risco

- **`risco_favorabilidade.py`** [incerteza] — mapa de confiança (P × σ), curva
  área-favorabilidade P10–P90, fração da favorabilidade alta vinda de SWIR
  estimado (`swir_estimado`, `origem`).

### Fase 4 — campo

- **`selecionar_alvos_campo.py`** [incerteza] — 10–20 locais por rodada com
  P ≈ 0,5 e σ alto, sem vizinhança repetida.
- **Rodada de campo 1** [incerteza] — litologia, pXRF de P₂O₅ e U, GNSS.
- **Assimilar rodada 1** [incerteza] — presenças e ausências verdadeiras;
  refazer ensemble; medir queda do σ.

### Fase 5 — priorização

- **Relatório de alvos** [incerteza] — P10 alto, SWIR medido, ≥ 2 hipóteses;
  ficha por alvo.

### Acervo — outros alvos (quando houver tempo)

- **`ocorrencias_<alvo>`** [acervo] — calcário/dolomito, Pb–Zn–F–Ba,
  talco/magnesita, ouro, com tabela de sinônimos.
- **Geoquímica multielementar no GEE** [acervo] — CaO, MgO, SiO₂, Al₂O₃, Pb, Zn,
  F, Ba, Mn, Au, Cr, Ni.
- **Avaliar aerolevantamentos 1023, 1027, 1115, 3067 e FIOL** [acervo] —
  cobertura gama fora do 1042.

### Manutenção

- **Corrigir `.copyProperties` em `mascaraNuvemASTER`** [infra] —
  `select_black_list_save_mosaicSem.js` e
  `scripts_mosaico_ASTER_samples_Landsat.js`.
- **Atualizar ou aposentar `train_swir_models.py`** [infra] — usa CSV antigo
  (`samples_v2`).
- **Definir o que versionar dos `registro_*.csv`** [infra] — hoje fora do git.

---

## Tarefas em processo

- **Mosaico único etapa 1 — sem_swir L3C2–L3C5** [mosaico] — 21/25 salvos; L3C2
  `RUNNING`, L3C3–L3C5 `READY` em solkanCengine.
- **Mosaico único etapa 2 — sem_swir_estimado** [mosaico] — 21/25 (L3C1
  concluído em 2026-10-16; faltam L3C2–L3C5).
- **Maxent fosfato** [maxent] — script pronto (amostras, validação, mapa);
  aguarda o mosaico final escolhido.

---

## Tarefas feitas

- **Grid 5 × 5 com solape e área de estudo** [mosaico] — `grid_irece_5x5_solape100`,
  `data/area_estudo.shp`.
- **Mosaicos semestrais v4** [mosaico] — `ASTER_v4`, 548 assets.
- **Zero-cloud 2000–2008 e 2009–2020** [mosaico] — `ASTER_0cloud`, 42 assets.
- **Listas de cenas** [mosaico] — blacklist e zero-cloud revisadas.
- **Piloto de saturação VNIR × Landsat** [mosaico] — saturados nas zero-cloud são
  nuvem; marcados com `saturado` e excluídos.
- **Análises de histograma e normalização** [swir] — DN bom × ruim, LUT de
  casamento.
- **Pontos pareados zero-cloud** [swir] — 00–08 × 09–20.
- **Notebook de treino PLS × XGBoost × RF** [swir] — RF escolhido; teste temporal,
  mesma estação, TIR anomalia.
- **Normalização OLS por cena** [swir] — R² temporal −0,44 → +0,49.
- **SWIR estimado v1** [swir] — `ASTER_MOSAIC`, 50 assets.
- **Mosaico final v1** [mosaico] — `aster-mosaico-final`, 25/25, com `origem` e
  `swir_estimado`.
- **Scripts da cadeia v2 seco** [mosaico] — `_v2` prontos (quality v6, jul–out).
- **`build_mosaico_unico.py`** [mosaico] — etapas 1 e 2.
- **Etapa 3 do mosaico único (código)** [mosaico] — junta com_swir e sem_swir;
  SWIR re-estimado com retas entre os mosaicos + RF (2026-10-17).
- **Mosaico único etapa 1 — com_swir** [mosaico] — 25/25.
- **Acervo CPRM: catálogo, download, .grd → GeoTIFF** [acervo].
- **Acervo CPRM no GEE** [acervo] — 10 assets em `acervo_cprm/` (gama, mag,
  litologia, estruturas, dolinas, ocorrências, P₂O₅, recursos minerais).
- **Correção da leitura Morro do Chapéu** [acervo] — 4.077 ocorrências, 42 de
  fosfato.
- **Transferência de assets mapbiomas-brazil → mapbiomas-arida** [infra].
- **Ferramenta de cota de assets** [infra] — `listar_pastas_tabelas_gee.py`.
- **Análise dos documentos de incerteza e plano Irecê** [incerteza] —
  `docs/analise_fluxo_prospeccao_e_plano_irece.md`.
- **Skills e specs do Claude Code** [infra] — 10 skills em `.claude/skills/`.
