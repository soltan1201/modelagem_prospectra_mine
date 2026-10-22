# Quadro Kanri — Prospectra Irecê

Atualizado em 2026-10-22. Cada item é um cartão: o título em negrito vai no
cartão, o texto depois do travessão vai na descrição. Etiquetas sugeridas entre
colchetes: [mosaico], [swir], [acervo], [maxent], [incerteza], [infra].

---

## Tarefas a fazer

### Urgente (destrava o que está em processo)

- **Refazer contagem de cenas com SWIR** [mosaico] — `contar_imagens_semestre.py`
  corrigido (`com_swir` contava a B04 vazia de 2009–2020); rodar de novo, a
  célula [4] do `graficos_imagens_semestre.ipynb` e trocar `c4_disponibilidade_semestre.png`
  em `docs/figuras_contagem/` e no relatório do mosaico.
- **Mosaico único etapa 2 para L3C2–L3C5** [mosaico] — os bases sem_swir já
  estão 25/25; só completa `ASTER_mosaico_unico` (a etapa 3 não depende dela).
- **Conferir os 25 `sem_swir_estimado`** [mosaico] — cobertura, banda
  `normalizacao`, emendas entre células.

### Decisões

- **Mosaico de entrada da favorabilidade: v1 × único** [maxent] — fosfato:
  AUC por ponto (célula de fora) 0,83 no v1 × 0,62 no único com RF/GTB.
  Conferir se é assinatura de cena da colcha do v1 (AUC por célula, nulos).
- **Modelo da favorabilidade** [maxent] — RF/GTB ganham do Maxent em 8–9 de 10
  combinações (avaliar_maxent.ipynb, 24b); confirmar com nulos rodados com
  RF/GTB antes de fixar.
- **Viés residual do SWIR (~−9 %)** [swir] — testar RMA ou retas só com
  camada 1 na etapa 3.
- **Cadeia v2 seco (jul–out)** [mosaico] — exportar `ASTER_v5`,
  `ASTER_0cloud_v2`, `ASTER_MOSAIC_v2`, `aster-mosaico-final_v2` ou abandonar em
  favor do mosaico único.
- **Malha do modelo** [acervo] — como juntar covariáveis de ~500 m (gama 1042)
  com ASTER 30/90 m.
- **Extensão do modelo** [acervo] — área toda ou só a coberta pelo 1042.
- **Amostras R3 (P₂O₅ ≥ 5%)** [maxent] — treino ou só validação.
- **Fundo dos modelos** [maxent] — aleatório (atual), restrito aos carbonatos
  ou target-group (R5); o notebook já tem `FUNDO = 'target_group'`.

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

- **`stack_variaveis.py`** [acervo] — pilha por célula: eU, eTh, K, eU/eTh,
  eU/K, campo magnético, geologia (Fm. Salitre), SRTM/Copernicus DEM.
- **Revisar os pontos de fosfato** [maxent] — depósito × ocorrência, erro de
  posição, unidade geológica.

### Fase 1 — hipóteses

- **`docs/hipoteses_fosfato.md`** [incerteza] — H1 fosforito em carbonato (B08/B06,
  B13/B14, eU alto); H2 supergênico/cárstico; H3 controle estrutural. Variáveis
  e restrições de cada uma; `HIPOTESES` no Maxent.

### Fase 2 — ensemble

- **Nulos e AUC por célula para RF/GTB** [maxent] — `avaliar_maxent.ipynb`
  seção 5 hoje só roda o Maxent.
- **Favorabilidade dos outros alvos** [maxent] — carbonato, mvt, ouro,
  manganês: ligar em `ALVOS` do `exportar_favorabilidade_ml.py`.
- **Amostras/patches para deep learning** [maxent] — áreas com p ≥ 0,5 no RF e
  no GTB (`favorabilidade_ml`); antes: tirar manchas com borda de cena,
  conferir na litologia (Una/Salitre), rarefazer e excluir células de
  validação. Ver seção "Mapas de favorabilidade" do relatório de avaliação.
- **Conferir as manchas novas do sul e do centro** [maxent] — p alto sem pontos
  de treino nem de validação; checar no mosaico único e na geologia; candidatas
  ao campo.
- **Mapas RF/GTB no mosaico único** [maxent] — rodar `exportar_favorabilidade_ml.py`
  com o outro mosaico para comparar as manchas.
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

- **Correção da saturação com referência independente** [mosaico] — se voltar
  ao tema: usar a mesma área em cena com ganho normal (`GAIN_SETTING_B01/B02`)
  como verdade nos pixels saturados (`src/modeling/apontes.md`).

- **Corrigir `.copyProperties` em `mascaraNuvemASTER`** [infra] —
  `select_black_list_save_mosaicSem.js` e
  `scripts_mosaico_ASTER_samples_Landsat.js`.
- **Atualizar ou aposentar `train_swir_models.py`** [infra] — usa CSV antigo
  (`samples_v2`).
- **Definir o que versionar dos `registro_*.csv`** [infra] — hoje fora do git.

---

## Tarefas em processo

- **Mosaico único etapa 2 — sem_swir_estimado** [mosaico] — 21/25 (L3C1
  concluído em 2026-10-16; faltam L3C2–L3C5).
- **Relatórios LaTeX** [infra] — mosaico (`docs/relatorio_mosaico_aster/`),
  avaliação (`docs/relatorio_avaliacao_modelos/`) e acervo
  (`docs/relatorio_acervo_cprm/`) com figuras; o do mosaico ganhou a seção 5.1
  (saturação do VNIR, pseudo-saturação com Landsat). Falta completar as
  referências marcadas `[CONFERIR]` (Guha 2019, Gonzalez 2019, Chander 2013) e
  trocar o `c4` corrigido.

---

## Tarefas feitas

- **Grid 5 × 5 com solape e área de estudo** [mosaico] — `grid_irece_5x5_solape100`,
  `data/area_estudo.shp`.
- **Mosaicos semestrais v4** [mosaico] — `ASTER_v4`, 548 assets.
- **Zero-cloud 2000–2008 e 2009–2020** [mosaico] — `ASTER_0cloud`, 42 assets.
- **Listas de cenas** [mosaico] — blacklist e zero-cloud revisadas.
- **Piloto de saturação VNIR × Landsat** [mosaico] — pseudo-saturação (treino
  DN < 200, teste 200–253) em 4 cenas: correção pior que o valor cortado em 6
  de 8 combinações (melhor só na B02 de 2006, empate na B02 de 2010);
  saturados nas zero-cloud são nuvem → marcados com `saturado` e excluídos.
  Documentado no relatório do mosaico (seção 5.1).
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
- **Mosaico único etapa 1 — com_swir e sem_swir** [mosaico] — 25/25 cada.
- **Mosaico único etapa 3** [mosaico] — `ASTER_mosaico_unico_final` 25/25
  (2026-10-18/19).
- **Visualizadores dos mosaicos finais** [mosaico] —
  `visualizar_mosaico_unico_final.js` e `visualizar_comparar_mosaicos_finais.js`
  (v1 × único com cortina, 15 índices, Δ, acervo).
- **Pontos de referência unificados** [maxent] — `pontos_referencia` (1.103:
  fosfato 47 com `fonte` antigo+cprm/cprm e `papel`, carbonato 327, mvt 252,
  ouro 208, manganês 269); `preparar_pontos_referencia.py`.
- **Amostras dos dois mosaicos** [maxent] — `maxent_amostras_<mosaico>_bandas_indices`
  (31.860 cada: presença de todos os alvos + 10 mil de fundo, 36 features).
- **Maxent por alvo e mosaico** [maxent] — `maxent_fosfato.py` (menu alvo/mosaico)
  e `maxent_alvos.js` (prévia no Code Editor).
- **Avaliação Maxent × RF × GTB** [maxent] — `avaliar_maxent.ipynb`: célula de
  fora, ROC, Boyce, curva de sucesso, nulos, R2, p(outros alvos) como covariável.
- **Seleção de variáveis por SHAP** [maxent] — 21 de 36 (`FEATURES_DL`,
  `models/ranking_features_shap_20261021.csv`).
- **Export de bandas e índices para o Drive** [mosaico] —
  `features_process/exportar_mosaicos_indices.{js,py}` e versão Colab.
- **Acervo CPRM: catálogo, download, .grd → GeoTIFF** [acervo].
- **Acervo CPRM no GEE** [acervo] — 10 assets em `acervo_cprm/` (gama, mag,
  litologia, estruturas, dolinas, ocorrências, P₂O₅, recursos minerais).
- **Correção da leitura Morro do Chapéu** [acervo] — 4.077 ocorrências, 42 de
  fosfato.
- **Transferência de assets mapbiomas-brazil → mapbiomas-arida** [infra].
- **Ferramenta de cota de assets** [infra] — `listar_pastas_tabelas_gee.py`.
- **Análise dos documentos de incerteza e plano Irecê** [incerteza] —
  `docs/analise_fluxo_prospeccao_e_plano_irece.md`.
- **Favorabilidade RF/GTB do fosfato no GEE** [maxent] — `favorabilidade_ml`,
  mosaico v1, RF e GTB 25/25 cada (2026-10-22); visualizador com histogramas,
  curva área × limiar, concordância p ≥ limiar e legenda.
- **Contagem de cenas por semestre** [mosaico] — `contar_imagens_semestre.py` e
  `graficos_imagens_semestre.ipynb` (2000–2008 × 2009–2020).
- **Figuras do treino do SWIR a 450 dpi** [swir] — `train_swir_models.ipynb`
  salva `c<célula>_*.png` em `models/zero_cloud/figuras/`.
- **Mapa de localização** [infra] — `show_mosaics/mapa_localizacao.js` (biomas,
  estados, área de estudo).
- **Skills e specs do Claude Code** [infra] — 11 skills em `.claude/skills/`
  (nova: `modelagem-favorabilidade`, 2026-10-22).
