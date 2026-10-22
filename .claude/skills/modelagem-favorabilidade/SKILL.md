---
name: modelagem-favorabilidade
description: Use when the user works on mineral favorability modelling for Irecê — reference points (pontos_referencia), Maxent/RF/GTB samples and maps, the evaluation notebook (spatial validation, nulls, SHAP variable ranking), exporting favorabilidade_ml in GEE, its viewer, or the LaTeX report. Triggers on "maxent", "favorabilidade", "pontos_referencia", "avaliar_maxent", "SHAP", "FEATURES_DL", "RF/GTB", "favorabilidade_ml", "relatório da avaliação".
---

# Modelagem de favorabilidade (fosfato + 4 alvos)

Ver [SPEC.md](SPEC.md) para o contrato.

**O usuário roda tudo** (scripts, notebooks, tasks GEE). Claude escreve e
edita o código e confere só a sintaxe (`ast.parse`, `node --check`);
consultas GEE só leitura (`listAssets`, `size`) para conferir estado.

## Fluxo

| # | Script (`src/modeling/` salvo indicação) | Saída em `projects/mapbiomas-arida/mine/` |
|---|---|---|
| 1 | `preparar_pontos_referencia.py` | `pontos_referencia` + `src/Dados/referencias/pontos_referencia.csv` |
| 2 | `maxent_fosfato.py` etapa 1 (mosaico: único, v1 ou os dois) | `maxent_amostras_<mosaico>_bandas_indices` |
| 3 | `baixar_amostras_maxent.py` (ou export Drive em `export_points_mosaic.js`) | CSVs em `src/Dados/maxent/` |
| 4 | `avaliar_maxent.ipynb` | métricas, figuras `models/figuras_avaliacao_maxent/c<célula>_*.png`, `models/ranking_features_shap_*.csv` |
| 5 | `exportar_favorabilidade_ml.py` | `favorabilidade_ml/<modelo>_<alvo>_<mosaico>_<LxCy>` |
| 6 | `src/show_mosaics/visualizar_favorabilidade_ml.js` | histogramas, curva área × limiar |
| — | `maxent_alvos.js` | prévia Maxent no Code Editor |
| — | `docs/relatorio_avaliacao_modelos/` | relatório LaTeX |

`baixar_presenca_maxent.py` é legado (mosaico de prioridade montado na hora,
hoje = `aster-mosaico-final`).

## Pontos (`pontos_referencia`, 1.103)

- `alvo`: fosfato 47, carbonato 327, mvt (Ba/F/Pb/Zn) 252, ouro 208,
  manganes 269. `CODIGO_OBJ` único: o do `pontos_fosfato` antigo quando
  existe, senão `900000 + n`.
- fosfato: `fonte` = `antigo+cprm` (os 30 antigos estão todos na CPRM a
  ≤ 66 m) ou `cprm`; `papel` = `treino` (bacia de Irecê) ou `validacao` (R2).
- `subst` padronizada por `SUBST_ALVO`; mesmo alvo a ≤ 500 m = um ponto.

## Amostras

Um asset por mosaico com todos os alvos: presença = 20 pixels no buffer de
500 m (`alvo`, `fonte`, `papel`, `celula`); fundo = 10 mil pixels,
`alvo = 'fundo'`, `CODIGO_OBJ = -1`. 36 features de `calcular_features`
(14 bandas, 17 índices, `phosphorite_rbd`, 4 GLCM) + bandas de origem do
mosaico (só análise). Filtro por alvo/papel é feito na leitura.

## Variáveis — `FEATURES_DL` (21, do SHAP de 2026-10-21)

Mesma lista em dois lugares; mude os dois juntos:
`exportar_favorabilidade_ml.py` (`FEATURES`) e `maxent_alvos.js`
(`FEATURES_DL`, modo padrão). Origem: célula [34] do notebook.

## Notebook `avaliar_maxent.ipynb`

- Células numeradas na 1ª linha (`# [N]` código, `` `[N]` `` texto). O
  usuário pede mudanças pelo número; células novas recebem sufixo
  (`[13d]`, `[34a]`) para não renumerar.
- Maxent local = IWLR no sklearn (Fithian & Hastie 2013); RF e GTB com a
  mesma interface `fit(Xp, Xf)` / `predict(X)`; `MODELOS_CLS` (célula
  [13e]) alimenta 24b, 25, 26 e SHAP.
- Validação = célula da grid de fora (`oof_por_celula(..., Modelo=)`).
  Dentro da amostra só como contraste.
- Figuras com `salvar(fig, celula, nome)` → 350 dpi, prefixo `c<célula>_`.
- Edite com `nbformat` preservando as saídas que o usuário já tem.

## GEE: RF e GTB

`smileGradientTreeBoost` com no máximo ~80–90 árvores (mais fica pesado);
hoje 80, shrinkage 0,1. Parâmetros do notebook (`RF_PARAMS`, `GTB_PARAMS`)
espelham os do GEE: mude os dois juntos. Fundo sorteado = 3 × presença;
saída `probability` × 10000 Int16; treino com todos os pontos, export por
célula (grid sem solape). `p` ordena a área; 0,25 é o neutro (1:3).

## Relatório (`docs/relatorio_avaliacao_modelos/`)

`relatorio_avaliacao_modelos.tex` + `referencias.bib`; figuras em `image/`
(`\graphicspath`), nomes `c<célula>_*.png`. `lista_34a.tex` é gerado por
`gerar_lista_34a.py` a partir do CSV explicado — regenerar, não editar.
Texto colado de saída de célula precisa virar LaTeX (✔, `_`, ↑, ³⁺).
Compilar com `latexmk -pdf` num diretório temporário e copiar o PDF.

## Resultados de referência (2026-10-21)

AUC por ponto, célula de fora, fosfato: único Maxent 0,58 / RF 0,61 /
GTB 0,63; v1 0,70 / 0,83 / 0,84. Dentro da amostra ≥ 0,99 (RF/GTB) — o
sinal é em parte local. p(outros alvos) como covariável não melhora o
fosfato. Pendentes: nulos e AUC por célula com RF/GTB; v1 × único.
