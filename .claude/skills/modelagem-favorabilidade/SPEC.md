# Spec — modelagem-favorabilidade

## Objetivo
Manter o fluxo pontos de referência → amostras nos mosaicos finais →
avaliação local (Maxent × RF × GTB, SHAP) → mapas RF/GTB no GEE →
visualização e relatório, para fosfato e os outros quatro alvos.

## Pré-condições
- Mosaicos `ASTER_mosaico_unico_final` e `aster-mosaico-final` (25/25).
- Acervo CPRM no GEE (`acervo_cprm/ocorrencias_fosfato`,
  `recursos_minerais`) e `pontos_fosfato` antigo.
- Conta EE que escreve em `projects/mapbiomas-arida/mine`; os scripts
  novos inicializam com `project='mapbiomas-brazil'`.
- Notebook: pandas, numpy, scikit-learn, scipy, matplotlib; `shap` só na
  seção 9.

## Entradas
| Campo | Onde |
|---|---|
| alvo(s) | menu dos scripts; `ALVOS` (dict True/False) no export ML |
| mosaico(s) | menu: único, v1 ou os dois |
| variáveis | `FEATURES` / `FEATURES_DL` |
| parâmetros dos modelos | `MAXENT`, `MODELOS` (GEE); `RF_PARAMS`, `GTB_PARAMS` (notebook) |
| fundo | `FATOR_FUNDO` (GEE); `FUNDO` aleatório/target_group (notebook) |

## Saídas
- Assets: `pontos_referencia`, `maxent_amostras_<mosaico>_bandas_indices`,
  `maxent_<alvo>_<mosaico>` (Maxent), `favorabilidade_ml/*`.
- Locais: CSVs em `src/Dados/maxent/` e `src/Dados/referencias/`;
  `src/modeling/models/` (validação, ranking SHAP, figuras); registros
  `registro_favorabilidade_ml_<data>.csv`.
- Relatório LaTeX/PDF em `docs/relatorio_avaliacao_modelos/`.

## Efeitos colaterais
Tasks `Export.table/image.toAsset` reais — só o usuário dispara. Scripts
pulam assets existentes e tasks na fila.

## Modos de falha
- `getInfo` > 5.000 elementos: filtrar no servidor antes (ex.: `na_area`).
- Amostragem interativa pesada ("Too many concurrent aggregations"): salvar
  amostras como asset e ler a tabela.
- GTB com muitas árvores no GEE: lento/falha — manter ≤ 90.
- `PROBABILITY` não aceito por algum classificador: avisar o usuário.
- Nome de feature fora do mosaico: o export para antes das tasks.

## Não-objetivos
- Não decide mosaico final, fundo nem limiares de pseudo-rótulo sem o usuário.
- Não roda scripts, notebooks nem tasks.
