# Spec — merge-mosaico-final

## Objetivo
Editar/depurar o script que consolida os mosaicos semestrais em
`projects/mapbiomas-arida/MOSAICO_FINAL_ASTER_IRECE`.

## Pré-condições
- Todos os mosaicos semestrais relevantes já foram exportados pela etapa 1
  (skill `build-mosaico-semestral`) e existem como assets em
  `projects/mapbiomas-arida/mosaic_aster` e/ou `.../mosaic_aster_p2008`.
- `src/building_dataset/merge_mosaic_semestral.js` existe como referência.

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| Confirmação se deve unir `mosaic_aster` + `mosaic_aster_p2008` num único mosaico ou tratá-los separadamente | usuário | sim, se ambíguo |
| Ajustes de visualização (bandas RGB, min/max) | usuário | não |

## Saídas
- Script `.js` atualizado.
- Nenhuma execução real (Code Editor GEE é manual).

## Efeitos colaterais
Nenhum no filesystem além do arquivo `.js`. Export real de asset GEE só
acontece quando o usuário roda manualmente no Code Editor.

## Modos de falha e como reagir
- **Bandas divergentes entre as duas coleções**: reportar a lista de bandas
  de cada uma (`bandNames()`) antes de sugerir qualquer `.merge()` — não
  assumir compatibilidade.
- **Pedido implícito de "juntar tudo em um mosaico só"**: esclarecer se
  isso significa (a) unir as duas `ImageCollection`s antes do
  `qualityMosaic`, ou (b) manter dois mosaicos finais separados como está
  hoje. São comportamentos diferentes do pipeline.

## Não-objetivos
- Não recalcula a banda `quality` — ela já vem pronta de cada mosaico
  semestral (ver nota em readme_notas.md: "Não recalcula quality — usa a
  banda já salva em cada semestre").
