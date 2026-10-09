# Spec — build-mosaico-semestral

## Objetivo
Produzir/editar um script GEE (JS, Code Editor) que gera o mosaico de
qualidade ASTER de um semestre (`ano` + `mes`) e o exporta como asset INT16.

## Pré-condições
- O usuário informa `ano` e, implícita ou explicitamente, o semestre
  (1º = jan-jun, 2º = jul-dez).
- Repositório contém `src/building_dataset/select_black_list_save_mosaicSem.js`
  como template — se ausente, abortar e avisar (não recriar a lógica do
  zero por inferência).
- Execução real do script acontece fora do Claude Code, no
  code.earthengine.google.com — este skill não tem runtime para validar o
  resultado.

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| `ano` | usuário | sim |
| `mes` (1 ou outro) | usuário | sim |
| lista de IDs para blacklist | usuário (após inspeção visual no GEE) | não |
| lista de IDs 100% sem nuvem | usuário | não |
| ajuste de limiares | usuário (explícito) | não |

## Saídas
- Arquivo `.js` atualizado/criado com `ano`, `mes`, `blacklist`,
  `list_Cloud_zero` corretos.
- Diff que preserva toda a lógica de `radianceToTOA`,
  `mascaraNuvemASTER`, `addQualityASTER_v5_shadow_fix`,
  `converterPara16Bit` sem alteração não solicitada.
- Opcional: entradas atualizadas em `lista_imagens_aceites.txt`,
  `lista_imagens_Notaceites.txt`, `list_zero_cloud.txt`,
  `linksSemestres.txt`/`linksSemestres_corr.txt`, quando o usuário fornecer
  os dados de um export já concluído.

## Efeitos colaterais
- Nenhum fora do filesystem do repo. Nada é executado contra o projeto GEE
  `mapbiomas-arida` por este skill — apenas o usuário, manualmente, no
  Code Editor.

## Modos de falha e como reagir
- **Ano ambíguo quanto a SWIR (perto de 2008)**: perguntar ao usuário se o
  destino de export é `mosaic_aster` ou `mosaic_aster_p2008` em vez de
  assumir.
- **Pedido para "melhorar" a máscara de nuvem/sombra sem contexto**: tratar
  como mudança de parâmetro, não como bug — pedir confirmação explícita e
  mudar um limiar por vez.
- **Template ausente ou renomeado**: parar e perguntar, não reconstruir a
  lógica de quality score por conta própria (é um algoritmo calibrado
  empiricamente, não trivial de re-derivar).

## Não-objetivos
- Não roda `earthengine` CLI nem GEE Python API — este pipeline específico
  (etapa 1) é só JS de Code Editor.
- Não decide sozinho se uma cena deve entrar na blacklist — isso exige
  inspeção visual humana no mapa.
