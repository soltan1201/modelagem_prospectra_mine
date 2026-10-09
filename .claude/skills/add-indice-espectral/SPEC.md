# Spec — add-indice-espectral

## Objetivo
Adicionar um novo índice espectral/mineral (razão/expressão de bandas) ao
catálogo em `calculo_spectrals_index.js`, com sua respectiva entrada de
visualização, respeitando a convenção de nomes correta.

## Pré-condições
- Fórmula do índice e bandas ASTER envolvidas foram confirmadas com o
  usuário (ou citam literatura, como o índice Al-OH de Rowan & Mars 2003
  já referenciado no `README.md`).
- Usuário indicou (ou foi perguntado) se o índice é para a convenção
  pré-mosaico (`B01-B09`), pós-mosaico (`AST_*nm`, prefixo `IDX_`), ou
  ambas.

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| Fórmula do índice | usuário | sim |
| Nome do índice | usuário (ou sugerido seguindo convenção) | sim |
| Convenção alvo (pré/pós-mosaico) | usuário | sim |
| Parâmetros de visualização (min/max/paleta) | usuário ou inferido de índices análogos | não |

## Saídas
- Nova função/expressão adicionada em `calculo_spectrals_index.js` dentro
  da função certa (`adicionarIndicesGeo`, `calcularIndicesGeologia`, ou
  `addMineralIndices`).
- Nova entrada em `visualizacao_indices` com `band`, `min`, `max`,
  `colorPalette`, `description`.

## Efeitos colaterais
- Nenhum fora do arquivo `.js` editado — é biblioteca de referência, não
  chamada automaticamente por nenhum script de export.
- Se o índice for incorporado ao score de qualidade
  (`addQualityASTER_*`), isso propaga para os scripts da skill
  `build-mosaico-semestral` — tratar como mudança de pipeline, não como
  edição isolada.

## Modos de falha e como reagir
- **Pedido para "rodar" `calculo_spectrals_index.js` como está**: avisar
  que o arquivo tem trechos com variáveis não definidas (não é
  executável de ponta a ponta) antes de tentar rodar.
- **Nome de banda ASTER incerto/inventado**: parar e checar contra
  `bandas_ASTER.txt` / `readme_notas.md` antes de escrever a expressão.
- **Usuário pede para also corrigir o bug de `copyProperties` sem ponto**:
  aplicar a correção (`.copyProperties(...)`) apenas mediante pedido
  explícito — por padrão, apenas apontar o problema.

## Não-objetivos
- Não decide sozinho os limites de visualização (`min`/`max`) sem dado de
  referência — usar valores de índices análogos já presentes no
  dicionário como ponto de partida, mas sinalizar que são estimativas.
- Não corrige bugs pré-existentes no arquivo sem que o usuário peça.
