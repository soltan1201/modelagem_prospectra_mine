# Spec — triagem-bibliografica

## Objetivo
Rodar/estender `src/revision_blibiografic/filtrar_artigos.py`, que conta
ocorrências de palavras-chave por categoria em título/resumo/keywords de
um export Scopus e filtra artigos com hit em todas as 3 categorias.

## Pré-condições
- Um DataFrame `df_scopus` carregado (de um CSV/export Scopus) com colunas
  `Title`, `Abstract`, `Author Keywords`, `Index Keywords` — **não
  fornecido pelo script**, deve ser confirmado/adicionado antes de
  executar como script standalone.
- `pandas` disponível (confirmado no ambiente).

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| CSV/export Scopus de origem | usuário | sim |
| Novas palavras-chave a adicionar | usuário | não |
| Critério de filtro (AND entre categorias, hoje fixo) | usuário, se quiser mudar | não |

## Saídas
- `df_filtered`: subconjunto de `df_scopus` com `Mine>0 & DeepLearning>0 &
  RemoteSensing>0`.
- Print no console: total original, total filtrado, primeiros 10 títulos
  selecionados.

## Efeitos colaterais
Nenhum arquivo é escrito pelo script como está (`df_filtered` fica em
memória) — se o usuário quiser persistir o resultado, é preciso adicionar
um `to_csv`/`to_excel` explícito (confirmar caminho de saída antes).

## Modos de falha e como reagir
- **`NameError: df_scopus`**: o script foi rodado sem o load do CSV —
  perguntar o caminho do export Scopus e adicionar o `pd.read_csv`
  correspondente, não assumir um caminho.
- **Colunas ausentes** (`Title`/`Abstract`/`Author Keywords`/`Index
  Keywords`): conferir o nome exato das colunas no CSV real do Scopus
  (pode variar por tipo de export) antes de rodar `count_keywords`.
- **Keyword nova gerando muitos falsos positivos**: se o total filtrado
  disparar de forma suspeita após adicionar uma keyword curta, verificar
  substring matching antes de aceitar o resultado.

## Não-objetivos
- Não implementa exportação do resultado nem deduplicação de artigos —
  fora do escopo do script original.
- Não decide os critérios de inclusão/exclusão bibliográfica (isso é
  metodológico, do usuário/orientador da revisão).
