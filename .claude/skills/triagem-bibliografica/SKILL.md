---
name: triagem-bibliografica
description: Use when the user wants to filter or extend the Scopus literature-screening script by keyword category (mining/ore, deep learning, remote sensing), add new keywords, or debug why an article was or wasn't selected. Triggers on "filtrar_artigos", "triagem bibliográfica", "revisão Scopus", "palavras-chave mineração deep learning".
---

# Triagem Bibliográfica (Scopus)

Componente paralelo ao pipeline geoespacial: filtra um export Scopus para
manter apenas artigos com pelo menos 1 ocorrência em cada uma de 3
categorias de palavra-chave — mineração/prospecção, deep learning, e
sensoriamento remoto.

Ver [SPEC.md](SPEC.md) para o contrato formal.

## Arquivo de referência

[src/revision_blibiografic/filtrar_artigos.py](../../../src/revision_blibiografic/filtrar_artigos.py)

## Pré-condição: `df_scopus` não é carregado pelo script

O arquivo **assume que `df_scopus` já existe** no ambiente (é um snippet
para colar num notebook/REPL depois de um `pd.read_csv` do export Scopus,
não há esse `read_csv` no arquivo). Antes de "rodar" isso, confirme com o
usuário qual é o CSV de origem e adicione o carregamento
(`df_scopus = pd.read_csv(...)`) se for executar como script `.py`
standalone.

## Como o matching funciona (importante para não introduzir falsos positivos)

`count_keywords` faz **correspondência por substring**
(`kw in text_clean`), não por palavra inteira, depois de lowercase +
remover pontuação. Isso é intencional (`'mineralization'` casa com
`'mineral'`), mas significa que uma keyword curta pode aparecer dentro de
outras palavras sem relação — ao adicionar uma nova keyword em
`lstwordRS`/`lstwordDL`/`lstwordMine`, avalie o risco de substring
espúria (ex.: siglas de 3 letras) e avise o usuário se for o caso.

## Fluxo de trabalho

1. Confirmar/carregar `df_scopus` com as colunas esperadas: `Title`,
   `Abstract`, `Author Keywords`, `Index Keywords`.
2. Para adicionar uma keyword nova: colocar em minúsculas na lista certa
   (`lstwordMine`, `lstwordDL`, ou `lstwordRS`) — as listas já seguem essa
   convenção.
3. O filtro final exige `Mine > 0 AND DeepLearning > 0 AND RemoteSensing
   > 0` simultaneamente — se o usuário quiser relaxar (ex.: OR entre
   categorias, ou pesos diferentes), isso é mudança de critério de
   triagem, não just um ajuste de keyword; confirme antes de mudar a
   lógica de `df_filtered`.
4. Nota de ordem de execução: `import re` está no meio do arquivo (depois
   de `count_keywords`, que já usa `re`) — funciona em células de notebook
   executadas em ordem, mas quebra se rodado como script `.py` puro do
   jeito que está. Mover o import para o topo se for rodar como script.

## O que este skill NÃO faz

Não decide sozinho novas categorias de triagem (além de mineração/DL/RS)
nem pesos/critérios de corte — isso é escopo de revisão bibliográfica do
usuário/orientador, não uma decisão técnica do código.
