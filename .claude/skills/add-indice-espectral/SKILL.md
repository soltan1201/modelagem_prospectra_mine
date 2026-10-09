---
name: add-indice-espectral
description: Use when the user wants to add a new spectral/mineral index (band ratio or expression) to the ASTER feature library, register its visualization params, or reconcile index naming between the pre-mosaic and post-mosaic band conventions. Triggers on "novo índice espectral", "índice mineral", "calculo_spectrals_index", "visualizacao_indices", "IDX_".
---

# Add Índice Espectral

Ajuda a adicionar um novo índice espectral/mineral (razão de bandas ou
expressão) de forma consistente com o resto do projeto.

Ver [SPEC.md](SPEC.md) para o contrato formal.

## Arquivo de referência

[src/features_process/calculo_spectrals_index.js](../../../src/features_process/calculo_spectrals_index.js)

**Atenção**: este arquivo é uma **biblioteca de referência/snippets**, não
um script executável de ponta a ponta. As linhas finais (~386-445) chamam
funções e usam variáveis nunca definidas ali (`mosaico_final`, `imgClean`,
`imgComMineral`, `imgComTextura`) — são exemplos de uso colados de outro
contexto, não código que roda como está. Trate como catálogo de funções a
reutilizar, não como pipeline pronto.

## Duas convenções de nomes coexistem — escolha a certa

1. **Pré-mosaico, bandas cruas** (`B01`...`B09`): funções
   `adicionarIndicesGeo` / `adicionarIndices` / `addIndicesASTER`. Nomes
   sem prefixo: `NDVI`, `AlOH_Clay`, `Carbonate`, `Ferric_Fe`, `SAVI`,
   `NDWI`.
2. **Pós-mosaico, bandas renomeadas** (`AST_Green_556nm`,
   `AST_SWIR2_2167nm_AlOH`, etc. — ver dicionário `bandas_nomes_novos` nos
   scripts de amostragem): função `calcularIndicesGeologia`. Nomes com
   prefixo `IDX_`: `IDX_AlOH_Clay`, `IDX_Carbonate`, `IDX_Ferric_Fe`,
   `IDX_NDVI`, `IDX_SAVI`.

Há ainda `addMineralIndices` com uma terceira família de nomes em
minúsculas (`ferric_iron_2_1`, `sultan_r`, `abrams_g`, `kaolinite_7_5`,
etc.) — composições clássicas de sensoriamento mineral (Sultan, Abrams,
razões de banda TIR). Ao adicionar um índice novo, pergunte ao usuário em
qual convenção ele deve entrar (pode ser mais de uma) antes de escrever
código.

## Checklist ao adicionar um índice

1. Confirmar a fórmula com o usuário e a(s) banda(s) ASTER envolvidas
   contra [bandas_ASTER.txt](../../../bandas_ASTER.txt) e a tabela em
   [readme_notas.md](../../../src/building_dataset/readme_notas.md) — não
   inventar comprimento de onda/banda.
2. Adicionar a expressão na função certa (`img.expression(...)` seguindo o
   padrão já usado, com `.rename('<Nome>')`).
3. Adicionar entrada correspondente no dicionário `visualizacao_indices`
   (linha ~3): `band`, `min`, `max`, `colorPalette` (array de cores GEE
   válidas), `description` (PT-BR, curta), `gamma` opcional.
4. Se o índice for aplicado antes do quality mosaic (ex.: como parte do
   score de qualidade), avisar que isso muda `addQualityASTER_*` nos
   scripts da skill `build-mosaico-semestral` — mudança de comportamento
   do pipeline, não deste arquivo isolado; peça confirmação antes.

## O que este skill NÃO faz

Não corrige silenciosamente o bug de sintaxe presente em duas cópias da
função `mascaraNuvemASTER` (em
`select_black_list_save_mosaicSem.js` e
`scripts_mosaico_ASTER_samples_Landsat.js`): falta um `.` antes de
`copyProperties(...)` no final da função
(`.updateMask(clear)\n    copyProperties(...)` em vez de
`.updateMask(clear)\n    .copyProperties(...)`), o que quebra a execução
no Code Editor. Isso deve ser **sinalizado ao usuário como um bug real**
antes de qualquer correção, não tratado como estilo do código nem
corrigido sem aviso.
