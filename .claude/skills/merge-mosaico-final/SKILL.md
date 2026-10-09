---
name: merge-mosaico-final
description: Use when the user wants to merge all semestral ASTER mosaic assets into the final multi-temporal mosaic (MOSAICO_FINAL_ASTER_IRECE), check band/property consistency between the mosaic_aster and mosaic_aster_p2008 collections, or adjust the merge script. Triggers on "mosaico final", "merge_mosaic_semestral", "unir mosaicos semestrais", "MOSAICO_FINAL_ASTER_IRECE".
---

# Merge Mosaico Final (ASTER)

Etapa 2 do pipeline (ver
[src/building_dataset/readme_notas.md](../../../src/building_dataset/readme_notas.md)):
combina **todos** os mosaicos semestrais já exportados (etapa 1) num único
asset multitemporal, escolhendo pixel a pixel a melhor observação via
`qualityMosaic('quality')`.

Ver [SPEC.md](SPEC.md) para o contrato formal.

## Arquivo de referência

[src/building_dataset/merge_mosaic_semestral.js](../../../src/building_dataset/merge_mosaic_semestral.js)

## Pontos de atenção específicos deste script

1. **Duas coleções de origem**: `mosaic_aster` (anos com SWIR válido) e
   `mosaic_aster_p2008` (SWIR quebrado, pós-2008). O script atual **carrega
   e imprime as duas separadamente, mas NÃO as une antes do
   `qualityMosaic`** — cada uma gera seu próprio mosaico (`mosaico` e
   `mosaico_p2008`). Se o usuário pedir "um único mosaico final" combinando
   as duas eras, isso exige unir as `ImageCollection`s (`.merge()`) — bandas
   precisam bater primeiro (ver ponto 2). Não assuma isso já está feito;
   confirme a intenção antes de mudar o comportamento do script.
2. **Antes de unir**, sempre rode e compare:
   `colecao.first().bandNames()` vs `colecao_p2.first().bandNames()`. Se as
   bandas SWIR (B04-B09) estiverem ausentes em `mosaic_aster_p2008`, uma
   união direta via `.merge()` do GEE ainda funciona (bandas ausentes viram
   mascaradas), mas avise o usuário do efeito: pixels pós-2008 não
   competirão em SWIR no quality mosaic.
3. **Fator de escala fixo em `ESCALA = 10000`** só para visualização
   (`divide(ESCALA)`) — o asset exportado continua INT16 bruto.
4. **Export usa `pyramidingPolicy` com `'.default': 'sample'` e
   `'quality': 'max'`** — não mude para `'mean'` na banda quality (perderia
   o sentido de "melhor observação" nos overviews).

## O que este skill NÃO faz

Não decide sozinho por unir as duas coleções (`mosaic_aster` +
`mosaic_aster_p2008`) num só `qualityMosaic` — isso é uma mudança de
comportamento do pipeline, não um ajuste de parâmetro; peça confirmação
explícita do usuário antes de implementar. Não executa o script (Code
Editor GEE, sem runtime local).
