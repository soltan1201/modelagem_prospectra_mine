---
name: build-mosaico-semestral
description: Use when the user wants to build, update, or debug a semestral ASTER quality-mosaic script for the Irecê (BA) pipeline — creating a new ano/semestre script, tuning cloud/shadow thresholds, managing the image blacklist/zero-cloud lists, or fixing the quality-score function. Triggers on "mosaico semestral", "novo semestre ASTER", "blacklist de imagens", "quality score ASTER", "select_black_list_save_mosaicSem".
---

# Build Mosaico Semestral (ASTER)

Gera/ajusta o script GEE (JavaScript, Code Editor) que produz um mosaico de
qualidade para **um semestre específico** de imagens ASTER L1T sobre o bloco
de Irecê (BA), como etapa 1 do pipeline descrito em
[src/building_dataset/readme_notas.md](../../../src/building_dataset/readme_notas.md).

Ver [SPEC.md](SPEC.md) para o contrato formal (entradas/saídas/pré-condições).

## Arquivo de referência

[src/building_dataset/select_black_list_save_mosaicSem.js](../../../src/building_dataset/select_black_list_save_mosaicSem.js)
é o script canônico. Ao criar um script para um novo ano/semestre, copie este
arquivo como base — não reescreva a lógica do zero.

## Fluxo de trabalho

1. **Definir `ano` e `mes`** (linha ~33-35). `mes = 1` → semestre `'1'`
   (jan-jun); qualquer outro valor → semestre `'2'` (jul-dez). O período é
   sempre 6 meses a partir de `data_inicio`.
2. **Escolher a linha de qualidade certa pelo período**:
   - Anos com SWIR funcional (**≤ 2008**): usar `addQualityASTER_v5_shadow_fix`
     (não usa SWIR no score, mas SWIR ainda é exportado) ou a variante
     `addQualityASTER_v2`/`addQualityASTER_CC0` conforme o caso de cena
     100% limpa (`CLOUDCOVER == 0` → pula a máscara de nuvem).
   - Anos **pós-2008** (SWIR do ASTER quebrado): o destino de export muda
     para `projects/mapbiomas-arida/mosaic_aster_p2008/` (ver linha ~579)
     em vez de `mosaic_aster/`. Confirme com o usuário qual coleção de
     destino aplica-se ao ano pedido antes de exportar.
3. **Atualizar a blacklist** (`blacklist`, ~linha 39): IDs de cenas ruins
   (nuvem, artefato) vão aqui — cole o `system:index` completo
   (`ASTER/AST_L1T_003/<timestamp>`), o script extrai o sufixo sozinho.
   Cenas 100% sem nuvem entram em `list_Cloud_zero` (~linha 109), não na
   blacklist.
4. **Não altere os limiares** (`lim_visBright`, `lim_swir1`,
   `ajuste_termico`, `raio_borda_nuvem`, `altura_nuvem_estimada`,
   `lim_escuro_green`, `lim_escuro_nir`, `raio_borda_sombra`) a menos que o
   usuário peça explicitamente — eles foram calibrados empiricamente para o
   semiárido baiano. Se for ajustar, mude um de cada vez e peça para o
   usuário validar visualmente no Code Editor antes do export (este script
   só roda no Code Editor GEE — não há runtime JS local para testar).
5. **Registrar o resultado** nos três arquivos de rastreamento (o script em
   si não escreve neles — é trabalho manual do usuário, mas você pode
   preencher a entrada quando ele informar o resultado):
   - `lista_imagens_aceites.txt` / `lista_imagens_Notaceites.txt` — lista
     Python-like `ano = [...]` com os IDs aceitos/rejeitados.
   - `list_zero_cloud.txt` — mesma estrutura, chave `ano_semestre`.
   - `linksSemestres.txt` / `linksSemestres_corr.txt` — linha
     `ano_semestre = <link do Code Editor>` (o `_corr` é a versão
     corrigida/revisada do script salvo).
6. **Confirmar o nome do asset exportado**:
   `ASTER_QualityMosaic_<ano>_semestre_<N>_INT16`, escala 15 m,
   `pyramidingPolicy: {'.default': 'sample'}`.

## Limitações conhecidas (não tente "consertar" sem avisar o usuário)

Ver [aspectos_iniciais.txt](../../../src/building_dataset/aspectos_iniciais.txt)
para a crítica original: uso de radiância vs reflectância TOA, ausência de
correção atmosférica (DOS), e deslocamento espacial entre cenas de datas
diferentes (paralaxe em relevo). Essas são decisões de projeto conhecidas,
não bugs — não as "corrija" silenciosamente.

## O que este skill NÃO faz

Não executa o script (é JS de Code Editor GEE, sem runtime local). A ação de
rodar/exportar é sempre manual no navegador; este skill só edita/gera o
código e ajuda a interpretar prints/erros que o usuário colar de volta.
