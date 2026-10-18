# O que falta para o mosaico final (mosaico único ASTER)

Estado conferido em 2026-10-18 nos assets do GEE (só leitura).
Script: `src/building_dataset/build_mosaico_unico.py`.

## Onde estamos (2026-10-18, conferido nos assets)

| Item | Estado |
|---|---|
| Base `com_swir` (`mine/ASTER_mosaico_unico_base`) | 25/25 |
| Base `sem_swir` | 23/25 — **faltam L3C4 e L3C5** (tasks canceladas; nada na fila) |
| Mosaico final (`mine/ASTER_mosaico_unico_final`) | 0/25 |

## Passos

O script agora confere **nos assets** (não nos CSV de registro) o que falta e
lança tudo no rodízio caatinga01–05 + solkan1201, **2 tasks por conta**, todas
no projeto EE `mapbiomas-brazil`.

### 1. Primeira rodada — 2 bases + 23 finais

```
cd src/building_dataset
python build_mosaico_unico.py
```

| Pergunta | Resposta |
|---|---|
| O que fazer | **1** — revisar nos assets o que falta e lançar tudo |
| Modo | primeiro **1 — só mostrar** (tabela célula × com_swir / sem_swir / final) |

Rodar de novo com **Modo: 2 — exportar o que falta**:

- [ ] Tabela mostra `FALTA` só em L3C4/L3C5 (sem_swir) e nos 25 finais
- [ ] Disparou 2 bases (`ASTER_base_sem_swir_L3C4`, `_L3C5`, mesmos meses das
      outras sem_swir: ano todo) e 23 finais (as retas de cada célula são
      calculadas antes — ~25 s por célula, ~10 min)
- [ ] Gerou `retas_mosaico_unico_<data>.csv` e `registro_mosaico_unico_revisao_<data>.csv`
- [ ] Nenhuma célula com `sem_retas` / `erro_inicio` no registro

### 2. Segunda rodada — finais de L3C4 e L3C5

- [ ] Esperar `ASTER_base_sem_swir_L3C4` e `_L3C5` em `COMPLETED`
- [ ] Rodar de novo a opção 1 → exportar. Ela vê o que já está salvo ou na
      fila e lança só os 2 finais que faltam.

### 3. Conferir o resultado

- [ ] Opção 1 → só mostrar: as 25 linhas com `ok ok ok`
- [ ] Bandas: B01–B14, `quality`, `data`, `camada`, `swir_estimado`,
      `normalizacao`
- [ ] Inspeção visual (emendas entre células, SWIR estimado × medido)

## Pendências depois do mosaico

- [ ] **Visualizador**: `src/show_mosaics/visualizar_mosaico_unico.js` ainda não
      lê `ASTER_mosaico_unico_final` (adicionar a coleção e as camadas
      `swir_estimado` / `normalizacao`).
- [ ] **Viés residual do SWIR estimado** (~−9 % em L2C2, +3/+11 % em L1C1 no
      teste): decidir se aplica correção por célula pós-RF (reta SWIR estimado →
      medido nos pixels com |ΔNDVI| ≤ 0,05).
- [ ] **Maxent**: apontar `ASSET_MOSAICO` de `src/modeling/maxent_fosfato.py`
      para a nova coleção e refazer amostras / validação / mapa.

## Se algo der errado

| Sintoma | Causa provável | O que fazer |
|---|---|---|
| `sem_base` em células que não são L3C4/L3C5 | base não lida pela conta | conferir se a conta do rodízio lê `mine/ASTER_mosaico_unico_base` |
| `sem_retas` | célula sem pares com_swir × sem_swir com NDVI parecido | avisar; dá para cair nas retas `pareados` para essa célula |
| "Too many concurrent aggregations" no conferir | muitas consultas ao mesmo tempo | rodar de novo (as retas são uma consulta por vez) |
| Conta fora do rodízio (aviso) | solkanCengine / solkanGeodatin | normal — não usam o projeto mapbiomas-brazil |
