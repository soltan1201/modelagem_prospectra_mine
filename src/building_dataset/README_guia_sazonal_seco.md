# Guia de comandos — mosaico sazonal seco e SWIR estimado

Roteiro de execução, na ordem, do mosaico sazonal seco ASTER (períodos A e B) até o
SWIR estimado do período B. Cada etapa diz o comando, a conta usada, o que conferir
antes de seguir e onde fica o resultado.

- **A**: 2000-03 a 2008-03, VNIR + SWIR + TIR (SWIR medido).
- **B**: 2009–2020, VNIR + TIR (sem SWIR; é estimado na etapa 6).

Todos os scripts são interativos: perguntam períodos, células e modo ao rodar.
Rode sempre o **dry-run** primeiro quando o script oferecer.

| Etapa | Script | Conta | Saída |
|---|---|---|---|
| 1 | `building_dataset/build_mosaico_sazonal_seco.py` | `superconta` | `mapbiomas-brazil/.../CAATINGA/ASTER_sazonal` |
| 2 | `tools/transferir_assets_brazil_para_arida.py` | `caatinga04` (copia) / `superconta` (apaga) | `mapbiomas-arida/mine/ASTER_sazonal` |
| 3 | `show_mosaics/visualizar_mosaico_sazonal_seco.js` | Code Editor | — |
| 4 | `building_dataset/pontos_sazonal_seco_pareados.py` | `caatinga04` | `mapbiomas-arida/mine/points_sazonal_seco_A` / `_B` e CSVs |
| 5 | `modeling/train_swir_models_sazonal.ipynb` | local | `modeling/models/sazonal_seco/` |
| 6 | `modeling/estimar_swir_rf_v2.py` | conta ativa | `mapbiomas-arida/mine/ASTER_MOSAIC_v2` |

Por que contas diferentes: a `superconta` escreve em `mapbiomas-brazil` mas **não lê
`mapbiomas-arida`**; a `caatinga04` lê `mapbiomas-brazil` e lê/escreve em
`mapbiomas-arida`. Por isso o mosaico é exportado temporariamente no brazil e copiado.

Tasks de qualquer etapa: `python tools/revisarTasksXconta.py` (filtre pelo nome, ex.:
`ASTER_sazonal_`, `points_sazonal_`, `ASTER_MOSAIC_sazonal`).

---

## Etapa 1 — Mosaicos sazonais secos (A e B)

```bash
cd src/building_dataset
python build_mosaico_sazonal_seco.py
```

1. Períodos: `A e B`, `só A` ou `só B`.
2. Células: Enter para todas, ou ids (ex.: `6` ou `1,2,7` ou `11-15`).
3. Modo: dry-run (só contagens de cenas) e depois exportar.
4. Sobrescrever: **não**, a menos que queira refazer um asset existente
   (pede `APAGAR` para confirmar; apaga o asset antigo).

O script pula o que já existe ou está processando e mostra `última task FALHOU: ...`
para os itens cuja task anterior falhou. Registro em
`registro_sazonal_seco_<data>.csv` (status, cenas, task_id por período × célula).

**Conferir antes de seguir:** 25 imagens A + 25 imagens B na origem. Uma célula
que falta é refeita rodando de novo só com aquele período e aquela célula — o script
não reexporta o que já existe.

> 2026-10-29: falta **B `L2C1`** (id 6). A task `KT3QISZAL42WM6F56FTCB5XI`
> (superconta, 2026-10-28 22:01) foi iniciada mas o asset não existe — ver o erro em
> `revisarTasksXconta.py` e rodar de novo: período `só B`, célula `6`, exportar.

## Etapa 2 — Copiar para mapbiomas-arida

```bash
cd src/tools
python transferir_assets_brazil_para_arida.py
```

`COLECOES` no topo do script deve ter só `'ASTER_sazonal': 'ASTER_sazonal'` ativa.

1. `conferir` — lista origem × destino (o que falta copiar). Não altera nada.
2. `copiar` — com a `caatinga04`; copia só as imagens que ainda não estão no
   destino (rodar de novo depois de refazer uma célula na etapa 1 copia só ela).
3. `apagar originais` — opcional, libera espaço da superconta. Só para imagens já
   conferidas no destino; pede `APAGAR`. **Não tem volta.** Não apague antes de
   terminar a etapa 1 em todas as células.

Registro em `registro_transferencia_copia_<data>.csv`.

## Etapa 3 — Ver os mosaicos (opcional, recomendado)

Colar `src/show_mosaics/visualizar_mosaico_sazonal_seco.js` no Code Editor.
`FONTE = 'arida'` (padrão) lê a cópia; `'brazil'` lê o export temporário.

O painel mostra quantas células têm A e B e quais faltam; a grid fica vermelha onde
faltam os dois e laranja onde falta um. Conferir `n_obs` baixo, emendas entre células
e a diferença B − A nas bandas comuns.

## Etapa 3b — Harmonizar com o MODIS (piloto)

Corrige as emendas entre cenas sem refazer o mosaico. Para cada cena ASTER ajusta uma
reta contra o MCD43A4 do mesmo dia (Terra: mesma hora e atmosfera) e aplica a reta
aos pixels do mosaico que vieram daquela cena.

```bash
cd src/building_dataset
python harmonizar_sazonal_modis.py      # conta caatinga04, projeto mapbiomas-brazil
```

1. Etapa `1` — mapa de cenas (por célula) → `mapbiomas-arida/mine/ASTER_sazonal_cenas`.
   Reconstrói o qualityMosaic e grava `cena_t` (cena vencedora) e `confere`
   (1 = igual ao asset salvo).
2. Etapa `2` — pares ASTER × MODIS: uma task `pares_modis_<P>_<ano>_<mês>` por mês
   (A: 40, B: 60), em batch → CSVs em `Drive/ASTER_sazonal_harmonizacao/`. Pode rodar
   junto com a 1. Rodar de novo reenvia só os meses sem task concluída ou rodando.
3. Baixe a pasta do Drive e ponha os CSVs (`A_2003_07.csv` …) em
   `src/Dados/ASTER_sazonal_harmonizacao/ASTER_sazonal_harmonizacao/`.
4. Etapa `3` — retas (local, segundos): `coef_modis_A.csv` / `_B.csv`. Pode refazer à
   vontade mudando os parâmetros de ajuste.
5. Etapa `5` — status: `confere` tem de dar ≈ 1 em todas as células (se não, a
   blacklist ou a zero-cloud mudou depois do export e o mapa não vale); ver quantas
   cenas têm reta válida e a faixa dos ganhos.
6. Etapa `4` — aplicar → `mapbiomas-arida/mine/ASTER_sazonal_harmonizado`, mesmos nomes
   dos originais, mais as bandas `cena_t` e `harmonizado` (1 = reta da cena,
   2 = mediana do período). Os originais não mudam.

Comece por 2–3 células com emenda forte, nos dois períodos. A saída fica em
reflectância de superfície (NBAR) × 1000, não TOA: **use A e B harmonizados juntos**
nas etapas 4–6, nunca um harmonizado com o outro original. O TIR não é corrigido.

## Etapa 4 — Pontos pareados

```bash
cd src/building_dataset
python pontos_sazonal_seco_pareados.py
```

Rodar as opções **uma por vez, esperando a task terminar entre elas**:

1. `etapa 1` — sorteia até 2000 pontos por célula no mosaico A (30 m, onde VNIR +
   SWIR + TIR têm dado) → asset `points_sazonal_seco_A`.
2. `etapa 2` — amostra o mosaico B nos mesmos pontos → `points_sazonal_seco_B`.
   Precisa da etapa 1 pronta. Célula sem mosaico B é pulada (os pontos do A dela ficam
   sem par). Por isso: **termine as etapas 1 e 2 do guia (todas as células B) antes**.
3. `status` — pontos por célula nos dois assets.
4. `baixar` — baixa as duas tabelas (cache por célula) e pareia pelo `pid` em
   `src/Dados/ASTER_Samples_sazonal_seco/`:
   - `points_sazonal_seco_A.csv`, `points_sazonal_seco_B.csv`
   - `points_sazonal_seco_pareados.csv` (bandas com sufixo `_A` / `_B`)

Para refazer uma etapa, apague o asset de pontos antes (o script recusa sobrescrever).

## Etapa 5 — Análise e hiperparâmetros (notebook)

Abrir `src/modeling/train_swir_models_sazonal.ipynb` e rodar tudo, em ordem. Precisa
dos CSVs da etapa 4. RF e XGBoost usam `RandomizedSearchCV`: demora (reduza
`N_ITER_RF` / `N_ITER_XGB` para um teste rápido).

| Seção | O que responde |
|---|---|
| 2 | dados, correlações, distribuição do SWIR; filtro `N_OBS_MIN` |
| 3–6 | PLS, XGBoost e RF treinados no A (VNIR + TIR → SWIR); comparação |
| 7 | salva modelos e `best_params.json` (usado pelo grupo `sazonal_seco`) |
| 8 | teste temporal: modelo A aplicado ao VNIR + TIR do B × SWIR medido no A |
| 9 | treino direto nos pares (B → SWIR do A), validado com célula de fora; salva `best_params_pareado.json` (usado pelo grupo `sazonal_pareado`) |
| 10 | normalização B → A (`ols`/`rma`/`qm` × `global`/`celula` × pontos estáveis) |

Decisão para a etapa 6, comparando o R² do SWIR fora da amostra:
- seção 9 (`r2_pareado_celula_fora`) melhor → grupo `sazonal_pareado`;
- seção 10 com normalização melhor → grupo `sazonal_seco` com reta `global` ou
  `celula` (o GEE só aplica `ols`; se `rma`/`qm` ganharem com folga, o script precisa
  desse método);
- na dúvida, rode os dois e compare no mapa.

## Etapa 6 — SWIR estimado no GEE

```bash
cd src/modeling
python estimar_swir_rf_v2.py
```

Usa a conta ativa (`CONTA = None`); depois da etapa 4 ela é a `caatinga04`.

1. Grupos: `sazonal_seco`, `sazonal_pareado` ou os dois (`comparar`).
2. Normalização (só `sazonal_seco`): `nenhuma`, `celula` ou `global`
   (retas OLS dos pontos pareados → `models/sazonal_seco/normalizacao_ols_celula_coef.csv`).
3. Células, treino (`global` ou `por célula`, com vizinhas se faltar ponto) e modo
   (dry-run e depois exportar).

Hiperparâmetros: `models/sazonal_seco/best_params.json` (`sazonal_seco`) e
`best_params_pareado.json` (`sazonal_pareado`); sem o arquivo, usa `RF_PARAMS`.

Saída: `ASTER_MOSAIC_<grupo>_2009_2020_<LxCy>` em `mapbiomas-arida/mine/ASTER_MOSAIC_v2`
(VNIR/TIR medidos do B + SWIR estimado + banda `normalizacao`). Registro em
`registro_swir_rf_<data>.csv`. A célula sem mosaico B sai como `sem_entrada`.
