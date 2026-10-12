# Modelagem — Prospectra (bloco Irecê, BA)

Diário dos testes de modelagem: o que se testa, por quê, o que saiu e o que vem
depois. Cada seção registra a decisão tomada e a evidência que a sustenta, para
que o caminho possa ser refeito ou contestado.

Convenções dos dados (iguais em todo o pipeline v4 / zero-cloud):

- VNIR e SWIR em reflectância TOA × 1000 (INT32 nos assets e CSV); os notebooks
  dividem por 1000.
- TIR (B10–B14) em DN, sem conversão para temperatura.
- `saturado = 1` marca pixel com VNIR em DN ≥ 254. Nas cenas zero-cloud esses
  pixels se mostraram nuvem (TIR frio), não solo claro — ver
  `building_dataset/saturacao_landsat_piloto.py`. Saem de todas as amostras.

---

## 1. Problema central: estimar o SWIR depois de 2008

O detector SWIR do ASTER (B04–B09) parou em 2008-S2. As bandas SWIR são as que
mais informam sobre mineralogia (argilas, carbonatos — feição de 2,33 µm na B08),
então de 2009 a 2020 só existe VNIR + TIR. A proposta é aprender a relação
VNIR + TIR → SWIR no período em que as três existem (2000–2008) e aplicar a
2009–2020.

A pergunta que decide se isso funciona não é "o modelo acerta no período de
treino?", e sim **"o modelo transfere para outras datas?"**. Todos os testes
abaixo giram em torno disso.

## 2. Dados de treino: pontos pareados zero-cloud

**O que:** `building_dataset/pontos_zero_cloud_pareados.py` sorteia 2.000 pontos
por célula da grid 5 × 5 no mosaico zero-cloud 2000–2008 e amostra o mosaico
zero-cloud 2009–2020 nos **mesmos pontos**. `building_dataset/baixar_pontos_pareados.py`
baixa as tabelas para `src/Dados/ASTER_Samples_0cloud/pareados/`.

| Arquivo | Pontos | Conteúdo |
|---|---|---|
| `points_zeroCC_00_08.csv` | 34.000 (17 células) | VNIR + SWIR + TIR + data da cena |
| `points_zeroCC_09_20.csv` | 30.631 | VNIR + TIR + data da cena |
| `points_zeroCC_pareados.csv` | 30.631 | junção por `pid`, sufixos `_00_08` / `_09_20` |

**Por que:**

- Os CSV antigos (`samples_v2_ASTER_QualityMosaic_*`) vinham de mosaicos de
  qualidade com nuvem residual, TIR estourado em INT16 (100% = 32767) e SWIR
  saturado em 8–17% dos pontos. Treinar com isso ensina ruído.
- Cenas zero-cloud foram inspecionadas uma a uma; entram intactas (sem
  composição pixel a pixel entre cenas), então VNIR, SWIR e TIR de cada ponto
  são da mesma cena e da mesma hora.
- O pareamento permite um **teste temporal** honesto: no mesmo lugar, entra o
  VNIR + TIR de 2009–2020 e compara-se o SWIR estimado com o SWIR medido em
  2000–2008.

**Limites conhecidos:** 8 células (L1C3–L1C5, L2C5, L3C4, L3C5, L4C4, L4C5) não
têm cena zero-cloud em 2000–2008 e ficam sem pontos. L5C4 e L5C5 têm pouca
cobertura em 2009–2020 (573 e 345 pontos pareados).

## 3. Modelos de regressão (`train_swir_models.ipynb`)

Um modelo por banda SWIR, features VNIR + TIR padronizadas, split 80/20 aleatório
por ponto.

| Modelo | Por que entra | Situação |
|---|---|---|
| **PLS** | Linha de base linear; lida com a colinearidade forte entre VNIR e entre TIR | mantido |
| **XGBoost** | Captura relações não lineares; rápido | mantido |
| **RF** | É o `ee.Classifier.smileRandomForest` do GEE — os hiperparâmetros achados aqui vão direto para a aplicação no GEE (`best_params.json → rf.<banda>.gee`) | **entrou** no lugar do GPR |
| ~~GPR~~ | Custo O(n³), exigia subamostrar a 2.000 pontos e não superou o XGBoost | **removido** |

**Resultado:**

| Modelo | R² mesmo período (teste) | R² teste temporal (VNIR+TIR 2009–2020) |
|---|---|---|
| PLS | 0,90–0,93 | −0,5 a −0,8 |
| XGBoost | 0,85–0,91 | −0,22 a −0,28 |
| **RF** | **0,94–0,96** (RMSE 0,006–0,015) | −0,29 a −0,65 |

O RF da execução completa (30 iterações do `RandomizedSearchCV`) escolheu os
mesmos hiperparâmetros para as 6 bandas — no GEE: `numberOfTrees` 108,
`variablesPerSplit` 5, `minLeafPopulation` 2, `bagFraction` 0,983, sem
`maxNodes`. As colunas de PLS/XGBoost e o teste temporal vêm da busca reduzida
(indicativos) e precisam ser atualizados com a execução completa.

Dentro do período os modelos são bons, e o **RF é o escolhido** para a produção
(seção 6.1). Na transferência para 2009–2020 todos ficam **piores que prever a
média** (R² negativo). O problema não é o algoritmo.

## 4. Diagnóstico: por que não transfere (Célula 12A)

**O que:** assinaturas por banda (B01–B14) de 4 grupos de 10 pontos aleatórios,
2000–2008 ao lado de 2009–2020, mesmos pontos, escalas iguais.

**O que mostrou:**

- Correlação do mesmo ponto entre os dois períodos: B01/B02 ≈ 0,6; B3N 0,37;
  TIR 0,21–0,28.
- O TIR sobe em média ~100 DN de um período para outro e varia muito mais em
  2009–2020.
- As cenas de 2000–2008 amostradas são quase todas do 1º semestre / julho; as de
  2009–2020, de julho a outubro (estação seca, superfície mais quente). Pontos
  de cenas de outubro já em 2000–2008 (L4C3, 2002-10) têm TIR alto; a cena de
  junho de 2014 (L5C5) tem o TIR mais baixo de 2009–2020.

**Conclusão:** o TIR em DN carrega sobretudo a **data/estação da cena**
(temperatura), não a superfície. O modelo aprende a relação com cenas frias e
erra nas quentes. A B3N também muda com a fenologia da vegetação.

### 4.1 SWIR estimado nos mesmos pontos (Célula 21A)

**O que:** os mesmos 4 grupos da Célula 12A; à direita, o SWIR estimado pelo RF
a partir do VNIR + TIR de 2009–2020 (pontilhado), ao lado do SWIR medido em
2000–2008 à esquerda. `MODELO_EST` troca o modelo.

**O que mostrou:**

- A **forma** do espectro estimado é plausível: pico na B04, queda até a B09,
  feição da B08 presente.
- O **nível** sai deslocado para cima na maioria dos pontos (ex.: L1C2_127, B04
  medido ≈ 0,23 em 2000–2008 → estimado ≈ 0,32 em 2009–2020).
- O deslocamento acompanha o TIR: pontos de cenas de ago–out (TIR alto) recebem
  SWIR mais alto; o ponto da cena de junho de 2014 (L5C5_618, TIR frio) recebe o
  SWIR mais baixo.

**Conclusão:** confirma a seção 4 — o RF usa o TIR como "termômetro" e o SWIR
estimado herda a temperatura da data. Reforça o teste 1 da seção 5 (TIR
normalizado por cena) como prioridade.

### 4.2 Mesma estação e TIR anomalia (seção 9 do notebook)

**Teto do teste temporal:** usando a própria banda medida de 2009–2020 como
"previsão" da de 2000–2008 (sem modelo), o R² é −0,6 a −0,7 no VNIR e −0,8 a
−1,1 no TIR. O teste temporal com todos os pares mede sobretudo a mudança entre
datas, não o modelo.

**Pares da mesma estação (±1 mês, ~12,9 mil pares):** o R² do RF/XGBoost vai de
−0,26/−0,44 para **+0,56 a +0,62** — acima do teto das bandas medidas (VNIR
−0,09 a 0,22; TIR 0,37–0,60). O modelo funciona quando a estação bate.

**TIR como anomalia da cena:** melhora o teste com todos os pares (R² médio RF
−0,44 → −0,22; XGBoost −0,26 → −0,12), mas não melhora na mesma estação (≈ 0,55)
e perde um pouco no teste 00-08 (0,95 → 0,93). Remove parte do efeito de
temperatura, não todo — a sazonalidade da vegetação (B3N) continua.

(Busca reduzida — números indicativos.)

**Conclusão:** o limite é a **sazonalidade**, não o algoritmo. Para a aplicação
em 2009–2020, priorizar cenas da mesma época do ano das de treino, ou dar ao
modelo a informação de época (próximo teste).

### 4.3 Normalização radiométrica 2009–2020 → 2000–2008 (seção 10)

**O que:** antes do modelo de SWIR, cada banda VNIR/TIR de 2009–2020 é levada à
escala de 2000–2008 com os pontos pareados. Métodos `ols` / `rma` (eixo maior
reduzido) / `qm` (quantile mapping) × agrupamento `global` / `mes` / `cena` ×
todos os pontos ou só estáveis (|ΔNDVI| ≤ 0,05, 27% dos pares). Ajuste e
avaliação com uma célula de fora por vez (fora da amostra).

**Resultado (busca reduzida, indicativo), R² médio do SWIR estimado pelo RF:**

| Configuração | Todos os pares | Mesma estação | Razão de desvio padrão (VNIR / TIR) |
|---|---|---|---|
| sem normalização | −0,44 | 0,58 | 1,01 / 0,90 |
| **`ols` por cena, todos os pontos** | **+0,49** | 0,50 | 0,74 / 0,85 |
| `rma` por cena, todos | +0,32 | 0,39 | 1,07 / 1,15 |
| `rma` por cena, só estáveis | −0,08 | **0,65** | 0,94 / 1,20 |
| `ols` global, todos | +0,24 | 0,12 | 0,53 / 0,31 |

**Leitura:**

- Por **cena** é o agrupamento que mais ganha: o deslocamento é da data, não do
  período. Global quase não ajuda e achata a imagem (`ols` global: desvio padrão
  cai a 31% no TIR).
- `ols` por cena dá o maior ganho com todos os pares (−0,44 → +0,49) mantendo
  74–85% do contraste — é o candidato.
- Ajustar só com pontos estáveis melhora a mesma estação (0,65) mas não os pares
  em geral: a maioria dos pontos mudou de verdade entre as datas, e a
  transformação estável não "corrige" essa mudança — o que é o comportamento
  certo para geologia (estável no tempo), mas não aparece nesta métrica.
- Nas assinaturas (Célula 28), com a normalização o nível do SWIR estimado desce
  para perto do medido em 2000–2008, sem perder a ordem entre os pontos.

**Para aplicar no GEE:** ganho e offset por banda **por cena** — exige a banda
`data` no mosaico (zero-cloud tem; semestral v4 e `ASTER_MOSAIC` não) e cenas
2009–2020 que cruzem células com cobertura 2000–2008; cena sem pares suficientes
usa o ajuste global.

## 5. Próximos testes (em ordem)

0. **Época como feature** — acrescentar mês / dia do ano (seno e cosseno) da
   cena às features. Por quê: as cenas de treino vão de janeiro a outubro; o
   modelo pode aprender a corrigir o efeito da estação em vez de confundi-lo
   com material. Exige a banda `data` no mosaico de aplicação.
1. ~~**TIR normalizado por cena**~~ (feito — seção 4.2) — subtrair a média (ou z-score) do TIR da cena
   de cada ponto, usando a coluna `data`. Por quê: remove o nível de temperatura
   da data e mantém o contraste entre materiais (emissividade), que é o que
   deveria informar o SWIR.
2. **Sem TIR / com razões** — treinar só com VNIR, ou com razões entre bandas
   (B02/B01, B3N/B02, B14/B10…). Por quê: razões cancelam efeitos multiplicativos
   de iluminação e temperatura; mede quanto o TIR ajuda ou atrapalha.
3. **Pares da mesma estação** — restringir o teste temporal a pontos cujas cenas
   dos dois períodos caem na mesma época do ano. Por quê: separa o erro do
   modelo da mudança sazonal real da superfície.
4. **Validação espacial** — trocar o split aleatório por "uma célula de fora".
   Por quê: pontos vizinhos são quase iguais; o split aleatório infla o R².
5. **Execução completa do notebook** — rodar com `N_ITER_XGB = 60` e
   `N_ITER_RF = 30` e salvar (Célula 22) para gerar
   `models/zero_cloud/best_params.json`; o `estimar_swir_rf.py` passa a ler os
   hiperparâmetros por banda de lá.

Critério para levar um modelo ao GEE: R² do teste temporal positivo e estável
entre células.

## 6. Aplicação no GEE (`estimate_swir_gee.ipynb` / `.py`)

Treina `smileGradientTreeBoost` direto no GEE e aplica ao mosaico. Ainda lê os
pontos antigos (`projects/mapbiomas-arida/mine/points`). Só vale atualizar
depois que a seção 5 achar uma configuração que transfere; aí o RF com os
parâmetros do `best_params.json` é o candidato natural.

### 6.1 Produção por célula (`estimar_swir_rf.py`)

RF (`smileRandomForest` em REGRESSION, um por banda SWIR) treinado em
`points_zeroCC_00_08` e aplicado por célula aos mosaicos 2009–2020 de dois
grupos, salvos juntos na ImageCollection
`…/CAATINGA/ASTER_MOSAIC` (conta superconta). RF e não GTB porque foi o melhor
modelo no `train_swir_models.ipynb` (R² teste 0,94–0,96) e os hiperparâmetros
de lá (`numberOfTrees` 108, `variablesPerSplit` 5, `minLeafPopulation` 2,
`bagFraction` 0,983) passam direto para o GEE.

| Grupo (`grupo`, e no nome) | Entrada |
|---|---|
| `zero_cloud` | `ASTER_0cloud/ASTER_zero_cloud_2009_2020_<LxCy>` |
| `semestral_v4` | semestrais v4 da célula 2009–2020 (`ASTER_v4`), compostos por `qualityMosaic('quality')` |

Nome: `ASTER_MOSAIC_<grupo>_2009_2020_<LxCy>`. Bandas B01–B14, com B04–B09
estimadas (`swir = 'estimado_rf'`). Por que os dois grupos na mesma coleção:
comparar, célula a célula, o mosaico só de cenas limpas com o mosaico de maior
cobertura (com nuvem filtrada por quality), com o mesmo modelo de SWIR.
**Normalização (v2, `swir_rf_v2_norm`):** antes do RF, o VNIR/TIR de 2009–2020
passa pela reta OLS por cena da seção 4.3 (ganho e offset por banda, escolhidos
pela banda `data` com `remap`; cena com < 200 pares usa a reta global). Retas em
`models/zero_cloud/normalizacao_ols_cena_coef.csv` (23 de 34 cenas com reta
própria). `NORMALIZACAO = {'zero_cloud': 'cena', 'semestral_v4': 'nenhuma'}` —
os semestrais v4 não têm a banda `data`. A saída guarda VNIR/TIR **medidos** e o
SWIR estimado, mais a banda `normalizacao` (0 nenhuma, 1 global, 2 cena).
Teste em L3C3: mediana B04 251 → 230, B08 101 → 85; 36% dos pixels com reta da
cena, 64% com a global.

Ressalvas: a reta global tem ganho baixo (B14 0,31, B3N 0,40) e achata o que cai
nela; o grupo `semestral_v4` continua sem normalização, com o viés da seção 4.1. Rodar primeiro em dry-run (confere entradas e a leitura dos pontos pela
superconta) e só então exportar.

### 6.2 Mosaico final para inferência (`building_dataset/build_mosaico_final.py`)

**O que:** junta, por célula, os quatro mosaicos numa imagem só, na
ImageCollection `…/CAATINGA/aster-mosaico-final` (`aster-mosaico-final_<LxCy>`;
o mosaico único é `ee.ImageCollection(…).mosaic()`). Prioridade por pixel:

| `origem` | Camada | SWIR | Por que nessa posição |
|---|---|---|---|
| 1 | zero-cloud 2000–2008 | medido | sem nuvem e SWIR real — o dado mais confiável |
| 2 | zero-cloud 2009–2020 (`ASTER_MOSAIC`) | estimado (RF) | sem nuvem; cobre células sem zero-cloud antigo |
| 3 | semestrais v4 2000–2008 (`qualityMosaic`) | medido | preenche gaps; pode ter nuvem residual, mas o SWIR é real |
| 4 | semestrais v4 2009–2020 (`ASTER_MOSAIC`) | estimado (RF) | último recurso: nuvem residual + SWIR estimado |

Regras: em cada camada o pixel só vale com as 14 bandas e sem saturação (2008-S2,
sem SWIR, não entra na camada 3); a camada de maior prioridade com pixel válido
fica por cima. Bandas B01–B14 (INT32) + `origem` + `swir_estimado` (1 nas origens
2 e 4). Por que guardar a origem: as inferências de materiais (Maxent, índices de
SWIR) podem ser filtradas ou ponderadas pela confiabilidade do SWIR — os índices
de argila/carbonato só são diretos onde `swir_estimado = 0`.

Por que a nuvem pesa mais que o SWIR na ordem 2 × 3: nuvem contamina VNIR, SWIR
e TIR ao mesmo tempo e não é detectável depois; o SWIR estimado tem erro
conhecido e marcado (`swir_estimado`).

**Teste (leitura, superconta):** L3C3 → camadas 1 e 3 (~18% origem 1); L1C4 →
só camada 3. Camadas 2 e 4 ainda vazias.

**Pendências antes de exportar:**

- Rodar `estimar_swir_rf.py` (dois grupos) — sem `ASTER_MOSAIC`, ficam buracos
  onde só existe 2009–2020.
- `ASTER_v4` (superconta) tem 2000–2008 só de 2004 e 2006–2008; a coleção v4 das
  contas caatinga (`projects/mapbiomas-arida/mosaic_aster_v4`) não é legível pela
  superconta — compartilhar ou completar o v4 na superconta.

## 7. Favorabilidade para fosfato — Maxent (`maxent_fosfato.py`)

**O que:** mapa contínuo 0–1 de semelhança com as ocorrências de fosfato
(`projects/mapbiomas-arida/mine/pontos_fosfato`, 29 pontos na grid: 12 depósitos,
17 ocorrências) com `ee.Classifier.amnhMaxent` sobre um mosaico final
(`ASSET_MOSAICO`, a preencher).

**Por que Maxent:** só existem presenças; não há amostras confiáveis de "não
fosfato". Maxent contrasta presença × fundo (pixels aleatórios da paisagem) e
devolve um gradiente, não uma classe.

**Escolhas e porquês:**

- Presença = até 20 pixels num buffer de 500 m de cada ponto — o erro de posição
  dos pontos (carta 1:250.000) é de 250–1.000 m, muito maior que o pixel.
- Só funções linear e quadrática e `betaMultiplier = 2` — com 29 pontos, funções
  mais flexíveis (hinge, threshold) decoram as ocorrências.
- Validação deixando uma célula de fora por vez — os pontos são agrupados
  (Lapão, Morro Pelado, Lagedão…); deixar vizinhos no treino infla a AUC.
- Amostras salvas em asset antes de validar/mapear — a amostragem interativa
  estoura o limite do GEE ("Too many concurrent aggregations").

**Features (`FEATURES_MODO`, padrão `bandas+indices`):** além das bandas B01–B14,
os índices de `features_process/calculo_spectrals_index.js` — vegetação (NDVI,
SAVI), ferro (Fe³⁺, óxido férrico, gossan, ferroso-silicato), argilas
(sericita/illita, alunita/caulinita, muscovita, caulinita), carbonato (dolomita,
B08/B06, B13/B14), sílica (B14/B12), composições Sultan e textura GLCM da B3N
(contraste, variância, entropia, homogeneidade). Por quê: razões entre bandas
isolam as feições de absorção (Al-OH a 2,20 µm, CO₃ a 2,33 µm, Fe³⁺ no VNIR) e
cancelam iluminação e brilho geral; a textura separa afloramento de solo/cultivo.
Índices repetidos no JS (ex.: `abrams_r` = `sultan_r`, `abrams_b` = Fe³⁺) entram
uma vez — Maxent com features duplicadas só divide o peso. Mosaico sem SWIR →
ficam só os índices de VNIR/TIR e textura (17 features contra 35).

**Cobertura (achado):** num buffer de 500 m, **nenhum** dos 29 pontos tem pixel
no mosaico zero-cloud 2000–2008 (as células dos pontos são justamente as sem cena
zero-cloud nesse período), e só 14 têm no zero-cloud 2009–2020. O mosaico de
entrada do Maxent precisa cobrir os pontos: o candidato natural é o
`aster-mosaico-final` (seção 6.2), que cobre a grid toda; usar `swir_estimado`
para checar se as presenças caem em SWIR medido ou estimado.

**Cuidados de interpretação:** o mapa mede semelhança espectral com o entorno
das ocorrências — rocha hospedeira (carbonatos da Fm. Salitre) e cobertura —,
não o fósforo. Calcário sem fosfato também pontua alto. O SWIR (B08, carbonato)
só existe no mosaico 2000–2008; para 2009–2020 dependeria do SWIR estimado.

**Resultado:** pendente. O teste com o mosaico zero-cloud 2000–2008
(`maxent_fosfato_amostras_teste_zero_cloud_2000_2008`, só bandas) não serve: sem
cobertura nos pontos não há amostras de presença. Refazer com o
`aster-mosaico-final` e `FEATURES_MODO = 'bandas+indices'`.

---

## Arquivos

| Arquivo | Papel |
|---|---|
| `train_swir_models.ipynb` | treino e comparação PLS × XGBoost × RF; teste temporal; diagnósticos (Células 12A e 21A) |
| `train_swir_models.py` | versão script antiga (CSV `samples_v2`, com GPR) — não atualizada |
| `estimate_swir_gee.ipynb` / `.py` | GTB no GEE com os pontos antigos — não atualizado (seção 6) |
| `estimar_swir_rf.py` | produção: SWIR estimado por RF, por célula, → `ASTER_MOSAIC` |
| `maxent_fosfato.py` | favorabilidade para fosfato (Maxent), bandas + índices espectrais + GLCM |
| `../building_dataset/build_mosaico_final.py` | mosaico final por prioridade (seção 6.2) → `aster-mosaico-final` |
| `models/zero_cloud/` | figuras, modelos e métricas da versão com pontos pareados |
| `apontes.md` | notas sobre a correção da saturação do VNIR com Landsat |
