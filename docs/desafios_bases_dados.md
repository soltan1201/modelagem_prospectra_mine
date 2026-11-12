# Desafios das bases de dados — nova área de Irecê (RF/GTB com geofísica)

Levantamento feito em 2026-10-26, com consultas só de leitura no GEE, para
planejar os modelos de probabilidade de fosfato (Random Forest e Gradient Tree
Boost) na nova área de Irecê. Os modelos juntam a geofísica do acervo CPRM, as
bandas ASTER e os índices espectrais.

```js
var area_IRECE = ee.Geometry.Rectangle([-42.31762892604435, -12.359936623152336,
                                        -41.2382221877631,  -10.5375816299746]);
```

A área também está em `src/modeling/area_reduzida/js/novas_areas_estudo.js`. Ela toca as células
L2C3, L2C4, L3C3 e L3C4 da grid. O mosaico `aster-mosaico-final` (v1) cobre a
área inteira.

## 1. O que existe na área

| Base | Asset (`projects/mapbiomas-arida/mine/…`) | Na área | Resolução | Observação |
|---|---|---|---|---|
| Ocorrências de fosfato | `pontos_referencia` (`alvo = fosfato`) | **13**, todas `papel = treino` | ponto | o acervo tem 14; a fusão de pontos a ≤ 500 m no `preparar_pontos_referencia.py` junta dois |
| Ocorrências do acervo | `acervo_cprm/ocorrencias_fosfato` | 14 | ponto | fonte das 13 acima |
| Outros alvos | `pontos_referencia` | carbonato 24, mvt 36, manganês 1, ouro 1 | ponto | não entram como presença; podem servir de contraste |
| Gamaespectrometria | `acervo_cprm/1042_gama` | eU_ppm, eTh_ppm, K_pct, CT, U_Th, U_K, Th_K | ~493 m | cobre ~63–65% da área |
| Magnetometria | `acervo_cprm/1042_mag` | MAG_nT | ~223 m | cobertura parecida com a do gama |
| Mosaico ASTER v1 | `aster-mosaico-final` | 14 bandas, 18 índices, 4 texturas GLCM | 30 m | cobre a área toda |
| Litologia 1:250 mil | `acervo_cprm/lito_jacobina_250k` | 335 polígonos | vetor | também `lito_bahia` (escala menor) |
| Estruturas 1:250 mil | `acervo_cprm/estrut_jacobina_250k` | 1 027 feições | vetor | também `estrut_bahia` |
| Geoquímica de rocha | `acervo_cprm/geoquimica_p2o5` | 79 amostras | ponto | P₂O₅ de 0,01% a 40,27%, média 1,03% |
| Dolinas | `acervo_cprm/dolinas_morro_chapeu` | — | vetor | possível covariável de carste |
| Recursos minerais | `acervo_cprm/recursos_minerais` | — | ponto | consulta |

As 34 ocorrências de validação da R2 ficam fora desta área e fora da cobertura
geofísica. Por isso não servem para avaliar um modelo que use geofísica.

## 2. Primeiro sinal da geofísica nas ocorrências

Valores nas ocorrências cobertas, comparados com a distribuição da área:

| Variável | Nas ocorrências | Área: p10 | Área: p50 | Área: p90 | Leitura |
|---|---|---|---|---|---|
| eU (ppm) | 1,27 – 3,7 | — | 1,39 | 2,3 | a maioria acima da mediana, várias acima do p90 |
| U/Th | 0,11 – 0,29 | — | 0,14 | 0,25 | sobe um pouco |
| MAG (nT) | −13 a +4 | −141 | −15 | 93 | perto da mediana, sem sinal claro |

O eU alto é coerente com o fosforito, que concentra urânio. O campo magnético
bruto não separa as ocorrências; as derivadas (gradiente, *tilt*) ainda precisam
ser testadas.

## 3. Desafios

### 3.1 Poucas ocorrências: o n real é 11–13

São 13 ocorrências, e 2 ficam fora da geofísica. Um buffer de 250 m tem
cerca de 218 pixels de 30 m. Esses pixels não são observações independentes:
no gama, o buffer inteiro cabe num único pixel de ~500 m. Todos os pixels do
buffer recebem o mesmo valor de eU, eTh e K.

Consequências:

- a validação precisa deixar a ocorrência inteira de fora
  (*leave-one-occurrence-out*) ou um bloco espacial, nunca pixels soltos;
- cada ocorrência deve pesar igual no treino (peso 1/n_pixels do seu buffer),
  senão as ocorrências com buffer maior ou mais pixels válidos dominam;
- qualquer importância de variável calculada uma vez só, com 11–13 pontos, é
  instável.

### 3.2 Cobertura parcial da geofísica

O RF/GTB do GEE não aceita pixel sem dado. A área com gama e mag cobre ~64% do
retângulo. Há duas saídas:

1. mapear só a área com geofísica;
2. fazer dois modelos: ASTER puro na área toda e ASTER + geofísica na área coberta.

### 3.3 Resoluções diferentes (30 m × 220 m × 500 m)

**Decisão:** reamostrar a geofísica para 30 m, a grade do ASTER, repetindo o
mesmo valor em todos os pixels de 30 m que caem num pixel de 500 m (gama) ou de
220 m (mag). É a reamostragem pelo vizinho mais próximo, sem interpolação.

No GEE isso é o padrão: basta não chamar `.resample('bilinear')` e amostrar ou
exportar na escala de 30 m, na projeção do mosaico ASTER:

```js
var proj30 = ee.Image('projects/mapbiomas-arida/mine/aster-mosaico-final').projection();
var gama30 = ee.Image('projects/mapbiomas-arida/mine/acervo_cprm/1042_gama')
               .reproject({crs: proj30.crs(), scale: 30});   // vizinho mais próximo
```

Em Python (`ee`), o equivalente é o mesmo `reproject` sem `resample`. Para RF/GTB
nem o `reproject` é obrigatório: `sampleRegions` e o export em `scale: 30`
já leem o pixel de 500 m que contém cada pixel de 30 m.

O que essa escolha não resolve:

- a reamostragem não cria informação; os pixels do buffer continuam com o mesmo
  valor de gama (ver 3.1);
- o mapa de probabilidade vai mostrar blocos de ~500 m nas bordas dos pixels do
  gama, porque a árvore corta em degraus. É esperado; a bilinear suavizaria o
  desenho mas mudaria os valores medidos;
- a borda da cobertura geofísica fica em degraus de 500 m; o recorte da área de
  modelagem deve usar a máscara do gama reamostrado, não um polígono suavizado.

### 3.4 Variáveis redundantes

As razões do gama (U/Th, U/K, Th/K) derivam de eU, eTh e K; CT é a soma. As
bandas TIR (B10–B14) e vários índices também andam juntos. Sem um filtro de
redundância, a importância se divide entre variáveis parecidas e engana. No
SHAP da R2, por exemplo, `ferrous_sil_5_4` saiu por |r| = 0,92 com `sultan_r`.

### 3.5 Desenho do fundo

- sorteio na área com geofísica, a mais de 2 km das ocorrências;
- cerca de 3 × o número de pixels de presença, como no `exportar_favorabilidade_ml.py`
  (p ≈ 0,25 fica sendo o neutro);
- o mapa depende do sorteio; o ideal é repetir com sementes diferentes e tirar a
  média.
- **opção: fundo por cluster** (`FUNDO_MODO = 'cluster'` no notebook). Só o fundo é
  agrupado (k-means com K = 10 fixo, só o fundo com geofísica), em
  14 variáveis sem litologia: eU, eTh, U/Th, Th/K, eU_res, MAG, mag_res, mag_thg, Sultan
  (r, g, b), B10, quartz_14_12 e phosphorite_rbd. O "cluster fosfato" é o centro das
  ocorrências, e o raio do fosfato é a distância da ocorrência mais afastada desse centro.
  Entram no treino só os clusters sem ocorrência e com o centro além do raio (distância
  euclidiana nas variáveis padronizadas). A avaliação continua contra o fundo aleatório
  inteiro; contra o fundo distante, a AUC sairia inflada. O export refaz no GEE a atribuição
  de cada pixel de fundo ao centro mais próximo, com a padronização e os centros do
  notebook (`models/cluster_fundo_AR_<area>.json`), e acrescenta `_fcl` ao nome.
  Os conjuntos da área toda (sem gama/mag em ~1/3 do fundo) treinam com o fundo
  aleatório. Plano completo: `src/modeling/README.md`, seção 8.

### 3.6 Validação independente

Não há ocorrências de validação dentro da área com geofísica. As 79 amostras de
P₂O₅ em rocha são a única base que não entra no treino. Elas servem para conferir
o mapa: correlação de Spearman entre o p e o P₂O₅, e a fração das amostras com
P₂O₅ alto que caem no topo do mapa.

## 4. Análise recomendada para selecionar as variáveis

Em quatro etapas, nesta ordem:

| Etapa | O que faz | Por quê |
|---|---|---|
| 1. Filtro de redundância | Spearman entre as variáveis no fundo; agrupamento hierárquico; uma variável por grupo com \|r\| > 0,85 (ou VIF) | evita dividir a importância entre variáveis parecidas |
| 2. Ablação por grupo | validação *leave-one-occurrence-out* em quatro conjuntos: ASTER (bandas + índices); + gama; + gama + mag; + litologia e distância a estruturas | responde direto se a geofísica ajuda, e qual parte |
| 3. Ranking estável | SHAP fora da amostra (o mesmo do `avaliar_maxent.ipynb`) repetido em reamostragens das ocorrências; entra a variável que fica no topo em ≥ 70% das repetições | reduz o acaso de uma divisão só com 11–13 pontos |
| 4. Conferência independente | p do mapa contra o P₂O₅ das 79 amostras de rocha | é o único teste que não usa as ocorrências |

Métricas, como no restante do projeto: AP, índice de Boyce, AUC por ponto e
top10. Acurácia não serve com classes tão desbalanceadas.

O que evitar:

- importância Gini/MDI do RF dentro da amostra: favorece variáveis contínuas e
  lisas, como o gama de 500 m, que "decora" a posição dos poucos pontos;
- seleção por p-valor univariado;
- Boruta como critério único: com n = 13 tende a aceitar quase tudo ou quase
  nada. Serve só como segunda opinião.

## 5. Variáveis candidatas

| Grupo | Variáveis |
|---|---|
| ASTER | 14 bandas, 18 índices (`docs/indices_espectrais/`), 4 texturas GLCM, as 36 de `calcular_features` |
| Gama | eU, eTh, K, CT, U/Th, U/K, Th/K e `eU_res` (eU menos a média num raio de 5 km: anomalia local de urânio) |
| Mag | MAG, `mag_thg` (amplitude do gradiente horizontal, nT/km) e `mag_res` (MAG menos a média num raio de 5 km). A *tilt derivative* ficou de fora: precisa da derivada vertical, que pede FFT e não existe pronta no GEE |
| Geologia | `lito_<grupo>` = 1 dentro do grupo de unidades (prefixo do `SIGLA_UNID`: Nova América NP2sn, Jussara NP2sj, Gabriel NP2sg, Irecê NP2sr, Bebedouro NP2be, Morro do Chapéu M3N1mc, Caboclo MPc, Tombador MPt, coberturas Q) e distância, em km, a falhas, dobras e lineamentos (até 20 km) |

As derivadas do mag devem ser calculadas na resolução original (~220 m) e só
depois reamostradas para 30 m; derivar sobre a grade de 30 m em blocos gera
degraus artificiais.

## 6. Ocorrências, fundo e saída

- **Presença:** pixels de 30 m dentro do buffer de 250 m de cada uma das
  ocorrências com geofísica, com peso igual por ocorrência.
- **Fundo:** como em 3.5.
- **Modelos:** `smileRandomForest` e `smileGradientTreeBoost` (80 árvores,
  shrinkage 0,1), parâmetros espelhados no notebook.
- **Saída:** `probability` × 10000, numa coleção separada, por exemplo
  `projects/mapbiomas-arida/mine/favorabilidade_AR`.

## 7. Decisões

| Questão | Situação |
|---|---|
| Resolução da geofísica | **decidido:** reamostrar para 30 m pelo vizinho mais próximo (mesmo valor no pixel de 500 m) |
| Buffer das ocorrências | **decidido:** 250 m |
| Área do mapa: só com geofísica, ou também um modelo só com ASTER na área toda | a ablação testa as duas (`aster` e `aster_area_toda`); o export aceita qualquer conjunto |
| Incluir litologia e distância a estruturas | entram como grupo `geologia` na ablação. Cuidado: 10 das 13 ocorrências estão na Formação Nova América, então `lito_nova_america` tende a dominar e o mapa fica preso à unidade mapeada a 1:250 mil |

## 8. Implementação (prefixo `_AR`, áreas reduzidas)

| Passo | Arquivo (`src/modeling/`) | Saída |
|---|---|---|
| 1 | `amostras_AR.py` (task) e `amostras_AR.py --baixar` | asset `amostras_AR_irece`, `src/Dados/maxent/amostras_AR_irece.csv` e `_grupos.json` |
| 2 | `selecao_variaveis_AR.ipynb` | `models/variaveis_AR_irece.json`, CSVs `selecao_AR_irece_*`, figuras em `models/figuras_selecao_AR/` |
| 3 | `exportar_favorabilidade_AR.py` | coleção `favorabilidade_AR`, imagens `<RF\|GTB>_fosfato_AR_irece_<conjunto>` |

Detalhes que valem para os três:
- a mesma imagem de variáveis (`amostras_AR.imagem_variaveis()`) serve à amostra e ao mapa;
- cada ocorrência entra com 100 pixels sorteados (`N_POR_PONTO`), no notebook e no GEE;
- as 11 ocorrências com geofísica estão num trecho de ~8 × 10 km (vizinha mais próxima a
  0,6–3,2 km); agrupar por distância dava um fold só. Os folds agora são 5 grupos de
  ocorrências vizinhas (k-means das coordenadas), e o fundo vai em blocos de 10 km
  sorteados entre os folds. A ocorrência de teste sempre tem vizinhas no treino; a AUC
  mede "reconhece o distrito" e é otimista. O teste independente é a geoquímica;
- para ver o SHAP, o notebook precisa do pacote `shap` no ambiente.
