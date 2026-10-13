# Acervo CPRM — dados baixados para a área de Irecê

Inventário do que foi baixado do Acervo SIG do SGB/CPRM (plugin QGIS `acervo_cprm`) para o projeto de prospecção de **fosfato na Bacia de Irecê** (Grupo Una / Fm. Salitre) com mosaicos ASTER e MaxEnt. Para cada conjunto há o conteúdo, os números medidos nos próprios arquivos e uma avaliação de quanto ele ajuda no projeto. A coluna **Decisão** fica para preenchermos.

- **Área de estudo:** retângulo -44.46 / -14.79 / -39.20 / -9.42 (oeste/sul/leste/norte), cerca de 580 × 590 km.
- **Referência "Irecê":** (-41.86, -11.30). As distâncias abaixo são medidas a partir desse ponto.
- **Onde estão os dados:** `/mnt/Almacen/shpIRECE/` (1,4 GB, 86 pacotes, 469 shapefiles, 11 grids). O `manifesto.csv` lista o id, o link e os arquivos de cada pacote.
- **Data do levantamento:** 2026-10-12.

---

## 1. Achados que mudam o plano

| # | Achado | Consequência |
|---|---|---|
| 1 | **Os `.tif` de `942_1042-Geotif` são imagens coloridas (RGB 0–255), não valores físicos.** | **Resolvido:** os grids numéricos do pacote `1042-XYZ` foram convertidos para GeoTIFF float32 em `aerogeofisica/1042_grids/` (eU, eTh, K, CT, MAG, U/Th, U/K, Th/K). Ver seção 6. |
| 2 | **Há 41 ocorrências de fosfato na área (pontos únicos) e 13 delas ficam a até 63 km de Irecê**: a mina Irecê/Fazenda Juazeiro, mais Lapão, Faz. São Luiz, Faz. Três Irmãs, Nova Conquista, Lajedão I/II, Achado, Jacaré, Morro Pelado, Morrinhos, Lagedinho e Junco. | São os **pontos de presença** do MaxEnt. Vêm do `recmin` das cartas ao milionésimo (SC.24/SD.24) e se repetem na Geologia da Bahia. |
| 3 | **Nenhuma carta 1:100.000 cobre a Bacia de Irecê.** A folha 1:100k mais próxima fica a cerca de 93 km. | Na bacia, a melhor geologia é a Jacobina 1:250.000 (contém Irecê) e as cartas ao milionésimo. |
| 4 | **Há análises de P₂O₅ na área**: Agrominerais Irecê–Jaguarari (46 análises, máximo 16,5%, 4 com ≥5%) e Calcários da Bahia (288 análises, máximo 40,3%). | Servem de presença/ausência complementar e de validação. Os valores estão como texto com vírgula decimal (`'7,3'`) e precisam de conversão. |
| 5 | **Sistemas de coordenadas misturados**: 362 shapefiles em SIRGAS 2000 (EPSG:4674), 106 em WGS 84 (EPSG:4326) e 1 em UTM 24S (EPSG:31984). | Reprojetar tudo para um CRS único antes de juntar as camadas. |

---

## 2. Dicionário rápido dos temas guardados

| Tema | Geometria | Conteúdo | Uso no projeto |
|---|---|---|---|
| `lito` | polígono | unidade geológica: `SIGLA_UNID`, `NOME_UNIDA`, `HIERARQUIA`, idade (`IDADE_MAX/MIN`, era, período), `AMBSEDIMEN`, `SISTSEDIME`… | Mapear Fm. Salitre / Grupo Una (hospedeiros); máscara ou covariável categórica no MaxEnt; verdade de campo para validar os índices ASTER. |
| `recmin` / `recmin_mapa` | ponto | ocorrência mineral: `SUBS_PRINC` (ou `SUBST_PRIN`; fósforo = **`P`** ou "Fósforo"), `STATUS_ECO` (mina, garimpo, não explotado), `TOPONIMIA`, `ROCHA_HOSP`, `ROCHA_ENCA`, `MODELO_DEP` | Pontos de presença de fosfato; outros bens para checar o viés de amostragem. |
| `estrutura` / `estrut` | linha | falhas, fraturas, dobras, lineamentos, contatos | Covariável "distância a estruturas". |
| `afloramento` | ponto | estações de campo com descrição da rocha | Verdade de campo para a classificação ASTER. |
| `limite_da_folha` | polígono | contorno da folha | Recortar pela área e controlar a cobertura. |

---

## 3. Mapas geológicos (vetores)

Números = feições que cruzam a área de estudo.

### 3.1 `milionesimo/` — Cartas Geológicas ao Milionésimo (GIS Brasil 2004)
Folhas SC.23 (Rio São Francisco), SC.24 (Aracaju), SD.23 (Brasília) e SD.24 (Salvador). Juntas cobrem a área inteira.

| Tema | Feições na área |
|---|---|
| lito | 2.792 polígonos |
| recmin | 2.694 pontos |
| estrutura | 1.903 linhas |
| afloramento | 10.491 pontos |

- **Unidades-alvo presentes:** Salitre, Bebedouro, Irecê, Caboclo, Morro do Chapéu, Tombador, Sete Lagoas.
- **Fosfato:** é a principal fonte. A SC.24 tem a mina Irecê/Faz. Juazeiro e as 12 ocorrências vizinhas, todas a até 63 km.
- **Limitação:** escala 1:1M; os contatos são grosseiros para um pixel ASTER de 15–90 m.
- **Avaliação: AJUDA MUITO** (base contínua e pontos de presença).
- **Decisão:** ☐ usar ☐ descartar

### 3.2 `250k/` — Jacobina SC.24-Y-C (1:250.000)

| Tema | Feições na área |
|---|---|
| lito | 653 |
| estrutura | 1.979 |
| afloramento | 39 |

- **Contém Irecê.** Unidades presentes: Caboclo, Bebedouro, **Irecê** (25 polígonos), Morro do Chapéu, Tombador.
- **Limitação:** o pacote não traz `recmin`.
- **Avaliação: AJUDA MUITO.** É a geologia mais detalhada disponível sobre a bacia.
- **Decisão:** ☐ usar ☐ descartar

### 3.3 `100k/` — 56 Cartas Geológicas 1:100.000

| Tema | Feições na área |
|---|---|
| lito | 11.223 |
| recmin | 871 |
| estrutura | 58.323 |
| afloramento | 7.184 |

- Quando a folha tinha duas versões, ficou só a mais nova.
- **Cobertura:** as bordas da bacia — Chapada Diamantina, Jacobina, Campo Formoso, Sobradinho, Paramirim, Contas. **Nenhuma cobre a bacia em si.**
- **Unidades-alvo presentes:** Salitre (69 polígonos), Lagoa do Jacaré (51), Bebedouro, Caboclo, Morro do Chapéu, Tombador, Sete Lagoas, Serra de Santa Helena. As ocorrências de Salitre fora de Irecê são úteis para testar se o modelo generaliza.
- **Fosfato:** só 3 ocorrências (Serrote do Sobrado, Serrote Sacrabeto e Serrote da Batateira), todas no norte (Sobradinho), a cerca de 225–245 km.
- Fazem parte do conjunto, entre outras, as folhas Jacobina, Piritiba, Mirangaba, Lajedinho, Brejão da Caatinga, Tombador, Serra dos Brejinhos, Boquira, Caldeirão Grande, Mundo Novo, Campo Formoso, Itaetê, Ruy Barbosa e Senhor do Bonfim. A lista completa e as distâncias estão em `camadas_acervo_cprm_area_estudo.txt`.
- **Avaliação: AJUDA** para a validação regional e para as estruturas. Ainda falta filtrar pela área e juntar num arquivo só.
- **Decisão:** ☐ todas ☐ só as até ___ km de Irecê ☐ descartar

### 3.4 `50k/` — 11 Cartas Geológicas 1:50.000

| Tema | Feições na área |
|---|---|
| lito | 1.067 |
| recmin | 87 |
| estrutura | 4.242 |

- **Localização:** todas na região de Sobradinho/Remanso (SC.24-V-C, SC.23-X-D), a 140–230 km ao norte. Muita Salitre mapeada (72 polígonos).
- **Avaliação: TALVEZ.** Detalhe alto, mas longe do alvo. Pode servir como área de teste do modelo em Salitre fora de Irecê.
- **Decisão:** ☐ usar ☐ descartar

### 3.5 `regional/` — mapas sem código de folha

| id | Conjunto | O que tem de útil | Avaliação |
|---|---|---|---|
| 5984 | Geologia do Estado da Bahia (1:1M) | lito, recmin, estrut, afloramento contínuos para o estado. O recmin **repete os mesmos pontos** do milionésimo, sem toponímia. | AJUDA (substitui as 4 folhas do milionésimo num arquivo só) |
| 6691 | Projeto Morro do Chapéu (1995) | litologia de detalhe, com fácies e ambiente deposicional (inter-supramaré, submaré, planície de maré); **89 dolinas** (carste); 49 ocorrências (diamante, barita, diatomita — **nenhuma de fosfato**); 54 amostras com anomalias geoquímicas. **Atenção:** `MorChap_InfraEstrutura_ponto` entrou por engano (escolas e serviços) e pode ser apagado. | AJUDA localmente (vizinho de Irecê; o carste é relevante para o fosfato supergênico) |
| 6021 | Mapa Tectônico-Geocronológico da Bahia | domínios tectônicos, estruturas, recmin (1 ponto de fósforo a 13 km de Irecê) | TALVEZ (contexto) |
| 6714 | Geologia da Bacia do São Francisco | lito e estrutura regionais | TALVEZ (contexto) |
| 3594 | ARIM Serra de Jacobina | lito, recmin, estrutura (borda leste da bacia) | TALVEZ |
| 3395 | Ibitiara–Rio de Contas | lito, recmin, estrutura (Chapada, ao sul) | TALVEZ |
| 3398 | Barra–Oliveira dos Brejinhos | lito, recmin, estrutura (oeste) | TALVEZ |
| 3528 | ARIM Remanso–Sobradinho | lito, recmin com 2 ocorrências de fósforo (norte) | TALVEZ |

- **Decisão:** ☐ 5984 ☐ 6691 ☐ 6021 ☐ 6714 ☐ 3594 ☐ 3395 ☐ 3398 ☐ 3528

**Unidades-alvo nos regionais:** Tombador (232 polígonos), Caboclo (109), Morro do Chapéu (93), Salitre (70), Lagoa do Jacaré (65), Bebedouro (50), Sete Lagoas (48), Serra de Santa Helena (21), Irecê (14).

---

## 4. Ocorrências de fosfato na área (pontos de presença)

Há 88 registros, que se reduzem a **41 pontos únicos** porque as mesmas ocorrências aparecem no milionésimo e nos mapas regionais.

| Grupo | Pontos | Distância de Irecê | Status | Fonte |
|---|---|---|---|---|
| Bacia de Irecê | 13: Irecê/Faz. Juazeiro (**mina**), Faz. São Luiz, Faz. Três Irmãs, Lapão, Nova Conquista, Lajedão I e II, Achado, Jacaré, Morro Pelado, Morrinhos, Lagedinho, Junco | 5–63 km | 1 mina, 12 não explotados | `milionesimo/741_Aracaju_SC.24/sc24_recmin.shp` |
| Leste (~-39.7 a -40.0, -11.2 a -12.1) | ~22: fazendas Gameleira, Rio do Peixe, Pilão, Curralinho, Panelas, Retiro, Tomé… | 200–245 km | maioria garimpo / não explotado | SC.24 e SD.24 |
| Norte (Sobradinho) | 3: Serrote do Sobrado, Sacrabeto, Batateira | 224–247 km | não explotado | 100k, 3528, SC.24 |

**Como filtrar:** no milionésimo, selecionar `SUBST_PRIN = 'P'`. Nas folhas 100k e no ARIM, selecionar `SUBS_PRINC`/`SUBSTANCIA = 'Fósforo'`. O campo `ASSOC_MINE` também menciona apatita em outros 4 pontos, que ficaram de fora da contagem.

**Avaliação: ESSENCIAL.** Com só 13 pontos na bacia, é preciso cuidar da autocorrelação espacial e do viés de amostragem no MaxEnt.

---

## 5. Geoquímica (`geoquimica/`)

| id | Conjunto | Amostras | P / P₂O₅ | Avaliação |
|---|---|---|---|---|
| 7420 | **Agrominerais Irecê–Jaguarari (2016)** | 305 de rocha (300 na área, 80 a até 100 km) | **P₂O₅_pct**: 46 análises, máx. 16,5%, 4 com ≥5% (3 a até 100 km) | **AJUDA MUITO**: P₂O₅ medido dentro da bacia |
| 6719 | Calcários da Bahia (2021) | 2.245 de rocha (1.896 na área, 450 a até 100 km) | P₂O₅_pct: 288 análises, máx. 40,3%; quase todas baixas (só 2 com ≥1%) | AJUDA: carbonatos "negativos" (fundo) e validação |
| 4599 | Folha Morro do Chapéu PLGB (1997) | sedimento de corrente ~2.700, solo ~1.050, rocha ~40, concentrado de bateia ~260; shp dos pontos (1.098) | `P_ppm` (sedimento e solo), `P2O5_T_pct` (rocha), `FOSFATO` (bateia) | AJUDA: vizinha (44 km) |
| 4643 / 4646 | Geoquímica do Bambuí Etapas I e II – BA | sedimento ~21.000, solo ~4.100, rocha ~1.650, concentrado ~3.200; shp dos pontos | `P_ppm` | TALVEZ: Bambuí é correlato do Una; cobertura regional grande. Falta conferir quantos pontos caem na área. |

**Atenção nas planilhas:**
- Os valores vêm como texto com vírgula decimal (`'0,14'`) e às vezes com `<` (abaixo do limite de detecção).
- As coordenadas estão nas colunas `LATITUDE`/`LONGITUDE`.

- **Decisão:** ☐ 7420 ☐ 6719 ☐ 4599 ☐ 4643 ☐ 4646

---

## 6. Aerogeofísica — projeto 1042 Norte da Chapada Diamantina

- **Projeto:** 1042 Norte da Chapada Diamantina.
- **Extensão do grid:** -45.55 / -12.25 / -40.93 / -9.94, que contém Irecê.
- **Cobertura:** cerca de 27% da área de estudo. Não há dado a leste de -40.9° nem ao sul de -12.25°.
- **Resolução:** pixel de 0,0018° (~200 m), EPSG:4326.

| Arquivo | Produto |
|---|---|
| `1042_EU` | urânio equivalente |
| `1042_UTHRAZAO` | razão U/Th |
| `1042_UKRAZAO` | razão U/K |
| `1042_KPERC` | potássio (%) |
| `1042_ETH` | tório equivalente |
| `1042_THKRAZAO` | razão Th/K |
| `1042_CT` | contagem total |
| `1042_TERNARIO_RGB` | composição ternária K-Th-U |
| `1042_MAG` | campo magnético anômalo |
| `1042_1DV` | 1ª derivada vertical |
| `1042_SINAL` | sinal analítico |

**Os `.tif` de `942_1042-Geotif` são só imagens:** os 11 arquivos são imagens RGB de 3 bandas Byte (0–255), mapas renderizados em cor. Servem para inspeção visual, não como covariável.

**Grids numéricos (feito em 2026-10-12):** o pacote `1042-XYZ` (id 943) não traz dado de linha de voo. Ele traz **grids Geosoft `.grd`** com valores físicos, que estão em `aerogeofisica/943_1042-XYZ/`. O GDAL não lê esse formato; o script `geosoft_grd_para_tif.py` (nesta pasta) decodifica e grava GeoTIFF float32 em **`aerogeofisica/1042_grids/`**:

| Arquivo | Unidade | Pixel | p50 | p90 | máx. |
|---|---|---|---|---|---|
| `1042_eU_ppm.tif` | ppm | 0,00443° (~490 m) | 1,03 | 2,18 | 7,98 |
| `1042_eTh_ppm.tif` | ppm | idem | 5,77 | — | 39,4 |
| `1042_K_pct.tif` | % | idem | 0,39 | — | 3,79 |
| `1042_CT.tif` | contagem total | idem | 3,51 | — | 13,3 |
| `1042_U_Th.tif` | razão | idem | 0,18 | 0,43 | 3,66 |
| `1042_U_K.tif` | razão | idem | 2,43 | 5,19 | 40,1 |
| `1042_Th_K.tif` | razão | idem | 12,3 | — | 283 |
| `1042_MAG_nT.tif` | nT (campo anômalo) | 0,002° (~220 m) | -4,7 | — | 918 |

- **Grid:** EPSG:4326, NoData = -9999, extensão -45.55 / -12.25 / -40.93 / -9.94.
- **Conferência:** a máscara de dados coincide 100% com a das imagens RGB, então a orientação e o georreferenciamento estão certos.
- **Razões:** calculadas só onde K ≥ 0,1% (U/K, Th/K) ou eTh ≥ 0,5 ppm (U/Th); abaixo disso o denominador é ruído. Isso cortou cerca de 10–12 mil pixels.
- **Valores negativos pequenos:** eU, eTh e K têm valores um pouco abaixo de zero (até -0,7), comuns depois da correção gama. Por isso aparecem razões negativas. Convém truncar em 0 ou mascarar antes do MaxEnt.
- **Na mina de Irecê/Faz. Juazeiro (média 5×5 px):** eU 2,39 ppm (acima do p90 da área), U/Th 0,20 (perto da mediana), U/K 3,1 (entre p50 e p90). Lapão e Morro Pelado são parecidos. O eU destaca a bacia; as razões quase não destacam nessa resolução (~500 m).
- **Origem dos dados:** pelo `.doc` do pacote, os grids gama vêm do projeto BARMP e o MAG foi processado do XYZ da CPRM.
- **Cobertura:** só o oeste e o centro da área (~27%).
- **Avaliação: AJUDA.** O eU é a covariável geofísica mais promissora. Para o resto da área seriam precisos outros levantamentos.
- **Decisão:** ☐ usar eU ☐ usar U/Th ☐ usar U/K ☐ usar K/eTh/MAG ☐ procurar levantamentos para o leste/sul

---

## 7. Referências × covariáveis para os modelos

São dois papéis diferentes:

- **Referência** é o que o modelo tenta reproduzir ou o que usamos para avaliá-lo: pontos de presença, fundo/ausência e dados de validação. Entra como **rótulo**, nunca como preditor.
- **Covariável** é a variável explicativa, um raster contínuo ou categórico na malha do modelo, que o MaxEnt (ou outro modelo) usa para prever onde há fosfato.

**Regra:** a mesma informação não pode ser referência e covariável ao mesmo tempo. Se for, o modelo "aprende a resposta" e a validação fica otimista (vazamento). Exemplos de vazamento: usar "distância a ocorrências de fosfato" como covariável, ou interpolar o P₂O₅ das mesmas amostras usadas como presença.

Caminhos relativos a `/mnt/Almacen/shpIRECE/`.

### 7.1 Referências (rótulos e validação)

| # | Referência | Arquivo / filtro | Papel | Quantidade | Preparo |
|---|---|---|---|---|---|
| R1 | **Ocorrências de fosfato da bacia** | `milionesimo/741_Aracaju_SC.24/sc24_recmin.shp`, `SUBST_PRIN = 'P'` | **presença** (treino) | 13 pontos a até 63 km de Irecê (1 mina, 12 não explotados) | Tirar duplicatas da Geologia BA (5984); thinning espacial (1 ponto por célula da malha); com tão poucos pontos, avaliar por leave-one-out ou blocos espaciais. |
| R2 | Ocorrências de fosfato fora da bacia | mesmo filtro em SC.24/SD.24; `SUBS_PRINC = 'Fósforo'` em `100k/3516`, `100k/3551`, `regional/3528` | **validação independente** (generalização) | ~28 pontos (leste ~200–245 km; norte/Sobradinho ~225–247 km) | Não usar no treino; servem para testar se o modelo acerta fora de Irecê. Contexto geológico diferente: interpretar com cuidado. |
| R3 | **Amostras com P₂O₅ alto** | `geoquimica/7420…/projeto_agromin…xlsx`, coluna `P2O5_pct` ≥ 5 | presença complementar ou validação | 4 amostras na área (3 a até 100 km); máx. 16,5% | Converter vírgula decimal e `<`. Decidir se entram no treino (R1) ou só na validação. |
| R4 | Amostras com P₂O₅ baixo | `geoquimica/6719…/projeto_calcarios…xlsx` e 7420, `P2O5_pct` < 1 | **ausência / fundo informado** | ~330 análises na área (~95 a até 100 km) | Ausência verdadeira só se a amostra for do carbonato hospedeiro; útil para escolher o fundo do MaxEnt (fundo restrito aos carbonatos). |
| R5 | Outras ocorrências minerais (`recmin` sem P) | todas as folhas, `SUBST_PRIN ≠ 'P'` | **fundo com o mesmo viés de amostragem** (target-group background) | ~3,7–7,4 mil na área (com duplicatas entre séries) | Usar como pontos de fundo: corrigem o viés de "onde a CPRM foi a campo". |
| R6 | Afloramentos descritos | `*/afloramento.shp` (100k, 250k, regionais) | **verdade de campo** dos índices ASTER (litologia) | ~19 mil na área (sem contar o milionésimo) | Validar a classificação ASTER (carbonato, argilomineral, sílica), não o modelo de fosfato. |
| R7 | Unidades mapeadas (Salitre, Irecê, Bebedouro…) | `*/lito.shp`, `NOME_UNIDA` | **validação dos mapas ASTER** | Salitre: 69 (100k) + 72 (50k) + 70 (regional) polígonos; Irecê: 25 (250k) | Verdade de referência para o mapeamento litológico por sensoriamento remoto. Atenção: se a mesma lito entrar como covariável C3, ela não valida o modelo de fosfato. |

### 7.2 Covariáveis (preditores)

| # | Covariável | Fonte | Tipo | Cobertura | Preparo |
|---|---|---|---|---|---|
| C1 | **eU (ppm)** | `aerogeofisica/1042_grids/1042_eU_ppm.tif` | contínua | ~27% (oeste/centro, inclui Irecê) | Truncar negativos em 0; **mantém a resolução nativa (~490 m)**. Na mina: 2,39 ppm (> p90). **Principal covariável geofísica.** |
| C2 | Razões U/Th, U/K | `1042_grids/1042_U_Th.tif`, `1042_U_K.tif` | contínua | ~27%, com buracos onde K < 0,1% / eTh < 0,5 | Resolução nativa (~490 m). Destacam pouco na bacia (~500 m); testar com seleção de variáveis. |
| C3 | Unidade litológica / alvo | `250k/2617…/jacobina_lito.shp` (bacia) + `regional/5984…/bahia_lito.shp` (resto) | categórica (unidade-alvo sim/não, ou grupo litoestratigráfico) | 100% | Reclassificar `NOME_UNIDA` em classes (Salitre/Una carbonático, Bebedouro, Chapada/Tombador, embasamento…). Escala 1:250k–1:1M. |
| C4 | Distância ao contato da unidade-alvo | derivada de C3 | contínua | 100% | Distância euclidiana ao limite da Fm. Salitre / Grupo Una. |
| C5 | Distância a estruturas / densidade de lineamentos | `*/estrutura.shp` (250k na bacia; 5984 ou milionésimo no resto; 100k nas bordas) | contínua | 100% (detalhe desigual) | Juntar, reprojetar e calcular a distância e a densidade (km/km²). Cuidado: a densidade reflete o detalhe do mapeamento. |
| C6 | K %, eTh, contagem total | `1042_grids/1042_K_pct.tif`, `1042_eTh_ppm.tif`, `1042_CT.tif` | contínua | ~27% | Resolução nativa (~490 m). Contexto litológico (argilas, carbonato limpo × impuro). |
| C7 | Magnetometria | `1042_grids/1042_MAG_nT.tif` (~220 m) | contínua | ~27% | Resolução nativa (~220 m). Derivar a 1ª derivada vertical / sinal analítico para estruturas profundas. |
| C8 | Densidade de dolinas (carste) | `regional/6691…/Dolinas_area.shp` | contínua | só Morro do Chapéu | Cobertura local demais para ser covariável na área toda; manter como apoio. |
| C9 | P em sedimento de corrente / solo | `geoquimica/4599`, `4643`, `4646` (`P_ppm`) | contínua (interpolada) | regional, desigual | **Risco de vazamento** se R3/R4 vierem dos mesmos levantamentos; usar só se forem independentes. |
| — | Índices e bandas ASTER | mosaico final (`projects/mapbiomas-arida/mine/aster-mosaico-final`) | contínua | 100% | Fora do acervo CPRM, mas é a covariável central do projeto (carbonato, argila, sílica, TIR). |

### 7.3 Só contexto (nem referência nem covariável)

| Dado | Por quê |
|---|---|
| `geocronologia`, perfis geológicos, `medida_estrutural` | interpretação, sem valor preditivo direto |
| Domínios tectônicos (6021), Bacia do São Francisco (6714) | escala regional demais para a malha do modelo |
| `limite_da_folha` | só para recortar e checar cobertura |
| Imagens RGB `942_1042-Geotif` | renderização sem valor físico (substituídas por C1, C2, C6 e C7) |

### 7.4 Decisões

**Tomadas (2026-10-12):**
- **Resolução:** os dados que já vêm em ~500 m (grids gama do 1042) **ficam nessa resolução**. Não são reamostrados para a malha do ASTER. A magnetometria fica na resolução nativa dela (~220 m). Os vetores não têm resolução: sobem como vetores.
- **Extensão:** a extensão **de cada dado baixado** é a que vai para os assets do GEE. Nada é recortado pela área de estudo antes da exportação; o recorte, se preciso, fica para o GEE.
- **Exportação:** os dados vão para assets do GEE por um script próprio (seção 10). O script ainda não foi feito.

**Em aberto:**
- **Malha do modelo:** como combinar as covariáveis de ~500 m com as do ASTER (30/90 m) no momento da modelagem.
- **Extensão do modelo:** a área toda (onde a geofísica falta, o modelo usa só ASTER e geologia) ou só a área coberta pelo 1042?
- **R3:** as amostras com P₂O₅ ≥ 5% entram no treino ou só na validação?
- **Fundo:** aleatório, restrito aos carbonatos ou target-group (R5)?

## 8. O que dá para descartar

- `regional/6691…/MorChap_InfraEstrutura_ponto.*`: infraestrutura que entrou por engano no filtro de "estrut".
- Os `afloramento` do milionésimo (10 mil pontos de 1:1M), se não forem usados como verdade de campo.
- `aerogeofisica/942_1042-Geotif/` (imagens RGB), agora substituído por `1042_grids/`.

## 9. Arquivos de apoio nesta pasta

| Arquivo | O que é |
|---|---|
| `listar_QGIS.py` | classifica o catálogo do plugin contra a área (dentro / fora / sem caixa) |
| `inspecionar_pacotes.py` | lista o conteúdo dos zips remotos sem baixar |
| `baixar_acervo_irece.py` | baixa a seleção para `/mnt/Almacen/shpIRECE/` |
| `geosoft_grd_para_tif.py` | converte os grids Geosoft do 1042 para GeoTIFF e calcula as razões |
| `exportar_acervo_gee.py` | prepara os arquivos e os envia aos assets do GEE (seção 10) |
| `GUIA_GCLOUD_BUCKET_GEE.md` | guia de gcloud, bucket e envio com a conta `solkan.cengine17@gmail.com` |
| `camadas_acervo_cprm_area_estudo.txt` | catálogo completo da área, com distâncias e prioridades |
| `conteudo_pacotes_area_estudo.csv` | arquivo por arquivo de cada pacote |

---

## 10. Exportação para o GEE — `exportar_acervo_gee.py`

Regras: **resolução nativa** e **extensão do dado baixado**, sem recorte nem reamostragem. Destino: `projects/mapbiomas-arida/mine/acervo_cprm/`. Contas: as mesmas que escrevem em `mapbiomas-arida/mine` (todas, exceto a superconta).

### 10.1 Rasters → `ee.Image`

| Asset | Arquivo local (`/mnt/Almacen/shpIRECE/aerogeofisica/1042_grids/`) | Pixel | Tamanho | Extensão (O / S / L / N) | CRS |
|---|---|---|---|---|---|
| eU (ppm) | `1042_eU_ppm.tif` | 0,00443° (~490 m) | 1044 × 521 | -45.551 / -12.247 / -40.927 / -9.939 | EPSG:4326 |
| eTh (ppm) | `1042_eTh_ppm.tif` | idem | idem | idem | idem |
| K (%) | `1042_K_pct.tif` | idem | idem | idem | idem |
| Contagem total | `1042_CT.tif` | idem | idem | idem | idem |
| U/Th | `1042_U_Th.tif` | idem | idem | idem | idem |
| U/K | `1042_U_K.tif` | idem | idem | idem | idem |
| Th/K | `1042_Th_K.tif` | idem | idem | idem | idem |
| MAG (nT) | `1042_MAG_nT.tif` | 0,002° (~220 m) | 2168 × 1116 | -45.333 / -12.217 / -40.997 / -9.985 | EPSG:4326 |

- **Gama:** os 7 grids gama compartilham a mesma malha e podem subir como **uma imagem multibanda** (eU, eTh, K, CT, U_Th, U_K, Th_K). A MAG tem malha própria e sobe separada.
- **NoData:** -9999, que deve virar máscara no asset.
- **Negativos:** os valores negativos de eU/eTh/K **sobem como estão**. O corte em 0 é feito na modelagem, para o asset preservar o dado original.
- **Propriedades a gravar:** projeto 1042, fonte (SGB/CPRM, pacote id 943, grids BARMP), unidades e os limiares das razões (K ≥ 0,1%, eTh ≥ 0,5 ppm).

### 10.2 Vetores → `ee.FeatureCollection`

| Asset | Fonte | Feições (na área) | Preparo antes de subir |
|---|---|---|---|
| Ocorrências de fosfato (R1 + R2) | `sc24_recmin` / `sd24_recmin` (`SUBST_PRIN='P'`) + `Fósforo` das 100k/3528 | 41 pontos únicos | tirar duplicatas; campo `grupo` = bacia / leste / norte |
| Recursos minerais (R5) | todos os `recmin` | ~3,7–7,4 mil (com duplicatas) | juntar; tirar duplicatas por coordenada; manter `SUBST_PRIN`, `STATUS_ECO`, `TOPONIMIA`, `serie` |
| Geoquímica P₂O₅ (R3 + R4) | 7420 + 6719 (`.xlsx`) | 305 + 2.245 amostras | planilha → pontos (`LATITUDE`/`LONGITUDE`); `P2O5_pct` numérico (vírgula → ponto; `<` → limite) |
| Litologia (C3 / R7) | `250k/…/jacobina_lito` + `regional/5984…/bahia_lito` (+ 100k/50k se aprovados) | 653 + estado | reprojetar para EPSG:4326; manter `SIGLA_UNID`, `NOME_UNIDA`, `HIERARQUIA`, `serie` |
| Estruturas (C5) | `estrutura` de 250k + 5984 (+ 100k) | ~17 mil+ | idem; manter o tipo da estrutura |
| Dolinas (C8) | `regional/6691…/Dolinas_area` | 89 | idem |

- **CRS:** todos os vetores sobem em EPSG:4326, já que os originais misturam 4674, 4326 e 31984.
- **Extensão:** a do dado baixado (folha ou mapa inteiro), sem recorte pela área de estudo.
- **Pendências para o script:**
  - quais séries de `lito`/`estrutura` sobem (o milionésimo duplica a Geologia BA);
  - os nomes finais dos assets;
  - se as planilhas de Bambuí e Morro do Chapéu (C9) sobem.

### 10.3 Como rodar

```
cd src/dowloadQGIS
python exportar_acervo_gee.py
  Etapa 1 → prepara em /mnt/Almacen/shpIRECE/gee_upload/ (só local; dá para conferir no QGIS)
  Etapa 2 → pede bucket (gs://…), conta e se sobrescreve; sobe ao bucket e dispara as ingestões
  Etapa 3 → as duas
```

- **Bucket do Cloud Storage:** o GEE só ingere arquivos que estão em `gs://`. Precisa de um bucket em que a conta escolhida escreva e que o GEE consiga ler. O nome pode ficar fixo em `BUCKET` no topo do script.
- **Credenciais do Storage:** o script tenta primeiro as credenciais da própria conta do EE (a solkan1201 tem o escopo `devstorage`). Se não der, usa as do gcloud: `gcloud auth application-default login`. Hoje o gcloud está pedindo reautenticação.
- **Registro:** cada envio grava `registro_exportar_gee_<data>.csv` com o asset, o `gs://` e o id da task. A ingestão roda como task no GEE (ver `tools/revisarTasksXconta.py`).
- **Lito/estrutura 100k:** ficam de fora por padrão (`INCLUIR_100K = False`).
- **Duplicatas de fosfato:** ocorrências a menos de ~300 m contam como uma só (a mesma aparece como `P` numa série e `Fósforo` noutra). Lajedão I e II, a ~450 m uma da outra, continuam separadas.


---

## 11. Outros alvos além do fosfato — recomendações para tarefas futuras

Levantado em 2026-10-12 a partir dos dados já baixados.
- **Ocorrências:** todos os `recmin`/`ocormin`, sem duplicatas (~100 m + mesma substância), só os que caem na área de estudo: **4.180**, das quais só 42 são de fosfato.
- **Geoquímica:** conferidas as colunas das planilhas em `geoquimica/`.
- **Como contar de novo:** `ler_recmin()` + `sem_duplicatas()` de `exportar_acervo_gee.py`.

### 11.1 Referências disponíveis por material

| Material | Ocorrências na área | Até 150 km de Irecê | Minas / garimpos | Referência extra já baixada | Detecção pelo ASTER |
|---|---|---|---|---|---|
| **Calcário / dolomito** (corretivo agrícola, cimento) | 237 + mármore 124 + magnesita 45 | 54 + 26 + 0 | 37 / 72 (calcário) | **Calcários da Bahia (6719): 2.245 amostras de rocha com CaO, MgO, SiO₂, Al₂O₃, Fe₂O₃, K₂O, LOI** (+ 7420) | **Direta e forte:** CO₃ em ~2,33 µm (B08); calcita × dolomita no TIR (B13/B14). Mesmo hospedeiro do fosfato. |
| **Pb–Zn–F–Ba nos carbonatos do Salitre** | barita 218, fluorita 69, chumbo 49, zinco 9 | barita 50, chumbo 17, fluorita 6, zinco 3 | barita 21 minas | Pb, Zn, F, Ba, Mn, Cu em sedimento de corrente (4599 Morro do Chapéu; 4643/4646 Bambuí) e em rocha (6719, 7420) | **Indireta:** o mineral não aparece no espectro; o hospedeiro carbonático e as estruturas sim. Usa as covariáveis do fosfato. |
| **Talco / magnesita / vermiculita / amianto** (ultramáficas e magnesianas; Brumado, Campo Formoso) | 152 / 45 / 88 / 41 | ~0 | magnesita 28 minas; talco 8 | cromo: 102 ocorrências, 25 minas (Campo Formoso); Cr, Ni, MgO nas planilhas de rocha | **Direta:** Mg–OH em ~2,31–2,32 µm (B07/B08). Longe de Irecê. |
| **Ouro** (Jacobina, Chapada Diamantina) | 340 | 72 | 10 / 153 | ARIM Serra de Jacobina (3594); Au_ppm nas planilhas de rocha | **Indireta:** alteração hidrotermal (sericita/argila Al–OH em ~2,20 µm, B05/B06) e silicificação (índice de quartzo no TIR). |
| **Manganês** | 337 | 38 | 34 / 41 | Mn_ppm em sedimento e rocha | **Parcial:** óxidos de Mn escurecem o VNIR; sem feição diagnóstica no SWIR. |
| **Caulim / argila / pirofilita** | 26 / 19 / 8 | 3 / 1 / 2 | poucos | Al₂O₃ nas planilhas de rocha | **Direta e forte** (Al–OH, ~2,16–2,20 µm), mas há poucas referências. |
| **Urânio** (Lagoa Real / Caetité) | 30 | 0 | 1 mina | **eU do 1042** (mas Lagoa Real fica ao sul da cobertura do 1042) | ASTER não; gama sim. |
| Diamante | 254 | 161 | 1 / 230 | — | **Não:** aluvião / conglomerado do Tombador. |
| Quartzo / cristal de rocha | 552 | 175 | 9 / 298 | — | Fraca: veios pequenos para o pixel do TIR (90 m). |

### 11.2 Prioridade sugerida

1. **Calcário / dolomito: próximo alvo recomendado.**
   - Tem muitas referências: ocorrências e geoquímica com CaO/MgO para separar calcário de dolomito.
   - O ASTER detecta o carbonato diretamente.
   - Está na mesma bacia e no mesmo hospedeiro do fosfato.
   - Reaproveita o pipeline do MaxEnt quase inteiro, só trocando os pontos de presença.
2. **Pb–Zn–F–Ba (estilo MVT nos carbonatos do Salitre).**
   - Mesmo hospedeiro e mesmas covariáveis do fosfato, inclusive a distância a estruturas.
   - A geoquímica de sedimento de corrente ajuda a validar.
3. **Talco / magnesita** e **ouro**: só se o modelo for estendido além da Bacia de Irecê (Brumado / Campo Formoso / Jacobina).
4. **Caulim / argila**: espectralmente o melhor alvo do ASTER, mas há poucas referências. Pode servir para validar índices de argila mais do que como alvo de modelo.
5. **Não recomendados com ASTER:** diamante e quartzo/cristal. Os registros são quase todos garimpos, que refletem onde se garimpou, e o alvo não tem resposta espectral útil.

### 11.3 Tarefas futuras

- [x] **Leitura do Projeto Morro do Chapéu corrigida em `exportar_acervo_gee.py` (2026-10-12).**
  - Os nomes de campo agora são comparados sem acento e por prefixo (`DESCRI` pega `DESCRIÇÃ`; `RECURSO` pega o campo da folha Cabeça no Tempo).
  - Saíram da leitura `Mo_Chap_Amostragem_Rec_Min` (54 amostras de geoquímica, não ocorrências) e `P_Mo_Chapeu_OcorMin_ponto` (cópia exata, nas mesmas 49 coordenadas, de `Mo_Chap_OcorMin_Lito_RecMin`, mas com a descrição em branco).
  - Total na área: **4.077 ocorrências sem duplicatas**. Só 6 continuam sem substância, porque estão vazias na própria fonte (4 em Anagé, 2 em Barra–Oliveira dos Brejinhos). Fosfato continua com 42.
  - As contagens da tabela 11.1 foram feitas antes da correção, com 4.180 ocorrências. As diferenças são pequenas e não mudam a prioridade.
- [ ] **Asset de pontos por alvo:** gerar `ocorrencias_<alvo>` (calcário/dolomito, Pb–Zn–F–Ba, talco/magnesita, ouro) do mesmo jeito que `ocorrencias_fosfato`, com tabela de sinônimos de substância. Exemplos: `P`/Fósforo; `cc`/`ca`/Calcário; `mm`/Mármore; `Ba`/Barita; `F`/Fluorita; `Pb`; `Zn`; `tal`/Talco; `Mg`/Magnesita; `Au`/Ouro.
- [ ] **Geoquímica multielementar no GEE:** além do P₂O₅, subir CaO, MgO, SiO₂, Al₂O₃, Pb, Zn, F, Ba, Mn, Au, Cr, Ni de 6719/7420 (rocha) e o sedimento de corrente de 4599/4643/4646. Converter vírgula decimal e `<`.
- [ ] **Índices ASTER por alvo:** carbonato (B08), calcita × dolomita (TIR), Mg–OH (talco/magnesita), Al–OH (sericita/caulim), sílica (TIR). Validar com R6/R7 (afloramentos e litologia).
- [ ] **Fundo por alvo:** usar o target-group (os outros `recmin`) excluindo o alvo, como já previsto em R5.
- [ ] **Avaliar o viés:** a maioria dos registros de quartzo, diamante e gemas é de garimpo, então a distribuição reflete onde se garimpou. Considerar só minas e ocorrências descritas em campo, ou ponderar pelo status.
- [ ] **Dados do acervo ainda não baixados**, a avaliar conforme o alvo:
  - **Rochas Ornamentais** (6441): granito/mármore.
  - **Geodiversidade do Estado da Bahia** (2181, 1 GB): camada "Zonas de interesse mineral" e solos.
  - **Cartas de Anomalias FIOL** (32 folhas, ~5,5 GB): grids gama e magnetometria de alta resolução no sul da área, com cobertura de Brumado/Caetité, úteis para talco/magnesita e urânio.
  - **Aerolevantamentos 1023, 1027, 1115 e 3067:** cobertura gama do leste, do norte e do oeste da área, fora do 1042.

---

## 12. Passo a passo: credenciais do gcloud, bucket e envio ao GEE

> **Para a conta `solkan.cengine17@gmail.com` (`solkanCengine`, projeto `ee-solkancengine17`), seguir o guia completo [`GUIA_GCLOUD_BUCKET_GEE.md`](GUIA_GCLOUD_BUCKET_GEE.md).** Esta seção descreve a alternativa com a conta da Geodatin.

**Situação em 2026-10-12:**
- O gcloud está configurado com a conta `solkan@geodatin.com` e o projeto `geo-data-s`, mas as credenciais **expiraram** ("Reauthentication required").
- Essa é a mesma conta Google da conta EE `solkanGeodatin` (projeto `geo-data-s`), que escreve em `mapbiomas-arida/mine`.
- **Recomendação:** usar essa conta em tudo (gcloud, bucket e ingestão), para que o GEE leia o bucket sem precisar de permissão extra.

Os passos 1 e 2 abrem o navegador. Rodar no terminal; dentro do Claude Code, usar o prefixo `!`.

### Passo 1 — renovar o login do gcloud (CLI)
```bash
gcloud auth login solkan@geodatin.com
gcloud config set project geo-data-s
gcloud auth list                      # conferir: * solkan@geodatin.com
```

### Passo 2 — salvar as credenciais para o Python (Application Default Credentials)
O script usa a biblioteca `google-cloud-storage`. Quando as credenciais do EE não dão acesso ao Storage, ele usa estas:
```bash
gcloud auth application-default login
gcloud auth application-default set-quota-project geo-data-s
```
Elas ficam salvas em `~/.config/gcloud/application_default_credentials.json`. **É uma chave de acesso: não versionar nem compartilhar.**

### Passo 3 — conferir se o projeto pode ter bucket
```bash
gcloud billing projects describe geo-data-s     # billingEnabled: true
gcloud services enable storage.googleapis.com earthengine.googleapis.com
```
Sem faturamento ativo, o Cloud Storage não cria bucket. O custo aqui é mínimo: ~150 MB, que podem ser apagados depois da ingestão.

### Passo 4 — criar o bucket
O nome precisa ser único no mundo, só com letras minúsculas, números e hífen. Usar a região `US-CENTRAL1`, a mesma do Earth Engine.
```bash
gcloud storage buckets create gs://geodatin-acervo-cprm \
    --project=geo-data-s --location=US-CENTRAL1 --uniform-bucket-level-access
```

### Passo 5 — testar escrita e leitura
```bash
echo teste | gcloud storage cp - gs://geodatin-acervo-cprm/acervo_cprm/teste.txt
gcloud storage ls gs://geodatin-acervo-cprm/acervo_cprm/
gcloud storage rm gs://geodatin-acervo-cprm/acervo_cprm/teste.txt
```

### Passo 6 — (opcional) usar outra conta EE na ingestão
Se a ingestão for feita com outra conta (ex.: `caatinga04`), essa conta precisa **ler** o bucket:
```bash
gcloud storage buckets add-iam-policy-binding gs://geodatin-acervo-cprm \
    --member=user:<email-da-conta-caatinga04> --role=roles/storage.objectViewer
```

### Passo 7 — fixar o bucket no script (opcional)
Em `exportar_acervo_gee.py`:
```python
BUCKET = 'gs://geodatin-acervo-cprm/acervo_cprm'
```
Se ficar vazio, o script pergunta o bucket ao rodar.

### Passo 8 — preparar os arquivos e conferir
```bash
cd src/dowloadQGIS
python exportar_acervo_gee.py        # Etapa: 1
```
- Gera `/mnt/Almacen/shpIRECE/gee_upload/`: `1042_gama.tif`, `1042_mag.tif` e os vetores em `.zip`.
- Abrir no QGIS para conferir: as 7 bandas do gama, os pontos de fosfato com o campo `grupo`, e a litologia/estruturas em EPSG:4326.

### Passo 9 — enviar
```bash
python exportar_acervo_gee.py        # Etapa: 2
  Bucket ............ gs://geodatin-acervo-cprm/acervo_cprm
  Conta ............. 7 (solkanGeodatin)   ← mesma conta do gcloud/bucket
  Sobrescrever ...... n
```
- Copia cada arquivo para o bucket e dispara a ingestão em `projects/mapbiomas-arida/mine/acervo_cprm/`.
- O registro fica em `registro_exportar_gee_<data>.csv`.

### Passo 10 — acompanhar e limpar
- **Acompanhar as tasks:** `python ../tools/revisarTasksXconta.py` (conta `solkanGeodatin`, filtro vazio) ou a aba *Tasks* do Code Editor.
- **Conferir no GEE:** `ee.data.listAssets({'parent': 'projects/mapbiomas-arida/mine/acervo_cprm'})`.
- **Limpar o bucket**, depois que todas as tasks estiverem `COMPLETED`:
  ```bash
  gcloud storage rm -r gs://geodatin-acervo-cprm/acervo_cprm
  ```

### Problemas comuns

| Erro | Causa / solução |
|---|---|
| `Reauthentication required` / `RefreshError` | refazer os passos 1 e 2 |
| `403 … storage.buckets.create` | a conta não tem papel de *Storage Admin* no projeto, ou o faturamento está desligado (passo 3) |
| `409 … bucket already exists` | o nome já existe no mundo: escolher outro |
| Ingestão `FAILED: … cannot read gs://…` | a conta EE da ingestão não lê o bucket: usar a mesma conta do bucket, ou o passo 6 |
| `Asset already exists` | responder `s` em "sobrescrever", ou apagar o asset antes |
| Ingestão de tabela com erro de geometria | ver o shapefile em `gee_upload/` no QGIS (*Verificar validade*) e reenviar só ele |
