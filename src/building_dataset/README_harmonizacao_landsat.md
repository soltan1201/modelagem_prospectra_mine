# Harmonização do mosaico sazonal seco ASTER contra o Landsat

Scripts: `exportar_landsat_sazonal.py` (referência) e `harmonizar_sazonal_landsat.py`
(harmonização). Visualizador: `src/show_mosaics/visualizar_harmonizado_sazonal.js`
(`HARM = 'landsat_razao'` ou `'landsat_residuo'`; caixa "Landsat 32 dias à esquerda").

## Por que trocar o MODIS pelo Landsat

- **Referência contínua.** A mediana Landsat 32 dias dos meses secos (jun–out) do mesmo
  período não tem emendas entre cenas. O mosaico ASTER tem.
- **Resolução.** O Landsat tem 30 m; o MODIS tem 500 m. A comparação fica quase pixel a
  pixel, e cada pedaço de cena ganha a sua própria transformação, inclusive na coluna 1
  (L1C1 e L2C1), que antes dependia da ponte pelo MODIS mediano.
- **Por que a reta não bastava.** Nas cenas de 2000 ainda verdes, a reta do B3N contra o
  MODIS tinha R² ≈ 0 (`coef_rel_A.csv`). A relação com a referência não é linear. O
  casamento de histograma não presume forma nenhuma.

## Ideia

O mosaico é a soma de pedaços: cada pedaço é uma cena (`cena_t` do mapa de cenas) dentro
de uma célula. Cada pedaço recebe uma função valor → valor que leva o histograma dele ao
histograma do Landsat na **mesma área**. Como todos os pedaços vão para a mesma referência
contínua, os degraus entre cenas somem. A função é aplicada aos pixels de 15 m, então a
textura e o detalhe do ASTER continuam.

Para não destruir o sinal mineral, só a **banda guia** de cada grupo é casada com o Landsat.
As outras bandas seguem a guia:

| Grupo | Guia ↔ Landsat | Outras bandas |
|---|---|---|
| VNIR | B01 ↔ green, B02 ↔ red, B3N ↔ nir (cada uma casada) | — |
| SWIR (só A) | B04 ↔ swir1 | B05–B09 seguem o B04 |
| TIR | B13 ↔ temperatura de superfície Landsat (`st`) | B10, B11, B12, B14 seguem o B13 |

Por que o B04 é a guia do SWIR: é a banda mais próxima do swir1 (1,6 µm) e é a **fonte** do
crosstalk do ASTER, ou seja, a menos contaminada por ele (Tonooka e Iwasaki 2004).

Por que o TIR segue o B13: a diferença de TIR entre cenas é quase toda temperatura (data,
umidade, vento) e desloca as cinco bandas juntas. O sinal mineral (quartzo, carbonato)
está na forma da emissividade, isto é, nas razões entre as bandas.

## Passo a passo (`harmonizar_sazonal_landsat.py`)

### 1. Pixels estáveis ("pouca mudança")

Em Caatinga, o EVI responde melhor que o NDVI. O ASTER não tem banda azul, então usamos o
**EVI2** (Jiang et al. 2008):

```
EVI2 = 2,5 (NIR − R) / (NIR + 2,4 R + 1)
```

Para cada pixel de 60 m: d = EVI2_ASTER − EVI2_Landsat. **Estável** é o pixel com
|d − mediana(d) do pedaço| ≤ quantil `FRAC_ESTAVEL` (0,6) do pedaço. Usamos a mediana do
pedaço porque o deslocamento geral da cena é justamente o que queremos corrigir; o que se
descarta são os pixels que mudaram **mais que o resto da cena** (parcela verde, queimada,
lagoa). É a mesma lógica de normalização por pixels sem mudança do IR-MAD (Canty e Nielsen
2008), em versão simples. As curvas saem só dos estáveis e valem para todos os pixels.

### 2. VNIR: casamento de histograma por percentis

Percentis 1, 3, …, 99 do ASTER e do Landsat nos estáveis do pedaço →
`image.interpolate(perc_ASTER, perc_Landsat, 'clamp')`. Abaixo do 1º e acima do 99º
percentil, as caudas seguem a inclinação entre os percentis 5 e 95, para que anomalias
claras ou escuras não sejam cortadas.

A lógica é a do `src/show_mosaics/intercalibracao.js` (Planet): CDF da imagem → CDF da
referência. Lá, a tabela de consulta era feita com duas Random Forest; aqui ela é feita
com percentis e interpolação, o que é exato, determinístico e mais leve.

### 3. SWIR: a guia e as bandas que a seguem

O B04 é casado com o swir1 como no item 2, dando o B04'. Para B05–B09 há dois modos (o
script pergunta qual usar):

- **razao**: `Bk' = Bk × B04'/B04`. Todos os índices de razão (Al-OH (B05+B07)/B06,
  carbonato (B06+B09)/B08, …) ficam idênticos ao original. Não corrige diferenças aditivas
  entre cenas (atmosfera, crosstalk).
- **residuo**:
  `Bk' = [Bk − (a_p + b_p·B04) + (a_ref + b_ref·B04)] × B04'/B04`.
  - (a_p, b_p) é a reta Theil-Sen de Bk contra B04 nos pixels estáveis do pedaço.
  - (a_ref, b_ref) é a mediana dessas retas entre todos os pedaços do período.
  - O termo `Bk − (a_p + b_p·B04)` é o **resíduo**, o quanto o pixel foge do
    comportamento de fundo da cena. É o sinal mineral, e ele é preservado.
  - Risco: um pedaço dominado por uma litologia de verdade diferente tem uma relação de
    fundo real, que seria "corrigida". Para conferir, olhe se o índice Al-OH fica
    contínuo nas emendas e se as anomalias conhecidas sobrevivem.

### 4. TIR

Em cada pedaço, ajusta-se a reta quantil-quantil DN13 = a + b·st (percentis do B13 contra
percentis da `st` do Landsat). A referência do período é a mediana de (a, b) entre os
pedaços. O alvo do pedaço é `a_ref + b_ref·st_percentis`, e o B13 do pedaço é casado com
esse alvo. B10, B11, B12 e B14 seguem o B13, como no SWIR.

Pedaço sem `st` suficiente (menos de 50% dos estáveis): o TIR fica como no original.

**Hipótese a testar** (etapa 4 do script): o ajuste do TIR acompanha o do SWIR ou o desvio
de EVI2? Uma cena mais úmida tende a ser mais escura no SWIR e mais fria no TIR. A etapa 4
calcula a correlação de Spearman entre pedaços. Se |ρ| ≥ 0,5, dá para prever a correção do
TIR a partir do SWIR/EVI2 onde não há `st`.

### 5. Reservas

- Pedaço com menos de `N_MIN_PECA` (300) pixels estáveis usa a curva da célula inteira
  (todos os pedaços juntos).
- Pixel fora do mapa de cenas também usa a curva da célula.
- A banda `hm` marca o que foi usado: 1 = curva do pedaço, 2 = curva da célula.

## Como rodar

1. `python exportar_landsat_sazonal.py`
   - Gera `Landsat_sazonal/Landsat_sazonal_<P>_seco_mediana_<LxCy>`: uma task por período × célula, que rodam em paralelo.
   - Antes de exportar, confere se a composição vem em reflectância 0–1.
2. `python harmonizar_sazonal_landsat.py`
   1. **Etapa 1:** amostras de todas as células, baixadas direto como CSV, sem Drive.
      Ficam em `src/Dados/ASTER_sazonal_harmonizacao/landsat/amostras/`.
   2. **Etapa 2:** curvas, calculadas localmente e salvas em `landsat/curvas_<P>.json`.
   3. **Etapa 3:** aplicar, exportando para `ASTER_sazonal_hm_razao` ou
      `ASTER_sazonal_hm_residuo`. Para comparar os dois modos, rode duas vezes.
   4. **Etapa 4:** relatório por célula, ganhos por banda e o teste TIR × SWIR.
3. Visualizador:
   1. `HARM = 'landsat_razao'` (ou `'landsat_residuo'`).
   2. Marque a caixa do Landsat.
   3. Compare VNIR e EVI2, depois Al-OH e carbonato nas emendas.

## Validação

1. **Degrau nas emendas:** diferença de um lado e do outro das bordas entre cenas
   (faixa de cerca de 300 m), antes e depois, por banda e para EVI2, Al-OH e carbonato.
2. **Concordância com o Landsat:** R² por célula, no VNIR e no B04.
3. **Pontos de referência do fosfato:** os índices não podem perder o contraste.
4. **Teste final:** R² da estimativa de SWIR com RF/GTB, comparando o original com o
   harmonizado.

## Referências

- Canty, M. J.; Nielsen, A. A. (2008). Automatic radiometric normalization of
  multitemporal satellite imagery with the iteratively re-weighted MAD transformation.
  *Remote Sensing of Environment*, 112, 1025–1036.
  <https://www2.imm.dtu.dk/pubdb/edoc/imm5362.pdf>
- Helmer, E. H.; Ruefenacht, B. (2005). Cloud-free satellite image mosaics with regression
  trees and histogram matching. *Photogrammetric Engineering & Remote Sensing*, 71(9),
  1079–1089. <https://research.fs.usda.gov/treesearch/30086>
- Hewson, R. D.; Cudahy, T. J.; Mizuhiko, S.; Ueda, K.; Mauger, A. J. (2005). Seamless
  geological map generation using ASTER in the Broken Hill-Curnamona province of
  Australia. *Remote Sensing of Environment*, 99(1–2), 159–172. (A lista de autores veio
  de uma citação secundária; confira na página da revista.)
  <https://catalog.sarig.sa.gov.au/document/mesac20907>
- Geoscience Australia: mapas ASTER nacionais, com calibração cruzada das cenas por
  ganhos e offsets para uma resposta comum.
  <https://ecat.ga.gov.au/geonetwork/srv/api/records/c3e62161-0fe4-2149-e044-00144fdd4fa6>
- Tonooka, H.; Iwasaki, A. (2004). Improvement of ASTER/SWIR crosstalk correction.
  *Proc. SPIE* 5234. <https://doi.org/10.1117/12.511811>
- Jiang, Z.; Huete, A. R.; Didan, K.; Miura, T. (2008). Development of a two-band enhanced
  vegetation index without a blue band. *Remote Sensing of Environment*, 112, 3833–3845.
  (Citada de memória; não foi conferida em busca.)
