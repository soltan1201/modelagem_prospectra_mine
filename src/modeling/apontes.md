Pesquisei e não encontrei nenhum método publicado que faça exatamente isso, corrigir a saturação do VNIR do ASTER com Landsat. O trabalho mais próximo é o mosaico global de 15 m do ASTER (Gonzalez et al., 2019), que recupera B01/B02 saturadas a partir da B3N com polinômios de 1ª e 2ª ordem, usando MODIS como referência. Também há os ajustes entre sensores do tipo SBAF e as regressões ETM+ × ASTER banda a banda. A sua ideia é uma adaptação desses dois caminhos e é viável, mas tem um ponto crítico: medir a correção nos pixels saturados, onde não existe valor verdadeiro para comparar.

O que a coleção Landsat oferece
Bandas: LANDSAT/COMPOSITES/C02/T1_L2_32DAY tem blue, green, red, nir, swir1, swir2 e thermal, a 30 m.
Período: de 1984 a 2026, o que cobre todo o período ASTER (2000–2020).
Tipo de dado: é reflectância de superfície (Level-2). O ASTER do pipeline está em reflectância TOA, e a regressão vai ter que absorver essa diferença de atmosfera.
Composição: cada composto de 32 dias usa o pixel mais recente, não a mediana. A mediana do semestre seria feita sobre os ~6 compostos do período, e isso é fácil no GEE.
Sensores: o L7 é excluído depois de 2017 e o L8 antes de 2013-05. Em 2012 sobra praticamente só o L7 com SLC-off, que tem listras; a mediana atenua.
Proposta: correção por regressão ASTER ↔ Landsat
Pares por semestre e célula. Para cada cena ASTER (TOA × 1000, banda saturado), usar a mediana Landsat dos compostos de 32 dias do semestre. Melhor ainda é usar só os compostos a ±1 mês da data da cena, para reduzir a diferença de tempo.

Grade comum de 30 m. Agregar o ASTER de 15 para 30 m (reduceResolution, média 2×2) alinhado ao Landsat. Um pixel de 30 m conta como saturado se algum dos 4 sub-pixels estiver saturado. Agregar também reduz o erro de registro entre os dois sensores.

Amostras de treino. Usar só pixels não saturados e estáveis, descartando:

água;
mudanças de cobertura entre as datas (|NDVI_ASTER − NDVI_Landsat| alto).
Forçar amostras na faixa de brilho logo abaixo do teto (DN 200–253), que é a região que mais informa sobre o que existe acima dele.

Modelo por banda (B01 e B02). Prever ρ_B01 e ρ_B02 a partir de Landsat green, red, nir e da própria B3N do ASTER, que raramente satura. Testar em ordem de complexidade:

linear, no estilo SBAF;
polinômio de 2ª ordem, como no mosaico global;
smileRandomForest em modo regressão.
O random forest não extrapola além do maior valor visto no treino. Como os pixels saturados estão justamente acima da faixa de treino, a regressão linear ou polinomial tende a funcionar melhor ali.

Medição da melhora, com duas formas de validar:

Pseudo-saturação: baixar o limiar artificialmente, tratando DN ≥ 200 como "saturado". O modelo é treinado com DN < 200 e prevê a faixa 200–253, onde o valor real é conhecido. Isso dá RMSE, viés e R² da extrapolação, por banda, contra o valor cortado atual (baseline = DN 254 fixo).
Verdade independente pelo ganho: a propriedade GAIN_SETTING_B01/B02 diz o ganho de cada cena. Cenas em ganho normal do mesmo alvo não saturam. Comparar o valor corrigido com uma observação em ganho normal do mesmo pixel em outra data (descontando a sazonalidade) valida direto nos pixels saturados.
Aplicação. Nos pixels com saturado = 1:

substituir B01/B02 pelo valor previsto, com a restrição ρ_previsto ≥ ρ_cortado, já que o valor real só pode ser maior;
gravar uma banda corrigido (0/1) para rastrear o que foi alterado;
aplicar o valor de 30 m nos 4 sub-pixels de 15 m, ou redistribuir pela textura da B3N de 15 m numa versão mais refinada.
Métricas finais: RMSE/viés por banda antes × depois, recuperação da cauda do histograma e coerência dos índices que dependem de B01/B02 (Ferric_Fe = B02/B01, NDWI).

Riscos
Landsat também satura em alvo muito claro: TM e ETM+ são de 8 bits; o OLI, de 12 bits, quase não satura. Os compostos não trazem a banda de QA de saturação, então convém descartar no treino pixels Landsat no topo da distribuição.
Diferença de tempo: a mediana do semestre pode ser de outra condição de vegetação. Usar a janela de ±1 mês e o filtro de NDVI reduz isso.
Extrapolação: é a maior incerteza, e por isso o teste de pseudo-saturação é o que diz se vale a pena aplicar a correção.
Próximo passo sugerido
Um piloto em uma célula e um semestre com bastante saturação, por exemplo uma célula de calcário de Irecê em 2010-S2, que tem muitas cenas zero-cloud. O piloto mediria:

a fração de pixels saturados;
o erro da pseudo-saturação para os modelos linear, polinomial e RF;
uma comparação visual.
Posso montar esse piloto em JS ou em Python. Qual prefere?

Sobre o visualizar_zero_cloud.js: vi que você comentou cenas na lista (2000_1, 2001_1, 2001_2, 2006_1, 2007_2, 2015_1, 2019_1, 2020_1) e marcou 20100703125452 como "remover de vez". Quando fechar a lista, eu passo essas remoções para o list_zero_cloud.txt.


Sources:

Global 15-Meter Mosaic Derived from Simulated True-Color ASTER Imagery (Gonzalez et al., 2019) · registro DOAJ · MDPI
Landsat C2 T1 L2 32-Day Composite — catálogo GEE
Applications of Spectral Band Adjustment Factors (SBAF) for Cross-Calibration — USGS · pubs.usgs.gov · Chander et al. 2010
Cross-calibration of AM constellation sensors (Meyer & Chander) — USGS · SPIE
Xu & Zhang 2011 — ASTER × ETM+ NDVI/SAVI · resumo Worldwide Science
Cross-Calibration between ASTER and MODIS VNIR (Obata et al., 2017) · correção
Comparing ASTER and Landsat 7 ETM+ for Change Detection — ASPRS 2006
Global cross-calibration of Landsat spectral mixture models (arXiv)
Onboard electrical calibration of the ASTER VNIR — SPIE
AST L1T Quick Reference Guide