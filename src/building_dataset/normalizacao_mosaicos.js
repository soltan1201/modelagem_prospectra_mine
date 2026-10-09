// =================================================================
// --- LABORATÓRIO DE NORMALIZAÇÃO COM MOSAICO EM TEMPO REAL ---
// --- VERSÃO: ESTATÍSTICA ROBUSTA (MEDIANA/PERCENTIL + EPSILON) ---
// =================================================================

// 1. PARÂMETROS E ÁREA DE ESTUDO
var list_coord = [
    [-44.460483519818084,-14.790793773956535],
    [-39.203525512005584,-14.790793773956535],
    [-39.203525512005584,-9.416525285369607],
    [-44.460483519818084,-9.416525285369607],
    [-44.460483519818084,-14.790793773956535]
];
var area_estudo = ee.Geometry.Polygon(list_coord);

// 2. CONSTRUÇÃO DO GABARITO (MESTRE 2000-2008) EM MEMÓRIA
var id_asset_m1 = 'projects/mapbiomas-arida/mosaic_aster';
var colecao_m1 = ee.ImageCollection(id_asset_m1);

var mosaico_mestre = colecao_m1.qualityMosaic('quality').clip(area_estudo);

// 3. CARREGAR O SEMESTRE ALVO (RUIM)
var id_alvo = 'projects/mapbiomas-arida/mosaic_aster/ASTER_QualityMosaic_2003_semestre_2_INT16'; 
var img_alvo = ee.Image(id_alvo);

// 4. MOTOR DE NORMALIZAÇÃO ROBUSTA (Z-SCORE OTIMIZADO)
function normalizarImagem(imgTarget, imgRef, roi, scale) {
    var imgTargetFloat = imgTarget.divide(10000).float();
    var imgRefFloat = imgRef.divide(10000).float();
    
    // Passo A: Identificar quais bandas sobreviveram aos filtros (Prevenção de Nulls)
    var testTarget = imgTargetFloat.reduceRegion({
        reducer: ee.Reducer.first().unweighted(), geometry: roi, scale: scale, maxPixels: 1e10, bestEffort: true
    });
    var testRef = imgRefFloat.reduceRegion({
        reducer: ee.Reducer.first().unweighted(), geometry: roi, scale: scale, maxPixels: 1e10, bestEffort: true
    });
    
    var bandasValidasAlvo = ee.Dictionary(testTarget).keys();
    var bandasValidasRef = ee.Dictionary(testRef).keys();
    
    // Truque lógico para criar a interseção no Earth Engine
    var bandasExclusivasAlvo = bandasValidasAlvo.removeAll(bandasValidasRef);
    var bandasParaCorrigir = bandasValidasAlvo.removeAll(bandasExclusivasAlvo);
    
    // Passo B: Extrair Estatísticas Robustas Simultâneas (Percentis 5%, 50% e 95%)
    var pTarget = imgTargetFloat.select(bandasParaCorrigir).reduceRegion({
        reducer: ee.Reducer.percentile([5, 50, 95]), geometry: roi, scale: scale, maxPixels: 1e10, bestEffort: true
    });
    var pRef = imgRefFloat.select(bandasParaCorrigir).reduceRegion({
        reducer: ee.Reducer.percentile([5, 50, 95]), geometry: roi, scale: scale, maxPixels: 1e10, bestEffort: true
    });
    
    // Passo C: Aplicar a Transformação Afim
    var matchedBands = bandasParaCorrigir.map(function(bName) {
        bName = ee.String(bName);
        var targetBand = imgTargetFloat.select(bName);
        
        // Coleta Mediana (P50) e os extremos (P5 e P95)
        var tMedian = ee.Number(pTarget.get(bName.cat('_p50')));
        var tP5     = ee.Number(pTarget.get(bName.cat('_p5')));
        var tP95    = ee.Number(pTarget.get(bName.cat('_p95')));
        
        var rMedian = ee.Number(pRef.get(bName.cat('_p50')));
        var rP5     = ee.Number(pRef.get(bName.cat('_p5')));
        var rP95    = ee.Number(pRef.get(bName.cat('_p95')));
        
        // Pseudo-Desvio Padrão Robusto: Em uma distribuição normal, (P95 - P5) / 3.29 = Desvio Padrão
        var tStdRobust = tP95.subtract(tP5).divide(3.29);
        var rStdRobust = rP95.subtract(rP5).divide(3.29);
        
        // Proteção Numérica Epsilon: Evita divisão por zero em áreas perfeitamente homogêneas
        var eps = ee.Number(1e-6);
        var tStdSafe = tStdRobust.max(eps);
        
        // Fórmula Final: ((Alvo - Mediana_Alvo) / Std_Alvo_Seguro) * Std_Ref + Mediana_Ref
        var normalizado = targetBand.subtract(tMedian)
                                    .divide(tStdSafe)
                                    .multiply(rStdRobust)
                                    .add(rMedian)
                                    .rename([bName]);
                                    
        return normalizado;
    });
    
    // Remonta a imagem devolvendo as bandas originais intactas onde a correção não atuou
    var imgsCorrigidas = ee.ImageCollection(matchedBands).toBands().rename(bandasParaCorrigir);
    return imgTargetFloat.addBands(imgsCorrigidas, null, true); 
}

// =================================================================
// 5. EXECUÇÃO CLÍNICA NA TELA
// =================================================================
var visRGB = { bands: ['B3N', 'B02', 'B01'], min: 0.05, max: 0.35, gamma: 1.4 };
var visSWIR = { bands: ['B04', 'B06', 'B3N'], min: 0.1, max: 0.5, gamma: 1.2 };

if (typeof roi_ancora !== 'undefined') {
    
    var img_alvo_corrigida = normalizarImagem(img_alvo, mosaico_mestre, roi_ancora, 15);

    Map.centerObject(roi_ancora, 12); 
    Map.addLayer(mosaico_mestre.divide(10000).float(), visRGB, '1. GABARITO MESTRE (Ao Vivo)', false);
    Map.addLayer(img_alvo.divide(10000).float(), visRGB, '2. ALVO ORIGINAL', false);
    Map.addLayer(img_alvo_corrigida, visRGB, '3. ALVO NORMALIZADO ROBUSTO (RGB)', true);
    Map.addLayer(img_alvo_corrigida, visSWIR, '4. ALVO NORMALIZADO ROBUSTO (SWIR)', false);
    
    print('Normalização robusta (Mediana/Percentis + Epsilon) executada com sucesso.');
    
} else {
    Map.centerObject(area_estudo, 8);
    Map.addLayer(mosaico_mestre.divide(10000).float(), visRGB, '1. GABARITO (Procure a âncora aqui)', false);
    Map.addLayer(img_alvo.divide(10000).float(), visRGB, '2. ALVO ORIGINAL', true);
    print("PAUSA: Desenhe a roi_ancora estritamente sobre a mancha vegetal (para secar) ou sobre o afloramento (para calibrar). Aperte Run novamente.");
}

// =================================================================
// 6. PROVAS NUMÉRICAS E INTEGRIDADE ESPECTRAL (GRÁFICOS NATIVOS)
// =================================================================

if (typeof roi_ancora !== 'undefined') {
    
    // Escolhemos a banda B3N (Infravermelho Próximo) pois é ela que controla o "Vermelho"
    var bandaTeste = 'B3N';

    // 1. Preparar as 3 fatias de dados (Original, Referência e Normalizada)
    var imgOriginal = img_alvo.divide(10000).float().select([bandaTeste]).rename(['1_Alvo_Original']);
    var imgReferencia = mosaico_mestre.divide(10000).float().select([bandaTeste]).rename(['2_Gabarito_Mestre']);
    var imgNormalizada = img_alvo_corrigida.select([bandaTeste]).rename(['3_Alvo_Normalizado']);

    // Unir as 3 em uma única imagem multi-banda para o gerador de gráficos
    var imagemComparativa = ee.Image([imgOriginal, imgReferencia, imgNormalizada]);

    // =================================================================
    // PROVA 1: O HISTOGRAMA (Eficiência da Normalização Z-Score)
    // =================================================================
    var histograma = ui.Chart.image.histogram({
        image: imagemComparativa,
        region: roi_ancora,
        scale: 15,
        minBucketWidth: 0.005
    }).setOptions({
        title: 'Prova 1: Histograma da Banda B3N (Infravermelho) dentro da ROI',
        hAxis: {title: 'Reflectância', titleTextStyle: {italic: false, bold: true}},
        vAxis: {title: 'Contagem de Pixels', titleTextStyle: {italic: false, bold: true}},
        colors: ['red', 'blue', 'green'],
        interpolateNulls: true
    });
    
    print('GERANDO GRÁFICOS DE VALIDAÇÃO...');
    print(histograma);

    // =================================================================
    // PROVA 2: GRÁFICO DE DISPERSÃO (Prova de Integridade Espectral)
    // =================================================================
    // Sorteamos 500 pixels aleatórios dentro do polígono para provar a linearidade
    var amostras = ee.Image([imgOriginal.rename('Original'), imgNormalizada.rename('Normalizado')])
                     .sample({region: roi_ancora, scale: 15, numPixels: 500});

    var scatter = ui.Chart.feature.byFeature(amostras, 'Original', 'Normalizado')
        .setChartType('ScatterChart')
        .setOptions({
            title: 'Prova 2: Integridade Espectral (Alvo Original vs. Normalizado)',
            hAxis: {title: 'Reflectância Original (B3N)'},
            vAxis: {title: 'Reflectância Normalizada (B3N)'},
            colors: ['purple'],
            pointSize: 2,
            trendlines: { 0: { showR2: true, visibleInLegend: true, color: 'black' } }
        });
        
    print(scatter);
}