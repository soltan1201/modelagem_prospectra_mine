// ------- 1. PARÂMETROS E CONFIGURAÇÃO -------
var vis_raw = {
    mosaico_clean: {
        bands: ['B3N', 'B02', 'B01'], 
        min: 0.05, max: 0.3
    },
    mosaic_raw: {
        bands: ['B3N', 'B02', 'B01'], 
        min: 25, max: 150
    }
};

// Área de estudo
var list_coord = [
    [-44.460483519818084,-14.790793773956535],
    [-39.203525512005584,-14.790793773956535],
    [-39.203525512005584,-9.416525285369607], 
    [-44.460483519818084,-9.416525285369607], 
    [-44.460483519818084,-14.790793773956535]
];
var area_estudo = ee.Geometry.Polygon(list_coord);
Map.addLayer(area_estudo, {}, 'Área de Estudo');

// Bandas para processamento
var bandas_para_asset = [
    'B01', 'B02', 'B3N',           // VNIR 15m
    'B04', 'B05', 'B06', 'B07', 'B08', 'B09'  // SWIR 30m
];

// ------- 2. CARREGAR MOSAICOS SEMESTRAIS -------
var id_asset_mosaicos = 'projects/ee-solkancengine17/assets/prospectra/mosaic_aster';
var mosaicos_sem = ee.ImageCollection(id_asset_mosaicos);

print('📊 Total de mosaicos semestrais:', mosaicos_sem.size());

// 👉 Verificar propriedades disponíveis
var props = mosaicos_sem.first().propertyNames();
print('📋 Propriedades disponíveis:', props);

// ------- 3. FUNÇÃO PARA CRIAR TIMESTAMP A PARTIR DAS PROPRIEDADES -------
function adicionarTimestamp(img) {
    var year = ee.Number(img.get('year'));
    var semestre = ee.Number(img.get('semestre'));
    
    // Calcula o mês baseado no semestre (1 = Jan, 2 = Jul)
    var month = ee.Number(ee.Algorithms.If(
        semestre.eq(1), 1,  // Primeiro semestre = Janeiro
        7                    // Segundo semestre = Julho
    ));
    
    // Cria data (dia 1 do mês)
    var data = ee.Date.fromYMD(year, month, 1);
    
    // Adiciona timestamp como banda e propriedade
    var timestamp = data.millis();
    
    return img
        .set('system:time_start', timestamp)  // 👈 Adiciona propriedade temporal
        .set('timestamp', timestamp)
        .set('data_formatada', data.format('YYYY-MM-dd'));
}

// ------- 4. APLICAR TIMESTAMP A TODOS OS MOSAICOS -------
var mosaicos_com_data = mosaicos_sem.map(adicionarTimestamp);

print('✅ Timestamps adicionados aos mosaicos');
print('📅 Datas:', mosaicos_com_data.aggregate_array('data_formatada'));

// ------- 5. FUNÇÃO PARA CALCULAR QUALITY SCORE (MULTI-CRITÉRIO) -------
function calcularQualityScore(img) {
    // Converte para float (já que o asset é INT16)
    var fator_escala = ee.Number(img.get('SCALE_FACTOR')).divide(10000); // Se existir
    
    // Bandas necessárias (converte para reflectância 0-1)
    var green = img.select('B01').float().divide(10000);
    var red = img.select('B02').float().divide(10000);
    var nir = img.select('B3N').float().divide(10000);
    var swir1 = img.select('B04').float().divide(10000);
    
    // ==== CRITÉRIO 1: BRILHO (EVITA SOMBRAS) ====
    var score_brightness = nir.unitScale(0.05, 0.35).clamp(0, 1);
    
    // ==== CRITÉRIO 2: TEXTURA (EVITA BORDAS DE NUVEM) ====
    var nir_byte = nir.unitScale(0, 0.4).multiply(255).toByte();
    var glcm = nir_byte.glcmTexture({size: 1});
    var contrast = glcm.select('B3N_contrast');
    var score_texture = contrast.unitScale(0, 255).multiply(-1).add(1).clamp(0, 1);
    
    // ==== CRITÉRIO 3: VEGETAÇÃO (EVITA NUVENS) ====
    var ndvi = nir.subtract(red).divide(nir.add(red));
    var score_ndvi = ndvi.unitScale(0, 0.5).clamp(0, 1)
        .multiply(ndvi.multiply(-1).add(0.8).clamp(0, 1));
    
    // ==== CRITÉRIO 4: CONSISTÊNCIA SWIR ====
    var score_swir = swir1.gt(0.05).and(swir1.lt(0.5))
        .multiply(0.5).add(0.5);
    
    // ==== CRITÉRIO 5: VARIÂNCIA LOCAL ====
    var variance = nir.reduceNeighborhood({
        reducer: ee.Reducer.stdDev(),
        kernel: ee.Kernel.square(2, 'pixels')
    });
    var score_variance = variance.unitScale(0, 0.05).multiply(-1).add(1).clamp(0, 1);
    
    // ==== SCORE FINAL (PONDERADO) ====
    var quality = score_brightness.multiply(0.30)
        .add(score_texture.multiply(0.25))
        .add(score_ndvi.multiply(0.20))
        .add(score_swir.multiply(0.15))
        .add(score_variance.multiply(0.10))
        .rename('quality_score');
    
    return img.addBands(quality);
}



// ------- 6. FUNÇÃO PARA DETECTAR E MASCARAR PIXELS RUINS -------
function detectarPixelsRuins(img) {
    var green = img.select('B01').float().divide(10000);
    var nir = img.select('B3N').float().divide(10000);
    var swir1 = img.select('B04').float().divide(10000);
    
    // Detecta sombras residuais
    var shadow_residual = nir.lt(0.06).and(green.lt(0.08));
    
    // Detecta nuvens residuais
    var cloud_residual = green.gt(0.35).and(swir1.gt(0.4));
    
    // Detecta bordas de nuvem
    var nir_byte = nir.unitScale(0, 0.4).multiply(255).toByte();
    var glcm = nir_byte.glcmTexture({size: 1});
    var high_contrast = glcm.select('B3N_contrast').gt(180);
    
    // Detecta pixels saturados
    var saturated = green.gte(0.5).or(nir.gte(0.5)).or(swir1.gte(0.6));
    
    // Máscara de pixels ruins
    var bad_pixel = shadow_residual.or(cloud_residual).or(high_contrast).or(saturated);
    bad_pixel = bad_pixel.focal_max({radius: 30, units: 'meters'});
    
    return img.updateMask(bad_pixel.not());
}

// ------- 7. APLICAR QUALITY SCORE A TODOS OS MOSAICOS -------
// var mosaicos_com_quality = mosaicos_com_data.map(function(img) {
//     // Converte para float
//     var img_float = img.float();
    
//     // Calcula quality score
//     img_float = calcularQualityScore(img_float);
    
//     // Adiciona máscara de pixels ruins
//     img_float = detectarPixelsRuins(img_float);
    
//     // Adiciona banda de ano/semestre para rastreabilidade
//     var year_semestre = ee.Image.constant(
//         ee.Number(img.get('year')).multiply(10).add(ee.Number(img.get('semestre')))
//     ).rename('year_semestre').toInt16();
    
//     // Calcula qualidade média para metadados
//     var quality_mean = img_float.select('quality_score').reduceRegion({
//         reducer: ee.Reducer.mean(),
//         geometry: area_estudo,
//         scale: 90,
//         maxPixels: 1e9
//     }).get('quality_score');
    
//     return img_float
//         .addBands(year_semestre)
//         .set('quality_mean', quality_mean)
//         .set('year', img.get('year'))
//         .set('semestre', img.get('semestre'))
//         .set('bloco', img.get('bloco'));
// });

// ------- 5. APLICAR QUALITY SCORE (SEM MÁSCARA DE PIXELS RUINS) -------
var mosaicos_com_quality = mosaicos_com_data.map(function(img) {
    // Converte para float
    var img_float = img.float();
    
    // Calcula quality score
    img_float = calcularQualityScore(img_float);
    
    // 👉 REMOVIDO: detectarPixelsRuins() - Estava mascarando demais!
    // img_float = detectarPixelsRuins(img_float);
    
    // Adiciona banda de ano/semestre
    var year_semestre = ee.Image.constant(
        ee.Number(img.get('year')).multiply(10).add(ee.Number(img.get('semestre')))
    ).rename('year_semestre').toInt16();
    
    return img_float
        .addBands(year_semestre)
        .set('year', img.get('year'))
        .set('semestre', img.get('semestre'))
        .set('bloco', img.get('bloco'));
});

print('✅ Quality score calculado para todos os mosaicos');

// ------- 8. CRIAR MOSAICO FINAL (QUALITY MOSAIC MULTI-TEMPORAL) -------
var mosaico_final = mosaicos_com_quality.qualityMosaic('quality_score');

// Clip para área de estudo
mosaico_final = mosaico_final.clip(area_estudo);

print('🎯 Mosaico final criado!');

// ------- 9. PÓS-PROCESSAMENTO: CONVERTER PARA REFLECTÂNCIA -------
// Converte de volta para reflectância (0-1)
var bandas_reflectancia = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09'];

var reflectancia = mosaico_final.select(bandas_reflectancia)
    .divide(10000)
    .float()
    .rename(bandas_reflectancia);

// Mantém as outras bandas
var outras_bandas = mosaico_final.select(['quality_score', 'year_semestre']);

mosaico_final = reflectancia.addBands(outras_bandas);



// ------- 11. BANDA DE FREQUÊNCIA -------
var frequency = mosaicos_com_quality.select('quality_score')
    .reduce(ee.Reducer.count())
    .rename('frequency')
    .clip(area_estudo)
    .toInt16();

mosaico_final = mosaico_final.addBands(frequency);

// ------- 12. VISUALIZAÇÃO -------
Map.centerObject(area_estudo, 8);

// RGB Natural
Map.addLayer(mosaico_final, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35,
    gamma: 1.4
}, '🌟 MOSAICO FINAL - RGB Natural');

// Falsa-cor Geológica
Map.addLayer(mosaico_final, {
    bands: ['B04', 'B06', 'B3N'],
    min: 0.1, max: 0.5,
    gamma: 1.2
}, '🌟 MOSAICO FINAL - SWIR Geologia');

// Índice de Argila
Map.addLayer(mosaico_final.select('AlOH_Clay'), {
    min: 0.8, max: 1.4,
    palette: ['blue', 'cyan', 'green', 'yellow', 'red']
}, '🌟 Al-OH (Argilas)');

// Frequência de observações
Map.addLayer(mosaico_final.select('frequency'), {
    min: 1, max: mosaicos_sem.size(),
    palette: ['red', 'orange', 'yellow', 'green', 'darkgreen']
}, '📊 Frequência (Nº de mosaicos)');

// Quality Score
Map.addLayer(mosaico_final.select('quality_score'), {
    min: 0, max: 1,
    palette: ['red', 'yellow', 'green']
}, '📈 Quality Score Final');

// ------- 13. ANÁLISE DE CONTRIBUIÇÃO -------
var contribuicoes = mosaicos_com_quality.map(function(img) {
    return ee.Feature(null, {
        'ano': img.get('year'),
        'semestre': img.get('semestre'),
        'bloco': img.get('bloco'),
        'quality_medio': img.get('quality_mean')
    });
});

print('📅 Contribuição por mosaico:', contribuicoes);

// ------- 14. EXPORTAÇÃO -------
function converterPara16Bit(img) {
    var bandas = img.bandNames();
    
    var bandas_16bit = bandas.map(function(b) {
        b = ee.String(b);
        var banda = img.select(b);
        
        var fator = ee.Number(ee.Algorithms.If(
            b.equals('frequency').or(b.equals('year_semestre')),
            1,
            10000
        ));
        
        return banda.multiply(fator).round().toInt16().rename([b]);
    });
    
    return ee.ImageCollection(bandas_16bit).toBands()
        .rename(bandas)
        .copyProperties(img, img.propertyNames())
        .set('DATA_TYPE', 'INT16')
        .set('SCALE_FACTOR', 10000)
        .set('PROCESSING_DATE', ee.Date(Date.now()).format('YYYY-MM-dd'));
}

var mosaico_16bit = converterPara16Bit(mosaico_final);

// Exportar para Asset
Export.image.toAsset({
    image: mosaico_16bit,
    description: 'MOSAICO_FINAL_ASTER_MULTITEMPORAL',
    assetId: 'MOSAICO_FINAL_ASTER_MULTITEMPORAL',
    region: area_estudo,
    scale: 15,
    maxPixels: 1e13,
    pyramidingPolicy: {
        '.default': 'sample',
        'quality_score': 'max',
        'frequency': 'sum',
        'year_semestre': 'mode'
    }
});

print('✅ Processamento concluído!');
print('🎯 Mosaico final combina os melhores pixels de', mosaicos_sem.size(), 'mosaicos semestrais');