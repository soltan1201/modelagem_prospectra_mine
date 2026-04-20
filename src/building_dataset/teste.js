// ------- 1. PARÂMETROS E BLACKLIST -------
var list_coord = [
    [-44.460483519818084,-14.790793773956535],
    [-39.203525512005584,-14.790793773956535],
    [-39.203525512005584,-9.416525285369607], 
    [-44.460483519818084,-9.416525285369607], 
    [-44.460483519818084,-14.790793773956535]
];
var area_estudo = ee.Geometry.Polygon(list_coord);

var ano = 2007;
var mes = 1; 

// Blacklist
var blacklist = [
    "ASTER/AST_L1T_003/20070105131147",
    "ASTER/AST_L1T_003/20070105131156",
    // ... (sua lista completa)
];

var blacklist_clean = blacklist.map(function(id) {
    return id.split('/').pop();
});

var data_inicio = ee.Date.fromYMD(ano, mes, 1);
var data_fim = data_inicio.advance(6, 'month');

// ------- 2. FUNÇÕES DE PROCESSAMENTO -------

function radianceToTOA(imgRad) {
    imgRad = ee.Image(imgRad);  
    
    var ESUN = {
        'B01': 1848, 'B02': 1549,  'B3N': 1114,
        'B04': 225.4,'B05': 86.63, 'B06': 81.85,
        'B07': 74.85,'B08': 66.49, 'B09': 59.85
    };

    var solarElev = ee.Number(imgRad.get('SOLAR_ELEVATION'));
    var sunZen = ee.Number(90).subtract(solarElev);
    var cosz = sunZen.multiply(Math.PI / 180).cos();

    var doy = ee.Number(ee.Date(imgRad.get('system:time_start')).getRelative('day', 'year'));
    var d = ee.Number(1).subtract(
        ee.Number(0.01672).multiply(
            doy.multiply(2 * Math.PI / 365).cos()
        )
    );
    
    var optical = Object.keys(ESUN).map(function(b) {
        var esun = ee.Number(ESUN[b]);
        return ee.Image(imgRad).select(b)
                    .multiply(Math.PI)
                    .multiply(d.pow(2))
                    .divide(esun.multiply(cosz))
                    .rename(b)
                    .float();
    });

    return ee.Image(optical)
        .addBands(imgRad.select(['B10','B11','B12','B13','B14']).float())
        .copyProperties(imgRad, imgRad.propertyNames());
}

function mascaraNuvemASTER(imgToa) {
    imgToa = ee.Image(imgToa);
    var green = imgToa.select('B01');
    var red   = imgToa.select('B02');
    var nir   = imgToa.select('B3N');
    var swir1 = imgToa.select('B04');
    var tir10 = imgToa.select('B10');

    // Geometria solar
    var solar_el = ee.Number(imgToa.get('SOLAR_ELEVATION'));
    var solar_el_rad = solar_el.multiply(Math.PI / 180);
    var tan_solar_el = solar_el_rad.tan();
    var cloud_height = 2000;
    
    var shadow_dist = ee.Number(cloud_height)
        .divide(tan_solar_el.max(0.1))
        .max(200)
        .min(2000);

    var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
    var visBright = green.add(red).add(nir).divide(3);

    // Nuvens
    var cloud = visBright.gt(0.25)
        .and(swir1.gt(0.15))
        .and(ndvi.lt(0.4));

    var tirStats = tir10.reduceRegion({
        reducer: ee.Reducer.percentile([20]),
        geometry: imgToa.geometry(),
        scale: 90,
        maxPixels: 1e8,
        bestEffort: true
    });
    var p20 = ee.Number(tirStats.get('B10'));
    cloud = cloud.and(tir10.lt(p20));
    cloud = cloud.focal_max({radius: 120, units: 'meters'});

    // Sombras
    var dark_pixels = green.lt(0.14).and(nir.lt(0.175));
    var shadow = dark_pixels.and(
        cloud.focal_max({
            radius: ee.Number(shadow_dist).add(300),
            units: 'meters'
        })
    );
    shadow = shadow.focal_max({radius: 90, units: 'meters'});

    // Máscara
    var valid = imgToa.select(['B01','B02','B3N','B04']).reduce(ee.Reducer.min()).gt(0);
    var clear = valid.and(cloud.not()).and(shadow.not());

    return imgToa.updateMask(clear)
        .copyProperties(imgToa, imgToa.propertyNames());
}

// 👉 FUNÇÃO QUALITY (PARA QUALITYMOSAIC)
function addQualityASTER(img) {
    var b3n = img.select('B3N');
    
    // Converte Float (0-0.4) para Byte (0-255)
    var b3n_byte = b3n.unitScale(0, 0.4).multiply(255).toByte();
    
    // GLCM para detectar bordas de nuvem
    var glcm = b3n_byte.glcmTexture({size: 1});
    var contrast = glcm.select('B3N_contrast');
    
    // Qualidade: Alto NIR + Baixo Contraste
    var quality = b3n.unitScale(0, 0.35)
        .add(contrast.unitScale(0, 255).multiply(-1))
        .rename('quality');
    
    return img.addBands(quality);
}

// 👉 FUNÇÃO PARA CONVERTER PARA 16-BIT (INTEIRO)
function converterPara16Bit(img) {
    // Fatores de escala para cada banda (preserva precisão)
    var escala = {
        'B01': 10000,   // Reflectância 0-1 -> 0-10000
        'B02': 10000,
        'B3N': 10000,
        'B04': 10000,
        'B05': 10000,
        'B06': 10000,
        'B07': 10000,
        'B08': 10000,
        'B09': 10000,
        'B10': 10,      // Temperatura (0-400K -> 0-4000)
        'B11': 10,
        'B12': 10,
        'B13': 10,
        'B14': 10,
        'NDVI': 10000   // -1 a 1 -> -10000 a 10000
    };
    
    var bandas = img.bandNames();
    
    var bandas_16bit = bandas.map(function(b) {
        b = ee.String(b);
        var fator = ee.Number(escala[b] || 10000);  // Default 10000
        
        var banda = img.select(b)
            .multiply(fator)
            .round()
            .toInt16()
            .rename([b]);
        
        return banda;
    });
    
    return ee.ImageCollection(bandas_16bit).toBands()
        .rename(bandas)
        .copyProperties(img, img.propertyNames())
        .set('SCALE_FACTORS', escala)
        .set('DATA_TYPE', 'INT16');
}

// 👉 FUNÇÃO PARA CALCULAR ÍNDICES (ANTES DO MOSAICO)
function adicionarIndices(img) {
    var ndvi = img.normalizedDifference(['B3N', 'B02']).rename('NDVI');
    
    // Al-OH (Argilas) - se SWIR disponível
    var aloh = img.expression(
        "(b('B05') + b('B07')) / (2 * b('B06'))",
        {'B05': img.select('B05'), 'B06': img.select('B06'), 'B07': img.select('B07')}
    ).rename('AlOH_Clay');
    
    // Carbonato
    var carbonato = img.expression(
        "b('B08') / b('B06')",
        {'B06': img.select('B06'), 'B08': img.select('B08')}
    ).rename('Carbonate');
    
    // Ferric Iron
    var fe3 = img.expression(
        "b('B02') / b('B01')",
        {'B01': img.select('B01'), 'B02': img.select('B02')}
    ).rename('Ferric_Fe');
    
    return img.addBands([ndvi, aloh, carbonato, fe3]);
}

// ------- 3. PROCESSAMENTO DA COLEÇÃO -------
var colecao = ee.ImageCollection('ASTER/AST_L1T_003')
    .filterDate(data_inicio, data_fim)
    .filterBounds(area_estudo)
    .filter(ee.Filter.lt('CLOUDCOVER', 70))
    .filter(ee.Filter.inList('system:index', blacklist_clean).not())
    .map(radianceToTOA)
    .map(mascaraNuvemASTER)
    .map(adicionarIndices)
    .map(addQualityASTER);  // 👈 ADICIONA BANDA DE QUALIDADE

print("Total de imagens:", colecao.size());

// ------- 4. QUALITY MOSAIC (MELHOR QUE MEDIAN) -------
var mosaico_quality = colecao.qualityMosaic('quality').clip(area_estudo);

print("✅ Mosaico Quality criado:", mosaico_quality);

// ------- 5. VISUALIZAÇÃO -------
Map.centerObject(area_estudo, 8);

// RGB Natural
Map.addLayer(mosaico_quality, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35,
    gamma: 1.4
}, 'Quality Mosaic - RGB');

// Falsa-cor Geológica
Map.addLayer(mosaico_quality, {
    bands: ['B04', 'B06', 'B3N'],
    min: 0.1, max: 0.5,
    gamma: 1.2
}, 'SWIR Falsa-cor');

// Índice de Argila
Map.addLayer(mosaico_quality.select('AlOH_Clay'), {
    min: 0.8, max: 1.4,
    palette: ['blue', 'cyan', 'green', 'yellow', 'red']
}, 'Al-OH Clay Index');

// NDVI
Map.addLayer(mosaico_quality.select('NDVI'), {
    min: -0.2, max: 0.8,
    palette: ['brown', 'yellow', 'green']
}, 'NDVI');

// ------- 6. CONVERSÃO PARA 16-BIT E EXPORTAÇÃO -------
var mosaico_16bit = converterPara16Bit(mosaico_quality);

print("✅ Mosaico 16-bit:", mosaico_16bit);
print("Bandas disponíveis:", mosaico_16bit.bandNames());

// Lista de bandas para exportar
var bandas_para_asset = [
    'B01', 'B02', 'B3N',           // VNIR 15m
    'B04', 'B05', 'B06', 'B07', 'B08', 'B09',  // SWIR 30m
    'NDVI', 'AlOH_Clay', 'Carbonate', 'Ferric_Fe'  // Índices
];

// 👉 EXPORTAR PARA ASSET (Formato 16-bit)
Export.image.toAsset({
    image: mosaico_16bit.select(bandas_para_asset),
    description: 'ASTER_QualityMosaic_' + ano + '_' + mes + '_INT16',
    assetId: 'ASTER_QualityMosaic_' + ano + '_' + mes + '_INT16',
    region: area_estudo,
    scale: 15,  // Resolução do VNIR (SWIR será upsampled)
    maxPixels: 1e13,
    pyramidingPolicy: {
        '.default': 'sample',  // Melhor para dados contínuos
        'NDVI': 'sample',
        'AlOH_Clay': 'sample'
    }
});

// 👉 EXPORTAR PARA DRIVE (CSV com metadados)
var metadados = ee.FeatureCollection([
    ee.Feature(null, {
        'ano': ano,
        'mes': mes,
        'num_imagens': colecao.size(),
        'data_inicio': data_inicio.format('YYYY-MM-dd'),
        'data_fim': data_fim.format('YYYY-MM-dd'),
        'tipo_mosaico': 'qualityMosaic',
        'formato': 'INT16',
        'fator_escala_reflectancia': 10000,
        'fator_escala_termal': 10
    })
]);

Export.table.toDrive({
    collection: metadados,
    description: 'METADADOS_ASTER_' + ano + '_' + mes,
    folder: 'ASTER_Metadados',
    fileFormat: 'CSV'
});

print("✅ Exportações iniciadas!");
print("📊 Mosaico salvo como INT16 (economiza 50% de espaço)");

// ------- 7. VERIFICAÇÃO DE QUALIDADE -------
// Compara median vs qualityMosaic
var mosaico_median = colecao.median().clip(area_estudo);

Map.addLayer(mosaico_median, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35
}, 'Median Mosaic (Comparação)', false);

print("🔍 Compare as camadas 'Quality Mosaic' vs 'Median Mosaic'");
print("   O Quality Mosaic preserva bordas e evita ghosting!");