// ------- 1. PARÂMETROS E BLACKLIST -------
var list_coord = [
    [-44.460483519818084,-14.790793773956535],
    [-39.203525512005584,-14.790793773956535],
    [-39.203525512005584,-9.416525285369607], 
    [-44.460483519818084,-9.416525285369607], 
    [-44.460483519818084,-14.790793773956535]
];
var area_estudo = ee.Geometry.Polygon(list_coord);
Map.addLayer(area_estudo, {}, 'rectangulo');

var vis_raw = {
    mosaico_clean:  {
        bands: ['B3N', 'B02', 'B01'], 
        min: 0.05, max: 0.3
    },
    mosaic_raw: {
      bands: ['B3N', 'B02', 'B01'], 
      min: 25, max: 150
    }
}

var ano = 2007;
var mes = 6; 

// 🛑 COLE AQUI OS IDs DAS IMAGENS QUE VOCÊ QUER REMOVER
// Você encontra o ID clicando na imagem no mapa ou no Console
var blacklist = [
    // 'ASTER/AST_L1T_003/20210715132455', // Exemplo de ID (substitua pelos reais)
    "ASTER/AST_L1T_003/20070623130639",
    "ASTER/AST_L1T_003/20070623130648",
    "ASTER/AST_L1T_003/20070623130657",
    "ASTER/AST_L1T_003/20070911130548",
    "ASTER/AST_L1T_003/20070911130557",
    "ASTER/AST_L1T_003/20070911130708",
    "ASTER/AST_L1T_003/20070911130717",
    "ASTER/AST_L1T_003/20070911130725",
    "ASTER/AST_L1T_003/20070918131334",
];
// Divide pela barra e pega o último elemento
var blacklist_clean = blacklist.map(function(id) {
    return id.split('/').pop();
});

print('Blacklist limpa:', blacklist_clean);

var data_inicio = ee.Date.fromYMD(ano, mes, 1);
var data_fim = data_inicio.advance(6, 'month');

// ------- 2. FUNÇÕES DE PROCESSAMENTO -------

// -----------------------------
// 2) Radiância -> reflectância TOA
// apenas VNIR/SWIR
// -----------------------------
function radianceToTOA(imgRad) {
    // ASTER L1T já vem em RADIÂNCIA (W/m²/sr/µm)
    // Não precisa converter de DN! Os valores já são radiância.
    imgRad = ee.Image(imgRad);  
    
    var ESUN = {
        'B01': 1848, 'B02': 1549,  'B3N': 1114,
        'B04': 225.4,'B05': 86.63, 'B06': 81.85,
        'B07': 74.85,'B08': 66.49, 'B09': 59.85
    };

    var solarElev = ee.Number(imgRad.get('SOLAR_ELEVATION'));
    var sunZen = ee.Number(90).subtract(solarElev);
    var cosz = sunZen.multiply(Math.PI / 180).cos();

    var doy = ee.Number(ee.Date(imgRad.get('system:time_start')).getRelative('day', 'year'))//.add(1);
    var d = ee.Number(1).subtract(
        ee.Number(0.01672).multiply(
        doy.multiply(2 * Math.PI / 365).cos()
        )
    );

    var optical = Object.keys(ESUN).map(
        function(b) {
            var esun = ee.Number(ESUN[b]);
            return ee.Image(imgRad).select(b)
                        .multiply(Math.PI)
                        .multiply(d.pow(2))
                        .divide(esun.multiply(cosz))
                        .rename(b)
                        .float();
        }
    );

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

    // ===== GEOMETRIA SOLAR (CORREÇÃO DINÂMICA) =====
    var solar_az = ee.Number(imgToa.get('SOLAR_AZIMUTH'));
    var solar_el = ee.Number(imgToa.get('SOLAR_ELEVATION'));
    
    // Converte elevação para radianos e calcula tangente
    var solar_el_rad = solar_el.multiply(Math.PI / 180);
    var tan_solar_el = solar_el_rad.tan();
    
    // Estima distância da sombra baseado na elevação solar
    // Fórmula: distância = altura_nuvem / tan(elevação)
    // Altura típica de nuvens: 2000m (pode variar de 1000-4000m)
    var cloud_height = 2000; // metros
    
    var shadow_dist = ee.Number(cloud_height)
        .divide(tan_solar_el.max(0.1))  // Evita divisão por zero
        .max(200)   // Distância mínima: 200m
        .min(2000); // Distância máxima: 2000m

    var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
    // brilho visível
    var visBright = green.add(red).add(nir).divide(3);

    // nuvens: brilhantes, frias, pouco vegetadas
    var cloud = visBright.gt(0.25)
        .and(swir1.gt(0.15))
        .and(ndvi.lt(0.4));

    // frio relativo por percentil da própria cena
    var tirStats = tir10.reduceRegion({
        reducer: ee.Reducer.percentile([20]),
        geometry: imgToa.geometry(),
        scale: 90,
        maxPixels: 1e9,
        bestEffort: true
    });
    var p20 = ee.Number(tirStats.get('B10'));
    cloud = cloud.and(tir10.lt(p20));
    cloud = cloud.focal_max({radius: 120, units: 'meters'}).rename('cloud');

    // ===== DETECÇÃO DE SOMBRAS COM CORREÇÃO GEOMÉTRICA =====
    // print("show shadow_dist = ", shadow_dist);
    // sombras: escuras e próximas das nuvens
    var dark_pixels = green.lt(0.14).and(nir.lt(0.175));
    // Sombras estão na direção oposta ao sol
    // Usa a distância calculada dinamicamente baseada na elevação solar
    var shadow = dark_pixels.and(
        cloud.focal_max({
            radius: ee.Number(shadow_dist).add(300),
            units: 'meters'
        })
    );

    // Dilatação da sombra (suaviza bordas)
    shadow = shadow.focal_max({radius: 90, units: 'meters'}).rename('shadow');

    // máscara de pixels válidos
    var valid = imgToa.select(['B01','B02','B3N','B04']).reduce(ee.Reducer.min()).gt(0);
    var clear = valid.and(cloud.not()).and(shadow.not()).rename('clear_mask');

    return imgToa.addBands([ndvi])
                .updateMask(clear)
                copyProperties(imgToa, imgToa.propertyNames());
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


// Função para forçar tipo de forma mais robusta
function forcarTipoUniforme(img) {
    var bandas = ['B01','B02','B3N','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14'];
    
    var bandas_uniformes = bandas.map(function(b) {
        // Se a banda existe, converte para float com range fixo
        var banda = img.select(b);
        return banda.float().clamp(0, 1000);  // Range amplo para evitar overflow
    });
    
    return ee.Image(bandas_uniformes)
        .rename(bandas)
        .copyProperties(img, img.propertyNames());
}

// ------- 3. FILTRAGEM E APLICAÇÃO DA BLACKLIST -------
var total_nFilter = 214;
print(data_inicio, data_fim);
// ------- 3. PROCESSAMENTO DA COLEÇÃO -------
var colecao = ee.ImageCollection('ASTER/AST_L1T_003')
                    .filterDate(data_inicio, data_fim)
                    .filterBounds(area_estudo)
                    .filter(ee.Filter.lt('CLOUDCOVER', 70))
                    .filter(ee.Filter.inList('system:index', blacklist_clean).not())
                    .map(radianceToTOA)
                    .map(mascaraNuvemASTER)
                    .map(addQualityASTER);  // 👈 ADICIONA BANDA DE QUALIDADE

print("Total de imagens:", colecao.size());

// ------- 4. QUALITY MOSAIC (MELHOR QUE MEDIAN) -------
var mosaico_quality = colecao.qualityMosaic('quality').clip(area_estudo);

print("✅ Mosaico Quality criado:", mosaico_quality);
print(" total sem filtros ", total_nFilter);
print(colecao.limit(10));
// ------- 4. VISUALIZAÇÃO -------
// Map.centerObject(area_estudo, 8);

// Mosaico Resultante (Sem as imagens da Blacklist)
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

// ------- 6. CONVERSÃO PARA 16-BIT E EXPORTAÇÃO -------
var mosaico_16bit = ee.Image(converterPara16Bit(mosaico_quality));
print("✅ Mosaico 16-bit:", mosaico_16bit);
print("Bandas disponíveis:", mosaico_16bit.bandNames());
var semestre = '2';
mosaico_16bit = mosaico_16bit.set(
                  'year', ano,
                  'semestre', 2,
                  'bloco', 'Irece'
                );

// Lista de bandas para exportar
var bandas_para_asset = [
    'B01', 'B02', 'B3N',           // VNIR 15m
    'B04', 'B05', 'B06', 'B07', 'B08', 'B09',  // SWIR 30m   
];
var id_output = 'projects/ee-solkancengine17/assets/prospectra/mosaic_aster/';

var name_to_export = 'ASTER_QualityMosaic_' + ano + '_semestre_' + semestre + '_INT16';
// 👉 EXPORTAR PARA ASSET (Formato 16-bit)
Export.image.toAsset({
    image: mosaico_16bit.select(bandas_para_asset),
    description: name_to_export,
    assetId: id_output  + name_to_export,
    region: area_estudo,
    scale: 15,  // Resolução do VNIR (SWIR será upsampled)
    maxPixels: 1e13,
    pyramidingPolicy: {
        '.default': 'sample'  // Melhor para dados contínuos
    }
});

colecao.evaluate(function(coll) {
    coll.features.forEach(function(feat) {
        var id_simples = feat.id.split('/').pop(); // Pega apenas o final do ID
        // var data = new Date(feat.properties['system:time_start']).toISOString().split('T')[0];
        var img = ee.Image(feat.id);
        
        // Mostra o ID no console para facilitar a cópia para a blacklist
        print('"' + feat.id + '",'); 
        
        Map.addLayer(img, vis_raw.mosaic_raw, 'Inspecionar: ' + ' _ ' + id_simples , false);
    });
});