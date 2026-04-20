var vis_raw = {
    mosaico_clean:  {
        bands: ['B3N', 'B02', 'B01'], 
        min: 0.05, max: 0.3
    },
    mosaic_raw_vnir: {
        bands: ['B3N', 'B02', 'B01'],
        min: 30, max: 215,
        gamma: 1.15
    },
    mosaic_raw_swir: {
        bands: ['B04', 'B06', 'B09'],
        min: 25, max: 170,
        gamma: 1.2
    },
    layer_masc: {
        min: 0, max: 2,
        palette: ['black', 'red', 'blue']
    }
}

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
    print("numero PI ", Math.PI);
    print("date image ", ee.Date(imgRad.get('system:time_start')));
    print("number doy = ", doy);
    print("number d = ", d);
    print("lista de bandas ", Object.keys(ESUN));
    var optical = Object.keys(ESUN).map(
        function(b) {
            var esun = ee.Number(ESUN[b]);
            return ee.Image(imgRad).select(b)
                        .multiply(Math.PI)
                        .multiply(d.pow(2))
                        .divide(esun.multiply(cosz))
                        .rename(b);
        }
    );

    return ee.Image(optical)
                .addBands(imgRad.select(['B10','B11','B12','B13','B14']))
                .copyProperties(imgRad, imgRad.propertyNames());
}

// -----------------------------
// 3) Máscara de nuvem/sombra
// -----------------------------
function maskCloudShadowASTER(imgToa) {
  
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
        maxPixels: 1e8,
        bestEffort: true
    });
    var p20 = ee.Number(tirStats.get('B10'));
    cloud = cloud.and(tir10.lt(p20));
    cloud = cloud.focal_max({radius: 120, units: 'meters'}).rename('cloud');

    // ===== DETECÇÃO DE SOMBRAS COM CORREÇÃO GEOMÉTRICA =====
    print("show shadow_dist = ", shadow_dist);
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

    // var shadow = dark.and(cloud.focal_max({radius: 1500, units: 'meters'}));    
    // shadow = shadow.focal_max({radius: 90, units: 'meters'}).rename('shadow');
    // Map.addLayer(shadow.selfMask(), {min: 0, max: 1, palette: 'green'}, 'shadow');

    // máscara de pixels válidos
    var valid = imgToa.select(['B01','B02','B3N','B04']).reduce(ee.Reducer.min()).gt(0);
    var clear = valid.and(cloud.not()).and(shadow.not()).rename('clear_mask');

    return imgToa.addBands([ndvi])
                .updateMask(clear);
}

var asset_img = "ASTER/AST_L1T_003/20070105131232";
var img_rand = ee.Image(asset_img);
print("show metadata ", img_rand);
var imgToa = ee.Image(radianceToTOA(img_rand));
var imgClean = maskCloudShadowASTER(imgToa);
print("show all metadata and bands from image clean ", imgClean);

// Map.centerObject(img, 9);
// Map.addLayer(img, {bands:['B3N','B02','B01'], min: 33, max: 215}, 'ASTER raw');
Map.addLayer(img_rand, vis_raw.mosaic_raw_vnir, 'ASTER 2007 - VNIR', false);
Map.addLayer(img_rand, vis_raw.mosaic_raw_swir, 'ASTER 2007 - SWIR', false);
Map.addLayer(imgToa, {bands:['B3N','B02','B01'], min:0.03, max:0.35}, 'ASTER TOA');
Map.addLayer(imgClean, {bands:['B3N','B02','B01'], min:0.03, max:0.35}, 'ASTER sem nuvem', false);
Map.addLayer(imgClean.select('cloud').selfMask(), {palette:['red']}, 'nuvem'); //
Map.addLayer(imgClean.select('shadow').selfMask(), {palette:['blue']}, 'sombra');