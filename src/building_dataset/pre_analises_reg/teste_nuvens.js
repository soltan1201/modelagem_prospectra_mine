// =================================================================
// PAINEL DE CONTROLE DAS MÁSCARAS (
// =================================================================

//Proposta 02

// 1. Limiares de Nuvem
var lim_visBright = 0.30;      // => Vai para: visBright.gt(lim_visBright)
var lim_swir1 = 0.15;          // => Vai para: swir1.gt(lim_swir1) 
var ajuste_termico = 40;       // => Vai para: CLOUDCOVER / 2 + ajuste_termico
var raio_borda_nuvem = 150;    // => Vai para: cloud.focal_max({radius: raio_borda_nuvem})

// 2. Limiares de Sombra 
var altura_nuvem_estimada = 2000; // => Vai para: shadow_dist = ee.Number(altura_nuvem_estimada)
var lim_escuro_green = 0.18;   // => Vai para: green.lt(lim_escuro_green) 
var lim_escuro_nir = 0.20;    // => Vai para: nir.lt(lim_escuro_nir) 
var raio_borda_sombra = 90;    // => Vai para: shadow.focal_max({radius: raio_borda_sombra})


// =================================================================
// CÓDIGO DO LABORATÓRIO
// =================================================================

var vis_raw = {
    mosaico_clean:  { bands: ['B3N', 'B02', 'B01'], min: 0.05, max: 0.3 },
    mosaic_raw_vnir: { bands: ['B3N', 'B02', 'B01'], min: 30, max: 215, gamma: 1.15 },
    mosaic_raw_swir: { bands: ['B04', 'B06', 'B09'], min: 25, max: 170, gamma: 1.2 },
    layer_masc: { min: 0, max: 2, palette: ['black', 'red', 'blue'] }
}

// -----------------------------
// 2) Radiância -> reflectância TOA
// -----------------------------
function radianceToTOA(imgRad) {
    imgRad = ee.Image(imgRad);  
    var ESUN = { 'B01': 1848, 'B02': 1549,  'B3N': 1114, 'B04': 225.4,'B05': 86.63, 'B06': 81.85, 'B07': 74.85,'B08': 66.49, 'B09': 59.85 };
    var solarElev = ee.Number(imgRad.get('SOLAR_ELEVATION'));
    var sunZen = ee.Number(90).subtract(solarElev);
    var cosz = sunZen.multiply(Math.PI / 180).cos();
    var doy = ee.Number(ee.Date(imgRad.get('system:time_start')).getRelative('day', 'year'));
    var d = ee.Number(1).subtract(ee.Number(0.01672).multiply(doy.multiply(2 * Math.PI / 365).cos()));
    
    var optical = Object.keys(ESUN).map(function(b) {
        var esun = ee.Number(ESUN[b]);
        return ee.Image(imgRad).select(b).multiply(Math.PI).multiply(d.pow(2)).divide(esun.multiply(cosz)).rename(b);
    });
    return ee.Image(optical).addBands(imgRad.select(['B10','B11','B12','B13','B14'])).copyProperties(imgRad, imgRad.propertyNames());
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

    var solar_az = ee.Number(imgToa.get('SOLAR_AZIMUTH'));
    var solar_el = ee.Number(imgToa.get('SOLAR_ELEVATION'));
    var solar_el_rad = solar_el.multiply(Math.PI / 180);
    var tan_solar_el = solar_el_rad.tan();
    
    // APLICANDO VARIÁVEL 
    var shadow_dist = ee.Number(altura_nuvem_estimada).divide(tan_solar_el.max(0.1)).max(200).min(2000);

    var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
    var visBright = green.add(red).add(nir).divide(3);

    // APLICANDO VARIÁVEIS DO PAINEL (NDVI comentado)
    var cloud = visBright.gt(lim_visBright);
        //.and(swir1.gt(lim_swir1));
        //.and(ndvi.lt(0.4));
    
    var percentCloud = ee.Number(imgToa.get('CLOUDCOVER')).divide(2).add(ajuste_termico);
    
    print("Percentil térmico calculado = ", percentCloud);
    
    var tirStats = tir10.reduceRegion({
        reducer: ee.Reducer.percentile([percentCloud]),
        geometry: imgToa.geometry(),
        scale: 90,
        maxPixels: 1e8,
        bestEffort: true
    });
    var p20 = ee.Number(tirStats.get('B10'));
    cloud = cloud.and(tir10.lt(p20));
    
    // APLICANDO VARIÁVEL 
    cloud = cloud.focal_max({radius: raio_borda_nuvem, units: 'meters'}).rename('cloud');

    print("show shadow_dist = ", shadow_dist);
    
    // APLICANDO VARIÁVEIS 
    var dark_pixels = green.lt(lim_escuro_green).and(nir.lt(lim_escuro_nir));
    var shadow = dark_pixels.and(
        cloud.focal_max({
            radius: shadow_dist.add(300),
            units: 'meters'
        })
    );

    // APLICANDO VARIÁVEL 
    shadow = shadow.focal_max({radius: raio_borda_sombra, units: 'meters'}).rename('shadow');

    var valid = imgToa.select(['B01','B02','B3N','B04']).reduce(ee.Reducer.min()).gt(0);
    var clear = valid.and(cloud.not()).and(shadow.not()).rename('clear_mask');

    return imgToa.addBands([ndvi, cloud.rename('cloud'), shadow.rename('shadow'), visBright.rename('btb'), clear.rename('mask_clear')]);
}

var asset_img = "ASTER/AST_L1T_003/20001124133229"; // IMAGEM DE TESTE

var img_rand = ee.Image(asset_img);
var imgToa = ee.Image(radianceToTOA(img_rand));
var imgClean = maskCloudShadowASTER(imgToa);

print("CLOUDCOVER oficial do Satélite = ", img_rand.get('CLOUDCOVER'));
print("Elevação Solar = ", img_rand.get('SOLAR_ELEVATION'));
print("Ver todos os metadados originais = ", img_rand);

Map.addLayer(img_rand, vis_raw.mosaic_raw_vnir, 'ASTER 2007 - VNIR', false);
Map.addLayer(img_rand, vis_raw.mosaic_raw_swir, 'ASTER 2007 - SWIR', false);
Map.addLayer(imgToa, {bands:['B3N','B02','B01'], min:0.03, max:0.35}, 'ASTER TOA', false);

Map.addLayer(imgClean.updateMask(imgClean.select('mask_clear')), {bands:['B3N','B02','B01'], min:0.03, max:0.35}, 'ASTER sem nuvem (Fundo)');
Map.addLayer(imgClean.select('cloud').selfMask(), {palette:['red']}, 'nuvem (Máscara)', false); 
Map.addLayer(imgClean.select('shadow').selfMask(), {palette:['blue']}, 'sombra (Máscara)', false);