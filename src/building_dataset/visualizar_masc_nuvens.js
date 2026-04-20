var vis_raw = {
    mosaico_clean:  {
        bands: ['B3N', 'B02', 'B01'], 
        min: 0.05, max: 0.3
    },
    mosaic_raw: {
      bands: ['B3N', 'B02', 'B01'], 
      min: 25, max: 150
    },
    layer_masc: {
        min: 0, max: 2,
        palette: ['black', 'red', 'blue']
    }
}
function mascaraNuvemASTER(img) {
    img = ee.Image(img);
    var b10 = img.select('B10');
    var b01 = img.select('B01');
    var b3n = img.select('B3N');
    var b04 = img.select('B04');  // SWIR1 também ajuda
    
    // ===== DETECÇÃO DE NUVENS =====
    // Usa múltiplos critérios para maior precisão
    var cloud_bright = b01.gt(0.35);           // Muito brilhante no verde
    var cloud_swir = b04.gt(0.35);             // Brilhante no SWIR
    var cloud_cold = b10.lt(5.0);              // Muito fria no térmico
    
    // Combina critérios (nuvem precisa ser brilhante E fria)
    var cloud = cloud_bright.and(cloud_cold).and(cloud_swir);
    
    // Expande para cobrir bordas e nuvens finas
    cloud = cloud.focal_max({radius: 150, units: 'meters'});
    
    // ===== DETECÇÃO DE SOMBRAS =====
    // Sombras são escuras e próximas a nuvens
    var dark_vis = b01.lt(0.08);               // Escuro no visível
    var dark_nir = b3n.lt(0.05);               // Muito escuro no NIR
    
    var shadow = dark_nir.and(dark_vis)
                        .and(cloud.focal_max({radius: 800, units: 'meters'}));
    
    // ===== MÁSCARA FINAL =====
    var mask = cloud.or(shadow).not();
    
    // Remove pequenos artefatos (opcional)
    mask = mask.focal_min({radius: 30, units: 'meters'});
    
    return img.addBands(mask.rename('mask'));
}
function transfor_aster_TOA(img) {
    var ESUN_ASTER = {
        'B01': 1845.99, 'B02': 1555.74, 'B3N': 1119.47,
        'B04': 231.25, 'B05': 79.81, 'B06': 74.99, 
        'B07': 68.66, 'B08': 59.74, 'B09': 56.92
    };
    
    // Geometria solar com validação
    var solarElev = ee.Number(img.get('SOLAR_ELEVATION'));
    solarElev = solarElev.max(5).min(90);  // Evita ângulos inválidos
    
    var sunZenith = ee.Number(90).subtract(solarElev);
    var sunZenithRad = sunZenith.multiply(Math.PI / 180);
    
    // Distância Terra-Sol
    var doy = img.date().getRelative('day', 'year').add(1);
    var dist = ee.Number(1).subtract(
        ee.Number(0.01673).multiply(doy.multiply(2).multiply(Math.PI).divide(365).cos())
    );
    var d2 = dist.pow(2);
    
    // Processa todas as bandas ópticas
    var bandas_opticas = Object.keys(ESUN_ASTER);
    
    var bands = bandas_opticas.map(function(b) {
        var irr = ee.Number(ESUN_ASTER[b]);
        return img.select(b)
            .multiply(Math.PI)
            .multiply(d2)
            .divide(irr.multiply(sunZenithRad.cos()))
            .rename(b)
            .float()
            .clamp(0, 2.0);  // Reflectância TOA raramente > 2.0
    });
    
    var img_toa = ee.Image(bands);
    
    // Adiciona bandas térmicas e metadados
    return img_toa
        .addBands(img.select(['B10', 'B11', 'B12', 'B13', 'B14']))
        .copyProperties(img, img.propertyNames())
        .set('TOA_SCALE', 1.0)
        .set('TOA_OFFSET', 0.0);
}
var asset_img = "ASTER/AST_L1T_003/20210316131118";
var img_aster = ee.Image(asset_img);
print("know metadata of images Aster ", img_aster);
img_aster = ee.Image(transfor_aster_TOA(img_aster));
print("know metadata of images Aster TOA", img_aster);
// Uso:
var img_clean = mascaraNuvemASTER(img_aster);


Map.addLayer(img_aster, vis_raw.mosaico_clean, 'raster original TOA');

Map.addLayer(img_clean.select('mask'), vis_raw.layer_masc, 'Nuvem/Sombra (Vermelho=Nuvem, Azul=Sombra)');