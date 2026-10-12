

// 👉 FUNÇÃO QUALIDADE V6 - SOMBRAS COM DETECÇÃO ROBUSTA PARA TOA
var year = '2002';
//https://code.earthengine.google.com/b64f68aca9b02c573151b62708e627d4
function addQualityASTER_degrade(img) {
  var asset_mapbiomas = 'projects/mapbiomas-public/assets/brazil/lulc/collection10/mapbiomas_brazil_collection10_integration_v2';
  var colMapbiomas = ee.Image(asset_mapbiomas).select('classification_' + year)
                        .clip(img.geometry()).eq(33).focal_max(5);

  // === BANDAS ===
  var b1  = img.select('B01');
  var b2  = img.select('B02');
  var b3n = img.select('B3N');
  var b4  = img.select('B04');

  // === ÍNDICES BASE (simplificados) ===
  var ndvi = b3n.subtract(b2).divide(b3n.add(b2).add(1e-6));
  var brightness = b1.add(b2).add(b3n).add(b4).divide(4);
  var ndsi = b2.subtract(b4).divide(b2.add(b4).add(1e-6));

  // === SCORE BASE (normal) ===
  var bright_dev = brightness.subtract(0.25).abs();
  var q_bright = ee.Image(1).subtract(bright_dev.unitScale(0, 0.3));
  var q_ndsi = ee.Image(1).subtract(ndsi.unitScale(0.2, 0.6));
  var q_ndvi = ndvi.unitScale(-0.1, 0.6);

  var base_quality = q_bright.multiply(0.45)
                             .add(q_ndsi.multiply(0.30))
                             .add(q_ndvi.multiply(0.25))
                             .clamp(0, 1);

  // === DEGRADÊ DE QUALIDADE BASEADO EM B3N ===
  // Fórmula: quality_factor = (B3N / 0.16) * 0.6  [para B3N < 0.16]
  //          quality_factor = 1.0                   [para B3N >= 0.16]
  
  var quality_factor = b3n.unitScale(0, 0.16)  // Mapeia 0→0, 0.16→1
                          .multiply(0.6)        // Escala para 0→0, 0.16→0.6
                          .where(b3n.gte(0.16), 1.0);  // Acima de 0.16 = 1.0

  // === APLICAÇÃO DO DEGRADÊ + MAPBIOMAS ===
  // Multiplica qualidade base pelo fator de degradê
  var quality_with_shadow = base_quality.multiply(quality_factor);
  
  // Recupera solo exposto escuro (classe 33) que não é sombra
  quality_with_shadow = quality_with_shadow.where(
    colMapbiomas.and(b3n.lt(0.16)),  // Se é classe 33 E B3N < 0.16
    base_quality                      // Mantém qualidade base (sem penalidade)
  );

  return img.addBands(quality_with_shadow.rename('quality'));
           
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



var asset = "ASTER/AST_L1T_003/20020913130841";


var imgAster = ee.Image(asset);
imgAster = ee.Image.cat(radianceToTOA(imgAster));
imgAster = addQualityASTER_degrade(imgAster);

Map.addLayer(imgAster, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35,
    gamma: 1.4
}, 'Quality Mosaic - RGB');

Map.addLayer(imgAster.select('quality'), {
    min: 0, max: 1,
    palette: ['#8B0000', '#FF4500', '#FF8C00', '#FFD700', '#ADFF2F', '#006400'],
    opacity: 0.7
}, 'Quality (degradê)');

Map.addLayer(imgAster.select('shadow'),{}, 'shadow');

