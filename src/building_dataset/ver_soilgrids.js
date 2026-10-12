// ============================================================================
// SoilGrids 2.0 (ISRIC, 250 m) sobre o bloco de Irecê
// Mostra as variáveis de solo na grid 5x5 junto com os pontos de fosfato e
// imprime estatísticas perto dos pontos vs grid inteira.
// Assets: projects/soilgrids-isric/<variavel>_mean, bandas <var>_<prof>_mean
// ============================================================================

var grid    = ee.FeatureCollection('projects/mapbiomas-arida/mine/grid_irece_5x5_solape100');
var fosfato = ee.FeatureCollection('projects/mapbiomas-arida/mine/pontos_fosfato');
var area    = grid.geometry().bounds();

// Profundidade exibida: 0-5cm, 5-15cm, 15-30cm, 30-60cm, 60-100cm, 100-200cm
var PROF = '0-5cm';

// fator: divide o valor inteiro do SoilGrids para chegar na unidade indicada
var VARS = {
  phh2o:    {fator: 10,  unid: 'pH',           min: 5,   max: 8.5, pal: ['d7191c','fdae61','ffffbf','abd9e9','2c7bb6']},
  sand:     {fator: 10,  unid: 'areia %',      min: 20,  max: 80,  pal: ['f7fcb9','fec44f','d95f0e','993404']},
  clay:     {fator: 10,  unid: 'argila %',     min: 10,  max: 50,  pal: ['fff7ec','fdbb84','e34a33','7f0000']},
  silt:     {fator: 10,  unid: 'silte %',      min: 5,   max: 35,  pal: ['f7f4f9','c994c7','dd1c77','67001f']},
  cfvo:     {fator: 10,  unid: 'fragm. grossos %', min: 0, max: 30, pal: ['ffffff','bdbdbd','636363','000000']},
  cec:      {fator: 10,  unid: 'CTC cmol/kg',  min: 5,   max: 40,  pal: ['f7fbff','6baed6','08519c','08306b']},
  soc:      {fator: 10,  unid: 'C org g/kg',   min: 2,   max: 25,  pal: ['ffffe5','78c679','238443','004529']},
  bdod:     {fator: 100, unid: 'dens. g/cm3',  min: 1.1, max: 1.6, pal: ['fff5eb','fd8d3c','a63603']},
  nitrogen: {fator: 100, unid: 'N g/kg',       min: 0.3, max: 2.5, pal: ['fcfbfd','9e9ac8','54278f']},
  ocd:      {fator: 10,  unid: 'dens. C kg/m3', min: 5,  max: 50,  pal: ['ffffe5','78c679','004529']}
};

function camada(nome) {
  var p = VARS[nome];
  return ee.Image('projects/soilgrids-isric/' + nome + '_mean')
           .select(nome + '_' + PROF + '_mean')
           .divide(p.fator).rename(nome).clip(area);
}

var nomes = Object.keys(VARS);
var solo = ee.Image.cat(nomes.map(camada));

Map.centerObject(grid, 8);
nomes.forEach(function (n, i) {
  var p = VARS[n];
  // só pH, areia e fragmentos ligados de início
  Map.addLayer(solo.select(n), {min: p.min, max: p.max, palette: p.pal},
               n + ' (' + p.unid + ', ' + PROF + ')', ['phh2o', 'sand', 'cfvo'].indexOf(n) >= 0);
});
// estoque de carbono só existe para 0-30 cm
Map.addLayer(ee.Image('projects/soilgrids-isric/ocs_mean').clip(area),
             {min: 10, max: 60, palette: ['ffffe5', '78c679', '004529']}, 'ocs (t/ha, 0-30cm)', false);

Map.addLayer(ee.Image().paint(grid, 0, 1), {palette: '000000'}, 'grid 5x5');
Map.addLayer(fosfato, {color: 'ff00ff'}, 'pontos fosfato');

// --- Estatísticas: buffer 500 m dos pontos vs grid inteira ---
var reducer = ee.Reducer.median().combine(ee.Reducer.percentile([10, 90]), '', true);
var perto = solo.reduceRegion({reducer: reducer, geometry: fosfato.geometry().buffer(500),
                               scale: 250, maxPixels: 1e9});
var tudo  = solo.reduceRegion({reducer: reducer, geometry: area,
                               scale: 250, maxPixels: 1e9});
print('Mediana/p10/p90 a até 500 m dos pontos de fosfato', perto);
print('Mediana/p10/p90 na grid inteira', tudo);

// --- Histogramas comparando pontos vs fundo (pH e areia) ---
['phh2o', 'sand', 'cfvo'].forEach(function (n) {
  var pts = solo.select(n).sampleRegions({collection: fosfato, scale: 250, geometries: false})
                .map(function (f) { return f.set('grupo', 'fosfato'); });
  var fundo = solo.select(n).sample({region: area, scale: 250, numPixels: 2000, seed: 1})
                  .map(function (f) { return f.set('grupo', 'fundo'); });
  print(ui.Chart.feature.histogram(pts, n, 20).setOptions({title: n + ' nos pontos de fosfato'}));
  print(ui.Chart.feature.histogram(fundo, n, 20).setOptions({title: n + ' no fundo (grid)'}));
});
