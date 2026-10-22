// ============================================================================
// Mapa de localização — biomas, estados e área de estudo (bloco de Irecê, BA)
//
// Fontes (projects/mapbiomas-arida/)
//   biomas_2026   6 biomas IBGE (propriedade 'Bioma')
//   BR_UF_2025    27 UFs IBGE 2025 ('SIGLA_UF', 'NM_UF')
//   mine/grid_irece_5x5_solape100   grid 5×5 da área de estudo
// ============================================================================

var BIOMAS  = ee.FeatureCollection('projects/mapbiomas-arida/biomas_2026');
var ESTADOS = ee.FeatureCollection('projects/mapbiomas-arida/BR_UF_2025');
var GRID    = ee.FeatureCollection('projects/mapbiomas-arida/mine/grid_irece_5x5_solape100');
var AREA    = ee.Geometry.Rectangle([-44.460483519818084, -14.790793773956535,
                                     -39.203525512005584, -9.416525285369607]);

// cores dos biomas (padrão próximo ao do IBGE/MapBiomas)
var COR_BIOMA = {
  'Amazônia': '#1f8d49', 'Caatinga': '#F8FC68', 'Cerrado': '#d4a95a',
  'Mata Atlântica': '#5ba745', 'Pampa': '#b8d98f', 'Pantanal': '#8fd1e0'
};

var biomasPintados = BIOMAS.map(function (f) {
  return f.set('estilo', {fillColor: ee.String(ee.Dictionary(COR_BIOMA).get(f.get('Bioma'), 'cccccc')).cat('cc'),
                          color: '00000000', width: 0});
}).style({styleProperty: 'estilo'});
Map.addLayer(ee.Image.constant(1), {min:0, max: 1}, 'base');
Map.addLayer(biomasPintados, {}, 'Biomas (IBGE)');
Map.addLayer(ESTADOS.style({color: '444444', width: 1, fillColor: '00000000'}), {}, 'Estados (IBGE 2025)');
Map.addLayer(ESTADOS.filter(ee.Filter.eq('SIGLA_UF', 'BA')).style({color: '000000', width: 2.5, fillColor: '00000000'}),
             {}, 'Bahia');
Map.addLayer(ee.Image().paint(GRID, 0, 1), {palette: '555555'}, 'Grid 5×5', false);
Map.addLayer(ee.Image().paint(AREA, 0, 3), {palette: 'e31a1c'}, 'Área de estudo');

// siglas das UFs no centro de cada estado (camada de pontos com rótulo no Inspector)
Map.addLayer(ESTADOS.map(function (f) { return ee.Feature(f.geometry().centroid(1000), {uf: f.get('SIGLA_UF')}); })
               .style({color: '000000', pointSize: 2}), {}, 'Centro das UFs', false);

// ----------------------------------------------------------------------------
// Legenda
// ----------------------------------------------------------------------------
var legenda = ui.Panel({style: {position: 'bottom-left', padding: '8px 12px'}});
legenda.add(ui.Label('Localização', {fontWeight: 'bold', fontSize: '14px', margin: '0 0 6px 0'}));
function item(cor, texto, borda) {
  var caixa = ui.Label('', {backgroundColor: borda ? 'ffffff' : '#' + cor, padding: '8px', margin: '0 6px 4px 0',
                            border: borda ? ('2px solid #' + cor) : '1px solid #999999'});
  return ui.Panel([caixa, ui.Label(texto, {fontSize: '12px', margin: '0 0 4px 0'})], ui.Panel.Layout.flow('horizontal'));
}
Object.keys(COR_BIOMA).forEach(function (b) { legenda.add(item(COR_BIOMA[b], b)); });
legenda.add(item('000000', 'Bahia', true));
legenda.add(item('e31a1c', 'Área de estudo (Irecê)', true));
Map.add(legenda);

// área de estudo: tamanho e biomas que ela cruza
var info = ui.Label('calculando…', {fontSize: '11px', whiteSpace: 'pre-wrap'});
legenda.add(info);
ee.Dictionary({
  area_km2: AREA.area(100).divide(1e6),
  biomas: BIOMAS.filterBounds(AREA).aggregate_array('Bioma')
}).evaluate(function (r) {
  info.setValue('Área de estudo: ' + Math.round(r.area_km2).toLocaleString('pt-BR') + ' km²\nBiomas: ' +
                r.biomas.join(', '));
});

Map.setOptions('TERRAIN');
Map.centerObject(AREA, 5);
