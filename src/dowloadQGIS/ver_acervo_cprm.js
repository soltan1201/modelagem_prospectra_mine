// ============================================================================
// Acervo CPRM (SGB) — visualização dos assets em
//   projects/mapbiomas-arida/mine/acervo_cprm
// enviados por exportar_acervo_gee.py (ver ACERVO_CPRM_IRECE.md).
//
// Camadas
//   Aerogeofísica 1042 Norte da Chapada Diamantina (~490 m; MAG ~220 m):
//     eU, eTh, K, CT, razões U/Th, U/K, Th/K, ternário K-eTh-eU, MAG
//   Geologia: litologia Jacobina 1:250k e Bahia 1:1M, unidades do Grupo Una /
//     Fm. Salitre em destaque, estruturas, dolinas (Morro do Chapéu)
//   Referências: ocorrências de fosfato (por grupo), recursos minerais,
//     geoquímica de rocha (P2O5)
// Clique no mapa: valores do gama, unidade litológica e ocorrências próximas.
// ============================================================================

var PASTA = 'projects/mapbiomas-arida/mine/acervo_cprm/';

var gama     = ee.Image(PASTA + '1042_gama');   // eU_ppm eTh_ppm K_pct CT U_Th U_K Th_K
var mag      = ee.Image(PASTA + '1042_mag');    // MAG_nT
var litoJac  = ee.FeatureCollection(PASTA + 'lito_jacobina_250k');
var litoBa   = ee.FeatureCollection(PASTA + 'lito_bahia');
var estJac   = ee.FeatureCollection(PASTA + 'estrut_jacobina_250k');
var estBa    = ee.FeatureCollection(PASTA + 'estrut_bahia');
var dolinas  = ee.FeatureCollection(PASTA + 'dolinas_morro_chapeu');
var fosfato  = ee.FeatureCollection(PASTA + 'ocorrencias_fosfato');
var recmin   = ee.FeatureCollection(PASTA + 'recursos_minerais');
var geoq     = ee.FeatureCollection(PASTA + 'geoquimica_p2o5');
var grid     = ee.FeatureCollection('projects/mapbiomas-arida/mine/grid_irece_5x5_solape100');

// Área de estudo (mesmo retângulo de select_black_list_save_mosaicSem.js)
var area = ee.Geometry.Rectangle([-44.460483519818084, -14.790793773956535,
                                  -39.203525512005584, -9.416525285369607]);

// ----------------------------------------------------------------------------
// Parâmetros de visualização (faixas = percentis 2–98 % dos grids)
// ----------------------------------------------------------------------------
var PAL_GAMA = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];   // magma
var PAL_RAZAO = ['2c7bb6', 'abd9e9', 'ffffbf', 'fdae61', 'd7191c'];
var PAL_MAG = ['053061', '2166ac', '92c5de', 'f7f7f7', 'f4a582', 'b2182b', '67001f'];

var VIS = {
  eU_ppm:  {min: 0.2,  max: 3.1,  palette: PAL_GAMA,  nome: 'eU (ppm)'},
  eTh_ppm: {min: 0.4,  max: 19.4, palette: PAL_GAMA,  nome: 'eTh (ppm)'},
  K_pct:   {min: 0.07, max: 1.4,  palette: PAL_GAMA,  nome: 'K (%)'},
  CT:      {min: 1.2,  max: 7.1,  palette: PAL_GAMA,  nome: 'Contagem total'},
  U_Th:    {min: 0.06, max: 0.85, palette: PAL_RAZAO, nome: 'U/Th'},
  U_K:     {min: 0.7,  max: 8.7,  palette: PAL_RAZAO, nome: 'U/K'},
  Th_K:    {min: 1.6,  max: 62,   palette: PAL_RAZAO, nome: 'Th/K'}
};
var VIS_MAG = {min: -190, max: 103, palette: PAL_MAG};

// valores negativos pequenos (correção gama) ficam transparentes nas razões
var gamaPos = gama.updateMask(gama.select('eU_ppm').gte(0));

// ----------------------------------------------------------------------------
// Geologia: cor por unidade; Una / Salitre em destaque
// ----------------------------------------------------------------------------
var PAL_UNID = ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462', 'b3de69',
                'fccde5', 'd9d9d9', 'bc80bd', 'ccebc5', 'ffed6f', 'a6cee3', 'b2df8a',
                'fb9a99', 'cab2d6', 'e5d8bd', 'f1e2cc', 'c7e9c0', 'fdd0a2'];

function corPorUnidade(fc) {
  var unidades = fc.aggregate_array('NOME_UNIDA').distinct();
  var pal = ee.List(PAL_UNID);
  return fc.map(function (f) {
    var i = unidades.indexOf(f.get('NOME_UNIDA'));
    var cor = ee.String(pal.get(ee.Number(i).mod(pal.length())));
    return f.set('estilo', {color: '555555', width: 0.4, fillColor: cor.cat('b0')});
  }).style({styleProperty: 'estilo'});
}

// unidades do Grupo Una (Fm. Salitre e membros, Fm. Bebedouro) — hospedeiros do fosfato
var ALVO = ['Salitre', 'Irecê', 'Jussara', 'Nova América', 'Gabriel', 'Lapão', 'Bebedouro', 'Una'];
var filtroAlvo = ee.Filter.or.apply(null, ALVO.map(function (t) {
  return ee.Filter.stringContains('NOME_UNIDA', t);
}));
var alvoJac = litoJac.filter(filtroAlvo);
var alvoBa  = litoBa.filter(filtroAlvo);

// ----------------------------------------------------------------------------
// Pontos de referência
// ----------------------------------------------------------------------------
var CORES_GRUPO = {bacia_irece: 'ff0000', leste_outros: 'ff9900', norte: 'ffff00'};
var fosfatoEst = fosfato.map(function (f) {
  var cor = ee.Dictionary(CORES_GRUPO).get(f.get('grupo'), 'ffffff');
  return f.set('estilo', {color: '000000', fillColor: cor, pointSize: 7,
                          pointShape: 'star5', width: 1});
}).style({styleProperty: 'estilo'});

// P2O5 (%) da geoquímica de rocha: <1 | 1–5 | >=5
var geoqEst = geoq.map(function (f) {
  var p = ee.Number(f.get('p2o5_pct'));
  var cor = ee.Algorithms.If(p.gte(5), 'e31a1c', ee.Algorithms.If(p.gte(1), 'fd8d3c', '6baed6'));
  var tam = ee.Algorithms.If(p.gte(5), 7, ee.Algorithms.If(p.gte(1), 5, 3));
  return f.set('estilo', {color: '000000', fillColor: cor, pointSize: tam,
                          pointShape: 'circle', width: 0.5});
}).style({styleProperty: 'estilo'});

// ----------------------------------------------------------------------------
// Mapa (camadas desligadas por padrão: marcar no painel "Layers")
// ----------------------------------------------------------------------------
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);

// aerogeofísica
Map.addLayer(mag, VIS_MAG, '1042 · MAG (nT)', false);
Map.addLayer(gamaPos.select(['K_pct', 'eTh_ppm', 'eU_ppm']),
             {min: [0.07, 0.4, 0.2], max: [1.4, 19.4, 3.1]}, '1042 · Ternário K-eTh-eU', false);
['Th_K', 'CT', 'K_pct', 'eTh_ppm', 'U_K', 'U_Th'].forEach(function (b) {
  Map.addLayer(gamaPos.select(b), VIS[b], '1042 · ' + VIS[b].nome, false);
});
Map.addLayer(gamaPos.select('eU_ppm'), VIS.eU_ppm, '1042 · eU (ppm)', true, 0.85);

// geologia
Map.addLayer(corPorUnidade(litoBa), {}, 'Litologia Bahia 1:1M', false, 0.7);
Map.addLayer(corPorUnidade(litoJac), {}, 'Litologia Jacobina 1:250k', false, 0.8);
Map.addLayer(alvoBa.style({color: '00ffff', width: 1.5, fillColor: '00ffff30'}),
             {}, 'Una/Salitre · Bahia 1:1M', false);
Map.addLayer(alvoJac.style({color: '00ffff', width: 2, fillColor: '00ffff40'}),
             {}, 'Una/Salitre · Jacobina 1:250k', true);
Map.addLayer(estBa.style({color: 'ffffff', width: 0.8}), {}, 'Estruturas Bahia 1:1M', false);
Map.addLayer(estJac.style({color: 'ff00ff', width: 0.8}), {}, 'Estruturas Jacobina 1:250k', false);
Map.addLayer(dolinas.style({color: '3300ff', width: 1, fillColor: '3300ff80'}),
             {}, 'Dolinas (Morro do Chapéu)', false);

// referências e limites
Map.addLayer(ee.Image().paint(grid, 0, 1), {palette: 'ffffff'}, 'Grid 5x5', false);
Map.addLayer(ee.Image().paint(area, 0, 2), {palette: 'ffff00'}, 'Área de estudo', true);
Map.addLayer(recmin.style({color: '999999', pointSize: 2, width: 0.5}),
             {}, 'Recursos minerais (todos)', false);
Map.addLayer(geoqEst, {}, 'Geoquímica P2O5 (rocha)', false);
Map.addLayer(fosfatoEst, {}, 'Ocorrências de fosfato', true);

// ----------------------------------------------------------------------------
// Legenda
// ----------------------------------------------------------------------------
function barra(vis, titulo) {
  var img = ee.Image.pixelLonLat().select(0);
  var thumb = ui.Thumbnail({
    image: img, params: {bbox: [0, 0, 1, 0.1], dimensions: '160x10', format: 'png',
                         min: 0, max: 1, palette: vis.palette},
    style: {margin: '0 6px'}
  });
  return ui.Panel([
    ui.Label(titulo, {fontWeight: 'bold', fontSize: '12px', margin: '4px 6px 0 6px'}),
    thumb,
    ui.Panel([ui.Label(String(vis.min), {fontSize: '11px', margin: '0 6px'}),
              ui.Label(String(vis.max), {fontSize: '11px', margin: '0 0 0 100px'})],
             ui.Panel.Layout.flow('horizontal'))
  ]);
}
function item(cor, texto, forma) {
  return ui.Panel([
    ui.Label(forma || '■', {color: '#' + cor, fontSize: '16px', margin: '0 4px 0 6px'}),
    ui.Label(texto, {fontSize: '12px', margin: '3px 0 0 0'})
  ], ui.Panel.Layout.flow('horizontal'));
}

var legenda = ui.Panel({style: {position: 'bottom-left', padding: '6px', width: '230px'}});
legenda.add(ui.Label('Acervo CPRM · Irecê', {fontWeight: 'bold', fontSize: '14px'}));
legenda.add(barra(VIS.eU_ppm, 'eU (ppm) — 1042'));
legenda.add(ui.Label('Ocorrências de fosfato', {fontWeight: 'bold', fontSize: '12px', margin: '6px 6px 0 6px'}));
legenda.add(item(CORES_GRUPO.bacia_irece, 'Bacia de Irecê (≤100 km)', '★'));
legenda.add(item(CORES_GRUPO.leste_outros, 'Leste / outros', '★'));
legenda.add(item(CORES_GRUPO.norte, 'Norte (Sobradinho)', '★'));
legenda.add(ui.Label('Geoquímica P2O5 (%)', {fontWeight: 'bold', fontSize: '12px', margin: '6px 6px 0 6px'}));
legenda.add(item('e31a1c', '≥ 5', '●'));
legenda.add(item('fd8d3c', '1 – 5', '●'));
legenda.add(item('6baed6', '< 1', '●'));
legenda.add(item('00ffff', 'Grupo Una / Fm. Salitre'));
legenda.add(item('ffff00', 'Área de estudo', '□'));
Map.add(legenda);

// ----------------------------------------------------------------------------
// Inspeção por clique
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'top-right', width: '300px', maxHeight: '420px'}});
info.add(ui.Label('Clique no mapa para ver os valores.', {fontSize: '12px'}));
Map.add(info);
Map.style().set('cursor', 'crosshair');

Map.onClick(function (c) {
  var pt = ee.Geometry.Point(c.lon, c.lat);
  info.clear();
  info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4),
                    {fontWeight: 'bold', fontSize: '12px'}));
  info.add(ui.Label('carregando…', {fontSize: '11px'}));

  var valores = gama.addBands(mag).reduceRegion(ee.Reducer.first(), pt, 30);
  var unidJac = litoJac.filterBounds(pt).aggregate_array('NOME_UNIDA');
  var unidBa  = litoBa.filterBounds(pt).aggregate_array('NOME_UNIDA');
  var raio = pt.buffer(5000);
  var fos5 = fosfato.filterBounds(raio).aggregate_array('toponimia');
  var rec5 = recmin.filterBounds(raio).aggregate_array('subst');

  ee.Dictionary({v: valores, jac: unidJac, ba: unidBa, fos: fos5, rec: rec5})
    .evaluate(function (r) {
      info.remove(info.widgets().get(1));
      var v = r.v || {};
      var linhas = [
        'eU ' + fmt(v.eU_ppm) + ' ppm  ·  eTh ' + fmt(v.eTh_ppm) + ' ppm  ·  K ' + fmt(v.K_pct) + ' %',
        'U/Th ' + fmt(v.U_Th) + '  ·  U/K ' + fmt(v.U_K) + '  ·  Th/K ' + fmt(v.Th_K),
        'CT ' + fmt(v.CT) + '  ·  MAG ' + fmt(v.MAG_nT) + ' nT',
        'Unidade 1:250k: ' + (r.jac.length ? r.jac.join(', ') : '—'),
        'Unidade 1:1M: ' + (r.ba.length ? r.ba.join(', ') : '—'),
        'Fosfato a ≤5 km: ' + (r.fos.length ? r.fos.join(', ') : 'nenhum'),
        'Outras ocorrências a ≤5 km: ' + (r.rec.length ? r.rec.slice(0, 8).join(', ') +
                                         (r.rec.length > 8 ? '…' : '') : 'nenhuma')
      ];
      linhas.forEach(function (t) { info.add(ui.Label(t, {fontSize: '11px', margin: '2px 6px'})); });
    });
});

function fmt(x) { return (x === null || x === undefined) ? '—' : Number(x).toFixed(2); }

// ----------------------------------------------------------------------------
// Resumo no Console
// ----------------------------------------------------------------------------
print('Bandas 1042_gama:', gama.bandNames());
print('Fosfato por grupo:', fosfato.aggregate_histogram('grupo'));
print('Fosfato na Bacia de Irecê:', fosfato.filter(ee.Filter.eq('grupo', 'bacia_irece'))
                                           .aggregate_array('toponimia'));
print('Unidades-alvo (Jacobina 1:250k):', alvoJac.aggregate_array('NOME_UNIDA').distinct());
print('Amostras com P2O5 >= 5 %:', geoq.filter(ee.Filter.gte('p2o5_pct', 5)));
