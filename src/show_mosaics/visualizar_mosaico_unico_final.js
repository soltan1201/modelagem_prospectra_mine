// ============================================================================
// Mosaico único final ASTER (build_mosaico_unico.py etapa 3) — visualizador
//
// Fontes (projects/mapbiomas-arida/mine/)
//   ASTER_mosaico_unico_final  ASTER_mosaico_unico_<LxCy>  com_swir + sem_swir juntos
//   ASTER_mosaico_unico_base   ASTER_base_com_swir_<LxCy>  SWIR medido (referência)
//                              ASTER_base_sem_swir_<LxCy>  VNIR + TIR, sem SWIR
//   VNIR/SWIR TOA × 1000, TIR em DN; 'camada' 1 = zero-cloud, 2 = preenchimento;
//   'data' AAAAMMDD; 'quality' × 1000; 'swir_estimado' 0 medido / 1 estimado;
//   'normalizacao' 0 medido, 1 reta da célula, 2 reta da cena.
//
// Perguntas que o painel ajuda a responder
//   1. Onde o pixel final veio do com_swir (SWIR medido) e onde do sem_swir
//      (SWIR estimado)? Quanto de cada célula é estimado?
//   2. Há emendas visíveis entre medido e estimado, ou entre células?
//   3. O SWIR estimado no final concorda com o medido do com_swir base
//      onde os dois existem? (viés %, dispersão)
//   4. Os índices minerais (carbonato, Al-OH, Fe) ficam contínuos?
//
// Nada é calculado ao abrir: marque as camadas ou clique nos botões.
// Clique no mapa: valores do final e dos dois bases no ponto.
// ============================================================================

var MINE  = 'projects/mapbiomas-arida/mine/';
var COL_FIN = ee.ImageCollection(MINE + 'ASTER_mosaico_unico_final');
var FIN   = COL_FIN.mosaic();
var BASE  = ee.ImageCollection(MINE + 'ASTER_mosaico_unico_base');
var COM   = BASE.filter(ee.Filter.eq('tipo', 'com_swir')).mosaic();
var SEM   = BASE.filter(ee.Filter.eq('tipo', 'sem_swir')).mosaic();

var area = ee.Geometry.Rectangle([-44.460483519818084, -14.790793773956535,
                                  -39.203525512005584, -9.416525285369607]);
var grid = ee.FeatureCollection(MINE + 'grid_irece_5x5_solape100');

var VNIR = ['B01', 'B02', 'B3N'];
var SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09'];
var TIR  = ['B10', 'B11', 'B12', 'B13', 'B14'];
var BANDAS = VNIR.concat(SWIR, TIR);
// faixas por banda (VNIR/SWIR × 1000; TIR em DN) — as do mosaico v1
var FAIXA_BANDA = {
  B01: [69, 176], B02: [51, 204], B3N: [155, 332], B04: [167, 362], B05: [77, 232],
  B06: [75, 234], B07: [71, 228], B08: [45, 168], B09: [34, 115],
  B10: [1124, 1635], B11: [1228, 1803], B12: [1332, 1931], B13: [1580, 2203], B14: [1644, 2259]
};
var ESCALA_STATS = 300;   // m, estatísticas por célula

var PAL_VEG  = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];
var PAL_DIV  = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];
var PAL_ANO  = ['440154', '3b528b', '21918c', '5ec962', 'fde725'];
var PAL_MES  = ['313695', '4575b4', '74add1', 'abd9e9', 'e0f3f8', 'ffffbf',
                'fee090', 'fdae61', 'f46d43', 'd73027', 'a50026', '7f3b08'];
var PAL_IDX  = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];

// ----------------------------------------------------------------------------
// Derivados
// ----------------------------------------------------------------------------
function ndvi(img) { return img.normalizedDifference(['B3N', 'B02']).rename('ndvi'); }
function ano(img)  { return img.select('data').divide(10000).floor().rename('ano'); }
function mes(img)  { return img.select('data').divide(100).floor().mod(100).rename('mes'); }
function valido(img) { return img.select('B01').mask().gt(0); }
function razao(img, a, b) { return img.select(a).divide(img.select(b).max(1)); }

var EST_FIN = FIN.select('swir_estimado').eq(1);     // pixel do sem_swir
var FIN_COM = FIN.updateMask(EST_FIN.not());          // só pixels com SWIR medido
var FIN_SEM = FIN.updateMask(EST_FIN);                // só pixels com SWIR estimado

// viés relativo (%) do SWIR final estimado sobre o medido do com_swir base
function viesRel(b) {
  return FIN_SEM.select(b).subtract(COM.select(b)).divide(COM.select(b).max(1))
                .multiply(100).rename(b);
}

// índices minerais ASTER (Ninomiya / Rowan & Mars); faixas = percentis 2–98 do final a 1 km
var INDICES = {
  'carbonato (B06+B09)/B08': {f: function (i) { return i.select('B06').add(i.select('B09')).divide(i.select('B08').max(1)); },
                             min: 2.0, max: 2.43},
  'Al-OH (B05+B07)/B06':     {f: function (i) { return i.select('B05').add(i.select('B07')).divide(i.select('B06').max(1)); },
                             min: 1.9, max: 2.04},
  'Fe-OH/Mg-OH (B06+B09)/(B07+B08)': {f: function (i) {
                               return i.select('B06').add(i.select('B09')).divide(i.select('B07').add(i.select('B08')).max(1)); },
                             min: 0.84, max: 0.96},
  'ferro férrico B02/B01':   {f: function (i) { return razao(i, 'B02', 'B01'); }, min: 0.71, max: 1.36},
  'quartzo TIR B14/B12':     {f: function (i) { return razao(i, 'B14', 'B12'); }, min: 1.15, max: 1.25}
};

function rgb(img, b) {
  return [img.select(b), {bands: b, min: b.map(function (x) { return FAIXA_BANDA[x][0]; }),
                          max: b.map(function (x) { return FAIXA_BANDA[x][1]; })}];
}

var cel_fin = COL_FIN.aggregate_array('celula_nome');
var gridStatus = grid.map(function (f) {
  var ok = cel_fin.contains(f.get('nome'));
  return f.set('estilo', {color: ee.Algorithms.If(ok, 'ffffff', 'ff0000'),
                          width: ee.Algorithms.If(ok, 1, 3), fillColor: '00000000'});
});

// ----------------------------------------------------------------------------
// Catálogo de camadas
// ----------------------------------------------------------------------------
var CATALOGO = [
  {grupo: 'Mosaico final', itens: [
    {nome: 'final — VNIR (3N-2-1)',  criar: function () { return rgb(FIN, ['B3N', 'B02', 'B01']); }},
    {nome: 'final — SWIR (4-6-8)',   criar: function () { return rgb(FIN, ['B04', 'B06', 'B08']); }},
    {nome: 'final — SWIR (4-7-9)',   criar: function () { return rgb(FIN, ['B04', 'B07', 'B09']); }},
    {nome: 'final — TIR (13-12-10)', criar: function () { return rgb(FIN, ['B13', 'B12', 'B10']); }},
    {nome: 'final — NDVI',           criar: function () { return [ndvi(FIN), {min: 0.1, max: 0.7, palette: PAL_VEG}]; }}
  ]},
  {grupo: 'Final separado: com × sem SWIR', itens: [
    {nome: 'pixels com SWIR medido — SWIR (4-6-8)',   criar: function () { return rgb(FIN_COM, ['B04', 'B06', 'B08']); }},
    {nome: 'pixels com SWIR estimado — SWIR (4-6-8)', criar: function () { return rgb(FIN_SEM, ['B04', 'B06', 'B08']); }},
    {nome: 'pixels com SWIR medido — VNIR',           criar: function () { return rgb(FIN_COM, ['B3N', 'B02', 'B01']); }},
    {nome: 'pixels com SWIR estimado — VNIR',         criar: function () { return rgb(FIN_SEM, ['B3N', 'B02', 'B01']); }}
  ]},
  {grupo: 'Origem do pixel', itens: [
    {nome: 'swir_estimado (azul medido, laranja estimado)', criar: function () {
      return [FIN.select('swir_estimado'), {min: 0, max: 1, palette: ['2166ac', 'f46d43']}]; }},
    {nome: 'normalizacao (0 medido, 1 célula, 2 cena)', criar: function () {
      return [FIN.select('normalizacao'), {min: 0, max: 2, palette: ['2166ac', 'fee08b', '1a9850']}]; }},
    {nome: 'camada (1 zero-cloud, 2 preench.)', criar: function () {
      return [FIN.select('camada'), {min: 1, max: 2, palette: ['1a9850', 'fdae61']}]; }},
    {nome: 'ano da cena', criar: function () { return [ano(FIN), {min: 2000, max: 2020, palette: PAL_ANO}]; }},
    {nome: 'mês da cena', criar: function () { return [mes(FIN), {min: 1, max: 12, palette: PAL_MES}]; }},
    {nome: 'fora de jul–out', criar: function () {
      var m = mes(FIN); return [m.lt(7).or(m.gt(10)).selfMask(), {palette: ['ff00ff']}]; }},
    {nome: 'quality', criar: function () { return [FIN.select('quality'), {min: 1000, max: 5500, palette: PAL_VEG}]; }}
  ]},
  {grupo: 'Índices minerais (final)', itens: Object.keys(INDICES).map(function (k) {
    return {nome: k, criar: function () {
      return [INDICES[k].f(FIN), {min: INDICES[k].min, max: INDICES[k].max, palette: PAL_IDX}]; }};
  })},
  {grupo: 'Estimado (final) × medido (com_swir base)', itens: [
    {nome: 'com_swir base — SWIR medido (4-6-8)', criar: function () { return rgb(COM, ['B04', 'B06', 'B08']); }},
    {nome: 'sem_swir base — VNIR (3N-2-1)',       criar: function () { return rgb(SEM, ['B3N', 'B02', 'B01']); }},
    {nome: 'ΔNDVI (final estimado − com_swir)', criar: function () {
      return [ndvi(FIN_SEM).subtract(ndvi(COM)), {min: -0.3, max: 0.3, palette: PAL_DIV}]; }},
    {nome: 'ΔB04 % (estimado − medido)', criar: function () { return [viesRel('B04'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB06 % (estimado − medido)', criar: function () { return [viesRel('B06'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB08 % (estimado − medido)', criar: function () { return [viesRel('B08'), {min: -30, max: 30, palette: PAL_DIV}]; }}
  ]},
  {grupo: 'Referência', itens: [
    {nome: 'Grid (vermelho = falta final)', criar: function () { return [gridStatus.style({styleProperty: 'estilo'}), {}]; }},
    {nome: 'Área de estudo', criar: function () { return [ee.Image().paint(area, 0, 2), {palette: 'ffff00'}]; }}
  ]}
];

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var camadas = {};   // nome → ui.Map.Layer (criada na primeira marcação)

function alternar(item, ligado) {
  if (!camadas[item.nome]) {
    if (!ligado) { return; }
    var r = item.criar();
    camadas[item.nome] = ui.Map.Layer(r[0], r[1], item.nome, true);
    Map.layers().add(camadas[item.nome]);
  } else {
    camadas[item.nome].setShown(ligado);
  }
}

var painel = ui.Panel({style: {width: '340px', padding: '6px'}});
painel.add(ui.Label('Mosaico único final ASTER · Irecê', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('Etapa 3: com_swir (SWIR medido) + sem_swir (SWIR estimado por RF). ' +
                    'Marque o que quer ver.', {fontSize: '11px', color: '666666'}));
var lblCont = ui.Label('contando assets…', {fontSize: '11px'});
painel.add(lblCont);
ee.Dictionary({
  n: COL_FIN.size(),
  falta: grid.aggregate_array('nome').removeAll(cel_fin)
}).evaluate(function (r) {
  lblCont.setValue('Células no final: ' + r.n + '/25' +
                   (r.falta.length ? '\nFalta: ' + r.falta.join(', ') : ''));
  lblCont.style().set('whiteSpace', 'pre-wrap');
});

CATALOGO.forEach(function (g) {
  var corpo = ui.Panel({style: {shown: false, margin: '0 0 0 8px'}});
  var botao = ui.Button({label: '▸ ' + g.grupo, style: {stretch: 'horizontal', margin: '4px 0 0 0'}});
  botao.onClick(function () {
    var aberto = !corpo.style().get('shown');
    corpo.style().set('shown', aberto);
    botao.setLabel((aberto ? '▾ ' : '▸ ') + g.grupo);
  });
  g.itens.forEach(function (item) {
    corpo.add(ui.Checkbox({label: item.nome, value: false, style: {fontSize: '12px', margin: '1px 0'},
                           onChange: function (v) { alternar(item, v); }}));
  });
  painel.add(botao);
  painel.add(corpo);
});

// ----------------------------------------------------------------------------
// Estatísticas e gráficos (sob demanda)
// ----------------------------------------------------------------------------
var saida = ui.Panel();

function estatisticasCelulas() {
  var stk = ee.Image.cat([
    EST_FIN.rename('frac_est'),
    FIN.select('camada').eq(1).rename('c1'),
    FIN.select('normalizacao').eq(2).updateMask(EST_FIN).rename('norm_cena'),
    mes(FIN).gte(7).and(mes(FIN).lte(10)).rename('seco'),
    ndvi(FIN).rename('ndvi'),
    ee.Image.cat(['B04', 'B06', 'B08'].map(viesRel))
  ]);
  var fc = COL_FIN.map(function (img) {
    var g = grid.filter(ee.Filter.eq('nome', img.get('celula_nome'))).geometry();
    var a = stk.reduceRegion({reducer: ee.Reducer.mean(), geometry: g, scale: ESCALA_STATS,
                              maxPixels: 1e9, tileScale: 4});
    return ee.Feature(null, a).set('nome', img.get('celula_nome'),
                                   'n_reta', img.get('n_cenas_reta_propria'));
  }).sort('nome');
  var cols = [['nome', 'Célula', 0], ['frac_est', 'SWIR estimado %', 100], ['c1', 'Camada 1 %', 100],
              ['norm_cena', 'Reta da cena % (est.)', 100], ['n_reta', 'Cenas c/ reta', 0],
              ['seco', 'jul–out %', 100], ['ndvi', 'NDVI', 1],
              ['B04', 'viés B04 %', 1], ['B06', 'viés B06 %', 1], ['B08', 'viés B08 %', 1]];
  saida.clear();
  saida.add(ui.Label('Calculando 25 células a ' + ESCALA_STATS + ' m (pode levar alguns minutos)…',
                     {fontSize: '11px'}));
  fc.evaluate(function (res, erro) {
    saida.clear();
    if (erro) { saida.add(ui.Label('Erro: ' + erro, {color: 'red', fontSize: '11px'})); return; }
    var dt = [cols.map(function (c) { return c[1]; })];
    res.features.forEach(function (f) {
      var p = f.properties;
      dt.push(cols.map(function (c) {
        var v = p[c[0]];
        if (c[2] === 0) { return v; }
        return (v === null || v === undefined) ? null : Math.round(v * c[2] * 1000) / 1000;
      }));
    });
    saida.add(ui.Label('Por célula — viés: SWIR estimado do final × medido do com_swir base',
                       {fontSize: '11px', fontWeight: 'bold'}));
    saida.add(ui.Chart(dt, 'Table', {allowHtml: true, pageSize: 25}));
  });
}

function graficoDispersao() {
  var b = selBandaDisp.getValue();
  var pts = ee.Image.cat([COM.select(b).rename('medido'), FIN_SEM.select(b).rename('estimado')])
              .sample({region: area, scale: 120, numPixels: 20000, seed: 42, tileScale: 4})
              .limit(3000);
  saida.clear();
  saida.add(ui.Chart.feature.byFeature(pts, 'medido', 'estimado').setChartType('ScatterChart').setOptions({
    title: b + ' — final estimado × com_swir medido', pointSize: 1, legend: {position: 'none'},
    hAxis: {title: 'medido (TOA × 1000)'}, vAxis: {title: 'estimado (TOA × 1000)'},
    trendlines: {0: {type: 'linear', showR2: true, visibleInLegend: true, color: 'red'}}
  }));
}

function graficoHistogramas() {
  var b = selBandaDisp.getValue();
  var img = ee.Image.cat([FIN_COM.select(b).rename('medido'), FIN_SEM.select(b).rename('estimado')]);
  saida.clear();
  saida.add(ui.Chart.image.histogram({image: img, region: area, scale: 300, maxBuckets: 120, maxPixels: 1e9})
    .setOptions({title: b + ' no final: pixels medidos × estimados', hAxis: {title: 'TOA × 1000'},
                 vAxis: {title: 'pixels (300 m)'}}));
  var imgT = ee.Image.cat([ano(FIN).rename('ano')]);
  saida.add(ui.Chart.image.histogram({image: imgT, region: area, scale: 1000, minBucketWidth: 1, maxPixels: 1e9})
    .setOptions({title: 'Ano da cena por pixel (final)', hAxis: {title: 'ano', format: '####'},
                 vAxis: {title: 'pixels (1 km)'}}));
}

var selBandaDisp = ui.Select({items: SWIR, value: 'B08'});
painel.add(ui.Label('Análises', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
painel.add(ui.Button({label: 'Estatísticas por célula', onClick: estatisticasCelulas, style: {stretch: 'horizontal'}}));
painel.add(ui.Panel([ui.Label('Banda', {fontSize: '12px'}), selBandaDisp], ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Button({label: 'Dispersão estimado × medido', onClick: graficoDispersao, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'Histograma medido × estimado + anos', onClick: graficoHistogramas,
                      style: {stretch: 'horizontal'}}));
painel.add(saida);

// ----------------------------------------------------------------------------
// Inspeção por clique
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '340px', maxHeight: '380px'}});
info.add(ui.Label('Clique no mapa para ver os valores.', {fontSize: '11px'}));
Map.add(info);
Map.style().set('cursor', 'crosshair');

Map.onClick(function (c) {
  var pt = ee.Geometry.Point(c.lon, c.lat);
  info.clear();
  info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontWeight: 'bold', fontSize: '12px'}));
  info.add(ui.Label('carregando…', {fontSize: '11px'}));
  function ler(img) { return img.addBands(ndvi(img)).reduceRegion(ee.Reducer.first(), pt, 15); }
  ee.Dictionary({fin: ler(FIN), com: ler(COM), sem: ler(SEM)}).evaluate(function (r) {
    info.remove(info.widgets().get(1));
    function f(x, d) { return (x === null || x === undefined) ? '—' : Number(x).toFixed(d); }
    function bloco(titulo, v, bandas, extra) {
      v = v || {};
      if (v.B01 === null || v.B01 === undefined) {
        info.add(ui.Label(titulo + ': sem dado', {fontSize: '11px', fontWeight: 'bold'}));
        return;
      }
      info.add(ui.Label(titulo + ' — data ' + f(v.data, 0) + ' · camada ' + f(v.camada, 0) +
                        ' · NDVI ' + f(v.ndvi, 3) + (extra ? ' · ' + extra(v) : ''),
                        {fontSize: '11px', fontWeight: 'bold'}));
      info.add(ui.Label(bandas.map(function (b) { return b + ' ' + f(v[b], 0); }).join(' · '),
                        {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    }
    bloco('final', r.fin, BANDAS, function (v) {
      return (v.swir_estimado ? 'SWIR estimado' : 'SWIR medido') + ' · norm. ' + f(v.normalizacao, 0);
    });
    bloco('com_swir base', r.com, BANDAS);
    bloco('sem_swir base', r.sem, VNIR.concat(TIR));
  });
});

// ----------------------------------------------------------------------------
ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);
