// ============================================================================
// Mosaicos sazonais secos ASTER (build_mosaico_sazonal_seco.py) — visualizador
//
// Fontes (FONTE abaixo: 'arida' depois de copiar, 'brazil' = export temporário)
//   ASTER_sazonal  ASTER_sazonal_A_seco_quality_<LxCy>  2000-03 a 2008-03, VNIR+SWIR+TIR
//                  ASTER_sazonal_B_seco_quality_<LxCy>  2009–2020, VNIR+TIR (sem SWIR)
//   VNIR/SWIR TOA × 1000, TIR em DN; 'n_obs' = nº de cenas válidas no pixel;
//   propriedades 'periodo' (A/B), 'celula_nome', 'num_images', 'num_zero_cloud'.
//   Pontos pareados (pontos_sazonal_seco_pareados.py):
//   points_sazonal_seco_A / points_sazonal_seco_B em projects/mapbiomas-arida/mine.
//
// Perguntas que o painel ajuda a responder
//   1. Quais células têm mosaico A e B? Onde há poucas cenas (n_obs baixo)?
//   2. O B (2009–2020) fica radiometricamente parecido com o A nas bandas
//      comuns (VNIR, TIR)? Onde a diferença é grande (mudança real ou emenda)?
//   3. Há emendas visíveis entre células em cada período?
//   4. Onde caíram os pontos de treino do SWIR?
//
// Nada é calculado ao abrir: marque as camadas ou clique nos botões.
// Clique no mapa: valores do A e do B no ponto.
// ============================================================================

var FONTE = 'arida';   // 'arida' | 'brazil' | 'harmonizado' (harmonizar_sazonal_modis.py)
var PASTAS = {
  arida:  'projects/mapbiomas-arida/mine/ASTER_sazonal',
  harmonizado: 'projects/mapbiomas-arida/mine/ASTER_sazonal_harmonizado',
  brazil: 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA/ASTER_sazonal'
};
var MINE  = 'projects/mapbiomas-arida/mine/';
var COL   = ee.ImageCollection(PASTAS[FONTE]);
var COL_A = COL.filter(ee.Filter.eq('periodo', 'A'));
var COL_B = COL.filter(ee.Filter.eq('periodo', 'B'));
var A     = COL_A.mosaic();
var B     = COL_B.mosaic();
var PTS_A = MINE + 'points_sazonal_seco_A';
var PTS_B = MINE + 'points_sazonal_seco_B';

var area = ee.Geometry.Rectangle([-44.460483519818084, -14.790793773956535,
                                  -39.203525512005584, -9.416525285369607]);
var grid = ee.FeatureCollection(MINE + 'grid_irece_5x5_solape100');

var VNIR = ['B01', 'B02', 'B3N'];
var SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09'];
var TIR  = ['B10', 'B11', 'B12', 'B13', 'B14'];
var COMUNS = VNIR.concat(TIR);
// faixas por banda (VNIR/SWIR × 1000; TIR em DN) — as do mosaico v1
var FAIXA_BANDA = {
  B01: [69, 176], B02: [51, 204], B3N: [155, 332], B04: [167, 362], B05: [77, 232],
  B06: [75, 234], B07: [71, 228], B08: [45, 168], B09: [34, 115],
  B10: [1124, 1635], B11: [1228, 1803], B12: [1332, 1931], B13: [1580, 2203], B14: [1644, 2259]
};
var ESCALA_STATS = 300;   // m, estatísticas por célula

var PAL_VEG = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];
var PAL_DIV = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];
var PAL_OBS = ['d73027', 'fc8d59', 'fee08b', 'd9ef8b', '91cf60', '1a9850'];
var PAL_IDX = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];

// ----------------------------------------------------------------------------
// Derivados
// ----------------------------------------------------------------------------
function ndvi(img) { return img.normalizedDifference(['B3N', 'B02']).rename('ndvi'); }
function razao(img, a, b) { return img.select(a).divide(img.select(b).max(1)); }

// diferença relativa (%) B − A numa banda comum
function difRel(b) {
  return B.select(b).subtract(A.select(b)).divide(A.select(b).max(1)).multiply(100).rename(b);
}

// índices minerais ASTER (Ninomiya / Rowan & Mars); só o A tem SWIR
var INDICES = {
  'carbonato (B06+B09)/B08': {f: function (i) { return i.select('B06').add(i.select('B09')).divide(i.select('B08').max(1)); },
                             min: 2.0, max: 2.43},
  'Al-OH (B05+B07)/B06':     {f: function (i) { return i.select('B05').add(i.select('B07')).divide(i.select('B06').max(1)); },
                             min: 1.9, max: 2.04},
  'Fe-OH/Mg-OH (B06+B09)/(B07+B08)': {f: function (i) {
                               return i.select('B06').add(i.select('B09')).divide(i.select('B07').add(i.select('B08')).max(1)); },
                             min: 0.84, max: 0.96}
};
var INDICES_COMUNS = {
  'ferro férrico B02/B01': {f: function (i) { return razao(i, 'B02', 'B01'); }, min: 0.71, max: 1.36},
  'quartzo TIR B14/B12':   {f: function (i) { return razao(i, 'B14', 'B12'); }, min: 1.15, max: 1.25}
};

function rgb(img, b) {
  return [img.select(b), {bands: b, min: b.map(function (x) { return FAIXA_BANDA[x][0]; }),
                          max: b.map(function (x) { return FAIXA_BANDA[x][1]; })}];
}

// grid: vermelho = falta A e B, laranja = falta um dos dois, branco = completo
var cel_A = COL_A.aggregate_array('celula_nome');
var cel_B = COL_B.aggregate_array('celula_nome');
var gridStatus = grid.map(function (f) {
  var n = ee.Number(cel_A.contains(f.get('nome'))).add(ee.Number(cel_B.contains(f.get('nome'))));
  var cor = ee.List(['ff0000', 'ff9900', 'ffffff']).get(n);
  return f.set('estilo', {color: cor, width: ee.Algorithms.If(n.eq(2), 1, 3), fillColor: '00000000'});
});

function pontos(asset, cor) {
  return [ee.FeatureCollection(asset).style({color: cor, pointSize: 2, width: 0}), {}];
}

// ----------------------------------------------------------------------------
// Catálogo de camadas
// ----------------------------------------------------------------------------
function camadasPeriodo(rot, img, comSwir) {
  var itens = [
    {nome: rot + ' — VNIR (3N-2-1)',  criar: function () { return rgb(img, ['B3N', 'B02', 'B01']); }}
  ];
  if (comSwir) {
    itens.push({nome: rot + ' — SWIR (4-6-8)', criar: function () { return rgb(img, ['B04', 'B06', 'B08']); }});
    itens.push({nome: rot + ' — SWIR (4-7-9)', criar: function () { return rgb(img, ['B04', 'B07', 'B09']); }});
  }
  return itens.concat([
    {nome: rot + ' — TIR (13-12-10)', criar: function () { return rgb(img, ['B13', 'B12', 'B10']); }},
    {nome: rot + ' — NDVI',           criar: function () { return [ndvi(img), {min: 0.1, max: 0.7, palette: PAL_VEG}]; }},
    {nome: rot + ' — n_obs (cenas válidas)', criar: function () {
      return [img.select('n_obs'), {min: 1, max: comSwir ? 15 : 25, palette: PAL_OBS}]; }},
    {nome: rot + ' — n_obs = 1 (uma cena só)', criar: function () {
      return [img.select('n_obs').eq(1).selfMask(), {palette: ['ff00ff']}]; }}
  ]);
}

var CATALOGO = [
  {grupo: 'Período A (2000–2008, SWIR medido)', itens: camadasPeriodo('A', A, true)},
  {grupo: 'Período B (2009–2020, sem SWIR)',    itens: camadasPeriodo('B', B, false)},
  {grupo: 'B × A (bandas comuns)', itens: [
    {nome: 'ΔNDVI (B − A)', criar: function () {
      return [ndvi(B).subtract(ndvi(A)), {min: -0.3, max: 0.3, palette: PAL_DIV}]; }},
    {nome: 'ΔB01 % (B − A)', criar: function () { return [difRel('B01'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB02 % (B − A)', criar: function () { return [difRel('B02'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB3N % (B − A)', criar: function () { return [difRel('B3N'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB13 % (B − A)', criar: function () { return [difRel('B13'), {min: -10, max: 10, palette: PAL_DIV}]; }},
    {nome: 'só A (falta B)', criar: function () {
      return [A.select('B01').mask().and(B.select('B01').mask().not()).selfMask(), {palette: ['ff0000']}]; }},
    {nome: 'só B (falta A)', criar: function () {
      return [B.select('B01').mask().and(A.select('B01').mask().not()).selfMask(), {palette: ['0000ff']}]; }}
  ]},
  {grupo: 'Índices minerais', itens: Object.keys(INDICES).map(function (k) {
    return {nome: 'A — ' + k, criar: function () {
      return [INDICES[k].f(A), {min: INDICES[k].min, max: INDICES[k].max, palette: PAL_IDX}]; }};
  }).concat([].concat.apply([], Object.keys(INDICES_COMUNS).map(function (k) {
    return [['A', A], ['B', B]].map(function (p) {
      return {nome: p[0] + ' — ' + k, criar: function () {
        return [INDICES_COMUNS[k].f(p[1]), {min: INDICES_COMUNS[k].min, max: INDICES_COMUNS[k].max,
                                            palette: PAL_IDX}]; }};
    });
  })))},
  {grupo: 'Pontos pareados (treino SWIR)', itens: [
    {nome: 'pontos A (etapa 1)', criar: function () { return pontos(PTS_A, 'ffff00'); }},
    {nome: 'pontos B (etapa 2, com dado no B)', criar: function () { return pontos(PTS_B, '00ffff'); }}
  ]},
  {grupo: 'Referência', itens: [
    {nome: 'Grid (vermelho falta A e B, laranja falta um)', criar: function () {
      return [gridStatus.style({styleProperty: 'estilo'}), {}]; }},
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
painel.add(ui.Label('Mosaicos sazonais secos ASTER · Irecê', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('jun–out, qualityMosaic v5. A = 2000–2008 (SWIR medido), B = 2009–2020 ' +
                    '(SWIR a estimar). Fonte: ' + PASTAS[FONTE], {fontSize: '11px', color: '666666'}));
var lblCont = ui.Label('contando assets…', {fontSize: '11px'});
painel.add(lblCont);
ee.Dictionary({
  nA: COL_A.size(), nB: COL_B.size(),
  faltaA: grid.aggregate_array('nome').removeAll(cel_A),
  faltaB: grid.aggregate_array('nome').removeAll(cel_B)
}).evaluate(function (r, erro) {
  if (erro) { lblCont.setValue('Erro ao ler ' + PASTAS[FONTE] + ': ' + erro); return; }
  lblCont.setValue('Células A: ' + r.nA + '/25' + (r.faltaA.length ? ' — falta ' + r.faltaA.join(', ') : '') +
                   '\nCélulas B: ' + r.nB + '/25' + (r.faltaB.length ? ' — falta ' + r.faltaB.join(', ') : ''));
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
    A.select('n_obs').rename('obsA'), B.select('n_obs').rename('obsB'),
    ndvi(A).rename('ndviA'), ndvi(B).rename('ndviB'),
    A.select('B01').mask().rename('cobA'), B.select('B01').mask().rename('cobB'),
    difRel('B3N'), difRel('B13')
  ]).unmask(0, false);
  var fc = grid.map(function (f) {
    var nome = f.get('nome');
    var a = stk.reduceRegion({reducer: ee.Reducer.mean(), geometry: f.geometry(), scale: ESCALA_STATS,
                              maxPixels: 1e9, tileScale: 4});
    var imA = COL_A.filter(ee.Filter.eq('celula_nome', nome)).first();
    var imB = COL_B.filter(ee.Filter.eq('celula_nome', nome)).first();
    return ee.Feature(null, a).set('nome', nome,
      'cenasA', ee.Algorithms.If(imA, ee.Image(imA).get('num_images'), null),
      'cenasB', ee.Algorithms.If(imB, ee.Image(imB).get('num_images'), null));
  }).sort('nome');
  var cols = [['nome', 'Célula', 0], ['cenasA', 'Cenas A', 0], ['cenasB', 'Cenas B', 0],
              ['cobA', 'Cobertura A %', 100], ['cobB', 'Cobertura B %', 100],
              ['obsA', 'n_obs A', 1], ['obsB', 'n_obs B', 1],
              ['ndviA', 'NDVI A', 1], ['ndviB', 'NDVI B', 1],
              ['B3N', 'ΔB3N %', 1], ['B13', 'ΔB13 %', 1]];
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
    saida.add(ui.Label('Por célula — n_obs e NDVI com 0 fora da cobertura; Δ = (B − A)/A',
                       {fontSize: '11px', fontWeight: 'bold'}));
    saida.add(ui.Chart(dt, 'Table', {allowHtml: true, pageSize: 25}));
  });
}

function graficoDispersao() {
  var b = selBanda.getValue();
  var pts = ee.Image.cat([A.select(b).rename('A'), B.select(b).rename('B')])
              .sample({region: area, scale: 120, numPixels: 20000, seed: 42, tileScale: 4})
              .limit(3000);
  saida.clear();
  saida.add(ui.Chart.feature.byFeature(pts, 'A', 'B').setChartType('ScatterChart').setOptions({
    title: b + ' — B (2009–2020) × A (2000–2008)', pointSize: 1, legend: {position: 'none'},
    hAxis: {title: 'A'}, vAxis: {title: 'B'},
    trendlines: {0: {type: 'linear', showR2: true, visibleInLegend: true, color: 'red'}}
  }));
}

function graficoHistogramas() {
  var b = selBanda.getValue();
  var img = ee.Image.cat([A.select(b).rename('A'), B.select(b).rename('B')]);
  saida.clear();
  saida.add(ui.Chart.image.histogram({image: img, region: area, scale: 300, maxBuckets: 120, maxPixels: 1e9})
    .setOptions({title: b + ': A × B', hAxis: {title: TIR.indexOf(b) >= 0 ? 'DN' : 'TOA × 1000'},
                 vAxis: {title: 'pixels (300 m)'}}));
  var obs = ee.Image.cat([A.select('n_obs').rename('A'), B.select('n_obs').rename('B')]);
  saida.add(ui.Chart.image.histogram({image: obs, region: area, scale: 1000, minBucketWidth: 1, maxPixels: 1e9})
    .setOptions({title: 'n_obs (cenas válidas por pixel)', hAxis: {title: 'n_obs'},
                 vAxis: {title: 'pixels (1 km)'}}));
}

var selBanda = ui.Select({items: COMUNS, value: 'B3N'});
painel.add(ui.Label('Análises', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
painel.add(ui.Button({label: 'Estatísticas por célula', onClick: estatisticasCelulas, style: {stretch: 'horizontal'}}));
painel.add(ui.Panel([ui.Label('Banda comum', {fontSize: '12px'}), selBanda], ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Button({label: 'Dispersão B × A', onClick: graficoDispersao, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'Histograma A × B + n_obs', onClick: graficoHistogramas,
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
  var cel = grid.filterBounds(pt).aggregate_array('nome');
  ee.Dictionary({a: ler(A), b: ler(B), cel: cel}).evaluate(function (r) {
    info.remove(info.widgets().get(1));
    function f(x, d) { return (x === null || x === undefined) ? '—' : Number(x).toFixed(d); }
    info.add(ui.Label('Célula: ' + (r.cel.length ? r.cel.join(', ') : '—'), {fontSize: '11px'}));
    function bloco(titulo, v, bandas) {
      v = v || {};
      if (v.B01 === null || v.B01 === undefined) {
        info.add(ui.Label(titulo + ': sem dado', {fontSize: '11px', fontWeight: 'bold'}));
        return;
      }
      info.add(ui.Label(titulo + ' — n_obs ' + f(v.n_obs, 0) + ' · NDVI ' + f(v.ndvi, 3),
                        {fontSize: '11px', fontWeight: 'bold'}));
      info.add(ui.Label(bandas.map(function (b) { return b + ' ' + f(v[b], 0); }).join(' · '),
                        {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    }
    bloco('A (2000–2008)', r.a, VNIR.concat(SWIR, TIR));
    bloco('B (2009–2020)', r.b, COMUNS);
  });
});

// ----------------------------------------------------------------------------
ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);
