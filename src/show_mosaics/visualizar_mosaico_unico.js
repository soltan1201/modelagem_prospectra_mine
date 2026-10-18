// ============================================================================
// Mosaico único ASTER (build_mosaico_unico.py) — visualizador para decidir o
// produto final
//
// Fontes (projects/mapbiomas-arida/mine/)
//   ASTER_mosaico_unico_base  ASTER_base_com_swir_<LxCy>  cenas até 2008-S1, SWIR medido
//                             ASTER_base_sem_swir_<LxCy>  cenas 2008-S2 em diante, VNIR + TIR
//   ASTER_mosaico_unico       ASTER_sem_swir_estimado_<LxCy>  base sem_swir + SWIR por RF
//   aster-mosaico-final       mosaico v1 (origem 1–4), para comparação
//   VNIR/SWIR TOA × 1000, TIR em DN; 'camada' 1 = zero-cloud, 2 = preenchimento;
//   'data' AAAAMMDD; 'quality' × 1000; 'normalizacao' 0/1/2.
//
// Perguntas que o painel ajuda a responder
//   1. O SWIR estimado concorda com o medido onde a vegetação é parecida?
//      (Comparação → ΔB08, estatísticas por célula, dispersão)
//   2. Os dois mosaicos são da mesma estação? com_swir usa jul–out; o
//      sem_swir das rodadas atuais usa o ano todo (histograma de meses).
//   3. Quanto de cada célula vem de cena zero-cloud (camada 1)?
//   4. Qual regra de combinação dá o melhor produto? (Prévia combinada)
//
// Nada é calculado ao abrir: marque as camadas ou clique nos botões.
// Clique no mapa: valores dos três mosaicos no ponto.
// ============================================================================

var MINE  = 'projects/mapbiomas-arida/mine/';
var BASE  = ee.ImageCollection(MINE + 'ASTER_mosaico_unico_base');
var COM   = BASE.filter(ee.Filter.eq('tipo', 'com_swir')).mosaic();
var SEM   = BASE.filter(ee.Filter.eq('tipo', 'sem_swir')).mosaic();
var COL_EST = ee.ImageCollection(MINE + 'ASTER_mosaico_unico');
var EST   = COL_EST.mosaic();
var V1    = ee.ImageCollection(MINE + 'aster-mosaico-final').mosaic();

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
var LIMIAR_NDVI = 0.05;   // |ΔNDVI| padrão para "mesma condição de vegetação"

var PAL_VEG  = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];
var PAL_DIV  = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];
var PAL_ANO  = ['440154', '3b528b', '21918c', '5ec962', 'fde725'];
var PAL_MES  = ['313695', '4575b4', '74add1', 'abd9e9', 'e0f3f8', 'ffffbf',
                'fee090', 'fdae61', 'f46d43', 'd73027', 'a50026', '7f3b08'];

// ----------------------------------------------------------------------------
// Derivados
// ----------------------------------------------------------------------------
function ndvi(img) { return img.normalizedDifference(['B3N', 'B02']).rename('ndvi'); }
function ano(img)  { return img.select('data').divide(10000).floor().rename('ano'); }
function mes(img)  { return img.select('data').divide(100).floor().mod(100).rename('mes'); }
function valido(img) { return img.select('B01').mask().gt(0); }

var dNDVI = ndvi(EST).subtract(ndvi(COM)).rename('dNDVI');              // só onde os dois existem
function mesmaVeg(L) { return dNDVI.abs().lte(L); }
// viés relativo (%) do estimado sobre o medido, por banda SWIR
function viesRel(b) {
  return EST.select(b).subtract(COM.select(b)).divide(COM.select(b).max(1)).multiply(100).rename(b);
}

function rgb(img, b) {
  return [img.select(b), {bands: b, min: b.map(function (x) { return FAIXA_BANDA[x][0]; }),
                          max: b.map(function (x) { return FAIXA_BANDA[x][1]; })}];
}

// células com/sem o asset estimado
var cel_est = COL_EST.aggregate_array('celula_nome');
var gridStatus = grid.map(function (f) {
  var ok = cel_est.contains(f.get('nome'));
  return f.set('estilo', {color: ee.Algorithms.If(ok, 'ffffff', 'ff0000'),
                          width: ee.Algorithms.If(ok, 1, 3), fillColor: '00000000'});
});

// ----------------------------------------------------------------------------
// Catálogo de camadas
// ----------------------------------------------------------------------------
var CATALOGO = [
  {grupo: 'Composições', itens: [
    {nome: 'com_swir — VNIR (3N-2-1)',          criar: function () { return rgb(COM, ['B3N', 'B02', 'B01']); }},
    {nome: 'sem_swir — VNIR (3N-2-1)',          criar: function () { return rgb(SEM, ['B3N', 'B02', 'B01']); }},
    {nome: 'v1 final — VNIR (3N-2-1)',          criar: function () { return rgb(V1,  ['B3N', 'B02', 'B01']); }},
    {nome: 'com_swir — SWIR medido (4-6-8)',    criar: function () { return rgb(COM, ['B04', 'B06', 'B08']); }},
    {nome: 'estimado — SWIR RF (4-6-8)',        criar: function () { return rgb(EST, ['B04', 'B06', 'B08']); }},
    {nome: 'v1 final — SWIR (4-6-8)',           criar: function () { return rgb(V1,  ['B04', 'B06', 'B08']); }},
    {nome: 'com_swir — TIR (13-12-10)',         criar: function () { return rgb(COM, ['B13', 'B12', 'B10']); }},
    {nome: 'sem_swir — TIR (13-12-10)',         criar: function () { return rgb(SEM, ['B13', 'B12', 'B10']); }}
  ]},
  {grupo: 'Composição temporal', itens: [
    {nome: 'com_swir — camada (1 zero-cloud, 2 preench.)', criar: function () {
      return [COM.select('camada'), {min: 1, max: 2, palette: ['1a9850', 'fdae61']}]; }},
    {nome: 'sem_swir — camada (1 zero-cloud, 2 preench.)', criar: function () {
      return [SEM.select('camada'), {min: 1, max: 2, palette: ['1a9850', 'fdae61']}]; }},
    {nome: 'com_swir — ano da cena', criar: function () { return [ano(COM), {min: 2000, max: 2008, palette: PAL_ANO}]; }},
    {nome: 'sem_swir — ano da cena', criar: function () { return [ano(SEM), {min: 2008, max: 2020, palette: PAL_ANO}]; }},
    {nome: 'com_swir — mês da cena', criar: function () { return [mes(COM), {min: 1, max: 12, palette: PAL_MES}]; }},
    {nome: 'sem_swir — mês da cena', criar: function () { return [mes(SEM), {min: 1, max: 12, palette: PAL_MES}]; }},
    {nome: 'sem_swir — fora de jul–out', criar: function () {
      var m = mes(SEM); return [m.lt(7).or(m.gt(10)).selfMask(), {palette: ['ff00ff']}]; }},
    {nome: 'com_swir — quality', criar: function () { return [COM.select('quality'), {min: 1000, max: 5500, palette: PAL_VEG}]; }},
    {nome: 'sem_swir — quality', criar: function () { return [SEM.select('quality'), {min: 1000, max: 5500, palette: PAL_VEG}]; }},
    {nome: 'sem_swir — saturado', criar: function () { return [SEM.select('saturado').selfMask(), {palette: ['ff0000']}]; }},
    {nome: 'estimado — normalização (0/1/2)', criar: function () {
      return [EST.select('normalizacao'), {min: 0, max: 2, palette: ['d73027', 'fee08b', '1a9850']}]; }}
  ]},
  {grupo: 'Comparação estimado × medido', itens: [
    {nome: 'NDVI com_swir', criar: function () { return [ndvi(COM), {min: 0.1, max: 0.7, palette: PAL_VEG}]; }},
    {nome: 'NDVI sem_swir', criar: function () { return [ndvi(EST), {min: 0.1, max: 0.7, palette: PAL_VEG}]; }},
    {nome: 'ΔNDVI (sem − com)', criar: function () { return [dNDVI, {min: -0.3, max: 0.3, palette: PAL_DIV}]; }},
    {nome: 'ΔB04 % (estimado − medido)', criar: function () { return [viesRel('B04'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB06 % (estimado − medido)', criar: function () { return [viesRel('B06'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'ΔB08 % (estimado − medido)', criar: function () { return [viesRel('B08'), {min: -30, max: 30, palette: PAL_DIV}]; }},
    {nome: 'carbonate_8_6 medido', criar: function () {
      return [COM.select('B08').divide(COM.select('B06')), {min: 0.58, max: 0.74, palette: PAL_DIV.slice().reverse()}]; }},
    {nome: 'carbonate_8_6 estimado', criar: function () {
      return [EST.select('B08').divide(EST.select('B06')), {min: 0.58, max: 0.74, palette: PAL_DIV.slice().reverse()}]; }}
  ]},
  {grupo: 'Referência', itens: [
    {nome: 'Grid (vermelho = falta estimado)', criar: function () { return [gridStatus.style({styleProperty: 'estilo'}), {}]; }},
    {nome: 'Área de estudo', criar: function () { return [ee.Image().paint(area, 0, 2), {palette: 'ffff00'}]; }}
  ]}
];

// ----------------------------------------------------------------------------
// Prévia combinada: com_swir (medido) × sem_swir_estimado
// ----------------------------------------------------------------------------
var REGRAS = {
  'A · medido onde |ΔNDVI| ≤ limiar, senão estimado': function (L) {
    var ok = mesmaVeg(L).unmask(0);
    return valido(COM).and(valido(EST).not().or(ok));
  },
  'B · medido sempre que existir, estimado nas lacunas': function () { return valido(COM); },
  'C · estimado sempre, medido nas lacunas': function () { return valido(COM).and(valido(EST).not()); },
  'D · medido só na camada 1 (zero-cloud)': function () {
    return valido(COM).and(COM.select('camada').eq(1).unmask(0).or(valido(EST).not()));
  }
};
var BANDAS_PREVIA = BANDAS.concat(['data', 'camada']);

function combinado(regra, L) {
  var usarCom = REGRAS[regra](L);
  var est = EST.select(BANDAS_PREVIA).toInt32().addBands(ee.Image.constant(2).toInt32().rename('fonte'));
  var com = COM.select(BANDAS_PREVIA).toInt32().addBands(ee.Image.constant(1).toInt32().rename('fonte'))
               .updateMask(usarCom);
  return ee.ImageCollection([est, com]).mosaic();     // com_swir por cima onde a regra manda
}

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

var painel = ui.Panel({style: {width: '330px', padding: '6px'}});
painel.add(ui.Label('Mosaico único ASTER · Irecê', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('com_swir (medido, jul–out) × sem_swir_estimado (RF). Marque o que quer ver.',
                    {fontSize: '11px', color: '666666'}));
var lblCont = ui.Label('contando assets…', {fontSize: '11px'});
painel.add(lblCont);
ee.Dictionary({
  com: BASE.filter(ee.Filter.eq('tipo', 'com_swir')).size(),
  sem: BASE.filter(ee.Filter.eq('tipo', 'sem_swir')).size(),
  est: COL_EST.size(),
  falta: grid.aggregate_array('nome').removeAll(cel_est)
}).evaluate(function (r) {
  lblCont.setValue('Células: com_swir ' + r.com + '/25 · sem_swir ' + r.sem + '/25 · estimado ' +
                   r.est + '/25' + (r.falta.length ? '\nFalta estimado: ' + r.falta.join(', ') : ''));
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

// --- prévia combinada (redesenha quando muda regra ou limiar)
painel.add(ui.Label('Prévia combinada', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
var selRegra = ui.Select({items: Object.keys(REGRAS), value: Object.keys(REGRAS)[0],
                          style: {stretch: 'horizontal'}});
var sldLimiar = ui.Slider({min: 0.01, max: 0.3, step: 0.01, value: LIMIAR_NDVI,
                           style: {stretch: 'horizontal'}});
var chkFonte = ui.Checkbox('Fonte (azul medido, laranja estimado)', false);
var chkSwir  = ui.Checkbox('SWIR (4-6-8)', false);
var chkVnir  = ui.Checkbox('VNIR (3N-2-1)', false);
var chkMesma = ui.Checkbox('|ΔNDVI| ≤ limiar', false);
var camadasPrevia = [];

function desenharPrevia() {
  camadasPrevia.forEach(function (l) { Map.layers().remove(l); });
  camadasPrevia = [];
  var L = sldLimiar.getValue();
  var img = combinado(selRegra.getValue(), L);
  var pedidos = [
    [chkVnir,  function () { return rgb(img, ['B3N', 'B02', 'B01']).concat(['Prévia — VNIR']); }],
    [chkSwir,  function () { return rgb(img, ['B04', 'B06', 'B08']).concat(['Prévia — SWIR']); }],
    [chkFonte, function () { return [img.select('fonte'), {min: 1, max: 2, palette: ['2166ac', 'f46d43']}, 'Prévia — fonte']; }],
    [chkMesma, function () { return [mesmaVeg(L).selfMask(), {palette: ['00ffff']}, '|ΔNDVI| ≤ ' + L]; }]
  ];
  pedidos.forEach(function (p) {
    if (!p[0].getValue()) { return; }
    var r = p[1]();
    var l = ui.Map.Layer(r[0], r[1], r[2], true);
    camadasPrevia.push(l);
    Map.layers().add(l);
  });
}
[selRegra, sldLimiar, chkFonte, chkSwir, chkVnir, chkMesma].forEach(function (w) { w.onChange(desenharPrevia); });
painel.add(selRegra);
painel.add(ui.Panel([ui.Label('Limiar |ΔNDVI|', {fontSize: '12px'}), sldLimiar], ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Panel([chkVnir, chkSwir], ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Panel([chkFonte, chkMesma], ui.Panel.Layout.flow('vertical')));

// ----------------------------------------------------------------------------
// Estatísticas e gráficos (sob demanda)
// ----------------------------------------------------------------------------
var saida = ui.Panel();

function estatisticasCelulas() {
  var L = sldLimiar.getValue();
  var ok = mesmaVeg(L);
  var m = mes(SEM);
  var stk = ee.Image.cat([
    valido(COM).unmask(0).rename('cob_com'),
    valido(EST).unmask(0).rename('cob_est'),
    COM.select('camada').eq(1).rename('c1_com'),
    SEM.select('camada').eq(1).rename('c1_sem'),
    m.gte(7).and(m.lte(10)).rename('seco_sem'),
    ndvi(COM).rename('ndvi_com'),
    ndvi(EST).rename('ndvi_est'),
    ok.unmask(0).updateMask(valido(COM).and(valido(EST))).rename('frac_mesma'),
    ee.Image.cat(SWIR.map(viesRel)).updateMask(ok)
  ]);
  var par = ee.Image.cat([COM.select('B08'), EST.select('B08').rename('B08e')]).updateMask(ok);
  var fc = grid.map(function (f) {
    var g = f.geometry();
    var a = stk.reduceRegion({reducer: ee.Reducer.mean(), geometry: g, scale: ESCALA_STATS,
                              maxPixels: 1e9, tileScale: 4});
    var r = par.reduceRegion({reducer: ee.Reducer.pearsonsCorrelation(), geometry: g,
                              scale: ESCALA_STATS, maxPixels: 1e9, tileScale: 4});
    return ee.Feature(null, a).set('nome', f.get('nome'), 'r_B08', r.get('correlation'));
  });
  var cols = [['nome', 'Célula', 0], ['cob_com', 'Cob. com %', 100], ['cob_est', 'Cob. est %', 100],
              ['c1_com', 'Camada 1 com %', 100], ['c1_sem', 'Camada 1 sem %', 100],
              ['seco_sem', 'sem jul–out %', 100], ['ndvi_com', 'NDVI com', 1], ['ndvi_est', 'NDVI sem', 1],
              ['frac_mesma', '|ΔNDVI|≤L %', 100]]
    .concat(SWIR.map(function (b) { return [b, 'viés ' + b + ' %', 1]; }))
    .concat([['r_B08', 'r B08', 1]]);
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
    saida.add(ui.Label('Por célula — viés e r só onde |ΔNDVI| ≤ ' + L, {fontSize: '11px', fontWeight: 'bold'}));
    saida.add(ui.Chart(dt, 'Table', {allowHtml: true, pageSize: 25}));
  });
}

function graficoDispersao() {
  var b = selBandaDisp.getValue();
  var L = sldLimiar.getValue();
  var pts = ee.Image.cat([COM.select(b).rename('medido'), EST.select(b).rename('estimado')])
              .updateMask(mesmaVeg(L))
              .sample({region: area, scale: 120, numPixels: 20000, seed: 42, tileScale: 4})
              .limit(3000);
  saida.clear();
  saida.add(ui.Chart.feature.byFeature(pts, 'medido', 'estimado').setChartType('ScatterChart').setOptions({
    title: b + ' estimado × medido (|ΔNDVI| ≤ ' + L + ')', pointSize: 1, legend: {position: 'none'},
    hAxis: {title: 'medido (TOA × 1000)'}, vAxis: {title: 'estimado (TOA × 1000)'},
    trendlines: {0: {type: 'linear', showR2: true, visibleInLegend: true, color: 'red'}}
  }));
}

function graficoMeses() {
  var img = ee.Image.cat([mes(COM).rename('com_swir'), mes(SEM).rename('sem_swir')]);
  saida.clear();
  saida.add(ui.Chart.image.histogram({image: img, region: area, scale: 1000, minBucketWidth: 1, maxPixels: 1e9})
    .setOptions({title: 'Mês da cena por pixel', hAxis: {title: 'mês'}, vAxis: {title: 'pixels (1 km)'}}));
  var imgAno = ee.Image.cat([ano(COM).rename('com_swir'), ano(SEM).rename('sem_swir')]);
  saida.add(ui.Chart.image.histogram({image: imgAno, region: area, scale: 1000, minBucketWidth: 1, maxPixels: 1e9})
    .setOptions({title: 'Ano da cena por pixel', hAxis: {title: 'ano', format: '####'}, vAxis: {title: 'pixels (1 km)'}}));
}

var selBandaDisp = ui.Select({items: SWIR, value: 'B08'});
painel.add(ui.Label('Análises', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
painel.add(ui.Button({label: 'Estatísticas por célula', onClick: estatisticasCelulas, style: {stretch: 'horizontal'}}));
painel.add(ui.Panel([ui.Button({label: 'Dispersão estimado × medido', onClick: graficoDispersao}), selBandaDisp],
                    ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Button({label: 'Histograma de meses e anos', onClick: graficoMeses, style: {stretch: 'horizontal'}}));
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
  ee.Dictionary({com: ler(COM), sem: ler(SEM), est: ler(EST),
                 v1: V1.select(['origem', 'swir_estimado']).reduceRegion(ee.Reducer.first(), pt, 15)})
    .evaluate(function (r) {
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
      bloco('com_swir', r.com, BANDAS);
      bloco('sem_swir_estimado', r.est, BANDAS, function (v) { return 'norm. ' + f(v.normalizacao, 0); });
      if (!r.est || r.est.B01 === null || r.est.B01 === undefined) {
        bloco('sem_swir (base)', r.sem, VNIR.concat(TIR));
      }
      var v1 = r.v1 || {};
      info.add(ui.Label('v1 final: origem ' + f(v1.origem, 0) +
                        (v1.swir_estimado ? ' (SWIR estimado)' : ' (SWIR medido)'), {fontSize: '11px'}));
    });
});

// ----------------------------------------------------------------------------
ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);
