// ============================================================================
// Favorabilidade RF × GTB (exportar_favorabilidade_ml.py) — visualizador
//
// Fonte: projects/mapbiomas-arida/mine/favorabilidade_ml
//   uma imagem por <modelo>_<alvo>_<mosaico>_<LxCy>; banda 'probability' =
//   p × 10000 (Int16). Aqui p volta para 0–1.
//   Treino: todos os pontos do alvo + fundo 3 × presença (n_presenca,
//   n_fundo nas propriedades) → p NÃO é a chance real de ocorrência; vale
//   para ordenar a área. Com 1 presença : 3 fundo no treino, p ≈ 0,25 é o
//   "neutro" do modelo.
//
// Painel
//   Mapas: p do RF, p do GTB, diferença GTB − RF, concordância acima de um limiar.
//   Histogramas (um por modelo) e curva área × limiar: intervalo dos valores
//   de p e como a área cai quando o limiar sobe.
//   Estatísticas: mín., percentis, média, máx. por modelo.
//   p nos pontos: distribuição de p nas ocorrências × na paisagem.
// Clique no mapa: p dos dois modelos no ponto.
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
var COL = ee.ImageCollection(MINE + 'favorabilidade_ml');
var GRID = ee.FeatureCollection(MINE + 'grid_irece_5x5');
var PONTOS = ee.FeatureCollection(MINE + 'pontos_referencia');
var MOSAICOS = {
  'aster-mosaico-final': MINE + 'aster-mosaico-final',
  'ASTER_mosaico_unico_final': MINE + 'ASTER_mosaico_unico_final'
};
var MODELOS = ['RF', 'GTB'];
var COR = {RF: '#222222', GTB: '#17becf'};
var ESCALA_HIST = 300;     // m — histogramas e estatísticas (área toda)
var N_CLASSES = 50;        // classes do histograma em 0–1
var regiao = GRID.geometry().dissolve(100);

var PAL_P = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
var PAL_DIV = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];
var VIS_P   = {min: 0.009, max: 0.92, palette: PAL_P};     // mapas de p (RF e GTB)
var VIS_DIF = {min: -0.5, max: 0.5, palette: PAL_DIV};     // diferença GTB − RF
var COR_LIMIAR = ['555555', '17becf', 'ff00ff'];           // só RF, só GTB, os dois

// ----------------------------------------------------------------------------
// Dados
// ----------------------------------------------------------------------------
function mapaP(modelo, alvo, mosaico) {
  return COL.filter(ee.Filter.eq('modelo', modelo))
            .filter(ee.Filter.eq('alvo', alvo))
            .filter(ee.Filter.eq('mosaico', mosaico))
            .mosaic().divide(10000).rename(modelo);
}

function pilha(alvo, mosaico) {
  return ee.Image.cat(MODELOS.map(function (m) { return mapaP(m, alvo, mosaico); }));
}

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var painel = ui.Panel({style: {width: '380px', padding: '6px'}});
painel.add(ui.Label('Favorabilidade RF × GTB', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('p = probability / 10000. Treino com fundo 3 × presença: p ordena a área, ' +
                    'não é a chance real de ocorrência.', {fontSize: '11px', color: '666666'}));

var selAlvo = ui.Select({items: ['fosfato'], value: 'fosfato', style: {stretch: 'horizontal'}});
var selMos = ui.Select({items: Object.keys(MOSAICOS), value: 'aster-mosaico-final', style: {stretch: 'horizontal'}});
painel.add(ui.Panel([ui.Label('Alvo', {fontSize: '12px'}), selAlvo], ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Panel([ui.Label('Mosaico', {fontSize: '12px'}), selMos], ui.Panel.Layout.flow('horizontal')));
var lblInfo = ui.Label('', {fontSize: '11px', whiteSpace: 'pre-wrap'});
painel.add(lblInfo);

// alvos e mosaicos que existem na coleção
ee.Dictionary({
  alvos: COL.aggregate_array('alvo').distinct(),
  mosaicos: COL.aggregate_array('mosaico').distinct(),
  n: COL.aggregate_histogram('modelo'),
  params: COL.first().toDictionary(['n_presenca', 'n_fundo', 'features'])
}).evaluate(function (r) {
  selAlvo.items().reset(r.alvos);
  selMos.items().reset(r.mosaicos);
  if (r.alvos.indexOf(selAlvo.getValue()) < 0) { selAlvo.setValue(r.alvos[0], false); }
  if (r.mosaicos.indexOf(selMos.getValue()) < 0) { selMos.setValue(r.mosaicos[0], false); }
  lblInfo.setValue('Células na coleção: ' + JSON.stringify(r.n) + '\nTreino: presença ' + r.params.n_presenca +
                   ', fundo ' + r.params.n_fundo + ' · ' + r.params.features.split(',').length + ' variáveis');
  desenhar();
});

var sldLimiar = ui.Slider({min: 0, max: 1, step: 0.01, value: 0.5, style: {stretch: 'horizontal'}});
painel.add(ui.Label('Limiar para "área favorável"', {fontSize: '12px', fontWeight: 'bold', margin: '8px 0 0 8px'}));
painel.add(sldLimiar);
var lblArea = ui.Label('', {fontSize: '11px', whiteSpace: 'pre-wrap'});
painel.add(lblArea);

var saida = ui.Panel();
painel.add(ui.Label('Análises', {fontWeight: 'bold', margin: '10px 0 2px 0'}));
painel.add(ui.Button({label: 'Histogramas + curva área × limiar', onClick: histogramas, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'Estatísticas (percentis)', onClick: estatisticas, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'p nas ocorrências × paisagem', onClick: pNosPontos, style: {stretch: 'horizontal'}}));
painel.add(saida);

// ----------------------------------------------------------------------------
// Mapa
// ----------------------------------------------------------------------------
var atual = null;
var camLimiar = null;

function desenhar() {
  var alvo = selAlvo.getValue(), mos = selMos.getValue();
  var img = pilha(alvo, mos);
  atual = {alvo: alvo, mosaico: mos, img: img};
  Map.layers().reset();
  camLimiar = null;
  saida.clear();
  Map.addLayer(ee.Image.constant(1), {min:0, max: 1}, 'base');
  // só o contorno externo da grid (células dissolvidas), preto, sem preenchimento
  Map.addLayer(ee.FeatureCollection([ee.Feature(regiao)])
                 .style({color: '000000', width: 2, fillColor: '00000000'}), {}, 'contorno da área');
  var base = ee.ImageCollection(MOSAICOS[mos]).mosaic();
  Map.addLayer(base, {bands: ['B04', 'B06', 'B08'], min: 55, max: 362}, mos + ' — SWIR', false);
  MODELOS.forEach(function (m, k) {
    Map.addLayer(img.select(m), VIS_P, alvo + ' — p ' + m, k === 1);
  });
  Map.addLayer(img.select('GTB').subtract(img.select('RF')), VIS_DIF, 'diferença GTB − RF', false);
  Map.addLayer(PONTOS.filter(ee.Filter.eq('alvo', alvo)).map(function (f) {
    return f.set('estilo', {color: '000000', pointShape: 'star5', pointSize: 7, width: 1,
                            fillColor: ee.Algorithms.If(ee.String(f.get('papel')).equals('treino'), 'ff0000', '00ffff')});
  }).style({styleProperty: 'estilo'}), {}, 'pontos ' + alvo + ' (vermelho treino, ciano validação)');
  Map.addLayer(ee.Image().paint(GRID, 0, 1), {palette: 'ffffff'}, 'grid', false);
  atualizarLimiar();
}

function atualizarLimiar() {
  if (!atual) { return; }
  var t = sldLimiar.getValue();
  if (camLimiar) { Map.layers().remove(camLimiar); }
  var acima = atual.img.gte(t);
  var conc = acima.select('RF').multiply(1).add(acima.select('GTB').multiply(2)).selfMask();   // 1 só RF, 2 só GTB, 3 os dois
  camLimiar = ui.Map.Layer(conc, {min: 1, max: 3, palette: COR_LIMIAR},
                           'p ≥ ' + t + ' (cinza só RF, ciano só GTB, magenta os dois)', true);
  Map.layers().add(camLimiar);
  lblLegLimiar.setValue('p ≥ ' + t);
  lblArea.setValue('calculando a área acima de ' + t + '…');
  acima.addBands(acima.select('RF').and(acima.select('GTB')).rename('ambos'))
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: regiao, scale: ESCALA_HIST, maxPixels: 1e9, bestEffort: true})
    .evaluate(function (r) {
      if (!r) { return; }
      lblArea.setValue('Área com p ≥ ' + t + ':  RF ' + (100 * r.RF).toFixed(1) + '%  ·  GTB ' +
                       (100 * r.GTB).toFixed(1) + '%  ·  os dois ' + (100 * r.ambos).toFixed(1) + '%');
    });
}
sldLimiar.onChange(atualizarLimiar);
selAlvo.onChange(desenhar);
selMos.onChange(desenhar);

// ----------------------------------------------------------------------------
// Histogramas e curva área × limiar
// ----------------------------------------------------------------------------
function histogramas() {
  if (!atual) { return; }
  saida.clear();
  saida.add(ui.Label('Histograma de p em toda a grid (' + ESCALA_HIST + ' m). Pico perto de 0 = maior parte ' +
                     'da área desfavorável; cauda à direita = áreas que o modelo destaca.',
                     {fontSize: '11px', whiteSpace: 'pre-wrap'}));
  MODELOS.forEach(function (m) {
    saida.add(ui.Chart.image.histogram({image: atual.img.select(m), region: regiao, scale: ESCALA_HIST,
                                        minBucketWidth: 1 / N_CLASSES, maxPixels: 1e9})
      .setSeriesNames([m])
      .setOptions({title: 'Histograma de p — ' + m, colors: [COR[m]], legend: {position: 'none'},
                   hAxis: {title: 'p', viewWindow: {min: 0, max: 1}}, vAxis: {title: 'pixels (' + ESCALA_HIST + ' m)'}}));
  });
  // os dois juntos, frequência em escala log (mostra a cauda)
  saida.add(ui.Chart.image.histogram({image: atual.img, region: regiao, scale: ESCALA_HIST,
                                      minBucketWidth: 1 / N_CLASSES, maxPixels: 1e9})
    .setOptions({title: 'RF × GTB (frequência em log)', colors: [COR.RF, COR.GTB],
                 hAxis: {title: 'p', viewWindow: {min: 0, max: 1}}, vAxis: {title: 'pixels', logScale: true}}));

  // curva área × limiar: fração da área com p ≥ t (1 − CDF)
  var hist = atual.img.reduceRegion({reducer: ee.Reducer.fixedHistogram(0, 1, N_CLASSES), geometry: regiao,
                                     scale: ESCALA_HIST, maxPixels: 1e9, bestEffort: true});
  hist.evaluate(function (h) {
    var linhas = [['limiar p', 'RF', 'GTB']];
    var acum = {};
    MODELOS.forEach(function (m) {
      var cont = h[m].map(function (x) { return x[1]; });
      var total = cont.reduce(function (a, b) { return a + b; }, 0);
      var resto = total;
      acum[m] = cont.map(function (c) { var v = resto / total; resto -= c; return v; });
    });
    for (var i = 0; i < N_CLASSES; i++) {
      linhas.push([h.RF[i][0], 100 * acum.RF[i], 100 * acum.GTB[i]]);
    }
    saida.add(ui.Chart(linhas, 'LineChart', {
      title: 'Área com p ≥ limiar (%)', colors: [COR.RF, COR.GTB], lineWidth: 2,
      hAxis: {title: 'limiar p'}, vAxis: {title: '% da área', viewWindow: {min: 0, max: 100}}
    }));
    saida.add(ui.Label('Curva que cai rápido = o modelo concentra p alto em pouca área (mapa seletivo); ' +
                       'curva que cai devagar = p alto espalhado (mapa pouco discriminante).',
                       {fontSize: '11px', whiteSpace: 'pre-wrap'}));
  });
}

// ----------------------------------------------------------------------------
// Estatísticas
// ----------------------------------------------------------------------------
function estatisticas() {
  if (!atual) { return; }
  saida.clear();
  saida.add(ui.Label('calculando…', {fontSize: '11px'}));
  var PCT = [1, 5, 25, 50, 75, 95, 99];
  var red = ee.Reducer.percentile(PCT).combine(ee.Reducer.mean(), '', true).combine(ee.Reducer.minMax(), '', true);
  atual.img.reduceRegion({reducer: red, geometry: regiao, scale: ESCALA_HIST, maxPixels: 1e9, bestEffort: true})
    .evaluate(function (r) {
      saida.clear();
      var nomes = ['min'].concat(PCT.map(function (p) { return 'p' + p; }), ['mean', 'max']);
      var rot = {min: 'mínimo', mean: 'média', max: 'máximo'};
      var dt = [['estatística', 'RF', 'GTB']];
      nomes.forEach(function (n) {
        dt.push([rot[n] || ('percentil ' + n.slice(1)), r['RF_' + n], r['GTB_' + n]]);
      });
      saida.add(ui.Label('p em toda a grid (' + ESCALA_HIST + ' m)', {fontSize: '11px', fontWeight: 'bold'}));
      saida.add(ui.Chart(dt, 'Table', {}));
      saida.add(ui.Label('Intervalo útil ≈ p5–p95. Se o p95 do RF fica muito abaixo do GTB, o RF "achata" ' +
                         'os valores (média de árvores) e o GTB espalha mais.', {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    });
}

// ----------------------------------------------------------------------------
// p nas ocorrências × paisagem
// ----------------------------------------------------------------------------
function pNosPontos() {
  if (!atual) { return; }
  saida.clear();
  var pts = PONTOS.filter(ee.Filter.eq('alvo', atual.alvo)).filterBounds(regiao);
  var nosPontos = atual.img.reduceRegions({collection: pts, reducer: ee.Reducer.first(), scale: 30});
  var paisagem = atual.img.sample({region: regiao, scale: ESCALA_HIST, numPixels: 3000, seed: 42});
  MODELOS.forEach(function (m) {
    var pp = nosPontos.filter(ee.Filter.notNull([m]));
    var pf = paisagem;
    saida.add(ui.Chart.feature.histogram({features: pp, property: m, minBucketWidth: 0.05})
      .setOptions({title: m + ' — p nas ocorrências', colors: ['#d62728'], legend: {position: 'none'},
                   hAxis: {title: 'p', viewWindow: {min: 0, max: 1}}}));
    saida.add(ui.Chart.feature.histogram({features: pf, property: m, minBucketWidth: 0.05})
      .setOptions({title: m + ' — p na paisagem (3 mil pixels)', colors: [COR[m]], legend: {position: 'none'},
                   hAxis: {title: 'p', viewWindow: {min: 0, max: 1}}}));
  });
  ee.Dictionary({
    rf: nosPontos.aggregate_array('RF'), gtb: nosPontos.aggregate_array('GTB'),
    papel: nosPontos.aggregate_array('papel')
  }).evaluate(function (r) {
    function med(v) { v = v.filter(function (x) { return x !== null; }).sort(function (a, b) { return a - b; });
                      return v.length ? v[Math.floor(v.length / 2)] : null; }
    saida.add(ui.Label('Mediana de p nas ' + r.rf.length + ' ocorrências: RF ' + (med(r.rf) || 0).toFixed(3) +
                       ' · GTB ' + (med(r.gtb) || 0).toFixed(3) +
                       '\nAtenção: o modelo foi treinado com esses pontos — p alto neles é esperado e não ' +
                       'mede generalização (para isso: avaliar_maxent.ipynb, célula de fora).',
                       {fontSize: '11px', whiteSpace: 'pre-wrap'}));
  });
}

// ----------------------------------------------------------------------------
// Legenda (canto inferior esquerdo; para os prints do relatório)
// ----------------------------------------------------------------------------
var FS_LEG = '11px';

function barraCores(vis) {
  // gradiente horizontal com a mesma paleta da camada
  return ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '150x10', format: 'png',
             min: 0, max: 1, palette: vis.palette},
    style: {stretch: 'horizontal', margin: '0 8px', maxHeight: '12px'}
  });
}

function linhaEscala(vis) {
  return ui.Panel([
    ui.Label(String(vis.min), {fontSize: '10px', margin: '0 8px'}),
    ui.Label('', {stretch: 'horizontal'}),
    ui.Label(String(vis.max), {fontSize: '10px', margin: '0 8px'})
  ], ui.Panel.Layout.flow('horizontal'));
}

function itemCor(cor, texto, forma) {
  var simbolo = forma === 'estrela'
    ? ui.Label('★', {color: '#' + cor, fontSize: '15px', margin: '0 6px 0 4px', padding: '0'})
    : ui.Label('', {backgroundColor: '#' + cor, padding: '7px', margin: '2px 6px 2px 4px',
                    border: '1px solid #888888'});
  return ui.Panel([simbolo, ui.Label(texto, {fontSize: FS_LEG, margin: '3px 0'})],
                  ui.Panel.Layout.flow('horizontal'));
}

function titulo(t) { return ui.Label(t, {fontSize: FS_LEG, fontWeight: 'bold', margin: '6px 0 2px 4px'}); }

var lblLegLimiar = ui.Label('p ≥ —', {fontSize: FS_LEG, fontWeight: 'bold', margin: '6px 0 2px 4px'});
var legenda = ui.Panel({style: {position: 'bottom-left', padding: '6px 8px', width: '200px'}});
legenda.add(ui.Label('Favorabilidade', {fontSize: '13px', fontWeight: 'bold', margin: '0 0 2px 4px'}));
legenda.add(titulo('p (RF / GTB)'));
legenda.add(barraCores(VIS_P));
legenda.add(linhaEscala(VIS_P));
legenda.add(titulo('Diferença GTB − RF'));
legenda.add(barraCores(VIS_DIF));
legenda.add(linhaEscala(VIS_DIF));
legenda.add(lblLegLimiar);
legenda.add(itemCor(COR_LIMIAR[0], 'só RF'));
legenda.add(itemCor(COR_LIMIAR[1], 'só GTB'));
legenda.add(itemCor(COR_LIMIAR[2], 'RF e GTB'));
legenda.add(titulo('Ocorrências'));
legenda.add(itemCor('ff0000', 'treino (bacia)', 'estrela'));
legenda.add(itemCor('00ffff', 'validação', 'estrela'));
Map.add(legenda);

var chkLegenda = ui.Checkbox('Mostrar legenda no mapa', true);
chkLegenda.onChange(function (v) { legenda.style().set('shown', v); });
painel.add(chkLegenda);

// ----------------------------------------------------------------------------
// Clique
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '240px'}});
info.add(ui.Label('Clique no mapa para ver p.', {fontSize: '11px'}));
Map.add(info);
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  if (!atual) { return; }
  info.clear();
  info.add(ui.Label('carregando…', {fontSize: '11px'}));
  atual.img.reduceRegion(ee.Reducer.first(), ee.Geometry.Point(c.lon, c.lat), 30).evaluate(function (r) {
    info.clear();
    r = r || {};
    function f(x) { return (x === null || x === undefined) ? '—' : x.toFixed(3); }
    info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontSize: '11px', fontWeight: 'bold'}));
    info.add(ui.Label('p RF ' + f(r.RF) + '  ·  p GTB ' + f(r.GTB), {fontSize: '12px'}));
  });
});

ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);
