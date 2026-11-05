// ============================================================================
// Mosaico sazonal seco ASTER: ORIGINAL × HARMONIZADO com o MODIS — visualizador
//
// Fontes (projects/mapbiomas-arida/mine)
//   ASTER_sazonal              original (build_mosaico_sazonal_seco.py), TOA × 1000
//   ASTER_sazonal_harmonizado  harmonizar_sazonal_modis.py etapa 4 (tasks harm_<P>_<LxCy>):
//                              mesmas bandas em escala NBAR (superfície × 1000), TIR igual,
//                              + 'cena_t' (s desde 2000-01-01 da cena do pixel)
//                              + 'harmonizado' (1 = reta da cena, 2 = mediana do período)
//
// Tela dividida com os mapas ligados: esquerda = original, direita = harmonizado.
// Arraste a divisória para comparar. As faixas de cor da direita são as da esquerda
// passadas pela reta mediana (ganho × valor + offset), então uma cena "média" fica
// com a mesma cor nos dois lados e o que muda é a emenda entre cenas.
//
// O que olhar
//   1. VNIR / NDVI: o degrau entre cenas sumiu ou diminuiu?
//   2. Bordas (gradiente B3N): emenda = linha reta longa; no harmonizado deve enfraquecer.
//   3. Ganho aplicado (harm / orig no B3N): mostra a colcha das cenas e quanto cada uma
//      foi corrigida. Cena = cor aleatória por cena_t.
//   4. harmonizado = 2: onde entrou a mediana (cena sem reta válida ou confere = 0).
// Clique no mapa: valores original e harmonizado, cena e data do pixel.
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
var COL_ORIG = ee.ImageCollection(MINE + 'ASTER_sazonal');
// 2ª passada (harmonizar_sazonal_relativo.py): células ajustadas contra as vizinhas
// substituem as mesmas células do harmonizado. USAR_REL = false volta à 1ª passada.
var USAR_REL = true;
var COL_REL = ee.ImageCollection(MINE + 'ASTER_sazonal_harm_rel');
var COL_HARM1 = ee.ImageCollection(MINE + 'ASTER_sazonal_harmonizado');
// mosaic() usa só as bandas da 1ª imagem: as células da 1ª passada ganham um
// 'ajuste_rel' todo mascarado para que a banda das células da 2ª passada não se perca.
var COL_MODIS = USAR_REL
  ? COL_HARM1.filter(ee.Filter.inList('system:index', COL_REL.aggregate_array('system:index')).not())
             .map(function (img) {
               return img.addBands(ee.Image(0).toInt32().selfMask().rename('ajuste_rel')); })
             .merge(COL_REL)
  : COL_HARM1;
// Lado direito: 'modis' (1ª + 2ª passada) ou o casamento de histograma contra o Landsat
// (harmonizar_sazonal_landsat.py): 'landsat_razao' / 'landsat_residuo'. HARM é só a
// escolha inicial: o seletor "Harmonização" do painel troca sem recarregar o script.
var HARM = 'landsat_razao';
var COLS_HARM = {
  modis: COL_MODIS,
  landsat_razao: ee.ImageCollection(MINE + 'ASTER_sazonal_hm_razao'),
  landsat_residuo: ee.ImageCollection(MINE + 'ASTER_sazonal_hm_residuo')
};
var NOMES_HARM = {
  landsat_razao: 'Landsat — razão (B05–B09 e TIR × guia\'/guia)',
  landsat_residuo: 'Landsat — resíduo (relação de fundo por pedaço)',
  modis: 'MODIS MCD43A4 (1ª + 2ª passada)'
};
var ROTULO_HARM = {landsat_razao: 'HARMONIZADO — Landsat razão',
                   landsat_residuo: 'HARMONIZADO — Landsat resíduo',
                   modis: 'HARMONIZADO — MODIS (NBAR)'};
// Landsat de referência: asset salvo (exportar_landsat_sazonal.py) ou calculado na hora.
var LS_ASSET = true;
var grid = ee.FeatureCollection(MINE + 'grid_irece_5x5_solape100');
var T0 = ee.Date('2000-01-01T00:00:00Z');

// faixas por banda no original (TOA × 1000; TIR em DN)
var FAIXA = {
  B01: [69, 176], B02: [51, 204], B3N: [155, 332], B04: [167, 362], B05: [77, 232],
  B06: [75, 234], B07: [71, 228], B08: [45, 168], B09: [34, 115],
  B10: [1124, 1635], B11: [1228, 1803], B12: [1332, 1931], B13: [1580, 2203], B14: [1644, 2259]
};
// faixa do lado harmonizado = faixa original × ganho + offset, para as duas telas ficarem
// comparáveis. [ganho, offset × 1000]
// MODIS: retas medianas (coef_modis_<P>.csv, 2026-10-31)
var RETA_MODIS = {
  A: {B01: [1.164, -51.1], B02: [1.124, -35.0], B3N: [1.277, -14.3], B04: [1.199, -19.5], SW2: [1.177, 0]},
  B: {B01: [1.189, -49.7], B02: [1.163, -31.9], B3N: [1.267, -15.5]}
};
// Landsat: ganho mediano no percentil 50 (etapa 4 de harmonizar_sazonal_landsat.py,
// 2026-10-31); B05–B09 seguem o B04. TIR: B13 ≈ 1, sem correção.
var RETA_LANDSAT = {
  A: {B01: [0.694, 0], B02: [0.795, 0], B3N: [1.036, 0], B04: [0.935, 0], SW2: [0.935, 0]},
  B: {B01: [0.688, 0], B02: [0.847, 0], B3N: [1.087, 0]}
};
var SW2 = ['B05', 'B06', 'B07', 'B08', 'B09'];

// Referência Landsat (só para olhar): composições de 32 dias, reflectância de superfície,
// mesma janela e meses secos do mosaico ASTER (build_mosaico_sazonal_seco.py), mediana.
// Bandas renomeadas para as equivalentes do ASTER (green→B01, red→B02, nir→B3N,
// swir1→B04, swir2→SW2) e × 1000, na escala do harmonizado (NBAR × 1000).
var LS_DATAS = {A: ['2000-03-01', '2008-04-01'], B: ['2009-01-01', '2021-01-01']};
var MESES_SECO = [6, 7, 8, 9, 10];

function landsat(periodo) {
  if (LS_ASSET) {
    var bl = periodo === 'A' ? ['green', 'red', 'nir', 'swir1', 'swir2'] : ['green', 'red', 'nir'];
    var ba = periodo === 'A' ? ['B01', 'B02', 'B3N', 'B04', 'SW2'] : ['B01', 'B02', 'B3N'];
    // um asset por célula (exportar_landsat_sazonal.py): junta as células do período
    return ee.ImageCollection(MINE + 'Landsat_sazonal')
             .filter(ee.Filter.eq('periodo', periodo))
             .mosaic().select(bl, ba);
  }
  var d = LS_DATAS[periodo];
  return ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_32DAY')
    .filterDate(d[0], d[1])
    .filter(ee.Filter.calendarRange(MESES_SECO[0], MESES_SECO[MESES_SECO.length - 1], 'month'))
    .select(['green', 'red', 'nir', 'swir1', 'swir2'], ['B01', 'B02', 'B3N', 'B04', 'SW2'])
    .median()
    .multiply(1000);
}

function faixa(b, periodo, harm) {
  var f = FAIXA[b];
  if (!harm) return f;
  var r = (estado.harm === 'modis' ? RETA_MODIS : RETA_LANDSAT)[periodo][SW2.indexOf(b) >= 0 ? 'SW2' : b];
  if (!r) return f;                                   // TIR: sem correção
  return [f[0] * r[0] + r[1], f[1] * r[0] + r[1]];
}

function rgb(img, bandas, periodo, harm) {
  return [img.select(bandas), {
    bands: bandas,
    min: bandas.map(function (b) { return faixa(b, periodo, harm)[0]; }),
    max: bandas.map(function (b) { return faixa(b, periodo, harm)[1]; })}];
}

var PAL_VEG = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];
var PAL_DIV = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];
var PAL_GRAD = ['000000', '440154', '3b528b', '21918c', '5ec962', 'fde725'];

function ndvi(img) { return img.normalizedDifference(['B3N', 'B02']); }
// EVI2 (Jiang et al. 2008): EVI sem a banda azul (o ASTER não tem); entrada × 1000
function evi2(img) {
  var n = img.select('B3N').divide(1000), r = img.select('B02').divide(1000);
  return n.subtract(r).multiply(2.5).divide(n.add(r.multiply(2.4)).add(1));
}

// módulo do gradiente do B3N, relativo ao valor (adimensional): emendas viram linhas
function bordas(img) {
  var b = img.select('B3N').toFloat();
  var g = b.gradient();
  return g.select('x').hypot(g.select('y')).divide(b.max(1)).multiply(1e3);
}

// ----------------------------------------------------------------------------
// Camadas: cada uma devolve [imagem, vis] para um lado (harm = false/true)
// ----------------------------------------------------------------------------
var CAMADAS = {
  'VNIR (3N-2-1)':  {f: function (img, p, h) { return rgb(img, ['B3N', 'B02', 'B01'], p, h); }},
  'SWIR (4-6-8) — só A': {soA: true, f: function (img, p, h) { return rgb(img, ['B04', 'B06', 'B08'], p, h); }},
  'SWIR (4-7-9) — só A': {soA: true, f: function (img, p, h) { return rgb(img, ['B04', 'B07', 'B09'], p, h); }},
  'TIR (14-12-10)': {f: function (img, p, h) { return rgb(img, ['B14', 'B12', 'B10'], p, h); }},
  'NDVI':           {f: function (img) { return [ndvi(img), {min: 0.1, max: 0.7, palette: PAL_VEG}]; }},
  'EVI2':           {f: function (img) { return [evi2(img), {min: 0.05, max: 0.45, palette: PAL_VEG}]; }},
  'Bordas (gradiente B3N)': {f: function (img) { return [bordas(img), {min: 0, max: 15, palette: PAL_GRAD}]; }},
  'Al-OH (B05+B07)/B06 — só A': {soA: true, f: function (img) {
    return [img.select('B05').add(img.select('B07')).divide(img.select('B06').max(1)),
            {min: 1.9, max: 2.04, palette: ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf']}]; }},
  'carbonato (B06+B09)/B08 — só A': {soA: true, f: function (img) {
    return [img.select('B06').add(img.select('B09')).divide(img.select('B08').max(1)),
            {min: 2.0, max: 2.43, palette: ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf']}]; }}
};
// só no lado harmonizado (usam as duas imagens)
var EXTRAS = {
  'Ganho aplicado B3N (harm / orig)': function (orig, harm) {
    return [harm.select('B3N').divide(orig.select('B3N').max(1)), {min: 1.0, max: 1.5, palette: PAL_DIV}]; },
  'Cena do pixel (cor por cena_t)': function (orig, harm) {
    return [harm.select('cena_t').randomVisualizer(), {}]; },
  'harmonizado = 2 (mediana)': function (orig, harm) {
    return [harm.select('harmonizado').eq(2).selfMask(), {palette: ['ff00ff']}]; },
  'hm = 2 (curva da célula, Landsat)': function (orig, harm) {
    return [harm.select('hm').eq(2).selfMask(), {palette: ['ff00ff']}]; },
  'ajuste_rel (1 = cena, 2 = mediana)': function (orig, harm) {
    return [harm.select('ajuste_rel'), {min: 1, max: 2, palette: ['00ffff', 'ff00ff']}]; },
  'Diferença NDVI (harm − orig)': function (orig, harm) {
    return [ndvi(harm).subtract(ndvi(orig)), {min: -0.1, max: 0.1, palette: PAL_DIV}]; }
};

// ----------------------------------------------------------------------------
// Mapas ligados
// ----------------------------------------------------------------------------
var mapaOrig = ui.Map();
var mapaHarm = ui.Map();
ui.Map.Linker([mapaOrig, mapaHarm]);
[mapaOrig, mapaHarm].forEach(function (m) { m.setOptions('SATELLITE'); m.setControlVisibility({layerList: true}); });
mapaOrig.add(ui.Label('ORIGINAL (TOA)', {fontWeight: 'bold', position: 'bottom-left'}));
var lblHarm = ui.Label(ROTULO_HARM[HARM], {fontWeight: 'bold', position: 'bottom-right'});
mapaHarm.add(lblHarm);

var estado = {periodo: 'A', camada: 'VNIR (3N-2-1)', extra: '(nenhuma)', landsat: false,
              harm: HARM};

function imagens() {
  var p = estado.periodo;
  var harm = COLS_HARM[estado.harm].filter(ee.Filter.eq('periodo', p));
  // original só nas células que já têm harmonizado (compara a mesma área). Filtra por
  // celula_nome: merge() troca o system:index por '1_…' / '2_…' e não casaria mais.
  var orig = COL_ORIG.filter(ee.Filter.eq('periodo', p))
                     .filter(ee.Filter.inList('celula_nome', harm.aggregate_array('celula_nome')));
  return {orig: orig.mosaic(), harm: harm.mosaic(), colHarm: harm};
}

function contorno(colHarm) {
  var nomes = colHarm.aggregate_array('celula_nome');
  var com = grid.filter(ee.Filter.inList('nome', nomes));
  return {todas: grid.style({color: '888888', width: 1, fillColor: '00000000'}),
          com: com.style({color: 'ffff00', width: 2, fillColor: '00000000'})};
}

function atualizar() {
  var im = imagens();
  var c = CAMADAS[estado.camada];
  mapaOrig.layers().reset();
  mapaHarm.layers().reset();
  if (c.soA && estado.periodo !== 'A') {
    lblInfo.setValue('Camada só existe no período A (B não tem SWIR).');
  } else {
    var o = c.f(im.orig, estado.periodo, false);
    var h = c.f(im.harm, estado.periodo, true);
    mapaOrig.addLayer(o[0], o[1], estado.periodo + ' original — ' + estado.camada);
    mapaHarm.addLayer(h[0], h[1], estado.periodo + ' ' + estado.harm + ' — ' + estado.camada);
  }
  if (estado.extra !== '(nenhuma)') {
    var e = EXTRAS[estado.extra](im.orig, im.harm);
    mapaHarm.addLayer(e[0], e[1], estado.extra);
  }
  // Landsat no lado esquerdo (por cima do original): compara com o harmonizado da direita,
  // os dois em reflectância de superfície. Mostra NDVI se a camada escolhida for NDVI.
  if (estado.landsat) {
    var ls = landsat(estado.periodo);
    var vn = rgb(ls, ['B3N', 'B02', 'B01'], estado.periodo, true);
    var verNdvi = estado.camada === 'NDVI', verEvi = estado.camada === 'EVI2';
    mapaOrig.addLayer(vn[0], vn[1], 'Landsat 32 dias ' + estado.periodo + ' — VNIR (nir-red-green)',
                      !verNdvi && !verEvi);
    mapaOrig.addLayer(ndvi(ls), {min: 0.1, max: 0.7, palette: PAL_VEG},
                      'Landsat 32 dias ' + estado.periodo + ' — NDVI', verNdvi);
    mapaOrig.addLayer(evi2(ls), {min: 0.05, max: 0.45, palette: PAL_VEG},
                      'Landsat 32 dias ' + estado.periodo + ' — EVI2', verEvi);
  }
  var g = contorno(im.colHarm);
  [mapaOrig, mapaHarm].forEach(function (m) {
    m.addLayer(g.todas, {}, 'grid', true, 0.6);
    m.addLayer(g.com, {}, 'células harmonizadas');
  });
  im.colHarm.aggregate_array('celula_nome').evaluate(function (n, erro) {
    if (erro) { lblCel.setValue('Erro ao ler a coleção harmonizada (' + estado.harm + '): ' + erro); return; }
    lblCel.setValue(NOMES_HARM[estado.harm] + '\nPeríodo ' + estado.periodo + ': ' + n.length +
                    ' células harmonizadas' +
                    (n.length ? ' — ' + n.sort().join(', ') : ''));
  });
}

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var lblCel = ui.Label('', {fontSize: '11px', whiteSpace: 'pre-wrap'});
var lblInfo = ui.Label('Clique no mapa para ver os valores.', {fontSize: '11px', whiteSpace: 'pre'});

var selPeriodo = ui.Select({items: ['A', 'B'], value: 'A', onChange: function (v) {
  estado.periodo = v; atualizar(); }});
var selCamada = ui.Select({items: Object.keys(CAMADAS), value: estado.camada, onChange: function (v) {
  estado.camada = v; atualizar(); }});
var selExtra = ui.Select({items: ['(nenhuma)'].concat(Object.keys(EXTRAS)), value: '(nenhuma)',
  onChange: function (v) { estado.extra = v; atualizar(); }});
var selHarm = ui.Select({
  items: Object.keys(NOMES_HARM).map(function (k) { return {label: NOMES_HARM[k], value: k}; }),
  value: HARM, onChange: function (v) {
    estado.harm = v; lblHarm.setValue(ROTULO_HARM[v]);
    // camadas que só existem numa fonte: volta para '(nenhuma)'
    var soModis = ['harmonizado = 2 (mediana)', 'ajuste_rel (1 = cena, 2 = mediana)'];
    var soLandsat = ['hm = 2 (curva da célula, Landsat)'];
    if ((v === 'modis' ? soLandsat : soModis).indexOf(estado.extra) >= 0) {
      estado.extra = '(nenhuma)';
      selExtra.setValue('(nenhuma)', false);
    }
    atualizar(); }});
var chkLandsat = ui.Checkbox('Landsat 32 dias à esquerda (mesmos meses secos)', false,
  function (v) { estado.landsat = v; atualizar(); });

var btnCentro = ui.Button('Centralizar na área harmonizada', function () {
  var col = imagens().colHarm;
  var geom = grid.filter(ee.Filter.inList('nome', col.aggregate_array('celula_nome'))).geometry();
  mapaOrig.centerObject(geom);
});

var painel = ui.Panel({style: {width: '330px', position: 'top-left'}, widgets: [
  ui.Label('Original × harmonizado', {fontWeight: 'bold', fontSize: '15px'}),
  ui.Label('Harmonização (lado direito)'), selHarm,
  ui.Label('Período'), selPeriodo,
  ui.Label('Camada (dos dois lados)'), selCamada,
  ui.Label('Camada extra (só harmonizado)'), selExtra,
  chkLandsat,
  btnCentro, lblCel, lblInfo
]});
mapaOrig.add(painel);

// ----------------------------------------------------------------------------
// Clique: original × harmonizado no ponto
// ----------------------------------------------------------------------------
function aoClicar(c) {
  var pt = ee.Geometry.Point([c.lon, c.lat]);
  var im = imagens();
  lblInfo.setValue('Lendo…');
  var leitura = ee.Dictionary({
    orig: im.orig.reduceRegion(ee.Reducer.first(), pt, 15),
    harm: im.harm.reduceRegion(ee.Reducer.first(), pt, 15),
    ls: estado.landsat
      ? landsat(estado.periodo).reduceRegion(ee.Reducer.first(), pt, 30) : ee.Dictionary({})
  });
  leitura.evaluate(function (r, erro) {
    if (erro) { lblInfo.setValue('Erro: ' + erro); return; }
    var o = r.orig, h = r.harm;
    if (!h || h.B01 === null || h.B01 === undefined) {
      lblInfo.setValue('Sem harmonizado neste ponto.'); return; }
    var bandas = ['B01', 'B02', 'B3N'];
    if (estado.periodo === 'A') bandas = bandas.concat(['B04', 'B05', 'B06', 'B07', 'B08', 'B09']);
    var txt = 'banda   orig   harm   harm/orig\n';
    bandas.forEach(function (b) {
      txt += (b + '     ').slice(0, 6) + ('     ' + o[b]).slice(-6) + ' ' + ('     ' + h[b]).slice(-6) +
             '   ' + (o[b] ? (h[b] / o[b]).toFixed(3) : '-') + '\n';
    });
    if (h.harmonizado !== null && h.harmonizado !== undefined) {
      txt += 'harmonizado = ' + h.harmonizado + (h.harmonizado === 2 ? ' (mediana)' : ' (reta da cena)') + '\n';
    }
    if (h.hm !== null && h.hm !== undefined) {
      txt += 'hm = ' + h.hm + (h.hm === 2 ? ' (curva da célula)' : ' (curva do pedaço)') + '\n';
    }
    if (h.ajuste_rel !== null && h.ajuste_rel !== undefined) {
      txt += 'ajuste_rel = ' + h.ajuste_rel + (h.ajuste_rel === 2 ? ' (mediana)' : ' (reta da cena)') + '\n';
    }
    if (estado.landsat && r.ls && r.ls.B3N !== null && r.ls.B3N !== undefined) {
      var l = r.ls;
      txt += 'Landsat ×1000: B01 ' + l.B01.toFixed(0) + '  B02 ' + l.B02.toFixed(0) +
             '  B3N ' + l.B3N.toFixed(0) + '  NDVI ' + ((l.B3N - l.B02) / (l.B3N + l.B02)).toFixed(3) + '\n';
    }
    txt += 'n_obs = ' + o.n_obs + '\ncena_t = ' + h.cena_t;
    lblInfo.setValue(txt + '\ndata: …');
    T0.advance(h.cena_t, 'second').format('YYYY-MM-dd HH:mm').evaluate(function (d) {
      lblInfo.setValue(txt + '\ndata: ' + d);
    });
  });
}
mapaOrig.onClick(aoClicar);
mapaHarm.onClick(aoClicar);

ui.root.widgets().reset([ui.SplitPanel({firstPanel: mapaOrig, secondPanel: mapaHarm, wipe: true})]);
mapaOrig.setCenter(-41.8, -12.1, 8);
atualizar();
