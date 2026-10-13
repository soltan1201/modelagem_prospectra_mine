// ============================================================================
// Mosaico final ASTER + índices espectrais + Acervo CPRM — visualizador
//
// Nada é desenhado ao abrir: escolha no painel à esquerda o que ver. Cada
// caixa liga/desliga uma camada (a camada é criada na primeira vez que é
// marcada, então o mapa só calcula o que você pedir).
//
// Fontes
//   Mosaico final : projects/mapbiomas-arida/mine/aster-mosaico-final (25 células, mosaic())
//                   B01–B09 reflectância TOA × 1000; B10–B14 TIR em DN;
//                   'origem' 1–4; 'swir_estimado' 0/1  (build_mosaico_final.py)
//   Índices       : mesmas fórmulas de modeling/maxent_fosfato.py (INDICES +
//                   phosphorite_rbd pela emissividade TIR)
//   Acervo CPRM   : projects/mapbiomas-arida/mine/acervo_cprm (ver ver_acervo_cprm.js)
//
// Faixas de cor: percentis 2–98 % na área de estudo (500 m), medidos em 2026-10-12.
// Clique no mapa: bandas, índices, gama e unidade geológica no ponto.
// ============================================================================

var MOSAICO = ee.ImageCollection('projects/mapbiomas-arida/mine/aster-mosaico-final').mosaic();
var ACERVO  = 'projects/mapbiomas-arida/mine/acervo_cprm/';
var ESCALA_TOA = 1000;

var area = ee.Geometry.Rectangle([-44.460483519818084, -14.790793773956535,
                                  -39.203525512005584, -9.416525285369607]);
var grid = ee.FeatureCollection('projects/mapbiomas-arida/mine/grid_irece_5x5_solape100');

var BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
              'B10', 'B11', 'B12', 'B13', 'B14'];
// faixas por banda (VNIR/SWIR × 1000; TIR em DN)
var FAIXA_BANDA = {
  B01: [69, 176], B02: [51, 204], B3N: [155, 332], B04: [167, 362], B05: [77, 232],
  B06: [75, 234], B07: [71, 228], B08: [45, 168], B09: [34, 115],
  B10: [1124, 1635], B11: [1228, 1803], B12: [1332, 1931], B13: [1580, 2203], B14: [1644, 2259]
};

// ----------------------------------------------------------------------------
// Índices (iguais aos de maxent_fosfato.py)
// ----------------------------------------------------------------------------
function razao(a, b) { return MOSAICO.select(a).divide(MOSAICO.select(b)); }

var PAL_SEQ = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
var PAL_VEG = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];

var INDICES = {
  NDVI:                 {f: function () { return MOSAICO.normalizedDifference(['B3N', 'B02']); },
                         min: 0.12, max: 0.68, pal: PAL_VEG, grupo: 'Vegetação'},
  SAVI:                 {f: function () {
                           return MOSAICO.expression('1.5 * (n - r) / (n + r + 0.5)', {
                             n: MOSAICO.select('B3N').divide(ESCALA_TOA), r: MOSAICO.select('B02').divide(ESCALA_TOA)}); },
                         min: 0.08, max: 0.41, pal: PAL_VEG, grupo: 'Vegetação'},
  ferric_iron_2_1:      {f: function () { return razao('B02', 'B01'); }, min: 0.70, max: 1.32, grupo: 'Ferro'},
  ferric_ox_4_3:        {f: function () { return razao('B04', 'B3N'); }, min: 0.61, max: 1.67, grupo: 'Ferro'},
  gossan_4_2:           {f: function () { return razao('B04', 'B02'); }, min: 1.39, max: 3.48, grupo: 'Ferro'},
  ferrous_sil_5_4:      {f: function () { return razao('B05', 'B04'); }, min: 0.45, max: 0.67, grupo: 'Ferro'},
  sericite_illite_smectite: {f: function () { return MOSAICO.select('B05').add(MOSAICO.select('B07')).divide(MOSAICO.select('B06')); },
                         min: 1.88, max: 2.04, grupo: 'Argilas / Al-OH'},
  alunite_kaolinite_pyrophyllite: {f: function () { return MOSAICO.select('B04').add(MOSAICO.select('B06')).divide(MOSAICO.select('B05')); },
                         min: 2.48, max: 3.23, grupo: 'Argilas / Al-OH'},
  muscovite_7_6:        {f: function () { return razao('B07', 'B06'); }, min: 0.90, max: 1.01, grupo: 'Argilas / Al-OH'},
  kaolinite_7_5:        {f: function () { return razao('B07', 'B05'); }, min: 0.88, max: 1.01, grupo: 'Argilas / Al-OH'},
  dolomite:             {f: function () { return MOSAICO.select('B06').add(MOSAICO.select('B08')).divide(MOSAICO.select('B07')); },
                         min: 1.65, max: 1.82, grupo: 'Carbonatos'},
  carbonate_8_6:        {f: function () { return razao('B08', 'B06'); }, min: 0.58, max: 0.74, grupo: 'Carbonatos'},
  carbonate_13_14:      {f: function () { return razao('B13', 'B14'); }, min: 0.950, max: 0.983, grupo: 'Carbonatos'},
  quartz_14_12:         {f: function () { return razao('B14', 'B12'); }, min: 1.155, max: 1.249, grupo: 'Sílica / TIR'},
  phosphorite_rbd:      {f: function () { return fosforitoTIR(); }, min: 0.986, max: 1.020, grupo: 'Sílica / TIR'}
};

// Emissividade TIR (normalização de emissividade) → (ε11 + ε13) / (2 ε12)
function fosforitoTIR() {
  var UCC = {B10: 6.822e-3, B11: 6.780e-3, B12: 6.590e-3, B13: 5.693e-3, B14: 5.225e-3};
  var LAM = {B10: 8.291, B11: 8.634, B12: 9.075, B13: 10.657, B14: 11.318};
  var C1 = 1.191042e8, C2 = 1.4387752e4, EPS = 0.96;
  var L = {}, T = [];
  Object.keys(UCC).forEach(function (b) {
    L[b] = MOSAICO.select(b).subtract(1).max(1).multiply(UCC[b]);
    T.push(L[b].expression('c2 / (lam * log(eps * c1 / (pow(lam, 5) * L) + 1))',
                           {c1: C1, c2: C2, lam: LAM[b], eps: EPS, L: L[b]}));
  });
  var tmax = ee.Image.cat(T).reduce(ee.Reducer.max());
  function e(b) {
    return L[b].divide(tmax.expression('c1 / (pow(lam, 5) * (exp(c2 / (lam * t)) - 1))',
                                       {c1: C1, c2: C2, lam: LAM[b], t: tmax}));
  }
  return e('B11').add(e('B13')).divide(e('B12').multiply(2));
}

function indice(nome) { return INDICES[nome].f().toFloat().rename(nome); }

// ----------------------------------------------------------------------------
// Acervo CPRM
// ----------------------------------------------------------------------------
var gama    = ee.Image(ACERVO + '1042_gama');
var gamaPos = gama.updateMask(gama.select('eU_ppm').gte(0));
var litoJac = ee.FeatureCollection(ACERVO + 'lito_jacobina_250k');
var litoBa  = ee.FeatureCollection(ACERVO + 'lito_bahia');
var fosfato = ee.FeatureCollection(ACERVO + 'ocorrencias_fosfato');
var geoq    = ee.FeatureCollection(ACERVO + 'geoquimica_p2o5');

var ALVO = ['Salitre', 'Irecê', 'Jussara', 'Nova América', 'Gabriel', 'Lapão', 'Bebedouro', 'Una'];
var filtroAlvo = ee.Filter.or.apply(null, ALVO.map(function (t) {
  return ee.Filter.stringContains('NOME_UNIDA', t);
}));

var PAL_RAZAO = ['2c7bb6', 'abd9e9', 'ffffbf', 'fdae61', 'd7191c'];
var PAL_MAG = ['053061', '2166ac', '92c5de', 'f7f7f7', 'f4a582', 'b2182b', '67001f'];
var PAL_UNID = ['8dd3c7', 'ffffb3', 'bebada', 'fb8072', '80b1d3', 'fdb462', 'b3de69',
                'fccde5', 'd9d9d9', 'bc80bd', 'ccebc5', 'ffed6f', 'a6cee3', 'b2df8a'];

function corPorUnidade(fc) {
  var unidades = fc.aggregate_array('NOME_UNIDA').distinct();
  var pal = ee.List(PAL_UNID);
  return fc.map(function (f) {
    var cor = ee.String(pal.get(ee.Number(unidades.indexOf(f.get('NOME_UNIDA'))).mod(pal.length())));
    return f.set('estilo', {color: '555555', width: 0.4, fillColor: cor.cat('b0')});
  }).style({styleProperty: 'estilo'});
}

var CORES_GRUPO = {bacia_irece: 'ff0000', leste_outros: 'ff9900', norte: 'ffff00'};

// ----------------------------------------------------------------------------
// Catálogo de camadas: grupo → [{nome, criar: function → [objeto, vis]}]
// ----------------------------------------------------------------------------
function rgb(b, faixas) {
  return function () {
    return [MOSAICO.select(b), {bands: b, min: b.map(function (x) { return faixas[x][0]; }),
                                max: b.map(function (x) { return faixas[x][1]; })}];
  };
}

var CATALOGO = [
  {grupo: 'Mosaico final — composições', itens: [
    {nome: 'Cor falsa VNIR (3N-2-1)',        criar: rgb(['B3N', 'B02', 'B01'], FAIXA_BANDA)},
    {nome: 'SWIR (4-6-8)',                   criar: rgb(['B04', 'B06', 'B08'], FAIXA_BANDA)},
    {nome: 'SWIR-VNIR (4-3N-2)',             criar: rgb(['B04', 'B3N', 'B02'], FAIXA_BANDA)},
    {nome: 'Al-OH / carbonato (5-7-8)',      criar: rgb(['B05', 'B07', 'B08'], FAIXA_BANDA)},
    {nome: 'TIR (13-12-10)',                 criar: rgb(['B13', 'B12', 'B10'], FAIXA_BANDA)},
    {nome: 'Origem do pixel (1–4)', criar: function () {
      return [MOSAICO.select('origem'), {min: 1, max: 4, palette: ['1a9850', '91cf60', 'fee08b', 'd73027']}]; }},
    {nome: 'SWIR estimado (origens 2 e 4)', criar: function () {
      return [MOSAICO.select('swir_estimado').selfMask(), {palette: ['ff00ff']}]; }}
  ]},
  {grupo: 'Composições de índices', itens: [
    {nome: 'Sultan RGB (4/7, 4/1, 2/3·4/3)', criar: function () {
      return [ee.Image.cat([indiceSultanR(), razao('B04', 'B01'),
                            razao('B02', 'B3N').multiply(razao('B04', 'B3N'))]),
              {min: [1.51, 1.51, 0.11], max: [2.45, 2.98, 1.20]}]; }},
    {nome: 'Abrams RGB (4/7, 4/3, 2/1)', criar: function () {
      return [ee.Image.cat([indiceSultanR(), razao('B04', 'B3N'), razao('B02', 'B01')]),
              {min: [1.51, 0.61, 0.70], max: [2.45, 1.67, 1.32]}]; }},
    {nome: 'Carbonato · Al-OH · Fe³⁺', criar: function () {
      return [ee.Image.cat([indice('carbonate_8_6').multiply(-1), indice('alunite_kaolinite_pyrophyllite'),
                            indice('ferric_ox_4_3')]),
              {min: [-0.74, 2.48, 0.61], max: [-0.58, 3.23, 1.67]}]; }}
  ]},
  {grupo: 'Aerogeofísica 1042 (~490 m)', itens: [
    {nome: 'eU (ppm)',  criar: function () { return [gamaPos.select('eU_ppm'), {min: 0.2, max: 3.1, palette: PAL_SEQ}]; }},
    {nome: 'U/Th',      criar: function () { return [gamaPos.select('U_Th'), {min: 0.06, max: 0.85, palette: PAL_RAZAO}]; }},
    {nome: 'U/K',       criar: function () { return [gamaPos.select('U_K'), {min: 0.7, max: 8.7, palette: PAL_RAZAO}]; }},
    {nome: 'eTh (ppm)', criar: function () { return [gamaPos.select('eTh_ppm'), {min: 0.4, max: 19.4, palette: PAL_SEQ}]; }},
    {nome: 'K (%)',     criar: function () { return [gamaPos.select('K_pct'), {min: 0.07, max: 1.4, palette: PAL_SEQ}]; }},
    {nome: 'Ternário K-eTh-eU', criar: function () {
      return [gamaPos.select(['K_pct', 'eTh_ppm', 'eU_ppm']), {min: [0.07, 0.4, 0.2], max: [1.4, 19.4, 3.1]}]; }},
    {nome: 'MAG (nT, ~220 m)', criar: function () {
      return [ee.Image(ACERVO + '1042_mag'), {min: -190, max: 103, palette: PAL_MAG}]; }}
  ]},
  {grupo: 'Geologia (Acervo CPRM)', itens: [
    {nome: 'Una / Salitre — 1:250k', criar: function () {
      return [litoJac.filter(filtroAlvo).style({color: '00ffff', width: 2, fillColor: '00ffff40'}), {}]; }},
    {nome: 'Una / Salitre — 1:1M', criar: function () {
      return [litoBa.filter(filtroAlvo).style({color: '00ffff', width: 1.5, fillColor: '00ffff30'}), {}]; }},
    {nome: 'Litologia Jacobina 1:250k', criar: function () { return [corPorUnidade(litoJac), {}]; }},
    {nome: 'Litologia Bahia 1:1M',      criar: function () { return [corPorUnidade(litoBa), {}]; }},
    {nome: 'Estruturas 1:250k', criar: function () {
      return [ee.FeatureCollection(ACERVO + 'estrut_jacobina_250k').style({color: 'ff00ff', width: 0.8}), {}]; }},
    {nome: 'Estruturas 1:1M', criar: function () {
      return [ee.FeatureCollection(ACERVO + 'estrut_bahia').style({color: 'ffffff', width: 0.8}), {}]; }},
    {nome: 'Dolinas (Morro do Chapéu)', criar: function () {
      return [ee.FeatureCollection(ACERVO + 'dolinas_morro_chapeu').style({color: '3300ff', fillColor: '3300ff80'}), {}]; }}
  ]},
  {grupo: 'Referências', itens: [
    {nome: 'Ocorrências de fosfato', criar: function () {
      return [fosfato.map(function (f) {
        return f.set('estilo', {color: '000000', pointShape: 'star5', pointSize: 7, width: 1,
                                fillColor: ee.Dictionary(CORES_GRUPO).get(f.get('grupo'), 'ffffff')});
      }).style({styleProperty: 'estilo'}), {}]; }},
    {nome: 'Geoquímica P2O5 (rocha)', criar: function () {
      return [geoq.map(function (f) {
        var p = ee.Number(f.get('p2o5_pct'));
        return f.set('estilo', {color: '000000', width: 0.5, pointShape: 'circle',
          fillColor: ee.Algorithms.If(p.gte(5), 'e31a1c', ee.Algorithms.If(p.gte(1), 'fd8d3c', '6baed6')),
          pointSize: ee.Algorithms.If(p.gte(5), 7, ee.Algorithms.If(p.gte(1), 5, 3))});
      }).style({styleProperty: 'estilo'}), {}]; }},
    {nome: 'Recursos minerais (todos)', criar: function () {
      return [ee.FeatureCollection(ACERVO + 'recursos_minerais').style({color: '999999', pointSize: 2, width: 0.5}), {}]; }},
    {nome: 'Grid 5x5', criar: function () { return [ee.Image().paint(grid, 0, 1), {palette: 'ffffff'}]; }},
    {nome: 'Área de estudo', criar: function () { return [ee.Image().paint(area, 0, 2), {palette: 'ffff00'}]; }}
  ]}
];

function indiceSultanR() { return razao('B04', 'B07'); }

// índices individuais entram como um grupo por categoria
var gruposIdx = {};
Object.keys(INDICES).forEach(function (n) {
  var g = 'Índice — ' + INDICES[n].grupo;
  if (!gruposIdx[g]) { gruposIdx[g] = []; }
  gruposIdx[g].push({nome: n, criar: function () {
    var p = INDICES[n];
    return [indice(n), {min: p.min, max: p.max, palette: p.pal || PAL_SEQ}];
  }});
});
Object.keys(gruposIdx).forEach(function (g, i) {
  CATALOGO.splice(2 + i, 0, {grupo: g, itens: gruposIdx[g]});
});

// ----------------------------------------------------------------------------
// Painel de seleção (todas as camadas começam desligadas)
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

var painel = ui.Panel({style: {width: '300px', padding: '6px'}});
painel.add(ui.Label('Mosaico final ASTER · Irecê', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('Marque o que quer ver. Nada é calculado antes de marcar.',
                    {fontSize: '11px', color: '666666'}));

// banda individual
var selBanda = ui.Select({items: BANDAS, value: 'B3N'});
var chkBanda = ui.Checkbox('Banda individual', false);
var camadaBanda = null;
function desenharBanda() {
  if (camadaBanda) { Map.layers().remove(camadaBanda); camadaBanda = null; }
  if (!chkBanda.getValue()) { return; }
  var b = selBanda.getValue();
  camadaBanda = ui.Map.Layer(MOSAICO.select(b), {min: FAIXA_BANDA[b][0], max: FAIXA_BANDA[b][1]},
                             'Banda ' + b, true);
  Map.layers().add(camadaBanda);
}
chkBanda.onChange(desenharBanda);
selBanda.onChange(desenharBanda);
painel.add(ui.Label('Mosaico final — banda individual', {fontWeight: 'bold', margin: '10px 0 2px 0'}));
painel.add(ui.Panel([chkBanda, selBanda], ui.Panel.Layout.flow('horizontal')));

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

var botaoLimpar = ui.Button({label: 'Desligar tudo', style: {margin: '10px 0 0 0'}});
botaoLimpar.onClick(function () {
  Object.keys(camadas).forEach(function (n) { camadas[n].setShown(false); });
  chkBanda.setValue(false);
  painel.widgets().forEach(function (w) {
    if (w instanceof ui.Panel) {
      w.widgets().forEach(function (c) { if (c instanceof ui.Checkbox) { c.setValue(false, false); } });
    }
  });
});
painel.add(botaoLimpar);

// ----------------------------------------------------------------------------
// Inspeção por clique
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '320px', maxHeight: '360px'}});
info.add(ui.Label('Clique no mapa para ver os valores.', {fontSize: '11px'}));
Map.add(info);
Map.style().set('cursor', 'crosshair');

var nomesIdx = Object.keys(INDICES);
Map.onClick(function (c) {
  var pt = ee.Geometry.Point(c.lon, c.lat);
  info.clear();
  info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontWeight: 'bold', fontSize: '12px'}));
  info.add(ui.Label('carregando…', {fontSize: '11px'}));
  var tudo = MOSAICO.addBands(ee.Image.cat(nomesIdx.map(indice))).addBands(gama);
  ee.Dictionary({
    v: tudo.reduceRegion(ee.Reducer.first(), pt, 15),
    unid: litoJac.filterBounds(pt).aggregate_array('NOME_UNIDA').cat(
          litoBa.filterBounds(pt).aggregate_array('NOME_UNIDA'))
  }).evaluate(function (r) {
    info.remove(info.widgets().get(1));
    var v = r.v || {};
    function f(x, d) { return (x === null || x === undefined) ? '—' : Number(x).toFixed(d); }
    info.add(ui.Label('Origem ' + f(v.origem, 0) + (v.swir_estimado ? ' (SWIR estimado)' : ' (SWIR medido)'),
                      {fontSize: '11px'}));
    info.add(ui.Label(BANDAS.map(function (b) { return b + ' ' + f(v[b], 0); }).join(' · '),
                      {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    info.add(ui.Label(nomesIdx.map(function (n) { return n + ' ' + f(v[n], 3); }).join(' · '),
                      {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    info.add(ui.Label('eU ' + f(v.eU_ppm, 2) + ' ppm · U/Th ' + f(v.U_Th, 2) + ' · U/K ' + f(v.U_K, 2),
                      {fontSize: '11px'}));
    info.add(ui.Label('Unidade: ' + (r.unid.length ? r.unid.join(' | ') : '—'), {fontSize: '11px'}));
  });
});

// ----------------------------------------------------------------------------
ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 9);
