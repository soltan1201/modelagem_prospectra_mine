// ============================================================================
// Comparação dos dois mosaicos finais ASTER + índices de mineração + Acervo CPRM
//
// Tela dividida com cortina (wipe): à ESQUERDA o mosaico final v1, à DIREITA
// o mosaico único final. Os dois mapas andam juntos (zoom e posição).
//
// Fontes (projects/mapbiomas-arida/mine/)
//   v1    : aster-mosaico-final (build_mosaico_final.py)
//           'origem' 1–4; 'swir_estimado' 0/1 (sem banda 'data')
//   único : ASTER_mosaico_unico_final (build_mosaico_unico.py etapa 3)
//           'camada' 1 zero-cloud / 2 preenchimento; 'swir_estimado' 0/1;
//           'normalizacao' 0 medido, 1 reta da célula, 2 reta da cena
//   Nos dois: B01–B09 TOA × 1000; B10–B14 TIR em DN.
//   Índices  : mesmas fórmulas de modeling/maxent_fosfato.py (iguais às de
//              visualizar_mosaico_final_acervo.js)
//   Acervo   : projects/mapbiomas-arida/mine/acervo_cprm
//
// Cada linha do painel tem três caixas:
//   v1 → camada no mapa da esquerda; único → no da direita;
//   Δ  → (único − v1) no mapa da direita (só bandas/índices).
// As faixas de cor são as mesmas nos dois lados (percentis 2–98 % do v1,
// 500 m, medidos em 2026-10-12), então a cor compara direto.
// Nada é calculado ao abrir. Clique no mapa: valores dos dois mosaicos.
// ============================================================================

var MINE   = 'projects/mapbiomas-arida/mine/';
var ACERVO = MINE + 'acervo_cprm/';
var MOS = {
  v1:    ee.ImageCollection(MINE + 'aster-mosaico-final').mosaic(),
  unico: ee.ImageCollection(MINE + 'ASTER_mosaico_unico_final').mosaic()
};
var ESCALA_TOA = 1000;

var area = ee.Geometry.Rectangle([-44.460483519818084, -14.790793773956535,
                                  -39.203525512005584, -9.416525285369607]);
var grid = ee.FeatureCollection(MINE + 'grid_irece_5x5_solape100');

var BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
              'B10', 'B11', 'B12', 'B13', 'B14'];
// faixas por banda (VNIR/SWIR × 1000; TIR em DN)
var FAIXA_BANDA = {
  B01: [69, 176], B02: [51, 204], B3N: [155, 332], B04: [167, 362], B05: [77, 232],
  B06: [75, 234], B07: [71, 228], B08: [45, 168], B09: [34, 115],
  B10: [1124, 1635], B11: [1228, 1803], B12: [1332, 1931], B13: [1580, 2203], B14: [1644, 2259]
};

var PAL_SEQ = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
var PAL_VEG = ['8c510a', 'd8b365', 'f6e8c3', 'c7eae5', '5ab4ac', '01665e'];
var PAL_DIV = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];

// ----------------------------------------------------------------------------
// Índices (iguais aos de maxent_fosfato.py), em função do mosaico
// ----------------------------------------------------------------------------
function razao(img, a, b) { return img.select(a).divide(img.select(b)); }

var INDICES = {
  NDVI:                 {f: function (m) { return m.normalizedDifference(['B3N', 'B02']); },
                         min: 0.12, max: 0.68, pal: PAL_VEG, grupo: 'Vegetação'},
  SAVI:                 {f: function (m) {
                           return m.expression('1.5 * (n - r) / (n + r + 0.5)', {
                             n: m.select('B3N').divide(ESCALA_TOA), r: m.select('B02').divide(ESCALA_TOA)}); },
                         min: 0.08, max: 0.41, pal: PAL_VEG, grupo: 'Vegetação'},
  ferric_iron_2_1:      {f: function (m) { return razao(m, 'B02', 'B01'); }, min: 0.70, max: 1.32, grupo: 'Ferro'},
  ferric_ox_4_3:        {f: function (m) { return razao(m, 'B04', 'B3N'); }, min: 0.61, max: 1.67, grupo: 'Ferro'},
  gossan_4_2:           {f: function (m) { return razao(m, 'B04', 'B02'); }, min: 1.39, max: 3.48, grupo: 'Ferro'},
  ferrous_sil_5_4:      {f: function (m) { return razao(m, 'B05', 'B04'); }, min: 0.45, max: 0.67, grupo: 'Ferro'},
  sericite_illite_smectite: {f: function (m) { return m.select('B05').add(m.select('B07')).divide(m.select('B06')); },
                         min: 1.88, max: 2.04, grupo: 'Argilas / Al-OH'},
  alunite_kaolinite_pyrophyllite: {f: function (m) { return m.select('B04').add(m.select('B06')).divide(m.select('B05')); },
                         min: 2.48, max: 3.23, grupo: 'Argilas / Al-OH'},
  muscovite_7_6:        {f: function (m) { return razao(m, 'B07', 'B06'); }, min: 0.90, max: 1.01, grupo: 'Argilas / Al-OH'},
  kaolinite_7_5:        {f: function (m) { return razao(m, 'B07', 'B05'); }, min: 0.88, max: 1.01, grupo: 'Argilas / Al-OH'},
  dolomite:             {f: function (m) { return m.select('B06').add(m.select('B08')).divide(m.select('B07')); },
                         min: 1.65, max: 1.82, grupo: 'Carbonatos'},
  carbonate_8_6:        {f: function (m) { return razao(m, 'B08', 'B06'); }, min: 0.58, max: 0.74, grupo: 'Carbonatos'},
  carbonate_13_14:      {f: function (m) { return razao(m, 'B13', 'B14'); }, min: 0.950, max: 0.983, grupo: 'Carbonatos'},
  quartz_14_12:         {f: function (m) { return razao(m, 'B14', 'B12'); }, min: 1.155, max: 1.249, grupo: 'Sílica / TIR'},
  phosphorite_rbd:      {f: function (m) { return fosforitoTIR(m); }, min: 0.986, max: 1.020, grupo: 'Sílica / TIR'}
};

// Emissividade TIR (normalização de emissividade) → (ε11 + ε13) / (2 ε12)
function fosforitoTIR(m) {
  var UCC = {B10: 6.822e-3, B11: 6.780e-3, B12: 6.590e-3, B13: 5.693e-3, B14: 5.225e-3};
  var LAM = {B10: 8.291, B11: 8.634, B12: 9.075, B13: 10.657, B14: 11.318};
  var C1 = 1.191042e8, C2 = 1.4387752e4, EPS = 0.96;
  var L = {}, T = [];
  Object.keys(UCC).forEach(function (b) {
    L[b] = m.select(b).subtract(1).max(1).multiply(UCC[b]);
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

function indice(m, nome) { return INDICES[nome].f(m).toFloat().rename(nome); }

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
// Catálogo
//   Itens 'mosaico': criar(m) → [imagem, vis] para cada lado; delta(v1, un)
//   opcional → [imagem, vis] da diferença. Itens 'comum': criar() → [obj, vis],
//   vão para os dois mapas.
// ----------------------------------------------------------------------------
function rgb(b) {
  return function (m) {
    return [m.select(b), {bands: b, min: b.map(function (x) { return FAIXA_BANDA[x][0]; }),
                          max: b.map(function (x) { return FAIXA_BANDA[x][1]; })}];
  };
}
function sultanR(m) { return razao(m, 'B04', 'B07'); }

var CATALOGO = [
  {grupo: 'Composições', tipo: 'mosaico', itens: [
    {nome: 'Cor falsa VNIR (3N-2-1)',   criar: rgb(['B3N', 'B02', 'B01'])},
    {nome: 'SWIR (4-6-8)',              criar: rgb(['B04', 'B06', 'B08'])},
    {nome: 'SWIR-VNIR (4-3N-2)',        criar: rgb(['B04', 'B3N', 'B02'])},
    {nome: 'Al-OH / carbonato (5-7-8)', criar: rgb(['B05', 'B07', 'B08'])},
    {nome: 'TIR (13-12-10)',            criar: rgb(['B13', 'B12', 'B10'])}
  ]},
  {grupo: 'Origem do pixel', tipo: 'mosaico', itens: [
    {nome: 'SWIR estimado (magenta)', criar: function (m) {
      return [m.select('swir_estimado').selfMask(), {palette: ['ff00ff']}]; }},
    {nome: 'único: ano da cena', soUnico: true, criar: function (m) {
      return [m.select('data').divide(10000).floor(),
              {min: 2000, max: 2020, palette: ['440154', '3b528b', '21918c', '5ec962', 'fde725']}];}},
    {nome: 'v1: origem 1–4', soV1: true, criar: function (m) {
      return [m.select('origem'), {min: 1, max: 4, palette: ['1a9850', '91cf60', 'fee08b', 'd73027']}]; }},
    {nome: 'único: camada (1 zero-cloud, 2 preench.)', soUnico: true, criar: function (m) {
      return [m.select('camada'), {min: 1, max: 2, palette: ['1a9850', 'fdae61']}]; }},
    {nome: 'único: normalização (0/1/2)', soUnico: true, criar: function (m) {
      return [m.select('normalizacao'), {min: 0, max: 2, palette: ['2166ac', 'fee08b', '1a9850']}]; }}
  ]},
  {grupo: 'Composições de índices', tipo: 'mosaico', itens: [
    {nome: 'Sultan RGB (4/7, 4/1, 2/3·4/3)', criar: function (m) {
      return [ee.Image.cat([sultanR(m), razao(m, 'B04', 'B01'),
                            razao(m, 'B02', 'B3N').multiply(razao(m, 'B04', 'B3N'))]),
              {min: [1.51, 1.51, 0.11], max: [2.45, 2.98, 1.20]}]; }},
    {nome: 'Abrams RGB (4/7, 4/3, 2/1)', criar: function (m) {
      return [ee.Image.cat([sultanR(m), razao(m, 'B04', 'B3N'), razao(m, 'B02', 'B01')]),
              {min: [1.51, 0.61, 0.70], max: [2.45, 1.67, 1.32]}]; }},
    {nome: 'Carbonato · Al-OH · Fe³⁺', criar: function (m) {
      return [ee.Image.cat([indice(m, 'carbonate_8_6').multiply(-1), indice(m, 'alunite_kaolinite_pyrophyllite'),
                            indice(m, 'ferric_ox_4_3')]),
              {min: [-0.74, 2.48, 0.61], max: [-0.58, 3.23, 1.67]}]; }}
  ]},
  {grupo: 'Aerogeofísica 1042 (~490 m)', tipo: 'comum', itens: [
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
  {grupo: 'Geologia (Acervo CPRM)', tipo: 'comum', itens: [
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
  {grupo: 'Referências', tipo: 'comum', itens: [
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

// índices individuais entram como um grupo por categoria, com Δ (único − v1)
var gruposIdx = {};
Object.keys(INDICES).forEach(function (n) {
  var g = 'Índice — ' + INDICES[n].grupo;
  if (!gruposIdx[g]) { gruposIdx[g] = []; }
  var p = INDICES[n];
  var meia = (p.max - p.min) / 2;
  gruposIdx[g].push({
    nome: n,
    criar: function (m) { return [indice(m, n), {min: p.min, max: p.max, palette: p.pal || PAL_SEQ}]; },
    delta: function () {
      return [indice(MOS.unico, n).subtract(indice(MOS.v1, n)), {min: -meia, max: meia, palette: PAL_DIV}];
    }
  });
});
Object.keys(gruposIdx).forEach(function (g, i) {
  CATALOGO.splice(3 + i, 0, {grupo: g, tipo: 'mosaico', itens: gruposIdx[g]});
});

// ----------------------------------------------------------------------------
// Mapas lado a lado
// ----------------------------------------------------------------------------
var mapaV1 = ui.Map();
var mapaUn = ui.Map();
ui.Map.Linker([mapaV1, mapaUn]);
[mapaV1, mapaUn].forEach(function (m) { m.setOptions('HYBRID'); m.style().set('cursor', 'crosshair'); });
mapaV1.add(ui.Label('v1 · aster-mosaico-final', {position: 'top-left', fontWeight: 'bold'}));
mapaUn.add(ui.Label('único · ASTER_mosaico_unico_final', {position: 'top-right', fontWeight: 'bold'}));

var camadas = {};   // chave → [ui.Map.Layer, mapa] (criada na primeira marcação)

function alternar(chave, mapa, fabrica, nome, ligado) {
  if (!camadas[chave]) {
    if (!ligado) { return; }
    var r = fabrica();
    camadas[chave] = ui.Map.Layer(r[0], r[1], nome, true);
    mapa.layers().add(camadas[chave]);
  } else {
    camadas[chave].setShown(ligado);
  }
}

function caixa(rotulo, onChange) {
  return ui.Checkbox({label: rotulo, value: false, style: {fontSize: '11px', margin: '1px 4px 1px 0'},
                      onChange: onChange});
}

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var painel = ui.Panel({style: {width: '360px', padding: '6px'}});
painel.add(ui.Label('Mosaicos finais ASTER · Irecê', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('Esquerda v1, direita único (arraste a cortina). Em cada linha: ' +
                    'v1 → esquerda, único → direita, Δ (único − v1) → direita.',
                    {fontSize: '11px', color: '666666', whiteSpace: 'pre-wrap'}));

// banda individual
var selBanda = ui.Select({items: BANDAS, value: 'B3N'});
var camadasBanda = [];
function desenharBanda() {
  camadasBanda.forEach(function (x) { x[1].layers().remove(x[0]); });
  camadasBanda = [];
  var b = selBanda.getValue();
  var vis = {min: FAIXA_BANDA[b][0], max: FAIXA_BANDA[b][1]};
  var meia = (vis.max - vis.min) / 4;
  [[chkBandaV1, mapaV1, MOS.v1.select(b), vis, 'v1 ' + b],
   [chkBandaUn, mapaUn, MOS.unico.select(b), vis, 'único ' + b],
   [chkBandaD, mapaUn, MOS.unico.select(b).subtract(MOS.v1.select(b)),
    {min: -meia, max: meia, palette: PAL_DIV}, 'Δ ' + b]
  ].forEach(function (p) {
    if (!p[0].getValue()) { return; }
    var l = ui.Map.Layer(p[2], p[3], p[4], true);
    p[1].layers().add(l);
    camadasBanda.push([l, p[1]]);
  });
}
var chkBandaV1 = caixa('v1', desenharBanda);
var chkBandaUn = caixa('único', desenharBanda);
var chkBandaD  = caixa('Δ', desenharBanda);
selBanda.onChange(desenharBanda);
painel.add(ui.Label('Banda individual', {fontWeight: 'bold', margin: '10px 0 2px 0'}));
painel.add(ui.Panel([selBanda, chkBandaV1, chkBandaUn, chkBandaD], ui.Panel.Layout.flow('horizontal')));

var todasCaixas = [chkBandaV1, chkBandaUn, chkBandaD];

CATALOGO.forEach(function (g) {
  var corpo = ui.Panel({style: {shown: false, margin: '0 0 0 8px'}});
  var botao = ui.Button({label: '▸ ' + g.grupo, style: {stretch: 'horizontal', margin: '4px 0 0 0'}});
  botao.onClick(function () {
    var aberto = !corpo.style().get('shown');
    corpo.style().set('shown', aberto);
    botao.setLabel((aberto ? '▾ ' : '▸ ') + g.grupo);
  });
  g.itens.forEach(function (item) {
    var linha = ui.Panel({layout: ui.Panel.Layout.flow('horizontal'), style: {margin: '0'}});
    var cxs = [];
    if (g.tipo === 'comum') {
      cxs.push(caixa('ambos', function (v) {
        alternar('c|' + item.nome + '|v1', mapaV1, item.criar, item.nome, v);
        alternar('c|' + item.nome + '|un', mapaUn, item.criar, item.nome, v);
      }));
    } else {
      if (!item.soUnico) {
        cxs.push(caixa('v1', function (v) {
          alternar(item.nome + '|v1', mapaV1, function () { return item.criar(MOS.v1); }, 'v1 ' + item.nome, v);
        }));
      }
      if (!item.soV1) {
        cxs.push(caixa('único', function (v) {
          alternar(item.nome + '|un', mapaUn, function () { return item.criar(MOS.unico); }, 'único ' + item.nome, v);
        }));
      }
      if (item.delta) {
        cxs.push(caixa('Δ', function (v) {
          alternar(item.nome + '|d', mapaUn, item.delta, 'Δ ' + item.nome, v);
        }));
      }
    }
    cxs.forEach(function (c) { linha.add(c); todasCaixas.push(c); });
    linha.add(ui.Label(item.nome, {fontSize: '12px', margin: '2px 0 0 4px'}));
    corpo.add(linha);
  });
  painel.add(botao);
  painel.add(corpo);
});

var botaoLimpar = ui.Button({label: 'Desligar tudo', style: {margin: '10px 0 0 0'}});
botaoLimpar.onClick(function () {
  Object.keys(camadas).forEach(function (k) { camadas[k].setShown(false); });
  todasCaixas.forEach(function (c) { c.setValue(false, false); });
  desenharBanda();
});
painel.add(botaoLimpar);

// ----------------------------------------------------------------------------
// Dispersão v1 × único para um índice ou banda (sob demanda)
// ----------------------------------------------------------------------------
var nomesIdx = Object.keys(INDICES);
var selDisp = ui.Select({items: nomesIdx.concat(BANDAS), value: 'carbonate_8_6'});
var saida = ui.Panel();
function valorDe(m, n) { return INDICES[n] ? indice(m, n) : m.select(n).toFloat(); }
function graficoDispersao() {
  var n = selDisp.getValue();
  var pts = ee.Image.cat([valorDe(MOS.v1, n).rename('v1'), valorDe(MOS.unico, n).rename('unico')])
              .sample({region: area, scale: 120, numPixels: 20000, seed: 42, tileScale: 4})
              .limit(3000);
  saida.clear();
  saida.add(ui.Chart.feature.byFeature(pts, 'v1', 'unico').setChartType('ScatterChart').setOptions({
    title: n + ' — único × v1', pointSize: 1, legend: {position: 'none'},
    hAxis: {title: 'v1'}, vAxis: {title: 'único'},
    trendlines: {0: {type: 'linear', showR2: true, visibleInLegend: true, color: 'red'}}
  }));
}
painel.add(ui.Label('Análise', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
painel.add(ui.Panel([selDisp, ui.Button({label: 'Dispersão único × v1', onClick: graficoDispersao})],
                    ui.Panel.Layout.flow('horizontal')));
painel.add(saida);

// ----------------------------------------------------------------------------
// Inspeção por clique (nos dois mapas)
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '360px', maxHeight: '400px'}});
info.add(ui.Label('Clique no mapa para ver os valores.', {fontSize: '11px'}));
mapaUn.add(info);

function aoClicar(c) {
  var pt = ee.Geometry.Point(c.lon, c.lat);
  info.clear();
  info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontWeight: 'bold', fontSize: '12px'}));
  info.add(ui.Label('carregando…', {fontSize: '11px'}));
  function lerEm(m) {
    return m.addBands(ee.Image.cat(nomesIdx.map(function (n) { return indice(m, n); })))
            .reduceRegion(ee.Reducer.first(), pt, 15);
  }
  ee.Dictionary({
    v1: lerEm(MOS.v1), un: lerEm(MOS.unico),
    g: gama.reduceRegion(ee.Reducer.first(), pt, 15),
    unid: litoJac.filterBounds(pt).aggregate_array('NOME_UNIDA').cat(
          litoBa.filterBounds(pt).aggregate_array('NOME_UNIDA'))
  }).evaluate(function (r) {
    info.remove(info.widgets().get(1));
    var a = r.v1 || {}, b = r.un || {}, g = r.g || {};
    function f(x, d) { return (x === null || x === undefined) ? '—' : Number(x).toFixed(d); }
    info.add(ui.Label('v1: origem ' + f(a.origem, 0) +
                      (a.swir_estimado ? ' (SWIR estimado)' : ' (SWIR medido)'), {fontSize: '11px', fontWeight: 'bold'}));
    info.add(ui.Label('único: data ' + f(b.data, 0) + ' · camada ' + f(b.camada, 0) + ' · norm. ' +
                      f(b.normalizacao, 0) + (b.swir_estimado ? ' (SWIR estimado)' : ' (SWIR medido)'),
                      {fontSize: '11px', fontWeight: 'bold'}));
    info.add(ui.Label('banda: v1 / único', {fontSize: '11px', color: '666666'}));
    info.add(ui.Label(BANDAS.map(function (x) { return x + ' ' + f(a[x], 0) + '/' + f(b[x], 0); }).join(' · '),
                      {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    info.add(ui.Label(nomesIdx.map(function (n) { return n + ' ' + f(a[n], 3) + '/' + f(b[n], 3); }).join(' · '),
                      {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    info.add(ui.Label('eU ' + f(g.eU_ppm, 2) + ' ppm · U/Th ' + f(g.U_Th, 2) + ' · U/K ' + f(g.U_K, 2),
                      {fontSize: '11px'}));
    info.add(ui.Label('Unidade: ' + (r.unid.length ? r.unid.join(' | ') : '—'), {fontSize: '11px'}));
  });
}
mapaV1.onClick(aoClicar);
mapaUn.onClick(aoClicar);

// ----------------------------------------------------------------------------
var divisao = ui.SplitPanel({firstPanel: mapaV1, secondPanel: mapaUn, wipe: true,
                             style: {stretch: 'both'}});
ui.root.widgets().reset([painel, divisao]);
ui.root.setLayout(ui.Panel.Layout.flow('horizontal'));
mapaV1.centerObject(ee.Geometry.Point(-41.86, -11.30), 9);
