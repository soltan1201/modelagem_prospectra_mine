// ============================================================================
// Maxent por alvo mineral — versão Code Editor de maxent_fosfato.py, para ver
// o mapa antes de rodar o Python e salvar.
//
// Features: FEATURES_DL do SHAP (padrão, 21) ou as 36 do script (bandas,
// índices, phosphorite_rbd, GLCM da B3N); mesmo
// Maxent (linear + quadrática, betaMultiplier 2) e mesmos pontos
// (pontos_referencia, filtrado por 'alvo' e, opcionalmente, papel = treino).
//
// Diferença na amostragem (para caber no modo interativo):
//   presença  AMOSTRAS_POR_PONTO pontos aleatórios no buffer de BUFFER_M de
//             cada ocorrência (randomPoints + sampleRegions numa chamada só)
//   fundo     N_FUNDO pontos aleatórios na grid
//   O Python sorteia pixels com sample() por ponto/célula; os números saem
//   parecidos, não idênticos.
//
// Uso: escolha alvo, mosaico e features no painel → "Treinar e mostrar".
// Os botões de export criam as tasks (Run na aba Tasks), com os mesmos
// nomes de asset do Python: maxent_<alvo>_amostras_<mosaico>_<modo> e
// maxent_<alvo>_<mosaico> (bandas probability e complemento = 1 − p).
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
var ASSET_PONTOS = MINE + 'pontos_referencia';
var GRID = ee.FeatureCollection(MINE + 'grid_irece_5x5_solape100');
var MOSAICOS = {
  'ASTER_mosaico_unico_final': MINE + 'ASTER_mosaico_unico_final',
  'aster-mosaico-final': MINE + 'aster-mosaico-final'
};
var ALVOS = ['fosfato', 'carbonato', 'mvt', 'ouro', 'manganes'];
var MODOS = ['FEATURES_DL', 'bandas+indices', 'bandas', 'indices'];
// FEATURES_DL do SHAP (avaliar_maxent.ipynb, célula 34; ranking_features_shap_20261021.csv):
// 21 de 36, em ordem de importância, sem pares com |r| > 0,9 — a mesma lista do
// exportar_favorabilidade_ml.py
var FEATURES_DL = [
  'quartz_14_12', 'B10', 'ferric_iron_2_1', 'glcm_contrast', 'B09', 'sultan_g',
  'phosphorite_rbd', 'gossan_4_2', 'sultan_r', 'ferric_ox_4_3', 'SAVI', 'B01', 'B3N',
  'carbonate_13_14', 'kaolinite_7_5', 'carbonate_8_6', 'glcm_hom', 'dolomite',
  'sericite_illite_smectite', 'muscovite_7_6', 'glcm_ent'
];

var BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
              'B10', 'B11', 'B12', 'B13', 'B14'];
var BUFFER_M = 500;
var AMOSTRAS_POR_PONTO = 20;
var N_FUNDO = 10000;
var ESCALA = 30;
var SEED = 42;
var ESCALA_TOA = 1000;

var MAXENT = {
  autoFeature: false, linear: true, quadratic: true,
  product: false, threshold: false, hinge: false,
  betaMultiplier: 2.0, addSamplesToBackground: true, seed: SEED
};

var regiao = GRID.geometry().dissolve(100);

// ----------------------------------------------------------------------------
// Features (iguais às de maxent_fosfato.py)
// ----------------------------------------------------------------------------
function razao(a, b) { return function (i) { return i.select(a).divide(i.select(b)); }; }

var INDICES = {
  NDVI: function (i) { return i.normalizedDifference(['B3N', 'B02']); },
  SAVI: function (i) {
    return i.expression('1.5 * (n - r) / (n + r + 0.5)',
      {n: i.select('B3N').divide(ESCALA_TOA), r: i.select('B02').divide(ESCALA_TOA)}); },
  ferric_iron_2_1: razao('B02', 'B01'),
  ferric_ox_4_3:   razao('B04', 'B3N'),
  gossan_4_2:      razao('B04', 'B02'),
  ferrous_sil_5_4: razao('B05', 'B04'),
  sericite_illite_smectite: function (i) { return i.select('B05').add(i.select('B07')).divide(i.select('B06')); },
  alunite_kaolinite_pyrophyllite: function (i) { return i.select('B04').add(i.select('B06')).divide(i.select('B05')); },
  muscovite_7_6:   razao('B07', 'B06'),
  kaolinite_7_5:   razao('B07', 'B05'),
  dolomite:        function (i) { return i.select('B06').add(i.select('B08')).divide(i.select('B07')); },
  carbonate_8_6:   razao('B08', 'B06'),
  carbonate_13_14: razao('B13', 'B14'),
  quartz_14_12:    razao('B14', 'B12'),
  sultan_r:        razao('B04', 'B07'),
  sultan_g:        razao('B04', 'B01'),
  sultan_b:        function (i) { return i.select('B02').divide(i.select('B3N'))
                                           .multiply(i.select('B04').divide(i.select('B3N'))); }
};
// GLCM da B3N: nome no GEE → nome da feature
var TEXTURAS = {contrast: 'contrast', 'var': 'var', ent: 'ent', idm: 'hom'};

// Emissividade TIR (normalização, ε_max 0,96) → (ε11 + ε13) / (2 ε12)
function fosforitoTIR(img) {
  var UCC = {B10: 6.822e-3, B11: 6.780e-3, B12: 6.590e-3, B13: 5.693e-3, B14: 5.225e-3};
  var LAM = {B10: 8.291, B11: 8.634, B12: 9.075, B13: 10.657, B14: 11.318};
  var C1 = 1.191042e8, C2 = 1.4387752e4, EPS = 0.96;
  var L = {}, T = [];
  Object.keys(UCC).forEach(function (b) {
    L[b] = img.select(b).subtract(1).max(1).multiply(UCC[b]);
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

// Os dois mosaicos finais têm B01–B14; por isso não há checagem de bandas aqui.
function calcularFeatures(img, modo) {
  var r = (modo === 'FEATURES_DL')
    ? {img: featuresBrutas(img, 'bandas+indices').img.select(FEATURES_DL), nomes: FEATURES_DL}
    : featuresBrutas(img, modo);
  // pixel só entra com todas as features escolhidas válidas
  return {img: r.img.updateMask(r.img.mask().reduce(ee.Reducer.min())), nomes: r.nomes};
}

function featuresBrutas(img, modo) {
  img = img.select(BANDAS).toFloat();
  var partes = [], nomes = [];
  if (modo.indexOf('bandas') >= 0) { partes.push(img); nomes = nomes.concat(BANDAS); }
  if (modo.indexOf('indices') >= 0) {
    Object.keys(INDICES).forEach(function (n) {
      partes.push(INDICES[n](img).toFloat().rename(n));
      nomes.push(n);
    });
    partes.push(fosforitoTIR(img).toFloat().rename('phosphorite_rbd'));
    nomes.push('phosphorite_rbd');
    var b3n = img.select('B3N').unitScale(0, 0.35 * ESCALA_TOA).multiply(255).clamp(0, 255).toByte();
    var glcm = b3n.glcmTexture({size: 1});
    var de = Object.keys(TEXTURAS).map(function (t) { return 'B3N_' + t; });
    var para = Object.keys(TEXTURAS).map(function (t) { return 'glcm_' + TEXTURAS[t]; });
    partes.push(glcm.select(de, para).toFloat());
    nomes = nomes.concat(para);
  }
  return {img: ee.Image.cat(partes).toFloat(), nomes: nomes};
}

// ----------------------------------------------------------------------------
// Amostras e modelo
// ----------------------------------------------------------------------------
function lerPontos(alvo, soTreino) {
  var p = ee.FeatureCollection(ASSET_PONTOS).filter(ee.Filter.eq('alvo', alvo)).filterBounds(regiao);
  return soTreino ? p.filter(ee.Filter.eq('papel', 'treino')) : p;
}

function amostrasPresenca(feats, pontos) {
  var aleatorios = pontos.map(function (f) {
    var props = {classe: 1, CODIGO_OBJ: f.get('CODIGO_OBJ')};
    return ee.FeatureCollection.randomPoints(f.geometry().buffer(BUFFER_M), AMOSTRAS_POR_PONTO, SEED)
             .map(function (s) { return s.set(props); });
  }).flatten();
  return feats.sampleRegions({collection: aleatorios, scale: ESCALA, geometries: true, tileScale: 4});
}

function amostrasFundo(feats) {
  var aleatorios = ee.FeatureCollection.randomPoints(regiao, N_FUNDO, SEED)
                     .map(function (s) { return s.set('classe', 0); });
  return feats.sampleRegions({collection: aleatorios, scale: ESCALA, geometries: true, tileScale: 4});
}

function treinar(amostras, nomes) {
  return ee.Classifier.amnhMaxent(MAXENT).setOutputMode('PROBABILITY')
           .train({features: amostras, classProperty: 'classe', inputProperties: nomes});
}

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var PAL_P = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
var CORES_PAPEL = {treino: 'ff0000', validacao: '00ffff'};

var selAlvo = ui.Select({items: ALVOS, value: 'fosfato', style: {stretch: 'horizontal'}});
var selMos  = ui.Select({items: Object.keys(MOSAICOS), value: 'ASTER_mosaico_unico_final',
                         style: {stretch: 'horizontal'}});
var selModo = ui.Select({items: MODOS, value: 'FEATURES_DL', style: {stretch: 'horizontal'}});
var chkTreino = ui.Checkbox('Só papel = treino (tira a validação R2 do fosfato)', false);
var saida = ui.Panel();
var atual = null;   // {alvo, mosaico, modo, feats, nomes, amostras, mapa}

var painel = ui.Panel({style: {width: '340px', padding: '6px'}});
painel.add(ui.Label('Maxent por alvo · Irecê', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('Prévia do maxent_fosfato.py. Nada é salvo até clicar em exportar.',
                    {fontSize: '11px', color: '666666'}));
painel.add(ui.Label('Alvo', {fontSize: '12px', margin: '6px 0 0 8px'}));
painel.add(selAlvo);
painel.add(ui.Label('Mosaico', {fontSize: '12px', margin: '6px 0 0 8px'}));
painel.add(selMos);
painel.add(ui.Label('Features', {fontSize: '12px', margin: '6px 0 0 8px'}));
painel.add(selModo);
painel.add(chkTreino);

function nomeMosaico() { return selMos.getValue(); }

function treinarEMostrar() {
  Map.layers().reset();
  saida.clear();
  var alvo = selAlvo.getValue(), modo = selModo.getValue();
  var mosaico = ee.ImageCollection(MOSAICOS[nomeMosaico()]).mosaic();
  var f = calcularFeatures(mosaico, modo);
  var pontos = lerPontos(alvo, chkTreino.getValue());
  var pres = amostrasPresenca(f.img, pontos);
  var fundo = amostrasFundo(f.img);
  var amostras = pres.merge(fundo);
  var clf = treinar(amostras, f.nomes);
  var p = f.img.classify(clf).select('probability').toFloat();
  var mapa = p.addBands(ee.Image(1).subtract(p).rename('complemento'));
  atual = {alvo: alvo, modo: modo, mosaico: nomeMosaico(), nomes: f.nomes, amostras: amostras,
           mapa: mapa, pres: pres, fundo: fundo, clf: clf};

  Map.addLayer(mosaico.select(['B04', 'B06', 'B08']), {min: [167, 75, 45], max: [362, 234, 168]},
               'Mosaico SWIR (4-6-8)', false);
  Map.addLayer(p, {min: 0, max: 1, palette: PAL_P}, alvo + ' — probability');
  Map.addLayer(mapa.select('complemento'), {min: 0, max: 1, palette: PAL_P}, alvo + ' — complemento (1 − p)', false);
  Map.addLayer(pres, {color: 'ffffff'}, 'Amostras de presença', false);
  Map.addLayer(pontos.map(function (x) {
    return x.set('estilo', {color: '000000', pointShape: 'star5', pointSize: 7, width: 1,
                            fillColor: ee.Dictionary(CORES_PAPEL).get(x.get('papel'), 'ffffff')});
  }).style({styleProperty: 'estilo'}), {}, 'Pontos ' + alvo + ' (vermelho treino, ciano validação)');
  Map.addLayer(ee.Image().paint(GRID, 0, 1), {palette: 'ffffff'}, 'Grid', false);

  saida.add(ui.Label('Contando pontos e amostras…', {fontSize: '11px'}));
  ee.Dictionary({
    pontos: pontos.size(),
    fonte: pontos.aggregate_histogram('fonte'),
    papel: pontos.aggregate_histogram('papel'),
    n_pres: pres.size(), n_fundo: fundo.size()
  }).evaluate(function (r, erro) {
    saida.clear();
    if (erro) { saida.add(ui.Label('Erro: ' + erro, {color: 'red', fontSize: '11px'})); return; }
    saida.add(ui.Label(alvo + ': ' + r.pontos + ' pontos na grid · fonte ' + JSON.stringify(r.fonte) +
                       ' · papel ' + JSON.stringify(r.papel) + '\nAmostras: presença ' + r.n_pres +
                       ', fundo ' + r.n_fundo + ' · ' + f.nomes.length + ' features',
                       {fontSize: '11px', whiteSpace: 'pre-wrap'}));
  });
}

function histogramaProb() {
  if (!atual) { return; }
  var cls = atual.amostras.classify(atual.clf);
  saida.add(ui.Chart.feature.histogram({features: cls.filter(ee.Filter.eq('classe', 1)),
                                        property: 'probability', minBucketWidth: 0.05})
    .setOptions({title: 'Probability nas amostras de presença (treino)', legend: {position: 'none'}}));
  saida.add(ui.Chart.feature.histogram({features: cls.filter(ee.Filter.eq('classe', 0)),
                                        property: 'probability', minBucketWidth: 0.05})
    .setOptions({title: 'Probability no fundo', legend: {position: 'none'}}));
}

function exportarAmostras() {
  if (!atual) { return; }
  var nome = 'maxent_' + atual.alvo + '_amostras_' + atual.mosaico + '_' + atual.modo.replace('+', '_');
  Export.table.toAsset({collection: atual.amostras.map(function (x) { return x.set('features_modo', atual.modo); }),
                        description: nome, assetId: MINE + nome});
  saida.add(ui.Label('Task criada: ' + nome + ' (Run na aba Tasks)', {fontSize: '11px'}));
}

function exportarMapa() {
  if (!atual) { return; }
  var nome = 'maxent_' + atual.alvo + '_' + atual.mosaico;
  Export.image.toAsset({
    image: atual.mapa.set({mosaico: MOSAICOS[atual.mosaico], pontos: ASSET_PONTOS, alvo: atual.alvo,
                           papeis: chkTreino.getValue() ? 'treino' : 'todos', bandas: atual.nomes.join(','),
                           buffer_m: BUFFER_M, amostras_por_ponto: AMOSTRAS_POR_PONTO, n_fundo: N_FUNDO,
                           origem: 'maxent_alvos.js'}),
    description: nome, assetId: MINE + nome, region: regiao, scale: ESCALA, crs: 'EPSG:4326',
    maxPixels: 1e13, pyramidingPolicy: {'.default': 'mean'}});
  saida.add(ui.Label('Task criada: ' + nome + ' (Run na aba Tasks)', {fontSize: '11px'}));
}

painel.add(ui.Button({label: 'Treinar e mostrar', onClick: treinarEMostrar, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'Histograma da probability', onClick: histogramaProb, style: {stretch: 'horizontal'}}));
painel.add(ui.Label('Salvar', {fontWeight: 'bold', margin: '12px 0 2px 0'}));
painel.add(ui.Button({label: 'Exportar amostras (asset)', onClick: exportarAmostras, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'Exportar mapa (asset)', onClick: exportarMapa, style: {stretch: 'horizontal'}}));
painel.add(saida);

// ----------------------------------------------------------------------------
// Clique: probability no ponto
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '260px'}});
info.add(ui.Label('Treine e clique no mapa.', {fontSize: '11px'}));
Map.add(info);
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  if (!atual) { return; }
  info.clear();
  info.add(ui.Label('carregando…', {fontSize: '11px'}));
  atual.mapa.reduceRegion(ee.Reducer.first(), ee.Geometry.Point(c.lon, c.lat), ESCALA)
    .evaluate(function (r) {
      info.clear();
      r = r || {};
      var v = r.probability;
      info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontSize: '11px', fontWeight: 'bold'}));
      info.add(ui.Label(atual.alvo + ': p = ' + (v === null || v === undefined ? '—' : v.toFixed(3)),
                        {fontSize: '12px'}));
    });
});

ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);
