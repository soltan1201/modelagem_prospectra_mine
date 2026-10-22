// ============================================================================
// Baixar os mosaicos finais ASTER + índices espectrais (calculo_spectrals_index.js)
// para o Google Drive, em GeoTIFF, por célula da grid 5×5.
//
// Mosaicos (projects/mapbiomas-arida/mine/)
//   ASTER_mosaico_unico_final   etapa 3 do build_mosaico_unico.py
//   aster-mosaico-final         mosaico final v1 (build_mosaico_final.py)
//   B01–B09 reflectância TOA × 1000; B10–B14 TIR em DN.
//
// Índices — os de calculo_spectrals_index.js:
//   addMineralIndices    ferric_iron_2_1, ferric_ox_4_3, gossan_4_2,
//                        ferrous_sil_5_4, sericite_illite_smectite,
//                        alunite_kaolinite_pyrophyllite, muscovite_7_6,
//                        kaolinite_7_5, dolomite, carbonate_13_14,
//                        quartz_14_12, sultan_r/g/b, abrams_g/b
//   adicionarIndicesGeo  NDVI, AlOH_Clay, Carbonate (B08/B06), SAVI
//   texturas (GLCM B3N)  glcm_contrast, glcm_var, glcm_ent, glcm_hom
//   + phosphorite_rbd    (ε11 + ε13) / (2·ε12), feição do apatita em ~9 µm
//                        (Guha et al. 2019); emissividade por normalização
//                        (ε_max 0,96) sobre a radiância TOA — igual ao
//                        maxent_fosfato.py (não está em calculo_spectrals_index.js)
//   Repetidos no JS original entram uma vez: abrams_r = sultan_r (B04/B07) e
//   Ferric_Fe = ferric_iron_2_1 (B02/B01).
//   SAVI: o JS original usa os valores × 1000 direto com L = 0,5, o que anula
//   o fator de solo; aqui a reflectância é dividida por 1000 antes (como no
//   maxent_fosfato.py). GLCM: B3N (× 1000) 0–350 → 0–255 em byte.
//
// Arquivos (um por célula e por grupo, para não ficarem enormes)
//   <mosaico>_<LxCy>_bandas   B01–B14 + bandas de origem do mosaico, Int32
//   <mosaico>_<LxCy>_indices  índices (+ texturas, se marcadas), Float32
//
// Tamanho aproximado por célula (≈ 115 × 118 km): a 30 m, bandas ≈ 0,6–1 GB e
// índices ≈ 1,2 GB; a 90 m, ≈ 1/9 disso. O Drive divide arquivos grandes em
// pedaços (fileDimensions).
//
// Uso: escolha no painel → "Criar tasks" → Run em cada task na aba Tasks.
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
var MOSAICOS = {
  'ASTER_mosaico_unico_final': {
    col: ee.ImageCollection(MINE + 'ASTER_mosaico_unico_final'),
    origem: ['quality', 'data', 'camada', 'swir_estimado', 'normalizacao']
  },
  'aster-mosaico-final': {
    col: ee.ImageCollection(MINE + 'aster-mosaico-final'),
    origem: ['origem', 'swir_estimado']
  }
};
var GRID = ee.FeatureCollection(MINE + 'grid_irece_5x5');   // sem solape: arquivos não se repetem
var BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
              'B10', 'B11', 'B12', 'B13', 'B14'];
var ESCALA_TOA = 1000;

// ----------------------------------------------------------------------------
// Índices (calculo_spectrals_index.js)
// ----------------------------------------------------------------------------
function indicesMinerais(img) {
  var b = {};
  BANDAS.forEach(function (n) { b[n] = img.select(n).toFloat(); });
  return ee.Image.cat([
    b.B02.divide(b.B01).rename('ferric_iron_2_1'),
    b.B04.divide(b.B3N).rename('ferric_ox_4_3'),
    b.B04.divide(b.B02).rename('gossan_4_2'),
    b.B05.divide(b.B04).rename('ferrous_sil_5_4'),
    b.B05.add(b.B07).divide(b.B06).rename('sericite_illite_smectite'),
    b.B04.add(b.B06).divide(b.B05).rename('alunite_kaolinite_pyrophyllite'),
    b.B07.divide(b.B06).rename('muscovite_7_6'),
    b.B07.divide(b.B05).rename('kaolinite_7_5'),
    b.B06.add(b.B08).divide(b.B07).rename('dolomite'),
    b.B13.divide(b.B14).rename('carbonate_13_14'),
    b.B14.divide(b.B12).rename('quartz_14_12'),
    b.B04.divide(b.B07).rename('sultan_r'),
    b.B04.divide(b.B01).rename('sultan_g'),
    b.B02.divide(b.B3N).multiply(b.B04.divide(b.B3N)).rename('sultan_b'),
    b.B04.divide(b.B3N).rename('abrams_g'),
    b.B02.divide(b.B01).rename('abrams_b')
  ]);
}

function indicesGeo(img) {
  var n = img.select('B3N').toFloat().divide(ESCALA_TOA);
  var r = img.select('B02').toFloat().divide(ESCALA_TOA);
  return ee.Image.cat([
    img.normalizedDifference(['B3N', 'B02']).rename('NDVI'),
    img.expression("(b('B05') + b('B07')) / (2 * b('B06'))").rename('AlOH_Clay'),
    img.expression("b('B08') / b('B06')").rename('Carbonate'),
    n.subtract(r).multiply(1.5).divide(n.add(r).add(0.5)).rename('SAVI')
  ]).toFloat();
}

// Emissividade TIR (normalização de emissividade) → (ε11 + ε13) / (2 ε12)
function fosforitoTIR(img) {
  var UCC = {B10: 6.822e-3, B11: 6.780e-3, B12: 6.590e-3, B13: 5.693e-3, B14: 5.225e-3};
  var LAM = {B10: 8.291, B11: 8.634, B12: 9.075, B13: 10.657, B14: 11.318};
  var C1 = 1.191042e8, C2 = 1.4387752e4, EPS = 0.96;
  var L = {}, T = [];
  Object.keys(UCC).forEach(function (b) {
    L[b] = img.select(b).toFloat().subtract(1).max(1).multiply(UCC[b]);
    T.push(L[b].expression('c2 / (lam * log(eps * c1 / (pow(lam, 5) * L) + 1))',
                           {c1: C1, c2: C2, lam: LAM[b], eps: EPS, L: L[b]}));
  });
  var tmax = ee.Image.cat(T).reduce(ee.Reducer.max());
  function e(b) {
    return L[b].divide(tmax.expression('c1 / (pow(lam, 5) * (exp(c2 / (lam * t)) - 1))',
                                       {c1: C1, c2: C2, lam: LAM[b], t: tmax}));
  }
  return e('B11').add(e('B13')).divide(e('B12').multiply(2)).toFloat().rename('phosphorite_rbd');
}

function texturas(img) {
  var b3n = img.select('B3N').unitScale(0, 0.35 * ESCALA_TOA).multiply(255).clamp(0, 255).toByte();
  return b3n.glcmTexture({size: 1})
            .select(['B3N_contrast', 'B3N_var', 'B3N_ent', 'B3N_idm'],
                    ['glcm_contrast', 'glcm_var', 'glcm_ent', 'glcm_hom']).toFloat();
}

function imagemIndices(img, comTextura) {
  var out = indicesGeo(img).addBands(indicesMinerais(img)).addBands(fosforitoTIR(img));
  return comTextura ? out.addBands(texturas(img)) : out;
}

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var painel = ui.Panel({style: {width: '340px', padding: '6px'}});
painel.add(ui.Label('Baixar mosaicos + índices', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('GeoTIFF no Drive, um arquivo por célula e grupo. As tasks precisam de Run na aba Tasks.',
                    {fontSize: '11px', color: '666666'}));

function rotulo(t) { return ui.Label(t, {fontSize: '12px', fontWeight: 'bold', margin: '8px 0 0 8px'}); }

painel.add(rotulo('Mosaicos'));
var chkMos = {};
Object.keys(MOSAICOS).forEach(function (m) {
  chkMos[m] = ui.Checkbox(m, true);
  painel.add(chkMos[m]);
});

painel.add(rotulo('O que baixar'));
var chkBandas  = ui.Checkbox('Bandas B01–B14 + bandas de origem (Int32)', true);
var chkIndices = ui.Checkbox('Índices espectrais (Float32)', true);
var chkTextura = ui.Checkbox('… incluir texturas GLCM nos índices', false);
[chkBandas, chkIndices, chkTextura].forEach(function (w) { painel.add(w); });

painel.add(rotulo('Células'));
var txtCelulas = ui.Textbox({placeholder: 'vazio = todas; ex.: L2C2, L3C1', style: {stretch: 'horizontal'}});
painel.add(txtCelulas);

painel.add(rotulo('Escala (m) e pasta no Drive'));
var selEscala = ui.Select({items: ['15', '30', '90'], value: '30'});
var txtPasta  = ui.Textbox({value: 'aster_mosaicos_indices', style: {stretch: 'horizontal'}});
painel.add(ui.Panel([selEscala, txtPasta], ui.Panel.Layout.flow('horizontal')));

var saida = ui.Label('', {fontSize: '11px', whiteSpace: 'pre-wrap'});

function celulasEscolhidas() {
  var t = txtCelulas.getValue().trim();
  if (!t) { return GRID; }
  var nomes = t.split(/[\s,;]+/).filter(function (x) { return x; }).map(function (x) { return x.toUpperCase(); });
  return GRID.filter(ee.Filter.inList('nome', nomes));
}

function criarTasks() {
  var mosaicos = Object.keys(MOSAICOS).filter(function (m) { return chkMos[m].getValue(); });
  if (!mosaicos.length || (!chkBandas.getValue() && !chkIndices.getValue())) {
    saida.setValue('Marque ao menos um mosaico e um grupo.');
    return;
  }
  var escala = Number(selEscala.getValue());
  var pasta = txtPasta.getValue();
  saida.setValue('Lendo a grid…');
  celulasEscolhidas().aggregate_array('nome').evaluate(function (nomes) {
    var n = 0;
    nomes.sort().forEach(function (cel) {
      var geom = GRID.filter(ee.Filter.eq('nome', cel)).geometry();
      mosaicos.forEach(function (m) {
        var img = MOSAICOS[m].col.mosaic().clip(geom);
        var base = m + '_' + cel;
        if (chkBandas.getValue()) {
          Export.image.toDrive({
            image: img.select(BANDAS.concat(MOSAICOS[m].origem)).toInt32(),   // 'data' AAAAMMDD não cabe em Int16
            description: base + '_bandas', fileNamePrefix: base + '_bandas', folder: pasta,
            region: geom, scale: escala, crs: 'EPSG:4326', maxPixels: 1e13, fileFormat: 'GeoTIFF',
            formatOptions: {cloudOptimized: true}
          });
          n++;
        }
        if (chkIndices.getValue()) {
          Export.image.toDrive({
            image: imagemIndices(img, chkTextura.getValue()),
            description: base + '_indices', fileNamePrefix: base + '_indices', folder: pasta,
            region: geom, scale: escala, crs: 'EPSG:4326', maxPixels: 1e13, fileFormat: 'GeoTIFF',
            formatOptions: {cloudOptimized: true}
          });
          n++;
        }
      });
    });
    saida.setValue(n + ' tasks criadas (' + nomes.length + ' células × ' + mosaicos.length +
                   ' mosaico(s)) → Drive/' + pasta + ', ' + escala + ' m.\nClique em Run em cada uma na aba Tasks.');
  });
}

painel.add(ui.Button({label: 'Criar tasks', onClick: criarTasks, style: {stretch: 'horizontal', margin: '10px 0 0 0'}}));
painel.add(saida);

// ----------------------------------------------------------------------------
// Prévia (para conferir antes de baixar)
// ----------------------------------------------------------------------------
var FAIXA = {B01: [69, 176], B02: [51, 204], B3N: [155, 332], B04: [167, 362], B06: [75, 234], B08: [45, 168]};
function rgb(img, b) {
  return {bands: b, min: b.map(function (x) { return FAIXA[x][0]; }), max: b.map(function (x) { return FAIXA[x][1]; })};
}
var PAL = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
Object.keys(MOSAICOS).forEach(function (m) {
  var img = MOSAICOS[m].col.mosaic();
  var idx = imagemIndices(img, false);
  Map.addLayer(img, rgb(img, ['B3N', 'B02', 'B01']), m + ' — VNIR', false);
  Map.addLayer(img, rgb(img, ['B04', 'B06', 'B08']), m + ' — SWIR', false);
  Map.addLayer(idx.select(['sultan_r', 'sultan_g', 'sultan_b']), {min: [1.51, 1.51, 0.11], max: [2.45, 2.98, 1.20]},
               m + ' — Sultan RGB', false);
  Map.addLayer(idx.select('Carbonate'), {min: 0.58, max: 0.74, palette: PAL.slice().reverse()}, m + ' — Carbonate (B08/B06)', false);
  Map.addLayer(idx.select('AlOH_Clay'), {min: 0.94, max: 1.02, palette: PAL}, m + ' — AlOH_Clay', false);
  Map.addLayer(idx.select('phosphorite_rbd'), {min: 0.986, max: 1.020, palette: PAL}, m + ' — phosphorite_rbd', false);
});
Map.addLayer(ee.Image().paint(GRID, 0, 1), {palette: 'ffffff'}, 'Grid 5×5');
Map.addLayer(GRID.map(function (f) { return ee.Feature(f.geometry().centroid(100), {nome: f.get('nome')}); })
               .style({color: 'ffff00', pointSize: 3}), {}, 'Centro das células', false);

ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(GRID, 8);
