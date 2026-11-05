// ============================================================================
// Apresentação: mosaico final sazonal seco ORIGINAL × HARMONIZADO com o Landsat
//
// Os dois lados juntam A (2000–2008, SWIR medido) e B (2009–2020) na hora, com a mesma
// regra de build_mosaico_final_v2.py (montar_nuvem), sem exportar nada:
//   origem 1  A por cima;
//   origem 2  B onde falta A;
//   origem 3  B no lugar do A onde o A parece nuvem frente ao B:
//             A mais claro no B01 (B01_A − B01_B > DIF_NUVEM_B01) E
//             A mais frio no B13 (B13_B − B13_A > DIF_NUVEM_B13).
//
// Esquerda : mine/ASTER_sazonal (A e B), TOA × 1000. O B original não tem SWIR: nas
//            origens 2 e 3 o SWIR fica vazio.
// Direita  : mine/ASTER_sazonal_hm_razao (A e B), superfície × 1000. B com SWIR estimado
//            de ASTER_MOSAIC_v2/ASTER_MOSAIC_sazonal_seco_hm_razao_2009_2020_<LxCy>
//            (estimar_swir_rf_v2.py) nas células que já têm; no resto, B sem SWIR.
//
// Tela dividida com os mapas ligados; arraste a divisória. Seletor de composição no topo
// ("Origem do pixel" mostra de onde veio cada pixel).
// As faixas da direita são as da esquerda × ganho mediano TOA → superfície
// (etapa 4 de harmonizar_sazonal_landsat.py), para as cores ficarem comparáveis.
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
var DIF_NUVEM_B01 = 50;   // × 1000 (0,05)
var DIF_NUVEM_B13 = 60;   // DN (~3 °C)

var BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
              'B10', 'B11', 'B12', 'B13', 'B14'];
var SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09'];

// garante as 14 bandas (SWIR que falta entra todo mascarado), Int32
function completar(img) {
  img = ee.Image(img);
  var vazio = ee.Image.constant(SWIR.map(function () { return 0; })).rename(SWIR)
                .toInt32().updateMask(0);
  var comSwir = ee.Image(ee.Algorithms.If(img.bandNames().contains('B04'),
                                          img.select(BANDAS),
                                          img.addBands(vazio).select(BANDAS)));
  return comSwir.toInt32();
}

function periodo(col, p) {
  return completar(ee.ImageCollection(col).filter(ee.Filter.eq('periodo', p)).mosaic());
}

// mesma regra de montar_nuvem; 'origem' = 1, 2 ou 3
function integrar(a, b) {
  var nuvemA = a.select('B01').subtract(b.select('B01')).gt(DIF_NUVEM_B01)
                .and(b.select('B13').subtract(a.select('B13')).gt(DIF_NUVEM_B13))
                .unmask(0);
  function camada(img, origem) {
    return img.addBands(ee.Image.constant(origem).toInt32().rename('origem'))
              .updateMask(img.select('B01').mask());
  }
  // mosaic(): a última fica por cima
  return ee.ImageCollection([camada(b, 2),
                             camada(b.updateMask(nuvemA), 3),
                             camada(a.updateMask(nuvemA.not()), 1)]).mosaic();
}

var COL_ORIG = MINE + 'ASTER_sazonal';
var COL_HM = MINE + 'ASTER_sazonal_hm_razao';
var SWIR_EST = ee.ImageCollection(MINE + 'ASTER_MOSAIC_v2')
  .filter(ee.Filter.stringStartsWith('system:index', 'ASTER_MOSAIC_sazonal_seco_hm_razao_2009_2020_'))
  .map(completar);

// Lado esquerdo: duas camadas, ligar/desligar na lista de camadas.
//   v1   = mosaico final v1 (mine/aster-mosaico-final, build_mosaico_final.py:
//          zero-cloud + semestrais v4, TOA × 1000, origem 1–4) — por cima, ligado
//   seco = sazonal seco original integrado acima (A + B, TOA × 1000, origem 1–3)
// Desligue o v1 para ver o seco por baixo (ou mude V1_LIGADO).
var V1_LIGADO = true;
var V1 = ee.ImageCollection(MINE + 'aster-mosaico-final').mosaic();
var SECO = integrar(periodo(COL_ORIG, 'A'), periodo(COL_ORIG, 'B'));
// B harmonizado: SWIR estimado por cima onde já existe
var B_HM = ee.ImageCollection([periodo(COL_HM, 'B')]).merge(SWIR_EST).mosaic();
var HARM = integrar(periodo(COL_HM, 'A'), B_HM);
var grid = ee.FeatureCollection(MINE + 'grid_irece_5x5_solape100');

// Landsat de referência da harmonização (exportar_landsat_sazonal.py): mediana 32 dias dos
// meses secos, reflectância de superfície × 1000. A (2000–2008) por cima do B (2009–2020),
// como no ASTER. Fica ABAIXO do ASTER no mapa da direita: desligue o ASTER para vê-lo.
// Sempre em cor falsa nir-red-green (≈ ASTER 3N-2-1), faixas iguais às do harmonizado.
var LS_COL = ee.ImageCollection(MINE + 'Landsat_sazonal');
var LANDSAT = ee.ImageCollection([LS_COL.filter(ee.Filter.eq('periodo', 'B')).mosaic()
                                          .select(['nir', 'red', 'green']),
                                  LS_COL.filter(ee.Filter.eq('periodo', 'A')).mosaic()
                                          .select(['nir', 'red', 'green'])]).mosaic();

// faixas no original (TOA × 1000; TIR em DN)
var FAIXA = {
  B01: [69, 176], B02: [51, 204], B3N: [80, 330], B04: [167, 362], B05: [77, 232],
  B06: [75, 234], B07: [71, 228], B08: [45, 168], B09: [34, 115],
  B10: [1124, 1635], B12: [1332, 1931], B14: [1644, 2259]
};
// ganho mediano original → harmonizado (média de A e B); TIR sem correção
var GANHO = {B01: 0.69, B02: 0.82, B3N: 1.06, B04: 0.935, B05: 0.935, B06: 0.935,
             B07: 0.935, B08: 0.935, B09: 0.935};

var COMPOSICOES = {
  'VNIR cor falsa (3N-2-1)': ['B3N', 'B02', 'B01'],
  'SWIR (4-6-8)': ['B04', 'B06', 'B08'],
  'SWIR (4-7-9)': ['B04', 'B07', 'B09'],
  'TIR (14-12-10)': ['B14', 'B12', 'B10']
};

function vis(bandas, harm) {
  var g = function (b) { return harm && GANHO[b] ? GANHO[b] : 1; };
  return {bands: bandas,
          min: bandas.map(function (b) { return FAIXA[b][0] * g(b); }),
          max: bandas.map(function (b) { return FAIXA[b][1] * g(b); })};
}

var esq = ui.Map(), dir = ui.Map();
ui.Map.Linker([esq, dir]);
[esq, dir].forEach(function (m) {
  m.setOptions('SATELLITE');
  m.setControlVisibility({all: false, layerList: true, mapTypeControl: true,
                          zoomControl: true, scaleControl: true});
});
esq.add(ui.Label('Mosaico final v1 · seco original',
                 {fontWeight: 'bold', fontSize: '16px', padding: '4px 8px', position: 'bottom-left'}));
dir.add(ui.Label('Mosaico final seco — harmonizado com Landsat',
                 {fontWeight: 'bold', fontSize: '16px', padding: '4px 8px', position: 'bottom-right'}));

var ESTADOS = ee.FeatureCollection('projects/mapbiomas-arida/BR_UF_2025');
// vetores por cima do mosaico, nos dois lados (ligar/desligar na lista de camadas)
function vetores() {
  return [
    ui.Map.Layer(ESTADOS.style({color: 'ffff00', width: 1.5, fillColor: '00000000'}), {}, 'estados'),
    ui.Map.Layer(grid.style({color: 'ffffff', width: 1, fillColor: '00000000'}), {}, 'grid (células)'),
    ui.Map.Layer(ee.FeatureCollection([ee.Feature(grid.geometry().dissolve(100))])
                   .style({color: 'ff0000', width: 2.5, fillColor: '00000000'}), {}, 'limite da área')
  ];
}

var PAL_ORIGEM = ['1a9850', '4575b4', 'd73027'];   // 1 A · 2 B (falta A) · 3 B (nuvem no A)
// v1: 1 zero-cloud 00–08 · 2 zero-cloud 09–20 · 3 semestrais 00–08 · 4 semestrais 09–20
var PAL_V1 = ['1a9850', '4575b4', 'a6d96a', 'd73027'];

function camadaLandsat() {
  var f = vis(['B3N', 'B02', 'B01'], true);    // mesmas faixas do ASTER harmonizado
  return ui.Map.Layer(LANDSAT, {bands: ['nir', 'red', 'green'], min: f.min, max: f.max},
                      'Landsat referência (nir-red-green)');
}

function mostrar(nome) {
  if (nome === 'Origem do pixel') {
    var v = {min: 1, max: 3, palette: PAL_ORIGEM};
    esq.layers().reset([
                ui.Map.Layer(ee.Image.constant(1), {min: 0, max: 1}, 'base'),
                ui.Map.Layer(SECO.select('origem'), v, 'origem — seco original'),
                ui.Map.Layer(V1.select('origem'), {min: 1, max: 4, palette: PAL_V1},
                                     'origem — v1', V1_LIGADO)
            ].concat(vetores()));
    dir.layers().reset([
                ui.Map.Layer(ee.Image.constant(1), {min: 0, max: 1}, 'base'),
                camadaLandsat(),
                ui.Map.Layer(HARM.select('origem'), v, 'origem — harmonizado')
            ].concat(vetores()));
    return;
  }
  var b = COMPOSICOES[nome];
  esq.layers().reset([ui.Map.Layer(SECO, vis(b, false), 'seco original — ' + nome),
                      ui.Map.Layer(V1, vis(b, false), 'v1 — ' + nome, V1_LIGADO)].concat(vetores()));
  dir.layers().reset([camadaLandsat(),
                      ui.Map.Layer(HARM, vis(b, true), 'ASTER harmonizado — ' + nome)].concat(vetores()));
}

var sel = ui.Select({items: Object.keys(COMPOSICOES).concat(['Origem do pixel']), value: 'VNIR cor falsa (3N-2-1)',
                     onChange: mostrar, style: {position: 'top-left'}});
esq.add(sel);
esq.add(ui.Label('origem v1: verde escuro = zero-cloud 00–08 · azul = zero-cloud 09–20 · ' +
  'verde claro = semestrais 00–08 · vermelho = semestrais 09–20\n' +
  'origem seco: verde = A · azul = B (falta A) · vermelho = B (nuvem no A)',
  {fontSize: '11px', position: 'top-left', whiteSpace: 'pre-wrap', maxWidth: '320px'}));
dir.add(ui.Label('origem: verde = A · azul = B (falta A) · vermelho = B (nuvem no A)',
                 {fontSize: '11px', position: 'top-right'}));

ui.root.widgets().reset([ui.SplitPanel({firstPanel: esq, secondPanel: dir, wipe: true})]);
esq.setCenter(-41.8, -12.1, 9);
mostrar(sel.getValue());
