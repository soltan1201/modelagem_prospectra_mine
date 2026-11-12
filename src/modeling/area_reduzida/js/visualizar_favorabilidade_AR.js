// ============================================================================
// Favorabilidade de fosfato na área reduzida (AR) de Irecê — visualizador
//
// Fonte: exportar_favorabilidade_AR.py
//   favorabilidade_AR/<RF|GTB>_fosfato_AR_<ID_AR>_<conjunto>_fcl  banda
//     'probability' = p × 10000 (Int16, 0–10000). Aqui p volta para 0–1.
//     Treino: 100 px por ocorrência + fundo 3 × presença sorteado só nos
//     clusters de treino (fundo por cluster, _fcl) → p ordena a área; não é a
//     chance real de ocorrência (≈ 0,25 é o neutro com 1:3).
//   cluster_fundo_AR_<ID_AR>  bandas 'cluster' (0..K−1, k-means do notebook
//     selecao_variaveis_AR.ipynb, célula [10e]) e 'treino' (1 = cluster usado
//     como fundo de treino).
// Conjunto (CONJUNTO, abaixo): selecao (12 variáveis) ou aster_gama_mag_geo (todas as 42).
// Versões (VERSAO, abaixo; ID_AR = área + sufixo do mosaico, amostras_AR.py):
//   irece_shm  mosaico final sazonal harmonizado (aster-mosaico-final_sazonal_hm_razao),
//              12 variáveis de models/variaveis_AR_irece_shm.json   ← padrão
//   irece      mosaico v1 (aster-mosaico-final), 13 variáveis de variaveis_AR_irece.json
// Variáveis montadas aqui como em amostras_AR.imagem_variaveis(); ligar/desligar
//   no painel. Faixa de cores = percentis 2–98 na área (calculada ao abrir).
// Clique no mapa: p RF/GTB, cluster e valor das variáveis no ponto.
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
var ACERVO = MINE + 'acervo_cprm/';

// mesma ordem do JSON (ranking SHAP estável) — mude junto com variaveis_AR_<ID_AR>.json
var VERSOES = {
  irece_shm: {
    mosaico: 'aster-mosaico-final_sazonal_hm_razao',
    conjuntos: {
      selecao: ['dist_dobra_km', 'lito_nova_america', 'MAG_nT', 'eU_ppm', 'CT', 'B09', 'eTh_ppm',
                'U_K', 'mag_thg', 'phosphorite_rbd', 'ferric_ox_4_3', 'mag_res'],
      // todas as variáveis (20 ASTER + 8 gama + 3 mag + 11 geologia), ordem do JSON
      aster_gama_mag_geo: ['B01', 'B3N', 'B09', 'B10', 'SAVI', 'ferric_iron_2_1', 'ferric_ox_4_3', 'gossan_4_2',
                          'ferrous_sil_5_4', 'sericite_illite_smectite', 'kaolinite_7_5', 'dolomite', 'carbonate_8_6',
                          'carbonate_13_14', 'quartz_14_12', 'sultan_r', 'sultan_g', 'phosphorite_rbd', 'glcm_contrast',
                          'glcm_ent', 'eU_ppm', 'eTh_ppm', 'K_pct', 'CT', 'U_Th', 'U_K', 'Th_K', 'eU_res', 'MAG_nT',
                          'mag_thg', 'mag_res', 'lito_nova_america', 'lito_jussara', 'lito_gabriel', 'lito_irece',
                          'lito_bebedouro', 'lito_morro_chapeu', 'lito_caboclo', 'lito_tombador', 'lito_cobertura',
                          'dist_falha_km', 'dist_dobra_km']
    }
  },
  irece: {
    mosaico: 'aster-mosaico-final',
    conjuntos: {
      selecao: ['dist_dobra_km', 'gossan_4_2', 'lito_nova_america', 'B10', 'MAG_nT', 'quartz_14_12',
                'dist_falha_km', 'eU_ppm', 'CT', 'sultan_g', 'eTh_ppm', 'mag_thg', 'ferric_iron_2_1']
    }
  }
};
var VERSAO = 'irece_shm';
// conjunto = parte do nome do asset: favorabilidade_AR/<RF|GTB>_fosfato_AR_<VERSAO>_<CONJUNTO>_fcl
//   'selecao'             12 variáveis escolhidas (recomendado)
//   'aster_gama_mag_geo'  todas as 42 (só no irece_shm)
var CONJUNTO = 'aster_gama_mag_geo';

var CFG = VERSOES[VERSAO];
var NOME_MAPA = '_fosfato_AR_' + VERSAO + '_' + CONJUNTO + '_fcl';
var MAPAS = {
  RF:  ee.Image(MINE + 'favorabilidade_AR/RF' + NOME_MAPA),
  GTB: ee.Image(MINE + 'favorabilidade_AR/GTB' + NOME_MAPA)
};
var CLUSTER = ee.Image(MINE + 'cluster_fundo_AR_' + VERSAO);
var PONTOS = ee.FeatureCollection(MINE + 'pontos_referencia').filter(ee.Filter.eq('alvo', 'fosfato'));
var AREA = ee.Geometry.Rectangle([-42.31762892604435, -12.359936623152336, -41.2382221877631, -10.5375816299746]);
var MODELOS = ['RF', 'GTB'];
var COR = {RF: '#222222', GTB: '#17becf'};
var ESCALA_EST = 300;      // m — percentis, histogramas e área acima do limiar
var ESCALA_TOA = 1000;     // mosaico = reflectância × 1000

if (!CFG.conjuntos[CONJUNTO]) {
  throw new Error('Conjunto ' + CONJUNTO + ' não definido em VERSOES.' + VERSAO + '.conjuntos');
}
var SELECAO = CFG.conjuntos[CONJUNTO];
var GAMA = ['eU_ppm', 'eTh_ppm', 'K_pct', 'CT', 'U_Th', 'U_K', 'Th_K', 'eU_res'];
var MAG = ['MAG_nT', 'mag_thg', 'mag_res'];
// lito_<grupo> → prefixo do SIGLA_UNID (amostras_AR.LITO_GRUPOS)
var LITO = {
  lito_nova_america: 'NP2sn', lito_jussara: 'NP2sj', lito_gabriel: 'NP2sg', lito_irece: 'NP2sr',
  lito_bebedouro: 'NP2be', lito_morro_chapeu: 'M3N1mc', lito_caboclo: 'MPc', lito_tombador: 'MPt',
  lito_cobertura: 'Q'
};
var ESTRUT = {dist_falha_km: ['Falha'], dist_dobra_km: ['clinal', 'forme'], dist_lineamento_km: ['Lineamentos']};
function grupoDe(n) {
  if (GAMA.indexOf(n) >= 0) { return 'gama'; }
  if (MAG.indexOf(n) >= 0) { return 'mag'; }
  if (LITO[n] || ESTRUT[n]) { return 'geologia'; }
  return 'aster';
}
var GRUPO = {};
SELECAO.forEach(function (n) { GRUPO[n] = grupoDe(n); });
var DESCR = {
  dist_dobra_km: 'distância às dobras (km)', gossan_4_2: 'B04/B02 (gossan)',
  lito_nova_america: 'Fm. Nova América (NP2sn)', B09: 'SWIR B09', B10: 'TIR B10',
  MAG_nT: 'campo magnético (nT)', mag_res: 'MAG − média 5 km (nT)',
  quartz_14_12: 'B14/B12 (quartzo)', dist_falha_km: 'distância às falhas (km)', eU_ppm: 'eU (ppm)',
  CT: 'contagem total', sultan_g: 'B04/B01 (Sultan G)', eTh_ppm: 'eTh (ppm)', U_K: 'razão eU/K',
  mag_thg: 'gradiente horizontal (nT/km)', ferric_iron_2_1: 'B02/B01 (ferro férrico)',
  ferric_ox_4_3: 'B04/B3N (óxido férrico)', phosphorite_rbd: '(ε11+ε13)/(2·ε12) (fosforito)',
  B01: 'VNIR B01', B3N: 'VNIR B3N', SAVI: 'SAVI (B3N, B02)', ferrous_sil_5_4: 'B05/B04 (silicato ferroso)',
  sericite_illite_smectite: '(B05+B07)/B06 (sericita)', kaolinite_7_5: 'B07/B05 (caulinita)',
  dolomite: '(B06+B08)/B07 (dolomita)', carbonate_8_6: 'B08/B06 (carbonato)',
  carbonate_13_14: 'B13/B14 (carbonato TIR)', sultan_r: 'B04/B07 (Sultan R)',
  glcm_contrast: 'GLCM contraste B3N', glcm_ent: 'GLCM entropia B3N',
  K_pct: 'K (%)', U_Th: 'razão eU/eTh', Th_K: 'razão eTh/K', eU_res: 'eU − média 5 km (ppm)',
  lito_jussara: 'Fm. Salitre, Jussara (NP2sj)', lito_gabriel: 'Fm. Salitre, Gabriel (NP2sg)',
  lito_irece: 'Fm. Salitre, Irecê (NP2sr)', lito_bebedouro: 'Fm. Bebedouro (NP2be)',
  lito_morro_chapeu: 'Fm. Morro do Chapéu (M3N1mc)', lito_caboclo: 'Fm. Caboclo (MPc)',
  lito_tombador: 'Fm. Tombador (MPt)', lito_cobertura: 'coberturas (Q)',
  dist_lineamento_km: 'distância aos lineamentos (km)'
};

var PAL_P = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
var PAL_DIV = ['2166ac', '67a9cf', 'd1e5f0', 'f7f7f7', 'fddbc7', 'ef8a62', 'b2182b'];
var PAL_VAR = ['440154', '3b528b', '21918c', '5ec962', 'fde725'];
var PAL_DIST = ['fde725', '5ec962', '21918c', '3b528b', '440154'];   // perto = claro
var PAL_CL = ['e6194b', '3cb44b', 'ffe119', '4363d8', 'f58231', '911eb4', '46f0f0', 'f032e6', 'bcf60c', 'fabebe'];
var VIS_P   = {min: 0, max: 1, palette: PAL_P};
var VIS_DIF = {min: -0.5, max: 0.5, palette: PAL_DIV};
var COR_LIMIAR = ['555555', '17becf', 'ff00ff'];           // só RF, só GTB, os dois
var COR_CORTE = {RF: 'ff8c00', GTB: '17becf'};             // camada de cada modelo acima do corte

// ----------------------------------------------------------------------------
// Variáveis (iguais a amostras_AR.py; geofísica mascarada fora da cobertura)
// ----------------------------------------------------------------------------
// TIR: DN → radiância; emissividade por normalização (ε_max 0,96), como em
// models_ML_fosfato.emissividade_tir
var TIR_UCC = {B10: 6.822e-3, B11: 6.780e-3, B12: 6.590e-3, B13: 5.693e-3, B14: 5.225e-3};
var TIR_LAMBDA = {B10: 8.291, B11: 8.634, B12: 9.075, B13: 10.657, B14: 11.318};
var C1 = 1.191042e8, C2 = 1.4387752e4, EPS_MAX = 0.96;
function emissividade(img) {
  var L = {}, T = [];
  Object.keys(TIR_UCC).forEach(function (b) {
    L[b] = img.select(b).subtract(1).max(1).multiply(TIR_UCC[b]);
    T.push(L[b].expression('c2 / (lam * log(eps * c1 / (pow(lam, 5) * L) + 1))',
                           {c1: C1, c2: C2, lam: TIR_LAMBDA[b], eps: EPS_MAX, L: L[b]}));
  });
  var tmax = ee.Image.cat(T).reduce(ee.Reducer.max());
  var eps = {};
  Object.keys(TIR_UCC).forEach(function (b) {
    eps[b] = L[b].divide(tmax.expression('c1 / (pow(lam, 5) * (exp(c2 / (lam * t)) - 1))',
                                         {c1: C1, c2: C2, lam: TIR_LAMBDA[b], t: tmax}));
  });
  return eps;
}

function aster() {
  var img = ee.ImageCollection(MINE + CFG.mosaico).mosaic();
  var sat = ee.Image(ee.Algorithms.If(img.bandNames().contains('saturado'),
                                      img.select('saturado').eq(0), ee.Image(1)));
  img = img.updateMask(sat).toFloat();
  function razao(a, b) { return img.select(a).divide(img.select(b)); }
  // calculadas só as que a versão usa (o mosaico pode não ter todas as bandas)
  function soma2(a, b, c) { return img.select(a).add(img.select(c)).divide(img.select(b)); }   // (a + c) / b
  var f = {
    SAVI: function () {
      return img.expression('1.5 * (n - r) / (n + r + 0.5)',
                            {n: img.select('B3N').divide(ESCALA_TOA), r: img.select('B02').divide(ESCALA_TOA)});
    },
    ferrous_sil_5_4: function () { return razao('B05', 'B04'); },
    sericite_illite_smectite: function () { return soma2('B05', 'B06', 'B07'); },
    kaolinite_7_5: function () { return razao('B07', 'B05'); },
    dolomite: function () { return soma2('B06', 'B07', 'B08'); },
    carbonate_8_6: function () { return razao('B08', 'B06'); },
    carbonate_13_14: function () { return razao('B13', 'B14'); },
    sultan_r: function () { return razao('B04', 'B07'); },
    glcm_contrast: function () { return glcm().select('B3N_contrast'); },
    glcm_ent: function () { return glcm().select('B3N_ent'); },
    gossan_4_2: function () { return razao('B04', 'B02'); },
    quartz_14_12: function () { return razao('B14', 'B12'); },
    sultan_g: function () { return razao('B04', 'B01'); },
    ferric_iron_2_1: function () { return razao('B02', 'B01'); },
    ferric_ox_4_3: function () { return razao('B04', 'B3N'); },
    phosphorite_rbd: function () {
      var eps = emissividade(img);
      return eps.B11.add(eps.B13).divide(eps.B12.multiply(2));
    }
  };
  // GLCM da B3N como no calcular_features: refl. 0–0,35 → byte 0–255, janela 3 × 3
  function glcm() {
    return img.select('B3N').unitScale(0, 0.35 * ESCALA_TOA).multiply(255).clamp(0, 255).toByte()
              .glcmTexture({size: 1});
  }
  return ee.Image.cat(SELECAO.filter(function (n) { return GRUPO[n] === 'aster'; })
                             .map(function (n) {
                               return (/^B\d/.test(n) ? img.select(n) : f[n]()).toFloat().rename(n);
                             }));
}

function geofisica() {
  var gama = ee.Image(ACERVO + '1042_gama').select(GAMA.slice(0, 7));
  var mag = ee.Image(ACERVO + '1042_mag').select('MAG_nT');
  var pm = mag.projection(), pg = gama.projection();
  var euRes = gama.select('eU_ppm')
    .subtract(gama.select('eU_ppm').reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(5000, 'meters')))
    .reproject(pg).rename('eU_res');
  var grad = mag.gradient().reproject(pm);
  var thg = grad.select('x').hypot(grad.select('y')).multiply(1000).rename('mag_thg');
  // anomalia local: MAG − média num círculo de 5 km, na grade original (RAIO_RES_M)
  var res = mag.subtract(mag.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(5000, 'meters')))
               .reproject(pm).rename('mag_res');
  var geo = ee.Image.cat([gama, euRes, mag, thg, res]).toFloat();
  return geo.updateMask(geo.mask().reduce(ee.Reducer.min()));
}

function geologia() {
  var reg = AREA.buffer(20000);
  var lito = ee.FeatureCollection(ACERVO + 'lito_jacobina_250k').filterBounds(reg);
  var estr = ee.FeatureCollection(ACERVO + 'estrut_jacobina_250k').filterBounds(reg);
  function dist(trechos) {
    var f = ee.Filter.or.apply(null, trechos.map(function (t) { return ee.Filter.stringContains('TIPO', t); }));
    return estr.filter(f).distance(20000, 30).unmask(20000).divide(1000);
  }
  var partes = Object.keys(LITO).map(function (n) {
    return ee.Image(0).paint(lito.filter(ee.Filter.stringStartsWith('SIGLA_UNID', LITO[n])), 1).rename(n);
  });
  Object.keys(ESTRUT).forEach(function (n) { partes.push(dist(ESTRUT[n]).rename(n)); });
  return ee.Image.cat(partes).toFloat().clip(AREA);
}

var VARS = ee.Image.cat([aster(), geofisica(), geologia()]).select(SELECAO).clip(AREA);
var P = ee.Image.cat(MODELOS.map(function (m) { return MAPAS[m].divide(10000).rename(m); }));

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var painel = ui.Panel({style: {width: '380px', padding: '6px'}});
painel.add(ui.Label('Favorabilidade fosfato — AR ' + VERSAO, {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('RF × GTB, ' + SELECAO.length + ' variáveis (' + CONJUNTO + '), mosaico ' + CFG.mosaico + ', fundo por cluster (_fcl). p = probability / 10000; ' +
                    'ordena a área, não é a chance real de ocorrência.', {fontSize: '11px', color: '666666'}));
var lblInfo = ui.Label('', {fontSize: '11px', whiteSpace: 'pre-wrap'});
painel.add(lblInfo);

MAPAS.RF.toDictionary(['n_ocorrencias', 'n_presenca', 'n_fundo', 'n_fundo_elegivel', 'clusters_treino'])
  .evaluate(function (r) {
    if (!r) { return; }
    lblInfo.setValue('Treino: ' + r.n_ocorrencias + ' ocorrências, presença ' + r.n_presenca + ', fundo ' +
                     r.n_fundo + ' (de ' + r.n_fundo_elegivel + ' elegíveis, clusters ' + r.clusters_treino + ')');
  });

function secao(t) { painel.add(ui.Label(t, {fontWeight: 'bold', margin: '10px 0 2px 8px'})); }

// camadas fixas: chave → ui.Map.Layer
var camadas = {};
function chkCamada(chave, texto, ligada) {
  var c = ui.Checkbox(texto, ligada, function (v) { camadas[chave].setShown(v); }, false, {fontSize: '12px'});
  return c;
}

secao('Mapas');
var chkRF = chkCamada('RF', 'p RF', false);
var chkGTB = chkCamada('GTB', 'p GTB', true);
var chkDif = chkCamada('dif', 'diferença GTB − RF', false);
var chkLimRF = chkCamada('lim_RF', 'RF acima do corte', true);
var chkLimGTB = chkCamada('lim_GTB', 'GTB acima do corte', false);
var chkLim = chkCamada('limiar', 'concordância acima do corte de cada modelo', false);
var chkCl = chkCamada('cluster', 'clusters do fundo (k-means)', false);
var chkTr = chkCamada('treino', 'clusters de treino (fundo sorteado aqui)', false);
var chkPt = chkCamada('pontos', 'ocorrências de fosfato', true);
[chkRF, chkGTB, chkDif, chkLimRF, chkLimGTB, chkLim, chkCl, chkTr, chkPt].forEach(function (c) { painel.add(c); });

// um corte por modelo: o GTB tem p mais concentrado (pouca área acima do mesmo corte do RF)
painel.add(ui.Label('Limiar por modelo', {fontSize: '12px', margin: '6px 0 0 8px'}));
var sld = {};
MODELOS.forEach(function (m) {
  sld[m] = ui.Slider({min: 0, max: 1, step: 0.01, value: 0.5, style: {stretch: 'horizontal'}});
  painel.add(ui.Panel([ui.Label(m, {fontSize: '12px', color: COR[m], fontWeight: 'bold', width: '36px'}), sld[m]],
                      ui.Panel.Layout.flow('horizontal')));
});
var chkJuntos = ui.Checkbox('mover os dois juntos', false, null, false, {fontSize: '11px'});
painel.add(chkJuntos);
// corte do GTB que deixa a mesma fração da área que o RF deixa acima do corte dele
painel.add(ui.Button({label: 'GTB com a mesma área do RF', style: {stretch: 'horizontal'}, onClick: function () {
  chkJuntos.setValue(false);   // senão o RF iria junto
  var t = sld.RF.getValue();
  lblArea.setValue('procurando o corte do GTB…');
  P.select('RF').gte(t).reduceRegion({reducer: ee.Reducer.mean(), geometry: AREA, scale: ESCALA_EST,
                                      maxPixels: 1e9, bestEffort: true}).evaluate(function (r) {
    var pct = Math.min(99.9, Math.max(0.1, 100 * (1 - r.RF)));
    P.select('GTB').reduceRegion({reducer: ee.Reducer.percentile([pct]), geometry: AREA, scale: ESCALA_EST,
                                  maxPixels: 1e9, bestEffort: true}).evaluate(function (q) {
      var v = q[Object.keys(q)[0]];
      sld.GTB.setValue(Math.round(v * 100) / 100);   // dispara atualizarLimiar
    });
  });
}}));
var lblArea = ui.Label('', {fontSize: '11px', whiteSpace: 'pre-wrap'});
painel.add(lblArea);

secao('Variáveis do modelo (ordem do SHAP)');
var lblVars = ui.Label('calculando percentis 2–98…', {fontSize: '11px', color: '666666'});
painel.add(lblVars);
var painelVars = ui.Panel();
painel.add(painelVars);
var opacidade = ui.Slider({min: 0, max: 1, step: 0.05, value: 1, style: {stretch: 'horizontal'}});
opacidade.onChange(function (v) {
  SELECAO.forEach(function (n) { if (camadas[n]) { camadas[n].setOpacity(v); } });
});
painel.add(ui.Panel([ui.Label('opacidade', {fontSize: '11px'}), opacidade], ui.Panel.Layout.flow('horizontal')));
painel.add(ui.Button({label: 'Desligar todas as variáveis', style: {stretch: 'horizontal'}, onClick: function () {
  painelVars.widgets().forEach(function (w) { w.widgets().get(0).setValue(false); });
}}));

var saida = ui.Panel();
secao('Análises');
painel.add(ui.Button({label: 'Histogramas + área × limiar', onClick: histogramas, style: {stretch: 'horizontal'}}));
painel.add(ui.Button({label: 'p por cluster (média e área)', onClick: pPorCluster, style: {stretch: 'horizontal'}}));
painel.add(saida);

// ----------------------------------------------------------------------------
// Mapa
// ----------------------------------------------------------------------------
Map.addLayer(ee.FeatureCollection([ee.Feature(AREA)]).style({color: '000000', width: 2, fillColor: '00000000'}),
             {}, 'área AR');

function adicionar(chave, img, vis, nome, ligada) {
  camadas[chave] = ui.Map.Layer(img, vis, nome, ligada);
  Map.layers().add(camadas[chave]);
}

// variáveis embaixo dos mapas de p (camadas adicionadas depois ficam por cima)
var camVars = ui.Map.Layer(ee.Image(), {}, '(variáveis)', false);
Map.layers().add(camVars);

adicionar('cluster', CLUSTER.select('cluster'), {min: 0, max: 9, palette: PAL_CL}, 'cluster (0–9)', false);
adicionar('treino', CLUSTER.select('treino').selfMask(), {palette: ['ffffff']}, 'clusters de treino', false);
adicionar('RF', P.select('RF'), VIS_P, 'p RF', false);
adicionar('GTB', P.select('GTB'), VIS_P, 'p GTB', true);
adicionar('dif', P.select('GTB').subtract(P.select('RF')), VIS_DIF, 'diferença GTB − RF', false);
adicionar('limiar', ee.Image(), {}, 'p ≥ limiar', false);
adicionar('lim_GTB', ee.Image(), {}, 'GTB ≥ corte', false);
adicionar('lim_RF', ee.Image(), {}, 'RF ≥ corte', true);
adicionar('pontos', PONTOS.filterBounds(AREA).map(function (f) {
  return f.set('estilo', {color: '000000', pointShape: 'star5', pointSize: 8, width: 1, fillColor: 'ff0000'});
}).style({styleProperty: 'estilo'}), {}, 'ocorrências de fosfato', true);

function atualizarLimiar() {
  var t = {RF: sld.RF.getValue(), GTB: sld.GTB.getValue()};
  var acima = ee.Image.cat([P.select('RF').gte(t.RF), P.select('GTB').gte(t.GTB)]);
  var txt = t.RF === t.GTB ? 'p ≥ ' + t.RF : 'RF ≥ ' + t.RF + ' · GTB ≥ ' + t.GTB;
  var conc = acima.select('RF').add(acima.select('GTB').multiply(2)).selfMask();   // 1 só RF, 2 só GTB, 3 os dois
  camadas.limiar.setEeObject(conc);
  camadas.limiar.setVisParams({min: 1, max: 3, palette: COR_LIMIAR});
  camadas.limiar.setName(txt);
  MODELOS.forEach(function (m) {
    camadas['lim_' + m].setEeObject(acima.select(m).selfMask());
    camadas['lim_' + m].setVisParams({palette: [COR_CORTE[m]]});
    camadas['lim_' + m].setName(m + ' ≥ ' + t[m]);
  });
  lblLegLimiar.setValue(txt);
  lblArea.setValue('calculando a área acima do corte…');
  acima.addBands(acima.select('RF').and(acima.select('GTB')).rename('ambos'))
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: AREA, scale: ESCALA_EST, maxPixels: 1e9, bestEffort: true})
    .evaluate(function (r) {
      if (!r) { return; }
      lblArea.setValue('Área acima do corte (' + txt + '):\nRF ' + (100 * r.RF).toFixed(1) + '%  ·  GTB ' +
                       (100 * r.GTB).toFixed(1) + '%  ·  os dois ' + (100 * r.ambos).toFixed(1) + '%');
    });
}
var sincronizando = false;
MODELOS.forEach(function (m) {
  sld[m].onChange(function (v) {
    if (chkJuntos.getValue() && !sincronizando) {
      sincronizando = true;
      MODELOS.forEach(function (o) { if (o !== m) { sld[o].setValue(v); } });
      sincronizando = false;
    }
    atualizarLimiar();
  });
});

// uma camada por variável, criada depois dos percentis (faixa de cores)
var visVar = {};
VARS.reduceRegion({reducer: ee.Reducer.percentile([2, 98]), geometry: AREA, scale: ESCALA_EST,
                   maxPixels: 1e9, bestEffort: true}).evaluate(function (r) {
  lblVars.setValue('Faixa de cores = percentis 2–98 na área (' + ESCALA_EST + ' m). Geofísica só na cobertura.');
  var pos = Map.layers().indexOf(camVars);
  Map.layers().remove(camVars);
  SELECAO.forEach(function (n, k) {
    var vis;
    if (n.indexOf('lito_') === 0) {
      vis = {min: 0, max: 1, palette: ['ffffff', '8c510a']};
    } else {
      var lo = r[n + '_p2'], hi = r[n + '_p98'];
      if (lo === null || hi === null) { lo = 0; hi = 1; }
      vis = {min: Number(lo.toPrecision(3)), max: Number(hi.toPrecision(3)),
             palette: n.indexOf('dist_') === 0 ? PAL_DIST : PAL_VAR};
    }
    visVar[n] = vis;
    camadas[n] = ui.Map.Layer(VARS.select(n), vis, (k + 1) + '. ' + n, false, opacidade.getValue());
    Map.layers().insert(pos + k, camadas[n]);
    var chk = ui.Checkbox((k + 1) + '. ' + n, false, function (v) {
      camadas[n].setShown(v);
      if (v) { mostrarLegendaVar(n); }
    }, false, {fontSize: '12px', margin: '2px 4px'});
    var faixa = n.indexOf('lito_') === 0 ? '0 / 1' : vis.min + ' – ' + vis.max;
    painelVars.add(ui.Panel([chk, ui.Label(GRUPO[n] + ' · ' + DESCR[n] + ' · ' + faixa,
                                           {fontSize: '10px', color: '666666', margin: '5px 0'})],
                            ui.Panel.Layout.flow('horizontal', true)));
  });
});

// ----------------------------------------------------------------------------
// Análises
// ----------------------------------------------------------------------------
function histogramas() {
  saida.clear();
  saida.add(ui.Chart.image.histogram({image: P, region: AREA, scale: ESCALA_EST, minBucketWidth: 0.02, maxPixels: 1e9})
    .setOptions({title: 'Histograma de p (frequência em log)', colors: [COR.RF, COR.GTB],
                 hAxis: {title: 'p', viewWindow: {min: 0, max: 1}}, vAxis: {title: 'pixels', logScale: true}}));
  P.reduceRegion({reducer: ee.Reducer.fixedHistogram(0, 1, 50), geometry: AREA, scale: ESCALA_EST,
                  maxPixels: 1e9, bestEffort: true}).evaluate(function (h) {
    var acum = {};
    MODELOS.forEach(function (m) {
      var cont = h[m].map(function (x) { return x[1]; });
      var total = cont.reduce(function (a, b) { return a + b; }, 0), resto = total;
      acum[m] = cont.map(function (c) { var v = resto / total; resto -= c; return v; });
    });
    var linhas = [['limiar p', 'RF', 'GTB']];
    for (var i = 0; i < 50; i++) { linhas.push([h.RF[i][0], 100 * acum.RF[i], 100 * acum.GTB[i]]); }
    saida.add(ui.Chart(linhas, 'LineChart', {title: 'Área com p ≥ limiar (%)', colors: [COR.RF, COR.GTB],
      lineWidth: 2, hAxis: {title: 'limiar p'}, vAxis: {title: '% da área', viewWindow: {min: 0, max: 100}}}));
  });
}

function pPorCluster() {
  saida.clear();
  saida.add(ui.Label('calculando…', {fontSize: '11px'}));
  // RF e GTB: média e contagem; banda 2 = cluster
  var red = ee.Reducer.mean().combine(ee.Reducer.count(), '', true).repeat(2)
              .group({groupField: 2, groupName: 'cluster'});
  ee.Dictionary({
    g: P.addBands(CLUSTER.select('cluster')).reduceRegion({reducer: red, geometry: AREA, scale: ESCALA_EST,
                                                           maxPixels: 1e9, bestEffort: true}).get('groups'),
    tr: CLUSTER.get('clusters_treino'), d: CLUSTER.get('dist_fosfato'), raio: CLUSTER.get('raio')
  }).evaluate(function (r) {
    saida.clear();
    var treino = String(r.tr).split(',').map(Number), dist = String(r.d).split(',');
    var dt = [['cluster', 'dist. ao fosfato', 'treino', 'pixels (300 m)', 'p RF', 'p GTB']];
    r.g.forEach(function (x) {
      dt.push([x.cluster, Number(dist[x.cluster]), treino.indexOf(x.cluster) >= 0 ? 'sim' : '',
               x.count[0], x.mean[0], x.mean[1]]);
    });
    saida.add(ui.Label('raio do cluster de fosfato = ' + Number(r.raio).toFixed(2) +
                       ' (espaço padronizado). Treino = sem ocorrência e distância > raio.',
                       {fontSize: '11px', whiteSpace: 'pre-wrap'}));
    saida.add(ui.Chart(dt, 'Table', {}));
    saida.add(ui.Label('p médio baixo nos clusters de treino é esperado (são o fundo do modelo); ' +
                       'o que interessa é p alto fora deles.', {fontSize: '11px', whiteSpace: 'pre-wrap'}));
  });
}

// ----------------------------------------------------------------------------
// Legenda
// ----------------------------------------------------------------------------
var FS_LEG = '11px';
function barraCores(vis) {
  return ui.Thumbnail({image: ee.Image.pixelLonLat().select(0),
    params: {bbox: [0, 0, 1, 0.1], dimensions: '150x10', format: 'png', min: 0, max: 1, palette: vis.palette},
    style: {stretch: 'horizontal', margin: '0 8px', maxHeight: '12px'}});
}
function linhaEscala(vis) {
  return ui.Panel([ui.Label(String(vis.min), {fontSize: '10px', margin: '0 8px'}),
                   ui.Label('', {stretch: 'horizontal'}),
                   ui.Label(String(vis.max), {fontSize: '10px', margin: '0 8px'})], ui.Panel.Layout.flow('horizontal'));
}
function itemCor(cor, texto, forma) {
  var simbolo = forma === 'estrela'
    ? ui.Label('★', {color: '#' + cor, fontSize: '15px', margin: '0 6px 0 4px', padding: '0'})
    : ui.Label('', {backgroundColor: '#' + cor, padding: '7px', margin: '2px 6px 2px 4px', border: '1px solid #888888'});
  return ui.Panel([simbolo, ui.Label(texto, {fontSize: FS_LEG, margin: '3px 0'})], ui.Panel.Layout.flow('horizontal'));
}
function titulo(t) { return ui.Label(t, {fontSize: FS_LEG, fontWeight: 'bold', margin: '6px 0 2px 4px'}); }

var lblLegLimiar = titulo('p ≥ —');
var legenda = ui.Panel({style: {position: 'bottom-left', padding: '6px 8px', width: '210px'}});
legenda.add(ui.Label('Favorabilidade AR ' + VERSAO, {fontSize: '13px', fontWeight: 'bold', margin: '0 0 2px 4px'}));
legenda.add(titulo('p (RF / GTB)'));
legenda.add(barraCores(VIS_P));
legenda.add(linhaEscala(VIS_P));
legenda.add(titulo('Diferença GTB − RF'));
legenda.add(barraCores(VIS_DIF));
legenda.add(linhaEscala(VIS_DIF));
legenda.add(lblLegLimiar);
legenda.add(itemCor(COR_CORTE.RF, 'RF acima do corte'));
legenda.add(itemCor(COR_CORTE.GTB, 'GTB acima do corte'));
legenda.add(itemCor(COR_LIMIAR[0], 'só RF'));
legenda.add(itemCor(COR_LIMIAR[1], 'só GTB'));
legenda.add(itemCor(COR_LIMIAR[2], 'RF e GTB'));
legenda.add(itemCor('ff0000', 'ocorrência de fosfato', 'estrela'));
var legVar = ui.Panel();
legenda.add(legVar);
Map.add(legenda);

function mostrarLegendaVar(n) {
  legVar.clear();
  legVar.add(titulo(n));
  if (n.indexOf('lito_') === 0) {
    legVar.add(itemCor('8c510a', LITO[n]));
    legVar.add(itemCor('ffffff', 'outras unidades'));
    return;
  }
  legVar.add(barraCores(visVar[n]));
  legVar.add(linhaEscala(visVar[n]));
}

chkCl.onChange(function (v) {
  legVar.clear();
  if (!v) { return; }
  legVar.add(titulo('Cluster (* = treino)'));
  CLUSTER.get('clusters_treino').evaluate(function (tr) {
    var t = String(tr).split(',').map(Number);
    PAL_CL.forEach(function (c, k) { legVar.add(itemCor(c, 'cluster ' + k + (t.indexOf(k) >= 0 ? ' *' : ''))); });
  });
});

var chkLegenda = ui.Checkbox('Mostrar legenda no mapa', true);
chkLegenda.onChange(function (v) { legenda.style().set('shown', v); });
painel.add(chkLegenda);

// ----------------------------------------------------------------------------
// Clique
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '260px'}});
info.add(ui.Label('Clique no mapa para ver p e as variáveis.', {fontSize: '11px'}));
Map.add(info);
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  info.clear();
  info.add(ui.Label('carregando…', {fontSize: '11px'}));
  P.addBands(CLUSTER).addBands(VARS).reduceRegion(ee.Reducer.first(), ee.Geometry.Point(c.lon, c.lat), 30)
    .evaluate(function (r) {
      info.clear();
      r = r || {};
      function f(x, d) { return (x === null || x === undefined) ? '—' : x.toFixed(d); }
      info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontSize: '11px', fontWeight: 'bold'}));
      info.add(ui.Label('p RF ' + f(r.RF, 3) + '  ·  p GTB ' + f(r.GTB, 3), {fontSize: '12px'}));
      info.add(ui.Label('cluster ' + f(r.cluster, 0) + (r.treino === 1 ? ' (treino)' : ''), {fontSize: '11px'}));
      info.add(ui.Label(SELECAO.map(function (n) { return n + ': ' + f(r[n], 3); }).join('\n'),
                        {fontSize: '10px', whiteSpace: 'pre'}));
    });
});

ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(AREA, 9);
atualizarLimiar();
