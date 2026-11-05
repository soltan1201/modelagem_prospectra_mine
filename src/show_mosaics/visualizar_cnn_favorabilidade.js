// ============================================================================
// Favorabilidade da CNN (fosfato) — visualizador de todos os modelos
//
// Fonte: projects/mapbiomas-arida/mine/favorabilidade_ml, sem RF/GTB (antes: cnn_favorabilidade)
//   - upload à mão (juntar_mapas_gee.py + propriedades_mapas_gee.py): uma imagem
//     por rodada e passo, banda 'b1';
//   - bucket + ingestão (subir_predicao_gee.py): uma imagem por célula, banda
//     'probability'.
//   Nos dois casos o valor é p × 10000 (0 = sem dado) e as propriedades são
//   rodada, passo, preenchido (0/1), modelo (CNN_<backbone>), canais (A/B1), patch,
//   folds (media5 = média dos 5 fold*.pt, 1fold, final = final.pt da rodada _teste).
//   Aqui cada rodada + passo + preenchido + folds vira UM mapa (mosaic das imagens).
//
// Fundo: mosaico v1 (aster-mosaico-final) + os 15 índices espectrais do modelo A,
//   em grupos selecionáveis no painel (as mesmas fórmulas do
//   maxent_fosfato.INDICES / exportar_patches).
//
// Painel
//   - modelos encontrados na coleção (backbone, canais, folds, passo, bordas
//     preenchidas), com uma caixa para cada um entrar ou sair do consenso e da
//     % de área; filtros rápidos por canais (A/B1) e folds
//   - canais de entrada dos modelos A e B1
//   - índices espectrais, um seletor por grupo
//   - limiar: valor de p (0–10000, o mesmo para todos) ou top X% de cada modelo;
//     pixels acima do limiar em cada modelo, quantos dos marcados concordam e a %
//     da área acima do limiar em cada um
// Clique no mapa: p de todos os modelos no ponto.visualizar_cnn_favorabilidade
// ============================================================================

var MINE = 'projects/mapbiomas-arida/mine/';
// mapas da CNN subidos à mão no favorabilidade_ml, ao lado do RF/GTB
var COL = ee.ImageCollection(MINE + 'favorabilidade_ml')
          .filter(ee.Filter.neq( 'modelo', "RF"))
          .filter(ee.Filter.neq( 'modelo', "GTB"));
var MOSAICO_V1 = ee.ImageCollection(MINE + 'aster-mosaico-final');
var GRID = ee.FeatureCollection(MINE + 'grid_irece_5x5');
var PONTOS = ee.FeatureCollection(MINE + 'pontos_referencia').filter(ee.Filter.eq('alvo', 'fosfato'));
var regiao = GRID.geometry().dissolve(100);
var ESCALA_AREA = 300;      // m — % da área acima do limiar e estiramento dos índices
var ESCALA_TOA = 1000;      // mosaico = reflectância × 1000

// canais de entrada de cada modelo (exportar_patches.MODELOS)
var CANAIS = {
  A: ['ferric_iron_2_1', 'ferric_ox_4_3', 'gossan_4_2', 'ferrous_sil_5_4', 'SAVI',
      'sericite_illite_smectite', 'kaolinite_7_5', 'muscovite_7_6',
      'carbonate_8_6', 'dolomite', 'carbonate_13_14',
      'quartz_14_12', 'phosphorite_rbd', 'sultan_r', 'sultan_g', 'swir_estimado'],
  B1: ['log_B01', 'log_B02', 'log_B3N', 'log_B04', 'log_B05', 'log_B06', 'log_B07', 'log_B08', 'log_B09',
       'emis_B10', 'emis_B11', 'emis_B12', 'emis_B13', 'emis_B14', 'swir_estimado']
};
// propriedade folds: que modelos fizeram o mapa (propriedades_mapas_gee.py)
var ROTULO_FOLDS = {media5: 'média dos 5 folds', '1fold': '1 fold', final: 'final.pt (_teste)'};
var DESCRICAO = {
  A: '15 índices espectrais + swir_estimado (16 canais)',
  B1: 'log da reflectância B01–B09 + emissividade B10–B14 + swir_estimado (15 canais)'
};

var PAL_P = ['000004', '3b0f70', '8c2981', 'de4968', 'fe9f6d', 'fcfdbf'];
var VIS_P = {min: 0, max: 10000, palette: PAL_P};
// uma cor por modelo na camada do limiar (até 12; depois repete)
var CORES = ['e41a1c', '377eb8', '4daf4a', '984ea3', 'ff7f00', 'ffff33',
             'a65628', 'f781bf', '17becf', 'bcbd22', '00ff00', 'ff00ff'];
var PAL_CONSENSO = ['ffffb2', 'fecc5c', 'fd8d3c', 'f03b20', 'bd0026'];

// ----------------------------------------------------------------------------
// Índices espectrais no mosaico v1 (maxent_fosfato.INDICES)
// ----------------------------------------------------------------------------
var v1 = MOSAICO_V1.mosaic();
function razao(a, b) { return v1.select(a).divide(v1.select(b)); }

// emissividade (normalização, ε_max 0,96) para o phosphorite_rbd
var TIR_UCC = {B10: 6.822e-3, B11: 6.780e-3, B12: 6.590e-3, B13: 5.693e-3, B14: 5.225e-3};
var TIR_LAMBDA = {B10: 8.291, B11: 8.634, B12: 9.075, B13: 10.657, B14: 11.318};
var C1 = 1.191042e8, C2 = 1.4387752e4, EPS_MAX = 0.96;
function emissividade() {
  var L = {}, T = [];
  Object.keys(TIR_UCC).forEach(function (b) {
    L[b] = v1.select(b).subtract(1).max(1).multiply(TIR_UCC[b]);
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
var eps = emissividade();

// os 15 índices do modelo A, por grupo; paleta: baixo → alto (estiramento p2–p98 na área)
var PAL_IDX = ['2c7bb6', 'abd9e9', 'ffffbf', 'fdae61', 'd7191c'];
var GRUPOS = [
  {nome: 'Vegetação', indices: [
    ['SAVI', v1.expression('1.5 * (n - r) / (n + r + 0.5)',
                           {n: v1.select('B3N').divide(ESCALA_TOA), r: v1.select('B02').divide(ESCALA_TOA)}),
     'SAVI — vegetação (B3N, B02)']]},
  {nome: 'Óxidos de ferro (VNIR/SWIR)', indices: [
    ['ferric_iron_2_1', razao('B02', 'B01'), 'ferro férrico (B02/B01)'],
    ['ferric_ox_4_3', razao('B04', 'B3N'), 'óxido férrico (B04/B3N)'],
    ['gossan_4_2', razao('B04', 'B02'), 'gossan (B04/B02)'],
    ['ferrous_sil_5_4', razao('B05', 'B04'), 'silicato ferroso (B05/B04)']]},
  {nome: 'Argilominerais e micas (SWIR)', indices: [
    ['sericite_illite_smectite', v1.select('B05').add(v1.select('B07')).divide(v1.select('B06')),
     'sericita/ilita/esmectita ((B05+B07)/B06)'],
    ['kaolinite_7_5', razao('B07', 'B05'), 'caulinita (B07/B05)'],
    ['muscovite_7_6', razao('B07', 'B06'), 'muscovita (B07/B06)']]},
  {nome: 'Carbonatos (SWIR/TIR)', indices: [
    ['carbonate_8_6', razao('B08', 'B06'), 'carbonato SWIR (B08/B06)'],
    ['dolomite', v1.select('B06').add(v1.select('B08')).divide(v1.select('B07')), 'dolomita ((B06+B08)/B07)'],
    ['carbonate_13_14', razao('B13', 'B14'), 'carbonato TIR (B13/B14)']]},
  {nome: 'Sílica e fosfato (TIR)', indices: [
    ['quartz_14_12', razao('B14', 'B12'), 'quartzo (B14/B12)'],
    ['phosphorite_rbd', eps.B11.add(eps.B13).divide(eps.B12.multiply(2)), 'fosforito (ε11+ε13)/(2·ε12)']]},
  {nome: 'Sultan (litologia)', indices: [
    ['sultan_r', razao('B04', 'B07'), 'Sultan R (B04/B07)'],
    ['sultan_g', razao('B04', 'B01'), 'Sultan G (B04/B01)']]}
];
var INDICES = [];
GRUPOS.forEach(function (g) { INDICES = INDICES.concat(g.indices); });
var imgIndices = ee.Image.cat(INDICES.map(function (x) { return x[1].toFloat().rename(x[0]); }));
var sultan = ee.Image.cat([razao('B04', 'B07').rename('R'), razao('B04', 'B01').rename('G'),
                           razao('B02', 'B3N').multiply(razao('B04', 'B3N')).rename('B')]);

// ----------------------------------------------------------------------------
// Painel
// ----------------------------------------------------------------------------
var painel = ui.Panel({style: {width: '400px', padding: '6px'}});
var FS = '11px';
painel.add(ui.Label('Favorabilidade da CNN — fosfato', {fontWeight: 'bold', fontSize: '16px'}));
painel.add(ui.Label('Valor = p × 10000 (0 = sem dado). Pseudo-rótulos do RF/GTB: p ordena a área, ' +
                    'não é a chance real de ocorrência.', {fontSize: FS, color: '666666'}));

painel.add(ui.Label('Modelos na coleção', {fontWeight: 'bold', margin: '8px 8px 2px 8px'}));
painel.add(ui.Label('Marcados = entram no consenso e na % de área. Os filtros marcam de uma vez.',
                    {fontSize: FS, color: '666666'}));
var selCanais = ui.Select({items: ['todos', 'A', 'B1'], value: 'todos', style: {margin: '0 4px 0 8px'}});
var selFolds = ui.Select({items: [{label: 'todos', value: 'todos'}].concat(Object.keys(ROTULO_FOLDS).map(function (k) {
  return {label: ROTULO_FOLDS[k], value: k}; })), value: 'todos', style: {margin: '0 8px 0 4px'}});
painel.add(ui.Panel([ui.Label('Canais', {fontSize: FS, margin: '8px 0 0 8px'}), selCanais,
                     ui.Label('Folds', {fontSize: FS, margin: '8px 0 0 4px'}), selFolds],
                    ui.Panel.Layout.flow('horizontal')));
var pnlModelos = ui.Panel();
pnlModelos.add(ui.Label('carregando…', {fontSize: FS}));
painel.add(pnlModelos);

painel.add(ui.Label('Canais de entrada', {fontWeight: 'bold', margin: '8px 8px 2px 8px'}));
['A', 'B1'].forEach(function (m) {
  painel.add(ui.Label('Modelo ' + m + ': ' + DESCRICAO[m], {fontSize: FS, fontWeight: 'bold', margin: '2px 8px'}));
  painel.add(ui.Label(CANAIS[m].join(', '), {fontSize: FS, whiteSpace: 'pre-wrap', margin: '0 8px 4px 8px'}));
});

painel.add(ui.Label('Índices espectrais (mosaico v1)', {fontWeight: 'bold', margin: '8px 8px 2px 8px'}));
painel.add(ui.Label('Um seletor por grupo: o índice escolhido vai para o mapa (por cima do fundo, abaixo dos ' +
                    'modelos); "— nenhum —" tira. Cores: p2–p98 de cada índice na área.', {fontSize: FS, color: '666666'}));
var pnlIndices = ui.Panel();
painel.add(pnlIndices);

painel.add(ui.Label('Limiar', {fontWeight: 'bold', margin: '8px 8px 2px 8px'}));
painel.add(ui.Label('Valor de p: o mesmo corte (0–10000) para todos; perto de 10000 = só os pixels de maior p ' +
                    '(o máximo dos mapas atuais é ~9000–9350 e o 1% mais alto começa em ~6000–6900).\n' +
                    'Top %: cada modelo fica com os seus X% de maior p na área — compara modelos com escalas ' +
                    'diferentes.', {fontSize: FS, color: '666666', whiteSpace: 'pre-wrap'}));
var MODOS = {valor: 'valor de p (0–10000)', top: 'top % de cada modelo'};
var selModo = ui.Select({items: [{label: MODOS.valor, value: 'valor'}, {label: MODOS.top, value: 'top'}],
                         value: 'valor', style: {stretch: 'horizontal', margin: '0 8px'}});
painel.add(selModo);
var sldLimiar = ui.Slider({min: 0, max: 10000, step: 50, value: 6000, style: {stretch: 'horizontal'}});
painel.add(sldLimiar);
var lblArea = ui.Label('', {fontSize: FS, whiteSpace: 'pre-wrap'});
painel.add(lblArea);

// ----------------------------------------------------------------------------
// Modelos: agrupa as imagens por rodada + passo + preenchido
// ----------------------------------------------------------------------------
var mapas = [];        // [{nome, rodada, passo, preenchido, backbone, canais, n, img}]
var camLimiar = [], camConsenso = null;

function nomeMapa(m) {
  return m.rodada + ' · passo ' + m.passo + (m.preenchido ? ' · preenchido' : '') + ' · ' + m.folds;
}
function ativos() { return mapas.filter(function (m) { return m.ativo; }); }

ee.Dictionary({
  n: COL.size(),
  com_rodada: COL.filter(ee.Filter.notNull(['rodada', 'passo'])).size()
}).evaluate(function (r, erro) {
  pnlModelos.clear();
  if (erro || !r || r.n === 0) {
    pnlModelos.add(ui.Label(erro ? 'Erro: ' + erro : 'Coleção vazia: suba os mapas (juntar_mapas_gee.py + upload) ' +
                            'e grave as propriedades (propriedades_mapas_gee.py).', {fontSize: FS, color: 'b2182b'}));
    desenharFundo();
    return;
  }
  if (r.com_rodada !== r.n) {
    pnlModelos.add(ui.Label((r.n - r.com_rodada) + ' imagem(ns) sem as propriedades rodada/passo (upload sem ' +
                            'propriedades_mapas_gee.py) — ficam de fora.', {fontSize: FS, color: 'b2182b'}));
  }
  var grupos = {};
  // toList descarta a linha inteira se faltar uma coluna: padrão para as opcionais
  COL.filter(ee.Filter.notNull(['rodada', 'passo'])).map(function (im) {
    return im.set(ee.Dictionary({preenchido: -1, modelo: '', canais: '', folds: ''})
                    .combine(im.toDictionary(['rodada', 'passo', 'preenchido', 'modelo', 'canais', 'folds']), true));
  }).reduceColumns(
    ee.Reducer.toList(6), ['rodada', 'passo', 'preenchido', 'modelo', 'canais', 'folds']
  ).get('list').evaluate(function (linhas) {
    (linhas || []).forEach(function (l) {
      var chave = l[0] + '|' + l[1] + '|' + l[2] + '|' + l[5];
      if (!grupos[chave]) {
        grupos[chave] = {rodada: l[0], passo: l[1], preenchido: Math.max(l[2], 0), temPreenchido: l[2] !== -1,
                         backbone: String(l[3] || '').replace('CNN_', ''), canais: l[4] || '?', folds: l[5] || '?',
                         temFolds: l[5] !== '', n: 0, ativo: true};
      }
      grupos[chave].n += 1;
    });
    mapas = Object.keys(grupos).sort().map(function (k, i) {
      var m = grupos[k];
      m.nome = nomeMapa(m);
      m.cor = CORES[i % CORES.length];
      m.banda = 'm' + i;
      var sel = COL.filter(ee.Filter.eq('rodada', m.rodada)).filter(ee.Filter.eq('passo', m.passo));
      if (m.temPreenchido) { sel = sel.filter(ee.Filter.eq('preenchido', m.preenchido)); }
      if (m.temFolds) { sel = sel.filter(ee.Filter.eq('folds', m.folds)); }
      m.img = sel.map(function (im) { return im.select([0], ['p']); })
                 .mosaic().selfMask().rename(m.banda);
      return m;
    });
    listarModelos();
    desenhar();
  });
});

function listarModelos() {
  pnlModelos.clear();
  mapas.forEach(function (m) {
    var chk = ui.Checkbox({value: m.ativo, style: {margin: '3px 0 3px 8px'}});
    chk.onChange(function (v) { m.ativo = v; atualizarLimiar(); });
    pnlModelos.add(ui.Panel([
      chk,
      ui.Label('', {backgroundColor: '#' + m.cor, padding: '7px', margin: '3px 6px 3px 2px', border: '1px solid #888888'}),
      ui.Label(m.rodada + '  ·  ' + m.backbone + '  ·  ' + m.canais + '  ·  ' + (ROTULO_FOLDS[m.folds] || m.folds) +
               '  ·  passo ' + m.passo +
               ' (' + (30 * m.passo) + ' m)' + (m.preenchido ? '  ·  bordas preenchidas' : '') +
               '  ·  ' + m.n + (m.n === 1 ? ' imagem' : ' células'),
               {fontSize: FS, margin: '3px 0', width: '300px', whiteSpace: 'pre-wrap'})
    ], ui.Panel.Layout.flow('horizontal')));
  });
  pnlModelos.add(ui.Label(ativos().length + ' de ' + mapas.length + ' mapas marcados. A camada do limiar usa ' +
                          'estas cores; o consenso conta quantos dos marcados têm p ≥ limiar no pixel.',
                          {fontSize: FS, color: '666666'}));
}

function filtrar() {
  var c = selCanais.getValue(), f = selFolds.getValue();
  mapas.forEach(function (m) { m.ativo = (c === 'todos' || m.canais === c) && (f === 'todos' || m.folds === f); });
  listarModelos();
  atualizarLimiar();
}
selCanais.onChange(filtrar);
selFolds.onChange(filtrar);

// ----------------------------------------------------------------------------
// Mapa
// ----------------------------------------------------------------------------
function desenharFundo() {
  Map.layers().reset();
  Map.addLayer(ee.Image.constant(1), {min: 0, max: 1}, 'base branca', false);
  Map.addLayer(v1, {bands: ['B3N', 'B02', 'B01'], min: [155, 50, 68], max: [332, 202, 176]}, 'mosaico v1 — B3N/B02/B01 (falsa cor)', false);
  Map.addLayer(v1, {bands: ['B04', 'B06', 'B08'], min: 55, max: 362}, 'mosaico v1 — SWIR B04/B06/B08', true);
  Map.addLayer(sultan, {min: [1.5, 1.5, 0.11], max: [2.46, 2.96, 1.2]}, 'Sultan RGB (4/7, 4/1, 2/3N×4/3N)', false);   // p2–p98 na área
  // os índices escolhidos voltam ao mapa (logo acima do fundo)
  camIndices = {};
  INDICES.forEach(function (x) { if (marcados[x[0]]) { ligarIndice(x, true); } });
}

// índices: camada criada quando escolhida no seletor do grupo, removida quando trocada
var camIndices = {}, marcados = {}, faixa = {};
var POS_INDICES = 4;     // depois das 4 camadas do fundo (base, 2 RGB do v1, Sultan RGB), antes dos modelos

function ligarIndice(x, ligar) {
  var nome = x[0];
  marcados[nome] = ligar;
  if (camIndices[nome]) { Map.layers().remove(camIndices[nome]); delete camIndices[nome]; }
  if (!ligar) { return; }
  var vis = faixa[nome] ? {min: faixa[nome][0], max: faixa[nome][1], palette: PAL_IDX} : {palette: PAL_IDX};
  camIndices[nome] = ui.Map.Layer(imgIndices.select(nome), vis, 'índice ' + x[2], true);
  Map.layers().insert(Math.min(POS_INDICES, Map.layers().length()), camIndices[nome]);
}

// um ui.Select por grupo; a faixa p2–p98 do índice escolhido aparece embaixo
var lblFaixa = {};
function textoFaixa(nome) {
  return nome && faixa[nome] ? nome + ': ' + faixa[nome][0].toFixed(3) + ' – ' + faixa[nome][1].toFixed(3) : '';
}
GRUPOS.forEach(function (g) {
  var sel = ui.Select({
    items: [{label: '— nenhum —', value: ''}].concat(g.indices.map(function (x) { return {label: x[2], value: x[0]}; })),
    value: '', style: {stretch: 'horizontal', margin: '0 8px'}
  });
  var lbl = ui.Label('', {fontSize: '10px', color: '666666', margin: '0 8px 4px 16px'});
  lblFaixa[g.nome] = {lbl: lbl, sel: sel};
  sel.onChange(function (nome) {
    g.indices.forEach(function (x) { if (marcados[x[0]] && x[0] !== nome) { ligarIndice(x, false); } });
    g.indices.forEach(function (x) { if (x[0] === nome) { ligarIndice(x, true); } });
    lbl.setValue(textoFaixa(nome));
  });
  pnlIndices.add(ui.Label(g.nome, {fontSize: FS, fontWeight: 'bold', margin: '4px 8px 0 8px'}));
  pnlIndices.add(sel);
  pnlIndices.add(lbl);
});
function listarIndices() {
  GRUPOS.forEach(function (g) { var o = lblFaixa[g.nome]; o.lbl.setValue(textoFaixa(o.sel.getValue())); });
}

// estiramento p2–p98 de cada índice na área, uma vez (as camadas já ligadas são atualizadas)
imgIndices.reduceRegion({reducer: ee.Reducer.percentile([2, 98]), geometry: regiao, scale: ESCALA_AREA * 2,
                         maxPixels: 1e9, bestEffort: true})
  .evaluate(function (r) {
    if (!r) { return; }
    INDICES.forEach(function (x) {
      var lo = r[x[0] + '_p2'], hi = r[x[0] + '_p98'];
      if (lo !== null && lo !== undefined && hi !== null && hi !== undefined) {
        faixa[x[0]] = [lo, hi];
        if (camIndices[x[0]]) { camIndices[x[0]].setVisParams({min: lo, max: hi, palette: PAL_IDX}); }
      }
    });
    listarIndices();
  });

function desenhar() {
  desenharFundo();
  mapas.forEach(function (m, i) {
    Map.addLayer(m.img, VIS_P, 'p ' + m.nome, i === 0);
  });
  Map.addLayer(ee.FeatureCollection([ee.Feature(regiao)]).style({color: '000000', width: 2, fillColor: '00000000'}),
               {}, 'contorno da área');
  Map.addLayer(ee.Image().paint(GRID, 0, 1), {palette: 'ffffff'}, 'grid 5×5', false);
  Map.addLayer(PONTOS.map(function (f) {
    return f.set('estilo', {color: '000000', pointShape: 'star5', pointSize: 7, width: 1,
                            fillColor: ee.Algorithms.If(ee.String(f.get('papel')).equals('treino'), 'ff0000', '00ffff')});
  }).style({styleProperty: 'estilo'}), {}, 'ocorrências de fosfato (vermelho treino, ciano validação)');
  camLimiar = [];
  camConsenso = null;
  atualizarLimiar();
}

var pedido = 0;     // descarta respostas de limiares antigos (o slider dispara várias vezes)

function atualizarLimiar() {
  if (!mapas.length) { return; }
  camLimiar.forEach(function (c) { Map.layers().remove(c); });
  if (camConsenso) { Map.layers().remove(camConsenso); camConsenso = null; }
  camLimiar = [];
  var lista = ativos();
  if (!lista.length) {
    lblArea.setValue('Nenhum modelo marcado.');
    lblConsenso.setValue('Consenso (nenhum modelo marcado)');
    return;
  }
  var v = sldLimiar.getValue(), meu = ++pedido;
  if (selModo.getValue() === 'valor') {
    var lims = {};
    lista.forEach(function (m) { lims[m.banda] = v; });
    desenharLimiar(lista, lims, 'p ≥ ' + v);
    return;
  }
  // top X%: o percentil (100 − X) de cada modelo na área
  lblArea.setValue('calculando o percentil ' + (100 - v) + ' de cada modelo…');
  ee.Image.cat(lista.map(function (m) { return m.img; }))
    .reduceRegion({reducer: ee.Reducer.percentile([100 - v], ['lim']), geometry: regiao, scale: ESCALA_AREA,
                   maxPixels: 1e9, bestEffort: true})
    .evaluate(function (r) {
      if (!r || meu !== pedido) { return; }
      var lims = {};
      lista.forEach(function (m) {     // uma banda só: o GEE devolve o nome sem sufixo
        lims[m.banda] = r[m.banda + '_lim'] !== undefined ? r[m.banda + '_lim'] : r[m.banda];
      });
      desenharLimiar(lista, lims, 'top ' + v + '%');
    });
}

function desenharLimiar(lista, lims, rotulo) {
  camLimiar = lista.map(function (m) {
    var c = ui.Map.Layer(m.img.gte(lims[m.banda]).selfMask(), {palette: [m.cor]},
                         rotulo + ' — ' + m.nome + ' (p ≥ ' + Math.round(lims[m.banda]) + ')', false);
    Map.layers().add(c);
    return c;
  });
  var acima = ee.Image.cat(lista.map(function (m) { return m.img.unmask(0).gte(lims[m.banda]); }));
  var n = lista.length;
  camConsenso = ui.Map.Layer(acima.reduce(ee.Reducer.sum()).selfMask(),
                             {min: 1, max: Math.max(n, 2), palette: PAL_CONSENSO},
                             'consenso: nº de modelos com ' + rotulo + ' (1 a ' + n + ')', true);
  Map.layers().add(camConsenso);
  lblConsenso.setValue('Consenso (' + rotulo + '): 1 → ' + n + ' modelos');
  lblArea.setValue('calculando a área (' + rotulo + ')…');
  // % da área (só onde o modelo tem dado)
  var validos = ee.Image.cat(lista.map(function (m) { return m.img.mask().rename(m.banda); }));
  var meu = pedido;
  acima.updateMask(validos).reduceRegion({reducer: ee.Reducer.mean(), geometry: regiao, scale: ESCALA_AREA,
                                          maxPixels: 1e9, bestEffort: true})
    .evaluate(function (r) {
      if (!r || meu !== pedido) { return; }
      lblArea.setValue('Área com ' + rotulo + ':\n' + lista.map(function (m) {
        var a = r[m.banda];
        return '  ' + m.nome + ': ' + (a === null || a === undefined ? '—' : (100 * a).toFixed(2) + '%') +
               '  (p ≥ ' + Math.round(lims[m.banda]) + ')';
      }).join('\n'));
    });
}
selModo.onChange(function (modo) {
  if (modo === 'top') { sldLimiar.setMin(0.1); sldLimiar.setMax(20); sldLimiar.setStep(0.1); sldLimiar.setValue(1, false); }
  else { sldLimiar.setMin(0); sldLimiar.setMax(10000); sldLimiar.setStep(50); sldLimiar.setValue(6000, false); }
  atualizarLimiar();
});
sldLimiar.onChange(atualizarLimiar);

// ----------------------------------------------------------------------------
// Legenda
// ----------------------------------------------------------------------------
function barraCores(pal) {
  return ui.Thumbnail({image: ee.Image.pixelLonLat().select(0),
                       params: {bbox: [0, 0, 1, 0.1], dimensions: '150x10', format: 'png', min: 0, max: 1, palette: pal},
                       style: {stretch: 'horizontal', margin: '0 8px', maxHeight: '12px'}});
}
function escala(a, b) {
  return ui.Panel([ui.Label(a, {fontSize: '10px', margin: '0 8px'}), ui.Label('', {stretch: 'horizontal'}),
                   ui.Label(b, {fontSize: '10px', margin: '0 8px'})], ui.Panel.Layout.flow('horizontal'));
}
var lblConsenso = ui.Label('Consenso', {fontSize: FS, fontWeight: 'bold', margin: '6px 0 2px 4px'});
var legenda = ui.Panel({style: {position: 'bottom-left', padding: '6px 8px', width: '210px'}});
legenda.add(ui.Label('CNN — fosfato', {fontSize: '13px', fontWeight: 'bold', margin: '0 0 2px 4px'}));
legenda.add(ui.Label('p × 10000', {fontSize: FS, fontWeight: 'bold', margin: '6px 0 2px 4px'}));
legenda.add(barraCores(PAL_P));
legenda.add(escala('0', '10000'));
legenda.add(lblConsenso);
legenda.add(barraCores(PAL_CONSENSO));
legenda.add(escala('poucos', 'todos'));
legenda.add(ui.Label('Índices: baixo → alto (p2–p98)', {fontSize: FS, fontWeight: 'bold', margin: '6px 0 2px 4px'}));
legenda.add(barraCores(PAL_IDX));
Map.add(legenda);
var chkLegenda = ui.Checkbox('Mostrar legenda no mapa', true);
chkLegenda.onChange(function (v) { legenda.style().set('shown', v); });
painel.add(chkLegenda);

// ----------------------------------------------------------------------------
// Clique: p de todos os modelos no ponto
// ----------------------------------------------------------------------------
var info = ui.Panel({style: {position: 'bottom-right', width: '300px'}});
info.add(ui.Label('Clique no mapa para ver p de cada modelo.', {fontSize: FS}));
Map.add(info);
Map.style().set('cursor', 'crosshair');
Map.onClick(function (c) {
  if (!mapas.length) { return; }
  info.clear();
  info.add(ui.Label('carregando…', {fontSize: FS}));
  ee.Image.cat(mapas.map(function (m) { return m.img; })).addBands(imgIndices)
    .reduceRegion(ee.Reducer.first(), ee.Geometry.Point(c.lon, c.lat), 30)
    .evaluate(function (r) {
      info.clear();
      r = r || {};
      info.add(ui.Label('Lon ' + c.lon.toFixed(4) + '  Lat ' + c.lat.toFixed(4), {fontSize: FS, fontWeight: 'bold'}));
      mapas.forEach(function (m) {
        var v = r[m.banda];
        info.add(ui.Label(m.nome + ': ' + (v === null || v === undefined ? '—' : v), {fontSize: FS, color: '#' + m.cor}));
      });
      var ph = r.phosphorite_rbd;
      info.add(ui.Label('phosphorite_rbd: ' + (ph === null || ph === undefined ? '—' : ph.toFixed(4)),
                        {fontSize: FS, color: '666666'}));
    });
});

ui.root.insert(0, painel);
Map.setOptions('HYBRID');
Map.centerObject(ee.Geometry.Point(-41.86, -11.30), 8);
