// ============================================================================
// SCRIPT GEE - MAPEAMENTO ESPECTRAL E GEOFÍSICO - BACIA DE IRECÊ
// VERIFICA NÚMERO DE BANDAS PARA ESCOLHER PALETTE OU BANDS
// ============================================================================

var asset = 'projects/ee-solkancengine17/assets/prospectra/referencias';
var mapas_col = ee.ImageCollection(asset);

var asset_canarana = 'projects/ee-solkancengine17/assets/prospectra/canarana_gologia_misi_souto';
var feat_canarana = ee.FeatureCollection(asset_canarana);

var asset_geologia = 'projects/ee-solkancengine17/assets/prospectra/geologia_2004';
var feat_geologia = ee.FeatureCollection(asset_geologia);

// ==================== CONFIGURAÇÕES DE PALETAS (1 BANDA) ====================

var visParamsSingle = {
    'potassio': {min:0, max:2.5, palette: ['0000ff','00ffff','ffff00','ff0000'], name:'Potássio (%)'},
    'torio': {min:0, max:20, palette: ['440154','3b528b','21908c','5dc963','fde725'], name:'Tório (ppm)'},
    'uranio': {min:0, max:5, palette: ['440154','31688e','35b779','dce319','fde725'], name:'Urânio (ppm)'},
    'U_Th': {min:0, max:1.2, palette: ['2166ac','67a9cf','f7f7f7','ef8a62','b2182b'], name:'Razão U/Th'},
    'Th_K': {min:0, max:15, palette: ['762a83','9970ab','e7d4e8','faf0f0','1b7837'], name:'Razão Th/K'},
    'U_K': {min:0, max:3, palette: ['1a9641','a6d96a','ffffbf','fdae61','d7191c'], name:'Razão U/K'},
    'abrams': {min:-0.5, max:0.5, palette: ['313695','4575b4','74add1','fee090','f46d43','d73027','a50026'], name:'Índice Abrams'},
    'sultan': {min:-0.8, max:0.8, palette: ['053061','2166ac','92c5de','f7f7f7','f4a582','d6604d','b2182b'], name:'Índice Sultan'},
    'rf_class': {min:1, max:15, palette: ['a6cee3','1f78b4','b2df8a','33a02c','fb9a99','e31a23','fdbf6f','ff7f00','cab2d6','6a3d9a','ffff99','b15928','8dd3c7','ffffb3','fccde5'], name:'Classificação RF'},
    'rf_class_filtro': {min:1, max:15, palette: ['a6cee3','1f78b4','b2df8a','33a02c','fb9a99','e31a23','fdbf6f','ff7f00','cab2d6','6a3d9a','ffff99','b15928','8dd3c7','ffffb3','fccde5'], name:'RF com filtro'},
};

// Função cliente para obter parâmetros de uma banda (baseada no nome)
function getSingleParams(idStr) {
    var lower = idStr.toLowerCase();
    if (lower.indexOf('u_th') !== -1 || lower.indexOf('u/th') !== -1) return visParamsSingle['U_Th'];
    if (lower.indexOf('th_k') !== -1 || lower.indexOf('th/k') !== -1) return visParamsSingle['Th_K'];
    if (lower.indexOf('u_k') !== -1 || lower.indexOf('u/k') !== -1) return visParamsSingle['U_K'];
    if (lower.indexOf('ternario') !== -1 || lower.indexOf('ternary') !== -1) return null; // RGB
    if (lower.indexOf('abrams') !== -1) return visParamsSingle['abrams'];
    if (lower.indexOf('sultan') !== -1) return visParamsSingle['sultan'];
    if (lower.indexOf('rf') !== -1 || lower.indexOf('random_forest') !== -1) return visParamsSingle['rf_class'];
    if (lower.indexOf('filtro') !== -1) return visParamsSingle['rf_class_filtro'];
    if ((lower.indexOf('k') !== -1 || lower.indexOf('potassio') !== -1) && 
        lower.indexOf('th_k') === -1 && lower.indexOf('u_k') === -1) return visParamsSingle['potassio'];
    if ((lower.indexOf('th') !== -1 || lower.indexOf('torio') !== -1) && 
        lower.indexOf('u_th') === -1) return visParamsSingle['torio'];
    if ((lower.indexOf('u') !== -1 || lower.indexOf('uranio') !== -1) && 
        lower.indexOf('u_th') === -1 && lower.indexOf('u_k') === -1) return visParamsSingle['uranio'];
    return {min:0, max:255, palette: ['black','white'], name:idStr};
}

// ==================== ADICIONAR CAMADAS DE IMAGEM ====================

Map.centerObject(feat_geologia, 9);

var lst_map = mapas_col.reduceColumns(ee.Reducer.toList(), ['system:index']).get('list');
ee.List(lst_map).evaluate(function(lstIds) {
    for (var i = 0; i < lstIds.length; i++) {
        var id_img = lstIds[i];
        var idStr = String(id_img);
        var img = mapas_col.filter(ee.Filter.eq('system:index', id_img)).first();
        
        // Verifica o número de bandas da imagem (getInfo síncrono dentro do evaluate)
        var nBands = img.bandNames().size().getInfo();
        
        if (nBands > 1) {
        // Imagem com múltiplas bandas (RGB, ternário, etc.)
        var bandNames = img.bandNames().getInfo();
        // Usa as três primeiras bandas ou todas se menos de 3
        var bandsToUse = bandNames.slice(0, Math.min(3, bandNames.length));
        Map.addLayer(img, {min:0, max:255, bands: bandsToUse}, idStr);
        } else {
        // Imagem de uma banda - aplica paleta correspondente
        var params = getSingleParams(idStr);
        if (params) {
            Map.addLayer(img, {min: params.min, max: params.max, palette: params.palette}, params.name);
        } else {
            Map.addLayer(img, {min:0, max:255, palette: ['black','white']}, idStr);
        }
        }
    }
});

// ==================== ESTILIZAÇÃO DAS FEIÇÕES GEOLÓGICAS ====================

// Dicionário de cores (apenas para referência - usaremos uma função de estilo direta)
var unidadeColors = {
    'Unidade A': '#00aa00',
    'Unidade A1': '#55ff55',
    'Unidade B': '#ffaa00',
    'Unidade B1': '#ff8800',
    'Unidade C': '#aa6633',
    'Grupo Chapada Diamantina': '#888888',
    'Formação Caatinga': '#dddd99',
    'Cobertura TQ': '#ffff99'
};
var defaultColor = '#cccccc';

// Função de estilo que será aplicada diretamente no Map.addLayer
function estiloGeologia(feature) {
    // Tenta obter o nome da unidade de diferentes campos (ajuste conforme seus dados)
    var nome = feature.get('nome_unidade');
    if (!nome) nome = feature.get('unidade');
    if (!nome) nome = feature.get('classe');
    if (!nome) nome = 'default';
    
    var color = unidadeColors[nome] || defaultColor;
    return {color: color, fillColor: color, opacity: 0.7, width: 1.5};
}

// Adicionar camadas usando a função de estilo diretamente
Map.addLayer(feat_geologia, {style: estiloGeologia}, 'Geologia UFBA 2004');
Map.addLayer(feat_canarana, {style: estiloGeologia}, 'Geologia Misi & Souto');
print('Script executado. Verifica número de bandas automaticamente.');