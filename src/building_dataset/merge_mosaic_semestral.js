// ------- 1. PARÂMETROS -------
var list_coord = [
    [-44.460483519818084,-14.790793773956535],
    [-39.203525512005584,-14.790793773956535],
    [-39.203525512005584,-9.416525285369607],
    [-44.460483519818084,-9.416525285369607],
    [-44.460483519818084,-14.790793773956535]
];
var area_estudo = ee.Geometry.Polygon(list_coord);

var id_asset    = 'projects/mapbiomas-arida/mosaic_aster';
var id_asset_p2 = 'projects/mapbiomas-arida/mosaic_aster_p2008';
var nome_saida  = 'MOSAICO_FINAL_ASTER_IRECE';

// ------- 2. CARREGAR COLEÇÕES -------
var colecao = ee.ImageCollection(id_asset);
print('=== MOSAICO 1 (mosaic_aster) ===');
print('Total de mosaicos semestrais:', colecao.size());
print('Bandas disponíveis:', colecao.first().bandNames());
print('Propriedades:', colecao.first().propertyNames());

var colecao_p2 = ee.ImageCollection(id_asset_p2);
print('=== MOSAICO 2 (mosaic_aster_p2008) ===');
print('Total de mosaicos semestrais:', colecao_p2.size());
print('Bandas disponíveis:', colecao_p2.first().bandNames());
print('Propriedades:', colecao_p2.first().propertyNames());

// ------- 3. QUALITY MOSAIC USANDO BANDA EXISTENTE -------
var mosaico = colecao.qualityMosaic('quality').clip(area_estudo);
print('Mosaico 1 criado:', mosaico);
print('Bandas do mosaico 1:', mosaico.bandNames());

var mosaico_p2 = colecao_p2.qualityMosaic('quality').clip(area_estudo);
print('Mosaico 2 criado:', mosaico_p2);
print('Bandas do mosaico 2:', mosaico_p2.bandNames());

// ------- 4. VISUALIZAÇÃO -------
// Assets salvos em INT16 com fator de escala 10000 (exceto 'quality')
var ESCALA = 10000;

var vis_rgb    = mosaico.divide(ESCALA).float();
var vis_rgb_p2 = mosaico_p2.divide(ESCALA).float();

Map.centerObject(area_estudo, 8);
Map.addLayer(area_estudo, {color: 'red'}, 'Área de Estudo', false);

// --- Mosaico 1 (mosaic_aster) ---
Map.addLayer(vis_rgb, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35, gamma: 1.4
}, 'M1 - RGB Natural');

Map.addLayer(vis_rgb, {
    bands: ['B04', 'B06', 'B3N'],
    min: 0.1, max: 0.5, gamma: 1.2
}, 'M1 - SWIR Falsa-cor');

Map.addLayer(mosaico.select('quality'), {
    min: 0, max: 20000,
    palette: ['red', 'yellow', 'green'],
    opacity: 0.6
}, 'M1 - Quality Band');

// --- Mosaico 2 (mosaic_aster_p2008) ---
Map.addLayer(vis_rgb_p2, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35, gamma: 1.4
}, 'M2 - RGB Natural');

Map.addLayer(vis_rgb_p2, {
    bands: ['B04', 'B06', 'B3N'],
    min: 0.1, max: 0.5, gamma: 1.2
}, 'M2 - SWIR Falsa-cor');

Map.addLayer(mosaico_p2.select('quality'), {
    min: 0, max: 20000,
    palette: ['red', 'yellow', 'green'],
    opacity: 0.6
}, 'M2 - Quality Band');

// ------- 5. EXPORTAR PARA ASSET -------
Export.image.toAsset({
    image: mosaico,
    description: nome_saida,
    assetId: 'projects/mapbiomas-arida/' + nome_saida,
    region: area_estudo,
    scale: 15,
    maxPixels: 1e13,
    pyramidingPolicy: {
        '.default': 'sample',
        'quality': 'max'
    }
});
