// Download das FeatureCollections de pontos amostrados (asset GEE → CSV).
// Versão Code Editor de download_points_assets.py.
//
// Cria uma task Export.table.toDrive por tabela da pasta id_asset_pontos
// (um CSV por asset, com o mesmo nome). As tasks aparecem na aba Tasks
// e precisam ser iniciadas com "Run".

// ------- 1. PARÂMETROS -------
var id_asset_pontos = 'projects/mapbiomas-arida/mine/points';
var pasta_drive     = 'ASTER_Samples';
var incluir_geo     = false;   // true → mantém a coluna .geo no CSV

// ------- 2. LISTAR TABELAS DA PASTA -------
var assets = ee.data.listAssets(id_asset_pontos).assets || [];
var ids = assets
    .filter(function(a) { return a.type === 'TABLE'; })
    .map(function(a) { return a.id; });
print('Tabelas encontradas:', ids.length, ids);

// ------- 3. UMA TASK DE EXPORT POR TABELA -------
ids.forEach(function(id) {
    var nome = id.split('/').pop();
    var fc   = ee.FeatureCollection(id);

    // Colunas exportadas (sem propriedades system:*)
    var colunas = fc.first().propertyNames()
        .filter(ee.Filter.stringStartsWith('item', 'system:').not())
        .getInfo();
    if (incluir_geo) {
        colunas.push('.geo');
    }

    Export.table.toDrive({
        collection: fc,
        description: nome,
        folder: pasta_drive,
        fileNamePrefix: nome,
        fileFormat: 'CSV',
        selectors: colunas
    });
    print('Task criada:', nome);
});
