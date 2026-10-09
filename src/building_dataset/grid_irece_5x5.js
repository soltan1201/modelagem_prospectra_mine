// Grid 5 × 5 (25 células, com sobreposição opcional SOLAPE_M) sobre o retângulo do bloco Irecê (LIST_COORD de
// build_mosaicos_semestrais.py), para exportar os mosaicos semestrais em
// pedaços menores. O asset salvo aqui será a referência de recorte para
// todas as análises seguintes.
//
// Numeração das células: linha 1 = norte, coluna 1 = oeste;
//   id = (linha − 1) × 5 + coluna   → 1 no canto NW, 25 no canto SE
//   nome = 'L<linha>C<coluna>'       → ex.: L3C2
//
//   coluna →   1    2    3    4    5
//   linha 1    1    2    3    4    5     (norte)
//   linha 2    6    7    8    9   10
//   linha 3   11   12   13   14   15
//   linha 4   16   17   18   19   20
//   linha 5   21   22   23   24   25     (sul)
//
// Seções:
//   1. Parâmetros
//   2. Construção da grid
//   3. Visualização (clique numa célula mostra id/nome)
//   4. Exportação para asset

// ------- 1. PARÂMETROS -------
var OESTE = -44.460483519818084;
var LESTE = -39.203525512005584;
var SUL   = -14.790793773956535;
var NORTE = -9.416525285369607;

var N_COLUNAS = 5;
var N_LINHAS  = 5;

// Sobreposição entre células vizinhas (m). 0 → grid sem sobreposição.
// Cada borda INTERNA avança SOLAPE_M para dentro da vizinha (faixa comum de
// 2 × SOLAPE_M); as bordas externas ficam no limite do bloco.
var SOLAPE_M = 100;

var id_asset_grid = SOLAPE_M > 0
    ? 'projects/mapbiomas-arida/mine/grid_irece_5x5_solape' + SOLAPE_M
    : 'projects/mapbiomas-arida/mine/grid_irece_5x5';

// ------- 2. CONSTRUÇÃO DA GRID -------
// Coordenadas calculadas no cliente: bordas exatas, sem acúmulo de erro.
// geodesic false → lados seguem meridianos/paralelos, como o retângulo original.
var dx = (LESTE - OESTE) / N_COLUNAS;
var dy = (NORTE - SUL) / N_LINHAS;

// Sobreposição em graus. Longitude usa o cos da latitude mais ao sul (maior
// |lat|): garante pelo menos SOLAPE_M em todo o bloco.
var solape_lat = SOLAPE_M / 110574;
var solape_lon = SOLAPE_M / (111320 * Math.cos(Math.abs(SUL) * Math.PI / 180));

var celulas = [];
for (var linha = 1; linha <= N_LINHAS; linha++) {
    for (var coluna = 1; coluna <= N_COLUNAS; coluna++) {
        var x0 = OESTE + (coluna - 1) * dx;
        var x1 = (coluna === N_COLUNAS) ? LESTE : x0 + dx;
        var y1 = NORTE - (linha - 1) * dy;
        var y0 = (linha === N_LINHAS) ? SUL : y1 - dy;

        // expande só as bordas internas
        if (coluna > 1)         { x0 -= solape_lon; }
        if (coluna < N_COLUNAS) { x1 += solape_lon; }
        if (linha > 1)          { y1 += solape_lat; }
        if (linha < N_LINHAS)   { y0 -= solape_lat; }

        var geom = ee.Geometry.Rectangle([x0, y0, x1, y1], null, false);
        celulas.push(ee.Feature(geom, {
            id: (linha - 1) * N_COLUNAS + coluna,
            nome: 'L' + linha + 'C' + coluna,
            linha: linha,
            coluna: coluna,
            solape_m: SOLAPE_M,
            oeste: x0, leste: x1, sul: y0, norte: y1
        }));
    }
}

var grid = ee.FeatureCollection(celulas).map(function(f) {
    return f.set('area_km2', f.geometry().area(1).divide(1e6));
});

var area_total = ee.Geometry.Rectangle([OESTE, SUL, LESTE, NORTE], null, false);

print('Tamanho da célula (graus): ' + dx.toFixed(6) + ' × ' + dy.toFixed(6));
print('Grid (25 células):', grid);
print('Área total (km²):', area_total.area(1).divide(1e6));
print('Soma das células (km²) — maior que a total quando há sobreposição:',
      grid.aggregate_sum('area_km2'));
print('Sobreposição (m):', SOLAPE_M, '| asset:', id_asset_grid);

// ------- 3. VISUALIZAÇÃO -------
Map.setOptions('SATELLITE');
Map.centerObject(area_total, 7);

// Preenchimento em xadrez para distinguir células vizinhas
var xadrez = ee.Image().int().paint(
    grid.map(function(f) {
        return f.set('par', ee.Number(f.get('linha')).add(f.get('coluna')).mod(2));
    }), 'par');
Map.addLayer(xadrez, {min: 0, max: 1, palette: ['#2a78d6', '#e34948'], opacity: 0.25},
             'Células (xadrez)');
Map.addLayer(ee.Image().byte().paint(grid, 1, 2), {palette: ['#ffffff']}, 'Bordas da grid');
Map.addLayer(ee.Image().byte().paint(ee.FeatureCollection([ee.Feature(area_total)]), 1, 3),
             {palette: ['#ffcc00']}, 'Limite do bloco Irecê');
Map.addLayer(grid.map(function(f) { return ee.Feature(f.geometry().centroid(1)); }),
             {color: 'ffffff'}, 'Centro das células', false);

// Legenda com a numeração e informação da célula clicada
var info = ui.Label('Clique numa célula para ver id e nome');
var tabela = [];
for (var l = 1; l <= N_LINHAS; l++) {
    var nums = [];
    for (var c = 1; c <= N_COLUNAS; c++) {
        var n = (l - 1) * N_COLUNAS + c;
        nums.push((n < 10 ? ' ' : '') + n);
    }
    tabela.push(nums.join('  '));
}
var painel = ui.Panel({
    widgets: [ui.Label('Grid Irecê 5 × 5 — sobreposição ' + SOLAPE_M + ' m', {fontWeight: 'bold'}),
              ui.Label('N ↑   (linha 1 = norte, coluna 1 = oeste)', {fontSize: '11px'}),
              ui.Label(tabela.join('\n'), {whiteSpace: 'pre', fontFamily: 'monospace'}),
              info],
    style: {position: 'bottom-left', width: '300px'}
});
Map.add(painel);

Map.onClick(function(coords) {
    var ponto = ee.Geometry.Point([coords.lon, coords.lat]);
    grid.filterBounds(ponto).toList(1)
        .evaluate(function(lista) {
            if (!lista || lista.length === 0) {
                info.setValue('Fora da grid');
                return;
            }
            var d = lista[0].properties;
            info.setValue('Célula id ' + d.id + ' — ' + d.nome +
                          ' — ' + d.area_km2.toFixed(0) + ' km²');
        });
});

// ------- 4. EXPORTAÇÃO PARA ASSET -------
// Vai para a aba Tasks: clique em "Run" para salvar.
Export.table.toAsset({
    collection: grid,
    description: id_asset_grid.split('/').pop(),
    assetId: id_asset_grid
});

// Cópia opcional no Drive (Shapefile), para abrir em QGIS:
// Export.table.toDrive({collection: grid, description: 'grid_irece_5x5_shp',
//                       folder: 'ASTER_grid', fileFormat: 'SHP'});
