// =================================================================
// --- ANÁLISE DE HISTOGRAMAS: CENAS ASTER L1T "BOAS" vs "RUINS" ---
// =================================================================
// Compara pares de cenas (mesmo ano) onde a falsa-cor (B3N, B02, B01)
// fica "boa" numa e quase toda vermelha na outra.
//
// Para cada cena:
//   - metadados de aquisição (data, elevação solar, nuvem, GANHO por banda);
//   - histogramas por banda em DN bruto e em reflectância TOA;
//   - histograma sobreposto bom × ruim da banda B3N (NIR), por par;
//   - percentis (P2, P50, P98) por banda, em DN e TOA;
//   - camadas no mapa: falsa-cor em DN e em TOA.
//
// No fim: tabela-resumo única das 4 cenas (metadados, ganho, percentis DN/TOA,
// fração saturada DN ≥ 254, área comum do par) no Console e em CSV no Drive.
//
// Conversão usada (ASTER/AST_L1T_003 vem em DN, não em radiância):
//   radiância = (DN − 1) × GAIN_COEFFICIENT_Bxx   (DN = 0 → sem dado)
//   TOA       = π · L · d² / (ESUN · cos θz)

// ------- 1. PARÂMETROS -------
var pares = [
    {nome: '2003', bom: 'ASTER/AST_L1T_003/20030914131811', ruim: 'ASTER/AST_L1T_003/20031011130046'},
    {nome: '2001', bom: 'ASTER/AST_L1T_003/20010816131957', ruim: 'ASTER/AST_L1T_003/20010928130033'}
];

var BANDAS_VNIR = ['B01', 'B02', 'B3N'];
var ESCALA      = 60;      // m — histograma/percentis (15 m nativo fica lento numa cena inteira)
var BANDA_PAR   = 'B3N';   // banda do histograma sobreposto bom × ruim (NIR = "vermelho" da falsa-cor)
var DN_SAT      = 254;     // DN ≥ 254 = saturado (VNIR 8 bits)
var PASTA_DRIVE = 'ASTER_Histogramas';   // destino do CSV-resumo (Export.table.toDrive)

// Irradiância solar exoatmosférica (mesmos valores do radianceToTOA do pipeline)
var ESUN = {'B01': 1848, 'B02': 1549, 'B3N': 1114};

// ------- 2. FUNÇÕES -------

// DN → radiância → reflectância TOA, só VNIR
function dnParaTOA(img) {
    var elev = ee.Number(img.get('SOLAR_ELEVATION'));
    var cosz = ee.Number(90).subtract(elev).multiply(Math.PI / 180).cos();
    var doy  = ee.Date(img.get('system:time_start')).getRelative('day', 'year').add(1);
    // distância Terra-Sol (UA)
    var d    = ee.Number(1).subtract(
                   ee.Number(0.01672).multiply(doy.subtract(4).multiply(0.9856 * Math.PI / 180).cos()));

    var bandas = BANDAS_VNIR.map(function(b) {
        var dn    = img.select(b);
        var ganho = ee.Number(img.get('GAIN_COEFFICIENT_' + b));
        var rad   = dn.subtract(1).multiply(ganho);
        return rad.multiply(Math.PI).multiply(d.pow(2))
                  .divide(cosz.multiply(ESUN[b]))
                  .updateMask(dn.gt(0))
                  .rename(b);
    });
    return ee.Image.cat(bandas);
}

// Percentis P2/P50/P98 das bandas, na pegada da cena
function percentis(img, regiao) {
    return img.reduceRegion({
        reducer: ee.Reducer.percentile([2, 50, 98]),
        geometry: regiao, scale: ESCALA, maxPixels: 1e10, bestEffort: true
    });
}

function histograma(img, regiao, titulo, eixoX) {
    return ui.Chart.image.histogram({image: img, region: regiao, scale: ESCALA, maxPixels: 1e10})
        .setSeriesNames(BANDAS_VNIR)
        .setOptions({
            title: titulo,
            hAxis: {title: eixoX},
            vAxis: {title: 'Contagem de pixels'},
            colors: ['#2a78d6', '#e34948', '#7a4fbf'],   // B01, B02, B3N
            lineWidth: 2, pointSize: 0
        });
}

// Metadados que explicam diferença de brilho entre cenas (dicionário plano:
// data, sol, nuvem + todas as propriedades com "GAIN" no nome)
function metadados(img) {
    var chavesGanho = img.propertyNames().filter(ee.Filter.stringContains('item', 'GAIN'));
    return ee.Dictionary({
        data:            ee.Date(img.get('system:time_start')).format('YYYY-MM-dd HH:mm'),
        SOLAR_ELEVATION: img.get('SOLAR_ELEVATION'),
        SOLAR_AZIMUTH:   img.get('SOLAR_AZIMUTH'),
        CLOUDCOVER:      img.get('CLOUDCOVER')
    }).combine(img.toDictionary(chavesGanho));
}

// Acrescenta um prefixo às chaves de um dicionário (ex.: B01_p50 → TOA_B01_p50)
function prefixar(dic, prefixo) {
    dic = ee.Dictionary(dic);
    var chaves = dic.keys();
    return dic.rename(chaves, chaves.map(function(k) { return ee.String(prefixo).cat(k); }));
}

// Fração de pixels saturados (DN ≥ 254) por banda, na pegada da cena
function fracaoSaturada(img, regiao) {
    var vnir = img.select(BANDAS_VNIR);
    return vnir.gte(DN_SAT).updateMask(vnir.gt(0)).reduceRegion({
        reducer: ee.Reducer.mean(), geometry: regiao, scale: ESCALA, maxPixels: 1e10, bestEffort: true
    });
}

var visDN  = {bands: ['B3N', 'B02', 'B01'], min: 20, max: 160};
var visTOA = {bands: ['B3N', 'B02', 'B01'], min: 0.03, max: 0.35, gamma: 1.3};

// ------- 3. EXECUÇÃO POR PAR -------
var resumo = [];   // uma linha (Feature) por cena → tabela única no fim

pares.forEach(function(par) {
    print('══════════════ PAR ' + par.nome + ' ══════════════');

    var cenas = {bom: ee.Image(par.bom), ruim: ee.Image(par.ruim)};
    var toa   = {};

    // Área comum entre as duas cenas (se houver, dá para comparar o mesmo chão)
    var comum = cenas.bom.geometry().intersection(cenas.ruim.geometry(), 100);
    var areaComum = comum.area(100).divide(1e6);

    ['bom', 'ruim'].forEach(function(tipo) {
        var img    = cenas[tipo];
        var regiao = img.geometry();
        var rotulo = par.nome + ' ' + tipo.toUpperCase() + ' (' + img.id().getInfo() + ')';
        var dn     = img.select(BANDAS_VNIR).updateMask(img.select(BANDAS_VNIR).gt(0));
        toa[tipo]  = dnParaTOA(img);

        // Linha do resumo: metadados + ganho + percentis DN/TOA + saturação + área comum
        resumo.push(ee.Feature(null, metadados(img)
            .combine(prefixar(percentis(dn, regiao), 'DN_'))
            .combine(prefixar(percentis(toa[tipo], regiao), 'TOA_'))
            .combine(prefixar(fracaoSaturada(img, regiao), 'SAT_'))
            .combine({par: par.nome, tipo: tipo, cena: img.id(),
                      area_cena_km2: regiao.area(100).divide(1e6),
                      area_comum_par_km2: areaComum})));

        print(histograma(dn, regiao,        'DN bruto — ' + rotulo,        'DN'));
        print(histograma(toa[tipo], regiao, 'Reflectância TOA — ' + rotulo, 'Reflectância TOA'));

        Map.addLayer(dn,        visDN,  rotulo + ' — falsa-cor DN',  false);
        Map.addLayer(toa[tipo], visTOA, rotulo + ' — falsa-cor TOA', tipo === 'bom');
    });

    // Histograma sobreposto da banda NIR: bom × ruim (cada um na sua pegada)
    var regiaoPar = cenas.bom.geometry().union(cenas.ruim.geometry(), 100);
    var sobreposto = ee.Image.cat(
        toa.bom.select(BANDA_PAR).rename('bom'),
        toa.ruim.select(BANDA_PAR).rename('ruim')
    );
    print(ui.Chart.image.histogram({image: sobreposto, region: regiaoPar, scale: ESCALA, maxPixels: 1e10})
        .setSeriesNames(['bom', 'ruim'])
        .setOptions({
            title: 'Par ' + par.nome + ' — ' + BANDA_PAR + ' (TOA): bom × ruim',
            hAxis: {title: 'Reflectância TOA'},
            vAxis: {title: 'Contagem de pixels'},
            colors: ['#2a78d6', '#e34948'],
            lineWidth: 2, pointSize: 0
        }));
});

// ------- 4. RESUMO ÚNICO (4 cenas) -------
// Tabela no Console (botão de download CSV no canto do gráfico) + CSV no Drive.
var tabelaResumo = ee.FeatureCollection(resumo);
print('══════════════ RESUMO DAS CENAS ══════════════');
print(ui.Chart.feature.byFeature(tabelaResumo, 'cena')
    .setChartType('Table')
    .setOptions({title: 'Resumo: metadados, ganho, percentis DN/TOA, saturação, área comum'}));
print('Resumo (expandir para ver todas as propriedades):', tabelaResumo);

Export.table.toDrive({
    collection: tabelaResumo,
    description: 'resumo_cenas_histograma',
    folder: PASTA_DRIVE,
    fileNamePrefix: 'resumo_cenas_histograma',
    fileFormat: 'CSV'
});

Map.centerObject(ee.Image(pares[0].bom), 9);
