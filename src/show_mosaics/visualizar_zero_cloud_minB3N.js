// Mosaico das cenas zero-cloud com a regra "fica o pixel de MENOR B3N".
//
// Em cor falsa (B3N/B02/B01) o vermelho forte é NIR alto, ou seja, vegetação
// verde (estação chuvosa). Escolhendo, em cada pixel, a observação de menor
// B3N, o mosaico tende a mostrar a estação seca: menos vegetação e mais solo /
// rocha expostos — o que interessa para a geologia.
//
// Cuidado: B3N baixo também é sombra (nuvem, relevo), água e área queimada; a
// camada "data do pixel escolhido" e o clique no mapa ajudam a conferir de que
// cena veio cada pixel.
//
// Seleção: ano (ou todos) × semestre (ou os dois). Camadas para comparar:
//   - menor B3N (regra nova)
//   - mais recente por cima (mosaic comum, como no visualizar_zero_cloud.js)
//   - maior B3N (o oposto, para ver o efeito)
//
// A lista é a mesma do visualizar_zero_cloud.js (com as cenas comentadas lá
// já fora). Se a lista mudar lá, copie o bloco LISTA_ZERO_CLOUD para cá.

// ------- 1. PARÂMETROS -------
var id_grid = 'projects/mapbiomas-arida/mine/grid_irece_5x5_solape100';
var LIST_COORD = [[-44.460483519818084, -14.790793773956535], [-39.203525512005584, -14.790793773956535],
                  [-39.203525512005584, -9.416525285369607], [-44.460483519818084, -9.416525285369607],
                  [-44.460483519818084, -14.790793773956535]];
var DN_SATURADO = 254;
var ESUN = {B01: 1848, B02: 1549, B3N: 1114};
var VIS_TOA = {bands: ['B3N', 'B02', 'B01'], min: 0.03, max: 0.35, gamma: 1.2};

var LISTA_ZERO_CLOUD = {
    '2000_1': [
        '20000518132219',
        // '20000628131449'  poderia
    ],
    '2000_2': [
        '20000703133323',
        '20000703133332',
        '20000703133340',
        '20000710133923',
        '20000710133932',
        '20000710133941',
        '20000710133950',
        '20000710133958',
        '20000813132747',
        '20000813132805'
    ],
    '2001_1': [
        // '20010212133109',    poderia
        // '20010212133118',    poderia
        '20010212133127',
        '20010401132857',
        '20010401132906',
        '20010401132914',
        '20010505131643',
        '20010505131652',
        '20010521131542'
    ],
    '2001_2': [
        '20010706132700',
        '20010713133256',
        '20010713133304',
        '20010713133313',
        '20010814133213',
        '20010901132011',
        // '20010901132020'  poderia
    ],
    '2002_2': [
        '20021022131334',
        '20021022131343',
        '20021022131352',
        '20021022131401'
    ],
    '2003_1': [
        '20030306131938',
        '20030306131946',
        '20030306131955',
        '20030306132004'
    ],
    '2003_2': [
        '20030914131726'
    ],
    '2006_1': [
        // '20060517131818',   poderia
        '20060517131836',
        '20060517131845'
    ],
    '2006_2': [
        '20060708125414',
        '20060708125423',
        '20060711132349',
        '20060711132416',
        '20060711132425',
        '20060720131749',
        '20060720131758',
        '20060720131806',
        '20060720131815'
    ],
    '2007_1': [
        '20070121131310',
        '20070121131327'
    ],
    '2007_2': [
        // '20070709130553', 
        // '20070918131157'
    ],
    '2008_2': [
        '20081015130717',
        '20081015130726'
    ],
    '2009_1': [
        '20090101131917',
        '20090525132004',
        '20090525132013',
        '20090527130657',
        '20090619131328',
        '20090626131915'
    ],
    '2009_2': [
        '20090728131856',
        '20090728131905',
        '20090728131914',
        '20090728131923',
        '20090728131932',
        '20090728131941',
        '20090813131820',
        '20090925125935',
        '20090925125944',
        '20091018130534',
        '20091018130543',
        '20091018130711',
        '20091128125945'
    ],
    '2010_1': [
        '20100528131853',
        '20100528131902',
        '20100528131911',
        '20100528131919',
        '20100606131159',
        '20100606131234',
        '20100606131243',
        '20100606131252',
        '20100606131301',
        '20100606131310',
        '20100606131318',
        '20100606131327',
        '20100613131818',
        '20100613131836',
        '20100617125428',
        '20100617125437',
        '20100620132411',
        '20100620132420',
        '20100620132429',
        '20100620132437',
        '20100620132446',
        '20100620132455',
        '20100620132504',
        '20100620132513',
        '20100622131242',
        '20100622131309',
        '20100622131318'
    ],
    '2010_2': [
        // '20100703125452',  // remover de vez
        '20100708131154',
        '20100708131203',
        '20100708131211',
        '20100809131212',
        '20100809131221',
        '20100809131229',
        '20100809131238',
        '20100809131247',
        '20100809131256',
        '20100809131305',
        '20100816131814',
        '20100816131823',
        '20100816131831',
        '20100816131840',
        '20100816131849',
        '20100816131858',
        '20100816131907',
        '20100816131916',
        '20100825131205'
    ],
    '2011_1': [
        '20110609131121',
        '20110609131130',
        '20110609131139',
        '20110609131148',
        '20110609131157',
        '20110611130026'
    ],
    '2012_2': [
        '20120704131747',
        '20120704131756',
        '20120704131805',
        '20120704131814',
        '20120704131823',
        '20120704131832',
        '20120704131840',
        '20120704131849',
        '20120704131907',
        '20120704131916',
        '20120704131925',
        '20120713131155',
        '20120713131204',
        '20120713131213',
        '20120713131221',
        '20120713131257',
        '20120713131306'
    ],
    '2013_1': [
        '20130520131900',
        '20130520131918',
        '20130520131927',
        '20130529131315',
        '20130529131324'
    ],
    '2013_2': [
        '20130707131820',
        '20130707131829',
        '20130707131838',
        '20130716131321',
        '20130723131756',
        '20130723131805',
        '20130723131814',
        '20130723131823',
        '20130723131832',
        '20130723131841',
        '20130723131849',
        '20130723131858',
        '20130723131907',
        '20130723131916',
        '20130723131925',
        '20130909131819',
        '20130925131820',
        '20131004131143',
        '20131004131152',
        '20131004131219',
        '20131004131228',
        '20131004131237',
        '20131004131245',
        '20131006125946',
        '20131027131834',
        '20131027131842',
        '20131027131851',
        '20131027131900',
        '20131027131909',
        '20131027131918'
    ],
    '2014_1': [
        '20140603130057',
        '20140603130106',
        '20140603130115',
        '20140603130124'
    ],
    '2014_2': [
        '20140703131313',
        '20140703131322',
        '20140703131331',
        '20140703131340',
        '20140703131348',
        '20140710131822',
        '20140921131208',
        '20140921131216',
        '20140921131243',
        '20141014131844',
        '20141014131853',
        '20141014131902',
        '20141014131911',
        '20141014131920',
        '20141014131937',
        '20141014131946',
        '20141014131955',
        '20141014132004',
        '20141108131217'
    ],
    '2015_1': [
        '20150102131844',
        '20150102131853',
        '20150102131901',
        '20150102131910',
        // '20150104130549',    remover
        '20150526131925',
        '20150526131943',
        '20150526131952',
        '20150526132001',
        '20150526132010',
        '20150526132019'
    ],
    '2015_2': [
        '20150807131315',
        '20150807131323',
        '20150807131332',
        '20150807131341',
        '20150807131350',
        '20150807131359',
        '20150807131408',
        '20150807131417',
        '20150807131425',
        '20150807131434',
        '20150901130747',
        '20150901130756',
        '20150901130805',
        '20150908131231',
        '20150908131239',
        '20150908131248',
        '20150908131257',
        '20150908131341',
        '20150915131934',
        '20150915131943',
        '20150915131952',
        '20150915132001',
        '20151026131227',
        '20151026131245',
        '20151026131253',
        '20151026131302',
        '20151026131320',
        '20151026131329'
    ],
    '2016_1': [
        '20160514130615'
    ],
    '2016_2': [
        '20160706132505',
        '20160706132514',
        '20160708131252',
        '20160708131301',
        '20160708131310',
        '20160717130631',
        '20160717130639',
        '20160717130808',
        '20160809131240',
        '20160809131249',
        '20160809131258',
        '20160809131307',
        '20160809131316',
        '20160809131325',
        '20160809131333',
        '20160809131342',
        '20160809131351',
        '20160809131400',
        '20160919130709',
        '20161003131829',
        '20161019131816'
    ],
    '2017_1': [
        '20170107131921',
        '20170107131930',
        '20170305131145',
        '20170602130702',
        '20170602130711',
        '20170602130720'
    ],
    '2017_2': [
        '20170702131923',
        '20170702131932',
        '20170911132504',
        '20170911132513',
        '20170911132522',
        '20170911132531',
        '20170911132540'
    ],
    '2018_2': [
        '20180916131310',
        '20180916131328',
        '20180916131355',
        '20180916131404',
        '20180916131412',
        '20180916131430',
        '20181025131912'
    ],
    '2019_1': [
        // '20190304130639',  // remover
        // '20190304130647',   // remover
        // '20190304130705',   // remover
        '20190428131231',
        '20190428131306',
        '20190521131916',
        '20190521131951'
    ],
    '2019_2': [
        '20190717131229',
        '20190717131238',
        '20190717131247',
        '20190717131256',
        '20190717131305',
        '20190717131314',
        '20190717131322',
        '20190724131755',
        '20190724131804',
        '20190724131813',
        '20190903131246',
        '20190912130639',
        '20190912130648',
        '20190912130657'
    ],
    '2020_1': [
        // '20200111130031',    // remnover
        // '20200111130040',     // remover
        // '20200405131953'     // remover
    ],
    '2020_2': [
        '20200827131755',
        '20200827131831',
        '20200827131839',
        '20200928131846',
        '20200928131855',
        '20200928131904',
        '20200928131913',
        '20200928131922',
        '20201007131142',
        '20201007131150',
        '20201007131159',
        '20201007131208',
        '20201007131217',
        '20201007131252',
        '20201007131301',
        '20201009125946',
        '20201016130547'
    ]
};

// ------- 2. TOA (só VNIR, igual ao radiance_to_toa do pipeline) -------
// Bandas extras para conferir a origem do pixel: 'data' (AAAAMMDD) e 'mes'
function vnirToa(img) {
    img = ee.Image(img);
    var cosz = ee.Number(90).subtract(img.get('SOLAR_ELEVATION')).multiply(Math.PI / 180).cos();
    var data = ee.Date(img.get('system:time_start'));
    var doy = ee.Number(data.getRelative('day', 'year'));
    var d = ee.Number(1).subtract(ee.Number(0.01672).multiply(doy.multiply(2 * Math.PI / 365).cos()));
    var bandas = Object.keys(ESUN).map(function(b) {
        var dn = img.select(b);
        return dn.subtract(1).multiply(ee.Number(img.get('GAIN_COEFFICIENT_' + b)))
                 .multiply(Math.PI).multiply(d.pow(2)).divide(cosz.multiply(ESUN[b]))
                 .updateMask(dn.gt(0)).toFloat().rename(b);
    });
    var vnir = ee.Image.cat(bandas);
    var saturado = img.select(['B01', 'B02', 'B3N']).reduce(ee.Reducer.max())
                      .gte(DN_SATURADO).toByte().rename('saturado');
    var cod_data = ee.Image.constant(ee.Number.parse(data.format('YYYYMMdd'))).toInt32().rename('data');
    var mes = ee.Image.constant(data.get('month')).toByte().rename('mes');
    // prioridade = −B3N: o qualityMosaic pega o MAIOR valor, ou seja, o menor B3N
    var prioridade = vnir.select('B3N').multiply(-1).rename('prioridade');
    return ee.Image(vnir.addBands([saturado, cod_data, mes, prioridade])
                        .updateMask(vnir.select('B3N').mask())
                        .copyProperties(img, img.propertyNames()));
}

// ------- 3. MAPA -------
var grid = ee.FeatureCollection(id_grid);
var bloco = ee.Geometry.Polygon(LIST_COORD);
var semestres = Object.keys(LISTA_ZERO_CLOUD).sort();
var anos = [];
semestres.forEach(function(s) { var a = s.split('_')[0]; if (anos.indexOf(a) < 0) { anos.push(a); } });
var atual = null;   // mosaico de menor B3N da seleção (para o clique)

function contorno(fc, cor, largura) {
    return ee.Image().byte().paint(fc, 1, largura).visualize({palette: [cor]});
}

function idsDe(ano, sem) {
    var ids = [];
    semestres.forEach(function(chave) {
        var p = chave.split('_');
        if ((ano === 'todos' || p[0] === ano) && (sem === 'ambos' || p[1] === sem)) {
            ids = ids.concat(LISTA_ZERO_CLOUD[chave]);
        }
    });
    return ids;
}

function mostrar() {
    var ano = sel_ano.getValue(), sem = sel_sem.getValue();
    if (!ano || !sem) { return; }
    var ids = idsDe(ano, sem);
    var rotulo = (ano === 'todos' ? 'Todos os anos' : ano) +
                 (sem === 'ambos' ? ' (S1+S2)' : ' — ' + sem + 'º semestre');
    var col = ee.ImageCollection(ids.map(function(id) { return ee.Image('ASTER/AST_L1T_003/' + id); }))
                .sort('system:time_start').map(vnirToa);

    var menor = col.qualityMosaic('prioridade');
    var maior = col.qualityMosaic('B3N');
    var recente = col.mosaic();
    atual = menor;

    Map.layers().reset();
    Map.addLayer(recente, VIS_TOA, 'Mais recente por cima (mosaic comum)', false);
    Map.addLayer(maior, VIS_TOA, 'Maior B3N (oposto)', false);
    Map.addLayer(menor, VIS_TOA, 'MENOR B3N — ' + rotulo);
    Map.addLayer(menor.select('mes'), {min: 1, max: 12,
                 palette: ['#104281', '#2a78d6', '#7fb3e6', '#f0efec', '#f5c38a', '#e34948',
                           '#9e2a2a', '#e34948', '#f5c38a', '#f0efec', '#7fb3e6', '#2a78d6']},
                 'Mês do pixel escolhido', false);
    Map.addLayer(menor.select('saturado').selfMask(), {palette: ['#ff00ff']},
                 'VNIR saturado no pixel escolhido', false);
    Map.addLayer(contorno(grid, '#ffffff', 1), {}, 'Grid 5×5');
    Map.centerObject(bloco, 7);

    info.clear();
    info.add(ui.Label(rotulo + ': ' + ids.length + ' cenas', estilo_texto));
    info.add(ui.Label('Clique no mapa para ver de que cena veio o pixel.', estilo_texto));
    // de quais meses vêm os pixels escolhidos (amostra a 500 m)
    menor.select('mes').reduceRegion({
        reducer: ee.Reducer.frequencyHistogram(), geometry: bloco, scale: 500,
        maxPixels: 1e9, bestEffort: true
    }).get('mes').evaluate(function(h, erro) {
        if (erro || !h) { return; }
        var total = Object.keys(h).reduce(function(s, k) { return s + h[k]; }, 0);
        var linhas = [['Mês', '% dos pixels']];
        for (var m = 1; m <= 12; m++) { linhas.push([String(m), 100 * (h[m] || 0) / total]); }
        info.add(ui.Chart(linhas, 'ColumnChart', {
            title: 'Mês de origem dos pixels (menor B3N)', legend: {position: 'none'},
            vAxis: {title: '% dos pixels'}, colors: ['#2a78d6']}));
    });
}

// Clique: data, mês e valores do pixel escolhido
Map.onClick(function(c) {
    if (!atual) { return; }
    pixel.clear();
    pixel.add(ui.Label('Carregando…', estilo_texto));
    atual.select(['data', 'B01', 'B02', 'B3N', 'saturado'])
      .reduceRegion({reducer: ee.Reducer.first(), geometry: ee.Geometry.Point([c.lon, c.lat]), scale: 15})
      .evaluate(function(v) {
        pixel.clear();
        if (!v || v.data === null) { pixel.add(ui.Label('Sem dado.', estilo_texto)); return; }
        pixel.add(ui.Label('Cena de ' + v.data + (v.saturado ? '  — VNIR saturado' : ''), estilo_texto));
        pixel.add(ui.Label('B01 ' + v.B01.toFixed(3) + ' | B02 ' + v.B02.toFixed(3) +
                           ' | B3N ' + v.B3N.toFixed(3), estilo_texto));
    });
});

// ------- 4. INTERFACE -------
var estilo_titulo = {fontWeight: 'bold', margin: '10px 0 2px 0'};
var estilo_texto  = {margin: '1px 0', fontSize: '12px'};

function aoMudarAno(ano) {
    var n1 = idsDe(ano, '1').length, n2 = idsDe(ano, '2').length;
    var itens = [];
    if (n1 && n2) { itens.push({label: 'os dois (' + (n1 + n2) + ' cenas)', value: 'ambos'}); }
    if (n1) { itens.push({label: '1º semestre (' + n1 + ' cenas)', value: '1'}); }
    if (n2) { itens.push({label: '2º semestre (' + n2 + ' cenas)', value: '2'}); }
    sel_sem.items().reset(itens);
    sel_sem.setValue(itens[0].value, false);
    mostrar();
}

var sel_ano = ui.Select({
    items: [{label: 'todos os anos', value: 'todos'}].concat(anos.map(function(a) { return {label: a, value: a}; })),
    value: 'todos', onChange: aoMudarAno
});
var sel_sem = ui.Select({placeholder: 'Semestre', onChange: mostrar});
var info = ui.Panel();
var pixel = ui.Panel();

ui.root.insert(0, ui.Panel({
    widgets: [ui.Label('Zero-cloud — menor B3N', {fontSize: '18px', fontWeight: 'bold'}),
              ui.Label('Em cada pixel fica a cena de menor B3N (menos vegetação verde).',
                       {fontSize: '11px', color: '#666'}),
              ui.Panel([sel_ano, sel_sem], ui.Panel.Layout.flow('horizontal')),
              ui.Label('Seleção', estilo_titulo), info,
              ui.Label('Pixel clicado', estilo_titulo), pixel],
    style: {width: '380px'}
}));

Map.setOptions('SATELLITE');
Map.style().set('cursor', 'crosshair');
aoMudarAno('todos');
