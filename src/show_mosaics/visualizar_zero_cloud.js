// Visualização das cenas da lista zero-cloud (list_zero_cloud.txt):
// cenas inspecionadas por um analista e confirmadas sem nuvens, que entram
// intactas no mosaico (quality 5000).
//
// Escolha ano e semestre: cada cena vira uma camada própria no mapa (cor falsa
// B3N/B02/B01 em reflectância TOA com ganho, igual ao pipeline) e uma linha no
// painel com "ver" (mostra só aquela cena) e "tirar" (marca para sair da
// lista). As marcadas se acumulam entre semestres e aparecem num quadro de
// texto, já no formato do .txt, para copiar — o script não edita o arquivo.
// "Ver todas as cenas" (padrão ao abrir) mostra todas de uma vez.
//
// A lista abaixo foi gerada de list_zero_cloud.txt (324 cenas únicas em
// 33 semestres), agrupando cada cena pelo semestre da própria data —
// igual ao filtro de data do build. Se o arquivo mudar, gere de novo.

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
function vnirToa(img) {
    img = ee.Image(img);
    var cosz = ee.Number(90).subtract(img.get('SOLAR_ELEVATION')).multiply(Math.PI / 180).cos();
    var doy = ee.Number(ee.Date(img.get('system:time_start')).getRelative('day', 'year'));
    var d = ee.Number(1).subtract(ee.Number(0.01672).multiply(doy.multiply(2 * Math.PI / 365).cos()));
    var bandas = Object.keys(ESUN).map(function(b) {
        var dn = img.select(b);
        return dn.subtract(1).multiply(ee.Number(img.get('GAIN_COEFFICIENT_' + b)))
                 .multiply(Math.PI).multiply(d.pow(2)).divide(cosz.multiply(ESUN[b]))
                 .updateMask(dn.gt(0)).toFloat().rename(b);   // tipo fixo: mosaic exige coleção homogênea
    });
    var saturado = img.select(['B01', 'B02', 'B3N']).reduce(ee.Reducer.max())
                      .gte(DN_SATURADO).toByte().rename('saturado');
    // copyProperties devolve Element: ee.Image(...) para poder ir ao mapa
    return ee.Image(ee.Image.cat(bandas).addBands(saturado).copyProperties(img, img.propertyNames()));
}

// ------- 3. MAPA -------
var grid = ee.FeatureCollection(id_grid);
var bloco = ee.Geometry.Polygon(LIST_COORD);

function contorno(fc, cor, largura) {
    return ee.Image().byte().paint(fc, 1, largura).visualize({palette: [cor]});
}

function imagem(id) { return ee.Image('ASTER/AST_L1T_003/' + id); }

// Cenas marcadas para sair da lista (vale entre semestres)
var marcadas = [];
// camadas do semestre aberto: id → ui.Map.Layer
var camadas = {};

function atualizarMarcadas() {
    painel_marcadas.clear();
    if (!marcadas.length) {
        painel_marcadas.add(ui.Label('Nenhuma cena marcada.', estilo_texto));
        return;
    }
    marcadas.sort();
    painel_marcadas.add(ui.Label(marcadas.length + ' cena(s) para tirar da lista:', estilo_texto));
    painel_marcadas.add(ui.Textbox({
        value: marcadas.map(function(id) { return '"ASTER/AST_L1T_003/' + id + '",'; }).join('\n'),
        style: {width: '340px', fontSize: '11px', fontFamily: 'monospace'}
    }));
    painel_marcadas.add(ui.Button('Imprimir no Console', function() {
        print('Cenas para tirar de list_zero_cloud.txt:',
              marcadas.map(function(id) { return 'ASTER/AST_L1T_003/' + id; }));
    }));
}

function marcar(id, sim) {
    var i = marcadas.indexOf(id);
    if (sim && i < 0) { marcadas.push(id); }
    if (!sim && i >= 0) { marcadas.splice(i, 1); }
    if (camadas[id]) {
        camadas[id].setName((sim ? '✗ ' : '') + id);
    }
    atualizarMarcadas();
}

// Mostra só a cena `id` (as outras camadas do semestre ficam desligadas)
function soEssa(id) {
    Object.keys(camadas).forEach(function(k) { camadas[k].setShown(k === id); });
    Map.centerObject(imagem(id).geometry(), 10);
    mostrarMetadados(id);
}

// Modo "1 por 1": uma camada por cena (liga/desliga no painel de camadas do mapa)
// e uma linha por cena no painel lateral, com "ver" e "tirar da lista"
function abrirLista(ids, rotulo) {
    Map.layers().reset();
    camadas = {};
    Map.addLayer(contorno(grid, '#ffffff', 1), {}, 'Grid 5×5', false);
    ids.forEach(function(id, i) {
        var toa = vnirToa(imagem(id));
        var nome = (marcadas.indexOf(id) >= 0 ? '✗ ' : '') + id;
        camadas[id] = ui.Map.Layer(toa, VIS_TOA, nome, i === 0);
        Map.layers().add(camadas[id]);
    });
    Map.addLayer(contorno(ee.FeatureCollection(ids.map(function(id) {
        return ee.Feature(imagem(id).geometry(), {id: id}); })), '#ffcc00', 2), {}, 'Contorno das cenas');

    lista_cenas.clear();
    lista_cenas.add(ui.Label(rotulo + ': ' + ids.length +
                             ' cenas — "ver" mostra só a cena; marque as que saem da lista',
                             estilo_texto));
    ids.forEach(function(id, i) {
        var ver = ui.Button({label: 'ver', onClick: function() { soEssa(id); },
                             style: {padding: '0', margin: '0 4px 0 0'}});
        var tirar = ui.Checkbox({label: 'tirar', value: marcadas.indexOf(id) >= 0,
                                 onChange: function(v) { marcar(id, v); },
                                 style: {margin: '6px 0 0 4px'}});
        lista_cenas.add(ui.Panel([ui.Label((i + 1) + '. ' + id, {margin: '8px 6px 0 0', fontSize: '12px'}),
                                  ver, tirar], ui.Panel.Layout.flow('horizontal')));
    });
    if (ids.length) { soEssa(ids[0]); }
}

function mostrarMetadados(id) {
    var img = imagem(id);
    var toa = vnirToa(img);
    info.clear();
    info.add(ui.Label('Carregando metadados de ' + id + '…', estilo_texto));
    ee.Dictionary({
        data: ee.Date(img.get('system:time_start')).format('YYYY-MM-dd HH:mm'),
        CLOUDCOVER: img.get('CLOUDCOVER'),
        SOLAR_ELEVATION: img.get('SOLAR_ELEVATION'),
        tem_SWIR: img.propertyNames().contains('GAIN_COEFFICIENT_B04'),
        celulas: grid.filterBounds(img.geometry()).aggregate_array('nome'),
        saturado_pct: toa.select('saturado').reduceRegion({
            reducer: ee.Reducer.mean(), geometry: img.geometry(), scale: 90, maxPixels: 1e9,
            bestEffort: true}).get('saturado')
    }).evaluate(function(m, erro) {
        info.clear();
        if (erro) { info.add(ui.Label('Erro: ' + erro, {color: '#9e2a2a'})); return; }
        [['Cena', id], ['Data', m.data], ['CLOUDCOVER', m.CLOUDCOVER],
         ['Elevação solar', m.SOLAR_ELEVATION.toFixed(1) + '°'],
         ['SWIR', m.tem_SWIR ? 'sim' : 'não (a partir de 2008-S2)'],
         ['VNIR saturado', m.saturado_pct === null ? '—' : (100 * m.saturado_pct).toFixed(2) + ' %'],
         ['Células da grid', m.celulas.join(', ')]
        ].forEach(function(par) { info.add(ui.Label(par[0] + ': ' + par[1], estilo_texto)); });
    });
}

function todosIds() {
    var ids = [];
    Object.keys(LISTA_ZERO_CLOUD).sort().forEach(function(s) { ids = ids.concat(LISTA_ZERO_CLOUD[s]); });
    return ids;
}

// Modo "todas juntas": as cenas `ids` de uma vez (mosaic: a mais recente fica por cima)
function mostrarConjunto(ids, rotulo) {
    var col = ee.ImageCollection(ids.map(imagem)).sort('system:time_start');
    var n_cenas = col.map(function(i) {
        return i.select('B01').mask().gt(0).selfMask().rename('n');
    }).count().clip(bloco);

    Map.layers().reset();
    camadas = {};
    lista_cenas.clear();
    Map.addLayer(col.map(vnirToa).mosaic(), VIS_TOA, rotulo + ' — ' + ids.length + ' cenas zero-cloud');
    Map.addLayer(n_cenas, {min: 1, max: 15, palette: ['#f0efec', '#2a78d6', '#104281']},
                 'Nº de cenas zero-cloud por pixel', false);
    Map.addLayer(contorno(grid, '#ffffff', 1), {}, 'Grid 5×5');
    Map.centerObject(bloco, 7);

    info.clear();
    info.add(ui.Label(rotulo + ': ' + ids.length + ' cenas', estilo_texto));
    var por_celula = grid.map(function(f) { return f.set('n', col.filterBounds(f.geometry()).size()); });
    ee.Dictionary({nomes: por_celula.aggregate_array('nome'), n: por_celula.aggregate_array('n')})
      .evaluate(function(d, erro) {
        if (erro) { info.add(ui.Label('Erro: ' + erro, {color: '#9e2a2a'})); return; }
        var linhas = [['Célula', 'Cenas']];
        d.nomes.forEach(function(nm, i) { linhas.push([nm, d.n[i]]); });
        info.add(ui.Chart(linhas, 'ColumnChart', {
            title: 'Cenas zero-cloud por célula da grid', legend: {position: 'none'},
            hAxis: {slantedText: true}, vAxis: {title: 'cenas'}, colors: ['#2a78d6']}));
    });
}

// ------- 4. INTERFACE -------
var estilo_titulo = {fontWeight: 'bold', margin: '10px 0 2px 0'};
var estilo_texto  = {margin: '1px 0', fontSize: '12px'};
var semestres = Object.keys(LISTA_ZERO_CLOUD).sort();
var anos = [];
semestres.forEach(function(s) { var a = s.split('_')[0]; if (anos.indexOf(a) < 0) { anos.push(a); } });

var MAX_1POR1 = 80;   // acima disso o modo "1 por 1" vira "todas juntas" (camadas demais)

// cenas da seleção: ano ('todos' = todos) × semestre ('ambos' = 1 e 2)
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

function rotuloDe(ano, sem) {
    return (ano === 'todos' ? 'Todos os anos' : ano) +
           (sem === 'ambos' ? ' (S1+S2)' : ' — ' + sem + 'º semestre');
}

function abrir() {
    var ano = sel_ano.getValue(), sem = sel_sem.getValue();
    if (!ano || !sem) { return; }
    var ids = idsDe(ano, sem), rotulo = rotuloDe(ano, sem);
    if (sel_modo.getValue() === 'lista' && ids.length <= MAX_1POR1) {
        abrirLista(ids, rotulo);
    } else {
        mostrarConjunto(ids, rotulo);
        if (sel_modo.getValue() === 'lista') {
            lista_cenas.add(ui.Label(ids.length + ' cenas: muitas para "1 por 1" (máx. ' + MAX_1POR1 +
                                     ') — mostrando todas juntas. Escolha um ano ou semestre.',
                                     {fontSize: '12px', color: '#9e2a2a'}));
        }
    }
}

function contar(ano, sem) { return idsDe(ano, sem).length; }

var sel_ano = ui.Select({
    items: [{label: 'todos os anos', value: 'todos'}].concat(
        anos.map(function(a) { return {label: a, value: a}; })),
    placeholder: 'Ano',
    onChange: aoMudarAno
});

function aoMudarAno(ano) {
        var itens = [];
        var n1 = contar(ano, '1'), n2 = contar(ano, '2');
        if (n1 && n2) { itens.push({label: 'os dois (' + (n1 + n2) + ' cenas)', value: 'ambos'}); }
        if (n1) { itens.push({label: '1º semestre (' + n1 + ' cenas)', value: '1'}); }
        if (n2) { itens.push({label: '2º semestre (' + n2 + ' cenas)', value: '2'}); }
        sel_sem.items().reset(itens);
        sel_sem.setValue(itens[0].value, false);
        abrir();
}
var sel_sem = ui.Select({placeholder: 'Semestre', onChange: abrir});
var sel_modo = ui.Select({
    items: [{label: '1 por 1 (marcar cenas)', value: 'lista'},
            {label: 'todas juntas', value: 'juntas'}],
    value: 'lista', onChange: abrir
});

var bt_tudo = ui.Button({label: 'Ver todas as cenas (todos os anos)', style: {stretch: 'horizontal'},
                         onClick: function() {
                             sel_modo.setValue('juntas', false);
                             sel_ano.setValue('todos', false);
                             aoMudarAno('todos');         // semestre "os dois" + abrir
                         }});
var lista_cenas = ui.Panel();
var info = ui.Panel();
var painel_marcadas = ui.Panel();

var lateral = ui.Panel({
    widgets: [ui.Label('Cenas zero-cloud (inspecionadas)', {fontSize: '18px', fontWeight: 'bold'}),
              ui.Label('324 cenas em 33 semestres — list_zero_cloud.txt', {fontSize: '11px', color: '#666'}),
              ui.Panel([sel_ano, sel_sem], ui.Panel.Layout.flow('horizontal')),
              ui.Panel([ui.Label('Ver:'), sel_modo], ui.Panel.Layout.flow('horizontal')),
              bt_tudo,
              ui.Label('Cenas da seleção', estilo_titulo), lista_cenas,
              ui.Label('Cena selecionada', estilo_titulo), info,
              ui.Label('Marcadas para tirar da lista', estilo_titulo), painel_marcadas],
    style: {width: '400px'}
});
ui.root.insert(0, lateral);
atualizarMarcadas();

Map.setOptions('SATELLITE');
mostrarConjunto(todosIds(), 'Todos os anos');   // abre com todas; escolha ano/semestre/modo
