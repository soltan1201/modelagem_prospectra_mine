// Visualização e análise dos mosaicos semestrais ASTER salvos por célula
// (saída de build_mosaicos_semestrais.py, versão v4).
//
// As coleções guardam um asset por ano × semestre × célula da grid 5 × 5
// (propriedades year, semestre, celula_id, celula_nome). Há duas coleções —
// a das contas caatinga (mapbiomas-arida) e a da superconta (mapbiomas-brazil);
// dá para ver cada uma ou as duas juntas. Este script:
//   - mostra quantos mosaicos (células) foram salvos por ano e semestre,
//     de 2000 a 2020 (zero onde não há nada);
//   - deixa escolher a coleção, o ano (ou todos), o semestre (ou os dois) e a
//     célula (ou todas); com "todos"/"os dois" vira um mosaico multitemporal;
//   - junta as células escolhidas com qualityMosaic e mostra camadas,
//     estatísticas, histogramas e o perfil espectral ao clicar no mapa.
//
// Valores no asset = valor físico × fator (propriedade SCALE_FACTORS):
//   B01–B09  reflectância TOA × 1000
//   B10–B14  DN TIR × 1   (aqui convertido para temperatura de brilho, K)
//   quality  quality × 1000  (5 = cena zero-cloud; ≤ 2 = cena com máscara)
//   saturado 0/1 (VNIR com DN ≥ 254; só aparece em cena zero-cloud)
//
// Seções:
//   1. Parâmetros
//   2. Conversão para valores físicos
//   3. Camadas no mapa
//   4. Estatísticas e histogramas no painel
//   5. Perfil espectral ao clicar no mapa
//   6. Contagem de mosaicos por ano e semestre
//   7. Interface (seletores de ano, semestre e célula)

// ------- 1. PARÂMETROS -------
var COLECOES = {
    caatinga:   'projects/mapbiomas-arida/mosaic_aster_v4',
    superconta: 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA/ASTER_v4'
};
var FONTE_INICIAL = 'ambas';   // 'caatinga' | 'superconta' | 'ambas'
var VER_TUDO_AO_ABRIR = true;  // abre já com tudo: as duas coleções, todos os anos e células
var ANOS = [];
for (var a_ = 2000; a_ <= 2020; a_++) { ANOS.push(a_); }
var id_grid      = 'projects/mapbiomas-arida/mine/grid_irece_5x5_solape100';
var escala_estat = 150;    // m — estatísticas (15 m nativo é lento demais)
var escala_hist  = 300;    // m — histogramas

var VNIR = ['B01', 'B02', 'B3N'];
var SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09'];
var TIR  = ['B10', 'B11', 'B12', 'B13', 'B14'];

// Usado só se a imagem não tiver SCALE_FACTORS (fatores da v4)
var FATOR_PADRAO = {B01: 1000, B02: 1000, B3N: 1000, B04: 1000, B05: 1000,
                    B06: 1000, B07: 1000, B08: 1000, B09: 1000,
                    B10: 1, B11: 1, B12: 1, B13: 1, B14: 1, quality: 1000, saturado: 1};

var QUALITY_ZERO_CLOUD = 5;   // físico (5000 no asset)

// Coeficientes de conversão DN → radiância do TIR (fixos no ASTER L1T, W/m²/sr/µm por DN)
var UCC_TIR = {B10: 0.006822, B11: 0.006780, B12: 0.006590, B13: 0.005693, B14: 0.005225};
// Constantes de Planck (K1 em W/m²/sr/µm, K2 em K) — mesmas de analise_todas_bandas.py
var K1 = {B10: 3040.136402, B11: 2482.375199, B12: 1935.060183, B13: 866.468575, B14: 641.326517};
var K2 = {B10: 1735.337945, B11: 1666.398761, B12: 1585.420044, B13: 1350.069147, B14: 1271.221673};
// Comprimento de onda central (µm) para o perfil espectral
var LAMBDA = {B01: 0.556, B02: 0.661, B3N: 0.807, B04: 1.656, B05: 2.167,
              B06: 2.209, B07: 2.262, B08: 2.336, B09: 2.400};

var PALETA_DIV = ['#104281', '#2a78d6', '#f0efec', '#e34948', '#9e2a2a'];

var grid = ee.FeatureCollection(id_grid);

// Coleção da fonte escolhida; em 'ambas', o mesmo ano × semestre × célula
// salvo nas duas coleções entra uma vez só (vale o da superconta)
function colecaoDa(fonte) {
    if (fonte !== 'ambas') { return ee.ImageCollection(COLECOES[fonte]); }
    var sup = ee.ImageCollection(COLECOES.superconta);
    var caa = ee.ImageCollection(COLECOES.caatinga);
    var nomes_sup = sup.aggregate_array('system:index');
    return sup.merge(caa.filter(ee.Filter.inList('system:index', nomes_sup).not()));
}

var colecao = null;
// Estado da seleção atual (usado pelo clique no mapa)
var atual = {fisico: null, opticas: []};
// Índice cliente: "ano|semestre" → lista de nomes de célula salvos
var indice = {};

function presentes(lista, bandas) {
    return lista.filter(function(b) { return bandas.indexOf(b) >= 0; });
}

// ------- 2. CONVERSÃO PARA VALORES FÍSICOS -------
// toFloat antes de dividir: o asset é INT32 e INT ÷ INT no GEE pode truncar
function paraFisico(bruto, bandas, fatores) {
    function fator(b) { return fatores[b] || FATOR_PADRAO[b]; }
    var partes = presentes(VNIR.concat(SWIR), bandas).map(function(b) {
        return bruto.select(b).toFloat().divide(fator(b));
    });
    presentes(TIR, bandas).forEach(function(b) {
        var dn  = bruto.select(b).toFloat().divide(fator(b));
        var rad = dn.subtract(1).multiply(UCC_TIR[b]);
        partes.push(ee.Image(K2[b]).divide(ee.Image(K1[b]).divide(rad).add(1).log())
                      .updateMask(dn.gt(1)).rename(b));
    });
    ['quality', 'saturado'].forEach(function(b) {
        if (bandas.indexOf(b) >= 0) { partes.push(bruto.select(b).toFloat().divide(fator(b))); }
    });
    return ee.Image.cat(partes);
}

// ------- 3. CAMADAS NO MAPA -------
function adicionarCamadas(fisico, bandas, nome, celulas_fc) {
    Map.layers().reset();
    Map.addLayer(fisico, {bands: ['B3N', 'B02', 'B01'], min: 0.03, max: 0.35, gamma: 1.2},
                 nome + ' — cor falsa (B3N/B02/B01)');
    if (presentes(['B04', 'B06', 'B08'], bandas).length === 3) {
        Map.addLayer(fisico, {bands: ['B04', 'B06', 'B08'], min: 0.05, max: 0.45},
                     'SWIR (B04/B06/B08)', false);
    }
    if (bandas.indexOf('B13') >= 0) {
        Map.addLayer(fisico, {bands: ['B13'], min: 285, max: 325, palette: PALETA_DIV},
                     'Temperatura de brilho B13 (K)', false);
    }
    var savi = fisico.expression('1.5 * (n - r) / (n + r + 0.5)',
                                 {n: fisico.select('B3N'), r: fisico.select('B02')});
    Map.addLayer(savi, {min: -0.1, max: 0.6, palette: ['#a6611a', '#f5f5dc', '#1a9641']},
                 'SAVI', false);
    if (bandas.indexOf('quality') >= 0) {
        var quality = fisico.select('quality');
        Map.addLayer(quality, {min: 0, max: 2, palette: PALETA_DIV.slice().reverse()},
                     'quality cenas com máscara (0–2)', false);
        Map.addLayer(quality.gte(QUALITY_ZERO_CLOUD).selfMask(), {palette: ['#104281']},
                     'Pixel de cena zero-cloud', false);
    }
    if (bandas.indexOf('saturado') >= 0) {
        Map.addLayer(fisico.select('saturado').selfMask(), {palette: ['#ff00ff']},
                     'VNIR saturado (DN ≥ 254)', false);
    }
    Map.addLayer(ee.Image().byte().paint(grid, 1, 1), {palette: ['#ffffff']}, 'Grid 5×5');
    Map.addLayer(ee.Image().byte().paint(celulas_fc, 1, 3), {palette: ['#ffcc00']},
                 'Células salvas da seleção');
}

// ------- 4. ESTATÍSTICAS E HISTOGRAMAS NO PAINEL -------
function reduzir(img, reducer, regiao) {
    return img.reduceRegion({reducer: reducer, geometry: regiao, scale: escala_estat,
                             maxPixels: 1e10, bestEffort: true, tileScale: 4});
}

// Histograma de faixa FIXA (fixedHistogram) por banda. O histograma automático
// (ui.Chart.image.histogram) escolhe a largura das faixas pelo mínimo/máximo
// da região; poucos pixels extremos (ex.: bordas, saturação) deixam a faixa
// enorme e tudo cai numa barra só. Pixels fora de [minimo, maximo] são ignorados.
function histograma(img, regiao, titulo, eixo, minimo, maximo, n_faixas) {
    var bandas = img.bandNames();
    var painel = ui.Panel([ui.Label('Calculando ' + titulo + '…', estilo_texto)]);
    img.reduceRegion({
        reducer: ee.Reducer.fixedHistogram(minimo, maximo, n_faixas),
        geometry: regiao, scale: escala_hist, maxPixels: 1e9, bestEffort: true, tileScale: 4
    }).evaluate(function(h, erro) {
        painel.clear();
        if (erro) { painel.add(ui.Label(titulo + ': erro ' + erro, {color: '#9e2a2a'})); return; }
        var nomes = Object.keys(h).filter(function(b) { return h[b] !== null; }).sort();
        if (!nomes.length) { painel.add(ui.Label(titulo + ': sem pixels na faixa', estilo_texto)); return; }
        // h[b] = [[início da faixa, contagem], ...] — mesmas faixas em todas as bandas
        var linhas = [[eixo].concat(nomes)];
        for (var i = 0; i < n_faixas; i++) {
            var linha = [h[nomes[0]][i][0] + (maximo - minimo) / n_faixas / 2];
            nomes.forEach(function(b) { linha.push(h[b][i][1]); });
            linhas.push(linha);
        }
        painel.add(ui.Chart(linhas, 'LineChart', {
            title: titulo, hAxis: {title: eixo}, vAxis: {title: 'pixels'},
            legend: {position: 'top'}, lineWidth: 2, pointSize: 0
        }));
    });
    return painel;
}

function mostrarAnalise(fisico, bandas, info, regiao) {
    var vnir = presentes(VNIR, bandas);
    var swir = presentes(SWIR, bandas);
    var tir  = presentes(TIR, bandas);

    resultados.clear();
    resultados.add(ui.Label('Seleção', estilo_titulo));
    resultados.add(ui.Label('Assets juntados: ' + info.n_assets, estilo_texto));
    resultados.add(ui.Label('Células: ' + info.celulas.join(', '), estilo_texto));
    resultados.add(ui.Label('Cenas usadas (soma das células): ' + info.num_images +
                            ' | zero-cloud: ' + info.num_zero_cloud, estilo_texto));
    resultados.add(ui.Label('versao: ' + info.versao, estilo_texto));

    var grupos = [];
    if (vnir.length) { grupos.push(fisico.select(vnir[0]).mask().rename('VNIR')); }
    if (swir.length) { grupos.push(fisico.select(swir[0]).mask().rename('SWIR')); }
    if (tir.length)  { grupos.push(fisico.select(tir[0]).mask().rename('TIR')); }
    var cobertura = ee.Image.cat(grupos).unmask(0);

    var opticas = fisico.select(vnir.concat(swir));
    var q = fisico.select('quality');
    var estat = ee.Dictionary({
        cobertura: reduzir(cobertura, ee.Reducer.mean(), regiao),
        fora_faixa: reduzir(opticas.lt(0).or(opticas.gt(1)), ee.Reducer.mean(), regiao),
        mediana: reduzir(fisico, ee.Reducer.median(), regiao),
        zero_cloud: reduzir(q.gte(QUALITY_ZERO_CLOUD).rename('zero_cloud'), ee.Reducer.mean(), regiao)
    });

    var carregando = ui.Label('Calculando estatísticas…', estilo_texto);
    resultados.add(carregando);
    estat.evaluate(function(e, erro) {
        resultados.remove(carregando);
        if (erro) {
            resultados.add(ui.Label('Erro: ' + erro, {color: '#9e2a2a'}));
            return;
        }
        resultados.add(ui.Label('Fração das células com dado', estilo_titulo));
        Object.keys(e.cobertura).forEach(function(k) {
            resultados.add(ui.Label(k + ': ' + (100 * e.cobertura[k]).toFixed(1) + ' %', estilo_texto));
        });
        var zc = e.zero_cloud.zero_cloud;
        resultados.add(ui.Label('Pixels vindos de cena zero-cloud: ' +
                                (zc === null ? '—' : (100 * zc).toFixed(1) + ' %'), estilo_texto));
        resultados.add(ui.Label('Mediana por banda (reflectância / K / quality)', estilo_titulo));
        bandas.forEach(function(b) {
            var v = e.mediana[b];
            resultados.add(ui.Label(b + ': ' + (v === null || v === undefined ? '—' : v.toFixed(4)),
                                    estilo_texto));
        });
        resultados.add(ui.Label('Reflectância fora de [0, 1] (fração de pixels)', estilo_titulo));
        Object.keys(e.fora_faixa).forEach(function(k) {
            var v = e.fora_faixa[k];
            resultados.add(ui.Label(k + ': ' + (v === null ? '—' : (100 * v).toFixed(2) + ' %'),
                                    estilo_texto));
        });
    });

    if (vnir.length) {
        resultados.add(histograma(fisico.select(vnir), regiao, 'Histograma VNIR',
                                  'reflectância TOA', 0, 0.6, 120));
    }
    if (swir.length) {
        resultados.add(histograma(fisico.select(swir), regiao, 'Histograma SWIR',
                                  'reflectância TOA', 0, 0.6, 120));
    }
    if (tir.length) {
        resultados.add(histograma(fisico.select(tir), regiao, 'Histograma TIR',
                                  'temperatura de brilho (K)', 270, 340, 140));
    }
    if (bandas.indexOf('quality') >= 0) {
        resultados.add(histograma(fisico.select('quality'), regiao, 'Histograma do quality',
                                  'quality (5 = zero-cloud)', 0, 5.2, 104));
    }
}

// Junta as células da seleção (qualityMosaic) e monta mapa + painel
function carregar() {
    var ano = sel_ano.getValue(), sem = sel_sem.getValue(), cel = sel_cel.getValue();
    if (!ano || !sem || !cel) { return; }

    var filtrada = colecao;
    if (ano !== 'todos') { filtrada = filtrada.filter(ee.Filter.eq('year', Number(ano))); }
    if (sem !== 'ambos') { filtrada = filtrada.filter(ee.Filter.eq('semestre', String(sem))); }
    if (cel !== 'todas') { filtrada = filtrada.filter(ee.Filter.eq('celula_nome', cel)); }
    var nome = (ano === 'todos' ? 'todos os anos' : ano) +
               (sem === 'ambos' ? ' (S1+S2)' : '-S' + sem) +
               (cel === 'todas' ? ' — todas as células' : ' — ' + cel);

    resultados.clear();
    resultados.add(ui.Label('Carregando ' + nome + '…', estilo_texto));
    var primeira = filtrada.first();
    ee.Dictionary({
        bandas: primeira.bandNames(),
        fatores: primeira.get('SCALE_FACTORS'),
        versao: primeira.get('versao'),
        celulas: filtrada.aggregate_array('celula_nome').distinct().sort(),
        n_assets: filtrada.size(),
        num_images: filtrada.aggregate_sum('num_images'),
        num_zero_cloud: filtrada.aggregate_sum('num_zero_cloud')
    }).evaluate(function(info, erro) {
        if (erro) {
            resultados.clear();
            resultados.add(ui.Label('Erro: ' + erro, {color: '#9e2a2a'}));
            return;
        }
        var bandas = info.bandas;
        var fatores = info.fatores || {};
        if (!info.fatores) { print('Aviso: sem SCALE_FACTORS — usando fatores padrão (v4).'); }

        // células vizinhas compartilham a grade de pixels: qualityMosaic junta sem emendas
        var bruto  = filtrada.qualityMosaic('quality');
        var fisico = paraFisico(bruto, bandas, fatores);
        atual.fisico  = fisico;
        atual.opticas = presentes(VNIR.concat(SWIR), bandas);

        var celulas_fc = grid.filter(ee.Filter.inList('nome', info.celulas));
        var regiao = celulas_fc.geometry();
        adicionarCamadas(fisico, bandas, nome, celulas_fc);
        Map.centerObject(celulas_fc, cel === 'todas' ? 7 : 9);
        mostrarAnalise(fisico, bandas, info, regiao);
        painel_pixel.clear();
        painel_pixel.add(ui.Label('Clique no mapa para ver o perfil espectral do pixel'));
    });
}

// ------- 5. PERFIL ESPECTRAL AO CLICAR -------
var painel_pixel = ui.Panel({style: {position: 'bottom-right', width: '400px'}});
painel_pixel.add(ui.Label('Escolha ano, semestre e célula; depois clique no mapa'));
Map.add(painel_pixel);
Map.style().set('cursor', 'crosshair');

Map.onClick(function(coords) {
    if (!atual.fisico) { return; }
    var ponto = ee.Geometry.Point([coords.lon, coords.lat]);
    painel_pixel.clear();
    painel_pixel.add(ui.Label('Carregando ' + coords.lon.toFixed(5) + ', ' + coords.lat.toFixed(5) + '…'));

    atual.fisico.reduceRegion({reducer: ee.Reducer.first(), geometry: ponto, scale: 15})
      .evaluate(function(v) {
        painel_pixel.clear();
        if (!v || v.B01 === null || v.B01 === undefined) {
            painel_pixel.add(ui.Label('Pixel sem dado.'));
            return;
        }
        var linhas = [['λ (µm)', 'reflectância']];
        atual.opticas.forEach(function(b) {
            if (v[b] !== null && v[b] !== undefined) { linhas.push([LAMBDA[b], v[b]]); }
        });
        painel_pixel.add(ui.Chart(linhas, 'LineChart', {
            title: 'Perfil VNIR+SWIR', legend: {position: 'none'},
            hAxis: {title: 'λ (µm)'}, vAxis: {title: 'reflectância TOA'},
            pointSize: 5, colors: ['#2a78d6']
        }));
        if (v.quality !== undefined && v.quality !== null) {
            painel_pixel.add(ui.Label('quality: ' + v.quality.toFixed(3) +
                (v.quality >= QUALITY_ZERO_CLOUD ? '  (cena zero-cloud)' : '  (cena com máscara)') +
                (v.saturado === 1 ? '  — VNIR SATURADO' : '')));
        }
        if (v.B13 !== undefined) {
            painel_pixel.add(ui.Label('B13: ' + (v.B13 === null ? '—' : v.B13.toFixed(1) + ' K') +
                '   B14: ' + (v.B14 === null || v.B14 === undefined ? '—' : v.B14.toFixed(1) + ' K')));
        }
    });
});

// ------- 6. CONTAGEM DE MOSAICOS POR ANO E SEMESTRE -------
// Tabela (ano × S1/S2) e gráfico de barras com o número de células salvas
function mostrarContagem() {
    var anos = ANOS;
    var tabela = [['Ano', 'S1', 'S2', 'Total']];
    var barras = [['Ano', 'S1', 'S2']];
    var total = 0;
    anos.forEach(function(a) {
        var n1 = (indice[a + '|1'] || []).length;
        var n2 = (indice[a + '|2'] || []).length;
        tabela.push([String(a), n1, n2, n1 + n2]);
        barras.push([String(a), n1, n2]);
        total += n1 + n2;
    });
    painel_contagem.clear();
    painel_contagem.add(ui.Label('Mosaicos salvos por ano e semestre (células, máx. 25) — ' +
                                 sel_fonte.getValue(), estilo_titulo));
    painel_contagem.add(ui.Label('Total: ' + total + ' assets', estilo_texto));
    painel_contagem.add(ui.Chart(tabela, 'Table', {allowHtml: true}));
    painel_contagem.add(ui.Chart(barras, 'ColumnChart', {
        title: 'Células salvas por semestre', isStacked: false, legend: {position: 'top'},
        hAxis: {title: 'ano', slantedText: true}, vAxis: {title: 'células', viewWindow: {min: 0, max: 25}},
        colors: ['#2a78d6', '#e34948']
    }));
}

// ------- 7. INTERFACE -------
var estilo_titulo = {fontWeight: 'bold', margin: '10px 0 2px 0'};
var estilo_texto  = {margin: '1px 0', fontSize: '12px'};

// células salvas para o ano/semestre escolhidos ('todos'/'ambos' = união)
function celulasDe(ano, sem) {
    var anos = ano === 'todos' ? ANOS : [ano];
    var sems = sem === 'ambos' ? ['1', '2'] : [sem];
    var cels = [];
    anos.forEach(function(a) {
        sems.forEach(function(s) {
            (indice[a + '|' + s] || []).forEach(function(c) {
                if (cels.indexOf(c) < 0) { cels.push(c); }
            });
        });
    });
    return cels.sort();
}

function contar(ano, sem) {
    var anos = ano === 'todos' ? ANOS : [ano];
    var n = 0;
    anos.forEach(function(a) { n += (indice[a + '|' + sem] || []).length; });
    return n;
}

function atualizarSemestres() {
    var ano = sel_ano.getValue();
    var n1 = contar(ano, '1'), n2 = contar(ano, '2');
    var itens = [];
    if (n1 && n2) { itens.push({label: 'os dois (' + (n1 + n2) + ' assets)', value: 'ambos'}); }
    if (n1) { itens.push({label: '1º semestre (' + n1 + ' assets)', value: '1'}); }
    if (n2) { itens.push({label: '2º semestre (' + n2 + ' assets)', value: '2'}); }
    sel_sem.items().reset(itens);
    sel_sem.setValue(itens.length ? itens[itens.length > 1 && ano !== 'todos' ? 1 : 0].value : null, false);
    atualizarCelulas();
}

function atualizarCelulas() {
    var cels = celulasDe(sel_ano.getValue(), sel_sem.getValue());
    var itens = [{label: 'todas (' + cels.length + ')', value: 'todas'}].concat(
        cels.map(function(c) { return {label: c, value: c}; }));
    sel_cel.items().reset(itens);
    sel_cel.setValue('todas', false);
    carregar();
}

var sel_fonte = ui.Select({
    items: [{label: 'caatinga + superconta', value: 'ambas'},
            {label: 'caatinga (mapbiomas-arida)', value: 'caatinga'},
            {label: 'superconta (mapbiomas-brazil)', value: 'superconta'}],
    value: FONTE_INICIAL,
    onChange: function(f) { carregarFonte(f); }
});
var sel_ano = ui.Select({placeholder: 'Ano', onChange: atualizarSemestres});
var sel_sem = ui.Select({placeholder: 'Semestre', onChange: atualizarCelulas});
var sel_cel = ui.Select({placeholder: 'Célula', onChange: carregar});
[sel_ano, sel_sem, sel_cel].forEach(function(s) { s.setDisabled(true); });

var painel_contagem = ui.Panel([ui.Label('Contando mosaicos salvos…', {fontSize: '12px'})]);
var resultados = ui.Panel();

var lateral = ui.Panel({
    widgets: [ui.Label('Mosaicos semestrais ASTER (por célula)', {fontSize: '18px', fontWeight: 'bold'}),
              ui.Panel([ui.Label('Coleção:'), sel_fonte], ui.Panel.Layout.flow('horizontal')),
              ui.Button({label: 'Ver todas as imagens (as duas coleções, todos os anos)',
                         onClick: verTudo, style: {stretch: 'horizontal'}}),
              ui.Panel([sel_ano, sel_sem, sel_cel], ui.Panel.Layout.flow('horizontal')),
              painel_contagem, resultados],
    style: {width: '440px'}
});
ui.root.insert(0, lateral);

Map.centerObject(grid, 7);
Map.setOptions('SATELLITE');
Map.addLayer(ee.Image().byte().paint(grid, 1, 1), {palette: ['#ffffff']}, 'Grid 5×5');

// Lê ano/semestre/célula de todos os assets da fonte e monta contagem e seletores
// verTudo: ao terminar de ler, seleciona todos os anos / os dois semestres /
// todas as células e desenha o mosaico de tudo
function carregarFonte(fonte, verTudo) {
    colecao = colecaoDa(fonte);
    indice = {};
    [sel_ano, sel_sem, sel_cel].forEach(function(s) { s.setDisabled(true); });
    painel_contagem.clear();
    painel_contagem.add(ui.Label('Contando mosaicos salvos…', estilo_texto));
    resultados.clear();

    ee.Dictionary({
        anos: colecao.aggregate_array('year'),
        sems: colecao.aggregate_array('semestre'),
        cels: colecao.aggregate_array('celula_nome')
    }).evaluate(function(d, erro) {
        if (erro) {
            painel_contagem.clear();
            painel_contagem.add(ui.Label('Erro ao listar (sem acesso à coleção?): ' + erro,
                                         {color: '#9e2a2a'}));
            return;
        }
        for (var i = 0; i < d.anos.length; i++) {
            var chave = d.anos[i] + '|' + d.sems[i];
            (indice[chave] = indice[chave] || []).push(d.cels[i]);
        }
        Object.keys(indice).forEach(function(k) { indice[k].sort(); });
        print('Mosaicos (' + fonte + '):', d.anos.length, indice);
        mostrarContagem();

        var anos = ANOS.filter(function(a) { return indice[a + '|1'] || indice[a + '|2']; });
        if (!anos.length) {
            sel_ano.items().reset([]);
            sel_ano.setPlaceholder('Coleção vazia');
            return;
        }
        sel_ano.items().reset([{label: 'todos os anos', value: 'todos'}].concat(
            anos.map(function(a) { return {label: String(a), value: String(a)}; })));
        [sel_ano, sel_sem, sel_cel].forEach(function(s) { s.setDisabled(false); });
        sel_ano.setPlaceholder('Ano');
        if (verTudo) {
            sel_ano.setValue('todos');   // dispara semestres → células → carregar
        }
    });
}

// Botão: as duas coleções juntas, todos os anos, os dois semestres, todas as células
function verTudo() {
    if (sel_fonte.getValue() !== 'ambas') {
        sel_fonte.setValue('ambas', false);
    }
    carregarFonte('ambas', true);
}

carregarFonte(FONTE_INICIAL, VER_TUDO_AO_ABRIR);
