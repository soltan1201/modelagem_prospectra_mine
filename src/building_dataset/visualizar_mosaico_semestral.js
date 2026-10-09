// Visualização e análise dos mosaicos semestrais ASTER salvos
// (saída de build_mosaicos_semestrais.py).
//
// Lista as imagens da coleção id_colecao num seletor; ao escolher uma,
// carrega as camadas no mapa e mostra metadados, estatísticas e
// histogramas no painel lateral.
//
// Valores no asset = valor físico × fator. O fator vem da propriedade
// SCALE_FACTORS da imagem; se ela não existir, usa FATOR_PADRAO:
//   B01–B09  reflectância TOA × 10000
//   B10–B14  DN TIR × 10   (aqui convertido para temperatura de brilho, K)
//   quality  quality × 10000  (2,5 = cena zero-cloud; ≤ 2 = cena com máscara)
//
// Seções:
//   1. Parâmetros
//   2. Conversão para valores físicos
//   3. Camadas no mapa
//   4. Estatísticas e histogramas no painel
//   5. Perfil espectral ao clicar no mapa
//   6. Interface (seletor de imagem)

// ------- 1. PARÂMETROS -------
var id_colecao   = 'projects/mapbiomas-arida/mosaic_aster_v2';
var escala_estat = 150;    // m — estatísticas da cena inteira (15 m nativo é lento demais)
var escala_hist  = 300;    // m — histogramas

var VNIR = ['B01', 'B02', 'B3N'];
var SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09'];
var TIR  = ['B10', 'B11', 'B12', 'B13', 'B14'];

var FATOR_PADRAO = {B01: 10000, B02: 10000, B3N: 10000, B04: 10000, B05: 10000,
                    B06: 10000, B07: 10000, B08: 10000, B09: 10000,
                    B10: 10, B11: 10, B12: 10, B13: 10, B14: 10, quality: 10000};

// Coeficientes de conversão DN → radiância do TIR (fixos no ASTER L1T, W/m²/sr/µm por DN)
var UCC_TIR = {B10: 0.006822, B11: 0.006780, B12: 0.006590, B13: 0.005693, B14: 0.005225};
// Constantes de Planck (K1 em W/m²/sr/µm, K2 em K) — mesmas de analise_todas_bandas.py
var K1 = {B10: 3040.136402, B11: 2482.375199, B12: 1935.060183, B13: 866.468575, B14: 641.326517};
var K2 = {B10: 1735.337945, B11: 1666.398761, B12: 1585.420044, B13: 1350.069147, B14: 1271.221673};
// Comprimento de onda central (µm) para o perfil espectral
var LAMBDA = {B01: 0.556, B02: 0.661, B3N: 0.807, B04: 1.656, B05: 2.167,
              B06: 2.209, B07: 2.262, B08: 2.336, B09: 2.400};

var PALETA_DIV = ['#104281', '#2a78d6', '#f0efec', '#e34948', '#9e2a2a'];

// Estado da imagem selecionada (usado pelo clique no mapa)
var atual = {fisico: null, opticas: []};

function presentes(lista, bandas) {
    return lista.filter(function(b) { return bandas.indexOf(b) >= 0; });
}

// ------- 2. CONVERSÃO PARA VALORES FÍSICOS -------
function paraFisico(bruto, bandas, fatores) {
    var opticas = presentes(VNIR.concat(SWIR), bandas);
    var tir     = presentes(TIR, bandas);
    var partes  = opticas.map(function(b) {
        return bruto.select(b).divide(fatores[b] || FATOR_PADRAO[b]);
    });
    tir.forEach(function(b) {
        var dn  = bruto.select(b).divide(fatores[b] || FATOR_PADRAO[b]);
        var rad = dn.subtract(1).multiply(UCC_TIR[b]);
        partes.push(ee.Image(K2[b]).divide(ee.Image(K1[b]).divide(rad).add(1).log())
                      .updateMask(dn.gt(1)).rename(b));
    });
    if (bandas.indexOf('quality') >= 0) {
        partes.push(bruto.select('quality').divide(fatores.quality || FATOR_PADRAO.quality));
    }
    return ee.Image.cat(partes);
}

// ------- 3. CAMADAS NO MAPA -------
function adicionarCamadas(fisico, bandas, nome) {
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
        Map.addLayer(quality, {min: 0, max: 2.5, palette: PALETA_DIV.slice().reverse()},
                     'quality (0–2,5)', false);
        Map.addLayer(quality.eq(2.5).selfMask(), {palette: ['#104281']},
                     'Pixel de cena zero-cloud', false);
    }
}

// ------- 4. ESTATÍSTICAS E HISTOGRAMAS NO PAINEL -------
function reduzir(img, reducer, regiao) {
    return img.reduceRegion({reducer: reducer, geometry: regiao, scale: escala_estat,
                             maxPixels: 1e10, bestEffort: true, tileScale: 4});
}

function histograma(img, regiao, titulo, eixo, janela) {
    var opcoes = {title: titulo, hAxis: {title: eixo}};
    if (janela) { opcoes.hAxis.viewWindow = janela; }
    return ui.Chart.image.histogram({image: img, region: regiao, scale: escala_hist,
                                     maxPixels: 1e9, maxBuckets: 128}).setOptions(opcoes);
}

function mostrarAnalise(fisico, bandas, props, regiao) {
    var vnir = presentes(VNIR, bandas);
    var swir = presentes(SWIR, bandas);
    var tir  = presentes(TIR, bandas);

    resultados.clear();
    resultados.add(ui.Label('Metadados', estilo_titulo));
    ['year', 'semestre', 'num_images', 'num_zero_cloud', 'versao', 'DATA_TYPE'].forEach(function(k) {
        if (props[k] !== undefined) {
            resultados.add(ui.Label(k + ': ' + props[k], estilo_texto));
        }
    });
    resultados.add(ui.Label('Bandas: ' + bandas.join(', '), estilo_texto));

    // Cobertura por grupo de bandas (média da máscara 0/1)
    var grupos = [];
    if (vnir.length) { grupos.push(fisico.select(vnir[0]).mask().rename('VNIR')); }
    if (swir.length) { grupos.push(fisico.select(swir[0]).mask().rename('SWIR')); }
    if (tir.length)  { grupos.push(fisico.select(tir[0]).mask().rename('TIR')); }
    var cobertura = ee.Image.cat(grupos).unmask(0);

    var opticas = fisico.select(vnir.concat(swir));
    var estat = ee.Dictionary({
        cobertura: reduzir(cobertura, ee.Reducer.mean(), regiao),
        fora_faixa: reduzir(opticas.lt(0).or(opticas.gt(1)), ee.Reducer.mean(), regiao),
        mediana: reduzir(fisico, ee.Reducer.median(), regiao)
    });

    var carregando = ui.Label('Calculando estatísticas…', estilo_texto);
    resultados.add(carregando);
    estat.evaluate(function(e, erro) {
        resultados.remove(carregando);
        if (erro) {
            resultados.add(ui.Label('Erro: ' + erro, {color: '#9e2a2a'}));
            return;
        }
        resultados.add(ui.Label('Fração da área com dado', estilo_titulo));
        Object.keys(e.cobertura).forEach(function(k) {
            resultados.add(ui.Label(k + ': ' + (100 * e.cobertura[k]).toFixed(1) + ' %', estilo_texto));
        });
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
                                  'reflectância TOA', {min: 0, max: 0.6}));
    }
    if (swir.length) {
        resultados.add(histograma(fisico.select(swir), regiao, 'Histograma SWIR',
                                  'reflectância TOA', {min: 0, max: 0.6}));
    }
    if (tir.length) {
        resultados.add(histograma(fisico.select(tir), regiao, 'Histograma TIR',
                                  'temperatura de brilho (K)'));
    }
    if (bandas.indexOf('quality') >= 0) {
        resultados.add(histograma(fisico.select('quality'), regiao, 'Histograma do quality',
                                  'quality'));
    }
}

// Carrega a imagem escolhida: lê bandas e propriedades, depois monta tudo
function carregar(id) {
    var bruto = ee.Image(id_colecao + '/' + id);
    resultados.clear();
    resultados.add(ui.Label('Carregando ' + id + '…', estilo_texto));
    ee.Dictionary({bandas: bruto.bandNames(), props: bruto.toDictionary()})
      .evaluate(function(info, erro) {
        if (erro) {
            resultados.clear();
            resultados.add(ui.Label('Erro: ' + erro, {color: '#9e2a2a'}));
            return;
        }
        var bandas  = info.bandas;
        var fatores = info.props.SCALE_FACTORS || {};
        if (!info.props.SCALE_FACTORS) {
            print('Aviso: ' + id + ' não tem SCALE_FACTORS — usando fatores padrão.');
        }
        var fisico = paraFisico(bruto, bandas, fatores);
        atual.fisico  = fisico;
        atual.opticas = presentes(VNIR.concat(SWIR), bandas);

        adicionarCamadas(fisico, bandas, id);
        mostrarAnalise(fisico, bandas, info.props, bruto.geometry());
        painel_pixel.clear();
        painel_pixel.add(ui.Label('Clique no mapa para ver o perfil espectral do pixel'));
    });
}

// ------- 5. PERFIL ESPECTRAL AO CLICAR -------
var painel_pixel = ui.Panel({style: {position: 'bottom-right', width: '400px'}});
painel_pixel.add(ui.Label('Selecione uma imagem e clique no mapa'));
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
                (v.quality === 2.5 ? '  (cena zero-cloud)' : '  (cena com máscara)')));
        }
        if (v.B13 !== undefined) {
            painel_pixel.add(ui.Label('B13: ' + (v.B13 === null ? '—' : v.B13.toFixed(1) + ' K') +
                '   B14: ' + (v.B14 === null || v.B14 === undefined ? '—' : v.B14.toFixed(1) + ' K')));
        }
    });
});

// ------- 6. INTERFACE -------
var estilo_titulo = {fontWeight: 'bold', margin: '10px 0 2px 0'};
var estilo_texto  = {margin: '1px 0', fontSize: '12px'};

var seletor = ui.Select({placeholder: 'Carregando lista…', onChange: carregar});
seletor.setDisabled(true);
var resultados = ui.Panel();

var lateral = ui.Panel({
    widgets: [ui.Label('Mosaicos semestrais ASTER', {fontSize: '18px', fontWeight: 'bold'}),
              ui.Label(id_colecao, {fontSize: '11px', color: '#666'}),
              seletor, resultados],
    style: {width: '420px'}
});
ui.root.insert(0, lateral);

var colecao = ee.ImageCollection(id_colecao);
Map.centerObject(colecao.geometry(), 7);
Map.setOptions('SATELLITE');

colecao.aggregate_array('system:index').evaluate(function(ids, erro) {
    if (erro) {
        seletor.setPlaceholder('Erro ao listar: ' + erro);
        return;
    }
    ids.sort();
    print('Imagens em ' + id_colecao + ':', ids.length, ids);
    seletor.items().reset(ids);
    seletor.setPlaceholder('Escolha um mosaico (' + ids.length + ')');
    seletor.setDisabled(false);
});
