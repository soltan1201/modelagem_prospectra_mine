// ------- 1. PARÂMETROS E BLACKLIST -------
var list_coord = [
    [-44.460483519818084,-14.790793773956535],
    [-39.203525512005584,-14.790793773956535],
    [-39.203525512005584,-9.416525285369607], 
    [-44.460483519818084,-9.416525285369607], 
    [-44.460483519818084,-14.790793773956535]
];
var area_estudo = ee.Geometry.Polygon(list_coord);
Map.addLayer(area_estudo, {}, 'rectangulo');


// 1. Limiares de Nuvem
var lim_visBright = 0.30;      // => Vai para: visBright.gt(lim_visBright)
var lim_swir1 = 0.15;          // => Vai para: swir1.gt(lim_swir1) 
var ajuste_termico = 40;       // => Vai para: CLOUDCOVER / 2 + ajuste_termico
var raio_borda_nuvem = 180;    // => Vai para: cloud.focal_max({radius: raio_borda_nuvem})

// 2. Limiares de Sombra 
var altura_nuvem_estimada = 2000; // => Vai para: shadow_dist = ee.Number(altura_nuvem_estimada)
var lim_escuro_green = 0.18;   // => Vai para: green.lt(lim_escuro_green) 
var lim_escuro_nir = 0.20;    // => Vai para: nir.lt(lim_escuro_nir) 
var raio_borda_sombra = 90;    // => Vai para: shadow.focal_max({radius: raio_borda_sombra})

var vis_raw = {
    mosaico_clean:  { bands: ['B3N', 'B02', 'B01'], min: 0.05, max: 0.3 },
    mosaic_raw: {  bands: ['B3N', 'B02', 'B01'],  min: 25, max: 150  },
    mosaic_raw_vnir: { bands: ['B3N', 'B02', 'B01'], min: 30, max: 215, gamma: 1.15 },
    mosaic_raw_swir: { bands: ['B04', 'B06', 'B09'], min: 25, max: 170, gamma: 1.2 },
    layer_masc: { min: 0, max: 2, palette: ['black', 'red', 'blue'] }
}
var properti_CC = "CLOUDCOVER";
var ano = 2002;
var mes = 6; 
var semestre = '2';
// 🛑 COLE AQUI OS IDs DAS IMAGENS QUE VOCÊ QUER REMOVER
// Você encontra o ID clicando na imagem no mapa ou no Console
var blacklist = [
    // 'ASTER/AST_L1T_003/20210715132455', // Exemplo de ID (substitua pelos reais)
    "ASTER/AST_L1T_003/20020812130713",
    "ASTER/AST_L1T_003/20020812130722",
    "ASTER/AST_L1T_003/20020812130731",
    "ASTER/AST_L1T_003/20020812130740",
    "ASTER/AST_L1T_003/20020812130749",
    "ASTER/AST_L1T_003/20020911131942",
    "ASTER/AST_L1T_003/20020911131950",
    "ASTER/AST_L1T_003/20020913130712",
    "ASTER/AST_L1T_003/20020913130721",
    "ASTER/AST_L1T_003/20020913130730",
    "ASTER/AST_L1T_003/20020913130739",
    "ASTER/AST_L1T_003/20020913130748",
    "ASTER/AST_L1T_003/20020913130814",
    "ASTER/AST_L1T_003/20020913130823",
    "ASTER/AST_L1T_003/20020913130832",
    "ASTER/AST_L1T_003/20020915125535",
    "ASTER/AST_L1T_003/20020915125602",
    "ASTER/AST_L1T_003/20020915125611",
    "ASTER/AST_L1T_003/20020915125619",
    "ASTER/AST_L1T_003/20021022131259",
    "ASTER/AST_L1T_003/20021022131308",
    "ASTER/AST_L1T_003/20021022131316",
    "ASTER/AST_L1T_003/20021022131325",
    "ASTER/AST_L1T_003/20021022131410",
    "ASTER/AST_L1T_003/20021022131418",
    "ASTER/AST_L1T_003/20021022131427",
    "ASTER/AST_L1T_003/20020915125553",
];
// Extrai o ID final e adiciona prefixo "1_" para corresponder ao system:index do ASTER
var blacklist_clean = ee.List(blacklist.map(function(id) {
    return id.split('/').pop();  //'1_' +
}));

print('Blacklist limpa:', blacklist_clean);

var data_inicio = ee.Date.fromYMD(ano, mes, 1);
var data_fim = data_inicio.advance(6, 'month');

// ------- 2. FUNÇÕES DE PROCESSAMENTO -------

// -----------------------------
// 2) Radiância -> reflectância TOA
// apenas VNIR/SWIR
// -----------------------------
function radianceToTOA(imgRad) {
    // ASTER L1T já vem em RADIÂNCIA (W/m²/sr/µm)
    // Não precisa converter de DN! Os valores já são radiância.
    imgRad = ee.Image(imgRad);  
    
    var ESUN = {
        'B01': 1848, 'B02': 1549,  'B3N': 1114,
        'B04': 225.4,'B05': 86.63, 'B06': 81.85,
        'B07': 74.85,'B08': 66.49, 'B09': 59.85
    };

    var solarElev = ee.Number(imgRad.get('SOLAR_ELEVATION'));
    var sunZen = ee.Number(90).subtract(solarElev);
    var cosz = sunZen.multiply(Math.PI / 180).cos();

    var doy = ee.Number(ee.Date(imgRad.get('system:time_start')).getRelative('day', 'year'))//.add(1);
    var d = ee.Number(1).subtract(
        ee.Number(0.01672).multiply(
        doy.multiply(2 * Math.PI / 365).cos()
        )
    );

    var optical = Object.keys(ESUN).map(
        function(b) {
            var esun = ee.Number(ESUN[b]);
            return ee.Image(imgRad).select(b)
                        .multiply(Math.PI)
                        .multiply(d.pow(2))
                        .divide(esun.multiply(cosz))
                        .rename(b)
                        .float();
        }
    );

    return ee.Image(optical)
                .addBands(imgRad.select(['B10','B11','B12','B13','B14']).float())
                .copyProperties(imgRad, imgRad.propertyNames());
}

//https://code.earthengine.google.com/c7e78ccf5cc1de6273432f777b774377
function mascaraNuvemASTER(imgToa) {

    imgToa = ee.Image(imgToa);
    
    var green = imgToa.select('B01');
    var red   = imgToa.select('B02');
    var nir   = imgToa.select('B3N');
    var swir1 = imgToa.select('B04');
    var tir10 = imgToa.select('B10');

    // ===== GEOMETRIA SOLAR (CORREÇÃO DINÂMICA) =====
    var solar_az = ee.Number(imgToa.get('SOLAR_AZIMUTH'));
    var solar_el = ee.Number(imgToa.get('SOLAR_ELEVATION'));
    
    // Converte elevação para radianos e calcula tangente
    var solar_el_rad = solar_el.multiply(Math.PI / 180);
    var tan_solar_el = solar_el_rad.tan();
    
    // Estima distância da sombra baseado na elevação solar
    // Fórmula: distância = altura_nuvem / tan(elevação)
    // Altura típica de nuvens: 2000m (pode variar de 1000-4000m)
    var cloud_height = 2000; // metros
    
    var shadow_dist = ee.Number(cloud_height)
        .divide(tan_solar_el.max(0.1))  // Evita divisão por zero
        .max(200)   // Distância mínima: 200m
        .min(2000); // Distância máxima: 2000m

    var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
    // brilho visível
    var visBright = green.add(red).add(nir).divide(3);

    // nuvens: brilhantes, frias, pouco vegetadas
    var cloud = visBright.gt(0.25)
        // .and(swir1.gt(0.15))
        // .and(ndvi.lt(0.4));

    var percentCloud = ee.Number(imgToa.get('CLOUDCOVER')).divide(2).add(ajuste_termico);    
    // print("Percentil térmico calculado = ", percentCloud);

    // frio relativo por percentil da própria cena
    var tirStats = tir10.reduceRegion({
        reducer: ee.Reducer.percentile([20]),
        geometry: imgToa.geometry(),
        scale: 180,
        maxPixels: 1e10,
        bestEffort: true
    });
    var p20 = ee.Number(tirStats.get('B10'));
    cloud = cloud.and(tir10.lt(p20));
    // APLICANDO VARIÁVEL 
    var pmtro_focal = {radius: raio_borda_nuvem,  kernelType: 'square', units: 'meters'};
    cloud = cloud.focal_min(pmtro_focal)
                    .focal_max(pmtro_focal).rename('cloud');

    // ===== DETECÇÃO DE SOMBRAS COM CORREÇÃO GEOMÉTRICA =====
    // sombras: escuras e próximas das nuvens
    var dark_pixels = green.lt(lim_escuro_green).and(nir.lt(lim_escuro_nir));
    // Sombras estão na direção oposta ao sol
    // Usa a distância calculada dinamicamente baseada na elevação solar
    var shadow = dark_pixels.and(
        cloud.focal_max({
            radius: ee.Number(shadow_dist).add(300),
            units: 'meters'
        })
    );

    // Dilatação da sombra (suaviza bordas)
    pmtro_focal = {radius: raio_borda_sombra, kernelType: 'square', units: 'meters'}
    shadow = shadow.focal_max().rename('shadow');

    // máscara de pixels válidos
    var valid = imgToa.select(['B01','B02','B3N','B04']).reduce(ee.Reducer.min()).gt(0);
    var clear = valid.and(cloud.not()).and(shadow.not()).rename('clear_mask');

    return imgToa.addBands([ndvi])
                .updateMask(clear)
                copyProperties(imgToa, imgToa.propertyNames());
}

// 👉 FUNÇÃO QUALIDADE APRIMORADA (PARA QUALITYMOSAIC)
function addQualityASTER_v2(img) {
    
    // === 1. SELECIONAR BANDAS RELEVANTES ===
    var blue  = img.select('B01');   // Visível azul
    var green = img.select('B02');   // Visível verde  
    var nir   = img.select('B3N');   // NIR
    var swir1 = img.select('B04');   // SWIR
    var tir   = img.select('B10');   // TIR para temperatura
    
    // === 2. ÍNDICES ESPECTRAIS PARA QUALIDADE ===
    
    // NDVI: vegetação saudável tem assinatura estável
    var ndvi = nir.subtract(green).divide(nir.add(green).add(1e-6));
    
    // NDSI simplificado para detecção de nuvens (brilho + SWIR)
    var ndsi = green.subtract(swir1).divide(green.add(swir1).add(1e-6));
    
    // Brilho geral (nuvens são muito brilhantes em todas as bandas)
    var brightness = blue.add(green).add(nir).add(swir1).divide(4);
    
    // Razão NIR/SWIR: nuvens têm assinatura diferente de solo/vegetação
    var nirSwirRatio = nir.divide(swir1.add(1e-6));
    
    // === 3. MÉTRICAS DE TEXTURA LEVES (alternativa ao GLCM) ===
    
    // Desvio padrão local: nuvens têm alta variabilidade espacial
    var kernel = ee.Kernel.circle({radius: 30, units: 'meters'});
    var nirStd = nir.reduceNeighborhood(ee.Reducer.stdDev(), kernel);
    
    // === 4. CONSTRUIR SCORE DE QUALIDADE PIXEL A PIXEL ===
    
    // Componente 1: Preferir NDVI moderado-alto (evita água/nuvem)
    var scoreNdvi = ndvi.unitScale(-0.2, 0.7);  // -0.2 a 0.7 → 0 a 1
    
    // Componente 2: Penalizar brilho excessivo (nuvens)
    var scoreBright = ee.Image(1).subtract(
        brightness.unitScale(0.15, 0.5)  // >0.5 = muito brilhante
    );
    
    // Componente 3: Penalizar NDSI alto (assinatura de nuvem/gelo)
    var scoreNdsi = ee.Image(1).subtract(
        ndsi.unitScale(0.2, 0.6)  // NDSI >0.4 típico de nuvens
    );
    
    // Componente 4: Penalizar alta variabilidade local (bordas de nuvem)
    var scoreTexture = ee.Image(1).subtract(
        nirStd.unitScale(0.02, 0.1)  // std >0.1 = alta variabilidade
    );
    
    // Componente 5: Razão NIR/SWIR dentro de faixa "natural"
    var ratioDev = nirSwirRatio.subtract(2).abs();  // desvio do valor esperado ~2
    var scoreRatio = ee.Image(1).subtract(
        ratioDev.unitScale(0, 3)
    );
    
    // === 5. COMBINAR COMPONENTES COM PESOS ===
    var quality = scoreNdvi.multiply(0.25)
        .add(scoreBright.multiply(0.25))
        .add(scoreNdsi.multiply(0.20))
        .add(scoreTexture.multiply(0.15))
        .add(scoreRatio.multiply(0.15));
    
    // Normalizar para 0-1
    quality = quality.unitScale(0, 1).clamp(0, 1);
    
    // === 6. BÔNUS INTELIGENTE PARA CLOUDCOVER ===
    var cc = ee.Number(img.get(properti_CC));
    
    // Bônus escalonado: quanto menor o CC, maior o bônus
    var cloudBonus = ee.Image(1).subtract(
        cc.divide(70).min(1)  // CC=0 → +1, CC=70 → +0
    ).multiply(0.8);  // Peso máximo do bônus
    
    // Bônus extra para CC == 0 (cenas totalmente limpas)
    var perfectBonus = ee.Algorithms.If(cc.eq(0), 0.3, 0);
    
    quality = quality.add(cloudBonus).add(ee.Image.constant(perfectBonus));
    
    // === 7. MÁSCARA DE VALIDADE (opcional, mas recomendado) ===
    // Zera qualidade para pixels fora do range físico ou com mask original = 0
    var validPixel = img.select('B01').mask().and(
        brightness.gt(0.01).and(brightness.lt(0.8))  // Remove extremos
    );
    quality = quality.updateMask(validPixel);
    
    return img.addBands(quality.rename('quality'));  // redundância para compatibilidade
}

function addQualityASTER_CC0(img) {
  // === 1. SELEÇÃO DE BANDAS ===
  var base_qual  = img.select('B01').gt(0).multiply(25000);
  return img.addBands(base_qual.rename('quality'));
  
}


// 👉 FUNÇÃO QUALITY (PARA QUALITYMOSAIC)
// 👉 FUNÇÃO QUALIDADE LIGHT + SOMBRAS (PERFORMANCE + PRECISÃO)
// 👉 FUNÇÃO QUALIDADE COM PENALIDADE RÍGIDA PARA SOMBRAS
// 👉 FUNÇÃO QUALIDADE V5 - CORREÇÃO DE SOMBRAS E ÁGUA
function addQualityASTER_v5_shadow_fix(img) {
  
  // === 1. SELEÇÃO DE BANDAS ===
  var b1  = img.select('B01');
  var b2  = img.select('B02');
  var b3n = img.select('B3N');
  var b4  = img.select('B04');

  // === 2. ÍNDICES BASE ===
  var ndvi = b3n.subtract(b2).divide(b3n.add(b2).add(1e-6));
  var brightness = b1.add(b2).add(b3n).add(b4).divide(4);
  var ndsi = b2.subtract(b4).divide(b2.add(b4).add(1e-6));

  // === 3. SCORE BASE (Mantendo lógica anti-nuvem) ===
  // Penaliza brilho excessivo (nuvens) e NDSI alto
  var q_bright = ee.Image(1).subtract(brightness.unitScale(0.15, 0.5));
  var q_ndsi   = ee.Image(1).subtract(ndsi.unitScale(0.2, 0.6));
  var q_ndvi   = ndvi.unitScale(-0.1, 0.6); 

  var base_quality = q_bright.multiply(0.4)
                             .add(q_ndsi.multiply(0.3))
                             .add(q_ndvi.multiply(0.3))
                             .clamp(0, 1);

  // === 4. DETECÇÃO RÍGIDA DE SOMBRAS (NOVO) ===
  // Lógica: Sombra é um pixel MUITO escuro que NÃO é água.
  
  // a) Identificar Água (Absorve forte no SWIR e NIR)
  // Água típica: B4 < 0.04 e B3N < 0.06
  var is_water = b4.lt(0.04).and(b3n.lt(0.06));
  
  // b) Identificar Escuridão Extrema (Limiar de Sombra)
  // Pixels com B3N < 0.10 são suspeitos de serem sombra (ou água)
  var is_dark = b3n.lt(0.10);
  
  // c) Máscara Final: É escuro E NÃO é água -> É Sombra
  var shadow_mask = is_dark.and(is_water.not());
  
  // Suavização leve para pegar bordas de sombra (penumbra)
  shadow_mask = shadow_mask.focal_max({radius: 30, units: 'meters'});

  // === 5. BÔNUS DE CLOUD COVER ===
  var cc = ee.Number(img.get('CLOUDCOVER'));
  // Bônus escalonado + bônus extra para cena 100% limpa
  var cloud_bonus = ee.Image(1).subtract(cc.divide(80).min(1)).multiply(0.8);
  var perfect_bonus = cc.eq(0).multiply(0.3); 

  // === 6. APLICAÇÃO DO TETO DE QUALIDADE (O PULO DO GATO) ===
  
  // Calcula qualidade potencial (base + bônus)
  var potential_quality = base_quality.add(cloud_bonus).add(perfect_bonus);
  
  // 🔒 CRÍTICO: Se for sombra, FORÇA qualidade para 0.1 (valor baixo).
  // Isso anula qualquer bônus e garante que sombra PERCA para pixel limpo.
  var final_quality = ee.Image(0.1) 
                          .where(shadow_mask.not(), potential_quality)
                          .clamp(0, 2.0);

  // === 7. MÁSCARA DE VALIDADE ===
  var valid = b1.mask().and(brightness.gt(0.005));
  
  return img.addBands(final_quality.rename('quality'));
}

// 👉 FUNÇÃO PARA CONVERTER PARA 16-BIT (INTEIRO)
function converterPara16Bit(img) {
    // Fatores de escala para cada banda (preserva precisão)
    var escala = {
        'B01': 10000,   // Reflectância 0-1 -> 0-10000
        'B02': 10000,
        'B3N': 10000,
        'B04': 10000,
        'B05': 10000,
        'B06': 10000,
        'B07': 10000,
        'B08': 10000,
        'B09': 10000,
        'B10': 10,      // Temperatura (0-400K -> 0-4000)
        'B11': 10,
        'B12': 10,
        'B13': 10,
        'B14': 10,
        'NDVI': 10000,   // -1 a 1 -> -10000 a 10000
        'quality': 10000
    };
    
    var bandas = img.bandNames();
    
    var bandas_16bit = bandas.map(function(b) {
        b = ee.String(b);
        var fator = ee.Number(escala[b] || 10000);  // Default 10000
        
        var banda = img.select(b)
            .multiply(fator)
            .round()
            .toInt16()
            .rename([b]);
        
        return banda;
    });
    
    return ee.ImageCollection(bandas_16bit).toBands()
        .rename(bandas)
        .copyProperties(img, img.propertyNames())
        .set('SCALE_FACTORS', escala)
        .set('DATA_TYPE', 'INT16');
}


// Função para forçar tipo de forma mais robusta
function forcarTipoUniforme(img) {
    var bandas = ['B01','B02','B3N','B04','B05','B06','B07','B08','B09','B10','B11','B12','B13','B14'];
    
    var bandas_uniformes = bandas.map(function(b) {
        // Se a banda existe, converte para float com range fixo
        var banda = img.select(b);
        return banda.float().clamp(0, 1000);  // Range amplo para evitar overflow
    });
    
    return ee.Image(bandas_uniformes)
        .rename(bandas)
        .copyProperties(img, img.propertyNames());
}

// ------- 3. FILTRAGEM E APLICAÇÃO DA BLACKLIST -------

print(data_inicio, data_fim);
// ------- 3. PROCESSAMENTO DA COLEÇÃO -------
var colecao_base = ee.ImageCollection('ASTER/AST_L1T_003')
                    .filterDate(data_inicio, data_fim)
                    .filterBounds(area_estudo)
                    .filter(ee.Filter.lt('CLOUDCOVER', 70));
                    
var total_nFilter = colecao_base.size();

colecao_base = colecao_base
                    .filter(ee.Filter.inList('system:index', blacklist_clean).not())
                    .map(radianceToTOA);

// Imagens sem nuvens: CLOUDCOVER == 0 → não aplica máscara, mas adiciona NDVI para homogeneidade
var colecao_cc0 = colecao_base.filter(ee.Filter.eq(properti_CC, 0))
                    .map(function(img) {
                        var ndvi = img.select('B3N').subtract(img.select('B02'))
                                      .divide(img.select('B3N').add(img.select('B02')))
                                      .rename('NDVI');
                        var base_qual  = img.select('B01').gt(0).multiply(2.5);
                        return img.addBands(ndvi).addBands(base_qual.rename('quality'));
                    });
// Imagens com nuvens: CLOUDCOVER >= 1 → aplica mascaraNuvemASTER (já adiciona NDVI)
var colecao_cc_nuvens = colecao_base.filter(ee.Filter.gte(properti_CC, 1))
                             .map(mascaraNuvemASTER)
                             .map(addQualityASTER_v5_shadow_fix);  // 👈 Nova função
                             
                             
print("imagens com nuvens ", colecao_base.filterBounds(point));

var colecao = colecao_cc0.merge(colecao_cc_nuvens);
var numero_img = colecao.size();
print("Total de imagens:", numero_img);

// ------- 4. QUALITY MOSAIC (MELHOR QUE MEDIAN) -------
var mosaico_quality = colecao.qualityMosaic('quality').clip(area_estudo);

print("✅ Mosaico Quality criado:", mosaico_quality);
print(" total sem filtros ", total_nFilter);
print(colecao.limit(10));
// ------- 4. VISUALIZAÇÃO -------
// Map.centerObject(area_estudo, 8);

// Mosaico Resultante (Sem as imagens da Blacklist)
// RGB Natural
// Map.addLayer(mosaico_quality, {
//     bands: ['B3N', 'B02', 'B01'],
//     min: 0.05, max: 0.35,
//     gamma: 1.4
// }, 'Quality Mosaic - RGB');

// // Falsa-cor Geológica
// Map.addLayer(mosaico_quality, {
//     bands: ['B04', 'B06', 'B3N'],
//     min: 0.1, max: 0.5,
//     gamma: 1.2
// }, 'SWIR Falsa-cor');

// ------- 6. CONVERSÃO PARA 16-BIT E EXPORTAÇÃO -------
var mosaico_16bit = ee.Image(converterPara16Bit(mosaico_quality));
print("✅ Mosaico 16-bit:", mosaico_16bit);
print("Bandas disponíveis:", mosaico_16bit.bandNames());
var show_all_img = false;

mosaico_16bit = mosaico_16bit.set(
                  'year', ano,
                  'semestre', semestre,
                  'bloco', 'Irece',
                  'num_images', numero_img
                );
Map.addLayer(mosaico_16bit, {
    bands: ['B3N', 'B02', 'B01'],
    gamma:  1.4,
    max: 4354, min: 804,
}, 'Quality Mosaic - 16B RGB');
// 22855.92 : 15940.08
Map.addLayer(mosaico_16bit.select('quality'), {
    min: 15950, max: 22000,  // Ajustado para o novo range (base 1 + bônus até ~2.1)
    palette: ['red', 'yellow', 'green'],
    opacity: 0.6
}, 'Quality Band');
// Lista de bandas para exportar
var bandas_para_asset = [
    'B01', 'B02', 'B3N',           // VNIR 15m
    'B04', 'B05', 'B06', 'B07', 'B08', 'B09',  // SWIR 30m   
];
var id_output = 'projects/mapbiomas-arida/mosaic_aster/';

var name_to_export = 'ASTER_QualityMosaic_' + ano + '_semestre_' + semestre + '_INT16';
// 👉 EXPORTAR PARA ASSET (Formato 16-bit)
Export.image.toAsset({
    image: mosaico_16bit.select(bandas_para_asset),
    description: name_to_export,
    assetId: id_output  + name_to_export,
    region: area_estudo,
    scale: 15,  // Resolução do VNIR (SWIR será upsampled)
    maxPixels: 1e13,
    pyramidingPolicy: {
        '.default': 'sample'  // Melhor para dados contínuos
    }
});

if (show_all_img){
    colecao.evaluate(function(coll) {
        coll.features.forEach(function(feat) {
            var id_simples = feat.id.split('/').pop(); // Pega apenas o final do ID
            // var data = new Date(feat.properties['system:time_start']).toISOString().split('T')[0];
            var img = ee.Image(feat.id);
            
            // Mostra o ID no console para facilitar a cópia para a blacklist
            print('"' + feat.id + '",'); 
            
            Map.addLayer(img, vis_raw.mosaic_raw, 'Inspecionar: ' + ' _ ' + id_simples , false);
        });
    });
}