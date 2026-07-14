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
var ano = 2009;
var mes = 1; 
var semestre = '2';
if (mes === 1){semestre = '1';}
// 🛑 COLE AQUI OS IDs DAS IMAGENS QUE VOCÊ QUER REMOVER
// Você encontra o ID clicando na imagem no mapa ou no Console
var blacklist = [
    // 'ASTER/AST_L1T_003/20210715132455', // Exemplo de ID (substitua pelos reais)
    "ASTER/AST_L1T_003/20090628130700",
    "ASTER/AST_L1T_003/20090628130625",
    "ASTER/AST_L1T_003/20090628130634",
    "ASTER/AST_L1T_003/20090621130047",
    "ASTER/AST_L1T_003/20090619131252",
    "ASTER/AST_L1T_003/20090619131243",
    "ASTER/AST_L1T_003/20090614125517",
    "ASTER/AST_L1T_003/20090612130740",
    "ASTER/AST_L1T_003/20090612130713",
    "ASTER/AST_L1T_003/20090612130731",
    "ASTER/AST_L1T_003/20090612130638",
    "ASTER/AST_L1T_003/20090612130620",
    "ASTER/AST_L1T_003/20090605130023",
    "ASTER/AST_L1T_003/20090603131304",
    "ASTER/AST_L1T_003/20090603131255",
    "ASTER/AST_L1T_003/20090527130631",
    "ASTER/AST_L1T_003/20090520130134",
    "ASTER/AST_L1T_003/20090520130116",
    "ASTER/AST_L1T_003/20090513125509",
    "ASTER/AST_L1T_003/20090513125500",
    "ASTER/AST_L1T_003/20090504130149",
    "ASTER/AST_L1T_003/20090502131336",
    "ASTER/AST_L1T_003/20090425130727",
    "ASTER/AST_L1T_003/20090425130736",
    "ASTER/AST_L1T_003/20090425130745",
    "ASTER/AST_L1T_003/20090418130033",
    "ASTER/AST_L1T_003/20090416131406",
    "ASTER/AST_L1T_003/20090416131330",
    "ASTER/AST_L1T_003/20090407131915",
    "ASTER/AST_L1T_003/20090310125539",
    "ASTER/AST_L1T_003/20090317130116",
    "ASTER/AST_L1T_003/20090317130142",
    "ASTER/AST_L1T_003/20090310125521",
    "ASTER/AST_L1T_003/20090301130155",
    "ASTER/AST_L1T_003/20090301130120",
    "ASTER/AST_L1T_003/20090301130129",
    "ASTER/AST_L1T_003/20090301130138",
    "ASTER/AST_L1T_003/20090301130147",
    "ASTER/AST_L1T_003/20090301130053",
    "ASTER/AST_L1T_003/20090213130027",
    "ASTER/AST_L1T_003/20090126131322",
    "ASTER/AST_L1T_003/20090121125458",
    "ASTER/AST_L1T_003/20090121125507",
    "ASTER/AST_L1T_003/20090117131856",
    "ASTER/AST_L1T_003/20090112130133",
    "ASTER/AST_L1T_003/20090112130057",
    "ASTER/AST_L1T_003/20090112130106",
    "ASTER/AST_L1T_003/20090112130115",
    "ASTER/AST_L1T_003/20090112130124",
    "ASTER/AST_L1T_003/20090110131252",
    "ASTER/AST_L1T_003/20090103130705",
    "ASTER/AST_L1T_003/20090103130647",
    "ASTER/AST_L1T_003/20090103130638",
    "ASTER/AST_L1T_003/20090310125503",
    "ASTER/AST_L1T_003/20090310125548",
    "ASTER/AST_L1T_003/20090326125533",
    "ASTER/AST_L1T_003/20090407131857",
    "ASTER/AST_L1T_003/20090614125450",
    "ASTER/AST_L1T_003/20090619131217",
    "ASTER/AST_L1T_003/20090630125437",
];
// Extrai o ID final e adiciona prefixo "1_" para corresponder ao system:index do ASTER
var blacklist_clean = ee.List(blacklist.map(function(id) {
    return id.split('/').pop();  //'1_' +
}));

print('Blacklist limpa:', blacklist_clean);

var list_Cloud_zero = [    
    "ASTER/AST_L1T_003/20090101131917",
    "ASTER/AST_L1T_003/20090525132004",
    "ASTER/AST_L1T_003/20090525132013",
    "ASTER/AST_L1T_003/20090527130657",
    "ASTER/AST_L1T_003/20090619131328",
    "ASTER/AST_L1T_003/20090619131328",
    "ASTER/AST_L1T_003/20090626131915",
]
var Cloudlist_clean = ee.List(list_Cloud_zero.map(function(id) {
    return id.split('/').pop();  //'1_' +
}));

print('Cloudlist limpa:', Cloudlist_clean);

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

// Calcula SAVI, NDWI partir de uma imagem ASTER TOA (reflectância)
// Bandas ASTER usadas: B01=Verde, B02=Vermelho, B3N=NIR, B04=SWIR1(1.6µm), B05=SWIR2(2.1µm)
function addIndicesASTER(img) {
    img = ee.Image(img);
    var green = img.select('B01');  // Verde  0.52–0.60µm
    var red   = img.select('B02');  // Vermelho 0.63–0.69µm
    var nir   = img.select('B3N');  // NIR  0.76–0.86µm
    var swir1 = img.select('B04');  // SWIR1 1.60–1.70µm
    var swir2 = img.select('B05');  // SWIR2 2.145–2.185µm

    // SAVI - Soil-Adjusted Vegetation Index
    // Fórmula: 1.5·(NIR − RED) / (NIR + RED + 0.5) | Vegetação esparsa, Caatinga, cerrado aberto
    var savi = nir.subtract(red).multiply(1.5)
                  .divide(nir.add(red).add(0.5))
                  .rename('SAVI');

    // NDWI - Normalized Difference Water Index (McFeeters, 1996)
    // Fórmula: (GREEN − NIR) / (GREEN + NIR) | Água livre; valores > 0 indicam superfície hídrica
    var ndwi = green.subtract(nir)
                    .divide(green.add(nir).add(1e-6))
                    .rename('NDWI');   

    return img.addBands([savi, ndwi]);
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

    return addIndicesASTER(imgToa)
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
    
    // SAVI - Soil-Adjusted Vegetation Index
    // Fórmula: 1.5·(NIR − RED) / (NIR + RED + 0.5) | Vegetação esparsa, Caatinga, cerrado aberto
    var savi = nir.subtract(green).multiply(1.5).divide(nir.add(green).add(0.5).add(1e-6));
    
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
    
    // Componente 1: Preferir SAVI moderado-alto (evita água/nuvem)
    var scoreNdvi = savi.unitScale(-0.2, 0.7);  // -0.2 a 0.7 → 0 a 1
    
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
  var base_qual  = img.select('B01').gt(0).multiply(22000);   // antes do 2009 25000
  return img.addBands(base_qual.rename('quality'));
  
}


// 👉 FUNÇÃO QUALIDADE V6 — VNIR + TIR (sem SWIR, compatível com ASTER pós-2008)
// Discriminadores principais:
//   TIR (B10): nuvens frias (~970-1014) vs superfície quente (~1160-1245)
//   Brilho VNIR: nuvens muito brilhantes (>0.35), sombras muito escuras (<0.09)
//   SAVI: alto em vegetação → prefere pixels vegetados
//   NDWI: muito negativo em vegetação (-0.40), próximo a 0 em nuvens (-0.12)
function addQualityASTER_v5_shadow_fix(img) {

  // === 1. SELEÇÃO DE BANDAS ===
  var b1  = img.select('B01');  // Verde
  var b2  = img.select('B02');  // Vermelho
  var b3n = img.select('B3N');  // NIR
  var b10 = img.select('B10');  // TIR — principal separador nuvem/superfície

  // === 2. ÍNDICES ESPECTRAIS ===
  var savi = b3n.subtract(b2).multiply(1.5).divide(b3n.add(b2).add(0.5).add(1e-6));
  var ndwi = b1.subtract(b3n).divide(b1.add(b3n).add(1e-6));
  var brightness = b1.add(b2).add(b3n).divide(3);

  // === 3. SCORES INDIVIDUAIS ===

  // TIR: superfície quente = alta qualidade; nuvem fria = qualidade zero
  // Intervalos medidos: nuvem B10~970-1014; sombra~1091; veg~1160-1191; solo~1245
  var tir_score = b10.unitScale(1050, 1300).clamp(0, 1);

  // Brilho: U-invertido — penaliza escuro (sombra) e brilhante (nuvem), premia moderado
  var bright_up   = brightness.unitScale(0.07, 0.16).clamp(0, 1);
  var bright_down = ee.Image(1).subtract(brightness.unitScale(0.22, 0.42).clamp(0, 1));
  var bright_score = bright_up.multiply(bright_down);

  // SAVI: vegetação densa → 1; solo/nuvem → menor pontuação
  var savi_score = savi.unitScale(-0.1, 0.5).clamp(0, 1);

  // NDWI: quanto mais negativo, mais vegetado; nuvens têm NDWI próximo a 0
  var ndwi_score = ndwi.multiply(-1).unitScale(-0.10, 0.45).clamp(0, 1);

  // === 4. QUALIDADE COMBINADA ===
  var base_quality = tir_score.multiply(0.45)
                              .add(bright_score.multiply(0.25))
                              .add(savi_score.multiply(0.15))
                              .add(ndwi_score.multiply(0.15))
                              .clamp(0, 1);

  // === 5. BÔNUS DE CLOUD COVER ===
  var cc = ee.Number(img.get('CLOUDCOVER'));
  var cloud_bonus  = ee.Image(1).subtract(cc.divide(80).min(1)).multiply(0.8);
  var perfect_bonus = cc.eq(0).multiply(0.3);

  var final_quality = base_quality.add(cloud_bonus).add(perfect_bonus).clamp(0, 2.0);

  // === 6. MÁSCARA DE VALIDADE ===
  var valid = b1.mask().and(brightness.gt(0.005));
  final_quality = final_quality.updateMask(valid);

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
        'SAVI': 10000,   // -1 a 1  -> -10000 a 10000
        'NDWI': 10000,   // -1 a 1  -> -10000 a 10000
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
print(" total sem filtros ", total_nFilter);
colecao_base = colecao_base
                    .filter(ee.Filter.inList('system:index', blacklist_clean).not())
                    .map(radianceToTOA);



// Imagens sem nuvens: CLOUDCOVER == 0 → não aplica máscara, mas adiciona SAVI/NDWI para homogeneidade
var colecao_cc0 = colecao_base.filter(ee.Filter.inList('system:index', Cloudlist_clean))
                    .map(function(img) {
                        var base_qual = img.select('B01').gt(0).multiply(2.5);
                        return addIndicesASTER(img).addBands(base_qual.rename('quality'));
                    });
// Imagens com nuvens: CLOUDCOVER >= 1 → aplica mascaraNuvemASTER (já adiciona NDVI)
var colecao_cc_nuvens = colecao_base.filter(ee.Filter.gte(properti_CC, 1))
                             .map(mascaraNuvemASTER)
                             .map(addQualityASTER_v5_shadow_fix);  // 👈 Nova função
                             

var colecao = colecao_cc0.merge(colecao_cc_nuvens);
var numero_img = colecao.size();
print("Total de imagens:", numero_img);

// ------- 4. QUALITY MOSAIC (MELHOR QUE MEDIAN) -------
var mosaico_quality = colecao.qualityMosaic('quality').clip(area_estudo);
mosaico_quality = mosaico_quality.updateMask(mosaico_quality.select('B01').gt(0));
print("✅ Mosaico Quality criado:", mosaico_quality);

print(colecao.limit(10));
// ------- 4. VISUALIZAÇÃO -------
// Map.centerObject(area_estudo, 8);

// Mosaico Resultante (Sem as imagens da Blacklist)
// RGB Natural
Map.addLayer(mosaico_quality, {
    bands: ['B3N', 'B02', 'B01'],
    min: 0.05, max: 0.35,
    gamma: 1.4
}, 'Quality Mosaic - RGB');
var vis_index = {
    ind_savi : {min: -0.04, max: 0.5, palette: ["#e3dc1c", "#51b254"]},
    ind_ndwi : {min: -0.6, max: 0.3, palette: ["#c8a96e", "#0e5faf"]},
}
Map.addLayer(mosaico_quality.select("SAVI"), vis_index.ind_savi, 'SAVI Mosaic', false);
Map.addLayer(mosaico_quality.select("NDWI"), vis_index.ind_ndwi, 'NDWI Mosaic', false);


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
    'B10','B11','B12','B13','B14', 'quality'
];
var id_output = 'projects/mapbiomas-arida/mosaic_aster_p2008/';

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