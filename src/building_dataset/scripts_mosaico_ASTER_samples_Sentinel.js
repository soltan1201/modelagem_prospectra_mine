// ------- 1. PARÂMETROS -------
var area_estudo = ee.FeatureCollection("projects/mapbiomas-arida/ESTADOS_IBGE_2024")
                    .filter(ee.Filter.eq('CD_UF', '29')).geometry();

var ano_central = 2021;
var ano_inicio = ano_central - 1;
var ano_fim = ano_central + 1;

var start_date = ee.Date.fromYMD(ano_inicio, 1, 1);
var end_date = ee.Date.fromYMD(ano_fim, 12, 31);
var periodo_str = ano_inicio + '_' + ano_fim;

var bandas_aster = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09'];

// ------- 2. DICIONÁRIO DE BANDAS (NOMES PADRONIZADOS) -------
var bandas_nomes_novos = [
    'AST_Green_556nm',
    'AST_Red_661nm', 
    'AST_NIR_807nm',
    'AST_SWIR1_1656nm',
    'AST_SWIR2_2167nm_AlOH',
    'AST_SWIR3_2209nm_Ref',
    'AST_SWIR4_2263nm_AlMgOH',
    'AST_SWIR5_2336nm_CO3',
    'AST_SWIR6_2395nm'
];

// ------- 3. FUNÇÕES ASTER (CORRIGIDAS) -------
var ESUN_ASTER = {
    'B01': 1845.99, 'B02': 1555.74, 'B3N': 1119.47,
    'B04': 231.25, 'B05': 79.81, 'B06': 74.99, 
    'B07': 68.66, 'B08': 59.74, 'B09': 56.92
};

function asterTOA(img) {
    var solarElev = ee.Number(img.get('SOLAR_ELEVATION'));
    var sunZenith = ee.Number(90).subtract(solarElev);
    var sunZenithRad = sunZenith.multiply(Math.PI / 180);
    
    var doy = img.date().getRelative('day', 'year').add(1);
    var dist = ee.Number(1).subtract(ee.Number(0.01673).multiply(doy.multiply(2).multiply(Math.PI).divide(365).cos()));
    var d2 = dist.pow(2);
    
    var bands = bandas_aster.map(function(b) {
        var irr = ee.Number(ESUN_ASTER[b]);
        var refl = img.select(b).multiply(Math.PI).multiply(d2)
                       .divide(irr.multiply(sunZenithRad.cos()))
                       .rename(b);
        return refl.float();
    });
    
    var termal = img.select('B10');
    return ee.Image(bands).addBands(termal).copyProperties(img, img.propertyNames());
}

function mascaraNuvemASTER(img) {
    var b10 = img.select('B10');
    var b01 = img.select('B01');
    
    // Ajuste empírico para semiárido baiano
    var cloud = b01.gt(0.28).and(b10.lt(7.2)); 
    cloud = cloud.focal_max({radius: 60, units: 'meters'});
    
    var b3n = img.select('B3N');
    var shadow = b3n.lt(0.06).and(cloud.focal_max(120));
    
    return img.updateMask(cloud.or(shadow).not());
}

// ------- 4. CORREÇÃO: GLCM COM CONVERSÃO PARA BYTE -------
function addQualityASTER(img) {
    var b3n = img.select('B3N');
    
    // 🔧 CORREÇÃO: Converte Float (0-0.4) para Byte (0-255)
    var b3n_byte = b3n.unitScale(0, 0.4).multiply(255).toByte();
    
    // Agora o GLCM funciona!
    var glcm = b3n_byte.glcmTexture({size: 1});
    var contrast = glcm.select('B3N_contrast');
    
    // Qualidade: Alto NIR + Baixo Contraste (sem bordas de nuvem)
    // Normalizamos ambos para 0-1
    var quality = b3n.unitScale(0, 0.35).add(contrast.unitScale(0, 255).multiply(-1)).rename('quality');
    
    return img.addBands(quality);
}

// ------- 5. PREPARAR LANDSAT AZUL (30m) -------
function prepararLandsatAzul() {

    // Função auxiliar para processar cada satélite
    function processarLandsat(colecao_id) {
        return ee.ImageCollection(colecao_id)
            .filterDate(start_date, end_date)
            .filterBounds(area_estudo)
            .map(function(img) {
                // Escala para Reflectância de Superfície (SR)
                var azul = img.select('SR_B2')
                    .multiply(0.0000275)
                    .add(-0.2);
                
                // Máscara de qualidade
                var qa = img.select('QA_PIXEL');
                
                // Bits de QA_PIXEL (Landsat Collection 2)
                var cloud = qa.bitwiseAnd(1 << 3).neq(0);  // Cloud
                var shadow = qa.bitwiseAnd(1 << 4).neq(0); // Cloud Shadow
                var snow = qa.bitwiseAnd(1 << 5).neq(0);   // Snow/Ice
                var water = qa.bitwiseAnd(1 << 7).neq(0);  // Water (opcional)
                
                // Máscara composta (mantém apenas pixels limpos)
                var mask = cloud.or(shadow).or(snow).not();
                
                // Saturação (valores > 1.0 são inválidos)
                var valid_range = azul.lt(1.0).and(azul.gt(-0.1));
                
                return azul.updateMask(mask.and(valid_range))
                           .rename('Blue_SR')
                           .copyProperties(img, ['system:time_start']);
            });
    }

    // Processa L7, L8 e L9
    var l8 = processarLandsat('LANDSAT/LC08/C02/T1_L2');
    var l9 = processarLandsat('LANDSAT/LC09/C02/T1_L2');

    // Combina todas as coleções
    var colecao_total = l7.merge(l8).merge(l9);
    
    // Estatísticas para debug
    var count = colecao_total.size();
    print('📊 Total de cenas Landsat Azul:', count);
    
    // Mediana temporal (reduz ruído e nuvens residuais)
    var mosaico_azul = colecao_total
                            .median()
                            .clip(area_estudo)
                            .rename('Blue_30m');

    
        // 🔧 PÓS-PROCESSAMENTO: Preenche gaps com interpolação espacial
    var kernel_suave = ee.Kernel.gaussian({
        radius: 60, 
        units: 'meters',
        sigma: 30
    });
    
    // Preenche pixels vazios com média dos vizinhos
    var azul_preenchido = mosaico_azul
        .unmask({
            sameFootprint: false,
            method: 'nearest'
        })
        .focal_mean({
            kernel: kernel_suave,
            iterations: 2
        })
        .rename('Blue_30m');
    
    return azul_preenchido;
}


// ------- FUNÇÃO SENTINEL-2 AZUL (10m NATIVO) -------
function prepararSentinel2Azul() {
    var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterDate(start_date, end_date)
        .filterBounds(area_estudo)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
        .map(function(img) {
            // Sentinel-2 SR já está em reflectância (escala 1-10000)
            var azul = img.select('B2')
                .divide(10000)
                .rename('Blue_S2');
            
            // Máscara de nuvens SCL
            var scl = img.select('SCL');
            var mask = scl.eq(4)  // Vegetation
                         .or(scl.eq(5))  // Not Vegetated
                         .or(scl.eq(6))  // Water
                         .or(scl.eq(7)); // Unclassified
            
            // Filtro adicional de nuvens
            var cloudProb = img.select('MSK_CLDPRB');
            mask = mask.and(cloudProb.lt(30));
            
            return azul.updateMask(mask)
                       .copyProperties(img, ['system:time_start']);
        });
    
    var mosaico_azul = s2.median()
        .clip(area_estudo)
        .rename('Blue_10m');
    
    return mosaico_azul;
}
// ------- 6. FUSÃO ASTER + LANDSAT (HPF) -------
// // Use esta versão no script principal (equilíbrio qualidade/performance)
// function fuseAzulComASTER_Otimizado(img_aster, img_azul_30m) {
//     var pan_15m = img_aster.select('AST_NIR_807nm');
    
//     // 1. Degradação suave (evita aliasing)
//     var pan_30m = pan_15m
//         .focal_mean({radius: 7.5, units: 'meters'})
//         .reduceResolution({
//             reducer: ee.Reducer.mean(),
//             maxPixels: 1024
//         })
//         .reproject({crs: pan_15m.projection(), scale: 30});
    
//     // 2. Injeção de detalhes
//     var fator_detalhe = pan_15m.divide(pan_30m.resample('bicubic').max(0.001));
//     fator_detalhe = fator_detalhe.clamp(0.7, 1.5);
    
//     // 3. Fusão
//     var azul_fusionado = img_azul_30m
//         .resample('bicubic')
//         .multiply(fator_detalhe)
//         .rename('FUS_Blue_15m');
    
//     // 4. Pós-processamento: filtro de borda preservando detalhes
//     var azul_filtrado = azul_fusionado
//         .focal_mean({radius: 7.5, units: 'meters'})
//         .multiply(0.2)
//         .add(azul_fusionado.multiply(0.8));
    
//     return azul_filtrado.clamp(0.01, 0.4);
// }
// ------- FUNÇÃO DEFINITIVA: FUSÃO AZUL + ASTER (CORRIGIDA) -------
function fuseAzulComASTER_Definitivo(img_aster, img_azul_30m) {
    var pan_15m = img_aster.select('AST_NIR_807nm');
    
    // 🔧 CORREÇÃO CRÍTICA: Forçar projeção do Azul para o mesmo CRS do ASTER
    var proj_aster = pan_15m.projection();
    
    // Reprojeta o Azul para o CRS do ASTER (UTM) ANTES de qualquer operação
    var azul_30m_projetado = img_azul_30m.reproject({
        crs: proj_aster,
        scale: 30
    });
    
    // Agora sim, degrada o PAN para 30m no mesmo grid
    var pan_30m = pan_15m
        .reduceResolution({
            reducer: ee.Reducer.mean(),
            maxPixels: 2048,
            bestEffort: true
        })
        .reproject({
            crs: proj_aster,
            scale: 30
        });
    
    // Fator de detalhe (evita divisão por zero)
    var pan_30m_resampled = pan_30m.resample('bicubic');
    var fator_detalhe = pan_15m.divide(pan_30m_resampled.max(0.001));
    
    // Limita o fator para evitar outliers extremos
    fator_detalhe = fator_detalhe.clamp(0.6, 1.8);
    
    // Aplica fusão por multiplicação
    var azul_15m = azul_30m_projetado
        .resample('bicubic')
        .multiply(fator_detalhe)
        .reproject({
            crs: proj_aster,
            scale: 15
        })
        .rename('FUS_Blue_15m');
    
    // Suavização adaptativa (preserva bordas)
    var azul_suave = azul_15m.focal_mean({
        radius: 7.5,
        units: 'meters'
    });
    
    // Mistura 80% fusionado + 20% suavizado (reduz ruído)
    var azul_final = azul_15m.multiply(0.8)
        .add(azul_suave.multiply(0.2))
        .clamp(0.01, 0.45)
        .rename('FUS_Blue_15m');
    
    return azul_final;
}

// ------- VERSÃO ALTERNATIVA (SEM REPROJECT, MAIS EFICIENTE) -------
// Esta versão evita múltiplos reproject e é mais rápida
function fuseAzulComASTER_Eficiente(img_aster, img_azul_30m) {
    var pan_15m = img_aster.select('AST_NIR_807nm');
    
    // 1. Degrada PAN usando kernel Gaussiano (simula PSF)
    var kernel = ee.Kernel.gaussian({
        radius: 15,
        units: 'meters',
        sigma: 7.5
    });
    
    var pan_blur = pan_15m.convolve(kernel);
    
    // 2. Reduz resolução espacialmente (sem mudar projeção ainda)
    var pan_30m_nativo = pan_blur
        .reduceResolution({
            reducer: ee.Reducer.mean(),
            maxPixels: 2048
        });
    
    // 3. Força projeção única para AMBAS as imagens
    var proj_final = pan_15m.projection();
    
    var pan_30m = pan_30m_nativo.reproject({
        crs: proj_final,
        scale: 30
    });
    
    var azul_30m = img_azul_30m.reproject({
        crs: proj_final,
        scale: 30
    });
    
    // 4. Fusão
    var detalhe = pan_15m.subtract(pan_30m.resample('bicubic'));
    var azul_15m = azul_30m.resample('bicubic')
        .add(detalhe.multiply(0.25))  // Peso conservador
        .reproject({
            crs: proj_final,
            scale: 15
        })
        .rename('FUS_Blue_15m');
    
    return azul_15m.clamp(0.01, 0.45);
}

// ------- VERSÃO ROBUSTA (COM FALLBACK PARA GRANDES ÁREAS) -------
function fuseAzulComASTER_Robusto(img_aster, img_azul_30m) {
    try {
        // Tenta método HPF
        return fuseAzulComASTER_Definitivo(img_aster, img_azul_30m);
    } catch (e) {
        // Fallback: Método simples de interpolação
        var proj = img_aster.select('AST_NIR_807nm').projection();
        return img_azul_30m
            .reproject({crs: proj, scale: 15})
            .resample('bicubic')
            .rename('FUS_Blue_15m')
            .clamp(0.01, 0.45);
    }
}

// ------- 7. CÁLCULO DE TEXTURAS NO MOSAICO FINAL (NÃO DURANTE O MOSAICO) -------
function calcularTexturasFinais(img) {
    var nir = img.select('AST_NIR_807nm');
    
    // Converte para Byte para GLCM
    var nir_byte = nir.unitScale(0, 0.35).multiply(255).toByte();
    
    var glcm = nir_byte.glcmTexture({size: 1});
    
    var contrast = glcm.select('AST_NIR_807nm_contrast')
        .unitScale(0, 255)
        .rename('GLCM_Contrast');
        
    var variance = glcm.select('AST_NIR_807nm_var')
        .unitScale(0, 5000)
        .rename('GLCM_Variance');
        
    var entropy = glcm.select('AST_NIR_807nm_ent')
        .unitScale(0, 10)
        .rename('GLCM_Entropy');
    
    return img.addBands([contrast, variance, entropy]);
}

// ------- 8. ÍNDICES GEOLÓGICOS -------
function calcularIndicesGeologia(img) {
    // Al-OH (Caulinita, Alunita)
    var aloh = img.expression(
        "(b('AST_SWIR2_2167nm_AlOH') + b('AST_SWIR4_2263nm_AlMgOH')) / (2 * b('AST_SWIR3_2209nm_Ref'))"
    ).rename('IDX_AlOH_Clay');
    
    // Carbonato/Mg-OH
    var carbonato = img.expression(
        "b('AST_SWIR5_2336nm_CO3') / b('AST_SWIR3_2209nm_Ref')"
    ).rename('IDX_Carbonate');
    
    // Ferric Iron (Fe3+)
    var fe3 = img.expression(
        "b('AST_Red_661nm') / b('AST_Green_556nm')"
    ).rename('IDX_Ferric_Fe');
    
    // NDVI
    var ndvi = img.normalizedDifference(['AST_NIR_807nm', 'AST_Red_661nm']).rename('IDX_NDVI');
    
    // SAVI (Melhor para solo exposto)
    var savi = img.expression(
        "1.5 * (b('AST_NIR_807nm') - b('AST_Red_661nm')) / (b('AST_NIR_807nm') + b('AST_Red_661nm') + 0.5)"
    ).rename('IDX_SAVI');
    
    return img.addBands([aloh, carbonato, fe3, ndvi, savi]);
}

// ------- EXECUÇÃO PRINCIPAL (CORRIGIDA) -------
print("🔄 Processando ASTER (3 anos)...");

var colecao_aster = ee.ImageCollection('ASTER/AST_L1T_003')
    .filterDate(start_date, end_date)
    .filterBounds(area_estudo)
    .select(bandas_aster.concat(['B10']))
    .map(asterTOA)
    .map(mascaraNuvemASTER)
    .map(addQualityASTER);  // Agora funciona!

var mosaico_aster = colecao_aster.qualityMosaic('quality').clip(area_estudo);

// Renomear bandas ASTER
mosaico_aster = mosaico_aster.select(bandas_aster, bandas_nomes_novos);
print("revisado as bandas do mosaio aster ", mosaico_aster);

// 🔧 FORÇAR PROJEÇÃO NO MOSAICO ASTER
var proj_utm = ee.Projection('EPSG:32724');  // UTM 24S (Bahia)
mosaico_aster = mosaico_aster.reproject({
    crs: proj_utm,
    scale: 15
});

print("🔄 Processando Landsat Azul (3 anos)...");
var mosaico_azul_10m = prepararSentinel2Azul();      // Versão Sentinel-2 10m
print("revisando a banda do azul extraida do Landsat ", mosaico_azul_10m);

// revisar ajustes para 10 metros

// 🔧 FORÇAR PROJEÇÃO NO AZUL LANDSAT
mosaico_azul_30m = mosaico_azul_30m.reproject({
    crs: proj_utm,
    scale: 30
});


print("⚡ Realizando Fusão Espacial...");
var azul_15m_final = fuseAzulComASTER_Eficiente(mosaico_aster, mosaico_azul_30m);
print( "revisando as banda do mosaico azul ", azul_15m_final)

// ------- 10. MONTAGEM DO CUBO FINAL -------
var cubo_final_15m = mosaico_aster
    .addBands(azul_15m_final)
    .addBands(mosaico_aster.select('AST_NIR_807nm').rename('Pan_15m'));

// Aplica índices e texturas APÓS o mosaico (mais eficiente)
cubo_final_15m = calcularIndicesGeologia(cubo_final_15m);
cubo_final_15m = calcularTexturasFinais(cubo_final_15m);

print('✅ Bandas Finais (15m):', cubo_final_15m.bandNames());

// ------- 11. VISUALIZAÇÃO -------
Map.centerObject(area_estudo, 9);

// RGB Natural (15m) - AGORA FUNCIONA!
Map.addLayer(cubo_final_15m, {
    bands: ['AST_Red_661nm', 'AST_Green_556nm', 'FUS_Blue_15m'],
    min: 0.05,
    max: 0.35,
    gamma: 1.4
}, 'RGB Natural 15m (ASTER + Landsat Blue)');

// RGB Falsa-Cor Geológica (SWIR)
Map.addLayer(cubo_final_15m, {
    bands: ['AST_SWIR1_1656nm', 'AST_SWIR2_2167nm_AlOH', 'AST_NIR_807nm'],
    min: 0.1,
    max: 0.45,
    gamma: 1.2
}, 'Falsa-Cor Geologia (SWIR1, AlOH, NIR)');

// Índice de Argila
Map.addLayer(cubo_final_15m.select('IDX_AlOH_Clay'), {
    min: 0.8,
    max: 1.4,
    palette: ['blue', 'green', 'yellow', 'red']
}, 'Índice Al-OH (Caulinita)');

// ------- 12. EXPORTAÇÃO -------
var bandas_para_exportar = [
    'AST_Green_556nm',
    'AST_Red_661nm',
    'FUS_Blue_15m',
    'AST_NIR_807nm',
    'AST_SWIR1_1656nm',
    'AST_SWIR2_2167nm_AlOH',
    'AST_SWIR3_2209nm_Ref',
    'AST_SWIR4_2263nm_AlMgOH',
    'AST_SWIR5_2336nm_CO3',
    'AST_SWIR6_2395nm',
    'IDX_AlOH_Clay',
    'IDX_Carbonate',
    'IDX_Ferric_Fe',
    'IDX_NDVI',
    'IDX_SAVI',
    'GLCM_Contrast',
    'GLCM_Variance',
    'GLCM_Entropy'
];

Export.image.toDrive({
    image: cubo_final_15m.select(bandas_para_exportar).toFloat(),
    description: 'ASTER_AZUL_FUSION_15m_' + periodo_str + '_FINAL',
    folder: 'Mosaicos_ASTER_Azul',
    region: area_estudo,
    scale: 15,
    maxPixels: 1e13,
    crs: 'EPSG:32724'
});

print("✅ Script concluído com sucesso!");