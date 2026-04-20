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
// ---------------------------------------------
//    2) Radiância -> reflectância TOA
//      apenas VNIR/SWIR
// --------------------------------------------
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
    // print("numero PI ", Math.PI);
    // print("date image ", ee.Date(imgRad.get('system:time_start')));
    // print("number doy = ", doy);
    // print("number d = ", d);
    // print("lista de bandas ", Object.keys(ESUN));
    var optical = Object.keys(ESUN).map(
        function(b) {
            var esun = ee.Number(ESUN[b]);
            return ee.Image(imgRad).select(b)
                        .multiply(Math.PI)
                        .multiply(d.pow(2))
                        .divide(esun.multiply(cosz))
                        .rename(b);
        }
    );

    return ee.Image(optical)
                .addBands(imgRad.select(['B10','B11','B12','B13','B14']))
                .copyProperties(imgRad, imgRad.propertyNames());
}


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
        .and(swir1.gt(0.15))
        .and(ndvi.lt(0.4));

    // frio relativo por percentil da própria cena
    var tirStats = tir10.reduceRegion({
        reducer: ee.Reducer.percentile([20]),
        geometry: imgToa.geometry(),
        scale: 90,
        maxPixels: 1e8,
        bestEffort: true
    });
    var p20 = ee.Number(tirStats.get('B10'));
    cloud = cloud.and(tir10.lt(p20));
    cloud = cloud.focal_max({radius: 120, units: 'meters'}).rename('cloud');

    // ===== DETECÇÃO DE SOMBRAS COM CORREÇÃO GEOMÉTRICA =====
    // print("show shadow_dist = ", shadow_dist);
    // sombras: escuras e próximas das nuvens
    var dark_pixels = green.lt(0.14).and(nir.lt(0.175));
    // Sombras estão na direção oposta ao sol
    // Usa a distância calculada dinamicamente baseada na elevação solar
    var shadow = dark_pixels.and(
        cloud.focal_max({
            radius: ee.Number(shadow_dist).add(300),
            units: 'meters'
        })
    );

    // Dilatação da sombra (suaviza bordas)
    shadow = shadow.focal_max({radius: 90, units: 'meters'}).rename('shadow');

    // var shadow = dark.and(cloud.focal_max({radius: 1500, units: 'meters'}));    
    // shadow = shadow.focal_max({radius: 90, units: 'meters'}).rename('shadow');
    // Map.addLayer(shadow.selfMask(), {min: 0, max: 1, palette: 'green'}, 'shadow');

    // máscara de pixels válidos
    var valid = imgToa.select(['B01','B02','B3N','B04']).reduce(ee.Reducer.min()).gt(0);
    var clear = valid.and(cloud.not()).and(shadow.not()).rename('clear_mask');

    return imgToa.addBands([ndvi])
                .updateMask(clear)
                copyProperties(imgToa, imgToa.propertyNames());
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
    var colecao_total = l8.merge(l9);
    
    // Estatísticas para debug
    var count = colecao_total.size();
    print('📊 Total de cenas Landsat Azul:', count);
    
    // Mediana temporal (reduz ruído e nuvens residuais)
    var mosaico_azul = colecao_total
                            .median()
                            .clip(area_estudo)
                            .rename('Blue_30m');

    
    // 🔧 PÓS-PROCESSAMENTO: Preenche gaps com interpolação espacial
    // var kernel_suave = ee.Kernel.gaussian({
    //     radius: 60, 
    //     units: 'meters',
    //     sigma: 30
    // });
    
    // // Preenche pixels vazios com média dos vizinhos
    // var azul_preenchido = mosaico_azul
    //     .unmask({
    //         sameFootprint: false,
    //         method: 'nearest'
    //     })
    //     .focal_mean({
    //         kernel: kernel_suave,
    //         iterations: 2
    //     })
    //     .rename('Blue_30m');
    
    return mosaico_azul; //azul_preenchido;
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

// ------- EXECUÇÃO PRINCIPAL (CORRIGIDA) -------
print("🔄 Processando ASTER (3 anos)...");

var colecao_aster = ee.ImageCollection('ASTER/AST_L1T_003')
                        .filterDate(start_date, end_date)
                        .filterBounds(area_estudo)
                        .select(bandas_aster.concat(['B10']))
                        .map(radianceToTOA)
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
var mosaico_azul_30m = prepararLandsatAzul();
print("revisando a banda do azul extraida do Landsat ", mosaico_azul_30m);

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