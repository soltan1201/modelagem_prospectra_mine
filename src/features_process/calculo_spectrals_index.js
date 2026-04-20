
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


// 👉 FUNÇÃO PARA CALCULAR ÍNDICES (ANTES DO MOSAICO)
function adicionarIndices(img) {
    var ndvi = img.normalizedDifference(['B3N', 'B02']).rename('NDVI');
    
    // Al-OH (Argilas) - se SWIR disponível
    var aloh = img.expression(
        "(b('B05') + b('B07')) / (2 * b('B06'))",
        {'B05': img.select('B05'), 'B06': img.select('B06'), 'B07': img.select('B07')}
    ).rename('AlOH_Clay');
    
    // Carbonato
    var carbonato = img.expression(
        "b('B08') / b('B06')",
        {'B06': img.select('B06'), 'B08': img.select('B08')}
    ).rename('Carbonate');
    
    // Ferric Iron
    var fe3 = img.expression(
        "b('B02') / b('B01')",
        {'B01': img.select('B01'), 'B02': img.select('B02')}
    ).rename('Ferric_Fe');
    
    return img.addBands([ndvi, aloh, carbonato, fe3]);
}


function addMineralIndices(img) {
    var b1 = img.select('B01');
    var b2 = img.select('B02');
    var b3 = img.select('B3N');
    var b4 = img.select('B04');
    var b5 = img.select('B05');
    var b6 = img.select('B06');
    var b7 = img.select('B07');
    var b8 = img.select('B08');
    var b9 = img.select('B09');
    var b12 = img.select('B12');
    var b13 = img.select('B13');
    var b14 = img.select('B14');

    var ferricIron = b2.divide(b1).rename('ferric_iron_2_1');
    var ferricOx   = b4.divide(b3).rename('ferric_ox_4_3');
    var gossan     = b4.divide(b2).rename('gossan_4_2');
    var ferrousSil = b5.divide(b4).rename('ferrous_sil_5_4');

    var sericite   = b5.add(b7).divide(b6).rename('sericite_illite_smectite');
    var aloh       = b4.add(b6).divide(b5).rename('alunite_kaolinite_pyrophyllite');
    var muscovite  = b7.divide(b6).rename('muscovite_7_6');
    var kaolinite  = b7.divide(b5).rename('kaolinite_7_5');
    var dolomite   = b6.add(b8).divide(b7).rename('dolomite');
    var carbonate  = b13.divide(b14).rename('carbonate_13_14');
    var quartz     = b14.divide(b12).rename('quartz_14_12');

    var sultan_r = b4.divide(b7).rename('sultan_r');
    var sultan_g = b4.divide(b1).rename('sultan_g');
    var sultan_b = b2.divide(b3).multiply(b4.divide(b3)).rename('sultan_b');

    var abrams_r = b4.divide(b7).rename('abrams_r');
    var abrams_g = b4.divide(b3).rename('abrams_g');
    var abrams_b = b2.divide(b1).rename('abrams_b');

    return img.addBands([
        ferricIron, ferricOx, gossan, ferrousSil,
        sericite, aloh, muscovite, kaolinite,
        dolomite, carbonate, quartz,
        sultan_r, sultan_g, sultan_b,
        abrams_r, abrams_g, abrams_b
    ]);
}

// ------- 10. ADICIONAR ÍNDICES AO MOSAICO FINAL -------
function adicionarIndicesGeo(img) {
    // NDVI
    var ndvi = img.normalizedDifference(['B3N', 'B02']).rename('NDVI');
    
    // Al-OH (Argilas)
    var aloh = img.expression(
        "(b('B05') + b('B07')) / (2 * b('B06'))"
    ).rename('AlOH_Clay');
    
    // Carbonato
    var carbonate = img.expression(
        "b('B08') / b('B06')"
    ).rename('Carbonate');
    
    // Ferric Iron
    var fe3 = img.expression(
        "b('B02') / b('B01')"
    ).rename('Ferric_Fe');
    
    // SAVI
    var savi = img.expression(
        "1.5 * (b('B3N') - b('B02')) / (b('B3N') + b('B02') + 0.5)"
    ).rename('SAVI');
    
    return img.addBands([ndvi, aloh, carbonate, fe3, savi]);
}

mosaico_final = adicionarIndicesGeo(mosaico_final);

var mineral = addMineralIndices(imgClean);

// Composições coloridas clássicas
Map.addLayer(
  mineral.select(['sultan_r','sultan_g','sultan_b']),
  {min: 0.6, max: 2.5},
  'Sultan RGB'
);

Map.addLayer(
  mineral.select(['abrams_r','abrams_g','abrams_b']),
  {min: 0.6, max: 2.5},
  'Abrams RGB'
);