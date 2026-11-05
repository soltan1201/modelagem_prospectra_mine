var asset_munico = 'projects/mapbiomas-arida/mine/maxent_amostras_ASTER_mosaico_unico_final_bandas_indices';
var asset_mfinal = 'projects/mapbiomas-arida/mine/maxent_amostras_aster-mosaico-final_bandas_indices';
// var asset = 'projects/mapbiomas-brazil/assets/CARBON/MAPPING/COLLECTION_00/CERRADO_FOREST/SAMPLES/IFN_covariaveis_2023';
var rois_unico = ee.FeatureCollection(asset_munico);
print("mostrar os metadados dos rois ", rois_unico.limit(5));
print('size rois', rois_unico.size());

var rois_final = ee.FeatureCollection(asset_mfinal);
print("mostrar os metadados dos rois rois_final", rois_final.limit(5));
print('size rois rois_final', rois_final.size());

Export.table.toDrive({
  collection: rois_unico, 
  description: 'maxent_amostras_ASTER_mosaico_unico_final_bandas_indices', 
  folder: 'rois_mine'
})

Export.table.toDrive({
  collection: rois_final, 
  description: 'maxent_amostras_aster-mosaico-final_bandas_indices', 
  folder: 'rois_mine'
})