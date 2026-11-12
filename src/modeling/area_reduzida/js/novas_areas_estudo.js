
var area_IRECE = ee.Geometry.Rectangle([-42.31762892604435,-12.359936623152336,
                                  -41.2382221877631,-10.5375816299746]);
var area_JuazeiroCoite = ee.Geometry.Rectangle([-40.8409898019856,-12.359143050738856,
                                  -39.3688218332356,-9.164831060671531]);

Map.addLayer(area_IRECE, {color: 'red'}, 'irece');
Map.addLayer(area_JuazeiroCoite, {color: 'green'}, 'Juazeiro');


Export.table.toDrive({
  collection: area_IRECE, 
  description: 'regiao_Juazeiro_coite', 
  folder: 'shps', 
  fileFormat: "SHP"
});


Export.table.toDrive({
  collection: area_JuazeiroCoite, 
  description: 'regiao_IRECE', 
  folder: 'shps', 
  fileFormat: "SHP"
});
