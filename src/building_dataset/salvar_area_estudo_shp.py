"""
Salva a área de estudo (retângulo usado em select_black_list_save_mosaicSem.js)
como shapefile em EPSG:4326.

Uso:
    python salvar_area_estudo_shp.py [pasta_saida]

Requer: geopandas, shapely
"""
import sys
from pathlib import Path

import geopandas as gpd
from shapely.geometry import Polygon

list_coord = [
    [-44.460483519818084, -14.790793773956535],
    [-39.203525512005584, -14.790793773956535],
    [-39.203525512005584, -9.416525285369607],
    [-44.460483519818084, -9.416525285369607],
    [-44.460483519818084, -14.790793773956535],
]

pasta_saida = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'data'
pasta_saida.mkdir(parents=True, exist_ok=True)
arquivo_shp = pasta_saida / 'area_estudo.shp'

area_estudo = gpd.GeoDataFrame(
    {'id': [1], 'nome': ['area_estudo']},
    geometry=[Polygon(list_coord)],
    crs='EPSG:4326',
)
area_estudo.to_file(arquivo_shp, driver='ESRI Shapefile')

print(f'Shapefile salvo em: {arquivo_shp}')
print(f'Bounds: {area_estudo.total_bounds}')
