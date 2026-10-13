import ee
import os
import sys
import collections
collections.Callable = collections.abc.Callable

from pathlib import Path
pathparent = str(Path(os.getcwd()).parents[0])
sys.path.append(pathparent)
from configure_account_projects_ee import get_current_account, get_project_from_account
projAccount = get_current_account()
print(f"projetos selecionado >>> {projAccount} <<<")

try:
    ee.Initialize(project= projAccount)
    print('The Earth Engine package initialized successfully!')
except ee.EEException as e:
    print('The Earth Engine package failed to initialize!')
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise

def GetPolygonsfromFolder(assetFolder, sufixo, lstBacias= [], lstYear= [], play_eliminar= False):
    getlistPtos = ee.data.getList(assetFolder)
    lst_path = []
    sizeFiles = len(getlistPtos)

    for cc, idAsset in enumerate(getlistPtos): 
        path_ = idAsset.get('id') 
        name =  path_.split("/")[-1]        
        
        if len(lstBacias) > 0 and len(lstYear) > 0:    
            idBacia = name.split('_')[0]
            nyear = int(name.split('_')[1])        
            if idBacia in lstBacias and nyear in lstYear:
                print(" --- passo nas condicionais --- ")
                print(f' {idBacia}    {nyear}    {name}'   )
                # print(path_)
                lst_path.append(path_)
        else:
            if sufixo in str(name): 
                lst_path.append(path_)        
     
            # print(name)
        # if str(name).startswith(sufixo): AMOSTRAS/col7/CAATINGA/classificationV
    cc = 0
    sizeFiles = len(lst_path)
    lst_erros = []
    for npath in lst_path:
        name = npath.split("/")[-1]
        if len(name) > 0:
            print(f"eliminando {cc}/{sizeFiles}:  {name}")
            print(npath)
            if play_eliminar:
                try:
                    ee.data.deleteAsset(npath)
                except ee.EEException as e:
                    # asset já apagado ou conta sem permissão de apagar: registra e segue
                    print(f"   ERRO ao eliminar {name}: {e}")
                    lst_erros.append(npath)

        cc += 1

    print(lstBacias)
    if lst_erros:
        print(f"{len(lst_erros)}/{sizeFiles} assets NÃO foram eliminados (ver ERRO acima)")


# asset = {'id': 'projects/mapbiomas-arida/mine/points'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/agosto2022'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/SHP/others'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/SHP/shp'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/SHP/shp2022'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/SHP/shp2023'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/SHP/shp2024'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/SHP/shp2025'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/ROIs/2021ROIs'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/ROIs/2022ROIs'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/ROIs/2023ROIs'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/shps'}
# asset = {'id': 'projects/mapbiomas-arida/ALERTAS/statist'}
asset = {'id': 'projects/mapbiomas-arida/Barragem/ROIS'}

lstbacias = []  # 7764
lst_years = []
eliminar_files = True
GetPolygonsfromFolder(asset, '', lstBacias= lstbacias, lstYear= lst_years, play_eliminar= eliminar_files)  # 

