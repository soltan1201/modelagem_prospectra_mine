import ee
import os
import sys
import collections
collections.Callable = collections.abc.Callable
from pathlib import Path
# pathparent = str(Path(os.getcwd()).parents[0])
# sys.path.append(pathparent)
# from configure_account_projects_ee import get_current_account, get_project_from_account
# projAccount = get_current_account()
# print(f"projetos selecionado >>> {projAccount} <<<")

try:
    ee.Initialize(project= 'ee-solkancengine17')
    print('The Earth Engine package initialized successfully!')
except ee.EEException as e:
    print('The Earth Engine package failed to initialize!')
except:
    print("Unexpected error:", sys.exc_info()[0])
    raise


asset = "projects/mapbiomas-arida/mosaic_aster_p2008"

imgCol = ee.ImageCollection(asset)

lst_id = imgCol.reduceColumns(ee.Reducer.toList(), ['system:index']).get('list').getInfo()
print(imgCol.aggregate_histogram('name_country').getInfo())
for cc, idss in enumerate(lst_id):    
    # id_bacia = idss.split("_")[2]
    path_ = str(asset + '/' + idss)    
    print ("... eliminando ❌ ... item 📍{} : {}  ▶️ ".format(cc, idss))    
    try:
        # ee.data.deleteAsset(path_)
        print(path_)
    except:
        print(" NAO EXISTE!")
