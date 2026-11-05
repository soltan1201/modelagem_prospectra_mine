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
# asset = {'id': 'projects/mapbiomas-arida/Barragem/ROIS'}
# asset = {'id': 'projects/mapbiomas-arida/statistics_col10'}
# asset = {'id': 'projects/mapbiomas-arida/statistics_col9'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/ROIs'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/ROIs2_Rej_SR'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/ROIs2_SR'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/ROIs_Rej_SR'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/ROIs_SR'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/patches'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas/patchesValid'}
asset = {'id': 'projects/mapbiomas-arida/Mapbiomas'}


lstbacias = []  # 7764
lst_years = []
eliminar_files = False
# GetPolygonsfromFolder(asset, 'grade_patch', lstBacias= lstbacias, lstYear= lst_years, play_eliminar= eliminar_files)  # 


# =============================================================================
# Limpeza do projects/mapbiomas-arida/mine para refazer a modelagem (2026-11-02)
#   Grupo A: saídas da modelagem RF/GTB/CNN/AR (todas do mosaico v1)
#   Grupo B: mosaicos intermediários/antigos (fora do v1, seco e harmonizado)
# Ficam: aster-mosaico-final, ASTER_sazonal, ASTER_sazonal_hm_razao,
#   aster-mosaico-final_sazonal_hm_razao, ASTER_MOSAIC_v2, Landsat_sazonal,
#   ASTER_sazonal_cenas, points_sazonal_seco_hm_razao_A/_B, pontos_referencia,
#   grid_irece_5x5*, acervo_cprm, amostras_AR_irece_shm.
# Rodar primeiro com APAGAR_MINE = False (só lista); depois True.
# =============================================================================
MINE = 'projects/mapbiomas-arida/mine/'

# pastas (FOLDER) e coleções (IMAGE_COLLECTION): apaga o conteúdo e depois o contêiner
LST_PASTAS = [
    # --- grupo A: modelagem
    'favorabilidade_ml',           # 57: RF 25 + GTB 25 + 7 CNN (v1)
    'favorabilidade_AR',           # 2: RF/GTB AR irece selecao (v1)
    'cnn_favorabilidade',          # 0 (vazia)
    # --- grupo B: mosaicos antigos e intermediários
    'ASTER_v4',                    # 548 semestrais v4 (entrada do v1)
    'ASTER_0cloud',                # 42 zero-cloud (entrada do v1)
    'ASTER_MOSAIC',                # 50 SWIR estimado do v4 (entrada do v1)
    'ASTER',                       # 27 QualityMosaic antigos
    'ASTER_mosaico_unico',         # 21
    'ASTER_mosaico_unico_base',    # 50
    'ASTER_mosaico_unico_final',   # 25
    'ASTER_sazonal_harmonizado',   # 49 harmonização MODIS (abandonada)
    'ASTER_sazonal_harm_rel',      # 2 2ª passada relativa (abandonada)
    'ASTER_sazonal_hm_residuo',    # 50 modo resíduo (ficou o razão)
    'points',                      # pasta, 20 tabelas samples_v2_*
    'swir_estimated',              # pasta vazia
    'harmonizacao_modis_coef',     # pasta vazia
]

# assets soltos (IMAGE e TABLE)
LST_IMAGES = [
    # --- grupo A: modelagem
    'cluster_fundo_AR_irece',                                    # IMAGE (v1)
    'amostras_AR_irece',                                         # TABLE (v1)
    'maxent_amostras_aster-mosaico-final_bandas_indices',        # TABLE (v1)
    'maxent_amostras_ASTER_mosaico_unico_final_bandas_indices',  # TABLE
    'maxent_fosfato_amostras_teste_zero_cloud_2000_2008',        # TABLE
    # --- grupo B: pontos antigos
    'points_zeroCC_00_08',                                       # TABLE (SWIR do v1)
    'points_zeroCC_09_20',                                       # TABLE (SWIR do v1)
    'pontos_fosfato',                                            # TABLE (-> pontos_referencia)
]

# nada daqui pode ser apagado; o script para se aparecer numa das listas
PROTEGIDOS = {
    'aster-mosaico-final', 'ASTER_sazonal', 'ASTER_sazonal_hm_razao',
    'aster-mosaico-final_sazonal_hm_razao', 'ASTER_MOSAIC_v2', 'Landsat_sazonal',
    'ASTER_sazonal_cenas', 'points_sazonal_seco_hm_razao_A',
    'points_sazonal_seco_hm_razao_B', 'pontos_referencia', 'grid_irece_5x5',
    'grid_irece_5x5_solape100', 'acervo_cprm', 'amostras_AR_irece_shm',
}

APAGAR_MINE = False   # False = só lista (rodada de teste); True = apaga


def listar_filhos(parent):
    """Todos os filhos de uma pasta/coleção, com paginação (ASTER_v4 tem 548)."""
    filhos, token = [], None
    while True:
        params = {'parent': parent, 'pageSize': 1000}
        if token:
            params['pageToken'] = token
        resp = ee.data.listAssets(params)
        filhos += resp.get('assets', [])
        token = resp.get('nextPageToken')
        if not token:
            return filhos


def apagar(path, tipo, apagar_de_fato, erros):
    print(f"   {'apagando' if apagar_de_fato else '[teste]'} {tipo:16s} {path}")
    if apagar_de_fato:
        try:
            ee.data.deleteAsset(path)
        except ee.EEException as e:
            print(f"      ERRO: {e}")
            erros.append(path)


def apagar_pasta(path, apagar_de_fato, erros):
    """Apaga o conteúdo (recursivo) e depois a própria pasta/coleção."""
    for f in listar_filhos(path):
        if f['type'] in ('FOLDER', 'IMAGE_COLLECTION'):
            apagar_pasta(f['name'], apagar_de_fato, erros)
        else:
            apagar(f['name'], f['type'], apagar_de_fato, erros)
    apagar(path, 'CONTÊINER', apagar_de_fato, erros)


def limpar_mine(apagar_de_fato=False):
    conflito = PROTEGIDOS & (set(LST_PASTAS) | set(LST_IMAGES))
    if conflito:
        sys.exit(f"PARADO: assets protegidos nas listas de apagar: {sorted(conflito)}")

    existentes = {a['name'].split('/')[-1]: a['type']
                  for a in listar_filhos(MINE.rstrip('/'))}
    erros, n = [], 0
    print(f"\n=== pastas/coleções ({len(LST_PASTAS)}) ===")
    for nome in LST_PASTAS:
        if nome not in existentes:
            print(f" -- {nome}: não existe, pulando")
            continue
        filhos = listar_filhos(MINE + nome)
        print(f" >> {nome} ({existentes[nome]}, {len(filhos)} itens)")
        n += len(filhos) + 1
        apagar_pasta(MINE + nome, apagar_de_fato, erros)

    print(f"\n=== imagens/tabelas ({len(LST_IMAGES)}) ===")
    for nome in LST_IMAGES:
        if nome not in existentes:
            print(f" -- {nome}: não existe, pulando")
            continue
        n += 1
        apagar(MINE + nome, existentes[nome], apagar_de_fato, erros)

    print(f"\n{n} assets {'apagados' if apagar_de_fato else 'seriam apagados (teste)'}"
          f"; {len(erros)} erros")
    for e in erros:
        print(f"   não apagado: {e}")


limpar_mine(apagar_de_fato=APAGAR_MINE)

