"""
Amostragem dos mosaicos ASTER semestrais.

Amostragem na grade de 15 m em que os mosaicos foram exportados
(Export.image.toAsset com scale 15: o asset guarda todas as bandas em 15 m —
SWIR 30 m e TIR 90 m já vêm reamostrados por vizinho mais próximo no export).

Para cada imagem:
  1. Coleta N_PONTOS aleatórios sobre todos os pixels válidos, a 15 m
  2. Exporta CSV com todas as bandas + quality (para filtrar depois)
"""

import ee
import os
import sys
import time
import logging

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# ------- PARÂMETROS -------
ASSET_ID       = 'projects/mapbiomas-arida/mosaic_aster_v2'
PASTA_DRIVE    = 'ASTER_Samples'
ASSET_POINTS   = 'projects/mapbiomas-arida/mine/points'
N_PONTOS       = 30000
ESCALA_AMOST   = 15    # metros — grade nativa do asset do mosaico
TILE_SCALE     = 4     # fator de particionamento interno do GEE
SEED           = 42
INCLUIR_COORDS = True  # True → adiciona colunas latitude/longitude no CSV

# Destino de exportação: 'drive' | 'asset' | 'both'
EXPORTAR_PARA  = 'both'

# ------- INICIALIZAÇÃO -------
from pathlib import Path
pathparent = str(Path(os.getcwd()).parents[0])
sys.path.append(pathparent)
from configure_account_projects_ee import get_current_account, get_project_from_account
projAccount = get_current_account()
print(f"projetos selecionado >>> {projAccount} <<<")
try:
    ee.Initialize(project=projAccount)
    log.info('Earth Engine inicializado com sucesso.')
except Exception as e:
    log.error(f'Erro de inicialização: {e}')
    raise

# ------- CARREGAR COLEÇÃO -------
colecao = ee.ImageCollection(ASSET_ID)
n_total = colecao.size().getInfo()
print(f'Total de imagens no asset: {n_total}')

lista_imgs = colecao.toList(n_total)


# ------- FUNÇÕES DE EXPORTAÇÃO -------

def exportar_para_drive(amostras, sys_index):
    """Exporta FeatureCollection como CSV no Google Drive."""
    task = ee.batch.Export.table.toDrive(
        collection=amostras,
        description=f'sample_v2_{sys_index}',
        folder=PASTA_DRIVE,
        fileNamePrefix=sys_index,
        fileFormat='CSV'
    )
    task.start()
    print(f'  [Drive] Task iniciada → {PASTA_DRIVE}/{sys_index}.csv')
    return task


def exportar_para_asset(amostras, sys_index):
    """Exporta FeatureCollection como asset GEE (FeatureCollection)."""
    asset_id = f'{ASSET_POINTS}/samples_v2_{sys_index}'
    task = ee.batch.Export.table.toAsset(
        collection=amostras,
        description=f'samples_v2_{sys_index}',
        assetId=asset_id
    )
    task.start()
    print(f'  [Asset] Task iniciada → {asset_id}')
    return task


# ------- PROCESSAR CADA IMAGEM -------
tasks = []

for i in range(n_total):
    img = ee.Image(lista_imgs.get(i))
    sys_index = img.get('system:index').getInfo()
    print(f'\n[{i+1}/{n_total}] {sys_index}')

    regiao = img.geometry()

    # Coleta N_PONTOS aleatórios sobre todos os pixels válidos (grade de 15 m)
    amostras = img.sample(
        region=regiao,
        scale=ESCALA_AMOST,
        numPixels=N_PONTOS,
        seed=SEED,
        geometries=INCLUIR_COORDS,
        tileScale=TILE_SCALE
    )

    # Exporta conforme destino configurado
    if EXPORTAR_PARA in ('drive', 'both'):
        t = exportar_para_drive(amostras, sys_index)
        tasks.append({'index': sys_index, 'destino': 'drive', 'task': t})

    if EXPORTAR_PARA in ('asset', 'both'):
        t = exportar_para_asset(amostras, sys_index)
        tasks.append({'index': sys_index, 'destino': 'asset', 'task': t})

    time.sleep(0.3)  # evita burst de requisições

# ------- RESUMO -------
print(f'\n{"=" * 55}')
print(f'Tasks iniciadas: {len(tasks)} / {n_total}')
print('Acompanhe em: https://code.earthengine.google.com/tasks')
print('=' * 55)

# ------- MONITORAMENTO OPCIONAL -------
# Descomente abaixo para aguardar todas as tasks e exibir o status final.

# import sys
#
# print('\nAguardando conclusão das tasks...')
# while True:
#     pendentes = [t for t in tasks if t['task'].status()['state'] in ('READY', 'RUNNING')]
#     if not pendentes:
#         break
#     print(f'  {len(pendentes)} task(s) em execução...')
#     time.sleep(60)
#
# print('\nStatus final:')
# for t in tasks:
#     estado = t['task'].status()['state']
#     print(f'  {t["index"]:50s}  {estado}')
