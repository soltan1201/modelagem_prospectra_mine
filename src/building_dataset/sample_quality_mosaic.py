"""
Amostragem dos mosaicos ASTER semestrais.

Todas as bandas são harmonizadas para 30 m antes da amostragem:
  - VNIR (B01, B02, B3N)  15 m → média de 4 pixels (2×2)
  - SWIR (B04–B09)        30 m → sem alteração
  - TIR  (B10–B14)        90 m → interpolação bilinear (1 pixel → 3×3 sub-pixels)

Para cada imagem:
  1. Harmoniza todas as bandas para 30 m (inclui banda quality)
  2. Coleta N_PONTOS aleatórios sobre todos os pixels válidos
  3. Exporta CSV com todas as bandas + quality (para filtrar depois)
"""

import ee
import time
import logging

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# ------- PARÂMETROS -------
ASSET_ID       = 'projects/mapbiomas-arida/mosaic_aster'
PASTA_DRIVE    = 'ASTER_Samples'
ASSET_POINTS   = 'projects/mapbiomas-arida/mine/points'
N_PONTOS       = 30000
ESCALA_AMOST   = 30    # metros — resolução após harmonização
TILE_SCALE     = 4     # fator de particionamento interno do GEE
SEED           = 42
INCLUIR_COORDS = True  # True → adiciona colunas latitude/longitude no CSV

# Destino de exportação: 'drive' | 'asset' | 'both'
EXPORTAR_PARA  = 'both'

# ------- INICIALIZAÇÃO -------
projAccount ='mapbiomas-caatinga-cloud02'
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


# ------- FUNÇÃO: HARMONIZAR TODAS AS BANDAS PARA 30m -------
def harmonizar_para_30m(img):
    """
    VNIR 15m → 30m : reduceResolution com média (agrupa 2×2 = 4 pixels)
    SWIR 30m        : sem alteração (resolução nativa)
    TIR  90m → 30m : resample bilinear (interpola 1 pixel em 3×3 sub-pixels)

    A máscara de qualidade deve ser aplicada ANTES desta função,
    assim a média do VNIR exclui automaticamente pixels inválidos.
    """
    # Projeção de referência: SWIR nativo a 30m
    proj_30m = img.select('B04').projection()

    # VNIR: agrega 4 pixels de 15m em 1 pixel de 30m por média
    vnir = (
        img.select(['B01', 'B02', 'B3N'])
           .reduceResolution(
               reducer=ee.Reducer.mean(),
               bestEffort=False,
               maxPixels=4        # exige exatamente 2×2 pixels por célula
           )
           .reproject(proj_30m)
    )

    # SWIR: resolução nativa, sem processamento adicional
    swir = img.select(['B04', 'B05', 'B06', 'B07', 'B08', 'B09'])

    # TIR: upsample de 90m → 30m por interpolação bilinear
    tir = (
        img.select(['B10', 'B11', 'B12', 'B13', 'B14'])
           .resample('bilinear')
           .reproject(proj_30m)
    )

    # quality: mesma resolução do VNIR (15m) → agrega para 30m por média
    quality = (
        img.select(['quality'])
           .reduceResolution(reducer=ee.Reducer.mean(), bestEffort=False, maxPixels=4)
           .reproject(proj_30m)
    )

    return ee.Image(
        vnir.addBands(swir)
            .addBands(tir)
            .addBands(quality)
            .copyProperties(img, img.propertyNames())
    )


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

    # Harmoniza para 30m (VNIR: média 2×2 | TIR: bilinear) — inclui banda quality
    img_30m = harmonizar_para_30m(img)

    # Coleta N_PONTOS aleatórios sobre todos os pixels válidos
    amostras = img_30m.sample(
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
