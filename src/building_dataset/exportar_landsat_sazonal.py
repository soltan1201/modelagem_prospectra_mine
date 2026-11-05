"""
Mosaico Landsat de referência para a harmonização dos mosaicos sazonais secos ASTER
(harmonizar_sazonal_landsat.py). Um asset por período × célula da grid 5 × 5 (tasks
pequenas, rodam em paralelo), na grade de 30 m ancorada no canto NW do bloco (2 × o
pixel ASTER de 15 m): as células encaixam pixel a pixel.

Mesma janela e mesmos meses do mosaico ASTER (build_mosaico_sazonal_seco.CONFIG e
MESES_SECO = jun–out), mediana:
  green, red, nir, swir1, swir2  LANDSAT/COMPOSITES/C02/T1_L2_32DAY (reflectância de
                                 superfície), × 1000 → mesma escala do mosaico ASTER
  st                             temperatura de superfície Collection 2 L2 (Landsat 5/7 no
                                 A; 5/7/8/9 no B), sem nuvem/sombra (QA_PIXEL), °C × 100
  n_ls, n_st                     nº de composições / cenas usadas no pixel
O período B só tem green, red, nir e st (o ASTER B não tem SWIR).

Saída: COL_LS/Landsat_sazonal_<P>_seco_mediana_<LxCy>, 30 m, INT16, propriedades
periodo e celula_nome (o visualizador junta as células com mosaic()).

Conta: a mesma da harmonização (caatinga04, projeto mapbiomas-brazil).

Uso: python exportar_landsat_sazonal.py
"""

import logging

import ee

import build_mosaico_sazonal_seco as s
import harmonizar_sazonal_modis as h
from build_mosaicos_semestrais_v2 import (CRS_EXPORT, NORTE_BLOCO, OESTE_BLOCO, RES_GRAUS,
                                          ler_grid, listar_destino, perguntar_celulas,
                                          perguntar_opcao, perguntar_sim_nao)
from gee_tools import switch_user

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
COL_LS    = f'{h.PASTA_MINE}/Landsat_sazonal'
COMPOSTO  = 'LANDSAT/COMPOSITES/C02/T1_L2_32DAY'
BANDAS_LS = {'A': ['green', 'red', 'nir', 'swir1', 'swir2'],
             'B': ['green', 'red', 'nir']}
# grade de 30 m = 2 × a grade de 15 m dos mosaicos ASTER, mesmo canto NW
CRS_TRANSFORM_30 = [2 * RES_GRAUS, 0, OESTE_BLOCO, 0, -2 * RES_GRAUS, NORTE_BLOCO]
ESCALA_REFL = s.ESCALA_REFL            # 1000

# Temperatura de superfície (Collection 2 L2): coleção → banda ST
COLS_ST = {'LANDSAT/LT05/C02/T1_L2': 'ST_B6',
           'LANDSAT/LE07/C02/T1_L2': 'ST_B6',
           'LANDSAT/LC08/C02/T1_L2': 'ST_B10',
           'LANDSAT/LC09/C02/T1_L2': 'ST_B10'}
ST_GANHO, ST_OFFSET = 0.00341802, 149.0          # DN → kelvin
BITS_RUINS = (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4)   # nuvem dilatada, cirrus, nuvem, sombra


def nome_ls(periodo: str, celula: str) -> str:
    return f'Landsat_sazonal_{periodo}_seco_mediana_{celula}'


# ═══════════════════════════════════════════════════════════════════════════
def refletancia(periodo: str, regiao) -> tuple[ee.Image, ee.Image]:
    cfg = s.CONFIG[periodo]
    col = (ee.ImageCollection(COMPOSTO)
             .filterDate(cfg['data_inicio'], cfg['data_fim'])
             .filter(s.filtro_meses(s.MESES_SECO))
             .filterBounds(regiao)
             .select(BANDAS_LS[periodo]))
    return col.median(), col.select('red').count().rename('n_ls')


def checar_escala(med: ee.Image, regiao) -> float:
    """A composição deve vir em reflectância 0–1. Se vier em DN escalado, para aqui."""
    v = (med.select('red').reduceRegion(ee.Reducer.median(), regiao.centroid(1).buffer(20000),
                                        300, bestEffort=True).get('red').getInfo())
    if v is None or not 0 < v < 1.5:
        raise SystemExit(f'Vermelho mediano da composição = {v}: não parece reflectância 0–1. '
                         'Confira a escala de LANDSAT/COMPOSITES/C02/T1_L2_32DAY antes de exportar.')
    return v


def temperatura(periodo: str, regiao) -> tuple[ee.Image, ee.Image]:
    cfg = s.CONFIG[periodo]
    col = None
    for cid, banda in COLS_ST.items():
        c = (ee.ImageCollection(cid)
               .filterDate(cfg['data_inicio'], cfg['data_fim'])
               .filter(s.filtro_meses(s.MESES_SECO))
               .filterBounds(regiao)
               .map(lambda i, b=banda: i.select(b).multiply(ST_GANHO).add(ST_OFFSET - 273.15)
                    .rename('st')
                    .updateMask(i.select(b).gt(0))
                    .updateMask(i.select('QA_PIXEL').bitwiseAnd(BITS_RUINS).eq(0))))
        col = c if col is None else col.merge(c)
    return col.median(), col.count().rename('n_st')


def mosaico_landsat(periodo: str, celula: str, regiao) -> ee.Image:
    cfg = s.CONFIG[periodo]
    refl, n_ls = refletancia(periodo, regiao)
    st, n_st = temperatura(periodo, regiao)
    img = (refl.multiply(ESCALA_REFL).round().toInt16()
              .addBands(st.multiply(100).round().toInt16())
              .addBands(n_ls.toInt16()).addBands(n_st.toInt16()))
    return img.clip(regiao).set({
        'periodo': periodo, 'celula_nome': celula,
        'data_inicio': cfg['data_inicio'], 'data_fim': cfg['data_fim'],
        'meses': ','.join(map(str, s.MESES_SECO)),
        'refletancia': f'{COMPOSTO} mediana, reflectância de superfície × {ESCALA_REFL}',
        'st': 'Collection 2 L2 ST (L5/L7/L8/L9), mediana sem nuvem, °C × 100',
        'script': 'exportar_landsat_sazonal.py'})


def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    switch_user(h.CONTA)
    ee.Initialize(project=h.PROJETO)
    log.info(f'Earth Engine com {h.CONTA} (projeto {h.PROJETO})')
    grid = ler_grid()

    print('\n══════ Mosaico Landsat de referência (meses secos), por célula ══════')
    periodos = perguntar_opcao('Períodos:', [('A e B', ['A', 'B']), ('só A', ['A']), ('só B', ['B'])])
    celulas = perguntar_celulas(grid)
    feitos = listar_destino(COL_LS, criar=True)
    ativas = h.tasks_ativas()
    fila = []
    for p in periodos:
        for cid in celulas:
            c = grid[cid]['nome']
            if nome_ls(p, c) in feitos:
                log.info(f'  {p} {c}: pulado (já existe; apague para refazer)')
            elif f'landsat_{p}_{c}' in ativas:
                log.info(f'  {p} {c}: pulado (processando)')
            else:
                fila.append((p, cid, c))
    if not fila or not perguntar_sim_nao(f'Exportar {len(fila)} células Landsat?', padrao=True):
        return
    # escala conferida uma vez por período, no centro da 1ª célula da fila
    for p in sorted({p for p, _, _ in fila}):
        cid0 = next(cid for q, cid, _ in fila if q == p)
        v = checar_escala(refletancia(p, grid[cid0]['geom'])[0], grid[cid0]['geom'])
        log.info(f'  {p}: vermelho mediano = {v:.3f} (reflectância 0–1, ok)')
    for p, cid, c in fila:
        geom = grid[cid]['geom']
        task = ee.batch.Export.image.toAsset(
            image=mosaico_landsat(p, c, geom), description=f'landsat_{p}_{c}',
            assetId=f'{COL_LS}/{nome_ls(p, c)}', region=geom, crs=CRS_EXPORT,
            crsTransform=CRS_TRANSFORM_30, maxPixels=1e13,
            pyramidingPolicy={'.default': 'mean', 'n_ls': 'mode', 'n_st': 'mode'})
        task.start()
        log.info(f'  task {task.id} → {COL_LS}/{nome_ls(p, c)}')


if __name__ == '__main__':
    main()
