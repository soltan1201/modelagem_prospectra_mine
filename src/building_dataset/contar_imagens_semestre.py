"""
Conta as cenas ASTER L1T (ASTER/AST_L1T_003) que cruzam a área de estudo, por
ano e semestre, e salva em CSV nesta pasta.

Colunas (uma linha por ano × semestre; semestre 1 = jan–jun, 2 = jul–dez)
  total          cenas que cruzam a área
  diurnas        sol acima do horizonte e com ganho VNIR (as noturnas só têm TIR)
  noturnas       total − diurnas
  com_swir       diurnas com GAIN_COEFFICIENT_B04 (SWIR real; parou em 2008)
  cc_lt_max      diurnas com CLOUDCOVER < CC_MAX (mesmo filtro dos mosaicos)
  apos_blacklist cc_lt_max fora da blacklist (lista_imagens_Notaceites.txt)
  zero_cloud     apos_blacklist que estão na list_zero_cloud.txt
  cc_media       CLOUDCOVER médio das diurnas
  cc_mediana     CLOUDCOVER mediano das diurnas
São os mesmos filtros de montar_semestre (build_mosaicos_semestrais.py); a
coluna apos_blacklist é o número de cenas que entra no mosaico semestral.

Com POR_CELULA = True conta também por célula da grid 5×5 (arquivo separado).

Saída: contagem_imagens_aster_semestre.csv (e _por_celula.csv) nesta pasta.

Uso: python contar_imagens_semestre.py
"""

import logging
from pathlib import Path

import ee
import pandas as pd

import build_mosaicos_semestrais as b

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

PROJETO_EE = 'mapbiomas-brazil'
ANO_INI, ANO_FIM = 2000, 2020      # período dos mosaicos do projeto
POR_CELULA = False                 # True: também por célula (25 × 42 linhas)
PASTA = Path(__file__).resolve().parent
SAIDA = PASTA / 'contagem_imagens_aster_semestre.csv'


def contar(area, ano, semestre, blacklist, zero_cloud):
    """Dicionário (servidor) com as contagens de um semestre numa área."""
    inicio = ee.Date.fromYMD(ano, 1 if semestre == 1 else 7, 1)
    base = ee.ImageCollection('ASTER/AST_L1T_003').filterDate(inicio, inicio.advance(6, 'month')).filterBounds(area)
    diurnas = base.filter(ee.Filter.And(ee.Filter.gt('SOLAR_ELEVATION', 0),
                                        ee.Filter.notNull(['GAIN_COEFFICIENT_B01'])))
    # a B04 existe em todas as cenas (vazia depois de 2008-S1); o SWIR real tem ganho
    com_swir = diurnas.filter(ee.Filter.notNull(['GAIN_COEFFICIENT_B04']))
    cc = diurnas.filter(ee.Filter.lt('CLOUDCOVER', b.CC_MAX))
    filtrada = cc.filter(ee.Filter.inList('system:index', blacklist).Not())
    zero = filtrada.filter(ee.Filter.inList('system:index', zero_cloud))
    return ee.Dictionary({
        'ano': ano, 'semestre': semestre,
        'total': base.size(), 'diurnas': diurnas.size(), 'com_swir': com_swir.size(),
        'cc_lt_max': cc.size(), 'apos_blacklist': filtrada.size(), 'zero_cloud': zero.size(),
        'cc_media': ee.Algorithms.If(diurnas.size().gt(0), diurnas.aggregate_mean('CLOUDCOVER'), None),
        'cc_mediana': ee.Algorithms.If(diurnas.size().gt(0),
                                       diurnas.reduceColumns(ee.Reducer.median(), ['CLOUDCOVER']).get('median'), None),
    })


def tabela(area, blacklist, zero_cloud, extra=None):
    semestres = [(a, s) for a in range(ANO_INI, ANO_FIM + 1) for s in (1, 2)]
    linhas = ee.List([contar(area, a, s, blacklist, zero_cloud) for a, s in semestres]).getInfo()
    df = pd.DataFrame(linhas)
    df['noturnas'] = df.total - df.diurnas
    for k, v in (extra or {}).items():
        df.insert(0, k, v)
    return df


def main():
    ee.Initialize(project=PROJETO_EE)
    ids_black = b.ler_ids('lista_imagens_Notaceites.txt')
    ids_zero = b.ler_ids('list_zero_cloud.txt') - ids_black      # blacklist vence, como nos mosaicos
    blacklist, zero_cloud = ee.List(sorted(ids_black)), ee.List(sorted(ids_zero))
    log.info(f'blacklist {len(ids_black)} cenas · zero-cloud {len(ids_zero)} cenas · {ANO_INI}–{ANO_FIM}')

    area = ee.Geometry.Polygon([b.LIST_COORD])
    colunas = ['ano', 'semestre', 'total', 'diurnas', 'noturnas', 'com_swir', 'cc_lt_max',
               'apos_blacklist', 'zero_cloud', 'cc_media', 'cc_mediana']
    df = tabela(area, blacklist, zero_cloud)[colunas]
    df.to_csv(SAIDA, index=False)
    log.info(f'área de estudo: {len(df)} semestres → {SAIDA}')
    print(df.to_string(index=False))
    print('\nTotais:', df[['total', 'diurnas', 'com_swir', 'cc_lt_max', 'apos_blacklist', 'zero_cloud']]
          .sum().to_dict())

    if POR_CELULA:
        grid = b.ler_grid()
        partes = []
        for cid, c in sorted(grid.items()):
            partes.append(tabela(c['geom'], blacklist, zero_cloud, {'celula_id': cid, 'celula': c['nome']}))
            log.info(f'  célula {c["nome"]}')
        dc = pd.concat(partes, ignore_index=True)[['celula_id', 'celula'] + colunas]
        arq = SAIDA.with_name(SAIDA.stem + '_por_celula.csv')
        dc.to_csv(arq, index=False)
        log.info(f'por célula: {len(dc)} linhas → {arq}')


if __name__ == '__main__':
    main()
