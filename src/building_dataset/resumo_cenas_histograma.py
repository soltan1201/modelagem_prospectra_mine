"""
Resumo das cenas ASTER L1T "boas" × "ruins" — versão Python de
analise_histogramas_cenas.js (seção 4, tabela-resumo).

Para cada cena dos pares, busca no GEE e grava numa linha do CSV:
  - data, elevação/azimute solar, CLOUDCOVER;
  - todas as propriedades com "GAIN" no nome (ganho de cada banda);
  - percentis P2/P50/P98 de B01, B02, B3N em DN e em reflectância TOA;
  - fração de pixels saturados (DN ≥ 254) por banda;
  - área da cena e área em comum com a outra cena do par.

Também salva a lista completa de propriedades de cada cena (JSON), para
conferir nomes de metadados que não estejam no CSV.

Uso:
  python resumo_cenas_histograma.py
  python resumo_cenas_histograma.py --project outro-projeto-gee

Saídas (src/Dados/histograma/):
  resumo_cenas_histograma.csv
  propriedades_cenas.json
"""

import argparse
import json
import logging
from pathlib import Path

import ee
import pandas as pd

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# ------- PARÂMETROS -------
PARES = [
    {'nome': '2003', 'bom': 'ASTER/AST_L1T_003/20030914131811', 'ruim': 'ASTER/AST_L1T_003/20031011130046'},
    {'nome': '2001', 'bom': 'ASTER/AST_L1T_003/20010816131957', 'ruim': 'ASTER/AST_L1T_003/20010928130033'},
]

BANDAS_VNIR = ['B01', 'B02', 'B3N']
ESCALA      = 60      # m — mesma escala do script JS
DN_SAT      = 254     # DN ≥ 254 = saturado (VNIR 8 bits)
ESUN        = {'B01': 1848, 'B02': 1549, 'B3N': 1114}   # mesmos valores do pipeline

PASTA_SAIDA = Path(__file__).resolve().parents[1] / 'Dados' / 'histograma'

parser = argparse.ArgumentParser(description='Resumo das cenas ASTER boas × ruins.')
parser.add_argument('--project', default='mapbiomas-caatinga-cloud02',
                    help='projeto GEE usado no ee.Initialize')
args = parser.parse_args()

# ------- INICIALIZAÇÃO -------
try:
    ee.Initialize(project= 'ee-solkancengine17')
    log.info(f'Earth Engine inicializado (projeto {args.project}).')
except Exception as e:
    log.error(f'Erro de inicialização: {e}')
    raise

PASTA_SAIDA.mkdir(parents=True, exist_ok=True)


# ------- FUNÇÕES (equivalentes às do script JS) -------
def dn_para_toa(img: ee.Image) -> ee.Image:
    """DN → radiância [(DN − 1) × GAIN_COEFFICIENT] → reflectância TOA, só VNIR."""
    elev = ee.Number(img.get('SOLAR_ELEVATION'))
    cosz = ee.Number(90).subtract(elev).multiply(3.141592653589793 / 180).cos()
    doy  = ee.Date(img.get('system:time_start')).getRelative('day', 'year').add(1)
    d    = ee.Number(1).subtract(
        ee.Number(0.01672).multiply(doy.subtract(4).multiply(0.9856 * 3.141592653589793 / 180).cos()))

    bandas = []
    for b in BANDAS_VNIR:
        dn = img.select(b)
        ganho = ee.Number(img.get(f'GAIN_COEFFICIENT_{b}'))
        toa = (dn.subtract(1).multiply(ganho)
                 .multiply(3.141592653589793).multiply(d.pow(2))
                 .divide(cosz.multiply(ESUN[b]))
                 .updateMask(dn.gt(0))
                 .rename(b))
        bandas.append(toa)
    return ee.Image.cat(bandas)


def reduzir(img: ee.Image, reducer: ee.Reducer, regiao: ee.Geometry) -> ee.Dictionary:
    return img.reduceRegion(reducer=reducer, geometry=regiao, scale=ESCALA,
                            maxPixels=1e10, bestEffort=True)


def prefixar(dic: ee.Dictionary, prefixo: str) -> ee.Dictionary:
    chaves = dic.keys()
    return dic.rename(chaves, chaves.map(lambda k: ee.String(prefixo).cat(k)))


def resumo_cena(img: ee.Image, par: str, tipo: str, area_comum: ee.Number) -> ee.Dictionary:
    regiao = img.geometry()
    vnir   = img.select(BANDAS_VNIR)
    dn     = vnir.updateMask(vnir.gt(0))
    pct    = ee.Reducer.percentile([2, 50, 98])

    chaves_ganho = img.propertyNames().filter(ee.Filter.stringContains('item', 'GAIN'))
    return (ee.Dictionary({
                'par': par, 'tipo': tipo, 'cena': img.id(),
                'data': ee.Date(img.get('system:time_start')).format('YYYY-MM-dd HH:mm'),
                'SOLAR_ELEVATION': img.get('SOLAR_ELEVATION'),
                'SOLAR_AZIMUTH': img.get('SOLAR_AZIMUTH'),
                'CLOUDCOVER': img.get('CLOUDCOVER'),
                'area_cena_km2': regiao.area(100).divide(1e6),
                'area_comum_par_km2': area_comum,
            })
            .combine(img.toDictionary(chaves_ganho))
            .combine(prefixar(reduzir(dn, pct, regiao), 'DN_'))
            .combine(prefixar(reduzir(dn_para_toa(img), pct, regiao), 'TOA_'))
            .combine(prefixar(reduzir(vnir.gte(DN_SAT).updateMask(vnir.gt(0)),
                                      ee.Reducer.mean(), regiao), 'SAT_')))


# ------- EXECUÇÃO -------
linhas, propriedades = [], {}
for par in PARES:
    cenas = {tipo: ee.Image(par[tipo]) for tipo in ('bom', 'ruim')}
    area_comum = (cenas['bom'].geometry()
                  .intersection(cenas['ruim'].geometry(), 100)
                  .area(100).divide(1e6))
    for tipo, img in cenas.items():
        log.info(f'Par {par["nome"]} — {tipo}: {par[tipo]}')
        linhas.append(resumo_cena(img, par['nome'], tipo, area_comum).getInfo())
        propriedades[par[tipo]] = img.toDictionary().getInfo()

# Ordem das colunas: identificação, metadados, ganho, DN, TOA, saturação, área
df = pd.DataFrame(linhas)
inicio = ['par', 'tipo', 'cena', 'data', 'SOLAR_ELEVATION', 'SOLAR_AZIMUTH', 'CLOUDCOVER']
ganho  = sorted(c for c in df.columns if 'GAIN' in c)
blocos = [sorted(c for c in df.columns if c.startswith(p)) for p in ('DN_', 'TOA_', 'SAT_')]
fim    = ['area_cena_km2', 'area_comum_par_km2']
df = df[inicio + ganho + sum(blocos, []) + fim]

saida_csv = PASTA_SAIDA / 'resumo_cenas_histograma.csv'
df.to_csv(saida_csv, index=False)
(PASTA_SAIDA / 'propriedades_cenas.json').write_text(
    json.dumps(propriedades, indent=2, ensure_ascii=False, default=str))

pd.set_option('display.width', 220, 'display.max_columns', None)
print(df.set_index(['par', 'tipo']).drop(columns='cena').T.to_string())
log.info(f'CSV  → {saida_csv}')
log.info(f'JSON → {PASTA_SAIDA / "propriedades_cenas.json"}')
