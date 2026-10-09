"""
Diagnóstico das 14 bandas ASTER (B01–B14) nas cenas boas × ruins.

Para cada cena e banda compara três caminhos do mesmo DN bruto:

  1. FÍSICO (correto)
       radiância L = (DN − 1) × GAIN_COEFFICIENT_Bxx
       VNIR/SWIR → reflectância TOA = π·L·d² / (ESUN·cos θz)
       TIR       → temperatura de brilho T = K2 / ln(K1/L + 1)   [K]
  2. PIPELINE ATUAL (select_black_list_save_mosaicSem.js, com os dois bugs)
       VNIR/SWIR → π·DN·d² / (ESUN·cos θz) × 10000   (DN tratado como radiância)
       TIR       → DN × 10000                         (escala[b] sempre caía em 10000)
       → toInt16 corta em 32767
  3. PIPELINE CORRIGIDO
       VNIR/SWIR → TOA física × 10000
       TIR       → só o fix do INT16 (DN × 10) e a proposta (temperatura × 10)

Métricas por banda: percentis P2/P50/P98 de DN, radiância, valor físico e
dos valores do pipeline; fração de DN saturado; fração cortada em 32767.

Uso:
  python analise_todas_bandas.py
  python analise_todas_bandas.py --project outro-projeto-gee

Saídas (src/Dados/histograma/):
  analise_todas_bandas.csv      — uma linha por cena × banda
  figuras/analise_todas_bandas.png
"""

import argparse
import logging
import math
from pathlib import Path

import ee
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# ------- PARÂMETROS -------
PARES = [
    {'nome': '2003', 'bom': 'ASTER/AST_L1T_003/20030914131811', 'ruim': 'ASTER/AST_L1T_003/20031011130046'},
    {'nome': '2001', 'bom': 'ASTER/AST_L1T_003/20010816131957', 'ruim': 'ASTER/AST_L1T_003/20010928130033'},
]

ESCALA    = 60        # m
INT16_MAX = 32767

# Irradiância solar exoatmosférica — mesmos valores do radianceToTOA do pipeline
ESUN = {'B01': 1848, 'B02': 1549, 'B3N': 1114,
        'B04': 225.4, 'B05': 86.63, 'B06': 81.85, 'B07': 74.85, 'B08': 66.49, 'B09': 59.85}

# Constantes de Planck por banda TIR (K1 em W/m²/sr/µm, K2 em K),
# a partir dos comprimentos de onda centrais 8.291, 8.634, 9.075, 10.657, 11.318 µm
K_TIR = {'B10': (3040.136402, 1735.337945), 'B11': (2482.375199, 1666.398761),
         'B12': (1935.060183, 1585.420044), 'B13': (866.468575, 1350.069147),
         'B14': (641.326517, 1271.221673)}

OPTICAS = list(ESUN)
TIR     = list(K_TIR)
BANDAS  = OPTICAS + TIR
DN_MAX  = {**dict.fromkeys(OPTICAS, 255), **dict.fromkeys(TIR, 4095)}   # 8 bits / 12 bits

PASTA   = Path(__file__).resolve().parents[1] / 'Dados' / 'histograma'
FIGURAS = PASTA / 'figuras'

parser = argparse.ArgumentParser(description='Diagnóstico das 14 bandas ASTER.')
parser.add_argument('--project', default='mapbiomas-caatinga-cloud02',
                    help='projeto GEE usado no ee.Initialize')
args = parser.parse_args()

try:
    ee.Initialize(project='ee-solkancengine17')
    log.info(f'Earth Engine inicializado (projeto {args.project}).')
except Exception as e:
    log.error(f'Erro de inicialização: {e}')
    raise

FIGURAS.mkdir(parents=True, exist_ok=True)


# ------- CAMADAS POR BANDA -------
def geometria_solar(img: ee.Image):
    elev = ee.Number(img.get('SOLAR_ELEVATION'))
    cosz = ee.Number(90).subtract(elev).multiply(math.pi / 180).cos()
    data = ee.Date(img.get('system:time_start'))
    # d físico (UA) e o d exatamente como o pipeline calcula hoje
    doy_1 = data.getRelative('day', 'year').add(1)
    d_fis = ee.Number(1).subtract(
        ee.Number(0.01672).multiply(doy_1.subtract(4).multiply(0.9856 * math.pi / 180).cos()))
    d_pipe = ee.Number(1).subtract(
        ee.Number(0.01672).multiply(data.getRelative('day', 'year').multiply(2 * math.pi / 365).cos()))
    return cosz, d_fis, d_pipe


def camadas(img: ee.Image) -> tuple[ee.Image, ee.Image]:
    """Imagem de valores contínuos e imagem de flags (0/1), uma banda por métrica."""
    cosz, d_fis, d_pipe = geometria_solar(img)
    continuas, flags = [], []
    for b in BANDAS:
        dn = img.select(b).toFloat()
        dn = dn.updateMask(dn.gt(0))
        rad = dn.subtract(1).multiply(ee.Number(img.get(f'GAIN_COEFFICIENT_{b}')))

        if b in OPTICAS:
            fisico = rad.multiply(math.pi).multiply(d_fis.pow(2)).divide(cosz.multiply(ESUN[b]))
            atual  = (dn.multiply(math.pi).multiply(d_pipe.pow(2))
                        .divide(cosz.multiply(ESUN[b])).multiply(10000))
            corrigido = fisico.multiply(10000)
            so_int16  = atual                      # o fix do INT16 não muda VNIR/SWIR
        else:
            k1, k2 = K_TIR[b]
            fisico    = ee.Image(k2).divide(ee.Image(k1).divide(rad).add(1).log())   # K
            atual     = dn.multiply(10000)
            so_int16  = dn.multiply(10)
            corrigido = fisico.multiply(10)

        continuas += [dn.rename(f'{b}__dn'), rad.rename(f'{b}__rad'),
                      fisico.rename(f'{b}__fisico'), atual.rename(f'{b}__atual'),
                      so_int16.rename(f'{b}__so_int16'), corrigido.rename(f'{b}__corrigido')]
        flags += [dn.gte(DN_MAX[b] - (1 if b in OPTICAS else 0)).rename(f'{b}__sat_dn'),
                  atual.gte(INT16_MAX).rename(f'{b}__corte_atual'),
                  so_int16.gte(INT16_MAX).rename(f'{b}__corte_so_int16'),
                  corrigido.gte(INT16_MAX).rename(f'{b}__corte_corrigido')]
    return ee.Image.cat(continuas), ee.Image.cat(flags)


def reduzir(img: ee.Image, reducer: ee.Reducer, regiao: ee.Geometry) -> ee.Dictionary:
    return img.reduceRegion(reducer=reducer, geometry=regiao, scale=ESCALA,
                            maxPixels=1e10, bestEffort=True, tileScale=4)


# ------- EXECUÇÃO -------
linhas = []
for par in PARES:
    for tipo in ('bom', 'ruim'):
        cena_id = par[tipo]
        log.info(f'Par {par["nome"]} — {tipo}: {cena_id}')
        img = ee.Image(cena_id)
        cont, flag = camadas(img)
        regiao = img.geometry()
        pct = reduzir(cont, ee.Reducer.percentile([2, 50, 98]), regiao).getInfo()
        frac = reduzir(flag, ee.Reducer.mean(), regiao).getInfo()
        ganhos = img.toDictionary([f'GAIN_SETTING_{b}' for b in BANDAS]).getInfo()

        for b in BANDAS:
            linha = {'par': par['nome'], 'tipo': tipo, 'cena': cena_id.split('/')[-1],
                     'banda': b, 'grupo': 'VNIR' if b in OPTICAS[:3] else 'SWIR' if b in OPTICAS else 'TIR',
                     'ganho': ganhos.get(f'GAIN_SETTING_{b}')}
            for metrica in ('dn', 'rad', 'fisico', 'atual', 'so_int16', 'corrigido'):
                for p in (2, 50, 98):
                    linha[f'{metrica}_p{p}'] = pct.get(f'{b}__{metrica}_p{p}')
            for f in ('sat_dn', 'corte_atual', 'corte_so_int16', 'corte_corrigido'):
                v = frac.get(f'{b}__{f}')
                linha[f'{f}_%'] = None if v is None else 100 * v
            linhas.append(linha)

df = pd.DataFrame(linhas)
df.to_csv(PASTA / 'analise_todas_bandas.csv', index=False)

# ------- RESUMO NO TERMINAL -------
pd.set_option('display.width', 220, 'display.max_columns', None, 'display.max_rows', 200)
resumo = (df.groupby('banda', sort=False)
            [['dn_p50', 'fisico_p50', 'atual_p50', 'sat_dn_%', 'corte_atual_%',
              'corte_so_int16_%', 'corte_corrigido_%']]
            .mean().round(3))
print('\nMédia das 4 cenas por banda (físico: reflectância TOA em VNIR/SWIR, temperatura K em TIR)')
print(resumo.to_string())

# ------- FIGURA -------
COR = {('2003', 'bom'): '#2a78d6', ('2003', 'ruim'): '#e34948',
       ('2001', 'bom'): '#104281', ('2001', 'ruim'): '#9e2a2a'}
fig, eixos = plt.subplots(1, 3, figsize=(18, 5), gridspec_kw={'width_ratios': [9, 5, 14]})

for (p, t), g in df.groupby(['par', 'tipo'], sort=False):
    g = g.set_index('banda')
    estilo = '-' if t == 'bom' else '--'
    eixos[0].plot(OPTICAS, g.loc[OPTICAS, 'fisico_p50'], estilo, marker='o', ms=6,
                  color=COR[(p, t)], lw=2, label=f'{p} {t}')
    eixos[1].plot(TIR, g.loc[TIR, 'fisico_p50'], estilo, marker='o', ms=6,
                  color=COR[(p, t)], lw=2, label=f'{p} {t}')
eixos[0].set_title('Reflectância TOA física (mediana)')
eixos[0].set_ylabel('Reflectância TOA')
eixos[1].set_title('Temperatura de brilho TIR (mediana)')
eixos[1].set_ylabel('K')
eixos[0].legend(fontsize=8)

x = np.arange(len(BANDAS))
largura = 0.27
for i, (col, rotulo, cor) in enumerate([
        ('corte_atual_%', 'pipeline atual', '#e34948'),
        ('corte_so_int16_%', 'só fix do INT16', '#f0a35e'),
        ('corte_corrigido_%', 'INT16 + ganho corrigidos', '#2a78d6')]):
    valores = df.groupby('banda', sort=False)[col].mean().reindex(BANDAS)
    eixos[2].bar(x + (i - 1) * largura, valores, largura, label=rotulo, color=cor)
eixos[2].set_xticks(x, BANDAS)
eixos[2].set_ylabel('% de pixels = 32767 (média das 4 cenas)')
eixos[2].set_title('Pixels cortados no INT16 por banda')
eixos[2].set_ylim(0, 105)
eixos[2].legend(fontsize=8)

for ax in eixos:
    ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIGURAS / 'analise_todas_bandas.png', dpi=110, bbox_inches='tight')

log.info(f'CSV    → {PASTA / "analise_todas_bandas.csv"}')
log.info(f'Figura → {FIGURAS / "analise_todas_bandas.png"}')
