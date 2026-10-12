"""
Piloto — correção da saturação de B01/B02 do ASTER com Landsat (32 dias).

As bandas B01/B02 do ASTER são gravadas em ganho alto (HGH) e saturam em
DN ≥ 254 sobre alvos claros: o valor fica cortado no teto. A ideia, adaptada de
Gonzalez et al. (2019, mosaico global ASTER de 15 m, que recupera B01/B02 a
partir da B3N com polinômios usando MODIS), tem três etapas:

  1. Relação nos pixels NÃO saturados: ajustar ρ_ASTER(B01/B02) a partir do
     Landsat (green, red, nir — mediana dos compostos de 32 dias ao redor da
     data da cena) e escolher a melhor forma (linear, polinômio de 2ª ordem,
     com/sem a B3N do ASTER).
  2. Correção dos saturados: aplicar a relação aos pixels com DN ≥ 254, com a
     restrição ρ_corrigido ≥ ρ_cortado (o real só pode ser maior).
  3. Ajuste entre bandas do ASTER: relação B01/B02 ← B3N (+ SWIR + TIR) nos não
     saturados; serve para refinar a estimativa do Landsat (média ponderada pelo
     erro de validação) e de plano B onde o Landsat não tem dado.

Como medir sem verdade nos saturados — PSEUDO-SATURAÇÃO: treina-se só com
pixels de DN < LIMIAR_PSEUDO e testa-se na faixa LIMIAR_PSEUDO ≤ DN < 254, onde
o valor real é conhecido. A comparação é com a linha de base "cortado" (o valor
que o pixel teria se saturasse no limiar).

Tudo em 30 m: o ASTER (15 m) é agregado em blocos 2×2 (média da reflectância,
máximo do DN para marcar saturação) na grade de 30 m da própria cena.
Modelos ajustados por cena (cada cena tem sua atmosfera; o Landsat é
reflectância de superfície e o ASTER é TOA).

Saídas (src/Dados/saturacao/):
  amostras_<cena>.csv, metricas_pseudo_saturacao.csv, correcao_saturados.csv,
  figuras/pseudo_<cena>_<banda>.png, figuras/saturados_<cena>.png

Uso:  python saturacao_landsat_piloto.py
"""

import logging
import math
import sys
from pathlib import Path

import ee
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

PASTA = Path(__file__).resolve().parent
sys.path.append(str(PASTA.parent))
from configure_account_projects_ee import get_current_account   # noqa: E402

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
# Cenas de teste: as zero-cloud com mais saturação em B01/B02 (2006: tem SWIR)
CENAS = ['20190304130639', '20200111130040', '20100703125452', '20060708125414']

JANELA_DIAS   = 32        # compostos Landsat com início a ± JANELA_DIAS da data da cena
ESCALA        = 30        # m — grade comum ASTER/Landsat
DN_SATURADO   = 254
LIMIAR_PSEUDO = 200       # DN: treino abaixo, teste entre o limiar e 253
N_POR_CLASSE  = {0: 4000, 1: 3000, 2: 3000}   # 0: DN < limiar | 1: pseudo | 2: saturado
SEMENTE       = 42

ESUN = {'B01': 1848, 'B02': 1549, 'B3N': 1114,
        'B04': 225.4, 'B05': 86.63, 'B06': 81.85, 'B07': 74.85, 'B08': 66.49, 'B09': 59.85}
SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']
TIR  = ['B10', 'B11', 'B12', 'B13', 'B14']
LS   = ['ls_blue', 'ls_green', 'ls_red', 'ls_nir', 'ls_swir1', 'ls_swir2']
ALVOS = ['B01', 'B02']

SAIDA   = PASTA.parent / 'Dados' / 'saturacao'
FIGURAS = SAIDA / 'figuras'


# ═══════════════════════════════════════════════════════════════════════════
# GEE — AMOSTRAS POR CENA
# ═══════════════════════════════════════════════════════════════════════════
def toa(img: ee.Image, bandas: list[str]) -> ee.Image:
    """DN → radiância (ganho) → reflectância TOA 0–1, sem máscara de saturação."""
    cosz = ee.Number(90).subtract(img.get('SOLAR_ELEVATION')).multiply(math.pi / 180).cos()
    doy = ee.Number(ee.Date(img.get('system:time_start')).getRelative('day', 'year'))
    d = ee.Number(1).subtract(ee.Number(0.01672).multiply(doy.multiply(2 * math.pi / 365).cos()))
    saida = []
    for b in bandas:
        dn = img.select(b)
        saida.append(dn.subtract(1).multiply(ee.Number(img.get(f'GAIN_COEFFICIENT_{b}')))
                       .multiply(math.pi).multiply(d.pow(2)).divide(cosz.multiply(ESUN[b]))
                       .updateMask(dn.gt(0)).toFloat().rename(b))
    return ee.Image.cat(saida)


def amostrar_cena(cena: str) -> pd.DataFrame:
    img = ee.Image(f'ASTER/AST_L1T_003/{cena}')
    info = ee.Dictionary({'tem_swir': img.propertyNames().contains('GAIN_COEFFICIENT_B04'),
                          'data': ee.Date(img.get('system:time_start')).format('YYYY-MM-dd')}).getInfo()
    opticas = ['B01', 'B02', 'B3N'] + (SWIR if info['tem_swir'] else [])
    pj30 = img.select('B01').projection().atScale(ESCALA)

    refl = toa(img, opticas)
    dn = img.select(['B01', 'B02', 'B3N']).rename(['DN_B01', 'DN_B02', 'DN_B3N'])
    # agregação 2×2 na grade de 30 m da cena: média da reflectância e do DN, máximo do DN
    media = (refl.addBands(dn.toFloat()).addBands(img.select(TIR).toFloat())
                 .reduceResolution(ee.Reducer.mean(), maxPixels=64).reproject(pj30))
    dn_max = (dn.reduceResolution(ee.Reducer.max(), maxPixels=64).reproject(pj30)
                .rename(['DNmax_B01', 'DNmax_B02', 'DNmax_B3N']))

    data = ee.Date(img.get('system:time_start'))
    landsat = (ee.ImageCollection('LANDSAT/COMPOSITES/C02/T1_L2_32DAY')
                 .filterDate(data.advance(-JANELA_DIAS, 'day'), data.advance(JANELA_DIAS, 'day'))
                 .select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'], LS)
                 .median())

    # classe pela banda mais saturada de B01/B02 (máximo DN do bloco 2×2)
    dmax = dn_max.select(['DNmax_B01', 'DNmax_B02']).reduce(ee.Reducer.max())
    classe = (ee.Image(0).where(dmax.gte(LIMIAR_PSEUDO), 1).where(dmax.gte(DN_SATURADO), 2)
                .rename('classe').toInt())
    pilha = (media.addBands(dn_max).addBands(landsat).addBands(classe)
                  .updateMask(media.select('B3N').mask())
                  .updateMask(landsat.select('ls_red').mask()))

    fc = pilha.stratifiedSample(
        numPoints=0, classBand='classe', region=img.geometry(), scale=ESCALA, projection=pj30,
        classValues=list(N_POR_CLASSE), classPoints=list(N_POR_CLASSE.values()),
        seed=SEMENTE, tileScale=4, geometries=False)
    df = ee.data.computeFeatures({'expression': fc, 'fileFormat': 'PANDAS_DATAFRAME'})
    df = df.drop(columns=[c for c in df.columns if c == 'geo'])
    df['cena'], df['data'], df['tem_swir'] = cena, info['data'], info['tem_swir']
    return df


# ═══════════════════════════════════════════════════════════════════════════
# MODELOS
# ═══════════════════════════════════════════════════════════════════════════
def candidatos(tem_swir: bool) -> dict[str, tuple[list[str], int]]:
    """nome → (variáveis, grau do polinômio)."""
    ls3 = ['ls_green', 'ls_red', 'ls_nir']
    intra = ['B3N'] + (SWIR if tem_swir else []) + ['B10', 'B13']
    return {
        # etapa 1 — Landsat
        'LS_lin':        (ls3, 1),
        'LS_poly2':      (ls3, 2),
        'LS+B3N_lin':    (ls3 + ['B3N'], 1),
        'LS+B3N_poly2':  (ls3 + ['B3N'], 2),
        # etapa 3 — só ASTER (entre bandas)
        'B3N_poly2':     (['B3N'], 2),
        'ASTER_lin':     (intra, 1),
        'ASTER_poly2':   (intra[:3], 2),     # B3N + 2 primeiras (evita explosão de termos)
    }


def ajustar(df_treino: pd.DataFrame, vars_: list[str], grau: int, alvo: str):
    modelo = make_pipeline(PolynomialFeatures(grau, include_bias=False), LinearRegression())
    return modelo.fit(df_treino[vars_].values, df_treino[alvo].values)


def metricas(obs, prev) -> dict:
    erro = prev - obs
    return {'rmse': float(np.sqrt(np.mean(erro ** 2))), 'vies': float(np.mean(erro)),
            'r2': float(1 - np.sum(erro ** 2) / np.sum((obs - obs.mean()) ** 2))}


def avaliar_cena(df: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    """Pseudo-saturação por banda + correção dos saturados reais."""
    cena = df['cena'].iloc[0]
    tem_swir = bool(df['tem_swir'].iloc[0])
    cands = candidatos(tem_swir)
    linhas_met, linhas_sat = [], []

    for alvo in ALVOS:
        dn_max, dn_med = df[f'DNmax_{alvo}'], df[f'DN_{alvo}']
        validos = df.dropna(subset=[alvo] + sum((v for v, _ in cands.values()), []))
        treino = validos[dn_max.loc[validos.index] < LIMIAR_PSEUDO]
        teste = validos[(dn_max.loc[validos.index] >= LIMIAR_PSEUDO) &
                        (dn_max.loc[validos.index] < DN_SATURADO)]
        saturados = validos[dn_max.loc[validos.index] >= DN_SATURADO]
        if len(treino) < 200 or len(teste) < 50:
            log.warning(f'{cena} {alvo}: poucas amostras (treino {len(treino)}, teste {len(teste)})')
            continue

        # metade do teste escolhe os pesos da combinação, a outra metade avalia
        rng = np.random.default_rng(SEMENTE)
        metade = rng.random(len(teste)) < 0.5
        teste_a, teste_b = teste[metade], teste[~metade]

        # linha de base: o pixel "saturaria" no limiar → valor cortado em DN = limiar
        fator = (np.minimum(teste_b[f'DN_{alvo}'], LIMIAR_PSEUDO - 1) - 1) / (teste_b[f'DN_{alvo}'] - 1)
        base = metricas(teste_b[alvo].values, (teste_b[alvo] * fator).values)
        linhas_met.append({'cena': cena, 'banda': alvo, 'modelo': 'cortado (linha de base)',
                           'n_treino': len(treino), 'n_teste': len(teste_b), **base})

        prev_a, prev_b, modelos = {}, {}, {}
        for nome, (vars_, grau) in cands.items():
            m = ajustar(treino, vars_, grau, alvo)
            modelos[nome] = (m, vars_)
            prev_a[nome] = m.predict(teste_a[vars_].values)
            prev_b[nome] = m.predict(teste_b[vars_].values)
            linhas_met.append({'cena': cena, 'banda': alvo, 'modelo': nome,
                               'n_treino': len(treino), 'n_teste': len(teste_b),
                               **metricas(teste_b[alvo].values, prev_b[nome])})

        # melhor Landsat e melhor intra-ASTER escolhidos na metade A; combinação 1/MSE
        mse_a = {n: np.mean((p - teste_a[alvo].values) ** 2) for n, p in prev_a.items()}
        melhor_ls = min((n for n in cands if n.startswith('LS')), key=mse_a.get)
        melhor_in = min((n for n in cands if not n.startswith('LS')), key=mse_a.get)
        w_ls = (1 / mse_a[melhor_ls]) / (1 / mse_a[melhor_ls] + 1 / mse_a[melhor_in])
        comb_b = w_ls * prev_b[melhor_ls] + (1 - w_ls) * prev_b[melhor_in]
        linhas_met.append({'cena': cena, 'banda': alvo,
                           'modelo': f'combinado ({melhor_ls} {w_ls:.2f} + {melhor_in} {1 - w_ls:.2f})',
                           'n_treino': len(treino), 'n_teste': len(teste_b),
                           **metricas(teste_b[alvo].values, comb_b)})
        figura_pseudo(cena, alvo, teste_b, teste_b[alvo] * fator, prev_b[melhor_ls], comb_b, melhor_ls)

        # etapa 2 — saturados reais: re-treina com TODOS os não saturados e aplica
        nao_sat = validos[dn_max.loc[validos.index] < DN_SATURADO]
        if len(saturados):
            m_ls = ajustar(nao_sat, modelos[melhor_ls][1], cands[melhor_ls][1], alvo)
            m_in = ajustar(nao_sat, modelos[melhor_in][1], cands[melhor_in][1], alvo)
            p = (w_ls * m_ls.predict(saturados[modelos[melhor_ls][1]].values) +
                 (1 - w_ls) * m_in.predict(saturados[modelos[melhor_in][1]].values))
            corrigido = np.maximum(p, saturados[alvo].values)      # real ≥ cortado
            linhas_sat.append({
                'cena': cena, 'banda': alvo, 'n_saturados': len(saturados),
                'cortado_mediana': float(np.median(saturados[alvo])),
                'corrigido_mediana': float(np.median(corrigido)),
                'aumento_mediano_%': float(np.median(corrigido / saturados[alvo] - 1) * 100),
                'restricao_ativa_%': float(np.mean(p < saturados[alvo].values) * 100),
                'modelo': f'{melhor_ls} {w_ls:.2f} + {melhor_in} {1 - w_ls:.2f}'})
            figura_saturados(cena, alvo, nao_sat, saturados, corrigido)
    return linhas_met, linhas_sat


# ═══════════════════════════════════════════════════════════════════════════
# FIGURAS
# ═══════════════════════════════════════════════════════════════════════════
def figura_pseudo(cena, alvo, teste, cortado, prev_ls, comb, nome_ls):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    obs = teste[alvo].values
    lim = [min(obs.min(), 0), max(obs.max(), comb.max()) * 1.05]
    ax.plot(lim, lim, color='#888888', lw=1)
    ax.scatter(obs, cortado, s=4, alpha=0.4, color='#9e2a2a', label='cortado no limiar')
    ax.scatter(obs, prev_ls, s=4, alpha=0.4, color='#f0a35e', label=nome_ls)
    ax.scatter(obs, comb, s=4, alpha=0.4, color='#2a78d6', label='combinado')
    ax.set(xlim=lim, ylim=lim, xlabel=f'{alvo} observado (TOA)', ylabel=f'{alvo} estimado',
           title=f'Pseudo-saturação {cena} — {alvo}\n(treino DN < {LIMIAR_PSEUDO}, teste {LIMIAR_PSEUDO}–253)')
    ax.legend(fontsize=8, markerscale=3)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURAS / f'pseudo_{cena}_{alvo}.png', dpi=110)
    plt.close(fig)


def figura_saturados(cena, alvo, nao_sat, saturados, corrigido):
    fig, ax = plt.subplots(figsize=(6.5, 4))
    faixa = np.linspace(0, max(corrigido.max(), nao_sat[alvo].max()) * 1.05, 80)
    ax.hist(nao_sat[alvo], bins=faixa, color='#888888', alpha=0.5, label='não saturados')
    ax.hist(saturados[alvo], bins=faixa, color='#9e2a2a', alpha=0.6, label='saturados (cortado)')
    ax.hist(corrigido, bins=faixa, color='#2a78d6', alpha=0.6, label='saturados (corrigido)')
    ax.set(xlabel=f'{alvo} reflectância TOA', ylabel='pixels (amostra)', title=f'{cena} — {alvo}')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURAS / f'saturados_{cena}_{alvo}.png', dpi=110)
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════
def main():
    ee.Initialize(project=get_current_account())
    FIGURAS.mkdir(parents=True, exist_ok=True)
    todas_met, todas_sat = [], []
    for cena in CENAS:
        arq = SAIDA / f'amostras_{cena}.csv'
        if arq.exists():
            df = pd.read_csv(arq)
            log.info(f'{cena}: amostras lidas de {arq.name} ({len(df)})')
        else:
            log.info(f'{cena}: amostrando no GEE…')
            df = amostrar_cena(cena)
            df.to_csv(arq, index=False)
            log.info(f'{cena}: {len(df)} amostras (classes {df["classe"].value_counts().to_dict()})')
        met, sat = avaliar_cena(df)
        todas_met += met
        todas_sat += sat

    met = pd.DataFrame(todas_met)
    sat = pd.DataFrame(todas_sat)
    met.to_csv(SAIDA / 'metricas_pseudo_saturacao.csv', index=False)
    sat.to_csv(SAIDA / 'correcao_saturados.csv', index=False)

    pd.set_option('display.width', 200, 'display.max_columns', None, 'display.max_colwidth', 60)
    print('\n══ Pseudo-saturação (teste DN', LIMIAR_PSEUDO, '–253; reflectância TOA 0–1) ══')
    print(met.round(4).to_string(index=False))
    print('\n══ Saturados reais: valor cortado × corrigido ══')
    print(sat.round(4).to_string(index=False))
    log.info(f'Saídas em {SAIDA}')


if __name__ == '__main__':
    main()
