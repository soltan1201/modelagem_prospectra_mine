"""
Análise dos histogramas exportados por analise_histogramas_cenas.js (GEE).

Entrada (src/Dados/histograma/), para cada ano em PARES:
  DN_bruto_<ano>_bom.csv, DN_bruto_<ano>_ruim.csv → histogramas DN de B01, B02, B3N
  par_<ano>_bomxRuim.csv                          → histograma TOA da B3N, bom × ruim

Para cada banda:
  1. Estatísticas (média, desvio, P2/P50/P98) e fração saturada (DN = 255)
  2. Ajuste da cena ruim para a boa por três métodos:
       - linear média/desvio  (z-score)
       - linear por percentis (P2→P2, P98→P98)
       - casamento de histograma (CDF ruim → CDF boa), saturados excluídos
  3. Figura: histogramas, CDFs e curvas de transferência
  4. LUT do casamento de histograma pronta para o GEE (image.interpolate)

Saídas (por ano): figuras/dn_<ano>_bom_ruim_ajuste.png, figuras/toa_b3n_par_<ano>.png
                  e lut_casamento_histograma_<ano>.js
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PASTA   = Path(__file__).resolve().parents[1] / 'Dados' / 'histograma'
FIGURAS = PASTA / 'figuras'
BANDAS  = ['B01', 'B02', 'B3N']
DN_SAT  = 254          # DN ≥ 254 = saturado (VNIR 8 bits; o GEE espalha o 255 no bin 254–255)
N_LUT   = 41           # nós da LUT exportada para o GEE
PARES   = ['2003', '2001']

COR_BOM, COR_RUIM, COR_AJUSTE = '#2a78d6', '#e34948', '#2b2b2b'

FIGURAS.mkdir(parents=True, exist_ok=True)


# ------- LEITURA -------
def carregar(nome: str) -> pd.DataFrame:
    """CSV de histograma do GEE: 'Band Value' = 'lo - hi', contagens com milhar ','."""
    df = pd.read_csv(PASTA / nome, thousands=',')
    limites = df['Band Value'].str.split(' - ', expand=True).astype(float)
    df['centro'] = limites.mean(axis=1)
    df['lo'] = limites[0]
    return df.drop(columns='Band Value').fillna(0)


# ------- ESTATÍSTICA DE HISTOGRAMA -------
def cdf(contagem: np.ndarray) -> np.ndarray:
    c = np.cumsum(contagem)
    return c / c[-1]


def quantil(x: np.ndarray, contagem: np.ndarray, p: float) -> float:
    return float(np.interp(p, cdf(contagem), x))


def estatisticas(x: np.ndarray, contagem: np.ndarray) -> dict:
    n = contagem.sum()
    media = (x * contagem).sum() / n
    return {
        'n': int(round(n)),
        'media': media,
        'desvio': np.sqrt((contagem * (x - media) ** 2).sum() / n),
        'p2': quantil(x, contagem, 0.02),
        'p50': quantil(x, contagem, 0.50),
        'p98': quantil(x, contagem, 0.98),
    }


# ------- MÉTODOS DE AJUSTE (ruim → bom) -------
def linear_media_desvio(s_ruim: dict, s_bom: dict):
    ganho = s_bom['desvio'] / s_ruim['desvio']
    return ganho, s_bom['media'] - ganho * s_ruim['media']


def linear_percentis(s_ruim: dict, s_bom: dict):
    ganho = (s_bom['p98'] - s_bom['p2']) / (s_ruim['p98'] - s_ruim['p2'])
    return ganho, s_bom['p2'] - ganho * s_ruim['p2']


def casamento_histograma(x_ruim, c_ruim, x_bom, c_bom):
    """Transferência T(x): valor ruim x → valor bom com o mesmo percentil."""
    cdf_ruim, cdf_bom = cdf(c_ruim), cdf(c_bom)
    return lambda x: np.interp(np.interp(x, x_ruim, cdf_ruim), cdf_bom, x_bom)


def reamostrar(x, contagem, transf, bordas):
    """Histograma da cena ruim depois de aplicar a transferência."""
    return np.histogram(transf(x), bins=bordas, weights=contagem)[0]


# ------- ANÁLISE DE UM PAR -------
def analisar_par(ano: str) -> pd.DataFrame:
    bom  = carregar(f'DN_bruto_{ano}_bom.csv')
    ruim = carregar(f'DN_bruto_{ano}_ruim.csv')
    linhas, luts = [], {}

    fig, eixos = plt.subplots(3, len(BANDAS), figsize=(15, 11))
    for j, banda in enumerate(BANDAS):
        # Saturados (DN = 255) não têm valor real: ficam fora das estatísticas e do ajuste
        b = bom[bom['lo'] < DN_SAT]
        r = ruim[ruim['lo'] < DN_SAT]
        x_b, c_b = b['centro'].values, b[banda].values
        x_r, c_r = r['centro'].values, r[banda].values

        s_b, s_r = estatisticas(x_b, c_b), estatisticas(x_r, c_r)
        sat_b = bom.loc[bom['lo'] >= DN_SAT, banda].sum() / bom[banda].sum()
        sat_r = ruim.loc[ruim['lo'] >= DN_SAT, banda].sum() / ruim[banda].sum()

        g_md, o_md = linear_media_desvio(s_r, s_b)
        g_pc, o_pc = linear_percentis(s_r, s_b)
        transf = casamento_histograma(x_r, c_r, x_b, c_b)

        linhas.append({
            'banda': banda,
            'media_bom': s_b['media'], 'media_ruim': s_r['media'],
            'desvio_bom': s_b['desvio'], 'desvio_ruim': s_r['desvio'],
            'p50_bom': s_b['p50'], 'p50_ruim': s_r['p50'],
            'sat_bom_%': 100 * sat_b, 'sat_ruim_%': 100 * sat_r,
            'ganho_media_desvio': g_md, 'offset_media_desvio': o_md,
            'ganho_percentis': g_pc, 'offset_percentis': o_pc,
        })

        # LUT para o GEE: nós nos percentis da cena ruim
        nos = np.array([quantil(x_r, c_r, p) for p in np.linspace(0.001, 0.999, N_LUT)])
        nos = np.unique(np.round(nos, 1))
        luts[banda] = (nos, np.round(transf(nos), 2))

        # Linha 1: histogramas (bom, ruim, ruim ajustado)
        bordas = np.append(b['lo'].values, b['lo'].values[-1] + 1)
        ajust = reamostrar(x_r, c_r, transf, bordas)
        ax = eixos[0, j]
        ax.plot(x_b, c_b / c_b.sum(), color=COR_BOM, lw=2, label='bom')
        ax.plot(x_r, c_r / c_r.sum(), color=COR_RUIM, lw=2, label='ruim')
        ax.plot(x_b, ajust / ajust.sum(), color=COR_AJUSTE, lw=1.5, ls='--', label='ruim ajustado')
        ax.set_title(f'{banda} — histograma (DN, sem saturados)')
        ax.set_xlabel('DN'); ax.set_ylabel('Fração de pixels')
        ax.text(0.98, 0.95, f'saturado: bom {100 * sat_b:.1f}% | ruim {100 * sat_r:.1f}%',
                transform=ax.transAxes, ha='right', va='top', fontsize=8, color='#555')

        # Linha 2: CDFs
        ax = eixos[1, j]
        ax.plot(x_b, cdf(c_b), color=COR_BOM, lw=2, label='bom')
        ax.plot(x_r, cdf(c_r), color=COR_RUIM, lw=2, label='ruim')
        ax.set_title(f'{banda} — distribuição acumulada')
        ax.set_xlabel('DN'); ax.set_ylabel('CDF')

        # Linha 3: curvas de transferência ruim → bom
        ax = eixos[2, j]
        xs = np.linspace(x_r.min(), x_r.max(), 200)
        ax.plot(xs, xs, color='#bbb', lw=1, ls=':', label='identidade')
        ax.plot(xs, transf(xs), color=COR_AJUSTE, lw=2, label='casamento de histograma')
        ax.plot(xs, g_md * xs + o_md, color=COR_BOM, lw=1.5, label=f'média/desvio  {g_md:.2f}x{o_md:+.1f}')
        ax.plot(xs, g_pc * xs + o_pc, color=COR_RUIM, lw=1.5, label=f'percentis  {g_pc:.2f}x{o_pc:+.1f}')
        ax.set_title(f'{banda} — transferência ruim → bom')
        ax.set_xlabel('DN ruim'); ax.set_ylabel('DN ajustado')
        ax.legend(fontsize=7, loc='upper left')

    for ax in eixos[:2, 0]:
        ax.legend(fontsize=8)
    for ax in eixos.flat:
        ax.grid(alpha=0.3)
    fig.suptitle(f'Par {ano} — cena ruim × cena boa, DN bruto (B01, B02, B3N)', fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGURAS / f'dn_{ano}_bom_ruim_ajuste.png', dpi=110, bbox_inches='tight')
    plt.close(fig)

    # Tabela-resumo
    tabela = pd.DataFrame(linhas).set_index('banda').round(3)
    print(f'\n══════════ PAR {ano} ══════════')
    print(tabela.T)

    # Razão NIR/vermelho (o "vermelho" da falsa-cor) nas medianas
    for nome, col in [('bom', 'p50_bom'), ('ruim', 'p50_ruim')]:
        print(f'B3N/B02 (mediana DN) {nome}: {tabela.loc["B3N", col] / tabela.loc["B02", col]:.2f}')

    # TOA B3N
    par = carregar(f'par_{ano}_bomxRuim.csv')
    for nome in ('bom', 'ruim'):
        s = estatisticas(par['centro'].values, par[nome].values)
        print(f'TOA B3N {nome:4s}: média={s["media"]:.3f}  P50={s["p50"]:.3f}  desvio={s["desvio"]:.3f}')

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(par['centro'], par['bom'] / par['bom'].sum(), color=COR_BOM, lw=2, label='bom')
    ax.plot(par['centro'], par['ruim'] / par['ruim'].sum(), color=COR_RUIM, lw=2, label='ruim')
    ax.set_title(f'Par {ano} — B3N em reflectância TOA'); ax.set_xlabel('Reflectância TOA')
    ax.set_ylabel('Fração de pixels'); ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURAS / f'toa_b3n_par_{ano}.png', dpi=110, bbox_inches='tight')
    plt.close(fig)

    # LUT para o GEE
    js = [f'// LUT do casamento de histograma (cena ruim → cena boa), par {ano}, DN bruto.',
          '// Gerado por analise_histogramas_csv.py. Uso:',
          '//   var ajustada = aplicarLUT(ee.Image(id_ruim));',
          'var LUT = {']
    for banda, (xs, ys) in luts.items():
        js.append(f"  '{banda}': {{x: {xs.tolist()}, y: {ys.tolist()}}},")
    js += ['};',
           'function aplicarLUT(img) {',
           '  return ee.Image.cat(Object.keys(LUT).map(function(b) {',
           '    var banda = img.select(b);',
           '    // saturados (DN = 255) ficam mascarados: não há valor real para ajustar',
           f'    return banda.updateMask(banda.lt({DN_SAT}))',
           "                .interpolate(LUT[b].x, LUT[b].y, 'clamp').rename(b);",
           '  }));',
           '}']
    (PASTA / f'lut_casamento_histograma_{ano}.js').write_text('\n'.join(js) + '\n')
    return tabela


pd.set_option('display.width', 200)
for ano in PARES:
    analisar_par(ano)
print(f'\nFiguras em {FIGURAS}')
