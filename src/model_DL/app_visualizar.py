"""
Visualizador dos mapas da CNN (Streamlit + Plotly): o mapa de p de uma rodada à
esquerda e, à direita, com zoom e arrasto sincronizados, um RGB artificial, um
canal/índice, outra rodada ou a diferença entre rodadas. Ocorrências de fosfato
por cima, com p da CNN × p_RF × p_GTB.

Roda no servidor (lê os arquivos locais, sem GEE) e abre no navegador do seu
computador.

Entradas (as mesmas variáveis dos outros scripts)
  RUNS_DIR     runs_cnn/<rodada>/mapas_passo<k>/p_<LxCy>.tif   (UInt16, p × 10000, 0 = sem dado)
  TILES_DIR    tiles_cnn/canais_<A|B1>/<LxCy>/<LxCy>_r<i>c<j>.tif   (float32, bandas com nome)
  PATCHES_DIR  patches_cnn/indice_patches.csv   (ocorrências: lon, lat, papel, p_rf, p_gtb)

Instalar (uma vez, no ambiente da CNN)
  conda activate py_DL
  pip install streamlit plotly

Rodar no servidor (dentro de um tmux)
  cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
  streamlit run app_visualizar.py --server.port 8501 --server.headless true

Abrir no seu computador
  ssh -p 2222 -N -L 8501:localhost:8501 superusuario@192.168.25.250
  → http://localhost:8501
  (ou direto http://192.168.25.250:8501 se a porta estiver liberada na rede)
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import rasterio
import streamlit as st
from PIL import Image
from plotly.subplots import make_subplots
from rasterio.enums import Resampling
from rasterio.merge import merge
from rasterio.transform import Affine, rowcol

RUNS_DIR = Path(os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn')).expanduser()
TILES_DIR = Path(os.environ.get('TILES_DIR', '/srv/almacen/db_images/tiles_cnn')).expanduser()
PATCHES_DIR = Path(os.environ.get('PATCHES_DIR', '/srv/almacen/db_images/patches_cnn')).expanduser()

AREA_TODA = 'Área toda'

# composições prontas por modelo de canais (nomes = descrições das bandas dos pedaços)
COMPOSICOES = {
    'B1': {
        'Falsa-cor VNIR (3N, 2, 1)': ('log_B3N', 'log_B02', 'log_B01'),
        'SWIR (4, 6, 8)': ('log_B04', 'log_B06', 'log_B08'),
        'SWIR minerais (5, 7, 9)': ('log_B05', 'log_B07', 'log_B09'),
        'Emissividade TIR (14, 12, 10)': ('emis_B14', 'emis_B12', 'emis_B10'),
    },
    'A': {
        'Carbonato · fosforito · ferro': ('carbonate_8_6', 'phosphorite_rbd', 'ferric_iron_2_1'),
        'Al-OH · dolomita · quartzo': ('sericite_illite_smectite', 'dolomite', 'quartz_14_12'),
        'Sultan R · Sultan G · SAVI': ('sultan_r', 'sultan_g', 'SAVI'),
        'Gossan · caulinita · carbonato TIR': ('gossan_4_2', 'kaolinite_7_5', 'carbonate_13_14'),
    },
}
PERSONALIZADO = 'Personalizado (escolher 3 canais)'


# ═══════════════════════════════════════════════════════════════════════════
# LEITURA (com cache: trocar de painel não relê o disco)
# ═══════════════════════════════════════════════════════════════════════════
def listar_rodadas():
    """{rodada: [passos com mapas]} — só rodadas com algum p_*.tif."""
    saida = {}
    if not RUNS_DIR.is_dir():
        return saida
    for run in sorted(p for p in RUNS_DIR.iterdir() if p.is_dir()):
        passos = sorted(int(d.name.removeprefix('mapas_passo')) for d in run.glob('mapas_passo*')
                        if d.is_dir() and d.name.removeprefix('mapas_passo').isdigit() and any(d.glob('p_*.tif')))
        if passos:
            saida[run.name] = passos
    return saida


def celula_de(arq):
    """p_<LxCy>[_<rodada>].tif → LxCy."""
    return re.match(r'p_(L\d+C\d+)', arq.name).group(1)


def arquivos_p(rodada, passo):
    return sorted((RUNS_DIR / rodada / f'mapas_passo{passo}').glob('p_L*C*.tif'))


def celulas(rodada, passo):
    return [celula_de(a) for a in arquivos_p(rodada, passo)]


@st.cache_data(show_spinner='Lendo mapa de p…', max_entries=32)
def ler_p(rodada, passo, area, max_lado):
    """p em [0, 1] (NaN sem dado), transform e bounds, reduzido para caber em max_lado px."""
    arqs = arquivos_p(rodada, passo)
    if area != AREA_TODA:
        arqs = [a for a in arqs if celula_de(a) == area]
    if not arqs:
        return None
    with rasterio.open(arqs[0]) as d:
        res = d.res
    fontes = [rasterio.open(a) for a in arqs]
    try:
        esq = min(f.bounds.left for f in fontes)
        dir_ = max(f.bounds.right for f in fontes)
        inf = min(f.bounds.bottom for f in fontes)
        sup = max(f.bounds.top for f in fontes)
        lado = max((dir_ - esq) / res[0], (sup - inf) / res[1])
        fator = max(1.0, lado / max_lado)
        mos, tr = merge(fontes, res=(res[0] * fator, res[1] * fator), nodata=0,
                        resampling=Resampling.average if fator > 1 else Resampling.nearest)
    finally:
        for f in fontes:
            f.close()
    p = mos[0].astype(np.float32)
    p[p == 0] = np.nan
    return p / 10000, tr, (esq, inf, dir_, sup)


def pasta_canais(modelo):
    return TILES_DIR / f'canais_{modelo}'


def modelos_com_canais():
    return [m for m in ('B1', 'A') if pasta_canais(m).is_dir()]


@st.cache_data(max_entries=4)
def nomes_canais(modelo):
    primeiro = next(pasta_canais(modelo).glob('*/*_r*c*.tif'), None)
    if primeiro is None:
        return []
    with rasterio.open(primeiro) as d:
        return [n or f'banda_{i + 1}' for i, n in enumerate(d.descriptions)]


@st.cache_data(show_spinner='Lendo canais dos pedaços…', max_entries=32)
def ler_canais(modelo, celula, canais, coefs, forma):
    """Canais pedidos (C, H, W) na mesma grade do mapa de p (média dos pixels de 30 m).
    coefs = tuple(transform) do mapa de p (o Affine não é hashável para o cache)."""
    transform = Affine(*coefs[:6])
    arqs = sorted((pasta_canais(modelo) / celula).glob(f'{celula}_r*c*.tif'))
    if not arqs:
        return None
    nomes = nomes_canais(modelo)
    indices = [nomes.index(c) + 1 for c in canais]
    a, e = transform.a, transform.e
    esq, sup = transform.c, transform.f
    limites = (esq, sup + e * forma[0], esq + a * forma[1], sup)
    fontes = [rasterio.open(x) for x in arqs]
    try:
        mos, _ = merge(fontes, bounds=limites, res=(a, -e), indexes=indices, nodata=np.nan,
                       resampling=Resampling.average)
    finally:
        for f in fontes:
            f.close()
    saida = np.full((len(canais), *forma), np.nan, np.float32)      # arredondamento de ±1 px
    h, w = min(forma[0], mos.shape[1]), min(forma[1], mos.shape[2])
    saida[:, :h, :w] = mos[:, :h, :w]
    return saida


@st.cache_data(max_entries=2)
def ocorrencias():
    """Uma linha por ocorrência de fosfato (o ponto, sem os centros deslocados)."""
    csv = PATCHES_DIR / 'indice_patches.csv'
    if not csv.exists():
        return pd.DataFrame(columns=['codigo_obj', 'papel', 'lon', 'lat', 'p_rf', 'p_gtb'])
    df = pd.read_csv(csv)
    oc = df[df.origem_amostra.astype(str).str.startswith('ocorrencia')]
    if 'jitter' in oc:
        oc = oc[oc.jitter == 0]
    cols = [c for c in ['codigo_obj', 'papel', 'origem_amostra', 'lon', 'lat', 'p_rf', 'p_gtb'] if c in oc]
    return oc[cols].drop_duplicates('codigo_obj' if 'codigo_obj' in oc else ['lon', 'lat']).reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════
# FIGURA
# ═══════════════════════════════════════════════════════════════════════════
def eixos(transform, forma):
    """Longitudes e latitudes dos centros dos pixels."""
    lon = transform.c + transform.a * (np.arange(forma[1]) + 0.5)
    lat = transform.f + transform.e * (np.arange(forma[0]) + 0.5)
    return lon, lat


def realce(x, pmin, pmax):
    """Corte nos percentis e escala 0–255 por canal; NaN → preto."""
    saida = np.zeros(x.shape, np.uint8)
    for i, c in enumerate(x):
        v = c[np.isfinite(c)]
        if not v.size:
            continue
        lo, hi = np.percentile(v, [pmin, pmax])
        s = np.clip((c - lo) / max(hi - lo, 1e-12), 0, 1)
        saida[i] = np.nan_to_num(s * 255).astype(np.uint8)
    return np.moveaxis(saida, 0, -1)


def heatmap(z, lon, lat, nome, escala, zmin, zmax, barra_x, zmid=None):
    return go.Heatmap(z=np.round(z, 4), x=lon, y=lat, colorscale=escala, zmin=zmin, zmax=zmax, zmid=zmid,
                      name=nome, colorbar=dict(title=nome, x=barra_x, len=0.9),
                      hovertemplate='lon %{x:.4f}<br>lat %{y:.4f}<br>' + nome + ' %{z:.3f}<extra></extra>')


def figura(p, lon, lat, painel, titulo_esq, titulo_dir, oc, altura):
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.07, subplot_titles=(titulo_esq, titulo_dir))
    fig.add_trace(heatmap(p, lon, lat, 'p', 'Viridis', 0, 1, 0.45), 1, 1)

    tipo = painel['tipo']
    if tipo == 'rgb':
        rgb = painel['rgb']
        dx, dy = lon[1] - lon[0], lat[0] - lat[1]
        fig.add_layout_image(dict(source=Image.fromarray(rgb), xref='x2', yref='y2',
                                  x=lon[0] - dx / 2, y=lat[0] + dy / 2,
                                  sizex=dx * rgb.shape[1], sizey=dy * rgb.shape[0],
                                  sizing='stretch', layer='below'))
        # traço invisível: dá extensão aos eixos; hover com os 3 canais só em mapas pequenos
        # (customdata vai inteiro para o navegador)
        valores = painel['valores'] if rgb.shape[0] * rgb.shape[1] <= 400_000 else None
        dica = ('<br>'.join(f'{c} %{{customdata[{i}]:.4g}}' for i, c in enumerate(painel['canais']))
                if valores is not None else 'canais: reduza a resolução para ver os valores')
        fig.add_trace(go.Heatmap(z=np.zeros(rgb.shape[:2]), x=lon, y=lat, opacity=0, showscale=False,
                                 customdata=valores,
                                 hovertemplate='lon %{x:.4f}<br>lat %{y:.4f}<br>' + dica + '<extra></extra>'),
                      1, 2)
    elif tipo == 'diferenca':
        d = np.abs(painel['z'][np.isfinite(painel['z'])])
        lim = float(np.percentile(d, 99)) if d.size else 1.0
        lim = lim or 1.0
        fig.add_trace(heatmap(painel['z'], lon, lat, 'Δp', 'RdBu_r', -lim, lim, 1.0, zmid=0), 1, 2)
    else:
        z = painel['z']
        v = z[np.isfinite(z)]
        zmin, zmax = (0, 1) if tipo == 'rodada' or not v.size else np.percentile(v, [2, 98])
        fig.add_trace(heatmap(z, lon, lat, painel['nome'], painel.get('escala', 'Viridis'),
                              zmin, zmax, 1.0), 1, 2)

    if oc is not None and len(oc):
        cores = {'treino': '#ff2d55', 'validacao': '#ffd60a'}
        for papel, g in oc.groupby(oc.papel.fillna('').astype(str)):
            texto = [f'{r.codigo_obj} · {papel}<br>p CNN {r.p_cnn:.3f}<br>p RF {r.p_rf:.3f} · p GTB {r.p_gtb:.3f}'
                     for r in g.itertuples()]
            for col in (1, 2):
                fig.add_trace(go.Scatter(x=g.lon, y=g.lat, mode='markers', name=f'ocorrência {papel or "?"}',
                                         legendgroup=papel, showlegend=col == 1, text=texto,
                                         hovertemplate='%{text}<extra></extra>',
                                         marker=dict(size=9, color=cores.get(papel, '#ffffff'), symbol='x',
                                                     line=dict(width=1, color='black'))), 1, col)

    latc = float(np.nanmean(lat))
    fig.update_xaxes(showgrid=False, range=[lon[0], lon[-1]])
    fig.update_yaxes(showgrid=False, range=[lat[-1], lat[0]])
    # painel da direita preso ao da esquerda (zoom e arrasto juntos); pixel quadrado em metros
    fig.update_xaxes(matches='x', row=1, col=2)
    fig.update_yaxes(matches='y', row=1, col=2)
    fig.update_yaxes(scaleanchor='x', scaleratio=1 / np.cos(np.radians(latc)), row=1, col=1)
    fig.update_layout(height=altura, margin=dict(l=10, r=10, t=40, b=10), uirevision='mapa',
                      legend=dict(orientation='h', y=-0.05), plot_bgcolor='#202020')
    return fig


# ═══════════════════════════════════════════════════════════════════════════
# INTERFACE
# ═══════════════════════════════════════════════════════════════════════════
def main():
    st.set_page_config(page_title='CNN fosfato — mapas', layout='wide')
    rodadas = listar_rodadas()
    if not rodadas:
        st.error(f'Nenhum mapa em {RUNS_DIR}/<rodada>/mapas_passo<k>/p_<LxCy>.tif — rode inferir_cnn.py primeiro.')
        return

    with st.sidebar:
        st.header('Mapa de p (esquerda)')
        rodada = st.selectbox('Rodada', list(rodadas))
        passo = st.selectbox('Passo', rodadas[rodada], index=len(rodadas[rodada]) - 1,
                             format_func=lambda k: f'{k} ({k * 30} m)')
        area = st.selectbox('Área', celulas(rodada, passo) + [AREA_TODA])
        max_lado = st.slider('Resolução máxima (px por lado)', 400, 2000, 1000, 100,
                             help='Mapas maiores são reduzidos pela média. Mais px = mais lento no navegador.')

        st.header('Painel da direita')
        opcoes = ['RGB artificial', 'Canal / índice', 'Outra rodada', 'Diferença entre rodadas']
        if area == AREA_TODA:
            opcoes = opcoes[2:]
            st.caption('RGB e canais só por célula (a área toda são 900 pedaços).')
        tipo = st.radio('Mostrar', opcoes)

        st.header('Análise')
        limiar = st.slider('Limiar de p', 0.0, 1.0, 0.5, 0.05)
        mascarar = st.checkbox('Esconder p abaixo do limiar', False)
        mostrar_oc = st.checkbox('Ocorrências de fosfato', True)
        altura = st.slider('Altura da figura (px)', 500, 1200, 750, 50)

    lido = ler_p(rodada, passo, area, max_lado)
    if lido is None:
        st.warning('Mapa não encontrado.')
        return
    p, tr, limites = lido
    lon, lat = eixos(tr, p.shape)
    p_vis = np.where(p >= limiar, p, np.nan) if mascarar else p

    # --- painel da direita ---------------------------------------------------
    painel, titulo_dir = None, ''
    if tipo in ('RGB artificial', 'Canal / índice'):
        modelos = modelos_com_canais()
        if not modelos:
            st.error(f'Sem pedaços de canais em {TILES_DIR} (baixar_tiles.py).')
            return
        with st.sidebar:
            modelo = st.selectbox('Canais do modelo', modelos,
                                  help='B1 = log das bandas + emissividade; A = índices espectrais')
            nomes = [n for n in nomes_canais(modelo) if n != 'swir_estimado'] + ['swir_estimado']
            if tipo == 'RGB artificial':
                comps = {k: v for k, v in COMPOSICOES.get(modelo, {}).items() if all(c in nomes for c in v)}
                escolha = st.selectbox('Composição', list(comps) + [PERSONALIZADO])
                if escolha == PERSONALIZADO:
                    canais = tuple(st.selectbox(cor, nomes, index=min(i, len(nomes) - 1))
                                   for i, cor in enumerate(['R', 'G', 'B']))
                else:
                    canais = comps[escolha]
                pmin, pmax = st.slider('Realce (percentis)', 0.0, 100.0, (2.0, 98.0), 0.5)
            else:
                canais = (st.selectbox('Canal', nomes),)
                escala = st.selectbox('Paleta', ['Viridis', 'Cividis', 'Turbo', 'RdBu_r', 'Greys'])
        x = ler_canais(modelo, area, canais, tuple(tr), p.shape)
        if x is None:
            st.error(f'Sem pedaços da célula {area} em {pasta_canais(modelo)}.')
            return
        if tipo == 'RGB artificial':
            painel = {'tipo': 'rgb', 'rgb': realce(x, pmin, pmax), 'canais': canais,
                      'valores': np.moveaxis(np.round(x, 4), 0, -1)}
            titulo_dir = f'{modelo}: R {canais[0]} · G {canais[1]} · B {canais[2]}'
        else:
            painel = {'tipo': 'canal', 'z': x[0], 'nome': canais[0], 'escala': escala}
            titulo_dir = f'{modelo}: {canais[0]}'
    else:
        outras = [r for r in rodadas if r != rodada and passo in rodadas[r]]
        if not outras:
            st.warning(f'Nenhuma outra rodada com mapas de passo {passo}.')
            return
        with st.sidebar:
            outra = st.selectbox('Outra rodada', outras)
        lido2 = ler_p(outra, passo, area, max_lado)
        if lido2 is None or lido2[0].shape != p.shape:
            st.warning(f'{outra} não tem esta área no passo {passo} (ou a grade é diferente).')
            return
        if tipo == 'Outra rodada':
            painel = {'tipo': 'rodada', 'z': lido2[0], 'nome': 'p'}
            titulo_dir = f'{outra} · passo {passo}'
        else:
            painel = {'tipo': 'diferenca', 'z': p - lido2[0]}
            titulo_dir = f'{rodada} − {outra}'

    # --- ocorrências na área, com p da CNN no pixel ---------------------------
    oc = None
    if mostrar_oc:
        todas = ocorrencias()
        dentro = todas[(todas.lon >= limites[0]) & (todas.lon <= limites[2])
                       & (todas.lat >= limites[1]) & (todas.lat <= limites[3])].copy()
        if len(dentro):
            lin, col = rowcol(tr, dentro.lon.values, dentro.lat.values)
            lin = np.clip(np.asarray(lin), 0, p.shape[0] - 1)
            col = np.clip(np.asarray(col), 0, p.shape[1] - 1)
            dentro['p_cnn'] = p[lin, col]
            oc = dentro

    # --- números ---------------------------------------------------------------
    validos = p[np.isfinite(p)]
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Área com p ≥ limiar', f'{(validos >= limiar).mean():.1%}' if validos.size else '—')
    c2.metric('p médio', f'{validos.mean():.3f}' if validos.size else '—')
    if oc is not None:
        acima = (oc.p_cnn >= limiar).sum()
        c3.metric('Ocorrências na área', len(oc))
        c4.metric('Ocorrências com p ≥ limiar', f'{acima} de {len(oc)}')

    titulo_esq = f'{rodada} · passo {passo} · {area}'
    fig = figura(p_vis, lon, lat, painel, titulo_esq, titulo_dir, oc, altura)
    try:
        st.plotly_chart(fig, width='stretch', config={'scrollZoom': True})
    except TypeError:                                   # Streamlit antigo
        st.plotly_chart(fig, use_container_width=True, config={'scrollZoom': True})
    st.caption('Roda do mouse = zoom; arrastar = mover; duplo clique = volta ao início. '
               'Os dois painéis andam juntos.')

    # --- detalhes ----------------------------------------------------------------
    esq, dir_ = st.columns(2)
    with esq:
        st.subheader('Distribuição de p')
        amostra = validos if validos.size <= 200_000 else np.random.default_rng(0).choice(validos, 200_000, False)
        h = go.Figure(go.Histogram(x=amostra, nbinsx=50, marker_color='#3b7dd8'))
        h.add_vline(x=limiar, line_dash='dash', line_color='red')
        h.update_layout(height=280, margin=dict(l=10, r=10, t=10, b=10), xaxis_title='p', yaxis_title='pixels')
        st.plotly_chart(h, config={'displayModeBar': False})
    with dir_:
        st.subheader('Ocorrências na área')
        if oc is not None:
            cols = [c for c in ['codigo_obj', 'papel', 'p_cnn', 'p_rf', 'p_gtb', 'lon', 'lat'] if c in oc]
            st.dataframe(oc[cols].sort_values('p_cnn', ascending=False).round(3), hide_index=True)
        else:
            st.caption('Nenhuma ocorrência nesta área (ou indice_patches.csv não encontrado).')


if __name__ == '__main__':
    main()
