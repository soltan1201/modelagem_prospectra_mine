"""
Harmonização dos mosaicos SAZONAIS SECOS ASTER por casamento de histograma contra o
mosaico Landsat dos mesmos meses (exportar_landsat_sazonal.py). Substitui as passadas
com o MODIS (harmonizar_sazonal_modis.py / _relativo.py): parte do ORIGINAL TOA.
Método e referências: README_harmonizacao_landsat.md.

Unidade de correção = PEDAÇO = cena (cena_t do mapa de cenas) dentro da célula. Todos os
pedaços são casados com a mesma referência contínua, então as emendas somem; a
transformação é uma função de valor → valor, aplicada a 15 m (a textura do ASTER fica).

  1) pixels estáveis ("pouca mudança"), por pedaço: EVI2 (sem azul, o ASTER não tem)
       d = EVI2_ASTER − EVI2_Landsat;  estável = |d − mediana(d) do pedaço| ≤ quantil
       FRAC_ESTAVEL do pedaço. As curvas saem SÓ dos estáveis e valem para todos os pixels.
  2) VNIR: casamento de histograma por percentis, banda a banda
       B01 ↔ green, B02 ↔ red, B3N ↔ nir
  3) SWIR (só A): banda guia B04 ↔ swir1 (casamento de histograma); B05–B09 seguem a guia
       razao   : Bk' = Bk × B04'/B04         (todas as razões entre bandas iguais ao original)
       residuo : Bk' = [Bk − (a_p + b_p·B04) + (a_ref + b_ref·B04)] × B04'/B04
                 (a_p, b_p) = Theil-Sen Bk × B04 nos estáveis do pedaço; (a_ref, b_ref) =
                 mediana entre os pedaços do período. Normaliza a relação de fundo entre
                 bandas (crosstalk, atmosfera) e mantém o resíduo (sinal mineral).
  4) TIR: guia B13 (DN) casada com a temperatura de superfície Landsat (st): por pedaço,
     reta quantil-quantil DN13 = a + b·st; alvo do pedaço = a_ref + b_ref·st_percentis
     (mediana das retas do período). B10, B11, B12, B14 seguem a guia como no SWIR.
     Pedaço sem st suficiente: TIR fica como no original.

Etapas (perguntadas ao rodar):
  1) amostras (GEE → CSV local, sem Drive): por célula, até N_POR_PECA pixels de 60 m
     por pedaço → SAIDA_HM/amostras/<P>_<LxCy>.csv
  2) curvas (local): estáveis, percentis, Theil-Sen, referência → SAIDA_HM/curvas_<P>.json
  3) aplicar: COL_HM[modo]/<mesmo nome do original>, bandas do original +
       cena_t + hm (1 = curva do pedaço, 2 = curva da célula, para pedaço pequeno ou
       pixel fora do mapa de cenas)
  4) relatório (local): estáveis por pedaço e o teste "o ajuste do TIR acompanha o do
     SWIR / o desvio de EVI2?" (correlação de Spearman entre pedaços)

Conta: caatinga04, projeto mapbiomas-brazil (igual às outras harmonizações).

Uso: python harmonizar_sazonal_landsat.py
"""

import json
import logging

import ee
import numpy as np
import pandas as pd

import build_mosaico_sazonal_seco as s
import exportar_landsat_sazonal as L
import harmonizar_sazonal_modis as h
from build_mosaicos_semestrais_v2 import (ler_grid, listar_destino, perguntar_celulas,
                                          perguntar_opcao, perguntar_sim_nao)
from gee_tools import switch_user

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
COL_HM = {'razao': f'{h.PASTA_MINE}/ASTER_sazonal_hm_razao',
          'residuo': f'{h.PASTA_MINE}/ASTER_sazonal_hm_residuo'}
SAIDA_HM = h.SAIDA / 'landsat'

ESC_AMOSTRA  = 60       # m — pixel da amostra (4 × 4 ASTER, 2 × 2 Landsat)
FRAC_MIN     = 0.9      # fração do pixel de 60 m coberta por ASTER válido
N_POR_PECA   = 3000     # pixels sorteados por pedaço
SEMENTE      = 42
FRAC_ESTAVEL = 0.6      # fração dos pixels de cada pedaço tratada como "pouca mudança"
N_MIN_PECA   = 300      # estáveis mínimos para o pedaço ter curva própria
FRAC_MIN_ST  = 0.5      # fração dos estáveis com st válido para corrigir o TIR
PERC         = list(range(1, 100, 2))          # percentis das curvas
I05, I95     = PERC.index(5), PERC.index(95)   # inclinação das caudas
CAUDA        = (1000.0, 5000.0)                # extensão das caudas abaixo/acima (unidades da banda)
N_THEIL      = 800      # pixels por pedaço no Theil-Sen (pares: N²/2)

PARES_VNIR  = {'B01': 'green', 'B02': 'red', 'B3N': 'nir'}
GUIA_SWIR   = ('B04', 'swir1')
SW_OUTRAS   = ['B05', 'B06', 'B07', 'B08', 'B09']
GUIA_TIR    = 'B13'
TIR_OUTRAS  = ['B10', 'B11', 'B12', 'B14']
ST_NULO     = -9999


def bandas_aster(periodo: str) -> list[str]:
    return s.CONFIG[periodo]['bandas_toa'] + s.TIR


def bandas_ls(periodo: str) -> list[str]:
    return L.BANDAS_LS[periodo] + ['st']


def landsat(periodo: str, celula: str) -> ee.Image:
    return ee.Image(f'{L.COL_LS}/{L.nome_ls(periodo, celula)}')


def evi2(nir, red):
    """EVI2 (Jiang et al. 2008) com reflectância × 1000 na entrada."""
    n, r = nir / s.ESCALA_REFL, red / s.ESCALA_REFL
    return 2.5 * (n - r) / (n + 2.4 * r + 1)


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 1 — AMOSTRAS
# ═══════════════════════════════════════════════════════════════════════════
def imagem_amostra(periodo: str, celula: str) -> tuple[ee.Image, object]:
    nome = s.nome_mosaico(periodo, celula)
    ast = ee.Image(f'{s.DESTINO_FINAL}/{nome}')
    cen = ee.Image(f'{h.COL_CENAS}/{nome}')
    ls = landsat(periodo, celula)
    proj = ls.select('red').projection().atScale(ESC_AMOSTRA)
    red = lambda im, r: im.reduceResolution(r, maxPixels=64).reproject(proj)

    ba = bandas_aster(periodo)
    a = ast.select(ba).toFloat().updateMask(cen.select('confere').eq(1))
    t = cen.select('cena_t').updateMask(a.select('B3N').mask())
    a60 = red(a, ee.Reducer.mean()).rename([f'a_{b}' for b in ba])
    frac = red(a.select('B3N').mask().toFloat(), ee.Reducer.mean())
    tmin, tmax = red(t, ee.Reducer.min()), red(t, ee.Reducer.max())

    bl = L.BANDAS_LS[periodo]
    l60 = red(ls.select(bl).toFloat(), ee.Reducer.mean()).rename([f'l_{b}' for b in bl])
    st60 = (red(ls.select('st').toFloat(), ee.Reducer.mean()).divide(100)
              .unmask(ST_NULO).rename('l_st'))
    img = (a60.addBands(l60).addBands(st60).addBands(tmin.toInt32().rename('cena_t'))
              .updateMask(tmin.eq(tmax).And(frac.gte(FRAC_MIN))))
    return img, proj


def amostras_celula(periodo: str, nome_cel: str, geom):
    img, proj = imagem_amostra(periodo, nome_cel)
    fc = img.stratifiedSample(numPoints=N_POR_PECA, classBand='cena_t', region=geom,
                              projection=proj, scale=ESC_AMOSTRA, seed=SEMENTE,
                              dropNulls=True, tileScale=16, geometries=False)
    colunas = (['cena_t'] + [f'a_{b}' for b in bandas_aster(periodo)]
               + [f'l_{b}' for b in L.BANDAS_LS[periodo]] + ['l_st'])
    # getDownloadURL (CSV): o computeFeatures paginado repetia pontos entre páginas
    df = pd.read_csv(fc.getDownloadURL(filetype='csv', selectors=colunas))
    df['celula'] = nome_cel
    return df[['celula'] + colunas]


def arq_amostra(periodo: str, celula: str):
    return SAIDA_HM / 'amostras' / f'{periodo}_{celula}.csv'


def etapa1(periodos, celulas, grid):
    (SAIDA_HM / 'amostras').mkdir(parents=True, exist_ok=True)
    originais = listar_destino(s.DESTINO_FINAL, criar=False)
    cenas = listar_destino(h.COL_CENAS, criar=False)
    ls_feitos = listar_destino(L.COL_LS, criar=False)
    for p in periodos:
        for cid in celulas:
            c = grid[cid]['nome']
            nome = s.nome_mosaico(p, c)
            arq = arq_amostra(p, c)
            if arq.exists():
                log.info(f'  {p} {c}: pulado (já existe {arq.name}; apague para refazer)')
            elif L.nome_ls(p, c) not in ls_feitos:
                log.info(f'  {p} {c}: pulado (sem Landsat — rode exportar_landsat_sazonal.py)')
            elif nome not in originais or nome not in cenas:
                log.info(f'  {p} {c}: pulado (sem mosaico original ou sem mapa de cenas)')
            else:
                try:
                    df = h.com_tentativas(lambda: amostras_celula(p, c, grid[cid]['geom']))
                except Exception as e:  # noqa: BLE001 — segue para a próxima célula
                    log.warning(f'  {p} {c}: erro {str(e)[:150]}')
                    continue
                df.to_csv(arq, index=False)
                log.info(f'  {p} {c}: {len(df)} pixels em {df["cena_t"].nunique()} pedaços → {arq.name}')


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 2 — CURVAS (local)
# ═══════════════════════════════════════════════════════════════════════════
def estaveis(df: pd.DataFrame) -> pd.DataFrame:
    d = evi2(df['a_B3N'], df['a_B02']) - evi2(df['l_nir'], df['l_red'])
    dev = (d - d.groupby(df['cena_t']).transform('median')).abs()
    lim = dev.groupby(df['cena_t']).transform(lambda v: v.quantile(FRAC_ESTAVEL))
    out = df.assign(d_evi2=d)
    return out[dev <= lim]


def curva(x: np.ndarray, y: np.ndarray) -> list[list[float]]:
    """Percentis do ASTER → percentis da referência (pixels estáveis do pedaço)."""
    return curva_de_quantis(np.percentile(x, PERC), np.percentile(y, PERC))


def theil_sen(x: np.ndarray, y: np.ndarray, rng) -> list[float]:
    """y = a + b·x por Theil-Sen (mediana das inclinações entre pares). → [a, b]"""
    if len(x) > N_THEIL:
        i = rng.choice(len(x), N_THEIL, replace=False)
        x, y = x[i], y[i]
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    ok = np.triu(dx != 0, 1)
    b = float(np.median(dy[ok] / dx[ok])) if ok.any() else 0.0
    a = float(np.median(y - b * x))
    return [round(a, 4), round(b, 6)]


def transformacao(g: pd.DataFrame, periodo: str, rng) -> dict:
    T = {'n': int(len(g)), 'd_evi2': round(float(g['d_evi2'].median()), 4), 'curvas': {}, 'rel': {}}
    for b, lb in PARES_VNIR.items():
        T['curvas'][b] = curva(g[f'a_{b}'].values, g[f'l_{lb}'].values)
    if periodo == 'A':
        gb, lb = GUIA_SWIR
        T['curvas'][gb] = curva(g[f'a_{gb}'].values, g[f'l_{lb}'].values)
        for k in SW_OUTRAS:
            T['rel'][k] = theil_sen(g[f'a_{gb}'].values, g[f'a_{k}'].values, rng)
    st = g[g['l_st'] > ST_NULO + 1]
    T['tir_ok'] = bool(len(st) >= FRAC_MIN_ST * len(g) and len(st) >= N_MIN_PECA)
    if T['tir_ok']:
        q_st = np.percentile(st['l_st'].values, PERC)
        q_13 = np.percentile(st[f'a_{GUIA_TIR}'].values, PERC)
        b, a = np.polyfit(q_st, q_13, 1)          # reta quantil-quantil DN13 = a + b·st
        T['st'] = {'a': float(a), 'b': float(b), 'q_st': [float(v) for v in q_st],
                   'q_13': [float(v) for v in q_13]}
        for k in TIR_OUTRAS:
            T['rel'][k] = theil_sen(st[f'a_{GUIA_TIR}'].values, st[f'a_{k}'].values, rng)
    return T


def fechar_tir(T: dict, ref: dict) -> None:
    """Curva do B13: percentis do pedaço → o que a reta de referência daria com o st do pedaço."""
    if not T.get('tir_ok'):
        return
    st = T.pop('st')
    alvo = ref['st_a'] + ref['st_b'] * np.array(st['q_st'])
    T['curvas'][GUIA_TIR] = curva_de_quantis(np.array(st['q_13']), alvo)


def curva_de_quantis(qx: np.ndarray, qy: np.ndarray) -> list[list[float]]:
    """[x, y] para ee.Image.interpolate: x estritamente crescente (exigido), y não
    decrescente, com caudas lineares (inclinação entre os percentis 5 e 95) para os
    valores fora da faixa dos estáveis — anomalias claras/escuras não são cortadas."""
    qx = np.maximum.accumulate(qx) + np.arange(len(qx)) * 1e-3
    qy = np.maximum.accumulate(qy)
    dx = qx[I95] - qx[I05]
    g = (qy[I95] - qy[I05]) / dx if dx > 0 else 1.0
    xs = [qx[0] - CAUDA[0]] + list(qx) + [qx[-1] + CAUDA[1]]
    ys = [qy[0] - g * CAUDA[0]] + list(qy) + [qy[-1] + g * CAUDA[1]]
    return [[round(float(v), 4) for v in xs], [round(float(v), 4) for v in ys]]


def arq_curvas(periodo: str):
    return SAIDA_HM / f'curvas_{periodo}.json'


def etapa2(periodos):
    rng = np.random.default_rng(SEMENTE)
    for p in periodos:
        arqs = sorted((SAIDA_HM / 'amostras').glob(f'{p}_*.csv'))
        if not arqs:
            log.warning(f'Período {p}: sem amostras — rode a etapa 1.')
            continue
        df = estaveis(pd.concat([pd.read_csv(a) for a in arqs], ignore_index=True))
        celulas = {}
        for c, gc in df.groupby('celula'):
            if len(gc) < N_MIN_PECA:
                log.warning(f'  {p} {c}: só {len(gc)} pixels estáveis — célula fica sem curva')
                continue
            pecas = {str(int(t)): transformacao(gt, p, rng)
                     for t, gt in gc.groupby('cena_t') if len(gt) >= N_MIN_PECA}
            celulas[c] = {'*': transformacao(gc, p, rng), 'pecas': pecas,
                          'pecas_pequenas': int(gc['cena_t'].nunique() - len(pecas))}
        # referência do período: mediana entre os pedaços com curva própria
        todas = [T for cel in celulas.values() for T in cel['pecas'].values()]
        ref = {'rel': {}}
        for k in (SW_OUTRAS if p == 'A' else []) + TIR_OUTRAS:
            v = np.array([T['rel'][k] for T in todas if k in T['rel']])
            if len(v):
                ref['rel'][k] = [round(float(np.median(v[:, 0])), 4), round(float(np.median(v[:, 1])), 6)]
        st = [T['st'] for T in todas if T.get('tir_ok')]
        if st:
            ref['st_a'] = float(np.median([x['a'] for x in st]))
            ref['st_b'] = float(np.median([x['b'] for x in st]))
        for cel in celulas.values():
            for T in [cel['*'], *cel['pecas'].values()]:
                if 'st_a' in ref:
                    fechar_tir(T, ref)
                else:
                    T['tir_ok'] = False
                    T.pop('st', None)
        n_p = sum(len(c['pecas']) for c in celulas.values())
        n_peq = sum(c['pecas_pequenas'] for c in celulas.values())
        arq_curvas(p).write_text(json.dumps({'periodo': p, 'ref': ref, 'celulas': celulas}))
        log.info(f'Período {p}: {len(celulas)} células, {n_p} pedaços com curva própria, '
                 f'{n_peq} pequenos (usam a da célula) → {arq_curvas(p)}')
        log.info(f'  referência: ' + ', '.join(f'{k} a={a:.1f} b={b:.3f}' for k, (a, b) in ref['rel'].items())
                 + (f'; TIR DN13 = {ref["st_a"]:.1f} + {ref["st_b"]:.2f}·st' if 'st_a' in ref else '; TIR sem st'))


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 3 — APLICAR
# ═══════════════════════════════════════════════════════════════════════════
def casar(img: ee.Image, banda: str, cv: list) -> ee.Image:
    return img.select(banda).toFloat().interpolate(cv[0], cv[1], 'clamp')


def seguir_guia(img, guia: str, guia2, outras: list[str], T: dict, ref: dict, modo: str):
    """Bandas que seguem a guia: razão (e, no modo residuo, relação de fundo normalizada)."""
    g = img.select(guia).toFloat()
    razao = guia2.divide(g.max(1))
    saida = []
    for k in outras:
        bk = img.select(k).toFloat()
        if modo == 'residuo' and k in T['rel'] and k in ref['rel']:
            a_p, b_p = T['rel'][k]
            a_r, b_r = ref['rel'][k]
            bk = bk.subtract(g.multiply(b_p).add(a_p)).add(g.multiply(b_r).add(a_r))
        saida.append(bk.multiply(razao).rename(k))
    return saida


def corrigir(img: ee.Image, T: dict, ref: dict, periodo: str, modo: str) -> ee.Image:
    out = [casar(img, b, T['curvas'][b]).rename(b) for b in PARES_VNIR]
    if periodo == 'A':
        gb = GUIA_SWIR[0]
        g2 = casar(img, gb, T['curvas'][gb])
        out += [g2.rename(gb)] + seguir_guia(img, gb, g2, SW_OUTRAS, T, ref, modo)
    if T.get('tir_ok') and GUIA_TIR in T['curvas']:
        t2 = casar(img, GUIA_TIR, T['curvas'][GUIA_TIR])
        out += [t2.rename(GUIA_TIR)] + seguir_guia(img, GUIA_TIR, t2, TIR_OUTRAS, T, ref, modo)
    else:
        out += [img.select(s.TIR).toFloat()]
    return ee.Image.cat(out).select(bandas_aster(periodo)).round().toInt32()


def harmonizar(periodo: str, celula: str, C: dict, ref: dict, modo: str) -> ee.Image:
    nome = s.nome_mosaico(periodo, celula)
    ast = ee.Image(f'{s.DESTINO_FINAL}/{nome}')
    t = ee.Image(f'{h.COL_CENAS}/{nome}').select('cena_t')
    extras = lambda im, v: (ast.select('n_obs').updateMask(im.select('B01').mask())
                              .addBands(t.updateMask(im.select('B01').mask()).rename('cena_t'))
                              .addBands(ee.Image.constant(v).toInt32().rename('hm')
                                          .updateMask(im.select('B01').mask())))
    base = corrigir(ast, C['*'], ref, periodo, modo)
    imgs = [base.addBands(extras(ast, 2))]
    for tv, T in C['pecas'].items():
        peca = ast.updateMask(t.eq(int(tv)))
        imgs.append(corrigir(peca, T, ref, periodo, modo).addBands(extras(peca, 1)))
    saida = ee.ImageCollection.fromImages(imgs).mosaic()   # pedaços por cima da curva da célula
    return (ee.Image(saida.copyProperties(ast, ast.propertyNames()))
              .set({'harmonizacao': 'harmonizar_sazonal_landsat.py: casamento de histograma por '
                                    'pedaço (cena × célula) contra o Landsat 32 dias mediano seco',
                    'modo': modo, 'landsat': f'{L.COL_LS}/{L.nome_ls(periodo, celula)}',
                    'SCALE_FACTORS': 'VNIR/SWIR na escala da reflectância de superfície Landsat x1000; '
                                     'TIR em DN normalizado pela st Landsat',
                    'hm': '1 = curva do pedaço, 2 = curva da célula',
                    'pecas': len(C['pecas'])}))


def etapa3(periodos, celulas, grid, modo: str):
    col = COL_HM[modo]
    feitos = listar_destino(col, criar=True)
    ativas = h.tasks_ativas()
    fila = []
    for p in periodos:
        if not arq_curvas(p).exists():
            log.warning(f'Período {p}: {arq_curvas(p)} não existe — rode a etapa 2.')
            continue
        dados = json.loads(arq_curvas(p).read_text())
        for cid in celulas:
            c = grid[cid]['nome']
            nome = s.nome_mosaico(p, c)
            desc = f'hm_{modo[:3]}_{p}_{c}'
            if c not in dados['celulas']:
                motivo = 'sem curvas'
            elif nome in feitos:
                motivo = 'já existe'
            elif desc in ativas:
                motivo = 'processando'
            else:
                fila.append((p, cid, nome, desc, dados))
                continue
            log.info(f'  {p} {c}: pulado ({motivo})')
    if not fila or not perguntar_sim_nao(f'Exportar {len(fila)} mosaicos ({modo})?', padrao=True):
        return
    for p, cid, nome, desc, dados in fila:
        c = grid[cid]['nome']
        img = harmonizar(p, c, dados['celulas'][c], dados['ref'], modo)
        h.exportar_imagem(img, desc, f'{col}/{nome}', grid[cid]['geom'],
                          {'.default': 'mean', 'cena_t': 'mode', 'hm': 'mode', 'n_obs': 'mode'})


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 4 — RELATÓRIO (local)
# ═══════════════════════════════════════════════════════════════════════════
def no_meio(cv: list) -> tuple[float, float]:
    """(x, y) no percentil 50 da curva (as caudas ocupam a 1ª e a última posição)."""
    i = 1 + PERC.index(51) if 51 in PERC else 1 + len(PERC) // 2
    return cv[0][i], cv[1][i]


def etapa4(periodos):
    for p in periodos:
        if not arq_curvas(p).exists():
            print(f'\nPeríodo {p}: {arq_curvas(p)} não existe.')
            continue
        dados = json.loads(arq_curvas(p).read_text())
        linhas = []
        for c, C in dados['celulas'].items():
            for tv, T in C['pecas'].items():
                lin = {'celula': c, 'cena_t': int(tv), 'n': T['n'], 'd_evi2': T['d_evi2']}
                for b in list(PARES_VNIR) + ([GUIA_SWIR[0]] if p == 'A' else []) + [GUIA_TIR]:
                    if b in T['curvas']:
                        x, y = no_meio(T['curvas'][b])
                        lin[f'ganho_{b}'] = y / x if x else np.nan
                        lin[f'desloc_{b}'] = y - x
                linhas.append(lin)
        df = pd.DataFrame(linhas)
        print(f'\n══ Período {p}: {len(df)} pedaços com curva própria em {len(dados["celulas"])} células ══')
        resumo = df.groupby('celula').agg(pecas=('cena_t', 'size'), estaveis=('n', 'sum'),
                                          d_evi2=('d_evi2', 'median'))
        print(resumo.to_string())
        cols = [c for c in df.columns if c.startswith('ganho_')]
        print('\nGanho no percentil 50 (curva / ASTER), por banda — mediana [p10–p90]:')
        for c in cols:
            print(f'  {c[6:]:4s} {df[c].median():.3f} [{df[c].quantile(.1):.3f}–{df[c].quantile(.9):.3f}]')
        if f'desloc_{GUIA_TIR}' in df:
            alvo = f'desloc_{GUIA_TIR}'
            expl = [c for c in ['d_evi2', f'ganho_{GUIA_SWIR[0]}', 'ganho_B3N'] if c in df]
            sub = df.dropna(subset=[alvo])
            print(f'\nTeste TIR × SWIR ({len(sub)} pedaços): Spearman do deslocamento do B13 com')
            for c in expl:
                print(f'  {c:12s} ρ = {sub[alvo].corr(sub[c], method="spearman"):+.2f}')
            print('  |ρ| ≥ 0,5: vale usar o ajuste do SWIR/EVI2 para prever o do TIR onde falta st.')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    print('\n══════ Harmonização dos mosaicos sazonais secos contra o Landsat ══════')
    etapa = perguntar_opcao('Etapa:', [
        ('amostras por pedaço (GEE → CSV local)', 1),
        ('curvas a partir das amostras (local)', 2),
        ('aplicar as curvas (exporta assets)', 3),
        ('relatório e teste TIR × SWIR (local)', 4),
    ])
    periodos = perguntar_opcao('Períodos:', [('A e B', ['A', 'B']), ('só A', ['A']), ('só B', ['B'])])
    if etapa == 2:
        return etapa2(periodos)
    if etapa == 4:
        return etapa4(periodos)
    switch_user(h.CONTA)
    ee.Initialize(project=h.PROJETO)
    log.info(f'Earth Engine com {h.CONTA} (projeto {h.PROJETO})')
    grid = ler_grid()
    celulas = perguntar_celulas(grid)
    if etapa == 1:
        etapa1(periodos, celulas, grid)
    else:
        modo = perguntar_opcao('Bandas que seguem a guia (SWIR B05–B09, TIR):', [
            ('razao — mantém todas as razões do original', 'razao'),
            ('residuo — normaliza a relação de fundo e mantém o resíduo', 'residuo'),
        ])
        etapa3(periodos, celulas, grid, modo)


if __name__ == '__main__':
    main()
