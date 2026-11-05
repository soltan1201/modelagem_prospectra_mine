"""
Harmonização radiométrica dos mosaicos SAZONAIS SECOS (build_mosaico_sazonal_seco.py)
contra o MODIS do mesmo dia — corrige as emendas entre cenas sem refazer a escolha
de pixels do qualityMosaic.

Ideia: o ASTER e o MODIS estão no mesmo satélite (Terra), então cada cena ASTER tem
um MODIS da mesma data, mesma hora, mesma atmosfera e mesmo sol. Para cada cena
ajusta-se, em pixels homogêneos a 500 m, uma reta MODIS = ganho × ASTER + offset
por banda; a reta da cena é aplicada aos pixels do mosaico que vieram dela.
Diferenças de atmosfera, calibração e iluminação entre datas somem; diferenças de
fenologia (parcela verde numa cena, seca na outra) não.

Referência: MODIS/061/MCD43A4 (NBAR, Terra + Aqua, janela de 16 dias centrada no
dia, ajustado para visada nadir — o ASTER é quase nadir). Só pixels com
inversão BRDF completa (qualidade obrigatória = 0).

  ASTER       MODIS             ajuste
  B01 verde   banda 4           ganho + offset
  B02 verm.   banda 1           ganho + offset
  B3N NIR     banda 2           ganho + offset
  B04         banda 6           ganho + offset           (só período A)
  B05–B09     banda 7 (com B05) UM ganho, sem offset     (só período A)
O ganho único sem offset em B05–B09 preserva as razões entre essas bandas, de onde
saem os índices Al-OH e carbonato. TIR (DN) não é corrigido: o MCD43A4 não tem
termal. Saída fica na escala do NBAR (reflectância de superfície × 1000), não mais
TOA — por isso A e B têm de passar os dois pela harmonização.

Etapas (perguntadas ao rodar; uma task por item):
  1) mapa de cenas: reconstrói o qualityMosaic com as mesmas cenas e listas
     (colecao_seca do build) e exporta, por período × célula, as bandas
       cena_t  : segundos desde 2000-01-01 do system:time_start da cena vencedora
       confere : 1 onde o mosaico reconstruído é IGUAL ao asset salvo (todas as
                 bandas), 0 onde não (ex.: blacklist mudou depois do export)
     → COL_CENAS/ASTER_sazonal_<P>_seco_quality_<LxCy>
  2) pares: uma task Export.table.toDrive por período × mês (roda em batch, em
     paralelo). Por cena, só TOA (sem a máscara de nuvem do build, que era o gargalo),
     sorteia até N_CLASSES × N_POR_CLASSE pixels de 500 m homogêneos estratificados
     pelo vermelho MODIS → Drive/PASTA_DRIVE/<P>_<ano>_<mês>.csv. Baixe os CSVs para
     PASTA_PONTOS.
  3) retas (local): lê os CSVs, tira pares com nuvem/sombra, ajusta a reta de cada
     cena em duas passadas → src/Dados/ASTER_sazonal_harmonizacao/coef_modis_<P>.csv
  4) aplicar: precisa de 1 e 3 prontos. Mosaico salvo × reta da cena do pixel.
     Cena com ajuste ruim (poucos pixels, R² baixo, ganho fora da faixa), sem pares ou
     pixel com confere = 0 usa a mediana das retas boas do período.
     → COL_HARM/ASTER_sazonal_<P>_seco_quality_<LxCy> (mesmos nomes do original),
       bandas do original + cena_t + harmonizado (1 = reta da cena, 2 = mediana)
  5) status: situação dos assets, fração de confere por célula e resumo das retas.

Os originais em DESTINO_FINAL não são alterados. Para refazer um item, apague o
asset dele antes (o script pula o que já existe).

Conta: caatinga04 (lê a coleção ASTER pública e lê/escreve em mapbiomas-arida),
rodando no projeto PROJETO (mapbiomas-brazil), não no projeto padrão dela.

Uso: python harmonizar_sazonal_modis.py
"""

import logging
import statistics
import time
from pathlib import Path

import ee
import numpy as np
import pandas as pd

import build_mosaico_sazonal_seco as s
from build_mosaicos_semestrais_v2 import (CRS_EXPORT, CRS_TRANSFORM, GRID_ASSET, ler_grid,
                                          ler_ids, listar_destino, perguntar_celulas,
                                          perguntar_opcao, perguntar_sim_nao)
from gee_tools import switch_user

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CONTA      = 'caatinga04'
PROJETO    = 'mapbiomas-brazil'   # projeto em que a conta roda (cota/billing), não o padrão dela
PASTA_MINE = s.DESTINO_FINAL.rsplit('/', 1)[0]            # projects/mapbiomas-arida/mine
COL_CENAS  = f'{PASTA_MINE}/ASTER_sazonal_cenas'
COL_HARM   = f'{PASTA_MINE}/ASTER_sazonal_harmonizado'
SAIDA      = s.PASTA.parent / 'Dados' / 'ASTER_sazonal_harmonizacao'
PASTA_PONTOS = SAIDA / 'ASTER_sazonal_harmonizacao'   # CSVs baixados de Drive/PASTA_DRIVE

MODIS       = 'MODIS/061/MCD43A4'
MODIS_ESC   = 0.0001
T0_MS       = 946684800000          # 2000-01-01T00:00:00Z em ms (origem de cena_t)

# chave do ajuste → (banda ASTER usada no ajuste, banda MODIS, bandas corrigidas, com offset)
AJUSTES = {
    'B01': ('B01', 4, ['B01'], True),
    'B02': ('B02', 1, ['B02'], True),
    'B3N': ('B3N', 2, ['B3N'], True),
    'B04': ('B04', 6, ['B04'], True),
    'SW2': ('B05', 7, ['B05', 'B06', 'B07', 'B08', 'B09'], False),
}
CHAVES = {'A': ['B01', 'B02', 'B3N', 'B04', 'SW2'], 'B': ['B01', 'B02', 'B3N']}

# Pixels de 500 m usados no ajuste
ESCALA_INTER = 90      # m — passo intermediário da agregação ASTER → MODIS
CV_MAX      = 0.15     # desvio/média do B02 e do B3N (pixels de 90 m) dentro do pixel MODIS
FRAC_MIN    = 0.9      # fração do pixel MODIS coberta por ASTER válido
K_RESIDUO   = 3.0      # 2ª passada descarta |resíduo| > K × RMSE da 1ª
# Amostra estratificada por cena (pixels de 500 m)
CLASSE_MODIS = 1       # banda MODIS que define as faixas (1 = vermelho)
FAIXA_CLASSE = (0, 10000)   # valores brutos do MCD43A4 (× 0,0001 = reflectância)
N_CLASSES    = 10      # faixas iguais de FAIXA_CLASSE
N_POR_CLASSE = 200     # pixels sorteados por faixa (faixa com menos pixels entra inteira)
SEMENTE      = 42
# Reta da cena aceita
N_MIN       = 50
R2_MIN      = 0.5
GANHO_FAIXA = (0.5, 2.0)
DIF_NUVEM   = 0.10     # x_B01 − y_B01 acima disso: nuvem no ASTER (TOA − NBAR no verde)
DIF_SOMBRA  = 0.10     # y_B3N − x_B3N acima disso: sombra no ASTER
PASTA_DRIVE = 'ASTER_sazonal_harmonizacao'          # pasta no Drive dos pares (etapa 2)
ESCALA_INT  = 1e6      # ganho/offset × 1e6 no remap (remap trabalha com inteiros)


def banda_modis(n: int) -> str:
    return f'Nadir_Reflectance_Band{n}'


def qa_modis(n: int) -> str:
    return f'BRDF_Albedo_Band_Mandatory_Quality_Band{n}'


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 1 — MAPA DE CENAS
# ═══════════════════════════════════════════════════════════════════════════
def listas(periodo: str) -> tuple[ee.List, ee.List]:
    """Blacklist e zero-cloud do período, com a mesma regra do build."""
    black = ler_ids(s.ARQ_BLACKLIST)
    zero = ler_ids(s.ARQ_ZERO_CLOUD) - black
    return (ee.List(s.ids_do_periodo(black, periodo)),
            ee.List(s.ids_do_periodo(zero, periodo)))


def cena_t(img) -> ee.Number:
    return ee.Number(img.get('system:time_start')).subtract(T0_MS).divide(1000).round()


def mapa_cenas(periodo: str, geom, salvo: ee.Image, black, zero) -> ee.Image:
    cfg = s.CONFIG[periodo]
    _, _, _, col = s.colecao_seca(periodo, geom, black, zero)
    col = col.map(lambda i: i.addBands(
        ee.Image.constant(cena_t(i)).toInt32().rename('cena_t')
          .updateMask(i.select('quality').mask())))
    m = col.qualityMosaic('quality')
    # mesma conversão de montar_seco: reflectância × 1000 e TIR, arredondados
    recon = (m.select(cfg['bandas_toa']).multiply(s.ESCALA_REFL)
              .addBands(m.select(s.TIR)).round().toInt32())
    confere = recon.eq(salvo.select(cfg['bandas_toa'] + s.TIR)).reduce(ee.Reducer.min())
    return (m.select('cena_t')
             .addBands(confere.rename('confere').toInt32())
             .clip(geom)
             .set({'periodo': periodo, 'cena_t': f'segundos desde 2000-01-01 (T0_MS={T0_MS})'}))


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 2 — RETA POR CENA CONTRA O MCD43A4
# ═══════════════════════════════════════════════════════════════════════════
def pares_cena(img, periodo: str, proj_modis, area) -> ee.FeatureCollection:
    """Pares ASTER × MODIS de uma cena em pixels de 500 m homogêneos, por amostragem
    estratificada: a banda MODIS CLASSE_MODIS (0–10000) é cortada em N_CLASSES faixas
    iguais em FAIXA_CLASSE e sorteiam-se até N_POR_CLASSE pixels em cada faixa. As
    faixas espalham os pontos por toda a amplitude de brilho (melhor para a reta) e o
    A reta é ajustada no Python (etapa 3)."""
    img = ee.Image(img)
    regiao = img.geometry().intersection(area, ee.ErrorMargin(1000))
    chaves = CHAVES[periodo]
    data = ee.Date(img.get('system:time_start'))
    nums = sorted({AJUSTES[k][1] for k in chaves} | {CLASSE_MODIS})
    nomes_m = [banda_modis(n) for n in nums] + [qa_modis(n) for n in nums]
    col_m = ee.ImageCollection(MODIS).filterDate(data.format('YYYY-MM-dd'), data.advance(1, 'day'))
    vazio = ee.Image.constant([0] * len(nomes_m)).rename(nomes_m).updateMask(0)
    modis = ee.Image(ee.Algorithms.If(col_m.size().gt(0), col_m.first(), vazio))

    # ASTER → grade MODIS em dois passos (15/30 m → 90 m → 500 m): direto de 15 m
    # para 500 m são ~1000 pixels por saída e a task estoura memória.
    # Heterogeneidade = desvio dos pixels de 90 m dentro do pixel MODIS.
    p15 = img.select('B01').projection()
    p90 = p15.atScale(ESCALA_INTER)
    para_90 = lambda im: (im.reduceResolution(ee.Reducer.mean(), maxPixels=64)
                            .reproject(p90))
    bandas = s.VNIR + (['B04', 'B05'] if periodo == 'A' else [])
    m90 = para_90(img.select(s.VNIR).setDefaultProjection(p15))
    if periodo == 'A':
        sw = img.select(['B04', 'B05'])
        m90 = m90.addBands(para_90(sw.setDefaultProjection(sw.select('B04').projection())))
    f90 = para_90(img.select('B3N').mask().setDefaultProjection(p15))
    media = m90.select(bandas).reduceResolution(ee.Reducer.mean(), maxPixels=64).reproject(proj_modis)
    desvio = (m90.select(['B02', 'B3N']).reduceResolution(ee.Reducer.stdDev(), maxPixels=64)
                 .reproject(proj_modis))
    cv = desvio.divide(media.select(['B02', 'B3N'])).reduce(ee.Reducer.max())
    frac = f90.reduceResolution(ee.Reducer.mean(), maxPixels=64).reproject(proj_modis)
    homog = cv.lt(CV_MAX).And(frac.gte(FRAC_MIN))

    qa_ok = ee.Image.cat([modis.select(qa_modis(n)).eq(0) for n in nums]).reduce(ee.Reducer.min())
    x = media.select([AJUSTES[k][0] for k in chaves], [f'x_{k}' for k in chaves])
    y = (modis.select([banda_modis(AJUSTES[k][1]) for k in chaves], [f'y_{k}' for k in chaves])
              .multiply(MODIS_ESC))
    classe = (modis.select(banda_modis(CLASSE_MODIS)).unitScale(*FAIXA_CLASSE)
                   .multiply(N_CLASSES).floor().clamp(0, N_CLASSES - 1).toInt().rename('classe'))
    meta = {'cena_t': cena_t(img), 'cena_id': img.get('system:index'),
            'data': data.format('YYYY-MM-dd'), 'CLOUDCOVER': img.get('CLOUDCOVER')}
    return (x.addBands(y).addBands(classe).updateMask(homog.And(qa_ok))
             .stratifiedSample(numPoints=N_POR_CLASSE, classBand='classe', region=regiao,
                               projection=proj_modis, seed=SEMENTE, dropNulls=True,
                               tileScale=16, geometries=False)
             .map(lambda f: f.set(meta)))


def colunas_pontos(periodo: str) -> list[str]:
    chaves = CHAVES[periodo]
    return (['cena_t', 'cena_id', 'data', 'CLOUDCOVER', 'classe']
            + [f'x_{k}' for k in chaves] + [f'y_{k}' for k in chaves])


def _ols(x: np.ndarray, y: np.ndarray, com_offset: bool) -> dict[str, float]:
    n = len(x)
    if n < 3:
        return {'n': n, 'g': np.nan, 'o': np.nan, 'r2': np.nan, 'rmse': np.nan}
    if com_offset:
        g, o = np.polyfit(x, y, 1)
    else:
        g, o = float(x @ y / (x @ x)), 0.0
    res = y - (g * x + o)
    return {'n': n, 'g': float(g), 'o': float(o),
            'r2': float(np.corrcoef(x, y)[0, 1] ** 2), 'rmse': float(np.sqrt(np.mean(res ** 2)))}


def ajustar(pts: pd.DataFrame, periodo: str) -> dict[str, float]:
    """Reta por chave em duas passadas: a 2ª descarta |resíduo| > K_RESIDUO × RMSE da 1ª
    (nuvem que escapou, borda, mudança dentro da janela do MCD43A4)."""
    # ASTER sem máscara de nuvem: par com ASTER muito mais claro que o MODIS no verde
    # (nuvem) ou muito mais escuro no NIR (sombra) sai antes do ajuste
    nuvem = pts['x_B01'] - pts['y_B01'] > DIF_NUVEM
    sombra = pts['y_B3N'] - pts['x_B3N'] > DIF_SOMBRA
    pts = pts[~(nuvem | sombra)]
    saida = {'n_pontos': len(pts), 'n_nuvem_sombra': int((nuvem | sombra).sum())}
    for k in CHAVES[periodo]:
        com_offset = AJUSTES[k][3]
        par = pts[[f'x_{k}', f'y_{k}']].dropna().to_numpy(dtype=float)
        x, y = par[:, 0], par[:, 1]
        r = _ols(x, y, com_offset)
        if r['n'] >= 3:
            manter = np.abs(y - (r['g'] * x + r['o'])) <= K_RESIDUO * r['rmse']
            r = _ols(x[manter], y[manter], com_offset)
        saida.update({f'{m}_{k}': v for m, v in r.items()})
    return saida


def colecao_coef(periodo: str):
    """Cenas secas do período que tocam a grid (as mesmas do mosaico) e a função que
    sorteia os pares de uma. Só TOA, SEM a máscara de nuvem/sombra do build: o
    focalMax da sombra (até 2,3 km a 15 m) era o que deixava o cálculo lento. Nuvem e
    sombra saem depois, nos pares (DIF_NUVEM / DIF_SOMBRA) e no corte de resíduo."""
    black, zero = listas(periodo)
    area = ee.FeatureCollection(GRID_ASSET).geometry().bounds()
    _, filtrada, _, _ = s.colecao_seca(periodo, area, black, zero)
    bandas_toa = s.CONFIG[periodo]['bandas_toa']
    col = filtrada.map(lambda i: s.radiance_to_toa(i, bandas_toa))
    proj_modis = ee.ImageCollection(MODIS).first().select(banda_modis(1)).projection()
    return col, lambda i: pares_cena(i, periodo, proj_modis, area)


def colunas(periodo: str) -> list[str]:
    return (['cena_t', 'cena_id', 'data', 'CLOUDCOVER', 'n_pontos', 'n_nuvem_sombra']
            + [f'{m}_{k}' for k in CHAVES[periodo] for m in ('n', 'g', 'o', 'r2', 'rmse')])


def arq_coef(periodo: str) -> Path:
    return SAIDA / f'coef_modis_{periodo}.csv'


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 3 — APLICAR
# ═══════════════════════════════════════════════════════════════════════════
def retas(periodo: str) -> tuple[pd.DataFrame, dict[str, tuple[float, float]]]:
    """Tabela de retas com a coluna `valida` e a mediana das retas válidas (reserva)."""
    df = ler_coef(periodo)
    chaves = CHAVES[periodo]
    valida = pd.Series(True, index=df.index)
    for k in chaves:
        valida &= ((df[f'n_{k}'] >= N_MIN) & (df[f'r2_{k}'] >= R2_MIN)
                   & df[f'g_{k}'].between(*GANHO_FAIXA))
    df['valida'] = valida
    if not valida.any():
        raise SystemExit(f'Período {periodo}: nenhuma cena com reta válida — confira a tabela.')
    reserva = {k: (statistics.median(df.loc[valida, f'g_{k}']),
                   statistics.median(df.loc[valida, f'o_{k}'])) for k in chaves}
    return df, reserva


def harmonizar(periodo: str, salvo: ee.Image, cenas: ee.Image,
               df: pd.DataFrame, reserva: dict) -> ee.Image:
    cfg = s.CONFIG[periodo]
    t = cenas.select('cena_t')
    de = [int(v) for v in df['cena_t']]
    usa_cena = t.remap(de, [1 if v else 2 for v in df['valida']], 2)
    usa_cena = usa_cena.where(cenas.select('confere').eq(0), 2).unmask(2)

    corrigidas = []
    for k in CHAVES[periodo]:
        _, _, bandas, _ = AJUSTES[k]
        g_res, o_res = reserva[k]
        g = [r[f'g_{k}'] if r['valida'] else g_res for _, r in df.iterrows()]
        o = [r[f'o_{k}'] if r['valida'] else o_res for _, r in df.iterrows()]
        img_g = (t.remap(de, [round(v * ESCALA_INT) for v in g], round(g_res * ESCALA_INT))
                  .where(usa_cena.eq(2), round(g_res * ESCALA_INT)).unmask(round(g_res * ESCALA_INT))
                  .divide(ESCALA_INT))
        img_o = (t.remap(de, [round(v * ESCALA_INT) for v in o], round(o_res * ESCALA_INT))
                  .where(usa_cena.eq(2), round(o_res * ESCALA_INT)).unmask(round(o_res * ESCALA_INT))
                  .divide(ESCALA_INT))
        # o offset foi ajustado em reflectância 0–1; o asset está × 1000
        corrigidas.append(salvo.select(bandas).multiply(img_g)
                               .add(img_o.multiply(s.ESCALA_REFL)).round().toInt32())
    saida = (ee.Image.cat(corrigidas).select(cfg['bandas_toa'])
              .addBands(salvo.select(s.TIR + ['n_obs']))
              .addBands(t.rename('cena_t'))
              .addBands(usa_cena.updateMask(salvo.select('B01').mask()).rename('harmonizado').toInt32()))
    return (ee.Image(saida.copyProperties(salvo, salvo.propertyNames()))
              .set({'harmonizacao': f'{MODIS} por cena; ganho+offset VNIR/B04, '
                                    'ganho único B05–B09; TIR sem correção',
                    'SCALE_FACTORS': 'reflectancia de superficie (NBAR) x1000; TIR em DN',
                    'harmonizado': '1 = reta da cena, 2 = mediana das retas do período',
                    'cenas_validas': int(df['valida'].sum()), 'cenas_total': len(df)}))


# ═══════════════════════════════════════════════════════════════════════════
# UTILITÁRIOS
# ═══════════════════════════════════════════════════════════════════════════
def existe(asset: str) -> bool:
    try:
        ee.data.getAsset(asset)
        return True
    except ee.EEException:
        return False


def tasks_ativas() -> set[str]:
    return {t['description'] for t in ee.data.getTaskList() if t['state'] in ('READY', 'RUNNING')}


def ler_coef(periodo: str) -> pd.DataFrame:
    arq = arq_coef(periodo)
    if not arq.exists():
        raise SystemExit(f'{arq} não existe — rode as etapas 2 e 3 antes.')
    return pd.read_csv(arq).drop_duplicates('cena_t').sort_values('cena_t')


def exportar_imagem(img, nome: str, asset: str, geom, piramide: dict) -> None:
    task = ee.batch.Export.image.toAsset(
        image=img, description=nome, assetId=asset, region=geom, crs=CRS_EXPORT,
        crsTransform=CRS_TRANSFORM, maxPixels=1e13, pyramidingPolicy=piramide)
    task.start()
    log.info(f'  task {task.id} → {asset}')


def itens(etapa: int, periodos, celulas, grid, colecao: str) -> list[tuple[str, int, str]]:
    """(período, célula, nome) a exportar: tem original, não existe e não está rodando."""
    originais = listar_destino(s.DESTINO_FINAL, criar=False)
    feitos = listar_destino(colecao, criar=True)
    mapas = listar_destino(COL_CENAS, criar=False) if etapa == 4 else {}
    ativas = tasks_ativas()
    fila = []
    for p in periodos:
        for cid in celulas:
            nome = s.nome_mosaico(p, grid[cid]['nome'])
            if nome not in originais:
                motivo = 'sem mosaico original'
            elif nome in feitos:
                motivo = 'já existe'
            elif f'{"cenas" if etapa == 1 else "harm"}_{p}_{grid[cid]["nome"]}' in ativas:
                motivo = 'processando'
            elif etapa == 4 and nome not in mapas:
                motivo = 'sem mapa de cenas (etapa 1)'
            else:
                fila.append((p, cid, nome))
                continue
            log.info(f'  {p} {grid[cid]["nome"]}: pulado ({motivo})')
    return fila


# ═══════════════════════════════════════════════════════════════════════════
# ETAPAS
# ═══════════════════════════════════════════════════════════════════════════
def etapa1(periodos, celulas, grid):
    fila = itens(1, periodos, celulas, grid, COL_CENAS)
    if not fila or not perguntar_sim_nao(f'Exportar {len(fila)} mapas de cenas?', padrao=True):
        return
    cache = {p: listas(p) for p in periodos}
    for p, cid, nome in fila:
        salvo = ee.Image(f'{s.DESTINO_FINAL}/{nome}')
        img = mapa_cenas(p, grid[cid]['geom'], salvo, *cache[p]).set(
            {'celula_id': cid, 'celula_nome': grid[cid]['nome']})
        exportar_imagem(img, f'cenas_{p}_{grid[cid]["nome"]}', f'{COL_CENAS}/{nome}',
                        grid[cid]['geom'], {'cena_t': 'mode', 'confere': 'min'})


def com_tentativas(func, tentativas: int = 6):
    """Repete em erro transitório (429, timeout), esperando 5, 10, 20… s."""
    for i in range(tentativas):
        try:
            return func()
        except ee.EEException as e:
            transitorio = any(m in str(e).lower() for m in (
                'too many', 'timeout', 'deadline', 'unavailable', 'internal', 'quota'))
            if not transitorio or i == tentativas - 1:
                raise
            time.sleep(5 * 2 ** i)


def meses(periodo: str) -> list[tuple[int, int]]:
    """(ano, mês) secos dentro das datas do período."""
    cfg = s.CONFIG[periodo]
    ini, fim = cfg['data_inicio'][:7], cfg['data_fim'][:7]
    return [(a, m) for a in range(int(ini[:4]), int(fim[:4]) + 1) for m in s.MESES_SECO
            if ini <= f'{a}-{m:02d}' < fim]


def tasks_por_estado() -> dict[str, str]:
    """Estado da última task de cada descrição (getTaskList vem da mais recente)."""
    estados = {}
    for tk in ee.data.getTaskList():
        estados.setdefault(tk.get('description', ''), tk['state'])
    return estados


def etapa2(periodos):
    """Pares ASTER × MODIS em BATCH: uma task Export.table.toDrive por período × mês
    (CSV em Drive/PASTA_DRIVE/<P>_<ano>_<mês>.csv). As tasks rodam em paralelo nos
    servidores do GEE, sem o limite de 5 min do getInfo. Rodar de novo reenvia só os
    meses sem task concluída nem em andamento."""
    estados = tasks_por_estado()
    for p in periodos:
        col, sortear = colecao_coef(p)
        novas, puladas = 0, 0
        for ano, mes in meses(p):
            nome = f'pares_modis_{p}_{ano}_{mes:02d}'
            if estados.get(nome) in ('READY', 'RUNNING', 'COMPLETED'):
                puladas += 1
                continue
            col_m = (col.filter(ee.Filter.calendarRange(ano, ano, 'year'))
                        .filter(ee.Filter.calendarRange(mes, mes, 'month')))
            task = ee.batch.Export.table.toDrive(
                collection=col_m.map(sortear).flatten(), description=nome,
                folder=PASTA_DRIVE, fileNamePrefix=f'{p}_{ano}_{mes:02d}',
                fileFormat='CSV', selectors=colunas_pontos(p))
            task.start()
            novas += 1
            log.info(f'  task {task.id} → Drive/{PASTA_DRIVE}/{p}_{ano}_{mes:02d}.csv')
        log.info(f'Período {p}: {novas} tasks novas, {puladas} meses já concluídos ou rodando.')
    print(f'\nQuando terminarem, baixe a pasta Drive/{PASTA_DRIVE} e ponha os CSVs em\n'
          f'  {PASTA_PONTOS}\ne rode a etapa 3.')


def etapa3_retas(periodos):
    """Lê os CSVs de pares baixados do Drive, ajusta a reta de cada cena e grava
    coef_modis_<P>.csv. Só local: pode rodar de novo à vontade (ex.: mudando os
    parâmetros de ajuste)."""
    for p in periodos:
        arqs = sorted(PASTA_PONTOS.glob(f'{p}_*.csv'))
        esperados = {f'{p}_{a}_{m:02d}' for a, m in meses(p)}
        faltam = sorted(esperados - {a.stem for a in arqs})
        if not arqs:
            log.warning(f'Período {p}: nenhum CSV {p}_*.csv em {PASTA_PONTOS}.')
            continue
        if faltam:
            log.warning(f'Período {p}: faltam {len(faltam)} meses ({", ".join(faltam[:6])}'
                        f'{"…" if len(faltam) > 6 else ""}) — as cenas deles usarão a mediana.')
        tabelas = [pd.read_csv(a) for a in arqs if a.stat().st_size > 0]
        pts = pd.concat([tb for tb in tabelas if len(tb)], ignore_index=True)   # mês sem cena = só cabeçalho
        linhas = []
        for ct, g in pts.groupby('cena_t'):
            meta = g.iloc[0][['cena_t', 'cena_id', 'data', 'CLOUDCOVER']].to_dict()
            linhas.append({**meta, **ajustar(g, p)})
        df = pd.DataFrame(linhas).reindex(columns=colunas(p)).sort_values('cena_t')
        df.to_csv(arq_coef(p), index=False)
        log.info(f'Período {p}: {len(df)} cenas com pontos, {len(pts)} pares → {arq_coef(p)}')


def etapa4(periodos, celulas, grid):
    fila = itens(4, periodos, celulas, grid, COL_HARM)
    if not fila:
        return
    tab = {p: retas(p) for p in {p for p, _, _ in fila}}
    for p, (df, reserva) in tab.items():
        log.info(f'Período {p}: {int(df["valida"].sum())}/{len(df)} cenas com reta válida; '
                 f'reserva (mediana): ' + ', '.join(f'{k} g={g:.3f} o={o:.4f}'
                                                   for k, (g, o) in reserva.items()))
    if not perguntar_sim_nao(f'Exportar {len(fila)} mosaicos harmonizados?', padrao=True):
        return
    for p, cid, nome in fila:
        salvo = ee.Image(f'{s.DESTINO_FINAL}/{nome}')
        cenas = ee.Image(f'{COL_CENAS}/{nome}')
        img = harmonizar(p, salvo, cenas, *tab[p])
        exportar_imagem(img, f'harm_{p}_{grid[cid]["nome"]}', f'{COL_HARM}/{nome}',
                        grid[cid]['geom'], {'.default': 'mean', 'cena_t': 'mode',
                                            'harmonizado': 'mode', 'n_obs': 'mode'})


def status(periodos, grid):
    cenas = listar_destino(COL_CENAS, criar=False)
    harm = listar_destino(COL_HARM, criar=False)
    originais = listar_destino(s.DESTINO_FINAL, criar=False)
    for p in periodos:
        nomes = {cid: s.nome_mosaico(p, grid[cid]['nome']) for cid in sorted(grid)}
        print(f'\nPeríodo {p}: originais {sum(n in originais for n in nomes.values())}, '
              f'mapas de cenas {sum(n in cenas for n in nomes.values())}, '
              f'harmonizados {sum(n in harm for n in nomes.values())} (de {len(grid)})')
        com_mapa = {cid: n for cid, n in nomes.items() if n in cenas}
        if com_mapa:
            # uma célula por vez: as 25 amostras num getInfo só dão
            # "Too many concurrent aggregations" (429)
            print('  confere (fração de pixels iguais ao original):')
            for cid, n in com_mapa.items():
                amostra = (ee.Image(f'{COL_CENAS}/{n}').select('confere')
                             .sample(region=grid[cid]['geom'], numPixels=2000, seed=1,
                                     dropNulls=True)
                             .aggregate_mean('confere'))
                try:
                    f = com_tentativas(amostra.getInfo)
                except ee.EEException as e:
                    print(f'    {grid[cid]["nome"]:6s} erro: {str(e)[:80]}')
                    continue
                aviso = '' if f is not None and f >= 0.99 else '   ← mosaico não bate: lista mudou?'
                print(f'    {grid[cid]["nome"]:6s} {f if f is None else round(f, 4)}{aviso}')
        if arq_coef(p).exists():
            df, reserva = retas(p)
            print(f'  retas: {int(df["valida"].sum())}/{len(df)} cenas válidas, '
                  f'mediana de {df["n_pontos"].median():.0f} pontos/cena ← {arq_coef(p)}')
            for k in CHAVES[p]:
                v = df[df['valida']]
                print(f'    {k}: ganho mediana {reserva[k][0]:.3f} '
                      f'[{v[f"g_{k}"].quantile(.1):.3f}–{v[f"g_{k}"].quantile(.9):.3f}], '
                      f'offset {reserva[k][1]:+.4f}, R² mediano {v[f"r2_{k}"].median():.3f}')
        else:
            print(f'  retas: {arq_coef(p)} ainda não existe')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    switch_user(CONTA)
    ee.Initialize(project=PROJETO)
    log.info(f'Earth Engine com {CONTA} (projeto {PROJETO})')
    grid = ler_grid()

    print('\n══════ Harmonização dos mosaicos sazonais secos com o MCD43A4 ══════')
    etapa = perguntar_opcao('Etapa:', [
        ('mapa de cenas (cena_t + confere) por célula', 1),
        ('pares ASTER × MODIS: tasks por mês → Drive', 2),
        ('retas por cena a partir dos CSVs baixados (local)', 3),
        ('aplicar as retas (precisa de 1 e 3)', 4),
        ('status', 5),
    ])
    periodos = perguntar_opcao('Períodos:', [('A e B', ['A', 'B']), ('só A', ['A']), ('só B', ['B'])])
    if etapa in (1, 4):
        celulas = perguntar_celulas(grid)
        (etapa1 if etapa == 1 else etapa4)(periodos, celulas, grid)
    elif etapa == 2:
        etapa2(periodos)
    elif etapa == 3:
        etapa3_retas(periodos)
    else:
        status(periodos, grid)


if __name__ == '__main__':
    main()
