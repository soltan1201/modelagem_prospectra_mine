"""
Mosaicos ASTER por célula da grid 5 × 5 (Irecê): dois mosaicos base separados
pela presença do SWIR, e o SWIR estimado para o mosaico que não tem.

╔══ PASSOS PARA TERMINAR (2026-10-18) — guia: docs/passos_mosaico_final.md ══╗
  python build_mosaico_unico.py → "1 — revisar nos assets o que falta e lançar
  tudo": mostra célula × com_swir / sem_swir / final (ok · fila · FALTA) lido nos
  ASSETS e lança o que falta no rodízio caatinga01–05 + solkan1201, 2 tasks por
  conta, projeto EE mapbiomas-brazil. Agora: 2 bases (sem_swir L3C4, L3C5) + 23
  finais; depois que as 2 bases terminarem, rodar de novo para os 2 finais.
╚════════════════════════════════════════════════════════════════════════════╝

ETAPA 1 — mosaicos base (DOIS assets por célula, cada um com quality)
  com_swir  só cenas COM SWIR (até 2008-S1: têm GAIN_COEFFICIENT_B04); o pixel
            só vale com VNIR + SWIR + TIR.
  sem_swir  só cenas SEM SWIR (a partir de 2008-S2); bandas VNIR + TIR.
  Os dois são montados com a mesma regra:
    camada 1  cenas ZERO-CLOUD (lista inspecionada), empilhadas pelo NDVI
              mediano da cena dentro da célula: a MAIS SECA (menor NDVI) fica
              por cima. Cena intacta (só TOA), exceto pixels saturados (VNIR
              DN ≥ 254 — nas cenas zero-cloud eram nuvem), que viram lacuna.
    camada 2  as DEMAIS cenas (CC < 77, diurnas, fora da blacklist e da lista
              zero-cloud), com máscara de nuvem/sombra e saturação, compostas
              por qualityMosaic com o quality v6 (vegetação penaliza). Entra SÓ
              nas lacunas que a camada 1 não cobriu.
  Meses: perguntado ao rodar (padrão jul–out, estação seca; ou o ano todo).
  Saída: DESTINO_BASE/ASTER_base_com_swir_<LxCy> e ASTER_base_sem_swir_<LxCy>, INT32:
    reflectância TOA × 1000 (B01–B09), TIR em DN (B10–B14), 'saturado',
    'quality' (× 1000; camada 1 = 5000–5500), 'data' (AAAAMMDD da cena do
    pixel), 'camada' (1 zero-cloud, 2 preenchimento).

ETAPA 2 — SWIR estimado no mosaico sem SWIR
  Lê ASTER_base_sem_swir_<LxCy> e estima B04–B09 com o Random Forest de
  modeling/estimar_swir_rf.py (mesmo treino ASSET_PONTOS, hiperparâmetros e
  normalização OLS por cena — a banda 'data' diz a cena de cada pixel; cena
  sem reta própria usa a reta global).
  Saída: DESTINO_FINAL/ASTER_sem_swir_estimado_<LxCy>, INT32: B01–B14 (B04–B09
    estimados), 'quality', 'data', 'camada', 'normalizacao' (1 global, 2 cena).

ETAPA 3 — mosaico único juntando com_swir e sem_swir (SWIR estimado de novo)
  Lê os dois bases da célula. O SWIR do sem_swir é estimado como no
  train_swir_models.ipynb: reta por banda VNIR/TIR levando o 2009–2020 à
  escala de 2000–2008, depois o RF (mesmo treino da etapa 2). Reta RMA
  (padrão; ganho σ_com/σ_sem) ou OLS; pares todos (padrão) ou só zero-cloud.
  Teste L2C2/L1C1: RMA leva os ganhos de 0,53–0,71 para 0,83–1,12, mas o viés
  do SWIR quase não muda (−9 % e +3/+11 %) — o resto vem do RF, não da reta. Retas:
    'mosaicos'  (padrão) pares dos DOIS mosaicos da célula, pixel a pixel
                (sem_swir → com_swir), só onde os dois existem e a vegetação
                é parecida (|ΔNDVI| ≤ LIMIAR_NDVI), ajustadas por cena do
                sem_swir (banda 'data'); cena com < MIN_PIX_CENA pares usa a
                reta da célula. Motivo: os pontos pareados só têm reta para
                24 cenas — 78 % dos pixels da etapa 2 caíram na reta global
                (ganho 0,24–0,58), que achata a imagem (SWIR 15–29 % abaixo
                do medido).
    'pareados'  as retas da etapa 2 (pontos pareados), para comparação.
  Pixel a pixel vence: zero-cloud (camada 1) sobre preenchimento; empate →
  menor NDVI (o mais seco). Se vence o com_swir, o SWIR é o medido.
  Saída: DESTINO_JUNTO/ASTER_mosaico_unico_<LxCy>, INT32: B01–B14, 'quality',
    'data', 'camada', 'swir_estimado' (0 medido, 1 estimado), 'normalizacao'
    (0 medido, 1 reta da célula/global, 2 reta da cena).

Contas: as de build_mosaicos_semestrais_v2 (todas que escrevem em
mapbiomas-arida/mine; rodízio a cada TASKS_POR_CONTA tasks). Não reenvia
assets já salvos nem com task ativa nas contas escolhidas.

Uso: python build_mosaico_unico.py
"""
import csv
import logging
import sys
from datetime import datetime
from pathlib import Path

import ee

import build_mosaicos_semestrais_v2 as b
import build_mosaicos_zero_cloud_v2 as z

sys.path.append(str(Path(__file__).resolve().parents[1] / 'modeling'))
import estimar_swir_rf as rf   # noqa: E402  (treino, params e normalização do RF)

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
MINE          = 'projects/mapbiomas-arida/mine'
DESTINO_BASE  = f'{MINE}/ASTER_mosaico_unico_base'
DESTINO_FINAL = f'{MINE}/ASTER_mosaico_unico'
PREFIXO_BASE  = 'ASTER_base_'                  # + com_swir_<LxCy> / sem_swir_<LxCy>
PREFIXO_FINAL = 'ASTER_sem_swir_estimado_'
DESTINO_JUNTO = f'{MINE}/ASTER_mosaico_unico_final'
PREFIXO_JUNTO = 'ASTER_mosaico_unico_'
VERSAO        = 'unico_v2'

# superconta: não lê nem escreve mapbiomas-arida → salva em mapbiomas-brazil (CAATINGA),
# coleções com all_users_can_read para as contas arida lerem. As bases são procuradas
# nos dois projetos (primeiro arida).
CAATINGA_BR = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA'

# Todas as contas inicializam o EE no projeto mapbiomas-brazil (cota de tasks do projeto,
# maior que a dos projetos pessoais); o destino dos assets não muda. Contas testadas em
# 2026-10-17 (inicializam no projeto e leem/escrevem mapbiomas-arida). Fora:
# solkanCengine (sem serviceusage no mapbiomas-brazil) e solkanGeodatin (não lê arida
# por esse projeto).
PROJETO_EE        = 'mapbiomas-brazil'
TASKS_POR_CONTA   = 2       # rodízio: troca de conta a cada N tasks
CONTAS_PROJETO_EE = ['caatinga01', 'caatinga02', 'caatinga03', 'caatinga04', 'caatinga05',
                     'solkan1201', 'superconta']
RAIZES      = {'arida': MINE, 'superconta': CAATINGA_BR}
COLECOES    = {1: 'ASTER_mosaico_unico_base', 2: 'ASTER_mosaico_unico', 3: 'ASTER_mosaico_unico_final'}

# etapa 3: retas sem_swir → com_swir nos pixels dos dois mosaicos
LIMIAR_NDVI   = 0.05        # |ΔNDVI| máximo para o par entrar na reta
MIN_PIX_CENA  = 200         # pares mínimos para a cena ter reta própria
ESCALA_AJUSTE = 90          # m (resolução do TIR)
METODO_RETA   = 'rma'       # 'ols' (y sobre x; ganho < 1 achata) ou 'rma' (ganho = σ_com / σ_sem)
PARES_RETA    = 'todos'     # 'todos' ou 'zero_cloud' (só pixels da camada 1 nos dois mosaicos)

PERIODO     = ('2000-01-01', '2021-01-01')
MESES_SECOS = b.MESES_SECOS          # (7, 10)

VNIR = ['B01', 'B02', 'B3N']
SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']
TIR  = b.TIR
EXTRAS = ['saturado', 'quality', 'data', 'camada']

# tipo de mosaico → bandas ópticas e se a cena precisa ter SWIR (GAIN_COEFFICIENT_B04)
TIPOS = {
    'com_swir': {'opticas': VNIR + SWIR + TIR, 'swir': True},
    'sem_swir': {'opticas': VNIR + TIR,        'swir': False},
}


def filtro_swir(com_swir: bool):
    # criado sob demanda: objetos ee não podem ser montados antes do ee.Initialize
    f = ee.Filter.notNull(['GAIN_COEFFICIENT_B04'])
    return f if com_swir else f.Not()


# ═══════════════════════════════════════════════════════════════════════════
# PREPARO DE UMA CENA
# ═══════════════════════════════════════════════════════════════════════════
def _data(img):
    d = ee.Number.parse(ee.Date(img.get('system:time_start')).format('YYYYMMdd'))
    return ee.Image.constant(d).toInt32().rename('data')


def preparar_zero_cloud(img, opticas: list[str]):
    """Camada 1: TOA intacta; só os pixels saturados viram lacuna. Máscara igual
    em todas as bandas (o mosaic() compõe banda a banda)."""
    img = ee.Image(img)
    t = ee.Image(b.radiance_to_toa(img, mascarar_saturacao=False))
    nir, red = t.select('B3N'), t.select('B02')
    savi = nir.subtract(red).multiply(1.5).divide(nir.add(red).add(0.5 * b.ESCALA_TOA).add(1e-6))
    quality = ee.Image(1).subtract(savi.unitScale(-0.1, 0.5).clamp(0, 1)) \
                .multiply(500).add(b.QUALITY_ZERO_CLOUD).rename('quality')
    valido = t.select(opticas).mask().reduce(ee.Reducer.min()).And(t.select('saturado').eq(0))
    return (ee.Image.cat([t.select(opticas), t.select('saturado'), quality, _data(img),
                          ee.Image.constant(1).rename('camada')])
              .toFloat().updateMask(valido))


def preparar_preenchimento(img, opticas: list[str]):
    """Camada 2: máscara de saturação + nuvem/sombra; quality v6."""
    img = ee.Image(img)
    t = ee.Image(b.radiance_to_toa(img))                     # mascara saturação
    q = ee.Image(b.add_quality_v5(b.mascara_nuvem(t)))       # v6 no módulo v2 (vegetação penaliza)
    valido = q.select(opticas + ['quality']).mask().reduce(ee.Reducer.min())
    return (ee.Image.cat([q.select(opticas), q.select('saturado'), q.select('quality'), _data(img),
                          ee.Image.constant(2).rename('camada')])
              .toFloat().updateMask(valido))


def nos_meses(col, meses):
    """Filtra a coleção pelos meses (None = ano todo)."""
    return col.filter(ee.Filter.calendarRange(meses[0], meses[1], 'month')) if meses else col


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 1 — MOSAICOS BASE
# ═══════════════════════════════════════════════════════════════════════════
def montar_base(tipo: str, celula: dict, ids_zero: list[str], ids_black: list[str], meses):
    geom, cfg = celula['geom'], TIPOS[tipo]
    opticas = cfg['opticas']
    zc = (nos_meses(ee.ImageCollection([ee.Image(f'ASTER/AST_L1T_003/{i}') for i in ids_zero])
                      .filterDate(*PERIODO), meses)
            .filterBounds(geom).filter(filtro_swir(cfg['swir']))
            .map(lambda i: z.ndvi_cena(i, geom))
            .sort('ndvi_cena', False))                       # última (mais seca) por cima
    demais = (nos_meses(ee.ImageCollection('ASTER/AST_L1T_003').filterDate(*PERIODO), meses)
                .filterBounds(geom).filter(filtro_swir(cfg['swir']))
                .filter(ee.Filter.lt('CLOUDCOVER', b.CC_MAX))
                .filter(ee.Filter.And(ee.Filter.gt('SOLAR_ELEVATION', 0),
                                      ee.Filter.notNull(['GAIN_COEFFICIENT_B01'])))
                .filter(ee.Filter.inList('system:index', ee.List(ids_black + ids_zero)).Not()))

    camada1 = zc.map(lambda i: preparar_zero_cloud(i, opticas)).mosaic()
    camada2 = demais.map(lambda i: preparar_preenchimento(i, opticas)).qualityMosaic('quality')
    base = ee.ImageCollection([camada2, camada1]).mosaic()   # camada 1 por cima; 2 só nas lacunas

    bandas = opticas + EXTRAS
    saida = ee.Image.cat([base.select(x).round().toInt32().rename(x) for x in bandas]).clip(geom)
    props = {
        'tipo': tipo, 'celula_id': celula['id'], 'celula_nome': celula['nome'],
        'grid_asset': b.GRID_ASSET, 'etapa': 'base', 'versao': VERSAO,
        'meses': f'{meses[0]}-{meses[1]}' if meses else '1-12',
        'periodo': f'{PERIODO[0][:4]}-{int(PERIODO[1][:4]) - 1}',
        'regra': 'camada 1 zero-cloud (menor NDVI da cena por cima) + camada 2 qualityMosaic v6 nas lacunas',
        'n_zero_cloud': zc.size(), 'n_preenchimento': demais.size(),
        'cenas_zero_cloud': zc.aggregate_array('system:index').join(','),
        'ndvi_cenas_zero_cloud': zc.aggregate_array('ndvi_cena'),
        'SCALE_FACTORS': {**{x: (b.ESCALA_TOA if x in VNIR + SWIR else 1) for x in opticas},
                          'quality': 1000, 'saturado': 1, 'data': 1, 'camada': 1},
        'DATA_TYPE': 'INT32',
    }
    return ee.Image(saida.set(props)), {'zero_cloud': zc.size(), 'preenchimento': demais.size()}


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 2 — SWIR ESTIMADO NO MOSAICO SEM SWIR
# ═══════════════════════════════════════════════════════════════════════════
def montar_estimado(celula: dict, bases: dict[str, str], modelos: dict, coefs: dict):
    origem = bases[f'{PREFIXO_BASE}sem_swir_{celula["nome"]}']
    base = ee.Image(origem)
    est = rf.estimar(base, modelos, coefs, 'cena')            # VNIR medido, SWIR estimado, TIR medido
    saida = ee.Image.cat([est.select(VNIR + SWIR + TIR),
                          base.select(['quality', 'data', 'camada']),
                          est.select('normalizacao').toInt32()]).updateMask(est.select('B01').mask())
    props = {
        'tipo': 'sem_swir_estimado', 'celula_id': celula['id'], 'celula_nome': celula['nome'],
        'grid_asset': b.GRID_ASSET, 'etapa': 'swir_estimado', 'versao': VERSAO, 'base': origem,
        'swir': 'estimado_rf', 'bandas_estimadas': ','.join(SWIR),
        'pontos_treino': rf.ASSET_PONTOS, 'normalizacao': 'cena (reta OLS por data; senao global)',
        'normalizacao_coefs': rf.ARQ_COEFS.name,
        'SCALE_FACTORS': {**dict.fromkeys(VNIR + SWIR, b.ESCALA_TOA), **dict.fromkeys(TIR, 1),
                          'quality': 1000, 'data': 1, 'camada': 1, 'normalizacao': 1},
        'DATA_TYPE': 'INT32',
    }
    return ee.Image(saida.set(props))


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 3 — MOSAICO ÚNICO (com_swir + sem_swir com SWIR estimado)
# ═══════════════════════════════════════════════════════════════════════════
def _ndvi(img):
    return img.normalizedDifference(['B3N', 'B02'])


def listar_bases() -> dict[str, str]:
    """{nome: caminho} dos mosaicos base nos dois projetos (arida primeiro). Projeto sem
    leitura para a conta ativa (superconta × arida) é ignorado."""
    achados = {}
    for raiz in (MINE, CAATINGA_BR):
        col = f'{raiz}/{COLECOES[1]}'
        for nome in b.listar_destino(col, criar=False):
            achados.setdefault(nome, f'{col}/{nome}')
    return achados


def bases_celula(nome: str, bases: dict[str, str]):
    com = ee.Image(bases[f'{PREFIXO_BASE}com_swir_{nome}'])
    sem = ee.Image(bases[f'{PREFIXO_BASE}sem_swir_{nome}'])
    return com, sem.updateMask(sem.select('saturado').eq(0))


def _red_rma():
    """Média e desvio dos dois lados do par (x = sem_swir, y = com_swir) — base da reta RMA.
    Montado sob demanda: objetos ee não podem ser criados antes do ee.Initialize."""
    return ee.Reducer.mean().repeat(2).combine(ee.Reducer.stdDev().repeat(2), sharedInputs=True)


def _reta(r: dict, metodo: str):
    """(offset, ganho) a partir da saída do redutor (linearFit ou mean/stdDev repetidos)."""
    if metodo == 'ols':
        return (r['offset'], r['scale']) if r.get('scale') is not None else None
    if not r.get('mean') or None in r['mean'] or None in r['stdDev'] or not r['stdDev'][0]:
        return None
    (mx, my), (sx, sy) = r['mean'], r['stdDev']
    ganho = sy / sx
    return (my - ganho * mx, ganho)


def retas_mosaicos(celula: dict, bases: dict[str, str],
                   metodo: str = METODO_RETA, pares: str = PARES_RETA) -> dict:
    """Retas por banda VNIR/TIR, sem_swir → com_swir (com = offset + ganho × sem), nos
    pixels onde os dois existem e |ΔNDVI| ≤ LIMIAR_NDVI (pares 'zero_cloud': só camada 1
    nos dois). Por cena do sem_swir (banda 'data', ≥ MIN_PIX_CENA pares) e da célula.
    metodo 'ols' = mínimos quadrados de com sobre sem (como o notebook); 'rma' = eixo
    maior reduzido, ganho σ_com/σ_sem — não encolhe o ganho pelo ruído dos dois lados.
    Mesmo formato de rf.coeficientes() (reta da célula em 'global') para usar
    rf.normalizar(..., 'cena').

    Calculadas aqui (getInfo, uma redução por vez) e não dentro da imagem exportada:
    reduções aninhadas estouram "Too many concurrent aggregations".
    """
    com, sem = bases_celula(celula['nome'], bases)
    estavel = _ndvi(sem).subtract(_ndvi(com)).abs().lte(LIMIAR_NDVI)
    if pares == 'zero_cloud':
        estavel = estavel.And(sem.select('camada').eq(1)).And(com.select('camada').eq(1))
    data = sem.select('data').toInt32()
    kw = dict(geometry=celula['geom'], scale=ESCALA_AJUSTE, maxPixels=1e10, tileScale=4)
    contagem = (ee.Image.cat([estavel.rename('n'), data]).updateMask(estavel)
                  .reduceRegion(ee.Reducer.count().group(1, 'data'), **kw).get('groups').getInfo())
    n_cena = {int(g['data']): g['count'] for g in contagem}
    coefs = {'global': {}, 'cenas': {}, 'n': {'celula': sum(n_cena.values()), **n_cena}}
    red = ee.Reducer.linearFit() if metodo == 'ols' else _red_rma()
    for x in rf.FEATURES:
        par = ee.Image.cat([sem.select(x).toFloat(), com.select(x).toFloat()]).updateMask(estavel)
        cel = _reta(par.reduceRegion(red, **kw).getInfo(), metodo)
        if cel is None:
            raise ValueError(f'{celula["nome"]}: sem pares estáveis com_swir × sem_swir ({pares})')
        coefs['global'][x] = cel
        grupos = par.addBands(data).reduceRegion(red.group(2, 'data'), **kw).get('groups').getInfo()
        for g in grupos:
            d, reta = int(g['data']), _reta(g, metodo)
            if reta and n_cena.get(d, 0) >= MIN_PIX_CENA:
                coefs['cenas'].setdefault(d, {})[x] = reta
    coefs['cenas'] = {d: v for d, v in coefs['cenas'].items() if len(v) == len(rf.FEATURES)}
    return coefs


def gravar_retas(retas_cel: dict[str, dict]) -> None:
    linhas = [{'celula': nome, 'data': d, 'n_pares': c['n'].get(d, 0),
               **{f'{x}_{p}': v[x][i] for x in rf.FEATURES for i, p in enumerate(('offset', 'ganho'))}}
              for nome, c in retas_cel.items()
              for d, v in [('celula', c['global']), *sorted(c['cenas'].items())]]
    if not linhas:
        return
    arq = b.PASTA / f'retas_mosaico_unico_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    log.info(f'Retas → {arq}')


def montar_junto(celula: dict, bases: dict[str, str], modelos: dict, coefs: dict, descricao_retas: str):
    nome = celula['nome']
    com, sem = bases_celula(nome, bases)

    # SWIR estimado do sem_swir: features normalizadas → RF; VNIR/TIR da saída são os medidos
    feats, flag = rf.normalizar(sem, coefs, 'cena')
    mascara = sem.select(rf.FEATURES).mask().reduce(ee.Reducer.min())
    feats = feats.updateMask(mascara)
    swir = [feats.classify(modelos[s], s).rename(s) for s in SWIR]
    est = ee.Image.cat([sem.select(VNIR), *swir, sem.select(TIR),
                        sem.select(['quality', 'data', 'camada']),
                        ee.Image.constant(1).rename('swir_estimado'), flag]).updateMask(mascara)
    med = ee.Image.cat([com.select(VNIR + SWIR + TIR + ['quality', 'data', 'camada']),
                        ee.Image.constant(0).rename('swir_estimado'),
                        ee.Image.constant(0).rename('normalizacao')])

    # com_swir vence se: só ele existe; ou camada menor (zero-cloud); ou mesma camada e NDVI ≤
    tem_com = com.select('B01').mask().gt(0)
    tem_sem = est.select('B01').mask().gt(0)
    c_com, c_sem = com.select('camada').unmask(9), est.select('camada').unmask(9)
    usar_com = tem_com.And(tem_sem.Not()
                           .Or(c_com.lt(c_sem))
                           .Or(c_com.eq(c_sem).And(_ndvi(com).lte(_ndvi(sem)).unmask(0))))
    bandas = VNIR + SWIR + TIR + ['quality', 'data', 'camada', 'swir_estimado', 'normalizacao']
    junto = ee.ImageCollection([est.select(bandas).round().toInt32(),
                                med.select(bandas).round().toInt32().updateMask(usar_com)]).mosaic()
    props = {
        'tipo': 'mosaico_unico', 'celula_id': celula['id'], 'celula_nome': nome,
        'grid_asset': b.GRID_ASSET, 'etapa': 'junto', 'versao': VERSAO,
        'base_com_swir': com.get('system:id'), 'base_sem_swir': sem.get('system:id'),
        'meses_com_swir': com.get('meses'), 'meses_sem_swir': sem.get('meses'),
        'regra': 'zero-cloud (camada 1) vence; empate -> menor NDVI; com_swir = SWIR medido',
        'swir': 'medido (com_swir) / estimado_rf (sem_swir)', 'pontos_treino': rf.ASSET_PONTOS,
        'normalizacao': descricao_retas,
        'n_cenas_reta_propria': len(coefs['cenas']),
        'SCALE_FACTORS': {**dict.fromkeys(VNIR + SWIR, b.ESCALA_TOA), **dict.fromkeys(TIR, 1),
                          'quality': 1000, 'data': 1, 'camada': 1, 'swir_estimado': 1, 'normalizacao': 1},
        'DATA_TYPE': 'INT32',
    }
    return ee.Image(junto.clip(celula['geom']).set(props))


# ═══════════════════════════════════════════════════════════════════════════
# CONTAS, ESTADO, REGISTRO
# ═══════════════════════════════════════════════════════════════════════════
def usar_conta(conta: str) -> None:
    """Troca as credenciais para `conta` e inicializa o EE no PROJETO_EE (não no da conta)."""
    b.switch_user(conta)
    ee.Initialize(project=PROJETO_EE)
    log.info(f'Earth Engine inicializado com {conta} (projeto {PROJETO_EE}).')


def estado(contas: list[str], destino: str, prefixo: str, outro: str) -> tuple[set, set]:
    """Assets já salvos no destino (ou na mesma coleção do outro projeto) e tasks ativas com o
    prefixo nas contas escolhidas. Coleção nova no mapbiomas-brazil recebe all_users_can_read
    (as contas arida só leem o que é público lá)."""
    ativas = set()
    for conta in contas:
        usar_conta(conta)
        ativas |= {t['description'] for t in ee.data.getTaskList()
                   if t['state'] in ('READY', 'RUNNING') and t.get('description', '').startswith(prefixo)}
    novo = not b.listar_destino(destino, criar=False)
    salvos = set(b.listar_destino(destino, criar=True)) | set(b.listar_destino(outro, criar=False))
    if novo and destino.startswith(CAATINGA_BR):
        try:
            ee.data.setAssetAcl(destino, {'all_users_can_read': True})
            log.info(f'{destino}: leitura pública (all_users_can_read) ligada.')
        except ee.EEException as e:
            log.warning(f'{destino}: não foi possível ligar a leitura pública ({e}).')
    if b.ARQ_CONTA_ATUAL.read_text().strip() != contas[0]:
        usar_conta(contas[0])
    return salvos, ativas


def gravar_registro(linhas: list[dict], etapa: str) -> None:
    if not linhas:
        return
    arq = b.PASTA / f'registro_mosaico_unico_{etapa}_{datetime.now():%Y%m%d_%H%M}.csv'
    campos = ['asset', 'celula_nome', 'celula_id', 'tipo', 'zero_cloud', 'preenchimento',
              'status', 'conta', 'task_id']
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        w.writeheader()
        w.writerows(linhas)
    log.info(f'Registro → {arq}')


def exportar(imagem, nome: str, destino: str, celula: dict, piramide: dict):
    task = ee.batch.Export.image.toAsset(
        image=imagem, description=nome, assetId=f'{destino}/{nome}', region=celula['geom'],
        crs=b.CRS_EXPORT, crsTransform=b.CRS_TRANSFORM, maxPixels=1e13, pyramidingPolicy=piramide)
    task.start()
    return task.id


# ═══════════════════════════════════════════════════════════════════════════
# REVISÃO PELOS ASSETS — o que falta salvar, e lança tudo no rodízio
# ═══════════════════════════════════════════════════════════════════════════
def _meses_do_tipo(tipo: str, bases: dict[str, str]):
    """Meses usados nas bases já salvas do tipo (para a que falta sair igual às outras)."""
    salvo = next((p for n, p in sorted(bases.items()) if n.startswith(f'{PREFIXO_BASE}{tipo}_')), None)
    if salvo is None:
        return MESES_SECOS if tipo == 'com_swir' else None
    meses = ee.Image(salvo).get('meses').getInfo() or '1-12'
    ini, fim = map(int, meses.split('-'))
    return None if (ini, fim) == (1, 12) else (ini, fim)


def revisar_e_lancar() -> None:
    """Confere nos ASSETS (não nos CSV de registro) o que falta por célula — base com_swir,
    base sem_swir e mosaico final (etapa 3) — e exporta o que falta no rodízio das contas
    arida (projeto EE PROJETO_EE, TASKS_POR_CONTA tasks por conta). Etapa 3 só para células
    com as duas bases salvas; as que dependem de base ainda na fila ficam para a próxima
    rodada desta opção."""
    contas = [c for c in CONTAS_PROJETO_EE if c != 'superconta']
    print(f'Rodízio: {", ".join(contas)} — {TASKS_POR_CONTA} tasks por conta, projeto EE {PROJETO_EE}.')
    usar_conta(contas[0])
    grid = b.ler_grid()
    destino_base, destino_final = f'{MINE}/{COLECOES[1]}', f'{MINE}/{COLECOES[3]}'
    bases = listar_bases()
    finais = (set(b.listar_destino(destino_final, criar=True))
              | set(b.listar_destino(f'{CAATINGA_BR}/{COLECOES[3]}', criar=False)))
    # no projeto PROJETO_EE a lista de tasks é a do projeto (mesma para todas as contas)
    fila = {t['description'] for t in ee.data.getTaskList() if t['state'] in ('READY', 'RUNNING')}

    def situacao(nome: str, salvos) -> str:
        return 'ok' if nome in salvos else ('fila' if nome in fila else 'FALTA')

    lancar_base, lancar_final, depois = [], [], []
    print(f'\n{"célula":7} {"com_swir":9} {"sem_swir":9} {"final":9}')
    for c in sorted(grid):
        nome = grid[c]['nome']
        s = {t: situacao(f'{PREFIXO_BASE}{t}_{nome}', bases) for t in ('com_swir', 'sem_swir')}
        s['final'] = situacao(f'{PREFIXO_JUNTO}{nome}', finais)
        print(f'{nome:7} {s["com_swir"]:9} {s["sem_swir"]:9} {s["final"]:9}')
        lancar_base += [(t, c) for t in ('com_swir', 'sem_swir') if s[t] == 'FALTA']
        if s['final'] == 'FALTA':
            (lancar_final if s['com_swir'] == s['sem_swir'] == 'ok' else depois).append(c)
    print(f'\nA lançar: {len(lancar_base)} bases, {len(lancar_final)} mosaicos finais.')
    if depois:
        print('Final depois que a base terminar (rodar esta opção de novo): '
              + ', '.join(grid[c]['nome'] for c in depois))
    if not lancar_base and not lancar_final:
        print('Nada falta salvar além do que já está na fila.')
        return
    if b.perguntar_opcao('Modo:', [('só mostrar (não exporta)', True), ('exportar o que falta', False)]):
        return

    registro, n_disparadas, ativa = [], 0, contas[0]

    def disparar(imagem, nome: str, destino: str, celula: dict, piramide: dict, linha: dict):
        nonlocal n_disparadas, ativa
        conta = contas[(n_disparadas // TASKS_POR_CONTA) % len(contas)]
        if conta != ativa:
            usar_conta(conta)
            ativa = conta
        try:
            tid = exportar(imagem, nome, destino, celula, piramide)
        except ee.EEException as e:
            registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
            return
        n_disparadas += 1
        log.info(f'  [{conta}] task {tid} → {destino}/{nome}')
        registro.append({**linha, 'status': 'exportando', 'conta': conta, 'task_id': tid})

    if lancar_base:
        ids_zero = z.ler_cenas()
        ids_black = sorted(b.ler_ids('lista_imagens_Notaceites.txt'))
        meses = {t: _meses_do_tipo(t, bases) for t in {t for t, _ in lancar_base}}
        for t, c in lancar_base:
            celula = {'id': c, **grid[c]}
            img, n = montar_base(t, celula, ids_zero, ids_black, meses[t])
            n = ee.Dictionary(n).getInfo()
            nome = f'{PREFIXO_BASE}{t}_{grid[c]["nome"]}'
            linha = {'asset': nome, 'celula_nome': grid[c]['nome'], 'celula_id': c, 'tipo': t, **n}
            if n['zero_cloud'] + n['preenchimento'] == 0:
                registro.append({**linha, 'status': 'sem_cenas'})
                continue
            disparar(img, nome, destino_base, celula, {'.default': 'sample'}, linha)

    if lancar_final:
        modelos = rf.treinar(ee.FeatureCollection(rf.ASSET_PONTOS), rf.params_por_banda())
        descricao = (f'mosaicos: {METODO_RETA.upper()} sem->com por cena (>= {MIN_PIX_CENA} pares '
                     f'{PARES_RETA}, |dNDVI| <= {LIMIAR_NDVI}, {ESCALA_AJUSTE} m), senao reta da celula')
        retas_cel = {}
        for c in lancar_final:
            celula = {'id': c, **grid[c]}
            nome = f'{PREFIXO_JUNTO}{grid[c]["nome"]}'
            linha = {'asset': nome, 'celula_nome': grid[c]['nome'], 'celula_id': c, 'tipo': 'mosaico_unico'}
            try:
                coefs = retas_mosaicos(celula, bases)
            except (ee.EEException, ValueError) as e:
                registro.append({**linha, 'status': f'sem_retas: {str(e)[:60]}'})
                continue
            retas_cel[grid[c]['nome']] = coefs
            disparar(montar_junto(celula, bases, modelos, coefs, descricao), nome, destino_final, celula,
                     {'.default': 'mean', 'data': 'mode', 'camada': 'mode',
                      'swir_estimado': 'mode', 'normalizacao': 'mode'}, linha)
        gravar_retas(retas_cel)
    print(f'{n_disparadas} tasks disparadas.')
    gravar_registro(registro, 'revisao')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    print('══════ Mosaicos ASTER com / sem SWIR por célula ══════')
    if b.perguntar_opcao('O que fazer:', [
        (('revisar nos assets o que falta e lançar tudo (rodízio arida, '
          f'{TASKS_POR_CONTA} tasks por conta) — recomendado'), True),
        ('escolher etapa, contas e células à mão', False),
    ]):
        revisar_e_lancar()
        return
    etapa = b.perguntar_opcao('Etapa:', [
        (f'1 — mosaicos base com_swir e sem_swir → {COLECOES[1]}', 1),
        (f'2 — SWIR estimado no mosaico sem_swir → {COLECOES[2]} (substituída pela 3)', 2),
        (f'3 — mosaico único com_swir + sem_swir (SWIR estimado) → {COLECOES[3]}', 3),
    ])
    projeto = b.perguntar_opcao('Contas:', [
        (f'rodízio das contas arida (projeto EE {PROJETO_EE}) → {MINE}', 'arida'),
        (f'superconta (só ela; não lê/escreve arida) → {CAATINGA_BR}', 'superconta'),
    ])
    contas = ['superconta'] if projeto == 'superconta' else b.perguntar_contas()
    fora = [c for c in contas if c not in CONTAS_PROJETO_EE]
    if fora:
        print(f'Fora do rodízio (não usam o projeto {PROJETO_EE}): {", ".join(fora)}')
        contas = [c for c in contas if c in CONTAS_PROJETO_EE]
    if not contas:
        raise SystemExit(f'Nenhuma conta escolhida usa o projeto {PROJETO_EE}.')
    if projeto == 'superconta' and etapa != 1:
        print('ATENÇÃO: a superconta não lê as bases em mapbiomas-arida; células com base lá '
              'saem como sem_base. Para a etapa 3 use o rodízio arida (lê as duas).')
    destino = f'{RAIZES[projeto]}/{COLECOES[etapa]}'
    outro = f'{RAIZES["arida" if projeto == "superconta" else "superconta"]}/{COLECOES[etapa]}'
    prefixo = {1: PREFIXO_BASE, 2: PREFIXO_FINAL, 3: PREFIXO_JUNTO}[etapa]
    salvos, ativas = estado(contas, destino, prefixo, outro)
    grid = b.ler_grid()

    tipos, meses, retas, metodo, pares = [], None, None, METODO_RETA, PARES_RETA
    if etapa == 3:
        retas = b.perguntar_opcao('Retas de normalização 2009–2020 → 2000–2008:', [
            ('pares dos dois mosaicos da célula, por cena (recomendado)', 'mosaicos'),
            ('pontos pareados (as mesmas da etapa 2)', 'pareados'),
        ])
        if retas == 'mosaicos':
            metodo = b.perguntar_opcao('Reta:', [
                ('RMA — ganho σ_com/σ_sem, não achata (recomendado)', 'rma'),
                ('OLS — como o notebook (ganho < 1)', 'ols'),
            ])
            pares = b.perguntar_opcao('Pares:', [
                (f'todos com |ΔNDVI| ≤ {LIMIAR_NDVI} (recomendado)', 'todos'),
                ('só zero-cloud nos dois (poucos pares; célula sem → sem_retas)', 'zero_cloud'),
            ])
    if etapa == 1:
        tipos = b.perguntar_opcao('Mosaicos:', [
            ('os dois (com_swir e sem_swir)', ['com_swir', 'sem_swir']),
            ('só com_swir', ['com_swir']),
            ('só sem_swir', ['sem_swir']),
        ])
        meses = b.perguntar_opcao('Meses das cenas:', [
            (f'estação seca ({MESES_SECOS[0]}–{MESES_SECOS[1]})', MESES_SECOS),
            ('o ano todo (a ordem por NDVI decide)', None),
        ])
    celulas = b.perguntar_celulas(grid)
    dry_run = b.perguntar_opcao('Modo:', [
        ('só conferir (dry-run, não exporta)', True),
        ('exportar', False),
    ])

    if etapa == 1:
        itens = {(t, c): f'{PREFIXO_BASE}{t}_{grid[c]["nome"]}' for t in tipos for c in celulas}
    elif etapa == 2:
        itens = {('sem_swir_estimado', c): f'{PREFIXO_FINAL}{grid[c]["nome"]}' for c in celulas}
    else:
        itens = {('mosaico_unico', c): f'{PREFIXO_JUNTO}{grid[c]["nome"]}' for c in celulas}
    pendentes = [k for k, nome in itens.items() if dry_run or (nome not in salvos and nome not in ativas)]
    print(f'{len(itens)} assets pedidos; {len(itens) - len(pendentes)} já salvos ou em processamento; '
          f'{len(pendentes)} a {"conferir" if dry_run else "exportar"}.')
    if not pendentes or not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Nada a fazer.')

    registro, disparadas, conta_ativa = [], 0, contas[0]

    def proxima_conta():
        nonlocal conta_ativa
        conta = contas[(disparadas // TASKS_POR_CONTA) % len(contas)]
        if conta != conta_ativa:
            usar_conta(conta)
            conta_ativa = conta
        return conta

    if etapa == 1:
        ids_zero = z.ler_cenas()
        ids_black = sorted(b.ler_ids('lista_imagens_Notaceites.txt'))
        montados = {(t, c): montar_base(t, {'id': c, **grid[c]}, ids_zero, ids_black, meses)
                    for t, c in pendentes}
        contagens = ee.List([ee.Dictionary(n) for _, n in montados.values()]).getInfo()
        for ((t, c), (img, _)), n in zip(montados.items(), contagens):
            nome = itens[(t, c)]
            linha = {'asset': nome, 'celula_nome': grid[c]['nome'], 'celula_id': c, 'tipo': t, **n}
            log.info(f'{nome}: {n["zero_cloud"]} cenas zero-cloud, {n["preenchimento"]} de preenchimento')
            if n['zero_cloud'] + n['preenchimento'] == 0:
                registro.append({**linha, 'status': 'sem_cenas'})
                continue
            if dry_run:
                registro.append({**linha, 'status': 'dry_run'})
                continue
            conta = proxima_conta()
            try:
                tid = exportar(img, nome, destino, {'id': c, **grid[c]}, {'.default': 'sample'})
            except ee.EEException as e:
                registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
                continue
            disparadas += 1
            log.info(f'  [{conta}] task {tid} → {destino}/{nome}')
            registro.append({**linha, 'status': 'exportando', 'conta': conta, 'task_id': tid})
    elif etapa == 3:
        bases = listar_bases()
        coefs_pareados = rf.coeficientes() if retas == 'pareados' else None
        descricao = (f'mosaicos: {metodo.upper()} sem->com por cena (>= {MIN_PIX_CENA} pares {pares}, '
                     f'|dNDVI| <= {LIMIAR_NDVI}, {ESCALA_AJUSTE} m), senao reta da celula'
                     if retas == 'mosaicos' else 'pareados: retas da etapa 2 (cena, senao global)')
        modelos = rf.treinar(ee.FeatureCollection(rf.ASSET_PONTOS), rf.params_por_banda())
        retas_cel = {}
        for t, c in pendentes:
            nome = itens[(t, c)]
            linha = {'asset': nome, 'celula_nome': grid[c]['nome'], 'celula_id': c, 'tipo': t}
            faltam = [x for x in ('com_swir', 'sem_swir')
                      if f'{PREFIXO_BASE}{x}_{grid[c]["nome"]}' not in bases]
            if faltam:
                log.info(f'{nome}: falta o mosaico base {" e ".join(faltam)} (rode a etapa 1) — pulado')
                registro.append({**linha, 'status': 'sem_base'})
                continue
            if retas == 'mosaicos':
                try:
                    coefs = retas_mosaicos({'id': c, **grid[c]}, bases, metodo, pares)
                except (ee.EEException, ValueError) as e:
                    log.info(f'{nome}: retas não calculadas ({str(e)[:80]}) — pulado')
                    registro.append({**linha, 'status': f'sem_retas: {str(e)[:60]}'})
                    continue
                retas_cel[grid[c]['nome']] = coefs
                log.info(f'{nome}: {coefs["n"]["celula"]} pares estáveis; '
                         f'{len(coefs["cenas"])} cenas com reta própria; ganho da célula '
                         + ', '.join(f'{x} {coefs["global"][x][1]:.2f}' for x in rf.FEATURES))
            else:
                coefs = coefs_pareados
            if dry_run:
                registro.append({**linha, 'status': 'dry_run'})
                continue
            conta = proxima_conta()
            try:
                tid = exportar(montar_junto({'id': c, **grid[c]}, bases, modelos, coefs, descricao), nome, destino,
                               {'id': c, **grid[c]},
                               {'.default': 'mean', 'data': 'mode', 'camada': 'mode',
                                'swir_estimado': 'mode', 'normalizacao': 'mode'})
            except ee.EEException as e:
                registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
                continue
            disparadas += 1
            log.info(f'  [{conta}] task {tid} → {destino}/{nome}')
            registro.append({**linha, 'status': 'exportando', 'conta': conta, 'task_id': tid})
    else:
        bases = listar_bases()
        params = rf.params_por_banda()
        coefs = rf.coeficientes()
        modelos = rf.treinar(ee.FeatureCollection(rf.ASSET_PONTOS), params)
        for t, c in pendentes:
            nome = itens[(t, c)]
            linha = {'asset': nome, 'celula_nome': grid[c]['nome'], 'celula_id': c, 'tipo': t}
            if f'{PREFIXO_BASE}sem_swir_{grid[c]["nome"]}' not in bases:
                log.info(f'{nome}: sem mosaico base sem_swir (rode a etapa 1) — pulado')
                registro.append({**linha, 'status': 'sem_base'})
                continue
            if dry_run:
                registro.append({**linha, 'status': 'dry_run'})
                continue
            conta = proxima_conta()
            try:
                tid = exportar(montar_estimado({'id': c, **grid[c]}, bases, modelos, coefs), nome, destino,
                               {'id': c, **grid[c]},
                               {'.default': 'mean', 'data': 'mode', 'camada': 'mode', 'normalizacao': 'mode'})
            except ee.EEException as e:
                registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
                continue
            disparadas += 1
            log.info(f'  [{conta}] task {tid} → {destino}/{nome}')
            registro.append({**linha, 'status': 'exportando', 'conta': conta, 'task_id': tid})

    if etapa == 3:
        gravar_retas(retas_cel)
    gravar_registro(registro, f'etapa{etapa}')


if __name__ == '__main__':
    main()
