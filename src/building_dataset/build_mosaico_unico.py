"""
Mosaicos ASTER por célula da grid 5 × 5 (Irecê): dois mosaicos base separados
pela presença do SWIR, e o SWIR estimado para o mosaico que não tem.

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

Como combinar os dois depois fica em aberto (ex.: com_swir onde o NDVI é
parecido, sem_swir_estimado no resto) — os dois assets guardam 'quality' e
'data' para isso.

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
VERSAO        = 'unico_v2'

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
def montar_estimado(celula: dict, modelos: dict, coefs: dict):
    origem = f'{DESTINO_BASE}/{PREFIXO_BASE}sem_swir_{celula["nome"]}'
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
# CONTAS, ESTADO, REGISTRO
# ═══════════════════════════════════════════════════════════════════════════
def estado(contas: list[str], destino: str, prefixo: str) -> tuple[set, set]:
    """Assets já salvos no destino e tasks ativas com o prefixo nas contas escolhidas."""
    ativas = set()
    for conta in contas:
        b.usar_conta(conta)
        ativas |= {t['description'] for t in ee.data.getTaskList()
                   if t['state'] in ('READY', 'RUNNING') and t.get('description', '').startswith(prefixo)}
    salvos = set(b.listar_destino(destino, criar=True))
    if b.ARQ_CONTA_ATUAL.read_text().strip() != contas[0]:
        b.usar_conta(contas[0])
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
def main():
    print('══════ Mosaicos ASTER com / sem SWIR por célula ══════')
    etapa = b.perguntar_opcao('Etapa:', [
        (f'1 — mosaicos base com_swir e sem_swir → {DESTINO_BASE}', 1),
        (f'2 — SWIR estimado no mosaico sem_swir → {DESTINO_FINAL}', 2),
    ])
    contas = b.perguntar_contas()
    destino, prefixo = (DESTINO_BASE, PREFIXO_BASE) if etapa == 1 else (DESTINO_FINAL, PREFIXO_FINAL)
    salvos, ativas = estado(contas, destino, prefixo)
    grid = b.ler_grid()

    tipos, meses = [], None
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
    else:
        itens = {('sem_swir_estimado', c): f'{PREFIXO_FINAL}{grid[c]["nome"]}' for c in celulas}
    pendentes = [k for k, nome in itens.items() if dry_run or (nome not in salvos and nome not in ativas)]
    print(f'{len(itens)} assets pedidos; {len(itens) - len(pendentes)} já salvos ou em processamento; '
          f'{len(pendentes)} a {"conferir" if dry_run else "exportar"}.')
    if not pendentes or not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Nada a fazer.')

    registro, disparadas, conta_ativa = [], 0, contas[0]

    def proxima_conta():
        nonlocal conta_ativa
        conta = contas[(disparadas // b.TASKS_POR_CONTA) % len(contas)]
        if conta != conta_ativa:
            b.usar_conta(conta)
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
                tid = exportar(img, nome, DESTINO_BASE, {'id': c, **grid[c]}, {'.default': 'sample'})
            except ee.EEException as e:
                registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
                continue
            disparadas += 1
            log.info(f'  [{conta}] task {tid} → {DESTINO_BASE}/{nome}')
            registro.append({**linha, 'status': 'exportando', 'conta': conta, 'task_id': tid})
    else:
        bases = set(b.listar_destino(DESTINO_BASE, criar=False))
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
                tid = exportar(montar_estimado({'id': c, **grid[c]}, modelos, coefs), nome, DESTINO_FINAL,
                               {'id': c, **grid[c]},
                               {'.default': 'mean', 'data': 'mode', 'camada': 'mode', 'normalizacao': 'mode'})
            except ee.EEException as e:
                registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
                continue
            disparadas += 1
            log.info(f'  [{conta}] task {tid} → {DESTINO_FINAL}/{nome}')
            registro.append({**linha, 'status': 'exportando', 'conta': conta, 'task_id': tid})

    gravar_registro(registro, f'etapa{etapa}')


if __name__ == '__main__':
    main()
