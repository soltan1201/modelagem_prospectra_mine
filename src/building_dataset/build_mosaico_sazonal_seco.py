"""
Mosaico SAZONAL SECO ASTER (Irecê, BA) — qualityMosaic de todas as cenas secas
de um período, por célula da grid 5 × 5.

Versão Python (earthengine-api) do candidato "seco | quality" de
mosaico_sazonal_homogeneo.js, que foi o melhor resultado na comparação visual
(quality × mediana × mediana_norm). A lógica por cena é a MESMA do JS:
  - estação seca = jun, jul, ago, set, out de TODOS os anos do período;
  - DN → radiância (GAIN_COEFFICIENT_Bxx) → reflectância TOA 0–1 nas bandas
    ópticas do período; TIR em DN; pixel com algum VNIR ≥ 254 é mascarado;
  - máscara de nuvem/sombra só com VNIR + TIR (vale para os dois períodos);
  - quality v5 (VNIR + TIR, SAVI alto soma pontos, bônus por CLOUDCOVER baixo);
  - cenas de list_zero_cloud.txt que caem em jun–out: sem máscara de nuvem e
    quality v5 + BONUS_ZERO_CLOUD, então sempre vencem as cenas com nuvem; entre
    duas zero-cloud decide o quality v5 do pixel (cena na blacklist sai mesmo
    estando na lista zero-cloud);
  - cenas de lista_imagens_Notaceites.txt (blacklist)
    sempre ficam de fora — mesmo com CLOUDCOVER < 60 há cenas que estragam o
    mosaico. Só os IDs da estação seca do período entram no grafo.

Períodos:
  A — 2000-03 a 2008-03 (SWIR medido): VNIR + SWIR + TIR; só cenas com B09.
  B — 2009 a 2020 (sem SWIR): VNIR + TIR; o SWIR é estimado depois pelo RF.

Saída (mesma escala dos mosaicos v4/v5 e do export do JS): reflectância × 1000,
TIR em DN, banda n_obs (nº de cenas válidas no pixel), INT32, 15 m em EPSG:4326
com o CRS_TRANSFORM ancorado no canto NW do bloco (células encaixam pixel a
pixel). Nome: ASTER_sazonal_<A|B>_seco_quality_<LxCy> em DESTINO (temporário:
mapbiomas-brazil .../CAATINGA/ASTER_sazonal); depois copiados para
DESTINO_FINAL (projects/mapbiomas-arida/mine/ASTER_sazonal), lido pelas etapas
seguintes. Tudo roda com a conta CONTA inicializada no projeto PROJETO
(mapbiomas-brazil), sem rodízio.

Uso (interativo — pergunta períodos, células e modo):
  python build_mosaico_sazonal_seco.py
"""

import csv
import logging
import math
import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import ee

PASTA = Path(__file__).resolve().parent
sys.path.append(str(PASTA))
from build_mosaicos_semestrais_v2 import (   # noqa: E402  (grade, contas e perguntas)
    CRS_EXPORT, CRS_TRANSFORM, GRID_ASSET, ler_grid, ler_ids, listar_destino,
    perguntar_celulas, perguntar_opcao, perguntar_sim_nao)
from gee_tools import switch_user   # noqa: E402  (src/ entra no sys.path pelo import acima)

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS (iguais a mosaico_sazonal_homogeneo.js)
# ═══════════════════════════════════════════════════════════════════════════
CONTA = 'superconta' # 'caatinga05'
PROJETO = 'mapbiomas-brazil'       # projeto em que a conta roda (não o padrão dela)
# Export temporário em mapbiomas-brazil (superconta escreve lá). Depois os
# assets são copiados para DESTINO_FINAL, que é onde as etapas seguintes leem
# (pontos_sazonal_seco_pareados, estimar_swir_rf_v2, build_mosaico_final_v2).
DESTINO = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA/ASTER_sazonal'
DESTINO_FINAL = 'projects/mapbiomas-arida/mine/ASTER_sazonal'
ESTACAO = 'seco'
CANDIDATO = 'quality'
MESES_SECO = [6, 7, 8, 9, 10]

VNIR = ['B01', 'B02', 'B3N']
SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']
TIR = ['B10', 'B11', 'B12', 'B13', 'B14']

CONFIG = {
    'A': {'data_inicio': '2000-03-01',
          'data_fim': '2008-04-01',          # SWIR do ASTER falhou em abril de 2008
          'bandas_toa': VNIR + SWIR},
    'B': {'data_inicio': '2009-01-01',
          'data_fim': '2021-01-01',
          'bandas_toa': VNIR},
}

CC_MAX = 60                  # CLOUDCOVER < 60
ARQ_BLACKLIST = 'lista_imagens_Notaceites.txt'
ARQ_ZERO_CLOUD = 'list_zero_cloud.txt'
BONUS_ZERO_CLOUD = 3.0       # quality v5 vai até 2,1: zero-cloud sempre vence
ESCALA_REFL = 1000           # reflectância 0–1 → 0–1000 no asset
DN_SATURADO = 254

# Nuvem / sombra (reflectância 0–1)
LIM_VIS_BRIGHT = 0.19
RAIO_BORDA_NUVEM = 180       # m
LIM_ESCURO_GREEN = 0.12
LIM_ESCURO_NIR = 0.17

ESUN = {'B01': 1848, 'B02': 1549, 'B3N': 1114,
        'B04': 225.4, 'B05': 86.63, 'B06': 81.85,
        'B07': 74.85, 'B08': 66.49, 'B09': 59.85}


# ═══════════════════════════════════════════════════════════════════════════
# FUNÇÕES POR CENA (port de mosaico_sazonal_homogeneo.js)
# ═══════════════════════════════════════════════════════════════════════════
def radiance_to_toa(img, bandas_toa: list[str]):
    """DN → radiância → reflectância TOA 0–1 em `bandas_toa`; TIR segue em DN.
    Pixel com algum VNIR saturado (DN ≥ 254) sai mascarado."""
    img = ee.Image(img)
    cosz = ee.Number(90).subtract(ee.Number(img.get('SOLAR_ELEVATION'))) \
                        .multiply(math.pi / 180).cos()
    doy = ee.Number(ee.Date(img.get('system:time_start')).getRelative('day', 'year'))
    d = ee.Number(1).subtract(ee.Number(0.01672).multiply(doy.multiply(2 * math.pi / 365).cos()))

    opticas = []
    for b in bandas_toa:
        dn = img.select(b)
        ganho = ee.Number(img.get(f'GAIN_COEFFICIENT_{b}'))
        opticas.append(dn.subtract(1).multiply(ganho)
                         .multiply(math.pi).multiply(d.pow(2))
                         .divide(cosz.multiply(ESUN[b]))
                         .updateMask(dn.gt(0))
                         .rename(b).toFloat())

    nao_saturado = img.select(VNIR).reduce(ee.Reducer.max()).lt(DN_SATURADO)
    tir = img.select(TIR).toFloat()
    saida = ee.Image.cat(opticas).addBands(tir.updateMask(tir.gt(0))).updateMask(nao_saturado)
    return ee.Image(saida.copyProperties(img, img.propertyNames()))


def mascara_nuvem(img):
    """Nuvem (brilho visível + TIR frio relativo) e sombra (escuro perto de nuvem),
    só com VNIR + TIR."""
    img = ee.Image(img)
    green, red, nir = img.select('B01'), img.select('B02'), img.select('B3N')
    tir10 = img.select('B10')

    tan_elev = ee.Number(img.get('SOLAR_ELEVATION')).multiply(math.pi / 180).tan()
    shadow_dist = ee.Number(2000).divide(tan_elev.max(0.1)).max(200).min(2000)

    vis_bright = green.add(red).add(nir).divide(3)
    p20 = ee.Number(tir10.reduceRegion(
        reducer=ee.Reducer.percentile([20]), geometry=img.geometry(),
        scale=180, maxPixels=1e10, bestEffort=True).get('B10'))
    cloud = vis_bright.gt(LIM_VIS_BRIGHT).And(tir10.lt(p20))
    cloud = (cloud.focalMin(radius=RAIO_BORDA_NUVEM, kernelType='square', units='meters')
                  .focalMax(radius=RAIO_BORDA_NUVEM, kernelType='square', units='meters'))

    dark = green.lt(LIM_ESCURO_GREEN).And(nir.lt(LIM_ESCURO_NIR))
    shadow = dark.And(cloud.focalMax(radius=shadow_dist.add(300), units='meters')).focalMax()

    valid = img.select(VNIR + ['B10']).reduce(ee.Reducer.min()).gt(0)
    clear = valid.And(cloud.Not()).And(shadow.Not())
    return ee.Image(img.updateMask(clear).copyProperties(img, img.propertyNames()))


def add_quality_v5(img):
    """Quality v5 do JS (VNIR + TIR), reflectância 0–1, sem teto."""
    img = ee.Image(img)
    b1, b2, b3n, b10 = img.select('B01'), img.select('B02'), img.select('B3N'), img.select('B10')

    savi = b3n.subtract(b2).multiply(1.5).divide(b3n.add(b2).add(0.5).add(1e-6))
    ndwi = b1.subtract(b3n).divide(b1.add(b3n).add(1e-6))
    brightness = b1.add(b2).add(b3n).divide(3)

    tir_score = b10.unitScale(1050, 1300).clamp(0, 1)
    bright_up = brightness.unitScale(0.053, 0.122).clamp(0, 1)
    bright_down = ee.Image(1).subtract(brightness.unitScale(0.168, 0.32).clamp(0, 1))
    savi_score = savi.unitScale(-0.1, 0.5).clamp(0, 1)
    ndwi_score = ndwi.multiply(-1).unitScale(-0.10, 0.45).clamp(0, 1)

    base = (tir_score.multiply(0.45)
            .add(bright_up.multiply(bright_down).multiply(0.25))
            .add(savi_score.multiply(0.15))
            .add(ndwi_score.multiply(0.15))
            .clamp(0, 1))

    cc = ee.Number(img.get('CLOUDCOVER'))
    cloud_bonus = ee.Number(1).subtract(cc.divide(80).min(1)).multiply(0.8)
    quality = base.add(cloud_bonus).add(cc.eq(0).multiply(0.3))
    return img.addBands(quality.rename('quality').toFloat())


# ═══════════════════════════════════════════════════════════════════════════
# MOSAICO DE UM PERÍODO × CÉLULA
# ═══════════════════════════════════════════════════════════════════════════
def filtro_meses(meses: list[int]):
    return ee.Filter.Or(*[ee.Filter.calendarRange(m, m, 'month') for m in meses])


def ids_do_periodo(ids: set[str], periodo: str) -> list[str]:
    """Só os IDs (AAAAMMDDhhmmss) dentro das datas do período e nos meses secos:
    deixa o grafo de cada task bem menor do que a blacklist inteira."""
    cfg = CONFIG[periodo]
    ini, fim = cfg['data_inicio'].replace('-', ''), cfg['data_fim'].replace('-', '')
    return sorted(i for i in ids if ini <= i[:8] < fim and int(i[4:6]) in MESES_SECO)


def colecao_seca(periodo: str, area, blacklist: ee.List, zero_cloud: ee.List):
    """Cenas secas do período, já em TOA, mascaradas e com a banda quality, na
    mesma ordem e com os mesmos filtros do mosaico. Também usada por
    harmonizar_sazonal_modis.py para reconstruir de qual cena veio cada pixel.
    Devolve (base, filtrada, cc0, col)."""
    cfg = CONFIG[periodo]
    base = (ee.ImageCollection('ASTER/AST_L1T_003')
              .filterDate(cfg['data_inicio'], cfg['data_fim'])
              .filterBounds(area)
              .filter(filtro_meses(MESES_SECO))
              .filter(ee.Filter.lt('CLOUDCOVER', CC_MAX))
              .filter(ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B3N'))
              .filter(ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B10')))
    if periodo == 'A':
        base = base.filter(ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B09'))
    filtrada = base.filter(ee.Filter.inList('system:index', blacklist).Not())

    na_lista = ee.Filter.inList('system:index', zero_cloud)
    toa = lambda i: radiance_to_toa(i, cfg['bandas_toa'])
    cc0 = (filtrada.filter(na_lista).map(toa).map(add_quality_v5)
                   .map(lambda i: i.addBands(i.select('quality').add(BONUS_ZERO_CLOUD), None, True)))
    nuvens = filtrada.filter(na_lista.Not()).map(toa).map(mascara_nuvem).map(add_quality_v5)
    return base, filtrada, cc0, cc0.merge(nuvens)


def montar_seco(periodo: str, area, blacklist: ee.List, zero_cloud: ee.List):
    cfg = CONFIG[periodo]
    bandas = cfg['bandas_toa'] + TIR
    base, filtrada, cc0, col = colecao_seca(periodo, area, blacklist, zero_cloud)
    obs = col.select('B01').count().rename('n_obs')
    mosaico = col.qualityMosaic('quality').select(bandas)

    saida = (mosaico.select(cfg['bandas_toa']).multiply(ESCALA_REFL)
                    .addBands(mosaico.select(TIR))
                    .round().toInt32()
                    .addBands(obs.toInt32())
                    .clip(area)
                    .set({
                        'periodo': periodo,
                        'anos': f"{cfg['data_inicio']}_{cfg['data_fim']}",
                        'estacao': ESTACAO,
                        'meses': ','.join(map(str, MESES_SECO)),
                        'candidato': CANDIDATO,
                        'quality': 'v5 (JS mosaico_sazonal_homogeneo)',
                        'cc_max': CC_MAX,
                        'blacklist': ARQ_BLACKLIST,
                        'num_images': filtrada.size(),
                        'num_zero_cloud': cc0.size(),
                        'bonus_zero_cloud': BONUS_ZERO_CLOUD,
                        'swir': 'medido' if periodo == 'A' else 'ausente (estimar com RF)',
                        'SCALE_FACTORS': f'reflectancia x{ESCALA_REFL}; TIR em DN',
                        'DATA_TYPE': 'INT32',
                    }))
    contagens = {'cenas_seco': base.size(), 'apos_blacklist': filtrada.size(),
                 'zero_cloud': cc0.size()}
    return ee.Image(saida), contagens


# ═══════════════════════════════════════════════════════════════════════════
# CONTAS, ESTADO E PLANO
# ═══════════════════════════════════════════════════════════════════════════
def nome_mosaico(periodo: str, celula: str) -> str:
    return f'ASTER_sazonal_{periodo}_{ESTACAO}_{CANDIDATO}_{celula}'


def estado_geral() -> tuple[dict, dict]:
    """Assets já em DESTINO e a última task ASTER_sazonal_*_seco_quality_* da
    conta, para não exportar duas vezes."""
    switch_user(CONTA)
    ee.Initialize(project=PROJETO)
    log.info(f'Earth Engine inicializado com {CONTA} (projeto {PROJETO}).')
    tasks = {}
    for t in ee.data.getTaskList():        # da mais recente para a mais antiga
        desc = t.get('description', '')
        if (desc.startswith('ASTER_sazonal_') and f'_{ESTACAO}_{CANDIDATO}_' in desc
                and desc not in tasks):
            tasks[desc] = t
    return listar_destino(DESTINO, criar=True), tasks


def planejar(args, grid, existentes, tasks) -> dict[tuple, str]:
    ativas = {d for d, t in tasks.items() if t['state'] in ('READY', 'RUNNING')}
    print(f'\n══════ Destino: {DESTINO} ══════')
    plano = {}
    for periodo in args.periodos:
        for cid in args.celulas:
            nome = nome_mosaico(periodo, grid[cid]['nome'])
            if nome in ativas:
                plano[(periodo, cid)] = 'processando'
            elif nome in existentes:
                plano[(periodo, cid)] = 'refazer' if args.sobrescrever else 'existe'
            else:
                plano[(periodo, cid)] = 'novo'
    for periodo in args.periodos:
        cats = [c for (p, _), c in plano.items() if p == periodo]
        resumo = ', '.join(f'{cats.count(c)} {c}' for c in ('novo', 'refazer', 'existe', 'processando')
                           if cats.count(c))
        print(f'  período {periodo}: {resumo}')
    falhas = [(n, t.get('error_message', '')) for n, t in tasks.items()
              if t['state'] == 'FAILED' and n not in existentes and n not in ativas]
    for n, erro in sorted(falhas):
        print(f'  última task FALHOU: {n}: {" ".join(erro.split())[:90]}')
    return plano


def coletar_opcoes(grid) -> SimpleNamespace:
    print('\n══════ Mosaico sazonal SECO (quality) ASTER — configuração ══════')
    periodos = perguntar_opcao('Períodos:', [
        ('A e B (2000–2008 e 2009–2020)', ['A', 'B']),
        ('só A (2000-03 a 2008-03, com SWIR)', ['A']),
        ('só B (2009–2020, sem SWIR)', ['B']),
    ])
    celulas = perguntar_celulas(grid)
    dry_run = perguntar_opcao('Modo:', [
        ('só contagens de cenas (dry-run, não exporta)', True),
        ('exportar os mosaicos', False),
    ])
    sobrescrever = False
    if not dry_run:
        sobrescrever = perguntar_sim_nao('Se o asset já existir, apagar e exportar de novo?')
        if sobrescrever:
            print('  ATENÇÃO: os assets existentes nesses itens serão APAGADOS (não tem volta).')
            sobrescrever = input('  Digite APAGAR para confirmar: ').strip() == 'APAGAR'
            if not sobrescrever:
                print('  Confirmação não recebida — assets existentes serão pulados.')
    return SimpleNamespace(periodos=periodos, celulas=celulas,
                           dry_run=dry_run, sobrescrever=sobrescrever)


# ═══════════════════════════════════════════════════════════════════════════
# EXECUÇÃO
# ═══════════════════════════════════════════════════════════════════════════
def main():
    existentes, tasks = estado_geral()
    grid = ler_grid()
    args = coletar_opcoes(grid)
    plano = planejar(args, grid, existentes, tasks)
    print(f'\nModo: {"dry-run (só contagens)" if args.dry_run else "EXPORTAR"} | '
          f'conta: {CONTA} | sobrescrever={args.sobrescrever}')
    if not perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Cancelado.')
    executar(args, grid, plano)


def executar(args, grid, plano: dict[tuple, str]):
    ids_black = ler_ids(ARQ_BLACKLIST)
    ids_zero = ler_ids(ARQ_ZERO_CLOUD)
    conflito = ids_black & ids_zero
    if conflito:
        log.warning(f'{len(conflito)} cena(s) nas duas listas — ficam excluídas: {sorted(conflito)}')
    ids_zero -= ids_black
    log.info(f'Blacklist: {len(ids_black)} cenas | zero-cloud: {len(ids_zero)} cenas')

    registro = [{'periodo': p, 'celula_id': cid, 'celula_nome': grid[cid]['nome'],
                 'asset': f'{DESTINO}/{nome_mosaico(p, grid[cid]["nome"])}', 'status': cat}
                for (p, cid), cat in plano.items()
                if cat in ('existe', 'processando') and not args.dry_run]
    pendentes = [k for k, c in plano.items() if args.dry_run or c in ('novo', 'refazer')]
    if not pendentes:
        log.info('Nada a exportar.')

    disparadas = 0
    for periodo in args.periodos:
        celulas = [cid for p, cid in pendentes if p == periodo]
        if not celulas:
            continue
        blacklist = ids_do_periodo(ids_black, periodo)
        zero = ids_do_periodo(ids_zero, periodo)
        log.info(f'Período {periodo}, estação seca: {len(blacklist)} cenas da blacklist, '
                 f'{len(zero)} zero-cloud')
        montados = {cid: montar_seco(periodo, grid[cid]['geom'], ee.List(blacklist), ee.List(zero))
                    for cid in celulas}
        contagens = ee.List([ee.Dictionary(c) for _, c in montados.values()]).getInfo()

        for (cid, (mosaico, _)), cont in zip(montados.items(), contagens):
            nome = nome_mosaico(periodo, grid[cid]['nome'])
            linha = {'periodo': periodo, 'celula_id': cid, 'celula_nome': grid[cid]['nome'],
                     'asset': f'{DESTINO}/{nome}', **cont}
            cat = plano[(periodo, cid)]
            log.info(f'{periodo} {grid[cid]["nome"]} [{cat}]: {cont}')

            if cont['apos_blacklist'] == 0:
                registro.append({**linha, 'status': 'sem_cenas'})
                continue
            if args.dry_run:
                registro.append({**linha, 'status': f'dry_run ({cat})'})
                continue

            mosaico = ee.Image(mosaico.set({'celula_id': cid, 'celula_nome': grid[cid]['nome'],
                                            'grid_asset': GRID_ASSET}))
            # monta a task antes de apagar: se a montagem falhar, o asset antigo fica intacto
            task = ee.batch.Export.image.toAsset(
                image=mosaico, description=nome, assetId=linha['asset'],
                region=grid[cid]['geom'], crs=CRS_EXPORT, crsTransform=CRS_TRANSFORM,
                maxPixels=1e13, pyramidingPolicy={'.default': 'mean'})
            if cat == 'refazer':
                log.warning(f'Apagando {linha["asset"]} (sobrescrever = sim).')
                ee.data.deleteAsset(linha['asset'])
            try:
                task.start()
            except ee.EEException as e:
                log.error(f'  {CONTA}: falhou ao iniciar {nome}: {e}')
                registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': CONTA})
                continue
            disparadas += 1
            log.info(f'  task {task.id} iniciada por {CONTA} ({disparadas}) → {nome}')
            registro.append({**linha, 'status': 'exportando', 'task_id': task.id, 'conta': CONTA})

    saida = PASTA / f'registro_sazonal_seco_{datetime.now():%Y%m%d_%H%M}.csv'
    campos = ['periodo', 'celula_id', 'celula_nome', 'status', 'cenas_seco', 'apos_blacklist',
              'zero_cloud',
              'asset', 'conta', 'task_id']
    with open(saida, 'w', newline='') as f:
        escritor = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        escritor.writeheader()
        escritor.writerows(registro)
    log.info(f'Registro → {saida}')


if __name__ == '__main__':
    main()
