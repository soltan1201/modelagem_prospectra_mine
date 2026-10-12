"""
Mosaicos só com as cenas zero-cloud (inspecionadas), por célula da grid 5 × 5,
em dois períodos:

  2000_2008 → cenas de 2000-01-01 a 2008-12-31, com SWIR (B04–B09)
  2009_2020 → cenas de 2009-01-01 a 2020-12-31, sem SWIR (detector desligado)

Regra de composição: mosaico comum — a cena mais recente fica por cima
(ImageCollection.sort por data + mosaic()). Foi a regra que se mostrou melhor
no visualizar_zero_cloud_minB3N.js (contra "menor B3N").

As cenas entram intactas, como no build semestral: só DN → radiância → TOA × 1000
(radiance_to_toa de build_mosaicos_semestrais.py, sem máscara de saturação),
TIR em DN e a banda 'saturado'. Uma banda 'data' (AAAAMMDD, INT32) guarda de
que cena veio cada pixel.

Máscara por pixel igual em todas as bandas: o mosaic() do GEE compõe banda a
banda; sem isso, uma cena sem SWIR (2008-S2) deixaria passar o SWIR de uma cena
mais antiga sob o VNIR dela. No período 2000_2008 o pixel só vale se tiver VNIR
e SWIR — na prática, as 2 cenas de 2008-S2 (sem SWIR) não entram nesse período.

Lista de cenas: as linhas ATIVAS (não comentadas) do LISTA_ZERO_CLOUD em
src/show_mosaics/visualizar_zero_cloud.js — onde a lista é revisada visualmente.
Se esse arquivo não existir, usa list_zero_cloud.txt.

Contas e destino (perguntado no início, como no build semestral):
  superconta             → DESTINOS['superconta']
  rodízio caatinga01–05  → DESTINOS['caatinga'] (troca a cada 2 tasks)
Nome: ASTER_zero_cloud_<periodo>_L<linha>C<coluna>, ex.: ASTER_zero_cloud_2000_2008_L3C2

Uso: python build_mosaicos_zero_cloud.py
"""

import csv
import logging
import re
from datetime import datetime
from types import SimpleNamespace

import ee

import build_mosaicos_semestrais as b

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
DESTINOS = {
    'superconta': 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA/ASTER_0cloud',
    'caatinga':   'projects/mapbiomas-arida/mosaic_aster_zero_cloud',
}
PERIODOS = {
    '2000_2008': {'inicio': '2000-01-01', 'fim': '2009-01-01', 'swir': True},
    '2009_2020': {'inicio': '2009-01-01', 'fim': '2021-01-01', 'swir': False},
}
ARQ_LISTA_JS = b.PASTA.parent / 'show_mosaics' / 'visualizar_zero_cloud.js'
PREFIXO = 'ASTER_zero_cloud_'

VNIR = ['B01', 'B02', 'B3N']
SWIR = ['B04', 'B05', 'B06', 'B07', 'B08', 'B09']
# SCALE_FACTORS gravado no asset: valor no asset = valor físico × fator
SCALE_FACTORS = {**dict.fromkeys(VNIR + SWIR, b.ESCALA_TOA), **dict.fromkeys(b.TIR, 1),
                 'saturado': 1, 'data': 1}


# ═══════════════════════════════════════════════════════════════════════════
# LISTA DE CENAS
# ═══════════════════════════════════════════════════════════════════════════
def ler_cenas() -> list[str]:
    """IDs ativos (linhas sem //) do visualizador; senão, o .txt. Sem blacklist."""
    if ARQ_LISTA_JS.exists():
        texto = ARQ_LISTA_JS.read_text()
        bloco = texto[texto.index('var LISTA_ZERO_CLOUD'):]
        bloco = bloco[:bloco.index('\n};')]
        ativas = [l.split('//')[0] for l in bloco.splitlines()]
        ids = set(re.findall(r"'(\d{14})'", '\n'.join(ativas)))
        origem = ARQ_LISTA_JS.name
    else:
        ids = b.ler_ids('list_zero_cloud.txt')
        origem = 'list_zero_cloud.txt'
    ids -= b.ler_ids('lista_imagens_Notaceites.txt')
    log.info(f'{len(ids)} cenas zero-cloud ativas (de {origem}).')
    return sorted(ids)


# ═══════════════════════════════════════════════════════════════════════════
# MOSAICO DE UMA CÉLULA × PERÍODO
# ═══════════════════════════════════════════════════════════════════════════
def preparar(img, com_swir: bool):
    """TOA × 1000 intacta + banda 'data'; máscara igual em todas as bandas."""
    img = ee.Image(img)
    t = ee.Image(b.radiance_to_toa(img, mascarar_saturacao=False))
    data = ee.Number.parse(ee.Date(img.get('system:time_start')).format('YYYYMMdd'))
    bandas = VNIR + (SWIR if com_swir else []) + b.TIR
    valido = t.select(bandas).mask().reduce(ee.Reducer.min())
    return (t.select(bandas + ['saturado'])
             .addBands(ee.Image.constant(data).toInt32().rename('data'))
             .updateMask(valido))


def montar(ids: list[str], periodo: str, celula: dict):
    cfg = PERIODOS[periodo]
    col = (ee.ImageCollection([ee.Image(f'ASTER/AST_L1T_003/{i}') for i in ids])
             .filterDate(cfg['inicio'], cfg['fim'])
             .filterBounds(celula['geom'])
             .sort('system:time_start'))           # mosaic(): a última (mais recente) fica por cima
    mosaico = col.map(lambda i: preparar(i, cfg['swir'])).mosaic().clip(celula['geom'])
    bandas = VNIR + (SWIR if cfg['swir'] else []) + b.TIR + ['saturado', 'data']
    # radiance_to_toa já entrega TOA × 1000: só arredonda para INT32
    mosaico = ee.Image.cat([mosaico.select(x).round().toInt32().rename(x) for x in bandas])
    mosaico = mosaico.set({
        'periodo': periodo, 'celula_id': celula['id'], 'celula_nome': celula['nome'],
        'num_cenas': col.size(), 'cenas': col.aggregate_array('system:index').join(','),
        'regra': 'mosaic (mais recente por cima)', 'grid_asset': b.GRID_ASSET,
        'SCALE_FACTORS': {k: SCALE_FACTORS[k] for k in bandas},
        'DATA_TYPE': 'INT32', 'versao': 'zero_cloud_v1',
    })
    return ee.Image(mosaico), col.size()


# ═══════════════════════════════════════════════════════════════════════════
# CONTAS, DESTINO, TASKS
# ═══════════════════════════════════════════════════════════════════════════
def destino_da_conta(conta: str) -> str:
    return DESTINOS['superconta'] if conta == 'superconta' else DESTINOS['caatinga']


def tasks_zero_cloud() -> dict[str, dict]:
    estados = {}
    for t in ee.data.getTaskList():
        desc = t.get('description', '')
        if desc.startswith(PREFIXO) and desc not in estados:
            estados[desc] = t
    return estados


def estado(contas: list[str]) -> tuple[dict, dict]:
    """Assets do destino e tasks das contas do grupo (lê na ordem do rodízio)."""
    existentes, tasks = {}, {}
    for conta in contas:
        b.usar_conta(conta)
        for desc, t in tasks_zero_cloud().items():
            if desc not in tasks or t.get('creation_timestamp_ms', 0) > tasks[desc].get('creation_timestamp_ms', 0):
                tasks[desc] = {**t, 'conta': conta}
    existentes = b.listar_destino(destino_da_conta(contas[0]), criar=True)
    if b.ARQ_CONTA_ATUAL.read_text().strip() != contas[0]:
        b.usar_conta(contas[0])
    return existentes, tasks


def coletar(grid: dict) -> SimpleNamespace:
    print('\n══════ Mosaicos zero-cloud por período e célula — configuração ══════')
    periodos = b.perguntar_opcao('Períodos:', [
        ('os dois (2000_2008 e 2009_2020)', list(PERIODOS)),
        ('só 2000_2008 (com SWIR)', ['2000_2008']),
        ('só 2009_2020 (sem SWIR)', ['2009_2020']),
    ])
    celulas = b.perguntar_celulas(grid)
    dry_run = b.perguntar_opcao('Modo:', [
        ('só contagens de cenas (dry-run, não exporta)', True),
        ('exportar os mosaicos', False),
    ])
    return SimpleNamespace(periodos=periodos, celulas=celulas, dry_run=dry_run)


# ═══════════════════════════════════════════════════════════════════════════
def main():
    contas = b.perguntar_contas()
    existentes, tasks = estado(contas)
    grid = b.ler_grid()
    args = coletar(grid)
    ids = ler_cenas()
    destino = destino_da_conta(contas[0])
    ativas = {d for d, t in tasks.items() if t['state'] in ('READY', 'RUNNING')}

    itens = [(p, c) for p in args.periodos for c in args.celulas]
    pendentes = [(p, c) for p, c in itens
                 if args.dry_run or (f'{PREFIXO}{p}_{grid[c]["nome"]}' not in existentes
                                     and f'{PREFIXO}{p}_{grid[c]["nome"]}' not in ativas)]
    print(f'\nDestino: {destino}')
    print(f'{len(itens)} itens pedidos; {len(itens) - len(pendentes)} já salvos ou em processamento; '
          f'{len(pendentes)} a {"contar" if args.dry_run else "exportar"}.')
    falhas = [n for n, t in tasks.items() if t['state'] == 'FAILED' and n not in existentes]
    if falhas:
        print(f'  última task falhou em: {", ".join(sorted(falhas))} (reenviados se pedidos)')
    if not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Cancelado.')

    montados = {(p, c): montar(ids, p, {'id': c, **grid[c]}) for p, c in pendentes}
    contagens = ee.List([n for _, n in montados.values()]).getInfo() if montados else []

    registro, disparadas, conta_ativa = [], 0, contas[0]
    for ((p, c), (mosaico, _)), n in zip(montados.items(), contagens):
        nome = f'{PREFIXO}{p}_{grid[c]["nome"]}'
        linha = {'periodo': p, 'celula_id': c, 'celula_nome': grid[c]['nome'], 'num_cenas': n}
        log.info(f'{nome}: {n} cenas')
        if n == 0:
            registro.append({**linha, 'status': 'sem_cenas'})
            continue
        if args.dry_run:
            registro.append({**linha, 'status': 'dry_run'})
            continue
        conta = contas[(disparadas // b.TASKS_POR_CONTA) % len(contas)]
        if conta != conta_ativa:
            b.usar_conta(conta)
            conta_ativa = conta
        asset = f'{destino_da_conta(conta)}/{nome}'
        task = ee.batch.Export.image.toAsset(
            image=mosaico, description=nome, assetId=asset, region=grid[c]['geom'],
            crs=b.CRS_EXPORT, crsTransform=b.CRS_TRANSFORM, maxPixels=1e13,
            pyramidingPolicy={'.default': 'sample'})
        try:
            task.start()
        except ee.EEException as e:
            log.error(f'  {conta}: falhou ao iniciar {nome}: {e}')
            registro.append({**linha, 'status': f'erro_inicio: {str(e)[:80]}', 'conta': conta})
            continue
        disparadas += 1
        log.info(f'  task {task.id} iniciada por {conta} → {asset}')
        registro.append({**linha, 'status': 'exportando', 'asset': asset, 'conta': conta,
                         'task_id': task.id})

    saida = b.PASTA / f'registro_zero_cloud_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(saida, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['periodo', 'celula_id', 'celula_nome', 'num_cenas',
                                          'status', 'asset', 'conta', 'task_id'], extrasaction='ignore')
        w.writeheader()
        w.writerows(registro)
    log.info(f'Registro → {saida}')


if __name__ == '__main__':
    main()
