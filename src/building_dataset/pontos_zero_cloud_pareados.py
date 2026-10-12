"""
Pontos pareados nos mosaicos zero-cloud: mesmas localizações em 2000–2008 (com
SWIR) e em 2009–2020 (sem SWIR).

Etapa 1 → ASSET_00_08
  Por célula da grid 5 × 5: mosaico zero-cloud 2000–2008 (montar() de
  build_mosaicos_zero_cloud.py: mosaic(), mais recente por cima, cenas intactas,
  TOA × 1000, TIR em DN). Pixel válido = VNIR + SWIR e não saturado (nas cenas
  zero-cloud o saturado se mostrou nuvem). Sorteia N_POR_CELULA pontos a 30 m
  dentro da célula encolhida 100 m (sem repetir pontos na faixa de
  sobreposição). Uma task toAsset com todas as células.

Etapa 2 → ASSET_09_20
  Só depois que ASSET_00_08 existir: amostra o mosaico zero-cloud 2009–2020 nos
  MESMOS pontos (sampleRegions, 30 m), célula por célula. Ponto sem dado em
  2009–2020 (fora das cenas, ou saturado) fica de fora.

Colunas: pid (igual nas duas tabelas, para juntar), celula_id, celula_nome,
bandas (B01, B02, B3N, [B04–B09], B10–B14; reflectância × 1000, TIR em DN) e
data (AAAAMMDD da cena que deu o pixel).

Conta: a ATIVA no momento (precisa escrever em projects/mapbiomas-arida — use
uma caatinga). O script não troca de conta.

Uso: python pontos_zero_cloud_pareados.py
"""

import logging
import sys

import ee

import build_mosaicos_semestrais as b
import build_mosaicos_zero_cloud as z

log = logging.getLogger(__name__)

sys.path.append(str(b.PASTA.parent))
from configure_account_projects_ee import get_current_account   # noqa: E402

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
ASSET_00_08  = 'projects/mapbiomas-arida/mine/points_zeroCC_00_08'
ASSET_09_20  = 'projects/mapbiomas-arida/mine/points_zeroCC_09_20'
N_POR_CELULA = 2000
ESCALA       = 30       # m — resolução nativa do SWIR
RECUO_M      = 100      # encolhe a célula: tira a faixa de sobreposição da grid
SEMENTE      = 42

BANDAS = {
    '2000_2008': z.VNIR + z.SWIR + b.TIR + ['data'],
    '2009_2020': z.VNIR + b.TIR + ['data'],
}
PROPS_PONTO = ['pid', 'celula_id', 'celula_nome']


def mosaico_valido(ids, periodo, cid, grid):
    """Mosaico zero-cloud do período na célula, sem pixels saturados."""
    img, n = z.montar(ids, periodo, {'id': cid, **grid[cid]})
    return img.updateMask(img.select('saturado').eq(0)).select(BANDAS[periodo]), n


def celulas_com_cenas(ids, periodo, grid) -> list[int]:
    """Células com ≥ 1 cena no período (uma chamada ao GEE). Célula sem cena daria
    um mosaic() sem bandas, e o select() quebra a task inteira."""
    cids = sorted(grid)
    n = ee.List([z.montar(ids, periodo, {'id': c, **grid[c]})[1] for c in cids]).getInfo()
    vazias = [grid[c]['nome'] for c, k in zip(cids, n) if k == 0]
    if vazias:
        log.info(f'{periodo}: sem cenas em {", ".join(vazias)} — células puladas.')
    return [c for c, k in zip(cids, n) if k > 0]


def existe(asset: str) -> bool:
    try:
        ee.data.getAsset(asset)
        return True
    except ee.EEException:
        return False


def tarefa(nome: str):
    for t in ee.data.getTaskList():
        if t.get('description') == nome:
            return t
    return None


# ═══════════════════════════════════════════════════════════════════════════
# ETAPAS
# ═══════════════════════════════════════════════════════════════════════════
def etapa1(ids, grid):
    if existe(ASSET_00_08):
        raise SystemExit(f'{ASSET_00_08} já existe — apague antes para refazer.')
    partes = []
    for cid in celulas_com_cenas(ids, '2000_2008', grid):
        cel = grid[cid]
        img, _ = mosaico_valido(ids, '2000_2008', cid, grid)
        # stratifiedSample sobre uma banda constante só nos pixels válidos: sorteia
        # N_POR_CELULA pontos COM dado (sample(numPixels) sorteia na área toda e
        # descarta os sem dado — célula com 17% de cobertura daria ~17% dos pontos)
        valido = ee.Image(1).rename('valido').updateMask(img.select('B01').mask()).toInt()
        pts = (img.addBands(valido)
                  .stratifiedSample(numPoints=N_POR_CELULA, classBand='valido',
                                    region=cel['geom'].buffer(-RECUO_M), scale=ESCALA,
                                    seed=SEMENTE + cid, geometries=True, dropNulls=True, tileScale=4)
                  .select(BANDAS['2000_2008']))
        partes.append(pts.map(lambda f, cid=cid, nome=cel['nome']: f.set({
            'pid': ee.String(nome).cat('_').cat(f.id()), 'celula_id': cid, 'celula_nome': nome})))
    fc = ee.FeatureCollection(partes).flatten()
    nome = ASSET_00_08.rsplit('/', 1)[-1]
    task = ee.batch.Export.table.toAsset(collection=fc, description=nome, assetId=ASSET_00_08)
    task.start()
    log.info(f'Etapa 1: task {task.id} → {ASSET_00_08} (até {N_POR_CELULA} pontos × {len(grid)} células)')


def etapa2(ids, grid):
    if not existe(ASSET_00_08):
        t = tarefa(ASSET_00_08.rsplit('/', 1)[-1])
        raise SystemExit(f'{ASSET_00_08} ainda não existe'
                         + (f' (task da etapa 1: {t["state"]})' if t else '') + ' — rode a etapa 2 depois.')
    if existe(ASSET_09_20):
        raise SystemExit(f'{ASSET_09_20} já existe — apague antes para refazer.')
    pontos = ee.FeatureCollection(ASSET_00_08).select(PROPS_PONTO)
    partes = []
    for cid in celulas_com_cenas(ids, '2009_2020', grid):
        img, _ = mosaico_valido(ids, '2009_2020', cid, grid)
        partes.append(img.sampleRegions(collection=pontos.filter(ee.Filter.eq('celula_id', cid)),
                                        properties=PROPS_PONTO, scale=ESCALA,
                                        geometries=True, tileScale=4))
    fc = ee.FeatureCollection(partes).flatten()
    nome = ASSET_09_20.rsplit('/', 1)[-1]
    task = ee.batch.Export.table.toAsset(collection=fc, description=nome, assetId=ASSET_09_20)
    task.start()
    log.info(f'Etapa 2: task {task.id} → {ASSET_09_20} (mesmos pontos de {ASSET_00_08})')


def status():
    for asset in (ASSET_00_08, ASSET_09_20):
        nome = asset.rsplit('/', 1)[-1]
        if existe(asset):
            fc = ee.FeatureCollection(asset)
            por_cel = fc.aggregate_histogram('celula_nome').getInfo()
            print(f'{nome}: {sum(por_cel.values())} pontos em {len(por_cel)} células')
            print('   ', ', '.join(f'{k} {v}' for k, v in sorted(por_cel.items())))
        else:
            t = tarefa(nome)
            print(f'{nome}: não existe' + (f' — task {t["state"]} {t.get("error_message", "")}' if t else ''))


# ═══════════════════════════════════════════════════════════════════════════
def main():
    conta = b.ARQ_CONTA_ATUAL.read_text().strip()
    ee.Initialize(project=get_current_account())
    log.info(f'Conta ativa: {conta}')
    acao = b.perguntar_opcao('O que fazer?', [
        (f'etapa 1 — pontos no mosaico 2000–2008 → {ASSET_00_08}', 'etapa1'),
        (f'etapa 2 — mesmos pontos no mosaico 2009–2020 → {ASSET_09_20}', 'etapa2'),
        ('status dos dois assets', 'status'),
    ])
    if acao == 'status':
        return status()
    grid = b.ler_grid()
    ids = z.ler_cenas()
    (etapa1 if acao == 'etapa1' else etapa2)(ids, grid)


if __name__ == '__main__':
    main()
