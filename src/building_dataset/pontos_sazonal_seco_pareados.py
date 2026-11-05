"""
Pontos pareados nos mosaicos SAZONAIS SECOS (build_mosaico_sazonal_seco.py):
mesmas localizações no período A (2000-03 a 2008-03, com SWIR) e no período B
(2009–2020, sem SWIR), para treinar e aplicar a estimativa do SWIR do B.

Mesmo desenho de pontos_zero_cloud_pareados.py, trocando a fonte: em vez de
montar o mosaico zero-cloud na hora, lê os assets já exportados
  <COL_FONTE>/ASTER_sazonal_A_seco_quality_<LxCy>
  <COL_FONTE>/ASTER_sazonal_B_seco_quality_<LxCy>

FONTE escolhe a coleção (FONTES): 'hm_razao' (padrão, harmonizado contra o Landsat,
harmonizar_sazonal_landsat.py, modo razão), 'hm_residuo' ou 'original'
(ASTER_sazonal, TOA). Cada fonte tem os seus assets de pontos
(points_sazonal_seco<SUFIXO>_A/_B) e a sua pasta local
(ASTER_Samples_sazonal_seco<SUFIXO>), para não misturar.

Etapa 1 → ASSET_A
  Por célula com asset A: sorteia N_POR_CELULA pontos a 30 m (resolução do
  SWIR) dentro da célula encolhida 100 m (sem repetir pontos na faixa de
  sobreposição), só onde VNIR + SWIR + TIR têm dado. Uma task toAsset.

Etapa 2 → ASSET_B
  Só depois que ASSET_A existir: amostra o asset B nos MESMOS pontos
  (sampleRegions, 30 m). Ponto sem dado no B fica de fora.

Baixar → src/Dados/ASTER_Samples_sazonal_seco<SUFIXO>/
  points_sazonal_seco_A.csv e _B.csv (cache por célula) e a junção pelo pid
  (points_sazonal_seco_pareados.csv; bandas com sufixo _A / _B).

Colunas: pid, celula_id, celula_nome, bandas (B01, B02, B3N, [B04–B09],
B10–B14; reflectância × 1000, TIR em DN) e n_obs (nº de cenas válidas no pixel).
Os mosaicos não têm banda de data (o pixel vem do qualityMosaic de várias
cenas), então não há normalização por cena: a diferença radiométrica B × A
tem de ser vista nos próprios pares.

Conta: CONTA abaixo (caatinga04, a mesma da cópia brazil → arida), no projeto
dela. Não usa a do build_mosaico_sazonal_seco.py: a superconta não lê
mapbiomas-arida. Tabelas na pasta-mãe de DESTINO_FINAL (mapbiomas-arida/mine).

Uso: python pontos_sazonal_seco_pareados.py
"""

import logging

import ee
import pandas as pd

import build_mosaico_sazonal_seco as s
from configure_account_projects_ee import get_project_from_account
from build_mosaicos_semestrais_v2 import ler_grid, perguntar_opcao
from gee_tools import switch_user

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CONTA        = 'caatinga04'   # lê e escreve em mapbiomas-arida/mine
PASTA_ASSETS = s.DESTINO_FINAL.rsplit('/', 1)[0]        # projects/mapbiomas-arida/mine
# fonte dos mosaicos → (coleção, sufixo dos pontos e da pasta local)
FONTES = {
    'hm_razao':   (f'{PASTA_ASSETS}/ASTER_sazonal_hm_razao', '_hm_razao'),
    'hm_residuo': (f'{PASTA_ASSETS}/ASTER_sazonal_hm_residuo', '_hm_residuo'),
    'original':   (s.DESTINO_FINAL, ''),
}
FONTE = 'hm_razao'
COL_FONTE, SUFIXO = FONTES[FONTE]
ASSET_A      = f'{PASTA_ASSETS}/points_sazonal_seco{SUFIXO}_A'
ASSET_B      = f'{PASTA_ASSETS}/points_sazonal_seco{SUFIXO}_B'
N_POR_CELULA = 2000
ESCALA       = 30       # m — resolução nativa do SWIR
RECUO_M      = 100      # encolhe a célula: tira a faixa de sobreposição da grid
SEMENTE      = 42

BANDAS = {
    'A': s.VNIR + s.SWIR + s.TIR + ['n_obs'],
    'B': s.VNIR + s.TIR + ['n_obs'],
}
TABELAS = {'A': ASSET_A, 'B': ASSET_B}
PROPS_PONTO = ['pid', 'celula_id', 'celula_nome']
CHAVES = PROPS_PONTO + ['longitude', 'latitude']
SAIDA = s.PASTA.parent / 'Dados' / f'ASTER_Samples_sazonal_seco{SUFIXO}'


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


def mosaicos(periodo: str, grid) -> dict[int, ee.Image]:
    """Assets do período por célula; célula sem asset é pulada (avisa)."""
    salvos = s.listar_destino(COL_FONTE, criar=False)
    imgs, faltam = {}, []
    for cid in sorted(grid):
        nome = s.nome_mosaico(periodo, grid[cid]['nome'])
        if nome in salvos:
            imgs[cid] = ee.Image(f'{COL_FONTE}/{nome}').select(BANDAS[periodo])
        else:
            faltam.append(grid[cid]['nome'])
    if faltam:
        log.warning(f'Período {periodo}: sem mosaico em {", ".join(faltam)} — células puladas.')
    return imgs


# ═══════════════════════════════════════════════════════════════════════════
# ETAPAS
# ═══════════════════════════════════════════════════════════════════════════
def etapa1(grid):
    if existe(ASSET_A):
        raise SystemExit(f'{ASSET_A} já existe — apague antes para refazer.')
    partes = []
    for cid, img in mosaicos('A', grid).items():
        cel = grid[cid]
        # stratifiedSample numa banda constante só onde todas as bandas têm dado:
        # sorteia N_POR_CELULA pontos COM dado, mesmo em célula pouco coberta
        valido = (ee.Image(1).rename('valido')
                    .updateMask(img.mask().reduce(ee.Reducer.min())).toInt())
        pts = (img.addBands(valido)
                  .stratifiedSample(numPoints=N_POR_CELULA, classBand='valido',
                                    region=cel['geom'].buffer(-RECUO_M), scale=ESCALA,
                                    seed=SEMENTE + cid, geometries=True, dropNulls=True, tileScale=4)
                  .select(BANDAS['A']))
        partes.append(pts.map(lambda f, cid=cid, nome=cel['nome']: f.set({
            'pid': ee.String(nome).cat('_').cat(f.id()), 'celula_id': cid, 'celula_nome': nome})))
    if not partes:
        raise SystemExit('Nenhum mosaico do período A exportado ainda.')
    fc = ee.FeatureCollection(partes).flatten()
    nome = ASSET_A.rsplit('/', 1)[-1]
    task = ee.batch.Export.table.toAsset(collection=fc, description=nome, assetId=ASSET_A)
    task.start()
    log.info(f'Etapa 1: task {task.id} → {ASSET_A} (até {N_POR_CELULA} pontos × {len(partes)} células)')


def etapa2(grid):
    if not existe(ASSET_A):
        t = tarefa(ASSET_A.rsplit('/', 1)[-1])
        raise SystemExit(f'{ASSET_A} ainda não existe'
                         + (f' (task da etapa 1: {t["state"]})' if t else '') + ' — rode a etapa 2 depois.')
    if existe(ASSET_B):
        raise SystemExit(f'{ASSET_B} já existe — apague antes para refazer.')
    pontos = ee.FeatureCollection(ASSET_A).select(PROPS_PONTO)
    partes = [img.sampleRegions(collection=pontos.filter(ee.Filter.eq('celula_id', cid)),
                                properties=PROPS_PONTO, scale=ESCALA,
                                geometries=True, tileScale=4)
              for cid, img in mosaicos('B', grid).items()]
    if not partes:
        raise SystemExit('Nenhum mosaico do período B exportado ainda.')
    fc = ee.FeatureCollection(partes).flatten()
    nome = ASSET_B.rsplit('/', 1)[-1]
    task = ee.batch.Export.table.toAsset(collection=fc, description=nome, assetId=ASSET_B)
    task.start()
    log.info(f'Etapa 2: task {task.id} → {ASSET_B} (mesmos pontos de {ASSET_A})')


def status():
    for asset in TABELAS.values():
        nome = asset.rsplit('/', 1)[-1]
        if existe(asset):
            por_cel = ee.FeatureCollection(asset).aggregate_histogram('celula_nome').getInfo()
            print(f'{nome}: {sum(por_cel.values())} pontos em {len(por_cel)} células')
            print('   ', ', '.join(f'{k} {v}' for k, v in sorted(por_cel.items())))
        else:
            t = tarefa(nome)
            print(f'{nome}: não existe' + (f' — task {t["state"]} {t.get("error_message", "")}' if t else ''))


# ═══════════════════════════════════════════════════════════════════════════
# BAIXAR E PAREAR (como baixar_pontos_pareados.py)
# ═══════════════════════════════════════════════════════════════════════════
def baixar_tabela(periodo: str) -> pd.DataFrame | None:
    asset = TABELAS[periodo]
    nome = asset.rsplit('/', 1)[-1]
    if not existe(asset):
        log.warning(f'{asset} não existe — pulado.')
        return None
    pasta = SAIDA / 'celulas' / nome
    pasta.mkdir(parents=True, exist_ok=True)
    fc = ee.FeatureCollection(asset)
    celulas = sorted(fc.aggregate_histogram('celula_nome').getInfo())
    colunas = CHAVES + BANDAS[periodo]

    for i, cel in enumerate(celulas, 1):
        arq = pasta / f'{cel}.csv'
        if arq.exists():
            continue
        # getDownloadURL (CSV): o computeFeatures paginado repetia pontos entre páginas
        sel = fc.filter(ee.Filter.eq('celula_nome', cel)).map(lambda f: f.set({
            'longitude': f.geometry().coordinates().get(0),
            'latitude': f.geometry().coordinates().get(1)}))
        df = pd.read_csv(sel.getDownloadURL(filetype='csv', selectors=colunas))
        df[colunas].to_csv(arq, index=False)
        log.info(f'{nome} [{i}/{len(celulas)}] {cel}: {len(df)} pontos')

    df = pd.concat([pd.read_csv(pasta / f'{c}.csv') for c in celulas], ignore_index=True)
    arq = SAIDA / f'points_sazonal_seco_{periodo}.csv'   # mesmo nome em todas as fontes
    df.to_csv(arq, index=False)
    log.info(f'{nome}: {len(df)} pontos em {len(celulas)} células → {arq}')
    return df


def baixar():
    SAIDA.mkdir(parents=True, exist_ok=True)
    tabelas = {per: baixar_tabela(per) for per in TABELAS}
    if any(df is None for df in tabelas.values()):
        log.info('Falta uma das tabelas — junção não feita.')
        return
    a, d = (df.rename(columns={c: f'{c}_{per}' for c in df.columns if c not in CHAVES})
            for per, df in tabelas.items())
    # longitude/latitude vêm da etapa 1; a etapa 2 amostra no mesmo ponto
    df = a.merge(d.drop(columns=['longitude', 'latitude', 'celula_id', 'celula_nome']), on='pid')
    arq = SAIDA / 'points_sazonal_seco_pareados.csv'
    df.to_csv(arq, index=False)
    log.info(f'Pareados: {len(df)} de {len(a)} pontos do A têm dado no B → {arq}')
    print(df.groupby('celula_nome').size().to_string())


# ═══════════════════════════════════════════════════════════════════════════
def main():
    switch_user(CONTA)
    projeto = get_project_from_account(CONTA)
    ee.Initialize(project=projeto)
    log.info(f'Earth Engine inicializado com {CONTA} (projeto {projeto}).')
    log.info(f'Fonte {FONTE}: {COL_FONTE} → pontos {ASSET_A} / _B, pasta {SAIDA}')
    acao = perguntar_opcao('O que fazer?', [
        (f'etapa 1 — pontos no mosaico seco A (2000–2008) → {ASSET_A}', 'etapa1'),
        (f'etapa 2 — mesmos pontos no mosaico seco B (2009–2020) → {ASSET_B}', 'etapa2'),
        ('status dos dois assets', 'status'),
        (f'baixar as duas tabelas e parear → {SAIDA}', 'baixar'),
    ])
    if acao == 'status':
        return status()
    if acao == 'baixar':
        return baixar()
    grid = ler_grid()
    (etapa1 if acao == 'etapa1' else etapa2)(grid)


if __name__ == '__main__':
    main()
