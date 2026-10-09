"""
Download das FeatureCollections de pontos amostrados (asset GEE → CSV local).

Lê todas as tabelas da pasta ASSET_POINTS (geradas por
sample_quality_mosaic.py) e salva um CSV por asset em PASTA_SAIDA —
mesma pasta que train_swir_models.py lê (data/samples/*.csv).

Para cada asset:
  1. Tenta download direto via getDownloadURL (rápido, um request)
  2. Se falhar (tabela grande demais p/ URL), baixa em blocos de
     TAMANHO_BLOCO features via getInfo e grava o CSV localmente
  3. Pula assets cujo CSV já existe (SOBRESCREVER = False)
"""
import os
import csv
import logging
from pathlib import Path
import sys
import ee
import requests

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# ------- PARÂMETROS -------
ASSET_POINTS   = 'projects/mapbiomas-arida/mine/points'
PASTA_SAIDA    = Path('data/samples')
SOBRESCREVER   = False   # True → baixa de novo mesmo se o CSV já existir
INCLUIR_GEO    = False   # True → mantém a coluna .geo (GeoJSON do ponto)
TAMANHO_BLOCO  = 5000    # features por getInfo no modo fallback
TIMEOUT_HTTP   = 600     # segundos

# ------- INICIALIZAÇÃO -------
from pathlib import Path
pathparent = str(Path(os.getcwd()).parents[0])
sys.path.append(pathparent)
from configure_account_projects_ee import get_current_account, get_project_from_account
projAccount = get_current_account()
print(f"projetos selecionado >>> {projAccount} <<<")
try:
    ee.Initialize(project=projAccount)
    log.info('Earth Engine inicializado com sucesso.')
except Exception as e:
    log.error(f'Erro de inicialização: {e}')
    raise

PASTA_SAIDA.mkdir(parents=True, exist_ok=True)


# ------- FUNÇÕES -------
def listar_tabelas(pasta: str) -> list[str]:
    """Lista os IDs de todas as FeatureCollections da pasta (com paginação)."""
    ids, token = [], None
    while True:
        params = {'parent': pasta}
        if token:
            params['pageToken'] = token
        resp = ee.data.listAssets(params)
        ids += [a['id'] for a in resp.get('assets', []) if a.get('type') == 'TABLE']
        token = resp.get('nextPageToken')
        if not token:
            return ids


def colunas_da_tabela(fc: ee.FeatureCollection) -> list[str]:
    """Nomes das propriedades (união sobre a coleção, ordem estável)."""
    nomes = fc.first().propertyNames().getInfo()
    return sorted(n for n in nomes if not n.startswith('system:'))


def baixar_por_url(fc: ee.FeatureCollection, colunas: list[str], destino: Path) -> None:
    seletores = colunas + (['.geo'] if INCLUIR_GEO else [])
    url = fc.getDownloadURL(filetype='csv', selectors=seletores)
    with requests.get(url, stream=True, timeout=TIMEOUT_HTTP) as r:
        r.raise_for_status()
        with open(destino, 'wb') as f:
            for pedaco in r.iter_content(chunk_size=1 << 20):
                f.write(pedaco)


def baixar_por_blocos(fc: ee.FeatureCollection, colunas: list[str], destino: Path) -> None:
    n = fc.size().getInfo()
    cabecalho = colunas + (['.geo'] if INCLUIR_GEO else [])
    lista = fc.toList(n)
    with open(destino, 'w', newline='') as f:
        escritor = csv.DictWriter(f, fieldnames=cabecalho, extrasaction='ignore')
        escritor.writeheader()
        for ini in range(0, n, TAMANHO_BLOCO):
            bloco = ee.List(lista.slice(ini, min(ini + TAMANHO_BLOCO, n))).getInfo()
            for feat in bloco:
                linha = dict(feat.get('properties', {}))
                if INCLUIR_GEO:
                    linha['.geo'] = feat.get('geometry')
                escritor.writerow(linha)
            log.info(f'    {min(ini + TAMANHO_BLOCO, n):,}/{n:,} features')


def baixar_tabela(asset_id: str) -> None:
    nome = asset_id.rsplit('/', 1)[-1]
    destino = PASTA_SAIDA / f'{nome}.csv'
    if destino.exists() and not SOBRESCREVER:
        log.info(f'  [skip] {destino} já existe')
        return

    fc = ee.FeatureCollection(asset_id)
    colunas = colunas_da_tabela(fc)
    tmp = destino.with_suffix('.csv.part')
    try:
        baixar_por_url(fc, colunas, tmp)
        modo = 'url'
    except Exception as exc:
        log.warning(f'  download direto falhou ({exc}); usando blocos via getInfo')
        baixar_por_blocos(fc, colunas, tmp)
        modo = 'blocos'
    tmp.replace(destino)
    log.info(f'  [ok:{modo}] {destino}')


# ------- EXECUÇÃO -------
tabelas = listar_tabelas(ASSET_POINTS)
log.info(f'Tabelas encontradas em {ASSET_POINTS}: {len(tabelas)}')

falhas = []
for i, asset_id in enumerate(tabelas, 1):
    log.info(f'[{i}/{len(tabelas)}] {asset_id}')
    try:
        baixar_tabela(asset_id)
    except Exception as exc:
        log.error(f'  falhou: {exc}')
        falhas.append(asset_id)

log.info(f'Concluído: {len(tabelas) - len(falhas)} ok, {len(falhas)} falhas.')
for a in falhas:
    log.info(f'  falha: {a}')
