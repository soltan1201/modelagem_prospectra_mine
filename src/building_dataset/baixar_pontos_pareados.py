"""
Baixa para o PC, em CSV, os pontos pareados dos mosaicos zero-cloud
(pontos_zero_cloud_pareados.py):

  ASSET_00_08 → points_zeroCC_00_08.csv   (VNIR + SWIR + TIR + data)
  ASSET_09_20 → points_zeroCC_09_20.csv   (VNIR + TIR + data)
  junção      → points_zeroCC_pareados.csv

Baixa célula por célula (cache: CSV de célula já baixado é pulado — apague a
pasta celulas/ para baixar de novo). A junção é pelo pid (inner join: só pontos
com dado nos dois períodos); as bandas ganham sufixo _00_08 / _09_20, e as
colunas pid, celula_id, celula_nome, longitude, latitude ficam sem sufixo.

Reflectância = valor / 1000; TIR em DN; data = AAAAMMDD da cena do pixel.

Saída: src/Dados/ASTER_Samples_0cloud/pareados/

Uso: python baixar_pontos_pareados.py
"""

import logging
import sys

import ee
import pandas as pd

import build_mosaicos_semestrais as b
import pontos_zero_cloud_pareados as p

log = logging.getLogger(__name__)

sys.path.append(str(b.PASTA.parent))
from configure_account_projects_ee import get_current_account   # noqa: E402

SAIDA = b.PASTA.parent / 'Dados' / 'ASTER_Samples_0cloud' / 'pareados'
TABELAS = {'2000_2008': p.ASSET_00_08, '2009_2020': p.ASSET_09_20}
SUFIXOS = {'2000_2008': '_00_08', '2009_2020': '_09_20'}
CHAVES = p.PROPS_PONTO + ['longitude', 'latitude']


def baixar(periodo: str) -> pd.DataFrame | None:
    asset = TABELAS[periodo]
    nome = asset.rsplit('/', 1)[-1]
    if not p.existe(asset):
        log.warning(f'{asset} não existe — pulado.')
        return None
    pasta = SAIDA / 'celulas' / nome
    pasta.mkdir(parents=True, exist_ok=True)
    fc = ee.FeatureCollection(asset)
    celulas = sorted(fc.aggregate_histogram('celula_nome').getInfo())
    colunas = CHAVES + p.BANDAS[periodo]

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
    arq = SAIDA / f'{nome}.csv'
    df.to_csv(arq, index=False)
    log.info(f'{nome}: {len(df)} pontos em {len(celulas)} células → {arq}')
    return df


def parear(tabelas: dict[str, pd.DataFrame]):
    renomeadas = [df.rename(columns={c: c + SUFIXOS[per] for c in df.columns if c not in CHAVES})
                  for per, df in tabelas.items()]
    a, d = renomeadas
    # longitude/latitude vêm da etapa 1; a etapa 2 amostra no mesmo ponto
    df = a.merge(d.drop(columns=['longitude', 'latitude', 'celula_id', 'celula_nome']), on='pid')
    arq = SAIDA / 'points_zeroCC_pareados.csv'
    df.to_csv(arq, index=False)
    log.info(f'Pareados: {len(df)} de {len(a)} pontos de 2000–2008 têm dado em 2009–2020 → {arq}')
    print(df.groupby('celula_nome').size().to_string())


def main():
    ee.Initialize(project=get_current_account())
    SAIDA.mkdir(parents=True, exist_ok=True)
    tabelas = {per: baixar(per) for per in TABELAS}
    if all(df is not None for df in tabelas.values()):
        parear(tabelas)
    else:
        log.info('Falta uma das tabelas — junção não feita.')


if __name__ == '__main__':
    main()
