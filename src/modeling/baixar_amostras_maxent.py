#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Baixa as amostras do Maxent (presença + fundo) salvas como asset pelo
maxent_fosfato.py e grava em CSV para análise local.

  Entrada : projects/mapbiomas-arida/mine/maxent_amostras_<mosaico>_<modo>
            (um por mosaico; propriedades: bandas/índices, extras do mosaico,
            'classe' 1 = presença / 0 = fundo, 'alvo', 'CODIGO_OBJ', 'fonte',
            'papel', 'celula'; geometria = ponto)
  Saída   : src/Dados/maxent/<nome_do_asset>.csv, com colunas lon/lat

Usa a conta ativa do GEE (~/.config/earthengine/current_user.txt). Baixa por
célula com getDownloadURL (computeFeatures paginado duplicava linhas) e junta
tudo no fim; as células já baixadas ficam em cache e não são baixadas de novo.

Uso:
  python baixar_amostras_maxent.py                      # asset padrão
  python baixar_amostras_maxent.py <asset> [<asset>...]  # outros assets
'''
import io
import sys
from pathlib import Path

import ee
import pandas as pd
import requests

RAIZ = Path(__file__).resolve().parents[1]
sys.path.append(str(RAIZ))
from configure_account_projects_ee import get_current_account   # noqa: E402

ASSETS = ['projects/mapbiomas-arida/mine/maxent_amostras_ASTER_mosaico_unico_final_bandas_indices',
          'projects/mapbiomas-arida/mine/maxent_amostras_aster-mosaico-final_bandas_indices']
PASTA  = RAIZ / 'Dados' / 'maxent'


def baixar_csv(fc, colunas):
    url = fc.getDownloadURL(filetype='csv', selectors=colunas)
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))


def baixar(asset):
    nome = asset.rstrip('/').split('/')[-1]
    cache = PASTA / 'celulas' / nome
    cache.mkdir(parents=True, exist_ok=True)

    fc = ee.FeatureCollection(asset).map(lambda f: f.set({
        'lon': f.geometry().coordinates().get(0),
        'lat': f.geometry().coordinates().get(1)}))
    props = ee.Feature(fc.first()).propertyNames().getInfo()
    colunas = sorted(p for p in props if not p.startswith('system:'))
    celulas = sorted(fc.aggregate_array('celula').distinct().getInfo())
    print(f'\n{nome}: {fc.size().getInfo()} amostras, {len(celulas)} células, {len(colunas)} colunas')

    partes = []
    for c in celulas:
        arq = cache / f'celula_{c}.csv'
        if arq.exists():
            print(f'  célula {c}: cache')
        else:
            df = baixar_csv(fc.filter(ee.Filter.eq('celula', c)), colunas)
            df.to_csv(arq, index=False)
            print(f'  célula {c}: {len(df)} linhas')
        partes.append(pd.read_csv(arq))

    df = pd.concat(partes, ignore_index=True)
    saida = PASTA / f'{nome}.csv'
    df.to_csv(saida, index=False)
    resumo = df['classe'].value_counts().to_dict() if 'classe' in df else {}
    print(f'  → {saida} ({len(df)} linhas; classe: {resumo})')


def main():
    ee.Initialize(project=get_current_account())
    PASTA.mkdir(parents=True, exist_ok=True)
    for asset in sys.argv[1:] or ASSETS:
        baixar(asset)


if __name__ == '__main__':
    main()
