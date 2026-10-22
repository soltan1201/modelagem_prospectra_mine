#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Amostras de PRESENÇA de um alvo mineral para o Maxent, baixadas direto em CSV.

Pontos: pontos_referencia (preparar_pontos_referencia.py), filtrado por 'alvo'
(menu: fosfato, carbonato, mvt, ouro, manganes) e, se PAPEIS não for None, por
'papel'.

O asset de amostras do teste (zero-cloud 2000–2008) saiu sem presença: esse
mosaico não tem nenhum pixel válido a até 500 m dos pontos de fosfato. Aqui
a presença é amostrada no mosaico com prioridade do build_mosaico_final.py,
montado na hora para a célula de cada ponto:

  origem 1  zero-cloud 2000–2008   (SWIR medido)
  origem 2  zero-cloud 2009–2020   (SWIR estimado RF)
  origem 3  semestrais v4 2000–2008 (SWIR medido)
  origem 4  semestrais v4 2009–2020 (SWIR estimado RF)

Mesma regra do maxent_fosfato.py: até AMOSTRAS_POR_PONTO pixels sorteados no
buffer de BUFFER_M m de cada ponto, escala 30 m, mesma semente. Bandas B01–B14
nas mesmas unidades das amostras de fundo (reflectância × 1000, TIR em DN),
mais 'origem' e 'swir_estimado' para saber de que camada veio cada pixel.

Saída (classe = 1, CODIGO_OBJ, alvo, fonte, papel, celula, lon, lat)
  CSV    src/Dados/maxent/maxent_<alvo>_presenca_mosaico_prioridade.csv
  asset  projects/mapbiomas-arida/mine/maxent_<alvo>_presenca_mosaico_prioridade
         (task toAsset com as mesmas amostras; pergunta antes de exportar)

Usa a conta ativa do GEE; a conta precisa ler as pastas ASTER_0cloud,
ASTER_MOSAIC e ASTER_v4 do mapbiomas-brazil (a pasta v4 do mapbiomas-arida é
usada se a conta tiver acesso).
'''
import io
import sys
from pathlib import Path

import ee
import pandas as pd
import requests

RAIZ = Path(__file__).resolve().parents[1]
sys.path.append(str(RAIZ))
sys.path.append(str(RAIZ / 'building_dataset'))
from configure_account_projects_ee import get_current_account   # noqa: E402
import build_mosaicos_semestrais as b                           # noqa: E402
import build_mosaico_final as mf                                # noqa: E402

ASSET_PONTOS       = 'projects/mapbiomas-arida/mine/pontos_referencia'
PASTA_ASSET        = 'projects/mapbiomas-arida/mine'
ALVOS              = ['fosfato', 'carbonato', 'mvt', 'ouro', 'manganes']
PAPEIS             = None    # None = todos; ['treino'] = sem os pontos de validação (R2 do fosfato)
BUFFER_M           = 500
AMOSTRAS_POR_PONTO = 20
ESCALA             = 30
SEED               = 42
PASTA_SAIDA        = RAIZ / 'Dados' / 'maxent'


def celula_dos_pontos(pontos, grid):
    """CODIGO_OBJ → id da primeira célula que contém o ponto (a grid tem solape)."""
    d = ee.Dictionary({str(cid): pontos.filterBounds(c['geom']).aggregate_array('CODIGO_OBJ')
                       for cid, c in grid.items()}).getInfo()
    cel = {}
    for cid in sorted(d, key=int):
        for cod in d[cid]:
            cel.setdefault(int(cod), int(cid))
    return cel


def main():
    ee.Initialize(project='mapbiomas-brazil')   # antes: get_current_account()
    try:
        ee.data.getAsset(ASSET_PONTOS)
    except ee.EEException:
        raise SystemExit(f'{ASSET_PONTOS} não existe — rode preparar_pontos_referencia.py e espere a task.')
    print('Alvo:\n' + '\n'.join(f'  {i}) {a}' for i, a in enumerate(ALVOS, 1)))
    alvo = ALVOS[int(input('Opção [1]: ').strip() or '1') - 1]
    saida = PASTA_SAIDA / f'maxent_{alvo}_presenca_mosaico_prioridade.csv'
    grid = b.ler_grid()
    pontos = ee.FeatureCollection(ASSET_PONTOS).filter(ee.Filter.eq('alvo', alvo))
    if PAPEIS:
        pontos = pontos.filter(ee.Filter.inList('papel', PAPEIS))
    cel_ponto = celula_dos_pontos(pontos, grid)
    print(f'{alvo}: {len(cel_ponto)} pontos dentro da grid, em {len(set(cel_ponto.values()))} células')

    existentes = set()
    for pasta in (mf.PASTA_ZERO_CLOUD, mf.PASTA_SWIR_RF):
        existentes |= {f'{pasta}/{n}' for n in b.listar_destino(pasta, criar=False)}
    v4 = ee.ImageCollection([ee.Image(i) for i in mf.ids_v4()])

    partes, fcs = [], []
    for cid in sorted(set(cel_ponto.values())):
        cel = {'id': cid, **grid[cid]}
        cam = mf.camadas(cel, v4, existentes)
        if cam[3][2].getInfo() == 0:
            del cam[3]
        mosaico = mf.montar(cel, cam)
        codigos = [k for k, v in cel_ponto.items() if v == cid]

        def por_ponto(f, cid=cid):
            props = {'classe': 1, 'CODIGO_OBJ': f.get('CODIGO_OBJ'), 'alvo': f.get('alvo'),
                     'fonte': f.get('fonte'), 'papel': f.get('papel'), 'celula': cid}
            return (mosaico.sample(region=f.geometry().buffer(BUFFER_M), scale=ESCALA,
                                   numPixels=AMOSTRAS_POR_PONTO * 5, seed=SEED, geometries=True)
                           .limit(AMOSTRAS_POR_PONTO)
                           .map(lambda s: s.set(props)))

        fc = (ee.FeatureCollection(pontos.filter(ee.Filter.inList('CODIGO_OBJ', codigos))
                                         .map(por_ponto)).flatten()
                .map(lambda f: f.set({'lon': f.geometry().coordinates().get(0),
                                      'lat': f.geometry().coordinates().get(1)})))
        colunas = mf.BANDAS + ['origem', 'swir_estimado', 'classe', 'CODIGO_OBJ', 'alvo', 'fonte', 'papel',
                               'celula', 'lon', 'lat']
        r = requests.get(fc.getDownloadURL(filetype='csv', selectors=colunas), timeout=600)
        r.raise_for_status()
        df = pd.read_csv(io.StringIO(r.text)) if r.text.strip().count('\n') else pd.DataFrame(columns=colunas)
        print(f'  célula {cid} ({cel["nome"]}): camadas {sorted(cam)}, {len(codigos)} pontos, '
              f'{len(df)} amostras' + (f', origem {df["origem"].value_counts().to_dict()}' if len(df) else ''))
        partes.append(df)
        fcs.append(fc)

    df = pd.concat([p for p in partes if len(p)], ignore_index=True)
    saida.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(saida, index=False)
    sem = sorted(set(cel_ponto) - set(df['CODIGO_OBJ'].astype(int)))
    print(f'\n→ {saida}: {len(df)} amostras de {df["CODIGO_OBJ"].nunique()} pontos'
          + (f'; pontos sem pixel válido: {sem}' if sem else ''))

    asset = f'{PASTA_ASSET}/{saida.stem}'
    if input(f'\nExportar as amostras para {asset}? [S/n]: ').strip().lower() == 'n':
        return
    try:
        ee.data.getAsset(asset)
        raise SystemExit(f'{asset} já existe — apague ou renomeie antes.')
    except ee.EEException:
        pass
    task = ee.batch.Export.table.toAsset(collection=ee.FeatureCollection(fcs).flatten(),
                                         description=saida.stem, assetId=asset)
    task.start()
    print(f'task {task.id} → {asset}')


if __name__ == '__main__':
    main()
