"""
Versão Drive do exportar_patches.py: cria a base de treino/validação/teste da
CNN do zero direto no Google Drive, por tasks do GEE (Export.table.toDrive) —
sem servidor, sem rclone e sem baixar nada na máquina que roda o script.

Sorteio: o MESMO do exportar_patches.py (importado, não copiado) — pontos de
fosfato (ocorrências + jitter) e pseudo-rótulos sorteados nos mapas RF e GTB
(favorabilidade_ml), divisão por célula (teste + 5 folds), descarte perto de
borda de outro grupo. Regras, canais e parâmetros: docstring do exportar_patches.py.

Patches: os mesmos 64 × 64 px de 30 m do computePixels (grade alinhada ao
mosaico, centro no pixel do ponto, todos os CANAIS em float32, sem dado =
NODATA), montados no GEE com neighborhoodToArray (um reduceRegion por centro) e exportados em TFRecord com
todos os metadados da amostra. As tasks vão para o Drive da conta autenticada
no Earth Engine (earthengine authenticate), na pasta --pasta-drive:

  <pasta>/patches_<LxCy>_<k>.tfrecord.gz   lote k da célula (padrão: célula inteira numa
                                           task, ~25 tasks; --lote N divide em N centros)
  <pasta>/indice_patches.csv               todas as amostras (sem X)
  <pasta>/divisao.csv                      divisao.json numa coluna 'json'

Depois, no Colab: tfrecord_para_npz.py converte para o formato do treino
(modelo_A/ e modelo_B1/ patches_<LxCy>.npz, indice_patches.csv, divisao.json).

Retomável: tasks com a mesma descrição já pendentes, rodando ou concluídas são
puladas (rode de novo depois de uma falha para reenviar só as que faltam).

Uso
  python exportar_patches_drive.py --projeto ee-solkancengine17
  python exportar_patches_drive.py --projeto ee-solkancengine17 --celulas L2C3,L3C3
  python exportar_patches_drive.py --so-sortear          (mostra as quantidades e sai)
"""

import argparse
import json
import math
import sys
from pathlib import Path

import ee
import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
sys.path.append(str(AQUI))
import exportar_patches as ep   # noqa: E402  (sorteio, canais, grade)

PASTA_DRIVE = 'patches_cnn_gee'     # pasta na raiz do Meu Drive (o GEE não cria subpastas)
PREFIXO     = 'cnn_patches'         # descrição das tasks: <PREFIXO>_<LxCy>_<k>
LOTE        = 0                     # centros por task; 0 = célula inteira (~430) numa task
TILE_SCALE  = 16


def imagem_vizinhanca():
    """Canais na grade de 30 m do mosaico e, em cada pixel, a janela PATCH × PATCH em volta
    (foco em PATCH // 2, como o grade_patch do computePixels: colunas col-32 … col+31)."""
    grade = [ep.RES, 0, ep.X0, 0, -ep.RES, ep.Y0]
    img = ep.imagem_canais().reproject(crs='EPSG:4326', crsTransform=grade)
    pesos = ee.List.repeat(ee.List.repeat(1, ep.PATCH), ep.PATCH)
    kernel = ee.Kernel.fixed(ep.PATCH, ep.PATCH, pesos, ep.PATCH // 2, ep.PATCH // 2)
    return img.neighborhoodToArray(kernel, defaultValue=ep.NODATA), grade


def amostrar(vizinhanca, grade, fc):
    """Janela de cada centro com um reduceRegion por feição. O sampleRegions junta todos os
    pontos de um tile num valor só e estoura ('Computed value is too large'): cada patch
    tem 31 × 64 × 64 valores."""
    def um(f):
        v = vizinhanca.reduceRegion(reducer=ee.Reducer.first(), geometry=f.geometry(), crs='EPSG:4326',
                                    crsTransform=grade, tileScale=TILE_SCALE)
        return f.set(v)
    return fc.map(um)


def para_ee(v):
    """Valor do DataFrame → tipo aceito nas propriedades do GEE (NaN vira -1)."""
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return -1.0 if math.isnan(v) else float(v)
    return v if isinstance(v, (int, str)) else str(v)


def colecao(df):
    return ee.FeatureCollection([ee.Feature(ee.Geometry.Point([r['lon'], r['lat']]),
                                            {k: para_ee(v) for k, v in r.items()})
                                 for r in df.to_dict('records')])


def tasks_existentes():
    """Descrições de tasks pendentes, rodando ou concluídas (não reenviar)."""
    return {t['description'] for t in ee.data.getTaskList()
            if t.get('state') in ('READY', 'RUNNING', 'COMPLETED')}


def main(lista=None):
    a = argparse.ArgumentParser(description='Base de patches da CNN direto no Google Drive (tasks do GEE)')
    a.add_argument('--projeto', default=ep.PROJETO_EE)
    a.add_argument('--pasta-drive', default=PASTA_DRIVE, help='pasta na raiz do Meu Drive')
    a.add_argument('--celulas', default='', help='ex.: L2C3,L3C3 (vazio = todas)')
    a.add_argument('--lote', type=int, default=LOTE,
                   help='centros por task (0 = célula inteira; use 100–200 se a task estourar memória)')
    a.add_argument('--so-sortear', action='store_true', help='só sorteia e mostra as quantidades')
    args = a.parse_args(lista)

    ee.Initialize(project=args.projeto)
    papeis = (ee.ImageCollection(ep.ASSET_FAV).filter(ee.Filter.eq('alvo', ep.ALVO))
              .aggregate_array('papeis').distinct().getInfo())
    if 'None' in papeis:
        ep.log.warning(f'{ep.ASSET_FAV} foi treinado com todos os pontos (papeis = None): o teste não é '
                       'independente. Ver docstring do exportar_patches.py.')

    grid = ep.ler_grid()
    df, divisao = ep.montar_amostras(grid, ep.imagem_p())
    print('\n' + pd.crosstab(df.origem_amostra, df.conjunto, margins=True).to_string())
    print(f'{len(df)} patches × {len(ep.CANAIS)} canais × {ep.PATCH} × {ep.PATCH} ≈ '
          f'{len(df) * len(ep.CANAIS) * ep.PATCH ** 2 * 4 / 1e9:.1f} GB em TFRecord (float32, antes do gzip)')
    if args.so_sortear:
        return

    pedidas = [c.strip().upper() for c in args.celulas.split(',') if c.strip()]
    df_exp = df[df.celula.isin(pedidas)] if pedidas else df
    feitas = tasks_existentes()
    pendentes = []
    for cel, g in df_exp.groupby('celula'):
        g = g.sort_values('id')
        lote = args.lote or len(g)
        for k, i in enumerate(range(0, len(g), lote)):
            nome = f'{PREFIXO}_{cel}_{k}'
            if nome not in feitas:
                pendentes.append((nome, cel, k, g.iloc[i:i + lote]))
    extras = [n for n in (f'{PREFIXO}_indice', f'{PREFIXO}_divisao') if n not in feitas]
    print(f'\n{len(pendentes)} tasks de patches + {len(extras)} de índice/divisão → '
          f'Drive da conta do EE, pasta "{args.pasta_drive}" (projeto {args.projeto})')
    if not pendentes and not extras:
        print('Nada a enviar: tudo já pendente, rodando ou concluído.')
        return
    if input('Enviar as tasks? [s/N]: ').strip().lower() != 's':
        return

    if f'{PREFIXO}_indice' in extras:
        ee.batch.Export.table.toDrive(collection=colecao(df), description=f'{PREFIXO}_indice',
                                      folder=args.pasta_drive, fileNamePrefix='indice_patches',
                                      fileFormat='CSV', selectors=list(df.columns)).start()
    if f'{PREFIXO}_divisao' in extras:
        fc = ee.FeatureCollection([ee.Feature(None, {'json': json.dumps(divisao, ensure_ascii=False)})])
        ee.batch.Export.table.toDrive(collection=fc, description=f'{PREFIXO}_divisao',
                                      folder=args.pasta_drive, fileNamePrefix='divisao',
                                      fileFormat='CSV', selectors=['json']).start()

    vizinhanca, grade = imagem_vizinhanca()
    colunas = list(df.columns)
    for nome, cel, k, g in pendentes:
        amostras = amostrar(vizinhanca, grade, colecao(g))
        ee.batch.Export.table.toDrive(collection=amostras, description=nome, folder=args.pasta_drive,
                                      fileNamePrefix=f'patches_{cel}_{k}', fileFormat='TFRecord',
                                      selectors=ep.CANAIS + colunas).start()
        print(f'  {nome}: {len(g)} centros')
    print(f'\n{len(pendentes)} tasks enviadas. Acompanhe em https://code.earthengine.google.com/tasks; '
          f'no fim, rode tfrecord_para_npz.py no Colab.')


if __name__ == '__main__':
    main()
