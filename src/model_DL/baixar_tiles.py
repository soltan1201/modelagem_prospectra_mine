"""
Baixa os canais de um modelo da CNN para o mosaico ASTER v1 inteiro, em
pedaços GeoTIFF, para a predição rodar depois sem GEE (inferir_cnn.py --tiles).

Cada célula da grid 5×5 é dividida em N × N pedaços (padrão 6 × 6 → ~19 × 20 km,
~640 × 660 px de 30 m). Cada pedaço vem com BORDA px a mais de cada lado (padrão
32 = meio patch de 64 px), para o pixel da beira ter o patch inteiro em volta —
sem emenda entre pedaços no mapa.

Canais: os mesmos do treino (exportar_patches.imagem_canais + MODELOS), uma
rodada por modelo:
  --modelo A   16 canais (15 índices + swir_estimado)
  --modelo B1  15 canais (9 log-bandas + 5 emissividades + swir_estimado)
Download por ee.data.computePixels (sem task; limite 48 MB por pedido — um
pedaço de ~720 × 720 px × 16 canais float32 tem ~33 MB).

Saída (TILES_DIR, padrão /srv/almacen/db_images/tiles_cnn)
  canais_<modelo>/<LxCy>/<LxCy>_r<i>c<j>.tif   GeoTIFF float32, EPSG:4326, grade de
      30 m do mosaico, NaN sem dado, deflate; bandas com o nome do canal; tags:
      modelo, celula, linha, coluna, borda, linha0/coluna0 (início do miolo na
      grade global)
  canais_<modelo>/indice_tiles.csv
Retomável: pedaço já gravado é pulado.

Espaço em disco antes de rodar:  df -h /srv/almacen
  (área toda ≈ 24 GB por modelo sem compressão; com deflate, menos)

Uso
  python baixar_tiles.py --modelo A
  python baixar_tiles.py --modelo B1 --celulas L2C3,L3C3
"""

import argparse
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import ee
import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
sys.path.append(str(AQUI))
import exportar_patches as ep   # noqa: E402  (canais, grade, grid, computePixels)

TILES_DIR = Path(os.environ.get('TILES_DIR', '/srv/almacen/db_images/tiles_cnn')).expanduser()
LIMITE_PEDIDO = 48 * 1024 ** 2      # bytes por computePixels


def limites_celula(c):
    """Colunas/linhas da grade de 30 m do mosaico que cobrem a célula [c0, c1) × [l0, l1)."""
    c0 = math.floor((c['oeste'] - ep.X0) / ep.RES)
    c1 = math.ceil((c['leste'] - ep.X0) / ep.RES)
    l0 = math.floor((ep.Y0 - c['norte']) / ep.RES)
    l1 = math.ceil((ep.Y0 - c['sul']) / ep.RES)
    return c0, c1, l0, l1


def pedacos(cel, n):
    """(i, j, linha0, linha1, coluna0, coluna1) do miolo de cada pedaço, na grade global."""
    c0, c1, l0, l1 = limites_celula(cel)
    cols = np.linspace(c0, c1, n + 1).round().astype(int)
    lins = np.linspace(l0, l1, n + 1).round().astype(int)
    return [(i, j, lins[i], lins[i + 1], cols[j], cols[j + 1]) for i in range(n) for j in range(n)]


def baixar(img, canais, lin0, col0, altura, largura):
    grade = {'dimensions': {'width': largura, 'height': altura},
             'affineTransform': {'scaleX': ep.RES, 'shearX': 0, 'translateX': ep.X0 + col0 * ep.RES,
                                 'shearY': 0, 'scaleY': -ep.RES, 'translateY': ep.Y0 - lin0 * ep.RES},
             'crsCode': 'EPSG:4326'}
    for t in range(ep.TENTATIVAS):
        try:
            arr = ee.data.computePixels({'expression': img, 'fileFormat': 'NUMPY_NDARRAY', 'grid': grade})
            x = np.stack([arr[c] for c in canais]).astype(np.float32)
            x[x == ep.NODATA] = np.nan
            return x
        except ee.EEException:
            if t == ep.TENTATIVAS - 1:
                raise
            time.sleep(2 ** t + np.random.rand())


def gravar_tif(arquivo, x, lin0, col0, canais, tags):
    """GeoTIFF float32 multibanda; (lin0, col0) = canto NW (com borda) na grade global."""
    import rasterio
    from rasterio.transform import from_origin
    tmp = arquivo.with_suffix('.tmp.tif')
    with rasterio.open(tmp, 'w', driver='GTiff', width=x.shape[2], height=x.shape[1], count=x.shape[0],
                       dtype='float32', crs='EPSG:4326', nodata=np.nan, compress='deflate', predictor=3,
                       tiled=True, blockxsize=256, blockysize=256,
                       transform=from_origin(ep.X0 + col0 * ep.RES, ep.Y0 - lin0 * ep.RES, ep.RES, ep.RES)) as d:
        d.write(x)
        for b, nome in enumerate(canais, 1):
            d.set_band_description(b, nome)
        d.update_tags(**{k: str(v) for k, v in tags.items()})
    tmp.replace(arquivo)


def main(lista=None):
    a = argparse.ArgumentParser(description='Canais do mosaico em pedaços GeoTIFF para a predição')
    a.add_argument('--modelo', required=True, choices=list(ep.MODELOS))
    a.add_argument('--n', type=int, default=6, help='pedaços por lado em cada célula')
    a.add_argument('--borda', type=int, default=ep.PATCH // 2, help='px extras de cada lado')
    a.add_argument('--celulas', default='', help='ex.: L2C3,L3C3 (vazio = todas)')
    a.add_argument('--threads', type=int, default=4)
    a.add_argument('--dir', default=str(TILES_DIR))
    a.add_argument('--projeto', default=ep.PROJETO_EE)
    args = a.parse_args(lista)

    ee.Initialize(project=args.projeto)
    canais = ep.MODELOS[args.modelo]
    img = ep.imagem_canais().select(canais)
    grid = ep.ler_grid()
    nomes = [c.strip().upper() for c in args.celulas.split(',') if c.strip()] or sorted(grid)
    raiz = Path(args.dir) / f'canais_{args.modelo}'
    raiz.mkdir(parents=True, exist_ok=True)

    tarefas = []
    for nome in nomes:
        (raiz / nome).mkdir(exist_ok=True)
        for i, j, l0, l1, c0, c1 in pedacos(grid[nome], args.n):
            arq = raiz / nome / f'{nome}_r{i}c{j}.tif'
            tarefas.append({'celula': nome, 'linha': i, 'coluna': j, 'arquivo': str(arq),
                            'lin0': l0 - args.borda, 'col0': c0 - args.borda,
                            'altura': l1 - l0 + 2 * args.borda, 'largura': c1 - c0 + 2 * args.borda,
                            'miolo_lin0': l0, 'miolo_col0': c0, 'miolo_lin1': l1, 'miolo_col1': c1})
    maior = max(t['altura'] * t['largura'] for t in tarefas) * len(canais) * 4
    print(f'modelo {args.modelo}: {len(canais)} canais · {len(nomes)} células × {args.n}×{args.n} = '
          f'{len(tarefas)} pedaços · maior pedido {maior / 1e6:.0f} MB · ≈ '
          f'{sum(t["altura"] * t["largura"] for t in tarefas) * len(canais) * 4 / 1e9:.1f} GB sem compressão')
    if maior > LIMITE_PEDIDO:
        raise SystemExit(f'Pedaço passa de 48 MB por pedido: aumente --n (hoje {args.n}).')
    pd.DataFrame(tarefas).to_csv(raiz / 'indice_tiles.csv', index=False)

    faltam = [t for t in tarefas if not Path(t['arquivo']).exists()]
    print(f'{len(tarefas) - len(faltam)} já baixados; faltam {len(faltam)} → {raiz}')

    def um(t):
        x = baixar(img, canais, t['lin0'], t['col0'], t['altura'], t['largura'])
        gravar_tif(Path(t['arquivo']), x, t['lin0'], t['col0'], canais,
                   {'modelo': args.modelo, 'celula': t['celula'], 'linha': t['linha'], 'coluna': t['coluna'],
                    'borda': args.borda, 'miolo_lin0': t['miolo_lin0'], 'miolo_col0': t['miolo_col0'],
                    'canais': ','.join(canais)})
        return t, float(np.isfinite(x[0]).mean())

    t0, feitos = time.time(), 0
    with ThreadPoolExecutor(args.threads) as ex:
        for f in as_completed([ex.submit(um, t) for t in faltam]):
            try:
                t, valido = f.result()
                feitos += 1
                print(f'  {Path(t["arquivo"]).name}  ({100 * valido:.0f}% com dado)  '
                      f'[{feitos}/{len(faltam)}, {time.time() - t0:.0f} s]')
            except Exception as e:   # pedaço perdido não para o resto; rode de novo para refazer
                print(f'  ERRO: {e}')
    print(f'Concluído: {feitos}/{len(faltam)} pedaços em {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
