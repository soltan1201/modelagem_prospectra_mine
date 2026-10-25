"""
Corta os GeoTIFF por célula de exportar_tiles_drive.py (tasks do GEE no Drive)
nos pedaços do baixar_tiles.py — mesma divisão N × N por célula, mesma borda,
nomes de banda e tags —, que é o que inferir_cnn.py --tiles lê:

  <destino>/canais_<modelo>/<LxCy>/<LxCy>_r<i>c<j>.tif
  <destino>/canais_<modelo>/indice_tiles.csv

A posição de cada arquivo na grade global sai da geotransformação dele (aceita
os vários arquivos que o GEE cria quando divide uma célula grande). Sem GEE.
Retomável: pedaço já gravado é pulado.

Uso (Colab, Drive montado)
  !python cortar_tiles_drive.py --modelo A --origem /content/drive/MyDrive/tiles_cnn_gee \
                                --destino /content/drive/MyDrive/prospectra/tiles_cnn
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
sys.path.append(str(AQUI))
import exportar_patches as ep   # noqa: E402  (canais, grade)
import baixar_tiles as bt       # noqa: E402  (pedaços, gravar_tif)


def posicao(fonte):
    """(linha0, coluna0) do canto NW do arquivo na grade global de 30 m."""
    tr = fonte.transform
    return round((ep.Y0 - tr.f) / ep.RES), round((tr.c - ep.X0) / ep.RES)


def ler_janela(fontes, lin0, col0, altura, largura, nb):
    """Janela na grade global montada a partir dos arquivos da célula (NaN onde nenhum cobre)."""
    from rasterio.windows import Window
    x = np.full((nb, altura, largura), np.nan, np.float32)
    for f in fontes:
        fl, fc = posicao(f)
        r0, r1 = max(lin0, fl), min(lin0 + altura, fl + f.height)
        c0, c1 = max(col0, fc), min(col0 + largura, fc + f.width)
        if r0 >= r1 or c0 >= c1:
            continue
        x[:, r0 - lin0:r1 - lin0, c0 - col0:c1 - col0] = f.read(
            window=Window(c0 - fc, r0 - fl, c1 - c0, r1 - r0)).astype(np.float32)
    x[x == ep.NODATA] = np.nan
    return x


def main(lista=None):
    import rasterio
    a = argparse.ArgumentParser(description='GeoTIFF por célula (GEE → Drive) → pedaços do inferir_cnn --tiles')
    a.add_argument('--modelo', required=True, choices=list(ep.MODELOS))
    a.add_argument('--origem', default='/content/drive/MyDrive/tiles_cnn_gee')
    a.add_argument('--destino', default=os.environ.get('TILES_DIR', '/content/drive/MyDrive/prospectra/tiles_cnn'))
    a.add_argument('--n', type=int, default=6, help='pedaços por lado em cada célula')
    a.add_argument('--borda', type=int, default=ep.PATCH // 2, help='a mesma do exportar_tiles_drive.py')
    a.add_argument('--celulas', default='', help='ex.: L2C3,L3C3 (vazio = todas)')
    args = a.parse_args(lista)

    canais = ep.MODELOS[args.modelo]
    por_celula = defaultdict(list)
    for arq in sorted(Path(args.origem).glob(f'canais_{args.modelo}_*.tif')):
        m = re.match(rf'canais_{args.modelo}_(L\d+C\d+)', arq.name)
        if m:
            por_celula[m.group(1)].append(arq)
    pedidas = [c.strip().upper() for c in args.celulas.split(',') if c.strip()]
    por_celula = {c: v for c, v in por_celula.items() if not pedidas or c in pedidas}
    if not por_celula:
        raise SystemExit(f'Nenhum canais_{args.modelo}_*.tif em {args.origem} — as tasks já terminaram?')

    raiz = Path(args.destino) / f'canais_{args.modelo}'
    b, linhas = args.borda, []
    for cel, arquivos in sorted(por_celula.items()):
        (raiz / cel).mkdir(parents=True, exist_ok=True)
        fontes = [rasterio.open(x) for x in arquivos]
        try:
            if fontes[0].count != len(canais):
                raise SystemExit(f'{arquivos[0].name}: {fontes[0].count} bandas, modelo {args.modelo} tem {len(canais)}')
            pos = [posicao(f) for f in fontes]
            l_min, c_min = min(p[0] for p in pos), min(p[1] for p in pos)
            l_max = max(p[0] + f.height for p, f in zip(pos, fontes))
            c_max = max(p[1] + f.width for p, f in zip(pos, fontes))
            limites = (c_min + b, c_max - b, l_min + b, l_max - b)      # miolo da célula (sem borda)
            feitos = 0
            for i, j, l0, l1, c0, c1 in bt.pedacos(None, args.n, limites):
                arq = raiz / cel / f'{cel}_r{i}c{j}.tif'
                t = {'celula': cel, 'linha': i, 'coluna': j, 'arquivo': str(arq),
                     'lin0': l0 - b, 'col0': c0 - b, 'altura': l1 - l0 + 2 * b, 'largura': c1 - c0 + 2 * b,
                     'miolo_lin0': l0, 'miolo_col0': c0, 'miolo_lin1': l1, 'miolo_col1': c1}
                linhas.append(t)
                if arq.exists():
                    continue
                x = ler_janela(fontes, t['lin0'], t['col0'], t['altura'], t['largura'], len(canais))
                bt.gravar_tif(arq, x, t['lin0'], t['col0'], canais,
                              {'modelo': args.modelo, 'celula': cel, 'linha': i, 'coluna': j, 'borda': b,
                               'miolo_lin0': l0, 'miolo_col0': c0, 'canais': ','.join(canais)})
                feitos += 1
        finally:
            for f in fontes:
                f.close()
        print(f'{cel}: {len(arquivos)} arquivo(s) · {c_max - c_min} × {l_max - l_min} px → '
              f'{feitos} pedaços novos de {args.n * args.n}')
    pd.DataFrame(linhas).to_csv(raiz / 'indice_tiles.csv', index=False)
    print(f'Concluído → {raiz}  (inferir_cnn.py --tiles {raiz})')


if __name__ == '__main__':
    main()
