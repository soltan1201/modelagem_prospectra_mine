"""
Junta as 25 células do mapa de p da CNN (p_<LxCy>.tif) num GeoTIFF só por rodada e
passo, comprimido, para subir à mão no GEE (Code Editor → Assets → NEW → GeoTIFF),
sem passar pelo bucket.

Entrada   <runs>/<rodada>/mapas_passo<k>[_preenchido]/p_<LxCy>_<rodada>.tif
          (a versão _preenchido, quando existe, no lugar da original)
Saída     <saida>/<rodada>_p<k>[_preenchido].tif   UInt16, p × 10000, 0 = sem dado,
          EPSG:4326, DEFLATE nível 9 + predictor 2, blocos 256
          <saida>/propriedades_upload.csv   asset e propriedades de cada arquivo,
          no mesmo formato do subir_predicao_gee.py, para preencher no upload

No upload: destino projects/mapbiomas-arida/mine/cnn_favorabilidade/<nome do arquivo>,
nodata 0, pyramiding MEAN; banda b1 = probability.

Uso
  python juntar_mapas_gee.py                          # todas as rodadas em RUNS_DIR
  python juntar_mapas_gee.py --run $RUNS/B1_regnety_64 --passo 3
  python juntar_mapas_gee.py --com-originais          # também as originais quando há _preenchido
"""

import argparse
import csv
import json
import os
import re
from pathlib import Path

COLECAO       = 'projects/mapbiomas-arida/mine/cnn_favorabilidade'
ALVO, MOSAICO = 'fosfato', 'aster-mosaico-final'


def listar(runs, passo, com_originais):
    """(run, passo, pasta) de cada mapas_passo<k>[_preenchido] com p_L*C*.tif."""
    alvos = []
    for run in runs:
        for pasta in sorted(run.glob('mapas_passo*')):
            m = re.fullmatch(r'mapas_passo(\d+)(_preenchido)?', pasta.name)
            if not (pasta.is_dir() and m and any(pasta.glob('p_L*C*.tif'))):
                continue
            k = int(m.group(1))
            if not m.group(2) and not com_originais and pasta.with_name(pasta.name + '_preenchido').is_dir():
                continue                                     # já vai a versão preenchida
            if passo is None or k == passo:
                alvos.append((run, k, pasta))
    return alvos


def juntar(arquivos, saida):
    """Mosaico das células num GeoTIFF só (UInt16, 0 = sem dado), comprimido para o upload."""
    import rasterio
    from rasterio.merge import merge
    fontes = [rasterio.open(a) for a in arquivos]
    try:
        mos, tr = merge(fontes, nodata=0)
        perfil = fontes[0].profile.copy()
    finally:
        for f in fontes:
            f.close()
    perfil.update(driver='GTiff', width=mos.shape[2], height=mos.shape[1], transform=tr, nodata=0,
                  dtype='uint16', compress='deflate', zlevel=9, predictor=2,
                  tiled=True, blockxsize=256, blockysize=256, BIGTIFF='IF_SAFER')
    with rasterio.open(saida, 'w', **perfil) as d:
        d.write(mos.astype('uint16'))
        d.set_band_description(1, 'probability')
        d.update_tags(escala='p x 10000', nodata='0')
    return mos.shape[1:]


def main(lista=None):
    a = argparse.ArgumentParser(description='Junta as células do mapa de p num GeoTIFF por rodada para o GEE')
    a.add_argument('--run', default='', help='pasta da rodada (padrão: todas em --runs)')
    a.add_argument('--runs', default=os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn'))
    a.add_argument('--passo', type=int, default=None, help='padrão: todos os mapas_passo<k>')
    a.add_argument('--saida', default='', help='padrão: <runs>/upload_gee')
    a.add_argument('--com-originais', action='store_true',
                   help='junta também mapas_passo<k> quando existe mapas_passo<k>_preenchido')
    a.add_argument('--sobrescrever', action='store_true')
    args = a.parse_args(lista)

    raiz = Path(args.runs).expanduser()
    runs = [Path(args.run)] if args.run else sorted(p for p in raiz.iterdir() if p.is_dir())
    alvos = listar(runs, args.passo, args.com_originais)
    if not alvos:
        raise SystemExit('Nenhuma pasta mapas_passo<k> com p_L*C*.tif — rode inferir_cnn.py --tiles antes')
    saida = Path(args.saida) if args.saida else raiz / 'upload_gee'
    saida.mkdir(parents=True, exist_ok=True)

    linhas = []
    for run, passo, pasta in alvos:
        cfg = json.loads((run / 'config.json').read_text()) if (run / 'config.json').exists() else {}
        preenchido = pasta.name.endswith('_preenchido')
        nome = f'{run.name}_p{passo}{"_preenchido" if preenchido else ""}'.replace('.', '_')
        destino = saida / f'{nome}.tif'
        arquivos = sorted(pasta.glob('p_L*C*.tif'))
        if destino.exists() and not args.sobrescrever:
            print(f'  {nome}: já existe — pulado')
        else:
            h, w = juntar(arquivos, destino)
            print(f'  {nome}: {len(arquivos)} células → {w}×{h} px, {destino.stat().st_size / 2**20:.1f} MB')
        linhas.append({'arquivo': destino.name, 'asset': f'{COLECAO}/{nome}',
                       'modelo': f'CNN_{cfg.get("backbone")}', 'alvo': ALVO, 'mosaico': MOSAICO,
                       'rodada': run.name, 'passo': passo, 'preenchido': int(preenchido),
                       'canais': cfg.get('modelo'), 'patch': cfg.get('patch'), 'escala_p': 10000,
                       'celulas': len(arquivos)})

    arq = saida / 'propriedades_upload.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    print(f'{len(linhas)} arquivos em {saida}\nPropriedades → {arq}')


if __name__ == '__main__':
    main()
