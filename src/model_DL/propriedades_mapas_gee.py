"""
Grava as propriedades dos mapas de p da CNN subidos à mão no GEE (Code Editor), a
partir do propriedades_upload.csv gerado pelo juntar_mapas_gee.py: uma linha por
asset, colunas = propriedades (modelo, alvo, mosaico, rodada, passo, preenchido,
canais, patch, escala_p, celulas). Assets que ainda não existem são listados como
"falta subir" e pulados.

O upload do Code Editor chama a banda de 'b1' (não 'probability'): no GEE use
.select(['b1'], ['probability']).

Conta EE: a ativa do servidor (current_user.txt), no projeto dela em
configure_account_projects_ee; precisa escrever em projects/mapbiomas-arida/mine.

Uso
  python propriedades_mapas_gee.py                    # <RUNS_DIR>/upload_gee/propriedades_upload.csv
  python propriedades_mapas_gee.py --conferir         # só mostra o que faria
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import ee

sys.path.append(str(Path(__file__).resolve().parents[1]))
from configure_account_projects_ee import get_current_account   # noqa: E402

INTEIROS = {'passo', 'preenchido', 'patch', 'escala_p', 'celulas'}


def propriedades(linha):
    """Colunas do CSV menos arquivo/asset; números como int, vazias fora."""
    props = {}
    for k, v in linha.items():
        if k in ('arquivo', 'asset') or v in ('', 'None', None):
            continue
        props[k] = int(v) if k in INTEIROS else v
    return props


def main(lista=None):
    a = argparse.ArgumentParser(description='Grava as propriedades dos mapas da CNN subidos à mão no GEE')
    raiz = Path(os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn'))
    a.add_argument('--csv', default=str(raiz / 'upload_gee' / 'propriedades_upload.csv'))
    a.add_argument('--projeto', default='', help='projeto EE (padrão: o da conta ativa)')
    a.add_argument('--conferir', action='store_true', help='não grava; só mostra')
    args = a.parse_args(lista)

    with open(args.csv, newline='') as f:
        linhas = list(csv.DictReader(f))
    projeto = args.projeto or get_current_account()
    ee.Initialize(project=projeto)
    print(f'EE: projeto {projeto} | {len(linhas)} linhas em {args.csv}')

    feitos, faltam = 0, []
    for linha in linhas:
        asset, props = linha['asset'], propriedades(linha)
        try:
            ee.data.getAsset(asset)
        except ee.EEException:
            print(f'  {asset}: falta subir')
            faltam.append(linha['arquivo'])
            continue
        if args.conferir:
            print(f'  {asset}: {props}')
            continue
        ee.data.setAssetProperties(asset, props)
        print(f'  {asset}: {len(props)} propriedades gravadas')
        feitos += 1

    print(f'Gravados {feitos}; faltam subir {len(faltam)}' + (f': {", ".join(faltam)}' if faltam else ''))


if __name__ == '__main__':
    main()
