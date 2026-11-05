"""
Preenche as propriedades dos mapas de p da CNN subidos à mão no Code Editor (que entram
sem propriedades) em projects/mapbiomas-arida/mine/favorabilidade_ml — outra coleção com
--colecao, ex.: cnn_favorabilidade —, para o visualizar_cnn_favorabilidade.js agrupar e
descrever os modelos.

Para cada imagem da coleção:
  1. pelo NOME do asset, na convenção do juntar_mapas_gee.py / subir_predicao_gee.py
       <rodada>_p<k>[_preenchido][_<LxCy>]      ex.: B1_regnety_64_p3_preenchido
     e a rodada na convenção do rodadas_treino.md
       <A|B1>_<backbone curto>_<patch>[_shm][_teste][_1fold]   (_shm = mosaico sazonal)
     → rodada, passo, preenchido, celula, canais (A/B1), modelo (CNN_<backbone>), patch,
       folds (media5 = média dos 5 fold*.pt; 1fold = rodada irmã _1fold; final = final.pt
       da rodada _teste), alvo, mosaico, escala_p
  2. se existir o propriedades_upload.csv do juntar_mapas_gee.py, a linha do asset
     completa/substitui (backbone exato do config.json, nº de células)

Por padrão só mexe nas imagens SEM a propriedade 'rodada' (--todas regrava todas).
Nomes fora da convenção são listados e pulados; imagens de outros modelos (RF/GTB,
quando os mapas foram subidos no favorabilidade_ml: --colecao) nunca são tocadas.

O upload do Code Editor chama a banda de 'b1' (não 'probability'); o visualizador lê a
1ª banda, qualquer que seja o nome.

Conta EE: a ativa do servidor (current_user.txt), no projeto dela em
configure_account_projects_ee; precisa escrever em projects/mapbiomas-arida/mine.

Uso
  python predict/propriedades_mapas_gee.py --conferir         # só mostra o que faria
  python predict/propriedades_mapas_gee.py                    # grava nas imagens sem propriedades
  python predict/propriedades_mapas_gee.py --todas            # regrava todas
  python predict/propriedades_mapas_gee.py --colecao projects/mapbiomas-arida/mine/cnn_favorabilidade
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import ee

sys.path.append(str(Path(__file__).resolve().parents[3]))   # src/
from configure_account_projects_ee import get_current_account   # noqa: E402

COLECAO       = 'projects/mapbiomas-arida/mine/favorabilidade_ml'   # os mapas da CNN foram subidos aqui
ALVO          = 'fosfato'

# Rodadas do mosaico sazonal harmonizado levam _shm no nome (rodadas_treino.md):
# <modelo>_<backbone curto>_<patch>_shm[_teste][_1fold]. Sem _shm = mosaico v1.
MOSAICO_V1, MOSAICO_SHM = 'aster-mosaico-final', 'aster-mosaico-final_sazonal_hm_razao'


def mosaico_da_rodada(rodada):
    return MOSAICO_SHM if '_shm' in rodada else MOSAICO_V1
# nome curto da rodada → --backbone do treinar_cnn.py (rodadas_treino.md)
BACKBONES = {'pequena': 'pequena', 'regnety': 'regnety_016', 'effv2s': 'tf_efficientnetv2_s',
             'resnet50': 'resnet50', 'coatnet0': 'coatnet_0_rw_224'}
INTEIROS = {'passo', 'preenchido', 'patch', 'escala_p', 'celulas'}
RE_ASSET  = re.compile(r'(?P<rodada>.+?)_p(?P<passo>\d+)(?P<pre>_preenchido)?(?:_(?P<cel>L\d+C\d+))?')
RE_RODADA = re.compile(r'(?P<canais>A|B1)_(?P<curto>[a-z0-9]+)_(?P<patch>\d+)(?P<sufixo>(?:_shm|_teste|_1fold)*)')


def props_do_nome(nome):
    """Propriedades tiradas do nome do asset; None se fora da convenção."""
    m = RE_ASSET.fullmatch(nome)
    if not m:
        return None
    props = {'rodada': m['rodada'], 'passo': int(m['passo']), 'preenchido': int(bool(m['pre'])),
             'alvo': ALVO, 'mosaico': mosaico_da_rodada(m['rodada']), 'escala_p': 10000}
    if m['cel']:
        props['celula'] = m['cel']
    r = RE_RODADA.fullmatch(m['rodada'])
    if r:
        props.update(canais=r['canais'], patch=int(r['patch']),
                     modelo=f'CNN_{BACKBONES.get(r["curto"], r["curto"])}', folds=folds_da_rodada(r['sufixo']))
    return props


def folds_da_rodada(sufixo):
    """Modelos que fizeram o mapa, pelo sufixo da rodada: _teste → final.pt, _1fold → um fold
    (inferir_cnn.py --folds auto), sem sufixo → média dos 5 fold*.pt."""
    if '_teste' in sufixo:
        return 'final'
    return '1fold' if '_1fold' in sufixo else 'media5'


def ler_csv(arq):
    """{asset: propriedades} do propriedades_upload.csv (vazio se não existir)."""
    if not Path(arq).exists():
        return {}
    with open(arq, newline='') as f:
        linhas = list(csv.DictReader(f))
    saida = {}
    for linha in linhas:
        saida[linha['asset']] = {k: (int(v) if k in INTEIROS else v) for k, v in linha.items()
                                 if k not in ('arquivo', 'asset') and v not in ('', 'None', None)}
    return saida


def listar_imagens(colecao):
    """[{id, properties}] de todas as imagens da coleção (com paginação)."""
    imagens, token = [], None
    while True:
        pedido = {'parent': colecao, 'pageSize': 1000}
        if token:
            pedido['pageToken'] = token
        r = ee.data.listAssets(pedido)
        imagens += r.get('assets', [])
        token = r.get('nextPageToken')
        if not token:
            return imagens


def main(lista=None):
    a = argparse.ArgumentParser(description='Preenche as propriedades dos mapas da CNN na coleção do GEE')
    raiz = Path(os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn'))
    a.add_argument('--colecao', default=COLECAO)
    a.add_argument('--csv', default=str(raiz / 'upload_gee' / 'propriedades_upload.csv'),
                   help='propriedades_upload.csv do juntar_mapas_gee.py (opcional)')
    a.add_argument('--projeto', default='', help='projeto EE (padrão: o da conta ativa)')
    a.add_argument('--todas', action='store_true', help='regrava também as que já têm rodada')
    a.add_argument('--conferir', action='store_true', help='não grava; só mostra')
    args = a.parse_args(lista)

    projeto = args.projeto or get_current_account()
    ee.Initialize(project=projeto)
    do_csv = ler_csv(args.csv)
    imagens = listar_imagens(args.colecao)
    print(f'EE: projeto {projeto} | {len(imagens)} imagens em {args.colecao} | '
          f'CSV: {len(do_csv)} linhas' + ('' if do_csv else f' (sem {args.csv}: só pelo nome)'))

    gravados, pulados, fora = 0, 0, []
    for im in imagens:
        asset = im['id']
        nome = asset.split('/')[-1]
        tem = ee.data.getAsset(asset).get('properties', {})
        if 'modelo' in tem and not str(tem['modelo']).startswith('CNN_'):
            continue                                  # RF/GTB do favorabilidade_ml: não mexe
        if not args.todas and 'rodada' in tem:
            pulados += 1
            continue
        props = props_do_nome(nome)
        if props is None:
            fora.append(nome)
            continue
        props.update(do_csv.get(asset, {}))
        if args.conferir:
            print(f'  {nome}: {props}')
            continue
        ee.data.setAssetProperties(asset, props)
        print(f'  {nome}: {len(props)} propriedades gravadas')
        gravados += 1

    print(f'Gravadas {gravados}; já tinham rodada {pulados} (--todas regrava)')
    if fora:
        print(f'Fora da convenção <rodada>_p<k>[_preenchido][_<LxCy>] ({len(fora)}), puladas: {", ".join(fora)}')


if __name__ == '__main__':
    main()
