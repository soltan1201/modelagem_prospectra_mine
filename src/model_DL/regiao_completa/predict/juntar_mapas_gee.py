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
  python predict/juntar_mapas_gee.py --listar                 # situação de cada mapa em RUNS_DIR
  python predict/juntar_mapas_gee.py                          # tabela + pergunta quais juntar
  python predict/juntar_mapas_gee.py --filtro _shm --sim      # junta os pendentes _shm sem perguntar
  python predict/juntar_mapas_gee.py --run $RUNS/B1_regnety_64 --passo 3
  python predict/juntar_mapas_gee.py --com-originais          # também as originais quando há _preenchido
"""

import argparse
import csv
import json
import os
import re
from pathlib import Path

COLECAO       = 'projects/mapbiomas-arida/mine/cnn_favorabilidade'
ALVO          = 'fosfato'

# Rodadas do mosaico sazonal harmonizado levam _shm no nome (rodadas_treino.md):
# <modelo>_<backbone curto>_<patch>_shm[_teste][_1fold]. Sem _shm = mosaico v1.
MOSAICO_V1, MOSAICO_SHM = 'aster-mosaico-final', 'aster-mosaico-final_sazonal_hm_razao'


def mosaico_da_rodada(rodada):
    return MOSAICO_SHM if '_shm' in rodada else MOSAICO_V1


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


def folds_usados(pasta):
    """media5 / final / fold<k>[+…] pelo config_inferencia.json do inferir_cnn.py (a versão
    _preenchido lê o da pasta original); 'media5' quando não registrado (padrão antigo)."""
    arq = pasta.with_name(pasta.name.removesuffix('_preenchido')) / 'config_inferencia.json'
    usados = json.loads(arq.read_text()).get('modelos_usados') if arq.exists() else None
    if not usados:
        return 'media5'
    if usados == ['final.pt']:
        return 'final'
    nomes = [u.removesuffix('.pt') for u in usados]
    return 'media5' if len(nomes) == 5 else '+'.join(nomes)


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


N_CELULAS = 25


def situacao(pasta, destino):
    """pendente (nunca juntado) / desatualizado (células mais novas que o .tif) /
    pronto; '(incompleto n/25)' quando faltam células."""
    arquivos = sorted(pasta.glob('p_L*C*.tif'))
    if not destino.exists():
        st = 'pendente'
    elif max(a.stat().st_mtime for a in arquivos) > destino.stat().st_mtime:
        st = 'desatualizado'
    else:
        st = 'pronto'
    return st, arquivos


def escolher(itens, todos):
    """Mostra a tabela e pergunta quais juntar; Enter = pendentes + desatualizados completos."""
    sugeridos = [i for i, it in enumerate(itens)
                 if it['st'] != 'pronto' and len(it['arquivos']) == N_CELULAS]
    print(f'\n  #  {"arquivo":46s} {"mosaico":8s} {"células":>8s}  situação')
    for i, it in enumerate(itens):
        n = len(it['arquivos'])
        inc = '' if n == N_CELULAS else '  (incompleto)'
        print(f'{i:3d}  {it["nome"]:46s} {"shm" if "_shm" in it["run"].name else "v1":8s} '
              f'{n:>4d}/{N_CELULAS}  {it["st"]}{inc}{"  ←" if i in sugeridos else ""}')
    if todos:
        return sugeridos
    r = input(f'\nJuntar quais? [Enter = os {len(sugeridos)} marcados com ←; ex.: 0,3,5; '
              f'"t" = todos; "n" = nenhum]: ').strip().lower()
    if not r:
        return sugeridos
    if r == 'n':
        return []
    if r == 't':
        return list(range(len(itens)))
    return [int(x) for x in r.replace(' ', '').split(',') if x]


def main(lista=None):
    a = argparse.ArgumentParser(description='Junta as células do mapa de p num GeoTIFF por rodada para o GEE')
    a.add_argument('--run', default='', help='pasta da rodada (padrão: todas em --runs)')
    a.add_argument('--runs', default=os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn'))
    a.add_argument('--passo', type=int, default=None, help='padrão: todos os mapas_passo<k>')
    a.add_argument('--saida', default='', help='padrão: <runs>/upload_gee')
    a.add_argument('--com-originais', action='store_true',
                   help='junta também mapas_passo<k> quando existe mapas_passo<k>_preenchido')
    a.add_argument('--filtro', default='',
                   help='só rodadas com este texto no nome (ex.: _shm)')
    a.add_argument('--listar', action='store_true', help='só mostra a situação de cada mapa e sai')
    a.add_argument('--sim', action='store_true',
                   help='sem pergunta: junta os pendentes e desatualizados completos')
    a.add_argument('--sobrescrever', action='store_true', help='refaz também os prontos escolhidos')
    args = a.parse_args(lista)

    raiz = Path(args.runs).expanduser()
    runs = [Path(args.run)] if args.run else sorted(p for p in raiz.iterdir() if p.is_dir())
    runs = [r for r in runs if args.filtro in r.name]
    alvos = listar(runs, args.passo, args.com_originais)
    if not alvos:
        raise SystemExit('Nenhuma pasta mapas_passo<k> com p_L*C*.tif — rode inferir_cnn.py --tiles antes')
    saida = Path(args.saida) if args.saida else raiz / 'upload_gee'
    saida.mkdir(parents=True, exist_ok=True)

    itens = []
    for run, passo, pasta in alvos:
        preenchido = pasta.name.endswith('_preenchido')
        nome = f'{run.name}_p{passo}{"_preenchido" if preenchido else ""}'.replace('.', '_')
        destino = saida / f'{nome}.tif'
        st, arquivos = situacao(pasta, destino)
        itens.append(dict(run=run, passo=passo, pasta=pasta, preenchido=preenchido,
                          nome=nome, destino=destino, st=st, arquivos=arquivos))

    escolhidos = escolher(itens, todos=args.sim or args.listar)
    if args.listar:
        return
    for i in escolhidos:
        it = itens[i]
        if it['st'] == 'pronto' and not args.sobrescrever:
            print(f'  {it["nome"]}: pronto — pulado (--sobrescrever para refazer)')
            continue
        if len(it['arquivos']) < N_CELULAS:
            print(f'  aviso: {it["nome"]} com {len(it["arquivos"])}/{N_CELULAS} células')
        h, w = juntar(it['arquivos'], it['destino'])
        print(f'  {it["nome"]}: {len(it["arquivos"])} células → {w}×{h} px, '
              f'{it["destino"].stat().st_size / 2**20:.1f} MB')

    # CSV com todos os .tif que existem em <saida>, não só os desta vez
    linhas = []
    for it in itens:
        if not it['destino'].exists():
            continue
        run = it['run']
        cfg = json.loads((run / 'config.json').read_text()) if (run / 'config.json').exists() else {}
        linhas.append({'arquivo': it['destino'].name, 'asset': f'{COLECAO}/{it["nome"]}',
                       'folds': folds_usados(it['pasta']),
                       'modelo': f'CNN_{cfg.get("backbone")}', 'alvo': ALVO, 'mosaico': mosaico_da_rodada(run.name),
                       'rodada': run.name, 'passo': it['passo'], 'preenchido': int(it['preenchido']),
                       'canais': cfg.get('modelo'), 'patch': cfg.get('patch'), 'escala_p': 10000,
                       'celulas': len(it['arquivos'])})
    if not linhas:
        print('Nenhum .tif juntado ainda.')
        return
    arq = saida / 'propriedades_upload.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    print(f'{len(linhas)} arquivos em {saida}\nPropriedades → {arq}')


if __name__ == '__main__':
    main()
