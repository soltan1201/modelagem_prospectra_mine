"""
Faixas sem p na beira das células do mapa da CNN: diagnóstico e preenchimento
sem prever de novo.

A janela deslizante não deveria deixar faixa nenhuma: os pedaços de canais têm
32 px de borda vindos do mosaico inteiro (inclusive da célula vizinha), então o
miolo de cada célula é previsto até a beira. Uma faixa sem p tem uma de duas
causas, e o diagnóstico separa as duas:

  entrada   os canais (tiles) também não têm dado ali — buraco do mosaico ASTER;
            não há o que prever
  predição  os canais têm dado e o p não — falha da inferência; o certo é prever
            só esses pedaços de novo (inferir_cnn.py pula os que existem: apague
            os pedaços afetados e rode de novo)

Diagnóstico (sempre), por célula: largura da faixa sem p em cada lado (N, S, O, L,
em px do mapa), quantos px sem p têm canais com dado (--tiles) e quantos a célula
vizinha cobre (as 8 vizinhas: L2C2 ← L1C2, L3C2, L2C1, L2C3 e diagonais).

Preenchimento (--preencher), sem tocar nos originais:
  1. vizinhas   p da célula vizinha no mesmo lugar (vizinho mais próximo na grade
                da célula), onde a célula não tem p
  2. --interpolar N   o que sobrar, até N px de um p válido, por interpolação
                (rasterio.fill.fillnodata) — só estético: p inventado, não previsto
  Saída: <mapas>_preenchido/p_<LxCy>_<rodada>.tif (mesmo formato: UInt16,
  p × 10000, 0 = sem dado) e relatorio_bordas.csv com a origem de cada px.

Uso (servidor; troque RODADA, MOD = A ou B1 conforme a rodada, e PASSO)
  cd /srv/almacen/projects/modelagem_prospectra_mine/src/model_DL
  RUNS=/srv/almacen/db_images/runs_cnn; TILES=/srv/almacen/db_images/tiles_cnn
  RODADA=B1_coatnet0_64; MOD=B1; PASSO=3

  # 1. diagnóstico (não grava mapas)
  python preencher_bordas.py --mapas $RUNS/$RODADA/mapas_passo$PASSO --tiles $TILES/canais_$MOD
  # 2. preencher só com as vizinhas
  python preencher_bordas.py --mapas $RUNS/$RODADA/mapas_passo$PASSO --preencher
  # 3. vizinhas + interpolação até 10 px
  python preencher_bordas.py --mapas $RUNS/$RODADA/mapas_passo$PASSO --preencher --interpolar 10

  # TODAS as rodadas com .pt e mapas (sem --mapas): RUNS_DIR e TILES_DIR, canais pelo config.json
  python preencher_bordas.py                                    (diagnóstico)
  python preencher_bordas.py --preencher --interpolar 10        (preenche todas)
  → resumo em $RUNS/relatorio_bordas_rodadas.csv
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

NODATA_CANAIS = -9999.0     # o mesmo NODATA dos pedaços de canais (exportar_patches.NODATA)


def celula_de(arq):
    m = re.match(r'p_L(\d+)C(\d+)', arq.name)
    return (int(m.group(1)), int(m.group(2))) if m else None


def faixa(valido, lado, minimo=0.05):
    """Linhas/colunas seguidas, a partir de um lado, com menos de `minimo` de px com p."""
    frac = {'N': valido.mean(1), 'S': valido.mean(1)[::-1],
            'O': valido.mean(0), 'L': valido.mean(0)[::-1]}[lado]
    ok = np.flatnonzero(frac >= minimo)
    return int(ok[0]) if len(ok) else len(frac)


def para_grade(fonte, banda, perfil, nearest=True):
    """Lê `banda` de um arquivo aberto reprojetado na grade de `perfil` (0 fora dele)."""
    from rasterio.warp import Resampling, reproject
    destino = np.zeros((perfil['height'], perfil['width']), banda.dtype)
    reproject(banda, destino, src_transform=fonte.transform, src_crs=fonte.crs, src_nodata=0,
              dst_transform=perfil['transform'], dst_crs=perfil['crs'], dst_nodata=0,
              resampling=Resampling.nearest if nearest else Resampling.max)
    return destino


def cobertura_entrada(pasta_celula, perfil):
    """1 onde o primeiro canal dos pedaços da célula tem dado, na grade do mapa."""
    import rasterio
    cob = np.zeros((perfil['height'], perfil['width']), np.uint8)
    for arq in sorted(pasta_celula.glob('*_r*c*.tif')):
        with rasterio.open(arq) as d:
            x = d.read(1)
            ok = (np.isfinite(x) & (x != NODATA_CANAIS)).astype(np.uint8)
            cob = np.maximum(cob, para_grade(d, ok, perfil, nearest=False))
    return cob.astype(bool)


def processar(pasta, tiles, preencher, interpolar):
    """Diagnóstico (e preenchimento) de uma pasta mapas_passo<k>; devolve a tabela por célula."""
    import rasterio
    arquivos = {celula_de(f): f for f in sorted(pasta.glob('p_L*C*.tif')) if celula_de(f)}
    if not arquivos:
        print(f'{pasta}: nenhum p_<LxCy>.tif — pulada')
        return None
    saida = pasta.with_name(pasta.name + '_preenchido')
    if preencher:
        saida.mkdir(exist_ok=True)

    linhas = []
    for (li, co), arq in sorted(arquivos.items()):
        nome = f'L{li}C{co}'
        with rasterio.open(arq) as d:
            p = d.read(1)
            perfil = d.profile.copy()
        valido = p > 0
        sem_p = ~valido
        r = {'celula': nome, 'px': p.size, 'sem_p': int(sem_p.sum()),
             **{f'faixa_{s}': faixa(valido, s) for s in 'NSOL'}}

        if tiles and (Path(tiles) / nome).is_dir():
            cob = cobertura_entrada(Path(tiles) / nome, perfil)
            r['sem_p_com_entrada'] = int((sem_p & cob).sum())      # > 0: falha da predição

        vizinhas = np.zeros_like(p)
        for dl in (-1, 0, 1):
            for dc in (-1, 0, 1):
                v = arquivos.get((li + dl, co + dc))
                if (dl or dc) and v is not None:
                    with rasterio.open(v) as dv:
                        q = para_grade(dv, dv.read(1), perfil)
                    vizinhas = np.where(vizinhas == 0, q, vizinhas)
        usa_viz = sem_p & (vizinhas > 0)
        r['preenchido_vizinhas'] = int(usa_viz.sum())

        if preencher:
            novo = np.where(usa_viz, vizinhas, p)
            r['preenchido_interpolacao'] = 0
            if interpolar:
                from rasterio.fill import fillnodata
                antes = novo == 0
                novo = fillnodata(novo.astype(np.float32), mask=(novo > 0).astype(np.uint8),
                                  max_search_distance=interpolar, smoothing_iterations=0)
                novo = np.clip(np.round(novo), 0, 10000).astype(np.uint16)
                r['preenchido_interpolacao'] = int((antes & (novo > 0)).sum())
            r['sem_p_final'] = int((novo == 0).sum())
            with rasterio.open(saida / arq.name, 'w', **perfil) as d:
                d.write(novo.astype(perfil['dtype']), 1)
                d.set_band_description(1, 'probability')
                d.update_tags(escala='p x 10000', nodata='0', preenchido='vizinhas'
                              + (f'+interpolacao_{interpolar}px' if interpolar else ''))
        linhas.append(r)
        print('  ' + ' · '.join(f'{k} {v}' for k, v in r.items()))

    tab = pd.DataFrame(linhas)
    destino = (saida if preencher else pasta) / 'relatorio_bordas.csv'
    tab.to_csv(destino, index=False)
    print(f'  → {destino}')
    if 'sem_p_com_entrada' in tab and tab.sem_p_com_entrada.sum():
        print('  ATENÇÃO: há px sem p onde os canais têm dado — falha da predição, não do mosaico: '
              'apague os pedaços dessas células e rode inferir_cnn.py de novo (só eles são refeitos).')
    return tab


def main(lista=None):
    import json
    import os
    a = argparse.ArgumentParser(description='Faixas sem p na beira das células: diagnóstico e preenchimento')
    a.add_argument('--mapas', default='', help='uma pasta mapas_passo<k>; vazio = todas as rodadas com .pt em --runs')
    a.add_argument('--runs', default=os.environ.get('RUNS_DIR', '/srv/almacen/db_images/runs_cnn'))
    a.add_argument('--tiles', default='', help='pasta canais_<modelo> (só com --mapas); diz se a faixa é '
                   'buraco da entrada. Em todas as rodadas: TILES_DIR/canais_<modelo do config.json>')
    a.add_argument('--tiles-dir', default=os.environ.get('TILES_DIR', '/srv/almacen/db_images/tiles_cnn'))
    a.add_argument('--preencher', action='store_true', help='grava <mapas>_preenchido/ com as vizinhas')
    a.add_argument('--interpolar', type=int, default=0, help='px máximos de interpolação no que sobrar (0 = não)')
    args = a.parse_args(lista)

    if args.mapas:
        if processar(Path(args.mapas), args.tiles, args.preencher, args.interpolar) is None:
            raise SystemExit(f'Nenhum p_<LxCy>.tif em {args.mapas}')
        return

    # todas as rodadas com modelo salvo (.pt) e alguma pasta mapas_passo<k>
    resumo = []
    for run in sorted(p for p in Path(args.runs).iterdir() if p.is_dir() and any(p.glob('*.pt'))):
        cfg = json.loads((run / 'config.json').read_text()) if (run / 'config.json').exists() else {}
        modelo = cfg.get('modelo') or run.name.split('_')[0]
        tiles = Path(args.tiles_dir) / f'canais_{modelo}'
        for pasta in sorted(run.glob('mapas_passo*')):
            if not pasta.is_dir() or not pasta.name.removeprefix('mapas_passo').isdigit():
                continue                                     # pula as _preenchido
            print(f'\n═══ {run.name} / {pasta.name} (canais_{modelo})')
            tab = processar(pasta, str(tiles) if tiles.is_dir() else '', args.preencher, args.interpolar)
            if tab is not None:
                resumo.append({'rodada': run.name, 'mapas': pasta.name, 'celulas': len(tab),
                               'sem_p': int(tab.sem_p.sum()),
                               'sem_p_com_entrada': int(tab.get('sem_p_com_entrada', pd.Series(dtype=int)).sum()),
                               'preenchido_vizinhas': int(tab.preenchido_vizinhas.sum()),
                               'sem_p_final': int(tab['sem_p_final'].sum()) if 'sem_p_final' in tab else None})
    if not resumo:
        raise SystemExit(f'Nenhuma rodada com .pt e mapas_passo<k> em {args.runs}')
    r = pd.DataFrame(resumo)
    r.to_csv(Path(args.runs) / 'relatorio_bordas_rodadas.csv', index=False)
    print('\n' + r.to_string(index=False))
    print(f'→ {Path(args.runs) / "relatorio_bordas_rodadas.csv"}')


if __name__ == '__main__':
    main()
