"""
Converte a base exportada por exportar_patches_drive.py (TFRecord das tasks do
GEE no Drive) para o formato do treino — o mesmo do exportar_patches.py:

  <destino>/modelo_A/patches_<LxCy>.npz    X (N, 16, 64, 64) float16 + metadados
  <destino>/modelo_B1/patches_<LxCy>.npz   X (N, 15, 64, 64) float16 + metadados
  <destino>/indice_patches.csv, divisao.json

Feito para o Colab (TensorFlow já vem instalado), com o Drive montado.
Retomável: célula com os dois .npz já gravados é pulada. Confere o número de
patches de cada célula contra o indice_patches.csv.

Uso (Colab)
  !python tfrecord_para_npz.py --origem /content/drive/MyDrive/patches_cnn_gee \
                               --destino /content/drive/MyDrive/prospectra/patches_cnn
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
sys.path.append(str(AQUI))
import exportar_patches as ep   # noqa: E402  (CANAIS, MODELOS, gravar_modelos)

TEXTO = ['id', 'origem_amostra', 'conjunto', 'celula', 'papel']


def ler_registros(arquivos):
    """Exemplos TFRecord → (X (N, len(CANAIS), PATCH, PATCH) float16, DataFrame de metadados)."""
    import tensorflow as tf
    comp = 'GZIP' if arquivos[0].name.endswith('.gz') else ''
    Xs, linhas = [], []
    for bruto in tf.data.TFRecordDataset([str(a) for a in arquivos], compression_type=comp):
        f = tf.train.Example.FromString(bruto.numpy()).features.feature
        Xs.append(np.stack([np.asarray(f[c].float_list.value, np.float32).reshape(ep.PATCH, ep.PATCH)
                            for c in ep.CANAIS]))
        linha = {}
        for k, v in f.items():
            if k in ep.CANAIS:
                continue
            tipo = v.WhichOneof('kind')
            valor = getattr(v, tipo).value[0] if getattr(v, tipo).value else None
            linha[k] = valor.decode() if isinstance(valor, bytes) else valor
        linhas.append(linha)
    X = np.stack(Xs)
    X[X == ep.NODATA] = np.nan
    meta = pd.DataFrame(linhas)
    for c in meta.columns:
        if c not in TEXTO:
            try:
                meta[c] = pd.to_numeric(meta[c])
            except (ValueError, TypeError):    # ex.: codigo_obj com texto
                meta[c] = meta[c].astype(str)
    for c in ('fold', 'uso_treino', 'em_magenta', 'jitter'):
        if c in meta:
            meta[c] = meta[c].astype(int)
    ordem = np.argsort(meta['id'].to_numpy())
    return X[ordem].astype(np.float16), meta.iloc[ordem].reset_index(drop=True)


def main(lista=None):
    a = argparse.ArgumentParser(description='TFRecord do GEE (exportar_patches_drive.py) → .npz do treino')
    a.add_argument('--origem', default='/content/drive/MyDrive/patches_cnn_gee')
    a.add_argument('--destino', default=os.environ.get('PATCHES_DIR', '/content/drive/MyDrive/prospectra/patches_cnn'))
    args = a.parse_args(lista)

    origem, destino = Path(args.origem), Path(args.destino)
    ep.PASTA_SAIDA = destino
    for m in ep.MODELOS:
        (destino / f'modelo_{m}').mkdir(parents=True, exist_ok=True)

    indice = None
    if (origem / 'indice_patches.csv').exists():
        indice = pd.read_csv(origem / 'indice_patches.csv')
        indice.drop(columns=['system:index', '.geo'], errors='ignore').to_csv(destino / 'indice_patches.csv',
                                                                            index=False)
    if (origem / 'divisao.csv').exists():
        divisao = json.loads(pd.read_csv(origem / 'divisao.csv')['json'].iloc[0])
        (destino / 'divisao.json').write_text(json.dumps(divisao, indent=2, ensure_ascii=False))

    por_celula = defaultdict(list)
    for arq in sorted(origem.glob('patches_*.tfrecord*')):
        m = re.match(r'patches_(L\d+C\d+)_\d+', arq.name)
        if m:
            por_celula[m.group(1)].append(arq)
    if not por_celula:
        raise SystemExit(f'Nenhum patches_*.tfrecord em {origem} — as tasks já terminaram?')

    for cel, arquivos in sorted(por_celula.items()):
        if all(ep.arquivo_modelo(m, cel).exists() for m in ep.MODELOS):
            print(f'{cel}: já convertido — pulado')
            continue
        X, meta = ler_registros(arquivos)
        esperado = int((indice.celula == cel).sum()) if indice is not None else len(meta)
        aviso = '' if len(meta) == esperado else f'  ATENÇÃO: índice tem {esperado} (falta task?)'
        ep.gravar_modelos(X, meta, cel)
        print(f'{cel}: {len(arquivos)} arquivo(s), {len(meta)} patches → modelo_A/ e modelo_B1/{aviso}')


if __name__ == '__main__':
    main()
