"""
Versão Drive do baixar_tiles.py: exporta os canais de um modelo da CNN para a
área toda por tasks do GEE (Export.image.toDrive) — sem servidor e sem rclone.

Uma task por célula da grid 5×5: GeoTIFF float32 na grade de 30 m do mosaico
(EPSG:4326, mesma grade do computePixels), com BORDA px a mais de cada lado
(padrão 32 = meio patch), sem dado = NODATA. Canais iguais aos do treino
(exportar_patches.imagem_canais + MODELOS):
  --modelo A   16 canais    --modelo B1  15 canais

Saída no Drive da conta autenticada no Earth Engine, pasta --pasta-drive:
  <pasta>/canais_<modelo>_<LxCy>[-<shard>].tif   (~1 GB por célula; o GEE pode
                                                  dividir em mais de um arquivo)

Depois, no Colab: cortar_tiles_drive.py corta cada célula nos 6 × 6 pedaços
com borda do baixar_tiles.py (mesmos nomes de banda e tags), em
TILES_DIR/canais_<modelo>/<LxCy>/ — o que inferir_cnn.py --tiles lê.

Retomável: tasks com a mesma descrição já pendentes, rodando ou concluídas são
puladas.

Uso
  python exportar_tiles_drive.py --modelo A --projeto ee-solkancengine17
  python exportar_tiles_drive.py --modelo B1 --projeto ee-solkancengine17 --celulas L2C3,L3C3
"""

import argparse
import sys
from pathlib import Path

import ee

AQUI = Path(__file__).resolve().parent
sys.path.append(str(AQUI))
import exportar_patches as ep   # noqa: E402  (canais, grade, grid)
import baixar_tiles as bt       # noqa: E402  (limites da célula na grade)

PASTA_DRIVE = 'tiles_cnn_gee'       # pasta na raiz do Meu Drive (o GEE não cria subpastas)
PREFIXO     = 'cnn_tiles'           # descrição das tasks: <PREFIXO>_<modelo>_<LxCy>


def main(lista=None):
    a = argparse.ArgumentParser(description='Canais da CNN por célula direto no Google Drive (tasks do GEE)')
    a.add_argument('--modelo', required=True, choices=list(ep.MODELOS))
    a.add_argument('--borda', type=int, default=ep.PATCH // 2, help='px extras de cada lado')
    a.add_argument('--celulas', default='', help='ex.: L2C3,L3C3 (vazio = todas)')
    a.add_argument('--pasta-drive', default=PASTA_DRIVE, help='pasta na raiz do Meu Drive')
    a.add_argument('--projeto', default=ep.PROJETO_EE)
    args = a.parse_args(lista)

    ee.Initialize(project=args.projeto)
    canais = ep.MODELOS[args.modelo]
    img = ep.imagem_canais().select(canais)
    grid = ep.ler_grid()
    nomes = [c.strip().upper() for c in args.celulas.split(',') if c.strip()] or sorted(grid)
    feitas = {t['description'] for t in ee.data.getTaskList()
              if t.get('state') in ('READY', 'RUNNING', 'COMPLETED')}
    pendentes = [n for n in nomes if f'{PREFIXO}_{args.modelo}_{n}' not in feitas]
    print(f'modelo {args.modelo}: {len(canais)} canais · {len(pendentes)} de {len(nomes)} células a enviar '
          f'→ Drive da conta do EE, pasta "{args.pasta_drive}" (projeto {args.projeto})')
    if not pendentes or input('Enviar as tasks? [s/N]: ').strip().lower() != 's':
        return

    grade = [ep.RES, 0, ep.X0, 0, -ep.RES, ep.Y0]
    eps = ep.RES / 100          # região 1 % de pixel para dentro: só os pixels [c0-b, c1+b), sem um a mais
    for nome in pendentes:
        c0, c1, l0, l1 = bt.limites_celula(grid[nome])
        b = args.borda
        regiao = ee.Geometry.Rectangle([ep.X0 + (c0 - b) * ep.RES + eps, ep.Y0 - (l1 + b) * ep.RES + eps,
                                        ep.X0 + (c1 + b) * ep.RES - eps, ep.Y0 - (l0 - b) * ep.RES - eps],
                                       proj='EPSG:4326', geodesic=False)
        ee.batch.Export.image.toDrive(image=img, description=f'{PREFIXO}_{args.modelo}_{nome}',
                                      folder=args.pasta_drive, fileNamePrefix=f'canais_{args.modelo}_{nome}',
                                      region=regiao, crs='EPSG:4326', crsTransform=grade,
                                      maxPixels=1e10, fileFormat='GeoTIFF').start()
        print(f'  {nome}: {c1 - c0 + 2 * b} × {l1 - l0 + 2 * b} px')
    print('Acompanhe em https://code.earthengine.google.com/tasks; no fim, rode cortar_tiles_drive.py no Colab.')


if __name__ == '__main__':
    main()
