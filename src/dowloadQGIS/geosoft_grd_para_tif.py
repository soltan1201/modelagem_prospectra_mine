"""
Converte os grids Geosoft (.grd binário) do aerolevantamento 1042 — Norte da
Chapada Diamantina — para GeoTIFF float32 e calcula as razões U/Th, U/K e Th/K.

O GDAL não lê o .grd binário da Geosoft. O formato é simples:
  cabeçalho de 512 bytes: ES, SF, NE, NV, KX (int32) e DE, DV, X0, Y0, ROT,
  ZBASE, ZMULT (float64). ES > 1024 = grid comprimido (tamanho do elemento =
  ES - 1024); KX = 1 = vetores são linhas, de baixo (Y0) para cima.
  compressão: int32 magic, tipo, nº de blocos, vetores por bloco; offsets dos
  blocos (int64); tamanhos (int32); cada bloco = 16 bytes de cabeçalho + zlib.
  X0, Y0 = centro da célula do canto inferior esquerdo. Valor nulo = -1e32.

Unidades (Norte_Chapada_Diamantina_1042.doc): CT (contagem total), K %,
eTh ppm, eU ppm, MAG nT. Gama a 0,004° (~500 m); MAG a 0,002° (~250 m).
Coordenadas geográficas WGS 84.

Entrada : /mnt/Almacen/shpIRECE/aerogeofisica/943_1042-XYZ/*.grd
Saída   : /mnt/Almacen/shpIRECE/aerogeofisica/1042_grids/1042_<produto>.tif

Uso: python geosoft_grd_para_tif.py
"""
import struct
import zlib
from pathlib import Path

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()

ORIGEM = Path('/mnt/Almacen/shpIRECE/aerogeofisica/943_1042-XYZ')
SAIDA = Path('/mnt/Almacen/shpIRECE/aerogeofisica/1042_grids')
NULO_GEOSOFT = -1e32
NODATA = -9999.0

GRIDS = {                       # arquivo → (produto, unidade)
    '1042ullwgs84.grd': ('eU_ppm', 'ppm'),
    '1042thllwgs84.grd': ('eTh_ppm', 'ppm'),
    '1042kllwgs84.grd': ('K_pct', '%'),
    '1042tcllwgs84.grd': ('CT', 'contagem total'),
    'nchapadamf.grd': ('MAG_nT', 'nT'),
}
# denominador mínimo para a razão: abaixo disso a razão explode no ruído
MIN_K, MIN_TH = 0.1, 0.5        # % e ppm


def ler_grd(arquivo: Path):
    """(matriz float64 com NaN, geotransform) do grid Geosoft."""
    b = arquivo.read_bytes()
    es, sf, ne, nv, kx = struct.unpack('<5i', b[:20])
    de, dv, x0, y0, rot, zbase, zmult = struct.unpack('<7d', b[20:76])
    if kx != 1 or rot != 0:
        raise ValueError(f'{arquivo.name}: KX={kx}, ROT={rot} não suportados')
    tam = es - 1024 if es > 1024 else es
    tipo = {(4, 2): '<f4', (8, 2): '<f8', (2, 1): '<i2', (4, 1): '<i4', (1, 1): '<i1',
            (2, 0): '<u2', (4, 0): '<u4', (1, 0): '<u1'}[(tam, sf)]

    if es > 1024:
        _, _, nblocos, _ = struct.unpack('<4i', b[512:528])
        offs = struct.unpack(f'<{nblocos}q', b[528:528 + 8 * nblocos])
        p = 528 + 8 * nblocos
        tams = struct.unpack(f'<{nblocos}i', b[p:p + 4 * nblocos])
        bruto = b''.join(zlib.decompress(b[o + 16:o + t]) for o, t in zip(offs, tams))
    else:
        bruto = b[512:512 + ne * nv * tam]

    z = np.frombuffer(bruto, dtype=tipo, count=ne * nv).astype('f8').reshape(nv, ne)
    nulo = z <= NULO_GEOSOFT / 10 if sf == 2 else np.zeros_like(z, bool)
    z = z / zmult + zbase
    z[nulo] = np.nan
    z = z[::-1]                                   # linha 0 = norte, como no GeoTIFF
    gt = (x0 - de / 2, de, 0.0, y0 + (nv - 1) * dv + dv / 2, 0.0, -dv)
    return z, gt


def gravar(z, gt, destino: Path, descricao: str):
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(str(destino), z.shape[1], z.shape[0], 1, gdal.GDT_Float32,
                    options=['COMPRESS=DEFLATE', 'PREDICTOR=3', 'TILED=YES'])
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    banda = ds.GetRasterBand(1)
    banda.WriteArray(np.where(np.isnan(z), NODATA, z).astype('f4'))
    banda.SetNoDataValue(NODATA)
    banda.SetDescription(descricao)
    banda.ComputeStatistics(False)
    ds = None


def main():
    SAIDA.mkdir(parents=True, exist_ok=True)
    grids = {}
    for arq, (produto, unidade) in GRIDS.items():
        z, gt = ler_grd(ORIGEM / arq)
        grids[produto] = (z, gt)
        gravar(z, gt, SAIDA / f'1042_{produto}.tif', f'{produto} ({unidade})')
        v = z[~np.isnan(z)]
        print(f'{produto:8s} {z.shape[1]}x{z.shape[0]}  pixel {gt[1]:.5f}°  '
              f'min {v.min():9.3f}  p50 {np.median(v):9.3f}  max {v.max():9.3f}  ({unidade})')

    # razões: mesma malha nos três canais gama
    (u, gt), (th, gt_th), (k, gt_k) = grids['eU_ppm'], grids['eTh_ppm'], grids['K_pct']
    assert u.shape == th.shape == k.shape and np.allclose(gt, gt_th) and np.allclose(gt, gt_k), \
        'canais gama em malhas diferentes'
    with np.errstate(divide='ignore', invalid='ignore'):
        razoes = {
            'U_Th': np.where(th >= MIN_TH, u / th, np.nan),
            'U_K': np.where(k >= MIN_K, u / k, np.nan),
            'Th_K': np.where(k >= MIN_K, th / k, np.nan),
        }
    for nome, z in razoes.items():
        gravar(z, gt, SAIDA / f'1042_{nome}.tif',
               f'razão {nome.replace("_", "/")} (K >= {MIN_K}%, eTh >= {MIN_TH} ppm)')
        v = z[~np.isnan(z)]
        print(f'{nome:8s} min {v.min():.3f}  p50 {np.median(v):.3f}  p99 {np.percentile(v, 99):.3f}  '
              f'max {v.max():.3f}  ({np.isnan(z).sum() - np.isnan(u).sum()} pixels cortados pelo denominador mínimo)')
    print(f'Saída: {SAIDA}')


if __name__ == '__main__':
    main()
