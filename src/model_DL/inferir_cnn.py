"""
Inferência da CNN de favorabilidade de fosfato sobre o mosaico ASTER v1, célula
a célula da grid 5×5, por janela deslizante: o patch centrado em cada pixel dá o
p daquele pixel (plano em plano_cnn.md).

Entrada
  --run   pasta de uma rodada de treinar_cnn.py. Usa final.pt se existir; senão
          a média dos fold<k>.pt (ensemble dos folds, cada um com a própria
          normalização).
  Canais calculados no GEE com as mesmas funções da exportação
  (exportar_patches.imagem_canais) e baixados por ee.data.computePixels em
  blocos de BLOCO × BLOCO px com borda de PATCH/2 px (sem tasks).

Janela deslizante
  --passo 3   um p a cada 3 pixels (90 m) — 9× mais rápido; use passo ímpar
              para o centro do pixel de saída cair no centro do patch
  --passo 1   um p por pixel (30 m), mapa final
  Patches de uma faixa de linhas são montados com unfold e previstos em lotes na
  GPU. Pixel central sem dado → sem predição (NaN).

Saída (--saida, padrão <run>/mapas_passo<k>/)
  p_<LxCy>.tif   GeoTIFF float32 (EPSG:4326, pixel = passo × 30 m), NaN sem dado;
                 com rasterio ou GDAL. Sem nenhum dos dois: p_<LxCy>.npz com a
                 matriz e a geotransformação.
  Para juntar: gdalbuildvrt mapa.vrt p_*.tif

Modo pedaços (--tiles, recomendado para o mapa inteiro)
  Lê os GeoTIFFs de canais de baixar_tiles.py (sem GEE), prevê pedaço a pedaço
  (--sub divide o pedaço se faltar memória; --lote controla a GPU) e grava
  <saida>/pedacos/<LxCy>/p_<LxCy>_r<i>c<j>.tif e, com os 36 prontos,
  <saida>/p_<LxCy>.tif — UInt16, p × 10000, 0 = sem dado (convenção do
  favorabilidade_ml), prontos para subir_predicao_gee.py.

Uso
  python inferir_cnn.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64 --tiles padrao --passo 1
  python inferir_cnn.py --run <pasta> --tiles padrao --passo 3 --celulas L2C3,L3C3
  python inferir_cnn.py --run <pasta> --passo 3 --celulas L2C3      # sem pedaços: baixa do GEE
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

try:
    import ee
except ImportError:   # --tiles não usa GEE: roda num ambiente só com torch/rasterio
    ee = None

AQUI = Path(__file__).resolve().parent
sys.path.append(str(AQUI))
import exportar_patches as ep   # noqa: E402  (canais, grade e computePixels)
import treinar_cnn as tc        # noqa: E402  (modelos e canal binário)

BLOCO = 512          # px de saída por bloco (sem a borda)


# ═══════════════════════════════════════════════════════════════════════════
# MODELOS
# ═══════════════════════════════════════════════════════════════════════════
def carregar_modelos(run, dispositivo):
    """Lista de (modelo, normalização) — final.pt ou todos os fold<k>.pt."""
    run = Path(run)
    arquivos = [run / 'final.pt'] if (run / 'final.pt').exists() else sorted(run.glob('fold*.pt'))
    if not arquivos:
        raise SystemExit(f'Nenhum final.pt ou fold*.pt em {run}')
    saida = []
    for a in arquivos:
        ck = torch.load(a, map_location='cpu', weights_only=False)
        args = ck['args']
        modelo = tc.criar_modelo(args['backbone'], len(ck['canais']), args['patch'], pretrained=False)
        modelo.load_state_dict(ck['estado'])
        modelo.to(dispositivo).eval()
        est = np.array([ck['est'][c] for c in ck['canais']], np.float32)       # (C, 4): lo, hi, mu, sd
        saida.append({'modelo': modelo, 'est': torch.from_numpy(est).to(dispositivo),
                      'bin': torch.tensor([c == tc.CANAL_BIN for c in ck['canais']], device=dispositivo)})
        print(f'  {a.name}: {args["backbone"]}, {len(ck["canais"])} canais, patch {args["patch"]} px')
    return saida, ck['canais'], args['patch'], args['modelo']


def normalizar(x, m):
    """x (C, H, W) bruto → padronizado como no treino (corte, média, desvio; binário cru; NaN → 0)."""
    lo, hi, mu, sd = (m['est'][:, k][:, None, None] for k in range(4))
    z = (torch.minimum(torch.maximum(x, lo), hi) - mu) / sd
    z = torch.where(m['bin'][:, None, None], x, z)
    return torch.nan_to_num(z, nan=0.0)


# ═══════════════════════════════════════════════════════════════════════════
# DADOS
# ═══════════════════════════════════════════════════════════════════════════
def limites_celula(c):
    """Colunas/linhas (grade de 30 m do mosaico) que cobrem a célula."""
    c0 = math.floor((c['oeste'] - ep.X0) / ep.RES)
    c1 = math.ceil((c['leste'] - ep.X0) / ep.RES)
    l0 = math.floor((ep.Y0 - c['norte']) / ep.RES)
    l1 = math.ceil((ep.Y0 - c['sul']) / ep.RES)
    return c0, c1, l0, l1


def baixar_bloco(img, canais, col, lin, largura, altura):
    """Canais (C, altura, largura) a partir da coluna/linha col/lin da grade de 30 m."""
    grade = {'dimensions': {'width': largura, 'height': altura},
             'affineTransform': {'scaleX': ep.RES, 'shearX': 0, 'translateX': ep.X0 + col * ep.RES,
                                 'shearY': 0, 'scaleY': -ep.RES, 'translateY': ep.Y0 - lin * ep.RES},
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


# ═══════════════════════════════════════════════════════════════════════════
# JANELA DESLIZANTE
# ═══════════════════════════════════════════════════════════════════════════
def prever_bloco(x, modelos, patch, passo, lote, dispositivo):
    """x (C, H + patch, W + patch) com borda → p (ceil(H/passo), ceil(W/passo))."""
    meio = patch // 2
    C, Ht, Wt = x.shape
    H, W = Ht - patch, Wt - patch
    linhas = list(range(0, H, passo))
    colunas = list(range(0, W, passo))
    xt = torch.from_numpy(x).to(dispositivo)
    centro_ok = torch.isfinite(xt[:, meio:meio + H, meio:meio + W]).any(0)       # pixel central com dado
    p = torch.zeros((len(linhas), len(colunas)), device=dispositivo)     # soma dos modelos
    cols_t = torch.tensor(colunas, device=dispositivo)
    for m in modelos:
        xn = normalizar(xt, m)
        for i, l in enumerate(linhas):
            faixa = xn[:, l:l + patch, :]                                     # (C, patch, Wt)
            jan = faixa.unfold(2, patch, 1)[:, :, cols_t].permute(2, 0, 1, 3)  # (n, C, patch, patch)
            saida = []
            for k in range(0, len(colunas), lote):
                with torch.no_grad(), tc.autocast(dispositivo, tc.tipo_amp(dispositivo)):   # bf16 se houver
                    saida.append(torch.sigmoid(m['modelo'](jan[k:k + lote].contiguous())).float())
            p[i] += torch.cat(saida)
    p = p / len(modelos)                                                  # média do ensemble
    ok = centro_ok[::passo, ::passo][:p.shape[0], :p.shape[1]]
    p[~ok] = float('nan')
    return p.cpu().numpy()


def inferir_celula(img, canais, modelos, patch, cel, passo, lote, dispositivo):
    c0, c1, l0, l1 = limites_celula(cel)
    W, H = c1 - c0, l1 - l0
    meio = patch // 2
    saida = np.full((math.ceil(H / passo), math.ceil(W / passo)), np.nan, np.float32)
    passo_bloco = BLOCO - BLOCO % passo                 # blocos múltiplos do passo
    for bl in range(0, H, passo_bloco):
        for bc in range(0, W, passo_bloco):
            h, w = min(passo_bloco, H - bl), min(passo_bloco, W - bc)
            x = baixar_bloco(img, canais, c0 + bc - meio, l0 + bl - meio, w + patch, h + patch)
            p = prever_bloco(x, modelos, patch, passo, lote, dispositivo)
            saida[bl // passo: bl // passo + p.shape[0], bc // passo: bc // passo + p.shape[1]] = p
    # geotransformação: pixel de saída = passo × 30 m, centrado no pixel do patch
    deslocamento = (passo - 1) / 2
    gt = (ep.X0 + (c0 - deslocamento) * ep.RES, passo * ep.RES, 0.0,
          ep.Y0 - (l0 - deslocamento) * ep.RES, 0.0, -passo * ep.RES)
    return saida, gt


def gravar(arquivo, p, gt):
    try:
        import rasterio
        from rasterio.transform import Affine
        with rasterio.open(arquivo, 'w', driver='GTiff', width=p.shape[1], height=p.shape[0], count=1,
                           dtype='float32', crs='EPSG:4326', transform=Affine.from_gdal(*gt),
                           nodata=np.nan, compress='deflate') as d:
            d.write(p, 1)
        return arquivo
    except ImportError:
        pass
    try:
        from osgeo import gdal, osr
        ds = gdal.GetDriverByName('GTiff').Create(str(arquivo), p.shape[1], p.shape[0], 1, gdal.GDT_Float32,
                                                  ['COMPRESS=DEFLATE'])
        ds.SetGeoTransform(gt)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        ds.SetProjection(srs.ExportToWkt())
        ds.GetRasterBand(1).WriteArray(p)
        ds.GetRasterBand(1).SetNoDataValue(float('nan'))
        ds = None
        return arquivo
    except ImportError:
        npz = Path(arquivo).with_suffix('.npz')
        np.savez_compressed(npz, p=p, geotransform=np.array(gt), crs='EPSG:4326')
        return npz


# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════
# MODO PEDAÇOS (--tiles): canais já baixados por baixar_tiles.py, sem GEE
# ═══════════════════════════════════════════════════════════════════════════
def inferir_tile(arquivo, canais, modelos, patch, passo, lote, dispositivo, sub=1):
    """Um GeoTIFF de canais (com borda) → p do miolo e a geotransformação da saída."""
    import rasterio
    with rasterio.open(arquivo) as d:
        nomes = list(d.descriptions)
        if nomes != list(canais):
            raise SystemExit(f'{arquivo.name}: canais {nomes} ≠ canais do modelo {canais} '
                             '(baixe os pedaços do mesmo modelo: baixar_tiles.py --modelo …)')
        x = d.read().astype(np.float32)
        borda = int(d.tags().get('borda', patch // 2))
        tr = d.transform
    corte = borda - patch // 2                 # borda maior que o meio patch (ex.: modelo de 32 px)
    if corte < 0:
        raise SystemExit(f'Borda dos pedaços ({borda} px) menor que meio patch ({patch // 2} px).')
    if corte:
        x = x[:, corte:x.shape[1] - corte, corte:x.shape[2] - corte]
    H, W = x.shape[1] - patch, x.shape[2] - patch
    if sub <= 1:
        p = prever_bloco(x, modelos, patch, passo, lote, dispositivo)
    else:                                      # subdivide o miolo se faltar memória
        p = np.full((math.ceil(H / passo), math.ceil(W / passo)), np.nan, np.float32)
        tam = math.ceil(math.ceil(H / sub) / passo) * passo
        tam_c = math.ceil(math.ceil(W / sub) / passo) * passo
        for r in range(0, H, tam):
            for c in range(0, W, tam_c):
                h, w = min(tam, H - r), min(tam_c, W - c)
                pb = prever_bloco(x[:, r:r + h + patch, c:c + w + patch], modelos, patch, passo, lote, dispositivo)
                p[r // passo: r // passo + pb.shape[0], c // passo: c // passo + pb.shape[1]] = pb
    res = tr.a
    x_miolo, y_miolo = tr.c + borda * res, tr.f - borda * res
    deslocamento = (passo - 1) / 2 * res
    gt = (x_miolo - deslocamento, passo * res, 0.0, y_miolo + deslocamento, 0.0, -passo * res)
    return p, gt


def gravar_uint16(arquivo, p, gt):
    """p × 10000 em UInt16 (0 = sem dado), mesma convenção do favorabilidade_ml no GEE."""
    import rasterio
    from rasterio.transform import Affine
    v = np.where(np.isfinite(p), np.clip(np.round(p * 10000), 1, 10000), 0).astype(np.uint16)
    with rasterio.open(arquivo, 'w', driver='GTiff', width=v.shape[1], height=v.shape[0], count=1,
                       dtype='uint16', crs='EPSG:4326', transform=Affine.from_gdal(*gt), nodata=0,
                       compress='deflate', tiled=True, blockxsize=256, blockysize=256) as d:
        d.write(v, 1)
        d.set_band_description(1, 'probability')
        d.update_tags(escala='p x 10000', nodata='0')


def juntar(arquivos, saida):
    """Junta os pedaços de uma célula num GeoTIFF só (UInt16, 0 = sem dado)."""
    import rasterio
    from rasterio.merge import merge
    fontes = [rasterio.open(a) for a in arquivos]
    try:
        mos, tr = merge(fontes, nodata=0)
        perfil = fontes[0].profile.copy()
    finally:
        for f in fontes:
            f.close()
    perfil.update(width=mos.shape[2], height=mos.shape[1], transform=tr)
    with rasterio.open(saida, 'w', **perfil) as d:
        d.write(mos)
        d.set_band_description(1, 'probability')
        d.update_tags(escala='p x 10000', nodata='0')


def rodar_tiles(args, modelos, canais, patch, modelo_dados, pasta, dispositivo):
    padrao = Path(os.environ.get('TILES_DIR', '/srv/almacen/db_images/tiles_cnn')) / f'canais_{modelo_dados}'
    raiz = padrao if args.tiles == 'padrao' else Path(args.tiles)
    celulas = sorted(p.name for p in raiz.iterdir() if p.is_dir())
    pedidas = [c.strip().upper() for c in args.celulas.split(',') if c.strip()]
    celulas = [c for c in celulas if not pedidas or c in pedidas]
    if not celulas:
        raise SystemExit(f'Nenhuma pasta de célula em {raiz} — rode baixar_tiles.py --modelo {modelo_dados}')
    print(f'pedaços: {raiz} · {len(celulas)} células')
    for cel in celulas:
        entradas = sorted((raiz / cel).glob(f'{cel}_r*c*.tif'))
        dest = pasta / 'pedacos' / cel
        dest.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        for arq in entradas:
            saida = dest / f'p_{arq.name}'
            if saida.exists():
                continue
            p, gt = inferir_tile(arq, canais, modelos, patch, args.passo, args.lote, dispositivo, args.sub)
            gravar_uint16(saida, p, gt)
        feitos = sorted(dest.glob('p_*.tif'))
        final = pasta / f'p_{cel}.tif'
        if len(feitos) == len(entradas) and entradas:
            juntar(feitos, final)
            print(f'{cel}: {len(feitos)} pedaços → {final.name} ({time.time() - t0:.0f} s)')
        else:
            print(f'{cel}: {len(feitos)}/{len(entradas)} pedaços (junção quando todos existirem)')


def main(lista=None):
    a = argparse.ArgumentParser(description='Inferência da CNN por janela deslizante')
    a.add_argument('--run', required=True)
    a.add_argument('--passo', type=int, default=3)
    a.add_argument('--celulas', default='', help='ex.: L2C3,L3C3 (vazio = todas)')
    a.add_argument('--lote', type=int, default=256)
    a.add_argument('--saida', default='')
    a.add_argument('--projeto', default=ep.PROJETO_EE)
    a.add_argument('--tiles', default='', help="pasta canais_<modelo> de baixar_tiles.py "
                   "('padrao' = $TILES_DIR/canais_<modelo>); vazio = baixa do GEE na hora")
    a.add_argument('--sub', type=int, default=1, help='subdivide cada pedaço em sub × sub (memória)')
    args = a.parse_args(lista)
    if args.passo % 2 == 0:
        print('Aviso: passo par desloca meio pixel o centro do pixel de saída; prefira 1, 3 ou 5.')

    dispositivo = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    modelos, canais, patch, modelo_dados = carregar_modelos(args.run, dispositivo)
    pasta = Path(args.saida) if args.saida else Path(args.run) / f'mapas_passo{args.passo}'
    pasta.mkdir(parents=True, exist_ok=True)
    if args.tiles:                              # pedaços locais: sem GEE
        (pasta / 'config_inferencia.json').write_text(json.dumps(
            {**vars(args), 'modelo': modelo_dados, 'canais': canais, 'patch': patch,
             'n_modelos': len(modelos)}, indent=2))
        return rodar_tiles(args, modelos, canais, patch, modelo_dados, pasta, dispositivo)
    if ee is None:
        raise SystemExit('Sem earthengine-api neste ambiente: use --tiles padrao (pedaços do baixar_tiles.py) '
                         'ou rode num ambiente com o pacote ee.')
    ee.Initialize(project=args.projeto)

    todos = ep.CANAIS
    faltam = [c for c in canais if c not in todos]
    if faltam:
        raise SystemExit(f'Canais do modelo fora da exportação: {faltam}')
    img = ep.imagem_canais().select(canais)
    grid = ep.ler_grid()
    nomes = [c.strip().upper() for c in args.celulas.split(',') if c.strip()] or sorted(grid)
    (pasta / 'config_inferencia.json').write_text(json.dumps(
        {**vars(args), 'modelo': modelo_dados, 'canais': canais, 'patch': patch,
         'n_modelos': len(modelos), 'celulas': nomes}, indent=2))
    print(f'{len(modelos)} modelo(s) · patch {patch} px · passo {args.passo} ({args.passo * 30} m) · '
          f'{len(nomes)} células · {dispositivo} → {pasta}')

    for nome in nomes:
        arq = pasta / f'p_{nome}.tif'
        if arq.exists() or arq.with_suffix('.npz').exists():
            print(f'{nome}: já existe — pulado')
            continue
        t0 = time.time()
        p, gt = inferir_celula(img, canais, modelos, patch, grid[nome], args.passo, args.lote, dispositivo)
        saida = gravar(arq, p, gt)
        print(f'{nome}: {p.shape[1]}×{p.shape[0]} px, p mediano {np.nanmedian(p):.3f}, '
              f'p ≥ 0,5 em {100 * np.nanmean(p >= 0.5):.1f}% → {Path(saida).name} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
