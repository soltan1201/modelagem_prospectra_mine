"""
Amostras de fosfato numa área reduzida (AR) com ASTER + geofísica do acervo CPRM +
geologia 1:250 mil, para a seleção de variáveis (selecao_variaveis_AR.ipynb) e o
mapa RF/GTB (exportar_favorabilidade_AR.py). Desafios e decisões:
docs/desafios_bases_dados.md.

Área: AREAS[AREA_NOME] (hoje só 'irece').

Variáveis (GRUPOS; imagem_variaveis() monta a mesma imagem para amostra e mapa)
  aster      36 de calcular_features do models_ML_fosfato.py sobre MOSAICO (padrão: o
             final sazonal harmonizado com o Landsat; v1 = aster-mosaico-final)
             (14 bandas, 17 índices, phosphorite_rbd, 4 GLCM) — 30 m
  gama       1042_gama: eU_ppm, eTh_ppm, K_pct, CT, U_Th, U_K, Th_K (~493 m) e
             eU_res = eU − média num círculo de RAIO_RES_M (anomalia local de U)
  mag        1042_mag: MAG_nT (~222 m), mag_thg = amplitude do gradiente horizontal
             (nT/km) e mag_res = MAG − média num círculo de RAIO_RES_M
  geologia   lito_<grupo> = 1 dentro das unidades de LITO_GRUPOS (prefixo do
             SIGLA_UNID do lito_jacobina_250k; fora de todas = 0 em todas) e
             dist_<tipo>_km = distância às estruturas de ESTRUTURAS
             (estrut_jacobina_250k), truncada em DIST_MAX_M
  Geofísica a 30 m pelo vizinho mais próximo: cada pixel de 30 m leva o valor do
  pixel de 500 m (gama) ou 222 m (mag) que o contém. As derivadas (eU_res,
  mag_thg, mag_res) são calculadas na grade original e só depois amostradas a
  30 m. Fora da cobertura geofísica: SEM_DADO e tem_geofisica = 0.

Amostras (tipo)
  presenca     TODOS os pixels de 30 m no buffer de BUFFER_M de cada ocorrência de
               fosfato da área (pontos_referencia), com CODIGO_OBJ, papel, fonte
  fundo        N_FUNDO pixels aleatórios da área a mais de DIST_FUNDO_M das
               ocorrências; 'parte' 0–9 só para baixar em pedaços
  geoquimica   amostras de rocha do geoquimica_p2o5 na área (p2o5_pct, p2o5_lt,
               amostra) — só conferência, nunca treino
  Pixel sem ASTER válido não entra (máscara do calcular_features).

Etapas
  python amostras_AR.py             task toAsset → PASTA/amostras_AR_<area>
  python amostras_AR.py --baixar    asset → src/Dados/maxent/amostras_AR_<area>.csv
                                    + amostras_AR_<area>_grupos.json (variáveis por grupo)
  python amostras_AR.py --baixar --refazer   baixa de novo, sem o cache das partes
                                    (CSV antigo sem a coluna id_amostra)

Nomes: ID_AR = <area> + SUFIXO do mosaico ('_shm' no harmonizado, '' no v1) em
todos os arquivos e assets (amostras_AR_<ID_AR>, CSV, _grupos.json; no notebook e no
exportar_favorabilidade_AR.py também), para não sobrescrever os do v1.

Conta: --conta (padrão caatinga02, de models_ML_fosfato.CONTAS), no projeto EE dela;
precisa escrever em projects/mapbiomas-arida/mine. Ao fim volta à conta que estava ativa.
"""

import argparse
import io
import json
import logging
import sys
from pathlib import Path

import ee

RAIZ = Path(__file__).resolve().parents[2]
sys.path.append(str(RAIZ))
sys.path.append(str(RAIZ / 'modeling' / 'regiao_completa'))   # models_ML_fosfato

import models_ML_fosfato as mx   # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
PASTA   = 'projects/mapbiomas-arida/mine'
ACERVO  = f'{PASTA}/acervo_cprm'
AREAS = {
    'irece': [-42.31762892604435, -12.359936623152336, -41.2382221877631, -10.5375816299746],
}
AREA_NOME = 'irece'
MOSAICO   = 'aster-mosaico-final_sazonal_hm_razao'   # v1: 'aster-mosaico-final'
SUFIXOS   = {'aster-mosaico-final': '', 'aster-mosaico-final_sazonal_hm_razao': '_shm'}
ID_AR     = AREA_NOME + SUFIXOS[MOSAICO]   # nome dos arquivos/assets (irece_shm)
ASSET_PONTOS = f'{PASTA}/pontos_referencia'
ALVO      = 'fosfato'

BUFFER_M     = 250      # buffer de cada ocorrência (todos os pixels de 30 m dentro)
DIST_FUNDO_M = 2000     # fundo longe das ocorrências
N_FUNDO      = 20000    # ~64% cai na cobertura geofísica
ESCALA       = 30
SEED         = 42
SEM_DADO     = -9999    # geofísica fora da cobertura
RAIO_RES_M   = 5000     # anomalia local (eU_res, mag_res)
DIST_MAX_M   = 20000    # distância às estruturas truncada
CSV_DIR      = RAIZ / 'Dados' / 'maxent'

GAMA_BANDAS = ['eU_ppm', 'eTh_ppm', 'K_pct', 'CT', 'U_Th', 'U_K', 'Th_K']
# grupo → prefixo do SIGLA_UNID (lito_jacobina_250k); 10 das 13 ocorrências estão em NP2sn*
LITO_GRUPOS = {
    'lito_nova_america':  'NP2sn',
    'lito_jussara':       'NP2sj',
    'lito_gabriel':       'NP2sg',
    'lito_irece':         'NP2sr',
    'lito_bebedouro':     'NP2be',
    'lito_morro_chapeu':  'M3N1mc',
    'lito_caboclo':       'MPc',
    'lito_tombador':      'MPt',
    'lito_cobertura':     'Q',
}
# banda → trechos do TIPO do estrut_jacobina_250k (qualquer um deles)
ESTRUTURAS = {
    'dist_falha_km':       ['Falha'],
    'dist_dobra_km':       ['clinal', 'forme'],          # anti/sinclinais, antiformes, sinformes
    'dist_lineamento_km':  ['Lineamentos'],
}
GRUPOS = {
    'gama':     GAMA_BANDAS + ['eU_res'],
    'mag':      ['MAG_nT', 'mag_thg', 'mag_res'],
    'geologia': list(LITO_GRUPOS) + list(ESTRUTURAS),
}   # + 'aster' (nomes do calcular_features), preenchido em imagem_variaveis()


def area():
    return ee.Geometry.Rectangle(AREAS[AREA_NOME])


def asset_amostras():
    return f'{PASTA}/amostras_AR_{ID_AR}'


# ═══════════════════════════════════════════════════════════════════════════
# VARIÁVEIS
# ═══════════════════════════════════════════════════════════════════════════
def anomalia_local(img, proj):
    """img − média no círculo de RAIO_RES_M, calculada na grade original (proj)."""
    media = img.reduceNeighborhood(ee.Reducer.mean(), ee.Kernel.circle(RAIO_RES_M, 'meters'))
    return img.subtract(media).reproject(proj)


def geofisica():
    """Gama + mag com derivadas na grade original; máscara = cobertura dos dois."""
    gama = ee.Image(f'{ACERVO}/1042_gama').select(GAMA_BANDAS)
    mag = ee.Image(f'{ACERVO}/1042_mag').select('MAG_nT')
    pg, pm = gama.projection(), mag.projection()
    eu_res = anomalia_local(gama.select('eU_ppm'), pg).rename('eU_res')
    # ee.Image.gradient dá a derivada por metro → × 1000 = nT/km (a escala não muda as árvores)
    grad = mag.gradient().reproject(pm)
    thg = grad.select('x').hypot(grad.select('y')).multiply(1000).rename('mag_thg')
    mag_res = anomalia_local(mag, pm).rename('mag_res')
    geo = ee.Image.cat([gama, eu_res, mag, thg, mag_res]).toFloat()
    return geo.updateMask(geo.mask().reduce(ee.Reducer.min()))


def geologia():
    """lito_<grupo> 0/1 e distância às estruturas (km), sem máscara na área."""
    reg = area().buffer(DIST_MAX_M)
    lito = ee.FeatureCollection(f'{ACERVO}/lito_jacobina_250k').filterBounds(reg)
    estr = ee.FeatureCollection(f'{ACERVO}/estrut_jacobina_250k').filterBounds(reg)
    partes = [ee.Image(0).paint(lito.filter(ee.Filter.stringStartsWith('SIGLA_UNID', pre)), 1)
                .toFloat().rename(nome) for nome, pre in LITO_GRUPOS.items()]
    for nome, trechos in ESTRUTURAS.items():
        filtro = ee.Filter.Or(*[ee.Filter.stringContains('TIPO', t) for t in trechos])
        partes.append(estr.filter(filtro).distance(DIST_MAX_M, 30).unmask(DIST_MAX_M).divide(1000)
                          .toFloat().rename(nome))
    return ee.Image.cat(partes)


def imagem_variaveis():
    """(imagem, grupos): ASTER (máscara do pixel) + geofísica com SEM_DADO fora da
    cobertura + tem_geofisica + geologia. grupos = {grupo: [bandas]}."""
    mx.ASSET_MOSAICO = mx.MOSAICOS[MOSAICO]
    aster, nomes, _ = mx.ler_mosaico()
    geo = geofisica()
    tem = geo.mask().reduce(ee.Reducer.min()).unmask(0).rename('tem_geofisica').toByte()
    img = ee.Image.cat([aster, geo.unmask(SEM_DADO), tem, geologia()])
    return img, {'aster': nomes, **GRUPOS}


# ═══════════════════════════════════════════════════════════════════════════
# AMOSTRAS
# ═══════════════════════════════════════════════════════════════════════════
def ocorrencias():
    return (ee.FeatureCollection(ASSET_PONTOS).filter(ee.Filter.eq('alvo', ALVO))
              .filterBounds(area()))


def amostras(img):
    pontos = ocorrencias()

    def presenca(f):
        props = {'tipo': 'presenca', 'classe': 1, 'CODIGO_OBJ': f.get('CODIGO_OBJ'),
                 'papel': f.get('papel'), 'fonte': f.get('fonte'), 'parte': -1}
        return (img.sample(region=f.geometry().buffer(BUFFER_M), scale=ESCALA, geometries=True)
                   .map(lambda s: s.set(props)))

    pres = ee.FeatureCollection(pontos.map(presenca)).flatten()

    longe = area().difference(pontos.map(lambda f: f.buffer(DIST_FUNDO_M)).union(10).geometry(), 10)
    fundo = (img.sample(region=longe, scale=ESCALA, numPixels=N_FUNDO, seed=SEED, geometries=True)
                .randomColumn('_r', SEED)
                .map(lambda s: s.set({'tipo': 'fundo', 'classe': 0, 'CODIGO_OBJ': -1,
                                      'parte': ee.Number(s.get('_r')).multiply(10).floor()})))

    geoq = (img.sampleRegions(
                collection=ee.FeatureCollection(f'{ACERVO}/geoquimica_p2o5').filterBounds(area())
                             .select(['amostra', 'p2o5_pct', 'p2o5_lt']),
                scale=ESCALA, geometries=True)
               .map(lambda s: s.set({'tipo': 'geoquimica', 'classe': -1, 'CODIGO_OBJ': -1, 'parte': -2})))
    return pres.merge(fundo).merge(geoq)


def exportar():
    asset = asset_amostras()
    try:
        ee.data.getAsset(asset)
        raise SystemExit(f'{asset} já existe — apague para refazer, ou rode --baixar.')
    except ee.EEException:
        pass
    img, grupos = imagem_variaveis()
    n = ocorrencias().size().getInfo()
    fc = amostras(img).map(lambda f: f.set({'area_ar': AREA_NOME, 'mosaico': MOSAICO,
                                            'buffer_m': BUFFER_M}))
    nome = asset.rsplit('/', 1)[-1]
    task = ee.batch.Export.table.toAsset(collection=fc, description=nome, assetId=asset)
    task.start()
    log.info(f'{AREA_NOME}: {n} ocorrências de {ALVO}, fundo {N_FUNDO}, '
             f'{sum(map(len, grupos.values()))} variáveis → task {task.id} → {asset}')


# ═══════════════════════════════════════════════════════════════════════════
# BAIXAR
# ═══════════════════════════════════════════════════════════════════════════
def baixar(refazer=False):
    import pandas as pd
    import requests

    asset = asset_amostras()
    _, grupos = imagem_variaveis()
    # id_amostra = system:index do asset: o exportar_favorabilidade_AR.py filtra o fundo por ele
    # quando o notebook escolhe o fundo por cluster (FUNDO_MODO = 'cluster')
    fc = ee.FeatureCollection(asset).map(lambda f: f.set({
        'lon': f.geometry().coordinates().get(0), 'lat': f.geometry().coordinates().get(1),
        'id_amostra': f.get('system:index')}))
    colunas = sorted(set().union(*[ee.Feature(fc.filter(ee.Filter.eq('tipo', t)).first())
                                   .propertyNames().getInfo() for t in ('presenca', 'fundo', 'geoquimica')])
                     - {'system:index', '_r'})
    partes = sorted(fc.aggregate_array('parte').distinct().getInfo())
    cache = CSV_DIR / 'celulas' / f'amostras_AR_{ID_AR}'
    cache.mkdir(parents=True, exist_ok=True)
    dfs = []
    for p in partes:
        arq = cache / f'parte_{int(p)}.csv'
        if refazer or not arq.exists():
            url = fc.filter(ee.Filter.eq('parte', p)).getDownloadURL(filetype='csv', selectors=colunas)
            r = requests.get(url, timeout=600)
            r.raise_for_status()
            arq.write_text(r.text)
        dfs.append(pd.read_csv(arq))
        log.info(f'  parte {int(p)}: {len(dfs[-1])} linhas')
    df = pd.concat(dfs, ignore_index=True)
    saida = CSV_DIR / f'amostras_AR_{ID_AR}.csv'
    df.to_csv(saida, index=False)
    meta = {'asset': asset, 'area': AREAS[AREA_NOME], 'mosaico': MOSAICO, 'buffer_m': BUFFER_M,
            'dist_fundo_m': DIST_FUNDO_M, 'sem_dado': SEM_DADO, 'grupos': grupos}
    (CSV_DIR / f'amostras_AR_{ID_AR}_grupos.json').write_text(json.dumps(meta, indent=1, ensure_ascii=False))
    log.info(f'→ {saida} ({len(df)} linhas; {df.tipo.value_counts().to_dict()})')


def main():
    a = argparse.ArgumentParser(description='Amostras de fosfato na área reduzida (ASTER + geofísica + geologia)')
    a.add_argument('--baixar', action='store_true', help='baixa o asset já exportado para CSV')
    a.add_argument('--refazer', action='store_true', help='com --baixar: ignora as partes já baixadas')
    a.add_argument('--conta', default=mx.CONTAS[0], choices=mx.CONTAS,
                   help='conta GEE (no projeto dela) que dispara a task / baixa')
    args = a.parse_args()
    original = mx.gee_tools.active_user()
    mx.ativar(args.conta)
    log.info(f'Mosaico {MOSAICO} → nomes com {ID_AR}')
    try:
        baixar(args.refazer) if args.baixar else exportar()
    finally:
        if original and original != mx.gee_tools.active_user():
            mx.gee_tools.switch_user(original)
            log.info(f'Conta de volta: {original}')


if __name__ == '__main__':
    main()
