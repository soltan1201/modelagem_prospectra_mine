"""
Mapas de probabilidade de fosfato com Random Forest e Gradient Tree Boosting no GEE
numa área reduzida (AR), com as variáveis escolhidas no selecao_variaveis_AR.ipynb.

Entrada
  Amostras  PASTA/amostras_AR_<area> (amostras_AR.py): tipo presenca / fundo /
            geoquimica. A geoquímica nunca entra no treino.
  Variáveis src/modeling/models/variaveis_AR_<area>.json (notebook, célula [19]):
            {'conjuntos': {nome: [variáveis]}, 'usa_geofisica': {nome: bool},
             'recomendado': nome, 'n_por_ponto': N, 'fundo_modo': 'aleatorio' | 'cluster'}
  Cluster   com fundo_modo = 'cluster': models/cluster_fundo_AR_<area>.json (célula
            [10e]) = k-means do fundo do notebook: 'variaveis', 'media', 'desvio'
            (padronização), 'centros' (K × variáveis, padronizados) e 'clusters_treino'.

Fundo de treino por cluster (o mesmo do notebook, [10a]–[10e])
  Cada pixel de fundo é padronizado, z = (x − média) / desvio, e vai para o centro
  mais próximo (distância euclidiana). Só os pixels dos clusters_treino (sem
  ocorrência e além do raio do fosfato) entram no treino; o nome do asset ganha
  '_fcl'. Só nos conjuntos com geofísica: os da área toda usam o fundo aleatório.

Treino (o mesmo desenho do notebook)
  presença = N_POR_PONTO pixels sorteados no buffer de cada ocorrência — cada
             ocorrência pesa igual (o GEE não tem peso por amostra)
  fundo    = FATOR_FUNDO × presença, sorteado entre os elegíveis (semente SEED + k)
  Conjunto com geofísica: presença e fundo só com tem_geofisica = 1, e o mapa sai
  só na cobertura geofísica. Sem geofísica: a área toda.
  N_SORTEIOS > 1: um classificador por sorteio de fundo e o mapa é a média.

Saída
  ImageCollection ASSET_SAIDA, uma imagem por <modelo>_fosfato_AR_<ID_AR>_<conjunto>[_fcl]
    (ID_AR = área + '_shm' no mosaico harmonizado, amostras_AR.ID_AR)
    banda 'probability' = p × 10000, Int16 em [0, 10000]
    propriedades modelo, alvo, area_ar, conjunto, features, n_presenca, n_fundo,
                 n_sorteios, parametros, usa_geofisica, fundo_modo, clusters_treino
  p ordena a área; com fundo 3 × presença, p ≈ 0,25 é o neutro.
  Opcional: ASSET_CLUSTER (imagem) com 'cluster' (0..K−1) e 'treino' (1 = cluster do
  fundo de treino), para ver no Code Editor onde o fundo foi excluído.

Contas: as tasks são repartidas entre CONTAS (caatinga02 e caatinga05, de
models_ML_fosfato.py), cada uma no projeto EE dela: um modelo por conta (RF na 1ª, GTB
na 2ª); a imagem de clusters vai com a 1ª. As duas precisam ler e escrever em
projects/mapbiomas-arida/mine. Ao fim volta à conta que estava ativa.

Uso: python exportar_favorabilidade_AR.py   (pergunta conjunto e modelos; mostra o
plano e pede confirmação). Pula assets que já existem e tasks na fila.
"""

import json
import logging
import sys
from pathlib import Path

import ee

sys.path.append(str(Path(__file__).resolve().parents[1] / 'regiao_completa'))
import amostras_AR as ar   # noqa: E402
import models_ML_fosfato as mx   # noqa: E402
from exportar_favorabilidade_ml import MODELOS, SEED, escolher, rotular_fundo, tasks_ativas   # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
PASTA         = ar.PASTA
ASSET_SAIDA   = f'{PASTA}/favorabilidade_AR'
ASSET_CLUSTER = f'{PASTA}/cluster_fundo_AR_{ar.ID_AR}'
ESCALA        = ar.ESCALA
FATOR_FUNDO   = 3
N_SORTEIOS    = 1          # > 1 = média de vários sorteios de fundo (custo × N)
ESCALA_P      = 10000
MODELS_DIR    = Path(__file__).resolve().parents[1] / 'models'   # src/modeling/models (compartilhada)
ARQ_VARIAVEIS = MODELS_DIR / f'variaveis_AR_{ar.ID_AR}.json'   # ID_AR = área + sufixo do mosaico
ARQ_CLUSTER   = MODELS_DIR / f'cluster_fundo_AR_{ar.ID_AR}.json'
CONTAS        = mx.CONTAS   # caatinga02 e caatinga05, cada uma no projeto dela


# ═══════════════════════════════════════════════════════════════════════════
def ler_variaveis():
    if not ARQ_VARIAVEIS.exists():
        raise SystemExit(f'{ARQ_VARIAVEIS} não existe — rode a célula [19] do selecao_variaveis_AR.ipynb.')
    return json.loads(ARQ_VARIAVEIS.read_text())


def ler_cluster():
    if not ARQ_CLUSTER.exists():
        raise SystemExit(f'{ARQ_CLUSTER.name} não existe — rode as células [2] a [10e] do '
                         'selecao_variaveis_AR.ipynb.')
    return json.loads(ARQ_CLUSTER.read_text())


def imagem_cluster(img, cl):
    """Banda 'cluster' (centro mais próximo) e 'treino' (1 = cluster do fundo de treino)."""
    z = (img.select(cl['variaveis']).subtract(ee.Image.constant(cl['media']))
            .divide(ee.Image.constant(cl['desvio'])))
    dist = [z.subtract(ee.Image.constant(c)).pow(2).reduce(ee.Reducer.sum()) for c in cl['centros']]
    menor, lab = dist[0], ee.Image(0)
    for k, d in enumerate(dist[1:], start=1):
        lab = lab.where(d.lt(menor), k)
        menor = menor.min(d)
    treino = lab.remap(cl['clusters_treino'], [1] * len(cl['clusters_treino']), 0)
    return (ee.Image.cat([lab.rename('cluster'), treino.rename('treino')]).toByte()
              .updateMask(img.select('tem_geofisica').eq(1)))


def presenca(fc, n_por_ponto):
    """N_POR_PONTO pixels sorteados por ocorrência (todas pesam igual)."""
    pres = fc.filter(ee.Filter.eq('tipo', 'presenca'))
    codigos = pres.aggregate_array('CODIGO_OBJ').distinct().getInfo()
    partes = [pres.filter(ee.Filter.eq('CODIGO_OBJ', c)).randomColumn('_s', SEED).sort('_s').limit(n_por_ponto)
              for c in codigos]
    return ee.FeatureCollection(partes).flatten(), len(codigos)


def amostras_treino(usa_geo, n_por_ponto, cl=None):
    fc = ee.FeatureCollection(ar.asset_amostras())
    if usa_geo:
        fc = fc.filter(ee.Filter.eq('tem_geofisica', 1))
    pres, n_pts = presenca(fc, n_por_ponto)
    n_pres = pres.size().getInfo()
    if n_pres == 0:
        raise SystemExit(f'Sem presença em {ar.asset_amostras()} (geofísica: {usa_geo}).')
    fundo = fc.filter(ee.Filter.eq('tipo', 'fundo'))
    if cl is not None:
        fundo = rotular_fundo(fundo, cl).filter(ee.Filter.inList('cluster', cl['clusters_treino']))
    n_eleg = fundo.size().getInfo()
    if cl is not None and n_eleg != cl['pixels_treino']:
        log.warning(f'fundo elegível no GEE ({n_eleg}) ≠ notebook ({cl["pixels_treino"]}): '
                    'o asset de amostras mudou depois do notebook?')
    n_fundo = min(FATOR_FUNDO * n_pres, n_eleg)
    sorteios = [fundo.randomColumn('_s', SEED + k).sort('_s').limit(n_fundo) for k in range(N_SORTEIOS)]
    return pres, sorteios, n_pts, n_pres, n_fundo, n_eleg


def mapa(nome_modelo, feats_img, features, pres, sorteios):
    construtor, params = MODELOS[nome_modelo]
    ps = []
    for k, fundo in enumerate(sorteios):
        treino = pres.merge(fundo).select(features + ['classe'])
        clf = (construtor(**{**params, 'seed': SEED + k}).setOutputMode('PROBABILITY')
                 .train(features=treino, classProperty='classe', inputProperties=features))
        ps.append(feats_img.classify(clf).select(0))
    p = ee.ImageCollection(ps).mean() if len(ps) > 1 else ps[0]
    return p.multiply(ESCALA_P).round().clamp(0, ESCALA_P).toInt16().rename('probability')


def exportar_cluster(img, cl, ja_ativas):
    nome = ASSET_CLUSTER.rsplit('/', 1)[-1]
    if nome in ja_ativas:
        log.info(f'  {nome}: na fila — pulado')
        return
    try:
        ee.data.getAsset(ASSET_CLUSTER)
        log.info(f'  {nome}: já existe — pulado')
        return
    except ee.EEException:
        pass
    imc = imagem_cluster(img, cl).set({
        'k': cl['k'], 'variaveis': ','.join(cl['variaveis']), 'raio': cl['raio'],
        'clusters_treino': ','.join(map(str, cl['clusters_treino'])),
        'dist_fosfato': ','.join(map(str, cl['dist_fosfato'])), 'area_ar': ar.AREA_NOME})
    task = ee.batch.Export.image.toAsset(
        image=imc, description=nome, assetId=ASSET_CLUSTER, region=ar.area(), scale=ESCALA,
        crs='EPSG:4326', maxPixels=1e13, pyramidingPolicy={'cluster': 'mode', 'treino': 'mode'})
    task.start()
    log.info(f'  {nome}: task {task.id}')


def exportar(conjunto, features, usa_geo, modelos, n_por_ponto, cl=None, com_cluster=False):
    """Modelos repartidos entre CONTAS; cada conta monta as amostras e dispara os seus mapas."""
    for k, (conta, mods) in enumerate(mx.repartir(modelos, CONTAS).items()):
        mx.ativar(conta)
        log.info(f'══ {conta}: {", ".join(mods)}' + (' + clusters' if com_cluster and k == 0 else ''))
        exportar_conta(conjunto, features, usa_geo, mods, n_por_ponto, cl, com_cluster and k == 0)


def exportar_conta(conjunto, features, usa_geo, modelos, n_por_ponto, cl=None, com_cluster=False):
    img, _ = ar.imagem_variaveis()
    faltam = [f for f in features if f not in img.bandNames().getInfo()]
    if faltam:
        raise SystemExit(f'Variáveis fora da imagem: {faltam}')
    feats_img = img.select(features)
    if usa_geo:
        feats_img = feats_img.updateMask(img.select('tem_geofisica').eq(1))

    try:
        ee.data.getAsset(ASSET_SAIDA)
    except ee.EEException:
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, ASSET_SAIDA)
        log.info(f'Coleção criada: {ASSET_SAIDA}')
    ja = {a['name'].rsplit('/', 1)[-1] for a in ee.data.listAssets({'parent': ASSET_SAIDA}).get('assets', [])}
    ativas = tasks_ativas()   # lista da conta/projeto ativo

    pres, sorteios, n_pts, n_pres, n_fundo, n_eleg = amostras_treino(usa_geo, n_por_ponto, cl)
    fundo_modo = 'cluster' if cl is not None else 'aleatorio'
    log.info(f'{conjunto}: {n_pts} ocorrências, presença {n_pres}, fundo {n_fundo} × {N_SORTEIOS} '
             f'(de {n_eleg} elegíveis, {fundo_modo})')
    for nome_modelo in modelos:
        nome = f'{nome_modelo}_fosfato_AR_{ar.ID_AR}_{conjunto}' + ('_fcl' if cl is not None else '')
        if nome in ja or nome in ativas:
            log.info(f'  {nome}: já existe/na fila — pulado')
            continue
        p = mapa(nome_modelo, feats_img, features, pres, sorteios).set({
            'modelo': nome_modelo, 'alvo': ar.ALVO, 'area_ar': ar.AREA_NOME, 'conjunto': conjunto,
            'features': ','.join(features), 'usa_geofisica': int(usa_geo), 'n_ocorrencias': n_pts,
            'n_presenca': n_pres, 'n_fundo': n_fundo, 'n_fundo_elegivel': n_eleg, 'n_sorteios': N_SORTEIOS,
            'n_por_ponto': n_por_ponto, 'fundo_modo': fundo_modo, 'buffer_m': ar.BUFFER_M, 'mosaico': ar.MOSAICO,
            'clusters_treino': ','.join(map(str, cl['clusters_treino'])) if cl is not None else '',
            'escala_p': ESCALA_P, 'amostras': ar.asset_amostras(),
            'parametros': str(MODELOS[nome_modelo][1])})
        task = ee.batch.Export.image.toAsset(
            image=p, description=nome, assetId=f'{ASSET_SAIDA}/{nome}', region=ar.area(),
            scale=ESCALA, crs='EPSG:4326', maxPixels=1e13, pyramidingPolicy={'probability': 'mean'})
        task.start()
        log.info(f'  {nome}: task {task.id}')
    if com_cluster:
        exportar_cluster(img, cl, ativas)


def main():
    original = mx.gee_tools.active_user()
    mx.ativar(CONTAS[0])
    v = ler_variaveis()
    nomes = list(v['conjuntos'])
    padrao = nomes.index(v['recomendado']) + 1 if v.get('recomendado') in nomes else 1
    conjunto = escolher('Conjunto de variáveis (variaveis_AR):', nomes, padrao)
    features, usa_geo = v['conjuntos'][conjunto], bool(v['usa_geofisica'][conjunto])
    cl = None
    if v.get('fundo_modo') == 'cluster' and usa_geo:
        cl = ler_cluster()
        log.info(f'fundo por cluster: K = {cl["k"]}, clusters de treino {cl["clusters_treino"]} '
                 f'({cl["pixels_treino"]} de {cl["pixels_fundo"]} px no notebook)')
    elif v.get('fundo_modo') == 'cluster':
        log.info(f'{conjunto}: área toda, sem geofísica para agrupar — fundo aleatório')
    nomes_mod = list(MODELOS)
    m = escolher('Modelos:', nomes_mod + ['os dois'], len(nomes_mod) + 1)
    modelos = nomes_mod if m == 'os dois' else [m]
    com_cluster = cl is not None and input(
        f'Exportar também a imagem de clusters do fundo ({ASSET_CLUSTER.rsplit("/", 1)[-1]})? [s/N]: '
    ).strip().lower() == 's'
    print(f'\n{conjunto}: {len(features)} variáveis ({", ".join(features)})'
          f'\nárea {ar.AREA_NOME} · {"só a cobertura geofísica" if usa_geo else "área toda"} · '
          f'{v["n_por_ponto"]} px por ocorrência · fundo {"cluster" if cl else "aleatorio"}, '
          f'{FATOR_FUNDO} × presença × {N_SORTEIOS} sorteio(s)'
          f'\n→ {len(modelos)} mapa(s) em {ASSET_SAIDA}' + (' + imagem de clusters' if com_cluster else ''))
    for conta, mods in mx.repartir(modelos, CONTAS).items():
        print(f'  {conta} ({mx.get_project_from_account(conta)}): {", ".join(mods)}')
    if input('Iniciar? [s/N]: ').strip().lower() != 's':
        return
    exportar(conjunto, features, usa_geo, modelos, int(v['n_por_ponto']), cl, com_cluster)
    if original and original != mx.gee_tools.active_user():
        mx.gee_tools.switch_user(original)
        log.info(f'Conta de volta: {original}')


if __name__ == '__main__':
    main()
