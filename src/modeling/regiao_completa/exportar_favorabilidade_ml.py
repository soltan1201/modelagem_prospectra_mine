"""
Mapas de favorabilidade com Random Forest e Gradient Tree Boosting no GEE,
treinados com as amostras do Maxent e salvos por célula da grid.

Treino (uma vez por modelo × alvo × mosaico × conjunto)
  Amostras: PASTA/maxent_amostras_<mosaico>_bandas_indices (models_ML_fosfato.py,
  etapa 1) — presença = todos os pontos do alvo (todos os papéis, ou PAPEIS);
  fundo = sorteio de FATOR_FUNDO × (nº de amostras de presença) pixels de fundo
  (semente SEED), entre os elegíveis. O GEE não tem class_weight; o sorteio deixa
  a proporção perto da usada na avaliação.
  Variáveis (menu "Conjunto"):
    conjuntos de models/variaveis_RC_<sufixo>.json (selecao_variaveis_RC.ipynb,
      célula [21]; sufixo = SUFIXOS[mosaico], ex. 'shm'); padrão = 'recomendado'
    'fixas' = FEATURES abaixo (21 do SHAP do v1, avaliar_maxent.ipynb [34]); nome
      do asset sem conjunto, como nos mapas do v1
  Os nomes são os de calcular_features do models_ML_fosfato.py — o mosaico é
  recalculado igual às amostras.
  Fundo por cluster (variaveis_RC fundo_modo = 'cluster'): models/cluster_fundo_RC_<sufixo>.json
  (notebook, célula [10e]) = k-means do fundo; cada pixel de fundo é padronizado,
  z = (x − média) / desvio, vai para o centro mais próximo e só os clusters_treino
  entram no treino. O nome do asset ganha '_fcl'.

Predição (uma task por célula, para cada task ficar pequena)
  O classificador é o mesmo em todas as células (treina com todos os pontos;
  a semente fixa garante o mesmo modelo). Grid sem solape: a predição é por
  pixel, não precisa de borda.

Saída
  ImageCollection ASSET_SAIDA, uma imagem por
    <modelo>_<alvo>_<mosaico>_<LxCy>                       (conjunto 'fixas')
    <modelo>_<alvo>_<mosaico>_<conjunto>[_fcl]_<LxCy>      (conjuntos do JSON)
    banda 'probability' = p × 10000, Int16
    propriedades modelo, alvo, mosaico, conjunto, fundo_modo, celula, features,
                 n_presenca, n_fundo, n_fundo_elegivel, parametros
  Atenção: exportar_patches.py e visualizar_favorabilidade_ml.js filtram por
  modelo/alvo/mosaico; com mais de um conjunto no mesmo mosaico, filtrar também
  por 'conjunto'.
  Mapa inteiro: ee.ImageCollection(ASSET_SAIDA).filter(modelo, alvo, mosaico).mosaic()
  Registro: registro_favorabilidade_ml_<data>.csv nesta pasta.

Alvos: ALVOS[alvo] = True entram; os outros ficam guardados para depois.

Pula assets que já existem e tasks com o mesmo nome na fila.

Uso: python exportar_favorabilidade_ml.py   (pergunta modelos, mosaicos e
células; mostra o plano e pede confirmação)
"""

import csv
import json
import logging
from datetime import datetime
from pathlib import Path

import ee

import models_ML_fosfato as mx

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CONTAS      = mx.CONTAS   # caatinga02 e caatinga05, cada uma no projeto dela; células repartidas
PASTA       = 'projects/mapbiomas-arida/mine'
ASSET_SAIDA = f'{PASTA}/favorabilidade_ml'
GRID_ASSET  = f'{PASTA}/grid_irece_5x5'          # sem solape
MOSAICOS    = mx.MOSAICOS                         # único, v1 e sazonal_hm_razao (padrão mx.MOSAICO_PADRAO)
ESCALA      = 30
SEED        = 42
FATOR_FUNDO = 3          # amostras de fundo = 3 × amostras de presença
PAPEIS      = None       # None = todos; ['treino'] = sem a validação R2 do fosfato
ESCALA_P    = 10000      # probability × 10000 → Int16
MODELS_DIR  = Path(__file__).resolve().parents[1] / 'models'   # src/modeling/models (compartilhada)
# nome curto do mosaico nos arquivos do selecao_variaveis_RC.ipynb (variaveis_RC_<sufixo>.json)
SUFIXOS     = {'aster-mosaico-final': 'v1', 'aster-mosaico-final_sazonal_hm_razao': 'shm',
               'ASTER_mosaico_unico_final': 'unico'}

# True = exporta agora; False = guardado para depois
ALVOS = {
    'fosfato':   True,
    'carbonato': False,
    'mvt':       False,
    'ouro':      False,
    'manganes':  False,
}

# conjunto 'fixas' — FEATURES_DL do SHAP (avaliar_maxent.ipynb, célula 34;
# ranking_features_shap_20261021.csv): 21 de 36, em ordem de importância, sem pares com |r| > 0,9
FEATURES = [
    'quartz_14_12', 'B10', 'ferric_iron_2_1', 'glcm_contrast', 'B09', 'sultan_g',
    'phosphorite_rbd', 'gossan_4_2', 'sultan_r', 'ferric_ox_4_3', 'SAVI', 'B01', 'B3N',
    'carbonate_13_14', 'kaolinite_7_5', 'carbonate_8_6', 'glcm_hom', 'dolomite',
    'sericite_illite_smectite', 'muscovite_7_6', 'glcm_ent',
]

# parâmetros próximos aos do notebook (RF_PARAMS / GTB_PARAMS)
#   RF : variablesPerSplit = √n (padrão do GEE)
#   GTB: maxNodes 16 folhas ≈ profundidade 4; 80 árvores (mais que ~90 fica pesado no
#        GEE) com shrinkage 0,1 — menos árvores pedem passo maior
MODELOS = {
    'RF':  (ee.Classifier.smileRandomForest,
            dict(numberOfTrees=300, minLeafPopulation=10, bagFraction=0.632, seed=SEED)),
    'GTB': (ee.Classifier.smileGradientTreeBoost,
            dict(numberOfTrees=80, shrinkage=0.1, samplingRate=0.7, maxNodes=16, seed=SEED)),
}


# ═══════════════════════════════════════════════════════════════════════════
# TREINO
# ═══════════════════════════════════════════════════════════════════════════
def asset_amostras(mosaico):
    return f'{PASTA}/maxent_amostras_{mosaico}_bandas_indices'


def ler_json(arq, celula):
    if not arq.exists():
        raise SystemExit(f'{arq.name} não existe — rode até a célula {celula} do selecao_variaveis_RC.ipynb.')
    return json.loads(arq.read_text())


def conjuntos_disponiveis(mosaico):
    """{nome: (features, cluster | None)}; 'fixas' sempre; os do JSON se existir."""
    conj = {'fixas': (FEATURES, None)}
    arq = MODELS_DIR / f'variaveis_RC_{SUFIXOS.get(mosaico, mosaico)}.json'
    if not arq.exists():
        return conj, 'fixas'
    v = json.loads(arq.read_text())
    cl = None
    if v.get('fundo_modo') == 'cluster':
        cl = ler_json(MODELS_DIR / f'cluster_fundo_RC_{SUFIXOS.get(mosaico, mosaico)}.json', '[10e]')
    if v.get('papeis') != PAPEIS:
        log.warning(f'{arq.name}: notebook com papeis = {v.get("papeis")}, aqui PAPEIS = {PAPEIS}')
    conj.update({k: (f, cl) for k, f in v['conjuntos'].items()})
    return conj, v.get('recomendado', 'fixas')


def rotular_fundo(fc, cl):
    """Propriedade 'cluster' = centro mais próximo no espaço padronizado do notebook
    (também usado pelo exportar_favorabilidade_AR.py)."""
    vs, med, dp, cen = cl['variaveis'], cl['media'], cl['desvio'], cl['centros']

    def rotulo(f):
        z = [f.getNumber(v).subtract(med[j]).divide(dp[j]) for j, v in enumerate(vs)]
        dist = []
        for c in cen:
            d = ee.Number(0)
            for j in range(len(vs)):
                d = d.add(z[j].subtract(c[j]).pow(2))
            dist.append(d)
        dist = ee.List(dist)
        return f.set('cluster', dist.indexOf(dist.reduce(ee.Reducer.min())))
    return fc.map(rotulo)


def amostras_treino(mosaico, alvo, features, cl=None):
    """Presença do alvo + fundo sorteado (FATOR_FUNDO × presença) entre os elegíveis
    (todos, ou só os clusters_treino de cl)."""
    fc = ee.FeatureCollection(asset_amostras(mosaico))
    pres = fc.filter(ee.Filter.eq('classe', 1)).filter(ee.Filter.eq('alvo', alvo))
    if PAPEIS:
        pres = pres.filter(ee.Filter.inList('papel', PAPEIS))
    n_pres = pres.size().getInfo()
    if n_pres == 0:
        raise SystemExit(f'Sem amostras de presença de {alvo} em {asset_amostras(mosaico)}.')
    fundo = fc.filter(ee.Filter.eq('classe', 0))
    if cl is not None:
        fundo = rotular_fundo(fundo, cl).filter(ee.Filter.inList('cluster', cl['clusters_treino']))
    n_eleg = fundo.size().getInfo()
    if cl is not None and n_eleg != cl['pixels_treino']:
        log.warning(f'fundo elegível no GEE ({n_eleg}) ≠ notebook ({cl["pixels_treino"]}): '
                    'o asset de amostras mudou depois do notebook?')
    n_fundo = min(FATOR_FUNDO * n_pres, n_eleg)
    fundo = fundo.randomColumn('_sorteio', SEED).sort('_sorteio').limit(n_fundo)
    return pres.merge(fundo).select(features + ['classe']), n_pres, n_fundo, n_eleg


def treinar(nome_modelo, amostras, features):
    construtor, params = MODELOS[nome_modelo]
    return (construtor(**params).setOutputMode('PROBABILITY')
              .train(features=amostras, classProperty='classe', inputProperties=features))


def features_mosaico(mosaico, features):
    """Mesmas features das amostras (calcular_features do models_ML_fosfato.py)."""
    mx.ASSET_MOSAICO = MOSAICOS[mosaico]
    feats, nomes, _ = mx.ler_mosaico()
    faltam = [f for f in features if f not in nomes]
    if faltam:
        raise SystemExit(f'Variáveis fora do mosaico {mosaico}: {faltam} (disponíveis: {nomes})')
    return feats.select(features)


def nome_asset(nome_modelo, alvo, mosaico, conjunto, cl, cel):
    meio = '' if conjunto == 'fixas' else f'_{conjunto}' + ('_fcl' if cl is not None else '')
    return f'{nome_modelo}_{alvo}_{mosaico}{meio}_{cel}'


# ═══════════════════════════════════════════════════════════════════════════
# EXPORT
# ═══════════════════════════════════════════════════════════════════════════
def garantir_colecao():
    try:
        ee.data.getAsset(ASSET_SAIDA)
    except ee.EEException:
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, ASSET_SAIDA)
        log.info(f'Coleção criada: {ASSET_SAIDA}')


def existentes() -> set[str]:
    return {a['name'].rsplit('/', 1)[-1] for a in ee.data.listAssets({'parent': ASSET_SAIDA}).get('assets', [])}


def tasks_ativas() -> set[str]:
    return {t['description'] for t in ee.data.getTaskList() if t['state'] in ('READY', 'RUNNING')}


def exportar(modelos, alvos, mosaicos, celulas, grid, conjunto='fixas'):
    """Células repartidas entre CONTAS (blocos contíguos); cada conta treina e dispara as suas."""
    registro = Path(__file__).resolve().parent / f'registro_favorabilidade_ml_{datetime.now():%Y%m%d_%H%M}.csv'
    linhas = []
    for conta, cels in mx.repartir(celulas, CONTAS).items():
        mx.ativar(conta)
        log.info(f'══ {conta}: {len(cels)} células ({", ".join(cels)})')
        exportar_conta(conta, modelos, alvos, mosaicos, cels, grid, linhas, conjunto)
    with open(registro, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['nome', 'conta', 'task', 'estado'])
        w.writeheader()
        w.writerows(linhas)
    n = sum(l['estado'] == 'iniciada' for l in linhas)
    log.info(f'{n} tasks iniciadas, {len(linhas) - n} puladas → {registro}')


def exportar_conta(conta, modelos, alvos, mosaicos, celulas, grid, linhas, conjunto='fixas'):
    garantir_colecao()
    ja, ativas = existentes(), tasks_ativas()   # tasks: lista da conta/projeto ativo
    for mosaico in mosaicos:
        conj, _ = conjuntos_disponiveis(mosaico)
        if conjunto not in conj:
            raise SystemExit(f'Conjunto {conjunto} não existe para {mosaico} (há: {list(conj)})')
        features, cl = conj[conjunto]
        fundo_modo = 'cluster' if cl is not None else 'aleatorio'
        feats = features_mosaico(mosaico, features)
        for alvo in alvos:
            amostras, n_pres, n_fundo, n_eleg = amostras_treino(mosaico, alvo, features, cl)
            log.info(f'{mosaico} · {alvo} · {conjunto}: presença {n_pres}, fundo {n_fundo} '
                     f'(de {n_eleg} elegíveis, {fundo_modo})')
            for nome_modelo in modelos:
                clf = treinar(nome_modelo, amostras, features)
                mapa = (feats.classify(clf).select(0).multiply(ESCALA_P).round().toInt16()
                             .rename('probability'))
                for cel in celulas:
                    nome = nome_asset(nome_modelo, alvo, mosaico, conjunto, cl, cel)
                    if nome in ja or nome in ativas:
                        linhas.append({'nome': nome, 'conta': conta, 'task': '', 'estado': 'pulado (existe/na fila)'})
                        continue
                    img = mapa.set({'modelo': nome_modelo, 'alvo': alvo, 'mosaico': mosaico, 'celula': cel,
                                    'conjunto': conjunto, 'fundo_modo': fundo_modo,
                                    'clusters_treino': ','.join(map(str, cl['clusters_treino'])) if cl else '',
                                    'features': ','.join(features), 'n_presenca': n_pres, 'n_fundo': n_fundo,
                                    'n_fundo_elegivel': n_eleg,
                                    'fator_fundo': FATOR_FUNDO, 'papeis': str(PAPEIS), 'escala_p': ESCALA_P,
                                    'amostras': asset_amostras(mosaico),
                                    'parametros': str(MODELOS[nome_modelo][1])})
                    task = ee.batch.Export.image.toAsset(
                        image=img, description=nome, assetId=f'{ASSET_SAIDA}/{nome}', region=grid[cel],
                        scale=ESCALA, crs='EPSG:4326', maxPixels=1e13,
                        pyramidingPolicy={'probability': 'mean'})
                    task.start()
                    linhas.append({'nome': nome, 'conta': conta, 'task': task.id, 'estado': 'iniciada'})
                    log.info(f'  {nome}: task {task.id}')


# ═══════════════════════════════════════════════════════════════════════════
def escolher(titulo, opcoes, padrao):
    print(f'\n{titulo}\n' + '\n'.join(f'  {i}) {o}' for i, o in enumerate(opcoes, 1)))
    return opcoes[int(input(f'Opção [{padrao}]: ').strip() or padrao) - 1]


def main():
    original = mx.gee_tools.active_user()
    mx.ativar(CONTAS[0])
    info = ee.FeatureCollection(GRID_ASSET).getInfo()['features']
    grid = {f['properties']['nome']: ee.Geometry(f['geometry']) for f in info}
    todas = sorted(grid)

    alvos = [a for a, ligado in ALVOS.items() if ligado]
    nomes_mod = list(MODELOS)
    m = escolher('Modelos:', nomes_mod + ['os dois'], len(nomes_mod) + 1)
    modelos = nomes_mod if m == 'os dois' else [m]
    nomes_mos = list(MOSAICOS)
    m = escolher('Mosaicos:', nomes_mos + ['todos'], nomes_mos.index(mx.MOSAICO_PADRAO) + 1)
    mosaicos = nomes_mos if m == 'todos' else [m]
    conj, recomendado = conjuntos_disponiveis(mosaicos[0])
    nomes_conj = list(conj)
    conjunto = escolher(f'Conjunto de variáveis ({mosaicos[0]}; "fixas" = 21 do v1):', nomes_conj,
                        nomes_conj.index(recomendado) + 1 if recomendado in nomes_conj else 1)
    features, cl = conj[conjunto]
    txt = input(f'\nCélulas (vazio = todas as {len(todas)}; ex.: L2C2, L3C1): ').strip()
    celulas = todas if not txt else [c.strip().upper() for c in txt.replace(';', ',').replace(' ', ',').split(',') if c.strip()]
    fora = [c for c in celulas if c not in grid]
    if fora:
        raise SystemExit(f'Células inexistentes: {fora} (válidas: {todas})')

    n = len(modelos) * len(alvos) * len(mosaicos) * len(celulas)
    print(f'\nAlvos {alvos} · modelos {modelos} · mosaicos {mosaicos} · {len(celulas)} células'
          f'\nconjunto {conjunto}: {len(features)} variáveis ({", ".join(features)})'
          f'\nfundo {"cluster " + str(cl["clusters_treino"]) if cl else "aleatorio"}, {FATOR_FUNDO} × presença · '
          f'papeis {PAPEIS} · {ESCALA} m'
          f'\nex.: {nome_asset(modelos[0], alvos[0], mosaicos[0], conjunto, cl, celulas[0])}'
          f'\n→ até {n} tasks em {ASSET_SAIDA}')
    if input('Iniciar? [s/N]: ').strip().lower() != 's':
        return
    for conta, cels in mx.repartir(celulas, CONTAS).items():
        print(f'  {conta} ({mx.get_project_from_account(conta)}): {len(cels)} células')
    exportar(modelos, alvos, mosaicos, celulas, grid, conjunto)
    if original and original != mx.gee_tools.active_user():
        mx.gee_tools.switch_user(original)
        log.info(f'Conta de volta: {original}')


if __name__ == '__main__':
    main()
