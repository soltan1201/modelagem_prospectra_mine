"""
Mapas de favorabilidade com Random Forest e Gradient Tree Boosting no GEE,
treinados com as amostras do Maxent e salvos por célula da grid.

Treino (uma vez por modelo × alvo × mosaico)
  Amostras: PASTA/maxent_amostras_<mosaico>_bandas_indices (maxent_fosfato.py,
  etapa 1) — presença = todos os pontos do alvo (todos os papéis, ou PAPEIS);
  fundo = sorteio de FATOR_FUNDO × (nº de amostras de presença) pixels de fundo
  (semente SEED). O GEE não tem class_weight; o sorteio deixa a proporção
  perto da usada na avaliação (avaliar_maxent.ipynb).
  Variáveis: FEATURES (lista abaixo) = FEATURES_DL do SHAP (avaliar_maxent.ipynb,
  célula 34). Os nomes são os de calcular_features
  do maxent_fosfato.py — o mosaico é recalculado igual às amostras.

Predição (uma task por célula, para cada task ficar pequena)
  O classificador é o mesmo em todas as células (treina com todos os pontos;
  a semente fixa garante o mesmo modelo). Grid sem solape: a predição é por
  pixel, não precisa de borda.

Saída
  ImageCollection ASSET_SAIDA, uma imagem por <modelo>_<alvo>_<mosaico>_<LxCy>
    banda 'probability' = p × 10000, Int16
    propriedades modelo, alvo, mosaico, celula, features, n_presenca, n_fundo,
                 parametros
  Mapa inteiro: ee.ImageCollection(ASSET_SAIDA).filter(modelo, alvo, mosaico).mosaic()
  Registro: registro_favorabilidade_ml_<data>.csv nesta pasta.

Alvos: ALVOS[alvo] = True entram; os outros ficam guardados para depois.

Pula assets que já existem e tasks com o mesmo nome na fila.

Uso: python exportar_favorabilidade_ml.py   (pergunta modelos, mosaicos e
células; mostra o plano e pede confirmação)
"""

import csv
import logging
from datetime import datetime
from pathlib import Path

import ee

import maxent_fosfato as mx

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
PROJETO_EE  = 'mapbiomas-brazil'
PASTA       = 'projects/mapbiomas-arida/mine'
ASSET_SAIDA = f'{PASTA}/favorabilidade_ml'
GRID_ASSET  = f'{PASTA}/grid_irece_5x5'          # sem solape
MOSAICOS    = mx.MOSAICOS                         # ASTER_mosaico_unico_final, aster-mosaico-final
ESCALA      = 30
SEED        = 42
FATOR_FUNDO = 3          # amostras de fundo = 3 × amostras de presença
PAPEIS      = None       # None = todos; ['treino'] = sem a validação R2 do fosfato
ESCALA_P    = 10000      # probability × 10000 → Int16

# True = exporta agora; False = guardado para depois
ALVOS = {
    'fosfato':   True,
    'carbonato': False,
    'mvt':       False,
    'ouro':      False,
    'manganes':  False,
}

# FEATURES_DL do SHAP (avaliar_maxent.ipynb, célula 34; ranking_features_shap_20261021.csv):
# 21 de 36, em ordem de importância, sem pares com |r| > 0,9
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


def amostras_treino(mosaico, alvo):
    """Presença do alvo + fundo sorteado (FATOR_FUNDO × presença)."""
    fc = ee.FeatureCollection(asset_amostras(mosaico))
    pres = fc.filter(ee.Filter.eq('classe', 1)).filter(ee.Filter.eq('alvo', alvo))
    if PAPEIS:
        pres = pres.filter(ee.Filter.inList('papel', PAPEIS))
    n_pres = pres.size().getInfo()
    if n_pres == 0:
        raise SystemExit(f'Sem amostras de presença de {alvo} em {asset_amostras(mosaico)}.')
    fundo = (fc.filter(ee.Filter.eq('classe', 0))
               .randomColumn('_sorteio', SEED).sort('_sorteio')
               .limit(FATOR_FUNDO * n_pres))
    n_fundo = fundo.size().getInfo()
    return pres.merge(fundo).select(FEATURES + ['classe']), n_pres, n_fundo


def treinar(nome_modelo, amostras):
    construtor, params = MODELOS[nome_modelo]
    return (construtor(**params).setOutputMode('PROBABILITY')
              .train(features=amostras, classProperty='classe', inputProperties=FEATURES))


def features_mosaico(mosaico):
    """Mesmas features das amostras (calcular_features do maxent_fosfato.py)."""
    mx.ASSET_MOSAICO = MOSAICOS[mosaico]
    feats, nomes, _ = mx.ler_mosaico()
    faltam = [f for f in FEATURES if f not in nomes]
    if faltam:
        raise SystemExit(f'FEATURES fora do mosaico {mosaico}: {faltam} (disponíveis: {nomes})')
    return feats.select(FEATURES)


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


def exportar(modelos, alvos, mosaicos, celulas, grid):
    garantir_colecao()
    ja, ativas = existentes(), tasks_ativas()
    registro = Path(__file__).resolve().parent / f'registro_favorabilidade_ml_{datetime.now():%Y%m%d_%H%M}.csv'
    linhas = []
    for mosaico in mosaicos:
        feats = features_mosaico(mosaico)
        for alvo in alvos:
            amostras, n_pres, n_fundo = amostras_treino(mosaico, alvo)
            log.info(f'{mosaico} · {alvo}: presença {n_pres}, fundo {n_fundo}')
            for nome_modelo in modelos:
                clf = treinar(nome_modelo, amostras)
                mapa = (feats.classify(clf).select(0).multiply(ESCALA_P).round().toInt16()
                             .rename('probability'))
                for cel in celulas:
                    nome = f'{nome_modelo}_{alvo}_{mosaico}_{cel}'
                    if nome in ja or nome in ativas:
                        linhas.append({'nome': nome, 'task': '', 'estado': 'pulado (existe/na fila)'})
                        continue
                    img = mapa.set({'modelo': nome_modelo, 'alvo': alvo, 'mosaico': mosaico, 'celula': cel,
                                    'features': ','.join(FEATURES), 'n_presenca': n_pres, 'n_fundo': n_fundo,
                                    'fator_fundo': FATOR_FUNDO, 'papeis': str(PAPEIS), 'escala_p': ESCALA_P,
                                    'amostras': asset_amostras(mosaico),
                                    'parametros': str(MODELOS[nome_modelo][1])})
                    task = ee.batch.Export.image.toAsset(
                        image=img, description=nome, assetId=f'{ASSET_SAIDA}/{nome}', region=grid[cel],
                        scale=ESCALA, crs='EPSG:4326', maxPixels=1e13,
                        pyramidingPolicy={'probability': 'mean'})
                    task.start()
                    linhas.append({'nome': nome, 'task': task.id, 'estado': 'iniciada'})
                    log.info(f'  {nome}: task {task.id}')
    with open(registro, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['nome', 'task', 'estado'])
        w.writeheader()
        w.writerows(linhas)
    n = sum(l['estado'] == 'iniciada' for l in linhas)
    log.info(f'{n} tasks iniciadas, {len(linhas) - n} puladas → {registro}')


# ═══════════════════════════════════════════════════════════════════════════
def escolher(titulo, opcoes, padrao):
    print(f'\n{titulo}\n' + '\n'.join(f'  {i}) {o}' for i, o in enumerate(opcoes, 1)))
    return opcoes[int(input(f'Opção [{padrao}]: ').strip() or padrao) - 1]


def main():
    ee.Initialize(project=PROJETO_EE)
    info = ee.FeatureCollection(GRID_ASSET).getInfo()['features']
    grid = {f['properties']['nome']: ee.Geometry(f['geometry']) for f in info}
    todas = sorted(grid)

    alvos = [a for a, ligado in ALVOS.items() if ligado]
    nomes_mod = list(MODELOS)
    m = escolher('Modelos:', nomes_mod + ['os dois'], len(nomes_mod) + 1)
    modelos = nomes_mod if m == 'os dois' else [m]
    nomes_mos = list(MOSAICOS)
    m = escolher('Mosaicos:', nomes_mos + ['os dois'], len(nomes_mos) + 1)
    mosaicos = nomes_mos if m == 'os dois' else [m]
    txt = input(f'\nCélulas (vazio = todas as {len(todas)}; ex.: L2C2, L3C1): ').strip()
    celulas = todas if not txt else [c.strip().upper() for c in txt.replace(';', ',').replace(' ', ',').split(',') if c.strip()]
    fora = [c for c in celulas if c not in grid]
    if fora:
        raise SystemExit(f'Células inexistentes: {fora} (válidas: {todas})')

    n = len(modelos) * len(alvos) * len(mosaicos) * len(celulas)
    print(f'\nAlvos {alvos} · modelos {modelos} · mosaicos {mosaicos} · {len(celulas)} células'
          f'\n{len(FEATURES)} features · fundo = {FATOR_FUNDO} × presença · {ESCALA} m'
          f'\n→ até {n} tasks em {ASSET_SAIDA}')
    if input('Iniciar? [s/N]: ').strip().lower() != 's':
        return
    exportar(modelos, alvos, mosaicos, celulas, grid)


if __name__ == '__main__':
    main()
