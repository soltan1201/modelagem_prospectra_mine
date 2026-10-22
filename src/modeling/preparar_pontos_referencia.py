"""
Junta os pontos de referência por alvo num asset só, para o Maxent de cada alvo
(maxent_fosfato.py, opção ALVO).

Saída
  ASSET_SAIDA  (FeatureCollection, pontos) — uma linha por ponto:
    alvo        'fosfato', 'carbonato', 'mvt', 'ouro', 'manganes'
    CODIGO_OBJ  id inteiro único (o do pontos_fosfato antigo quando existir;
                senão BASE_ID + n)
    fonte       só fosfato: 'antigo+cprm' (está nos dois) ou 'cprm'
                (só em acervo_cprm/ocorrencias_fosfato); demais alvos 'cprm'
    grupo       fosfato: bacia_irece / leste_outros / norte; demais ''
    papel       fosfato: 'treino' (bacia_irece) ou 'validacao' (fora da
                bacia, R2); demais 'treino'
    subst, subst_original, toponimia, status, serie, dist_irece
  CSV_SAIDA    a mesma tabela com lon/lat.

Regras
  - fosfato: acervo_cprm/ocorrencias_fosfato (50) + pontos_fosfato antigo
    (30). Antigo e CPRM são o mesmo ponto quando estão a ≤ DIST_PAR_M; aí
    fica a posição da CPRM e fonte = 'antigo+cprm'. Um antigo sem par entra
    com fonte = 'antigo'.
  - demais alvos: acervo_cprm/recursos_minerais, só na área (na_area = 1),
    'subst' padronizada por SUBST_ALVO (as grafias variam: 'Mn', 'Manganês',
    'Manganês (Mn)'…).
  - Pontos do mesmo alvo a ≤ DIST_DUP_M viram um só (séries diferentes do
    acervo repetem a ocorrência); fica o primeiro (antigo+cprm antes de cprm).

Uso: python preparar_pontos_referencia.py   (pergunta antes de exportar)
"""

import logging
import math
import re
from pathlib import Path

import ee
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
MINE              = 'projects/mapbiomas-arida/mine/'
ASSET_FOSFATO_OLD = MINE + 'pontos_fosfato'
ASSET_FOSFATO_CP  = MINE + 'acervo_cprm/ocorrencias_fosfato'
ASSET_RECMIN      = MINE + 'acervo_cprm/recursos_minerais'
ASSET_SAIDA       = MINE + 'pontos_referencia'
CSV_SAIDA         = Path(__file__).resolve().parents[1] / 'Dados' / 'referencias' / 'pontos_referencia.csv'

DIST_PAR_M = 250     # antigo × CPRM: mesmo ponto
DIST_DUP_M = 500     # mesmo alvo: duplicata (= BUFFER_M do Maxent)
BASE_ID    = 900000  # ids novos (não colidem com CODIGO_OBJ do SGB)

# subst padronizada (sem o sufixo ' (xx)', minúsculas) → alvo
SUBST_ALVO = {
    'carbonato': ['calcário', 'calcario', 'calcário dolomítico', 'calcita', 'dolomito',
                  'mármore', 'magnesita', 'ca'],
    'mvt':       ['ba', 'barita', 'bário', 'f', 'fluorita', 'flúor', 'pb', 'chumbo',
                  'zn', 'zinco', 'zinco(cobre+chumbo)'],
    'ouro':      ['au', 'ouro'],
    'manganes':  ['mn', 'manganês', 'manganes', 'mangans'],
}
ALVO_DE = {s: a for a, lista in SUBST_ALVO.items() for s in lista}


# ═══════════════════════════════════════════════════════════════════════════
def dist_m(a, b):
    (x1, y1), (x2, y2) = a, b
    return 6371000 * math.hypot(math.radians(x2 - x1) * math.cos(math.radians((y1 + y2) / 2)),
                                math.radians(y2 - y1))


def baixar(asset, filtro=None):
    fc = ee.FeatureCollection(asset)
    if filtro is not None:
        fc = fc.filter(filtro)
    return [(f['geometry']['coordinates'], f['properties']) for f in fc.getInfo()['features']]


def padronizar(subst):
    return re.sub(r'\s*\([^)]*\)\s*$', '', (subst or '').strip()).lower()


def tirar_duplicatas(linhas):
    """Mantém o primeiro de cada grupo a ≤ DIST_DUP_M (mesmo alvo)."""
    ficam = []
    for l in linhas:
        if all(o['alvo'] != l['alvo'] or dist_m((o['lon'], o['lat']), (l['lon'], l['lat'])) > DIST_DUP_M
               for o in ficam):
            ficam.append(l)
    return ficam


def pontos_fosfato():
    antigos = baixar(ASSET_FOSFATO_OLD)
    cprm = baixar(ASSET_FOSFATO_CP)
    usados, linhas = {}, []
    for c, p in antigos:
        i = min(range(len(cprm)), key=lambda k: dist_m(c, cprm[k][0]))
        if dist_m(c, cprm[i][0]) <= DIST_PAR_M:
            usados[i] = int(p['CODIGO_OBJ'])
        else:
            linhas.append({'lon': c[0], 'lat': c[1], 'CODIGO_OBJ': int(p['CODIGO_OBJ']), 'fonte': 'antigo',
                           'grupo': '', 'toponimia': p.get('TOPONIMIA', ''), 'status': p.get('STATUS_ECO', ''),
                           'serie': '', 'subst_original': p.get('SUBSTANCIA', ''), 'dist_irece': None})
    n_novo = 0
    for i, (c, p) in enumerate(cprm):
        if i in usados:
            cod, fonte = usados[i], 'antigo+cprm'
        else:
            cod, fonte = BASE_ID + n_novo, 'cprm'
            n_novo += 1
        linhas.append({'lon': c[0], 'lat': c[1], 'CODIGO_OBJ': cod, 'fonte': fonte, 'grupo': p['grupo'],
                       'toponimia': p['toponimia'], 'status': p['status'], 'serie': p['serie'],
                       'subst_original': p['subst'], 'dist_irece': p['dist_irece']})
    linhas.sort(key=lambda l: l['fonte'] != 'antigo+cprm')
    for l in linhas:
        l.update(alvo='fosfato', subst='fosforo',
                 papel='treino' if l['grupo'] == 'bacia_irece' else 'validacao')
    log.info(f'fosfato: {len(antigos)} antigos ({len(usados)} com par na CPRM a ≤ {DIST_PAR_M} m) + '
             f'{len(cprm)} CPRM')
    return linhas, BASE_ID + n_novo


def pontos_recmin(prox_id):
    linhas = []
    for c, p in baixar(ASSET_RECMIN, ee.Filter.eq('na_area', 1)):    # getInfo: máx. 5.000
        s = padronizar(p.get('subst'))
        if s not in ALVO_DE:
            continue
        linhas.append({'lon': c[0], 'lat': c[1], 'CODIGO_OBJ': prox_id, 'alvo': ALVO_DE[s], 'fonte': 'cprm',
                       'grupo': '', 'papel': 'treino', 'subst': s, 'subst_original': p['subst'],
                       'toponimia': p.get('toponimia', ''), 'status': p.get('status', ''),
                       'serie': p.get('serie', ''), 'dist_irece': p.get('dist_irece')})
        prox_id += 1
    return linhas


def para_fc(df):
    feats = [ee.Feature(ee.Geometry.Point([r.lon, r.lat]),
                        {k: (None if pd.isna(v) else v) for k, v in r._asdict().items()
                         if k not in ('Index', 'lon', 'lat')})
             for r in df.itertuples()]
    return ee.FeatureCollection(feats)


def main():
    ee.Initialize(project='mapbiomas-brazil')   # cota de tasks do projeto; destino continua em mapbiomas-arida
    fos, prox = pontos_fosfato()
    todas = fos + pontos_recmin(prox)
    antes = pd.Series([l['alvo'] for l in todas]).value_counts()
    df = pd.DataFrame(tirar_duplicatas(todas))
    resumo = pd.DataFrame({'antes': antes, f'sem duplicatas (≤ {DIST_DUP_M} m)': df['alvo'].value_counts()})
    print('\n' + resumo.to_string())
    print('\nfosfato por fonte × papel:')
    print(pd.crosstab(df.loc[df.alvo == 'fosfato', 'fonte'], df.loc[df.alvo == 'fosfato', 'papel']).to_string())

    CSV_SAIDA.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(CSV_SAIDA, index=False)
    log.info(f'CSV → {CSV_SAIDA}')

    if input(f'\nExportar {len(df)} pontos para {ASSET_SAIDA}? [s/N]: ').strip().lower() != 's':
        return
    try:
        ee.data.getAsset(ASSET_SAIDA)
        raise SystemExit(f'{ASSET_SAIDA} já existe — apague ou renomeie antes.')
    except ee.EEException:
        pass
    nome = ASSET_SAIDA.rsplit('/', 1)[-1]
    task = ee.batch.Export.table.toAsset(collection=para_fc(df), description=nome, assetId=ASSET_SAIDA)
    task.start()
    log.info(f'task {task.id} → {ASSET_SAIDA}')


if __name__ == '__main__':
    main()
