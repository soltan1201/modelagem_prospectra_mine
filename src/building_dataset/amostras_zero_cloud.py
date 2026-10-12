"""
Amostras de pontos nas cenas zero-cloud (inspecionadas), por cena e por ano,
para a estimativa das bandas SWIR.

  2000–2008 → cenas COM SWIR: colunas VNIR + SWIR + TIR → conjunto de TREINO
              (features VNIR+TIR, alvos B04–B09)
  2009–2020 → cenas SEM SWIR (detector desligado a partir de 2008-S2):
              colunas VNIR + TIR → conjunto de INFERÊNCIA (o modelo preenche o SWIR)

Cada cena entra intacta, igual aos builds: DN → radiância → TOA × 1000
(radiance_to_toa de build_mosaicos_semestrais.py, sem máscara de saturação), TIR
em DN. Uma cena sem SWIR dentro de 2000–2008 (as de 2008-S2) vai para o conjunto
de inferência — o que decide é ter ou não SWIR, não só o ano.

Pixels com VNIR saturado (DN ≥ 254) saem por padrão (REMOVER_SATURADOS): nas
cenas zero-cloud eles se mostraram nuvem (TIR frio), não solo claro — ver
saturacao_landsat_piloto.py.

Lista de cenas: as linhas ativas do visualizar_zero_cloud.js (mesma lista do
build_mosaicos_zero_cloud.py).

Destino (perguntado):
  1) local  → src/Dados/ASTER_Samples_0cloud/{treino_2000_2008,inferencia_2009_2020}/
              um CSV por cena (cache: cena já baixada é pulada) e um CSV por ano
  2) Drive  → uma task Export.table.toDrive por ano e conjunto, na conta ativa

Colunas: B01, B02, B3N, [B04–B09], B10–B14, longitude, latitude, cena, data,
ano, semestre, conjunto. Reflectância = valor / 1000; TIR em DN.

Uso: python amostras_zero_cloud.py
"""

import logging
import sys
from pathlib import Path

import ee
import pandas as pd

import build_mosaicos_semestrais as b
import build_mosaicos_zero_cloud as z

log = logging.getLogger(__name__)

sys.path.append(str(b.PASTA.parent))
from configure_account_projects_ee import get_current_account   # noqa: E402

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
N_POR_CENA        = 2000    # pontos aleatórios por cena (após máscaras)
ESCALA            = 30      # m — resolução nativa do SWIR (VNIR agregado, TIR reamostrado)
SEMENTE           = 42
REMOVER_SATURADOS = True
PASTA_DRIVE       = 'ASTER_Samples_0cloud'
SAIDA             = b.PASTA.parent / 'Dados' / 'ASTER_Samples_0cloud'
CONJUNTOS         = {True: 'treino_2000_2008', False: 'inferencia_2009_2020'}

VNIR, SWIR, TIR = z.VNIR, z.SWIR, b.TIR


# ═══════════════════════════════════════════════════════════════════════════
# AMOSTRAS DE UMA CENA
# ═══════════════════════════════════════════════════════════════════════════
def info_cenas(ids: list[str]) -> pd.DataFrame:
    """Data e presença de SWIR de cada cena (uma chamada ao GEE)."""
    col = ee.ImageCollection([ee.Image(f'ASTER/AST_L1T_003/{i}') for i in ids])
    d = ee.Dictionary({
        'id': col.aggregate_array('system:index'),
        'data': col.aggregate_array('system:time_start'),
        'swir': col.map(lambda i: i.set('tem', i.propertyNames().contains('GAIN_COEFFICIENT_B04')))
                   .aggregate_array('tem'),
    }).getInfo()
    df = pd.DataFrame(d)
    df['data'] = pd.to_datetime(df['data'], unit='ms')
    df['ano'] = df['data'].dt.year
    df['semestre'] = (df['data'].dt.month > 6).astype(int) + 1
    df['conjunto'] = df['swir'].map(CONJUNTOS)
    return df


def amostras_cena(cena: dict) -> ee.FeatureCollection:
    img = ee.Image(f'ASTER/AST_L1T_003/{cena["id"]}')
    bandas = VNIR + (SWIR if cena['swir'] else []) + TIR
    t = ee.Image(b.radiance_to_toa(img, mascarar_saturacao=False))
    valido = t.select(bandas).mask().reduce(ee.Reducer.min())
    if REMOVER_SATURADOS:
        valido = valido.And(t.select('saturado').Not())
    regiao = img.geometry().intersection(ee.Geometry.Polygon(b.LIST_COORD), 100)
    props = {'cena': cena['id'], 'data': cena['data'].strftime('%Y-%m-%d'),
             'ano': int(cena['ano']), 'semestre': int(cena['semestre']), 'conjunto': cena['conjunto']}
    return (t.select(bandas).updateMask(valido)
             .sample(region=regiao, scale=ESCALA, projection=img.select('B01').projection().atScale(ESCALA),
                     numPixels=N_POR_CENA, seed=SEMENTE, geometries=True, tileScale=4)
             .map(lambda f: f.set(props)))


def colunas(com_swir: bool) -> list[str]:
    return (VNIR + (SWIR if com_swir else []) + TIR +
            ['longitude', 'latitude', 'cena', 'data', 'ano', 'semestre', 'conjunto'])


# ═══════════════════════════════════════════════════════════════════════════
# DESTINOS
# ═══════════════════════════════════════════════════════════════════════════
def baixar_local(cenas: pd.DataFrame):
    for conj in cenas['conjunto'].unique():
        (SAIDA / conj / 'cenas').mkdir(parents=True, exist_ok=True)
    for i, cena in enumerate(cenas.to_dict('records'), 1):
        arq = SAIDA / cena['conjunto'] / 'cenas' / f'{cena["id"]}.csv'
        if arq.exists():
            continue
        fc = amostras_cena(cena)
        df = ee.data.computeFeatures({'expression': fc, 'fileFormat': 'PANDAS_DATAFRAME'})
        if len(df):
            coords = df['geo'].apply(lambda g: g['coordinates'])
            df['longitude'], df['latitude'] = coords.str[0], coords.str[1]
            df = df[[c for c in colunas(cena['swir']) if c in df.columns]]
        df.to_csv(arq, index=False)
        log.info(f'[{i}/{len(cenas)}] {cena["id"]} ({cena["conjunto"]}): {len(df)} pontos')

    # um CSV por ano e conjunto, juntando as cenas
    for (conj, ano), grupo in cenas.groupby(['conjunto', 'ano']):
        partes = [pd.read_csv(SAIDA / conj / 'cenas' / f'{c}.csv') for c in grupo['id']]
        partes = [p for p in partes if len(p)]
        if partes:
            arq = SAIDA / conj / f'amostras_0cloud_{ano}.csv'
            pd.concat(partes, ignore_index=True).to_csv(arq, index=False)
            log.info(f'{conj} {ano}: {sum(map(len, partes))} pontos → {arq.name}')


def exportar_drive(cenas: pd.DataFrame):
    for (conj, ano), grupo in cenas.groupby(['conjunto', 'ano']):
        fc = ee.FeatureCollection([amostras_cena(c) for c in grupo.to_dict('records')]).flatten()
        fc = fc.map(lambda f: f.set({'longitude': f.geometry().coordinates().get(0),
                                     'latitude': f.geometry().coordinates().get(1)}))
        nome = f'amostras_0cloud_{conj}_{ano}'
        task = ee.batch.Export.table.toDrive(
            collection=fc, description=nome, folder=PASTA_DRIVE, fileNamePrefix=nome,
            fileFormat='CSV', selectors=colunas(conj == CONJUNTOS[True]))
        task.start()
        log.info(f'task {task.id} → Drive/{PASTA_DRIVE}/{nome}.csv ({len(grupo)} cenas)')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    conta = b.ARQ_CONTA_ATUAL.read_text().strip()
    ee.Initialize(project=get_current_account())
    cenas = info_cenas(z.ler_cenas())

    print('\n══════ Amostras nas cenas zero-cloud ══════')
    print(cenas.groupby('conjunto').agg(cenas=('id', 'size'), anos=('ano', lambda a: f'{a.min()}–{a.max()}'))
               .to_string())
    print(f'{N_POR_CENA} pontos por cena, escala {ESCALA} m, '
          f'saturados {"removidos" if REMOVER_SATURADOS else "mantidos"}.')

    ano_ini = b.perguntar_int('Ano inicial', 2000, 2000, 2020)
    ano_fim = b.perguntar_int('Ano final', 2020, ano_ini, 2020)
    conjuntos = b.perguntar_opcao('Conjuntos:', [
        ('os dois (treino 2000–2008 com SWIR + inferência 2009–2020)', list(CONJUNTOS.values())),
        ('só treino (com SWIR)', [CONJUNTOS[True]]),
        ('só inferência (sem SWIR)', [CONJUNTOS[False]]),
    ])
    destino = b.perguntar_opcao('Destino:', [
        (f'baixar para o PC ({SAIDA})', 'local'),
        (f'tasks para o Drive da conta {conta} (pasta {PASTA_DRIVE})', 'drive'),
    ])
    cenas = cenas[cenas['ano'].between(ano_ini, ano_fim) & cenas['conjunto'].isin(conjuntos)]
    print(f'{len(cenas)} cenas selecionadas.')
    if not len(cenas) or not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Nada a fazer.')

    baixar_local(cenas) if destino == 'local' else exportar_drive(cenas)


if __name__ == '__main__':
    main()
