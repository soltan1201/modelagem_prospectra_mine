"""
Baixa os mosaicos finais ASTER + índices espectrais para o Google Drive, em
GeoTIFF, por célula da grid 5×5. Versão Python de exportar_mosaicos_indices.js
(mesmos índices, mesmos nomes de arquivo); aqui as tasks já saem iniciadas.

Mosaicos (projects/mapbiomas-arida/mine/)
  ASTER_mosaico_unico_final   etapa 3 do build_mosaico_unico.py
  aster-mosaico-final         mosaico final v1 (build_mosaico_final.py)
  B01–B09 reflectância TOA × 1000; B10–B14 TIR em DN.

Índices (21; os de features_process/calculo_spectrals_index.js + phosphorite_rbd)
  addMineralIndices    ferric_iron_2_1, ferric_ox_4_3, gossan_4_2,
                       ferrous_sil_5_4, sericite_illite_smectite,
                       alunite_kaolinite_pyrophyllite, muscovite_7_6,
                       kaolinite_7_5, dolomite, carbonate_13_14, quartz_14_12,
                       sultan_r/g/b, abrams_g/b
  adicionarIndicesGeo  NDVI, AlOH_Clay, Carbonate (B08/B06), SAVI
  phosphorite_rbd      (ε11 + ε13) / (2·ε12) — Guha et al. 2019; emissividade
                       por normalização (ε_max 0,96) sobre a radiância TOA
  texturas (opcional)  glcm_contrast, glcm_var, glcm_ent, glcm_hom (B3N)
  Repetidos no JS original entram uma vez: abrams_r = sultan_r, Ferric_Fe =
  ferric_iron_2_1. SAVI com reflectância (/1000), não com o valor × 1000.

Arquivos (um por célula e por grupo)
  <mosaico>_<LxCy>_bandas   B01–B14 + bandas de origem do mosaico, Int32
  <mosaico>_<LxCy>_indices  índices (+ texturas), Float32
  A 30 m: bandas ≈ 0,6–1 GB e índices ≈ 1,2 GB por célula; a 90 m ≈ 1/9.

Pula tasks com o mesmo nome que já estão na fila ou rodando.

Uso: python exportar_mosaicos_indices.py   (pergunta mosaicos, grupos, células,
escala e pasta; mostra o plano e pede confirmação antes de iniciar)
"""

import logging

import ee

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
PROJETO_EE = 'mapbiomas-brazil'
MINE = 'projects/mapbiomas-arida/mine/'
MOSAICOS = {
    'ASTER_mosaico_unico_final': {'asset': MINE + 'ASTER_mosaico_unico_final',
                                  'origem': ['quality', 'data', 'camada', 'swir_estimado', 'normalizacao']},
    'aster-mosaico-final':       {'asset': MINE + 'aster-mosaico-final',
                                  'origem': ['origem', 'swir_estimado']},
}
GRID_ASSET = MINE + 'grid_irece_5x5'      # sem solape: os arquivos não se repetem
BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
          'B10', 'B11', 'B12', 'B13', 'B14']
ESCALA_TOA = 1000
ESCALA_PADRAO = 30
PASTA_PADRAO = 'aster_mosaicos_indices'

# TIR: DN → radiância = (DN − 1) × UCC; λ central (µm)
TIR_UCC = {'B10': 6.822e-3, 'B11': 6.780e-3, 'B12': 6.590e-3, 'B13': 5.693e-3, 'B14': 5.225e-3}
TIR_LAMBDA = {'B10': 8.291, 'B11': 8.634, 'B12': 9.075, 'B13': 10.657, 'B14': 11.318}
EPS_MAX = 0.96
C1, C2 = 1.191042e8, 1.4387752e4


# ═══════════════════════════════════════════════════════════════════════════
# ÍNDICES
# ═══════════════════════════════════════════════════════════════════════════
def indices_minerais(img):
    b = {n: img.select(n).toFloat() for n in BANDAS}
    return ee.Image.cat([
        b['B02'].divide(b['B01']).rename('ferric_iron_2_1'),
        b['B04'].divide(b['B3N']).rename('ferric_ox_4_3'),
        b['B04'].divide(b['B02']).rename('gossan_4_2'),
        b['B05'].divide(b['B04']).rename('ferrous_sil_5_4'),
        b['B05'].add(b['B07']).divide(b['B06']).rename('sericite_illite_smectite'),
        b['B04'].add(b['B06']).divide(b['B05']).rename('alunite_kaolinite_pyrophyllite'),
        b['B07'].divide(b['B06']).rename('muscovite_7_6'),
        b['B07'].divide(b['B05']).rename('kaolinite_7_5'),
        b['B06'].add(b['B08']).divide(b['B07']).rename('dolomite'),
        b['B13'].divide(b['B14']).rename('carbonate_13_14'),
        b['B14'].divide(b['B12']).rename('quartz_14_12'),
        b['B04'].divide(b['B07']).rename('sultan_r'),
        b['B04'].divide(b['B01']).rename('sultan_g'),
        b['B02'].divide(b['B3N']).multiply(b['B04'].divide(b['B3N'])).rename('sultan_b'),
        b['B04'].divide(b['B3N']).rename('abrams_g'),
        b['B02'].divide(b['B01']).rename('abrams_b'),
    ])


def indices_geo(img):
    n = img.select('B3N').toFloat().divide(ESCALA_TOA)
    r = img.select('B02').toFloat().divide(ESCALA_TOA)
    return ee.Image.cat([
        img.normalizedDifference(['B3N', 'B02']).rename('NDVI'),
        img.expression("(b('B05') + b('B07')) / (2 * b('B06'))").rename('AlOH_Clay'),
        img.expression("b('B08') / b('B06')").rename('Carbonate'),
        n.subtract(r).multiply(1.5).divide(n.add(r).add(0.5)).rename('SAVI'),
    ]).toFloat()


def fosforito_tir(img):
    """(ε11 + ε13) / (2 ε12), emissividade pelo método de normalização."""
    L = {b: img.select(b).toFloat().subtract(1).max(1).multiply(u) for b, u in TIR_UCC.items()}
    T = [L[b].expression('c2 / (lam * log(eps * c1 / (pow(lam, 5) * L) + 1))',
                         {'c1': C1, 'c2': C2, 'lam': TIR_LAMBDA[b], 'eps': EPS_MAX, 'L': L[b]})
         for b in TIR_UCC]
    tmax = ee.Image.cat(T).reduce(ee.Reducer.max())

    def e(b):
        return L[b].divide(tmax.expression('c1 / (pow(lam, 5) * (exp(c2 / (lam * t)) - 1))',
                                           {'c1': C1, 'c2': C2, 'lam': TIR_LAMBDA[b], 't': tmax}))
    return e('B11').add(e('B13')).divide(e('B12').multiply(2)).toFloat().rename('phosphorite_rbd')


def texturas(img):
    b3n = img.select('B3N').unitScale(0, 0.35 * ESCALA_TOA).multiply(255).clamp(0, 255).toByte()
    return (b3n.glcmTexture(size=1)
               .select(['B3N_contrast', 'B3N_var', 'B3N_ent', 'B3N_idm'],
                       ['glcm_contrast', 'glcm_var', 'glcm_ent', 'glcm_hom']).toFloat())


def imagem_indices(img, com_textura=False):
    out = indices_geo(img).addBands(indices_minerais(img)).addBands(fosforito_tir(img))
    return out.addBands(texturas(img)) if com_textura else out


# ═══════════════════════════════════════════════════════════════════════════
# EXPORT
# ═══════════════════════════════════════════════════════════════════════════
def tasks_ativas() -> set[str]:
    return {t['description'] for t in ee.data.getTaskList() if t['state'] in ('READY', 'RUNNING')}


def planejar(mosaicos, grupos, celulas, com_textura):
    """Lista de (descrição, imagem, geometria) — nada é iniciado aqui."""
    grid = ee.FeatureCollection(GRID_ASSET)
    plano = []
    for cel in celulas:
        geom = grid.filter(ee.Filter.eq('nome', cel)).geometry()
        for m in mosaicos:
            img = ee.ImageCollection(MOSAICOS[m]['asset']).mosaic().clip(geom)
            base = f'{m}_{cel}'
            if 'bandas' in grupos:
                plano.append((f'{base}_bandas', img.select(BANDAS + MOSAICOS[m]['origem']).toInt32(), geom))
            if 'indices' in grupos:
                plano.append((f'{base}_indices', imagem_indices(img, com_textura), geom))
    return plano


def exportar(plano, escala, pasta):
    ativas = tasks_ativas()
    n = 0
    for desc, img, geom in plano:
        if desc in ativas:
            log.info(f'  {desc}: já na fila — pulado')
            continue
        task = ee.batch.Export.image.toDrive(
            image=img, description=desc, fileNamePrefix=desc, folder=pasta, region=geom,
            scale=escala, crs='EPSG:4326', maxPixels=1e13, fileFormat='GeoTIFF',
            formatOptions={'cloudOptimized': True})
        task.start()
        n += 1
        log.info(f'  {desc}: task {task.id}')
    log.info(f'{n} tasks iniciadas → Drive/{pasta}')


# ═══════════════════════════════════════════════════════════════════════════
def escolher(titulo, opcoes, padrao):
    print(f'\n{titulo}\n' + '\n'.join(f'  {i}) {o}' for i, o in enumerate(opcoes, 1)))
    return opcoes[int(input(f'Opção [{padrao}]: ').strip() or padrao) - 1]


def main():
    ee.Initialize(project=PROJETO_EE)
    todas = sorted(ee.FeatureCollection(GRID_ASSET).aggregate_array('nome').getInfo())

    nomes = list(MOSAICOS)
    m = escolher('Mosaicos:', nomes + ['os dois'], len(nomes) + 1)
    mosaicos = nomes if m == 'os dois' else [m]
    g = escolher('O que baixar:', ['bandas + índices', 'só bandas', 'só índices'], 1)
    grupos = {'bandas + índices': ['bandas', 'indices'], 'só bandas': ['bandas'], 'só índices': ['indices']}[g]
    com_textura = 'indices' in grupos and input('\nIncluir texturas GLCM nos índices? [s/N]: ').strip().lower() == 's'
    txt = input(f'\nCélulas (vazio = todas as {len(todas)}; ex.: L2C2, L3C1): ').strip()
    celulas = todas if not txt else [c.strip().upper() for c in txt.replace(';', ',').replace(' ', ',').split(',') if c.strip()]
    fora = [c for c in celulas if c not in todas]
    if fora:
        raise SystemExit(f'Células inexistentes: {fora} (válidas: {todas})')
    escala = int(input(f'Escala em m [{ESCALA_PADRAO}]: ').strip() or ESCALA_PADRAO)
    pasta = input(f'Pasta no Drive [{PASTA_PADRAO}]: ').strip() or PASTA_PADRAO

    plano = planejar(mosaicos, grupos, celulas, com_textura)
    print(f'\n{len(plano)} tasks: {len(celulas)} células × {len(mosaicos)} mosaico(s) × {grupos}'
          f'{" (+ texturas)" if com_textura else ""}, {escala} m → Drive/{pasta}')
    if input('Iniciar? [s/N]: ').strip().lower() != 's':
        return
    exportar(plano, escala, pasta)


if __name__ == '__main__':
    main()
