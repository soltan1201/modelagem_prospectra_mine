"""
Mosaicos semestrais ASTER (Irecê, BA) — todos os semestres de uma vez.

Versão Python (earthengine-api) de select_black_list_save_mosaicSem.js, com as
correções já aplicadas:
  - DN → radiância com GAIN_COEFFICIENT_Bxx antes do TOA (VNIR/SWIR);
  - pixel mascarado quando B01/B02/B3N saturam (DN ≥ 254);
  - TIR segue em DN, exportado ×10 (escala por banda);
  - exportação em INT32 (sem risco de estouro do teto 32767 do INT16);
  - limiares de brilho reescalados para o TOA com ganho;
  - ramos complementares: cena da lista zero-cloud → sem máscara, quality 2,5;
    demais → máscara de nuvem/sombra + quality v5 (máximo 2,0).

Listas lidas dos arquivos de rastreamento (mesma pasta deste script):
  lista_imagens_Notaceites.txt → cenas excluídas (blacklist)
  list_zero_cloud.txt          → cenas sem nuvem (sem máscara, quality 2,5)
Usa a união de todos os IDs de cada arquivo, ignorando o nome do bloco: o
filtro de data de cada semestre já seleciona as cenas certas. Cena presente
nas duas listas é excluída (a blacklist vence).

Versão 3 — exportação por célula da grid 5 × 5 (GRID_ASSET):
  um asset por ano × semestre × célula, nome
  ASTER_QualityMosaic_{ano}_semestre_{sem}_{L<linha>C<coluna>}_v3, com as
  propriedades celula_id / celula_nome. Cada célula usa só as cenas que a tocam
  e é recortada pela geometria da célula (grid com sobreposição de 100 m).
  Todas as células saem na MESMA grade de pixels (crs EPSG:4326 + CRS_TRANSFORM
  ancorado no canto NW do bloco), então encaixam pixel a pixel no mosaic().

Conta e destino:
  conta    → superconta (projeto mapbiomas-brazil); a conta ativa antes de
             rodar é restaurada no fim
  destino  → DESTINO (ImageCollection v3; todos os anos e células)
Antes de exportar, o script lista o que já está salvo no destino, consulta as
tasks da conta e mostra o plano para confirmação:
  novo         → exporta
  existe       → pula (ou apaga e refaz, se sobrescrever = sim)
  processando  → task READY/RUNNING na conta ou em TASKS_PROCESSANDO: nunca reenvia
Itens cuja última task falhou aparecem com o erro e são reenviados.
O limite de tasks por execução deixa o restante para a próxima rodada.

Uso (interativo — o script pergunta anos, semestres e modo):
  python build_mosaicos_semestrais.py
Enter aceita o valor padrão mostrado entre colchetes.
"""

import csv
import logging
import math
import re
import sys
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import ee

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

PASTA = Path(__file__).resolve().parent
sys.path.append(str(PASTA.parent))   # src/: configure_account_projects_ee, gee_tools
from configure_account_projects_ee import get_project_from_account
from gee_tools import switch_user

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CONTA_EXPORT = 'superconta'
DESTINO      = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA/ASTER_v3'
GRID_ASSET   = 'projects/mapbiomas-arida/mine/grid_irece_5x5_solape100'
VERSAO       = 'v3'
ARQ_CONTA_ATUAL = Path.home() / '.config' / 'earthengine' / 'current_user.txt'

# Tasks ainda em execução no GEE que NÃO devem ser reenviadas. O script também
# consulta sozinho as tasks READY/RUNNING da conta; esta lista é um reforço
# manual (ex.: tasks iniciadas por outra conta). Limpe quando terminarem.
TASKS_PROCESSANDO: list[str] = [
    # ex.: 'ASTER_QualityMosaic_2004_semestre_2_L3C2_v3',
]

LIST_COORD = [
    [-44.460483519818084, -14.790793773956535],
    [-39.203525512005584, -14.790793773956535],
    [-39.203525512005584, -9.416525285369607],
    [-44.460483519818084, -9.416525285369607],
    [-44.460483519818084, -14.790793773956535],
]
CC_MAX       = 70          # filtro de metadado CLOUDCOVER < 70

# Nuvem / sombra (limiares de brilho já reescalados para TOA com ganho)
LIM_VIS_BRIGHT       = 0.19   # era 0.25 com TOA sem ganho (÷1,31)
RAIO_BORDA_NUVEM     = 180    # m
ALTURA_NUVEM         = 2000   # m
LIM_ESCURO_GREEN     = 0.12   # era 0.18 (÷1,49)
LIM_ESCURO_NIR       = 0.17   # era 0.20 (÷1,17)
DN_SATURADO          = 254    # VNIR 8 bits: DN ≥ 254 = saturado

# Quality
QUALITY_ZERO_CLOUD   = 2.5    # cenas da lista zero-cloud (sempre vencem)
QUALITY_MAX_DEMAIS   = 2.0    # teto do quality v5

ESUN = {'B01': 1848, 'B02': 1549, 'B3N': 1114,
        'B04': 225.4, 'B05': 86.63, 'B06': 81.85,
        'B07': 74.85, 'B08': 66.49, 'B09': 59.85}
TIR = ['B10', 'B11', 'B12', 'B13', 'B14']

# Exportação INT32 (valores = físico × fator)
BANDAS_EXPORT = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
                 'B10', 'B11', 'B12', 'B13', 'B14', 'quality']
ESCALA_EXPORT_FATOR = {**dict.fromkeys(list(ESUN), 10000),   # reflectância 0–1 → 0–10000
                **dict.fromkeys(TIR, 10),              # TIR em DN (12 bits) → ×10
                'quality': 10000}                      # 0–2,5 → 0–25000
ESCALA_EXPORT = 15   # m (resolução VNIR)

# Grade de pixels comum a todas as células: EPSG:4326, pixel de 15 m em graus
# (mesma conversão que o GEE usa para scale=15 em EPSG:4326), origem no canto NW
# do bloco. Com ela as células vizinhas compartilham os mesmos pixels.
CRS_EXPORT     = 'EPSG:4326'
RES_GRAUS      = ESCALA_EXPORT / 111319.49079327357
OESTE_BLOCO    = min(x for x, _ in LIST_COORD)
NORTE_BLOCO    = max(y for _, y in LIST_COORD)
CRS_TRANSFORM  = [RES_GRAUS, 0, OESTE_BLOCO, 0, -RES_GRAUS, NORTE_BLOCO]
LIMITE_TASKS_PADRAO = 100


# ═══════════════════════════════════════════════════════════════════════════
# LISTAS DE IMAGENS
# ═══════════════════════════════════════════════════════════════════════════
def ler_ids(arquivo: str) -> set[str]:
    """Todos os timestamps ASTER (system:index) do arquivo, sem olhar o bloco."""
    return set(re.findall(r'AST_L1T_003/(\d{14})', (PASTA / arquivo).read_text()))


def ids_do_semestre(ids: set[str], ano: int, semestre: int) -> ee.List:
    """Só os IDs (AAAAMMDD...) do semestre: deixa o grafo de cada task bem menor
    do que mandar as ~2.500 cenas da blacklist inteira."""
    meses = range(1, 7) if semestre == 1 else range(7, 13)
    return ee.List(sorted(i for i in ids if int(i[:4]) == ano and int(i[4:6]) in meses))


# ═══════════════════════════════════════════════════════════════════════════
# FUNÇÕES DE PROCESSAMENTO (port de select_black_list_save_mosaicSem.js)
# ═══════════════════════════════════════════════════════════════════════════
def radiance_to_toa(img):
    """DN → radiância [(DN − 1) × ganho] → reflectância TOA (VNIR/SWIR); TIR segue em DN."""
    img = ee.Image(img)
    cosz = ee.Number(90).subtract(ee.Number(img.get('SOLAR_ELEVATION'))) \
                        .multiply(math.pi / 180).cos()
    # mesma fórmula de distância Terra-Sol do script JS
    doy = ee.Number(ee.Date(img.get('system:time_start')).getRelative('day', 'year'))
    d = ee.Number(1).subtract(ee.Number(0.01672).multiply(doy.multiply(2 * math.pi / 365).cos()))

    # A partir de 2008-S2 o SWIR do ASTER está desligado: B04–B09 vêm vazias e sem
    # GAIN_COEFFICIENT. Nesse caso as bandas SWIR saem totalmente mascaradas.
    tem_swir = img.propertyNames().contains('GAIN_COEFFICIENT_B04')

    def toa(b, esun):
        dn = img.select(b)
        ganho = ee.Number(img.get(f'GAIN_COEFFICIENT_{b}'))
        return (dn.subtract(1).multiply(ganho)
                  .multiply(math.pi).multiply(d.pow(2))
                  .divide(cosz.multiply(esun))
                  .updateMask(dn.gt(0))
                  .rename(b).toFloat())

    opticas = []
    for b, esun in ESUN.items():
        if b in ('B01', 'B02', 'B3N'):
            opticas.append(toa(b, esun))
        else:
            vazia = ee.Image.constant(0).toFloat().updateMask(0).rename(b)
            opticas.append(ee.Image(ee.Algorithms.If(tem_swir, toa(b, esun), vazia)))

    nao_saturado = img.select(['B01', 'B02', 'B3N']).reduce(ee.Reducer.max()).lt(DN_SATURADO)
    # copyProperties devolve Element no Python: ee.Image(...) devolve o tipo certo
    return ee.Image(ee.Image.cat(opticas)
              .addBands(img.select(TIR).toFloat())
              .updateMask(nao_saturado)
              .copyProperties(img, img.propertyNames())
              .set('TEM_SWIR', tem_swir))


def add_indices(img):
    """SAVI e NDWI (McFeeters) a partir do TOA."""
    img = ee.Image(img)
    green, red, nir = img.select('B01'), img.select('B02'), img.select('B3N')
    savi = nir.subtract(red).multiply(1.5).divide(nir.add(red).add(0.5)).rename('SAVI')
    ndwi = green.subtract(nir).divide(green.add(nir).add(1e-6)).rename('NDWI')
    return img.addBands([savi, ndwi])


def mascara_nuvem(img):
    """Nuvem (brilho visível + TIR frio relativo) e sombra (escuro perto de nuvem)."""
    img = ee.Image(img)
    green, red, nir = img.select('B01'), img.select('B02'), img.select('B3N')
    tir10 = img.select('B10')

    tan_elev = ee.Number(img.get('SOLAR_ELEVATION')).multiply(math.pi / 180).tan()
    shadow_dist = ee.Number(ALTURA_NUVEM).divide(tan_elev.max(0.1)).max(200).min(2000)

    vis_bright = green.add(red).add(nir).divide(3)
    p20 = ee.Number(tir10.reduceRegion(
        reducer=ee.Reducer.percentile([20]), geometry=img.geometry(),
        scale=180, maxPixels=1e10, bestEffort=True).get('B10'))
    cloud = vis_bright.gt(LIM_VIS_BRIGHT).And(tir10.lt(p20))
    cloud = (cloud.focalMin(radius=RAIO_BORDA_NUVEM, kernelType='square', units='meters')
                  .focalMax(radius=RAIO_BORDA_NUVEM, kernelType='square', units='meters')
                  .rename('cloud'))

    dark = green.lt(LIM_ESCURO_GREEN).And(nir.lt(LIM_ESCURO_NIR))
    shadow = dark.And(cloud.focalMax(radius=shadow_dist.add(300), units='meters'))
    shadow = shadow.focalMax().rename('shadow')

    # pixel válido: VNIR > 0 e, quando a cena tem SWIR, B04 > 0
    # (sem SWIR, exigir B04 apagaria a cena inteira)
    vnir_ok = img.select(['B01', 'B02', 'B3N']).reduce(ee.Reducer.min()).gt(0)
    swir_ok = ee.Image(ee.Algorithms.If(img.get('TEM_SWIR'), img.select('B04').gt(0), ee.Image(1)))
    valid = vnir_ok.And(swir_ok)
    clear = valid.And(cloud.Not()).And(shadow.Not())
    return ee.Image(add_indices(img).updateMask(clear).copyProperties(img, img.propertyNames()))


def add_quality_v5(img):
    """Quality v5 (VNIR + TIR), teto QUALITY_MAX_DEMAIS."""
    img = ee.Image(img)
    b1, b2, b3n, b10 = img.select('B01'), img.select('B02'), img.select('B3N'), img.select('B10')

    savi = b3n.subtract(b2).multiply(1.5).divide(b3n.add(b2).add(0.5).add(1e-6))
    ndwi = b1.subtract(b3n).divide(b1.add(b3n).add(1e-6))
    brightness = b1.add(b2).add(b3n).divide(3)

    tir_score = b10.unitScale(1050, 1300).clamp(0, 1)          # TIR em DN
    bright_up = brightness.unitScale(0.053, 0.122).clamp(0, 1)  # era (0.07, 0.16)
    bright_down = ee.Image(1).subtract(brightness.unitScale(0.168, 0.32).clamp(0, 1))  # era (0.22, 0.42)
    bright_score = bright_up.multiply(bright_down)
    savi_score = savi.unitScale(-0.1, 0.5).clamp(0, 1)
    ndwi_score = ndwi.multiply(-1).unitScale(-0.10, 0.45).clamp(0, 1)

    base = (tir_score.multiply(0.45).add(bright_score.multiply(0.25))
            .add(savi_score.multiply(0.15)).add(ndwi_score.multiply(0.15)).clamp(0, 1))

    cc = ee.Number(img.get('CLOUDCOVER'))
    cloud_bonus = ee.Number(1).subtract(cc.divide(80).min(1)).multiply(0.8)
    perfect_bonus = cc.eq(0).multiply(0.3)
    quality = base.add(cloud_bonus).add(perfect_bonus).clamp(0, QUALITY_MAX_DEMAIS)

    quality = quality.updateMask(b1.mask().And(brightness.gt(0.005)))
    return img.addBands(quality.rename('quality'))


def add_quality_zero_cloud(img):
    """Cena sem nuvem: sem máscara, quality fixo (sempre vence no qualityMosaic)."""
    img = ee.Image(img)
    quality = img.select('B01').gt(0).multiply(QUALITY_ZERO_CLOUD).rename('quality')
    return add_indices(img).addBands(quality)


def para_int32(img):
    """Escala cada banda pelo seu fator e converte para INT32."""
    bandas = [img.select(b).multiply(ESCALA_EXPORT_FATOR[b]).round().toInt32().rename(b)
              for b in BANDAS_EXPORT]
    return ee.Image(ee.Image.cat(bandas)
              .copyProperties(img, img.propertyNames())
              .set({'SCALE_FACTORS': ESCALA_EXPORT_FATOR, 'DATA_TYPE': 'INT32'}))


# ═══════════════════════════════════════════════════════════════════════════
# MOSAICO DE UM SEMESTRE
# ═══════════════════════════════════════════════════════════════════════════
def montar_semestre(ano: int, semestre: int, area, blacklist: ee.List, zero_cloud: ee.List):
    # semestre 1 = [ano-01-01, ano-07-01); semestre 2 = [ano-07-01, (ano+1)-01-01) — fim exclusivo
    inicio = ee.Date.fromYMD(ano, 1 if semestre == 1 else 7, 1)
    base = (ee.ImageCollection('ASTER/AST_L1T_003')
              .filterDate(inicio, inicio.advance(6, 'month'))
              .filterBounds(area)
              .filter(ee.Filter.lt('CLOUDCOVER', CC_MAX)))
    # cenas noturnas (sol abaixo do horizonte) só têm TIR: VNIR vazio e sem ganho
    diurnas = base.filter(ee.Filter.And(ee.Filter.gt('SOLAR_ELEVATION', 0),
                                        ee.Filter.notNull(['GAIN_COEFFICIENT_B01'])))
    filtrada = diurnas.filter(ee.Filter.inList('system:index', blacklist).Not())

    na_lista = ee.Filter.inList('system:index', zero_cloud)
    cc0 = filtrada.filter(na_lista).map(radiance_to_toa).map(add_quality_zero_cloud)
    nuvens = (filtrada.filter(na_lista.Not())
                .map(radiance_to_toa).map(mascara_nuvem).map(add_quality_v5))

    colecao = cc0.merge(nuvens)
    mosaico = colecao.qualityMosaic('quality').clip(area)
    mosaico = mosaico.updateMask(mosaico.select('B01').gt(0))
    mosaico = ee.Image(para_int32(mosaico).set({
        'year': ano, 'semestre': str(semestre), 'bloco': 'Irece',
        'num_images': colecao.size(),
        'num_zero_cloud': cc0.size(),
        'versao': VERSAO, 'processamento': 'ganho+saturacao+escala_por_banda+int32',
    }))
    contagens = {'total_cc_lt_70': base.size(), 'diurnas': diurnas.size(),
                 'apos_blacklist': filtrada.size(),
                 'zero_cloud': cc0.size(), 'com_mascara': nuvens.size()}
    return mosaico, contagens


# ═══════════════════════════════════════════════════════════════════════════
# CONTA E DESTINO
# ═══════════════════════════════════════════════════════════════════════════
def nome_mosaico(ano: int, semestre: int, celula: str) -> str:
    return f'ASTER_QualityMosaic_{ano}_semestre_{semestre}_{celula}_{VERSAO}'


def ler_grid() -> dict[int, dict]:
    """Células da grid: {id: {'nome', 'geom'}}. Falha cedo se a conta não lê o asset."""
    try:
        info = ee.FeatureCollection(GRID_ASSET).getInfo()
    except ee.EEException as e:
        raise SystemExit(f'Não foi possível ler a grid {GRID_ASSET} com {CONTA_EXPORT}: {e}\n'
                         'Compartilhe o asset com a conta ou copie a grid para o projeto dela.')
    grid = {int(f['properties']['id']): {'nome': f['properties']['nome'],
                                          'geom': ee.Geometry(f['geometry'], None, False)}
            for f in info['features']}
    log.info(f'Grid {GRID_ASSET}: {len(grid)} células.')
    return grid


def usar_conta(nome: str) -> None:
    """Troca as credenciais para `nome` e inicializa o EE no projeto da conta."""
    switch_user(nome)
    projeto = get_project_from_account(nome)
    ee.Initialize(project=projeto)
    log.info(f'Earth Engine inicializado com {nome} (projeto {projeto}).')


def listar_destino() -> dict[str, dict]:
    """Assets já salvos no destino: {nome: info do listAssets}. Cria a coleção se faltar."""
    try:
        info = ee.data.getAsset(DESTINO)
        log.info(f'Destino {DESTINO} existe (tipo {info["type"]}).')
    except ee.EEException:
        log.info(f'Destino {DESTINO} não existe — criando ImageCollection.')
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, DESTINO)
        return {}
    existentes, token = {}, None
    while True:
        params = {'parent': DESTINO, **({'pageToken': token} if token else {})}
        resp = ee.data.listAssets(params)
        for a in resp.get('assets', []):
            existentes[a['id'].rsplit('/', 1)[-1]] = a
        token = resp.get('nextPageToken')
        if not token:
            return existentes


def estado_tasks() -> dict[str, dict]:
    """Última task de cada descrição ASTER_QualityMosaic_*_v3 na conta (getTaskList vem
    da mais recente para a mais antiga)."""
    estados = {}
    for t in ee.data.getTaskList():
        desc = t.get('description', '')
        if (desc.startswith('ASTER_QualityMosaic_') and desc.endswith(f'_{VERSAO}')
                and desc not in estados):
            estados[desc] = t
    return estados


def planejar(args, grid: dict[int, dict], existentes: dict[str, dict],
             tasks: dict[str, dict]) -> dict[tuple, str]:
    """Classifica cada ano × semestre × célula: existe / processando / novo / refazer.
    Mostra o resumo por semestre."""
    ativas = {d for d, t in tasks.items() if t['state'] in ('READY', 'RUNNING')} | set(TASKS_PROCESSANDO)

    print(f'\n══════ Destino: {DESTINO} ══════')
    print(f'{len(existentes)} asset(s) já salvo(s); {len(ativas)} task(s) em processamento; '
          f'células pedidas: {len(args.celulas)}.')

    plano = {}
    for ano in range(args.anos[0], args.anos[1] + 1):
        for sem in args.semestres:
            for cid in args.celulas:
                nome = nome_mosaico(ano, sem, grid[cid]['nome'])
                if nome in ativas:
                    plano[(ano, sem, cid)] = 'processando'          # nunca reenviar
                elif nome in existentes:
                    plano[(ano, sem, cid)] = 'refazer' if args.sobrescrever else 'existe'
                else:
                    plano[(ano, sem, cid)] = 'novo'

    rotulos = {'novo': 'novas', 'refazer': 'a APAGAR e refazer',
               'existe': 'já salvas', 'processando': 'processando'}
    semestres = sorted({(a, s) for a, s, _ in plano})
    for a, s in semestres:
        cats = [c for (aa, ss, _), c in plano.items() if (aa, ss) == (a, s)]
        partes = [f'{cats.count(c)} {r}' for c, r in rotulos.items() if cats.count(c)]
        print(f'  {a}-S{s}: {", ".join(partes)}')
    total = {r: sum(c == k for c in plano.values()) for k, r in rotulos.items()}
    print('  TOTAL: ' + ', '.join(f'{n} {r}' for r, n in total.items() if n))

    falhas = [(n, t.get('error_message', '')) for n, t in tasks.items()
              if t['state'] == 'FAILED' and n not in existentes and n not in ativas]
    if falhas:
        print(f'  última task FALHOU em {len(falhas)} item(ns) (reenviados se pedidos):')
        for n, erro in sorted(falhas)[:20]:
            print(f'    {n}: {" ".join(erro.split())[:90]}')
        if len(falhas) > 20:
            print(f'    … e mais {len(falhas) - 20}')
    return plano


# ═══════════════════════════════════════════════════════════════════════════
# PERGUNTAS AO USUÁRIO
# ═══════════════════════════════════════════════════════════════════════════
def perguntar_int(texto: str, padrao: int, minimo: int, maximo: int) -> int:
    while True:
        resp = input(f'{texto} [{padrao}]: ').strip()
        if not resp:
            return padrao
        if resp.isdigit() and minimo <= int(resp) <= maximo:
            return int(resp)
        print(f'  → digite um número entre {minimo} e {maximo}.')


def perguntar_opcao(texto: str, opcoes: list[tuple[str, object]], padrao: int = 1):
    """Mostra opções numeradas e devolve o valor da escolhida."""
    print(texto)
    for i, (rotulo, _) in enumerate(opcoes, 1):
        print(f'  {i}) {rotulo}')
    escolha = perguntar_int('Opção', padrao, 1, len(opcoes))
    return opcoes[escolha - 1][1]


def perguntar_sim_nao(texto: str, padrao: bool = False) -> bool:
    sufixo = '[S/n]' if padrao else '[s/N]'
    while True:
        resp = input(f'{texto} {sufixo}: ').strip().lower()
        if not resp:
            return padrao
        if resp in ('s', 'sim', 'n', 'nao', 'não'):
            return resp.startswith('s')
        print('  → responda s ou n.')


def perguntar_celulas(grid: dict[int, dict]) -> list[int]:
    """'' ou 'todas' → todas; senão lista/intervalos de ids, ex.: 1,2,7 ou 11-15."""
    validos = sorted(grid)
    while True:
        resp = input(f'Células (ids {validos[0]}–{validos[-1]}; ex.: 1,2,7 ou 11-15) [todas]: ').strip()
        if resp in ('', 'todas'):
            return validos
        ids = set()
        try:
            for parte in resp.replace(' ', '').split(','):
                if '-' in parte:
                    a, b = map(int, parte.split('-'))
                    ids.update(range(a, b + 1))
                else:
                    ids.add(int(parte))
        except ValueError:
            print('  → formato inválido.')
            continue
        fora = ids - set(validos)
        if fora:
            print(f'  → ids inexistentes na grid: {sorted(fora)}')
            continue
        return sorted(ids)


def coletar_opcoes(grid: dict[int, dict]) -> SimpleNamespace:
    print('\n══════ Mosaicos semestrais ASTER v3 (por célula) — configuração ══════')
    ano_ini = perguntar_int('Ano inicial', 2000, 2000, 2021)
    ano_fim = perguntar_int('Ano final', max(ano_ini, 2020), ano_ini, 2021)

    semestres = perguntar_opcao('Semestres:', [
        ('os dois (jan–jun e jul–dez)', [1, 2]),
        ('só o 1º (jan–jun)', [1]),
        ('só o 2º (jul–dez)', [2]),
    ])
    celulas = perguntar_celulas(grid)

    dry_run = perguntar_opcao('Modo:', [
        ('só contagens de cenas (dry-run, não exporta)', True),
        ('exportar os mosaicos', False),
    ])

    sobrescrever, limite = False, None
    if not dry_run:
        limite = perguntar_int('Máximo de tasks a disparar nesta execução',
                               LIMITE_TASKS_PADRAO, 1, 3000)
        sobrescrever = perguntar_sim_nao('Se o asset já existir, apagar e exportar de novo?')
        if sobrescrever:
            print('  ATENÇÃO: os assets existentes nesses itens serão APAGADOS (não tem volta).')
            sobrescrever = input('  Digite APAGAR para confirmar: ').strip() == 'APAGAR'
            if not sobrescrever:
                print('  Confirmação não recebida — assets existentes serão pulados.')

    return SimpleNamespace(anos=[ano_ini, ano_fim], semestres=semestres, celulas=celulas,
                           sobrescrever=sobrescrever, dry_run=dry_run, limite=limite)


# ═══════════════════════════════════════════════════════════════════════════
# EXECUÇÃO
# ═══════════════════════════════════════════════════════════════════════════
def main():
    conta_original = ARQ_CONTA_ATUAL.read_text().strip()
    try:
        usar_conta(CONTA_EXPORT)
        grid = ler_grid()
        args = coletar_opcoes(grid)
        existentes = listar_destino()
        plano = planejar(args, grid, existentes, estado_tasks())
        print(f'\nModo: {"dry-run (só contagens)" if args.dry_run else "EXPORTAR"} | '
              f'conta: {CONTA_EXPORT} | sobrescrever={args.sobrescrever}'
              + ('' if args.dry_run else f' | limite de tasks={args.limite}'))
        if not perguntar_sim_nao('Continuar?', padrao=True):
            raise SystemExit('Cancelado.')
        executar(args, grid, plano)
    finally:
        if conta_original != CONTA_EXPORT:
            log.info(f'Voltando para a conta original: {conta_original}')
            switch_user(conta_original)


def executar(args, grid: dict[int, dict], plano: dict[tuple, str]):
    ids_blacklist = ler_ids('lista_imagens_Notaceites.txt')
    ids_zero = ler_ids('list_zero_cloud.txt')
    conflito = ids_blacklist & ids_zero
    if conflito:
        log.warning(f'{len(conflito)} cena(s) nas duas listas — ficam excluídas: {sorted(conflito)}')
    ids_zero -= ids_blacklist
    log.info(f'Blacklist: {len(ids_blacklist)} cenas | zero-cloud: {len(ids_zero)} cenas')

    def base_registro(ano, sem, cid):
        nome = nome_mosaico(ano, sem, grid[cid]['nome'])
        return {'ano': ano, 'semestre': sem, 'celula_id': cid, 'celula_nome': grid[cid]['nome'],
                'asset': f'{DESTINO}/{nome}'}

    registro = []
    for (ano, sem, cid), cat in plano.items():
        if cat in ('existe', 'processando') and not args.dry_run:
            registro.append({**base_registro(ano, sem, cid), 'status': cat})
    pendentes = [k for k, c in plano.items() if args.dry_run or c in ('novo', 'refazer')]
    if not pendentes:
        log.info('Nada a exportar.')

    disparadas = 0
    for ano, sem in sorted({(a, s) for a, s, _ in pendentes}):
        if not args.dry_run and disparadas >= args.limite:
            break
        blacklist_sem = ids_do_semestre(ids_blacklist, ano, sem)
        zero_sem = ids_do_semestre(ids_zero, ano, sem)
        celulas = [cid for a, s, cid in pendentes if (a, s) == (ano, sem)]

        # monta as células do semestre e busca as contagens numa única chamada
        montados = {cid: montar_semestre(ano, sem, grid[cid]['geom'], blacklist_sem, zero_sem)
                    for cid in celulas}
        contagens = ee.List([ee.Dictionary(c) for _, c in montados.values()]).getInfo()

        for (cid, (mosaico, _)), cont in zip(montados.items(), contagens):
            linha = {**base_registro(ano, sem, cid), **cont}
            nome = linha['asset'].rsplit('/', 1)[-1]
            cat = plano[(ano, sem, cid)]
            log.info(f'{ano}-S{sem} {grid[cid]["nome"]} [{cat}]: {cont}')

            if cont['apos_blacklist'] == 0:
                registro.append({**linha, 'status': 'sem_cenas'})
                continue
            if args.dry_run:
                registro.append({**linha, 'status': f'dry_run ({cat})'})
                continue
            if disparadas >= args.limite:
                registro.append({**linha, 'status': 'adiado (limite de tasks)'})
                continue

            mosaico = ee.Image(mosaico.set({'celula_id': cid, 'celula_nome': grid[cid]['nome'],
                                            'grid_asset': GRID_ASSET}))
            # Monta a task antes de apagar: se a montagem falhar, o asset antigo fica intacto
            task = ee.batch.Export.image.toAsset(
                image=mosaico, description=nome, assetId=linha['asset'],
                region=grid[cid]['geom'], crs=CRS_EXPORT, crsTransform=CRS_TRANSFORM,
                maxPixels=1e13, pyramidingPolicy={'.default': 'sample'})
            if cat == 'refazer':
                log.warning(f'Apagando {linha["asset"]} (sobrescrever = sim).')
                ee.data.deleteAsset(linha['asset'])
            task.start()
            disparadas += 1
            log.info(f'  task {task.id} iniciada ({disparadas}/{args.limite}) → {nome}')
            registro.append({**linha, 'status': 'exportando', 'task_id': task.id,
                             'conta': CONTA_EXPORT})

    if not args.dry_run and disparadas >= args.limite:
        log.warning(f'Limite de {args.limite} tasks atingido — rode de novo para o restante.')

    saida = PASTA / f'registro_mosaicos_{VERSAO}_{datetime.now():%Y%m%d_%H%M}.csv'
    campos = ['ano', 'semestre', 'celula_id', 'celula_nome', 'status', 'total_cc_lt_70',
              'diurnas', 'apos_blacklist', 'zero_cloud', 'com_mascara', 'asset', 'conta', 'task_id']
    with open(saida, 'w', newline='') as f:
        escritor = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        escritor.writeheader()
        escritor.writerows(registro)
    log.info(f'Registro → {saida}')


if __name__ == '__main__':
    main()
