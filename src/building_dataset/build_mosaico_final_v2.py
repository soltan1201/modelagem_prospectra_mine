"""
VERSÃO 2 (final_v2_seco) — mesma regra de prioridade do build_mosaico_final.py,
com as entradas da estação seca (jul–out):
  1  ASTER_0cloud_v2   (build_mosaicos_zero_cloud_v2.py)
  2  ASTER_MOSAIC_v2   zero_cloud 2009–2020 (estimar_swir_rf_v2.py)
  3  ASTER_v5          semestrais secos 2000–2008 (build_mosaicos_semestrais_v2.py)
  4  ASTER_MOSAIC_v2   semestral_v5 2009–2020 (estimar_swir_rf_v2.py)
Saída: projects/mapbiomas-arida/mine/aster-mosaico-final_v2 (sem superconta).

FONTE 'sazonal' (perguntada no início): só os mosaicos sazonais secos, sem
zero-cloud nem semestrais por cima (menos colcha de retalhos; o sazonal A já
inclui as cenas zero-cloud de jun–out com peso alto):
  1  ASTER_sazonal     ASTER_sazonal_A_seco_quality_<LxCy> 2000–2008, SWIR medido
                       (build_mosaico_sazonal_seco.py)
  2  ASTER_MOSAIC_v2   ASTER_MOSAIC_sazonal_seco_2009_2020_<LxCy>, SWIR estimado
                       (estimar_swir_rf_v2.py, grupo sazonal_seco)
Saída: projects/mapbiomas-arida/mine/aster-mosaico-final_sazonal, nome
aster-mosaico-final-sazonal_<LxCy>.

FONTE 'sazonal_hm': os mesmos dois mosaicos, harmonizados contra o Landsat
(harmonizar_sazonal_landsat.py; fonte = pontos_sazonal_seco_pareados.FONTE, padrão
'hm_razao'). Como A e B ficam na mesma escala (reflectância de superfície × 1000),
dá para comparar o pixel dos dois períodos e trocar o A quando ele parece nuvem:
  origem 1  A (SWIR medido) — regra geral, o A fica por cima
  origem 2  B (SWIR estimado) onde o A não tem pixel
  origem 3  B (SWIR estimado) no lugar do A com nuvem: A mais claro que o B no B01
            (> DIF_NUVEM_B01) E mais frio no B13 (> DIF_NUVEM_B13). Nuvem é clara e
            fria; solo claro ou mina é claro mas não frio (fica o A), e sombra é
            escura (não ganha por ser escura).
Entrada B: ASTER_MOSAIC_v2/ASTER_MOSAIC_sazonal_seco<SUFIXO>_2009_2020_<LxCy>
(estimar_swir_rf_v2.py com a mesma fonte). Saída:
mine/aster-mosaico-final_sazonal<SUFIXO>/aster-mosaico-final-sazonal<SUFIXO>_<LxCy>.

Texto original (v1):
Mosaico final ASTER para inferência de materiais: junta, por célula da grid
5 × 5, os quatro mosaicos disponíveis numa só imagem, com prioridade por pixel:

  origem 1  zero-cloud 2000–2008        — sem nuvem e SWIR MEDIDO (ganha de todos)
  origem 2  zero-cloud 2009–2020        — sem nuvem, SWIR ESTIMADO (RF)
  origem 3  semestrais v4 2000–2008     — preenche gaps; SWIR MEDIDO
  origem 4  semestrais v4 2009–2020     — último recurso; SWIR ESTIMADO (RF)

Fontes
  1  ASTER_0cloud/ASTER_zero_cloud_2000_2008_<LxCy>      (build_mosaicos_zero_cloud.py)
  2  ASTER_MOSAIC/ASTER_MOSAIC_zero_cloud_2009_2020_<LxCy>    (estimar_swir_rf.py)
  3  semestrais v4 2000–2008 da célula (PASTAS_V4), compostos por qualityMosaic
     ('quality'); só semestres/pixels com as 14 bandas (2008-S2 sem SWIR fica de fora)
  4  ASTER_MOSAIC/ASTER_MOSAIC_semestral_v5_2009_2020_<LxCy>  (estimar_swir_rf.py)

Em cada camada o pixel só vale com as 14 bandas válidas e sem saturação (VNIR
DN ≥ 254 — nas cenas zero-cloud era nuvem). A camada de maior prioridade com
pixel válido fica por cima. Camada que não existe para a célula é pulada (o
registro e a propriedade 'camadas' dizem quais entraram).

Saída (ImageCollection aster-mosaico-final)
  superconta      → DESTINOS['superconta'] (mapbiomas-brazil, CAATINGA)
  qualquer outra  → DESTINOS['outras']     (projects/mapbiomas-arida/mine)
  aster-mosaico-final_<LxCy>; a união das células é o mosaico final
  (ee.ImageCollection(destino).mosaic()). Com várias contas, as tasks são
  distribuídas em rodízio (TASKS_POR_CONTA por conta). Célula já salva em
  qualquer um dos dois destinos, ou com task ativa em alguma das contas, é pulada.
  Bandas: B01–B14 (INT32; reflectância × 1000, TIR em DN), 'origem' (1–4) e
  'swir_estimado' (1 nas origens 2 e 4).

Uso: python build_mosaico_final.py
"""

import csv
import logging
from datetime import datetime

import ee

import build_mosaico_sazonal_seco as sz
import build_mosaicos_semestrais_v2 as b
import pontos_sazonal_seco_pareados as ps

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CAATINGA  = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA'
DESTINOS  = {
    'superconta': 'projects/mapbiomas-arida/mine/aster-mosaico-final_v2',   # não usado (superconta fora)
    'outras':     'projects/mapbiomas-arida/mine/aster-mosaico-final_v2',
}
PREFIXO   = 'aster-mosaico-final_'
VERSAO    = 'final_v2_seco'

# mesmas contas de tools/revisarTasksXconta.py (aquele script não é importável)
# v2: só contas que escrevem em mapbiomas-arida/mine (sem superconta)
CONTAS = ['caatinga01', 'caatinga02', 'caatinga03', 'caatinga04', 'caatinga05',
          'solkan1201', 'solkanGeodatin', 'solkanCengine']
TASKS_POR_CONTA = 2          # rodízio: troca de conta a cada N tasks
# Projeto do GEE onde as tasks rodam. Com um nome: usa só a conta ativa
# (~/.config/earthengine/current_user.txt), sem perguntar contas nem rodízio.
# None: projeto de cada conta e menu de contas. A lista de tasks ativas é por
# projeto: tasks enviadas em outro projeto não aparecem e seriam reenviadas.
PROJETO_EE = 'mapbiomas-brazil'
# Só com PROJETO_EE: células mandadas para outras contas, no projeto PRÓPRIO de cada
# uma (não no PROJETO_EE), para salvar em paralelo. Conta sem leitura nas fontes
# sai e as células dela voltam para a conta ativa. {} desliga.
CELULAS_POR_CONTA = {
    'caatinga05': ['L2C5', 'L5C1', 'L5C2'],
    'caatinga02': ['L5C3', 'L5C4', 'L5C5'],
}
# PROJETO_EE = 'mapbiomas-caatinga-cloud04'
PASTA_ZERO_CLOUD = 'projects/mapbiomas-arida/mine/ASTER_0cloud_v2'
PASTA_SWIR_RF    = 'projects/mapbiomas-arida/mine/ASTER_MOSAIC_v2'
# semestrais v4: as duas coleções (superconta e rodízio caatinga) — nomes iguais
# nas duas valem uma vez, a primeira pasta tem preferência
PASTAS_V4 = ['projects/mapbiomas-arida/mine/ASTER_v5']      # v2: semestrais secos (jul–out)

BANDAS = ['B01', 'B02', 'B3N', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09',
          'B10', 'B11', 'B12', 'B13', 'B14']
CAMADAS = {   # origem → descrição
    1: 'zero_cloud_2000_2008 (SWIR medido)',
    2: 'zero_cloud_2009_2020 (SWIR estimado RF)',
    3: 'semestral_v5_2000_2008 seco (SWIR medido)',
    4: 'semestral_v5_2009_2020 seco (SWIR estimado RF)',
}
SWIR_ESTIMADO = {2, 4}
PASTAS_LEITURA = (PASTA_ZERO_CLOUD, PASTA_SWIR_RF)   # a conta precisa ler (e existentes vêm daqui)

# Fonte 'sazonal': substitui os parâmetros acima (usar_fonte)
PASTA_SAZONAL = sz.DESTINO_FINAL
B_SWIR_RF     = 'ASTER_MOSAIC_sazonal_seco_2009_2020_'      # + <LxCy>, em PASTA_SWIR_RF
# Fonte 'sazonal_hm': troca do A pelo B onde o A parece nuvem (montar_nuvem)
DIF_NUVEM_B01 = 50     # A mais claro que o B no B01, reflectância de superfície × 1000 (0,05)
DIF_NUVEM_B13 = 60     # E mais frio no B13, DN (DN13 ≈ 1070 + 22·°C → ~3 °C)
FONTES = {
    'v2': {},
    'sazonal': {
        'DESTINOS': dict.fromkeys(DESTINOS, 'projects/mapbiomas-arida/mine/aster-mosaico-final_sazonal'),
        'PREFIXO': 'aster-mosaico-final-sazonal_',
        'VERSAO': 'final_sazonal_seco',
        'CAMADAS': {1: 'sazonal_A seco quality 2000_2008 (SWIR medido)',
                    2: 'sazonal_B seco quality 2009_2020 (SWIR estimado RF)'},
        'SWIR_ESTIMADO': {2},
        'PASTAS_LEITURA': (PASTA_SAZONAL, PASTA_SWIR_RF),
    },
    'sazonal_hm': {
        'DESTINOS': dict.fromkeys(DESTINOS, f'projects/mapbiomas-arida/mine/aster-mosaico-final_sazonal{ps.SUFIXO}'),
        'PREFIXO': f'aster-mosaico-final-sazonal{ps.SUFIXO}_',
        'VERSAO': f'final_sazonal_seco{ps.SUFIXO}',
        'CAMADAS': {1: f'sazonal_A {ps.FONTE} 2000_2008 (SWIR medido)',
                    2: f'sazonal_B {ps.FONTE} 2009_2020 (SWIR estimado RF; onde falta A ou A com nuvem)'},
        'SWIR_ESTIMADO': {2, 3},
        'PASTA_SAZONAL': ps.COL_FONTE,
        'B_SWIR_RF': f'ASTER_MOSAIC_sazonal_seco{ps.SUFIXO}_2009_2020_',
        'PASTAS_LEITURA': (ps.COL_FONTE, PASTA_SWIR_RF),
    },
}
FONTE = 'v2'


def usar_fonte(nome: str) -> None:
    """Troca destino, prefixo, versão e camadas pelos da fonte escolhida."""
    global FONTE
    FONTE = nome
    globals().update(FONTES[nome])


# ═══════════════════════════════════════════════════════════════════════════
# CAMADAS DE UMA CÉLULA
# ═══════════════════════════════════════════════════════════════════════════
def valido(img):
    """14 bandas válidas e, se houver a banda, sem saturação."""
    img = ee.Image(img)
    m = img.select(BANDAS).mask().reduce(ee.Reducer.min())
    m = ee.Image(ee.Algorithms.If(img.bandNames().contains('saturado'),
                                  m.And(img.select('saturado').unmask(0).eq(0)), m))
    return img.updateMask(m)


def ids_v4() -> list[str]:
    """IDs dos semestrais v4 das PASTAS_V4 legíveis, sem nome repetido."""
    vistos, ids = set(), []
    for pasta in PASTAS_V4:
        try:
            ativos = b.listar_destino(pasta, criar=False)
        except ee.EEException as e:
            log.warning(f'{pasta}: sem acesso ({str(e)[:80]}) — fica de fora')
            continue
        for nome, info in ativos.items():
            if nome not in vistos:
                vistos.add(nome)
                ids.append(info['id'] if 'id' in info else f'{pasta}/{nome}')
    log.info(f'Semestrais v4: {len(ids)} imagens em {len(PASTAS_V4)} pastas')
    return ids


def camadas(cel: dict, v4: ee.ImageCollection, existentes: set[str]) -> dict[int, tuple]:
    """origem → (imagem, origem_ids) das camadas que existem para a célula."""
    nome = cel['nome']
    out = {}
    a1 = f'{PASTA_ZERO_CLOUD}/ASTER_zero_cloud_2000_2008_{nome}'
    a2 = f'{PASTA_SWIR_RF}/ASTER_MOSAIC_zero_cloud_2009_2020_{nome}'
    a4 = f'{PASTA_SWIR_RF}/ASTER_MOSAIC_semestral_v5_2009_2020_{nome}'
    if a1 in existentes:
        out[1] = (ee.Image(a1), a1)
    if a2 in existentes:
        out[2] = (ee.Image(a2), a2)
    sem = (v4.filter(ee.Filter.eq('celula_nome', nome))
             .filter(ee.Filter.rangeContains('year', 2000, 2008)))
    out[3] = (sem.map(valido).select(BANDAS + ['quality']).qualityMosaic('quality'),
              sem.aggregate_array('system:index').join(','), sem.size())
    if a4 in existentes:
        out[4] = (ee.Image(a4), a4)
    return out


def camadas_sazonal(cel: dict, existentes: set[str]) -> dict[int, tuple]:
    """origem → (imagem, asset) dos mosaicos sazonais secos que existem para a célula."""
    a1 = f'{PASTA_SAZONAL}/{sz.nome_mosaico("A", cel["nome"])}'
    a2 = f'{PASTA_SWIR_RF}/{B_SWIR_RF}{cel["nome"]}'
    return {o: (ee.Image(a), a) for o, a in ((1, a1), (2, a2)) if a in existentes}


def montar(cel: dict, cam: dict[int, tuple]):
    """Empilha da menor para a maior prioridade: mosaic() deixa a última por cima."""
    imgs = []
    for origem in sorted(cam, reverse=True):
        img = valido(cam[origem][0]).select(BANDAS).round().toInt32()
        imgs.append(img.addBands([ee.Image.constant(origem).toUint8().rename('origem'),
                                  ee.Image.constant(int(origem in SWIR_ESTIMADO)).toUint8()
                                    .rename('swir_estimado')])
                       .updateMask(img.select('B01').mask()))
    final = ee.ImageCollection(imgs).mosaic().clip(cel['geom'])
    return ee.Image.cat([final.select(BANDAS).toInt32(),
                         final.select(['origem', 'swir_estimado']).toUint8()])


def montar_nuvem(cel: dict, cam: dict[int, tuple]):
    """A por cima (origem 1); B onde falta A (2) ou onde o A parece nuvem frente ao B (3)."""
    if not (1 in cam and 2 in cam):
        return montar(cel, cam)
    a, bb = (valido(cam[o][0]).select(BANDAS).round().toInt32() for o in (1, 2))
    # só onde os dois têm pixel; no resto não há troca
    nuvem_a = (a.select('B01').subtract(bb.select('B01')).gt(DIF_NUVEM_B01)
                .And(bb.select('B13').subtract(a.select('B13')).gt(DIF_NUVEM_B13))
                .unmask(0))

    def camada(img, origem):
        return (img.addBands([ee.Image.constant(origem).toUint8().rename('origem'),
                              # SWIR do B (estimado) nas origens 2 e 3; SWIR_ESTIMADO {2, 4} é do v1
                              ee.Image.constant(int(origem in (2, 3))).toUint8()
                                .rename('swir_estimado')])
                   .updateMask(img.select('B01').mask()))
    # mosaic(): a última fica por cima
    final = ee.ImageCollection([camada(bb, 2),
                                camada(bb.updateMask(nuvem_a), 3),
                                camada(a.updateMask(nuvem_a.Not()), 1)]).mosaic().clip(cel['geom'])
    return ee.Image.cat([final.select(BANDAS).toInt32(),
                         final.select(['origem', 'swir_estimado']).toUint8()])


# ═══════════════════════════════════════════════════════════════════════════
# CONTAS
# ═══════════════════════════════════════════════════════════════════════════
def destino_da_conta(conta: str) -> str:
    return DESTINOS['superconta'] if conta == 'superconta' else DESTINOS['outras']


def ativar(conta: str) -> None:
    if b.ARQ_CONTA_ATUAL.read_text().strip() != conta:
        b.switch_user(conta)
    if PROJETO_EE and conta not in CELULAS_POR_CONTA:
        ee.Initialize(project=PROJETO_EE)
    else:
        ee.Initialize(project=b.get_project_from_account(conta))


def conta_da_celula(nome: str, contas: list[str], disparadas: int) -> str:
    """Com PROJETO_EE: conta de CELULAS_POR_CONTA (se ela ficou na lista) ou a ativa.
    Sem PROJETO_EE: rodízio de TASKS_POR_CONTA."""
    if PROJETO_EE:
        return next((c for c, cels in CELULAS_POR_CONTA.items() if nome in cels and c in contas),
                    contas[0])
    return contas[(disparadas // TASKS_POR_CONTA) % len(contas)]


def perguntar_contas() -> list[str]:
    """Uma conta ou várias (rodízio), da lista CONTAS."""
    print('\nContas disponíveis:')
    for i, c in enumerate(CONTAS, 1):
        print(f'  {i}) {c:15s} ({b.get_project_from_account(c)}) → {destino_da_conta(c)}')
    while True:
        resp = input('Contas para exportar (ex.: 4 ou 4,6,8) [1,2,3,4,5]: ').strip() or '1,2,3,4,5'
        try:
            escolha = [CONTAS[int(x) - 1] for x in resp.replace(' ', '').split(',')]
        except (ValueError, IndexError):
            print(f'  → números de 1 a {len(CONTAS)}, separados por vírgula.')
            continue
        # a superconta não lê mapbiomas-arida (mosaic_aster_v4, mine/aster-mosaico-final):
        # no rodízio, as tasks dela falhariam ao ler os semestrais v4 da arida
        if 'superconta' in escolha and len(set(escolha)) > 1:
            print('  → superconta não lê mapbiomas-arida: use a superconta sozinha ou só as outras contas.')
            continue
        return list(dict.fromkeys(escolha))


def estado_geral(contas: list[str]) -> tuple[set[str], set[str]]:
    """Células já salvas (nos dois destinos) e tasks ativas nas contas escolhidas.
    Cria o destino de cada conta se faltar. Termina na primeira conta."""
    salvos, ativas = set(), set()
    destinos_lidos = set()
    for conta in list(contas):
        try:
            ativar(conta)
        except Exception as e:
            raise SystemExit(f'{conta}: não foi possível inicializar o EE ({str(e)[:120]})')
        # a task roda com as permissões da conta: sem ler as fontes, ela falharia
        try:
            for pasta in PASTAS_LEITURA:
                ee.data.getAsset(pasta)
        except ee.EEException:
            log.warning(f'{conta}: sem leitura em {pasta} — conta fora do rodízio')
            contas.remove(conta)
            continue
        ativas |= {t['description'] for t in ee.data.getTaskList()
                   if t['state'] in ('READY', 'RUNNING')
                   and t.get('description', '').startswith(PREFIXO)}
        destino = destino_da_conta(conta)
        if destino not in destinos_lidos:
            salvos |= set(b.listar_destino(destino, criar=True))
            destinos_lidos.add(destino)
    if not contas:
        raise SystemExit('Nenhuma conta escolhida lê as fontes do mosaico.')
    ativar(contas[0])
    # o outro destino também conta como "já salvo", se a conta conseguir ler
    for destino in set(DESTINOS.values()) - destinos_lidos:
        try:
            salvos |= set(b.listar_destino(destino, criar=False))
        except ee.EEException as e:
            log.warning(f'{destino}: sem acesso ({str(e)[:80]}) — não conferido')
    return salvos, ativas


# ═══════════════════════════════════════════════════════════════════════════
def main():
    usar_fonte(b.perguntar_opcao('Fonte do mosaico final:', [
        ('v2 — zero-cloud v2 + semestrais v5 (4 camadas)', 'v2'),
        ('sazonal — sazonal seco A (SWIR medido) + B (SWIR estimado)', 'sazonal'),
        (f'sazonal_hm — A e B harmonizados ({ps.FONTE}); B no lugar do A com nuvem', 'sazonal_hm'),
    ]))
    if PROJETO_EE:
        ativa = b.ARQ_CONTA_ATUAL.read_text().strip()
        contas = [ativa] + [c for c in CELULAS_POR_CONTA if c != ativa]
        print(f'Conta GEE: {ativa} (projeto {PROJETO_EE})')
        for c, cels in CELULAS_POR_CONTA.items():
            print(f'  {c} ({b.get_project_from_account(c)}): {", ".join(cels)}')
    else:
        contas = perguntar_contas()
    salvos, ativas = estado_geral(contas)
    grid = b.ler_grid()

    existentes = set()
    for pasta in PASTAS_LEITURA:
        existentes |= {f'{pasta}/{n}' for n in b.listar_destino(pasta, criar=False)}
    v4 = ee.ImageCollection([ee.Image(i) for i in ids_v4()]) if FONTE == 'v2' else None

    print(f'\n══════ Mosaico final ASTER — fonte {FONTE} '
          f'(prioridade {" > ".join(map(str, CAMADAS))}) ══════')
    for c in contas:
        print(f'  conta {c} → {destino_da_conta(c)}')
    for o, d in CAMADAS.items():
        print(f'  origem {o}: {d}')
    celulas = b.perguntar_celulas(grid)
    dry_run = b.perguntar_opcao('Modo:', [
        ('só conferir camadas por célula (dry-run, não exporta)', True),
        ('exportar o mosaico final', False),
    ])
    nomes = {c: f'{PREFIXO}{grid[c]["nome"]}' for c in celulas}
    pendentes = [c for c in celulas if dry_run or (nomes[c] not in salvos and nomes[c] not in ativas)]
    print(f'{len(celulas)} células; {len(celulas) - len(pendentes)} já salvas ou em processamento; '
          f'{len(pendentes)} a {"conferir" if dry_run else "exportar"}.')
    if not pendentes or not b.perguntar_sim_nao('Continuar?', padrao=True):
        raise SystemExit('Nada a fazer.')

    if FONTE == 'v2':
        cams = {c: camadas({'id': c, **grid[c]}, v4, existentes) for c in pendentes}
        n_sem = ee.List([cams[c][3][2] for c in pendentes]).getInfo()
    else:   # sazonal / sazonal_hm: sem semestrais (coluna semestrais_2000_2008 fica 0)
        cams = {c: camadas_sazonal(grid[c], existentes) for c in pendentes}
        n_sem = [0] * len(pendentes)

    registro = []
    disparadas, conta_ativa = 0, contas[0]
    for c, n3 in zip(pendentes, n_sem):
        cel = {'id': c, **grid[c]}
        cam = {o: v for o, v in cams[c].items() if o != 3 or n3 > 0}
        usadas = sorted(cam)
        faltam = [o for o in CAMADAS if o not in cam]
        linha = {'celula_id': c, 'celula_nome': cel['nome'], 'camadas': ','.join(map(str, usadas)),
                 'faltam': ','.join(map(str, faltam)), 'semestrais_2000_2008': n3}
        log.info(f'{nomes[c]}: camadas {usadas}' + (f' (faltam {faltam})' if faltam else '')
                 + (f'; {n3} semestrais v4 2000–2008' if FONTE == 'v2' else ''))
        if not cam:
            registro.append({**linha, 'status': 'sem_camadas'})
            continue
        if dry_run:
            registro.append({**linha, 'status': 'dry_run'})
            continue

        final = (montar_nuvem if FONTE == 'sazonal_hm' else montar)(cel, cam).set({
            'celula_id': c, 'celula_nome': cel['nome'], 'camadas': ','.join(map(str, usadas)),
            'prioridade': ' > '.join(f'{o}={CAMADAS[o]}' for o in CAMADAS),
            **{f'origem_{o}': cam[o][1] if o in cam else '' for o in CAMADAS},
            'grid_asset': b.GRID_ASSET, 'DATA_TYPE': 'INT32',
            'SCALE_FACTORS': {**dict.fromkeys(BANDAS[:9], b.ESCALA_TOA), **dict.fromkeys(BANDAS[9:], 1),
                              'origem': 1, 'swir_estimado': 1},
            'versao': VERSAO,
            **({'regra_nuvem': f'B01_A - B01_B > {DIF_NUVEM_B01} e B13_B - B13_A > {DIF_NUVEM_B13}',
                'origem_3': 'B no lugar do A com nuvem'} if FONTE == 'sazonal_hm' else {}),
        })
        conta = conta_da_celula(cel['nome'], contas, disparadas)
        if conta != conta_ativa:
            ativar(conta)
            conta_ativa = conta
        destino = destino_da_conta(conta)
        task = ee.batch.Export.image.toAsset(
            image=ee.Image(final), description=nomes[c], assetId=f'{destino}/{nomes[c]}',
            region=cel['geom'], crs=b.CRS_EXPORT, crsTransform=b.CRS_TRANSFORM, maxPixels=1e13,
            pyramidingPolicy={'.default': 'mean', 'origem': 'mode', 'swir_estimado': 'mode'})
        try:
            task.start()
        except ee.EEException as e:
            log.error(f'{nomes[c]}: falhou ao iniciar: {e}')
            registro.append({**linha, 'conta': conta, 'status': f'erro_inicio: {str(e)[:80]}'})
            continue
        disparadas += 1
        log.info(f'  [{conta}] task {task.id} → {destino}/{nomes[c]}')
        registro.append({**linha, 'conta': conta, 'asset': f'{destino}/{nomes[c]}',
                         'status': 'exportando', 'task_id': task.id})

    if conta_ativa != contas[0]:      # deixa a conta ativa de antes (outros scripts usam ela)
        ativar(contas[0])
    arq = b.PASTA / f'registro_mosaico_final_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['celula_id', 'celula_nome', 'camadas', 'faltam',
                                          'semestrais_2000_2008', 'conta', 'asset', 'status', 'task_id'],
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(registro)
    log.info(f'Registro → {arq}')


if __name__ == '__main__':
    main()
