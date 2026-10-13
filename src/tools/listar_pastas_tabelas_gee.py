"""
Lista as pastas de um projeto do GEE com muitas FeatureCollections (TABLE)
pequenas, para analisar o que apagar. NÃO APAGA NADA.

Contexto (2026-10-12): o projeto mapbiomas-brazil (superconta) está com
199.189 de 200.000 assets (99,6%) e só ~10% do espaço usado: o limite que
estoura é a QUANTIDADE de assets, não o tamanho. Cada tabela, imagem, coleção
e pasta conta 1 — inclusive cada imagem dentro de uma ImageCollection.

O que faz
  1) mostra a cota do projeto (número de assets e espaço);
  2) percorre a árvore a partir da RAIZ escolhida (pastas em paralelo);
     coleções de imagem são só contadas (não abre imagem por imagem);
  3) por pasta: nº de tabelas, imagens, coleções (+ imagens dentro), subpastas,
     total de assets na subárvore e datas da tabela mais antiga e mais nova;
  4) opcional: consulta o tamanho de cada tabela (getAsset, 1 chamada por
     tabela) só nas N pastas com mais tabelas, e conta as "pequenas"
     (< LIMITE_PEQUENA_KB);
  5) grava dois CSV em src/tools/:
       pastas_<raiz>_<data>.csv   uma linha por pasta, ordenada por nº de tabelas
       tabelas_<raiz>_<data>.csv  uma linha por tabela (id, data, tamanho se consultado)

Leitura: no mapbiomas-brazil há pastas de várias equipes; antes de apagar,
confira de quem é a pasta. A data (updateTime) ajuda a achar resultados de
testes antigos.

Uso: python listar_pastas_tabelas_gee.py
"""
import csv
import statistics
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import ee

PASTA = Path(__file__).resolve().parent
sys.path.append(str(PASTA.parent))      # src/: configure_account_projects_ee, gee_tools
from configure_account_projects_ee import get_project_from_account  # noqa: E402
from gee_tools import switch_user                                    # noqa: E402

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
CONTAS = ['caatinga01', 'caatinga02', 'caatinga03', 'caatinga04', 'caatinga05',
          'solkan1201', 'solkanGeodatin', 'superconta', 'solkanCengine']
RAIZ_PADRAO = 'projects/mapbiomas-arida'
LIMITE_PEQUENA_KB = 100          # tabela "pequena" abaixo disso (quando o tamanho é consultado)
MOSTRAR = 30                     # pastas no resumo da tela
TRABALHADORES = 8                # chamadas em paralelo à API

ARQ_CONTA_ATUAL = Path.home() / '.config' / 'earthengine' / 'current_user.txt'


# ═══════════════════════════════════════════════════════════════════════════
# API
# ═══════════════════════════════════════════════════════════════════════════
def ativar(conta: str) -> None:
    if ARQ_CONTA_ATUAL.read_text().strip() != conta:
        switch_user(conta)
    ee.Initialize(project=get_project_from_account(conta))
    print(f'Earth Engine com {conta} ({get_project_from_account(conta)})')


def listar(pasta: str, basico: bool = False) -> list[dict]:
    """Filhos de `pasta`. `basico` = só id e tipo (rápido, para contar imagens)."""
    saida, token = [], None
    while True:
        params = {'parent': pasta, 'pageSize': 1000, **({'pageToken': token} if token else {})}
        if basico:
            params['view'] = 'BASIC'
        r = ee.data.listAssets(params)
        saida += r.get('assets', [])
        token = r.get('nextPageToken')
        if not token:
            return saida


def cota(raiz: str) -> None:
    projeto = raiz.split('/assets')[0] + '/assets/' if '/assets' in raiz else None
    if not projeto:
        return
    try:
        q = ee.data.getAssetRootQuota(projeto)
    except ee.EEException as e:
        print(f'Cota de {projeto}: indisponível ({str(e)[:80]})')
        return
    n, nl = int(q['asset_count']['usage']), int(q['asset_count']['limit'])
    s, sl = int(q['asset_size']['usage']), int(q['asset_size']['limit'])
    print(f'Cota {projeto}: {n:,} de {nl:,} assets ({100 * n / nl:.1f}%) | '
          f'{s / 2**40:.1f} de {sl / 2**40:.0f} TiB ({100 * s / sl:.1f}%)')


# ═══════════════════════════════════════════════════════════════════════════
# VARREDURA
# ═══════════════════════════════════════════════════════════════════════════
def varrer(raiz: str):
    """{pasta: estatísticas} e lista de tabelas, percorrendo a árvore em largura."""
    pastas = {}
    tabelas = []
    fila = [raiz]
    with ThreadPoolExecutor(TRABALHADORES) as ex:
        while fila:
            resultados = list(ex.map(lambda p: (p, _listar_seguro(p)), fila))
            fila = []
            colecoes = []
            for pasta, filhos in resultados:
                st = {'tabelas': 0, 'imagens': 0, 'colecoes': 0, 'imgs_em_colecoes': 0,
                      'subpastas': 0, 'datas_tab': [], 'erro': filhos if isinstance(filhos, str) else ''}
                if isinstance(filhos, str):
                    pastas[pasta] = st
                    continue
                for a in filhos:
                    t = a['type']
                    if t == 'TABLE':
                        st['tabelas'] += 1
                        st['datas_tab'].append(a.get('updateTime', '')[:10])
                        tabelas.append({'pasta': pasta, 'id': a['id'],
                                        'atualizado': a.get('updateTime', '')[:19], 'tamanho_kb': ''})
                    elif t == 'IMAGE':
                        st['imagens'] += 1
                    elif t == 'IMAGE_COLLECTION':
                        st['colecoes'] += 1
                        colecoes.append((pasta, a['id']))
                    elif t == 'FOLDER':
                        st['subpastas'] += 1
                        fila.append(a['id'])
                pastas[pasta] = st
            # imagens dentro das coleções também contam na cota
            for (pasta, col), n in zip(colecoes, ex.map(lambda c: _contar(c[1]), colecoes)):
                pastas[pasta]['imgs_em_colecoes'] += n
            print(f'  {len(pastas)} pastas lidas, {len(tabelas)} tabelas até agora…')
    return pastas, tabelas


def _listar_seguro(pasta):
    try:
        return listar(pasta)
    except ee.EEException as e:
        return f'sem acesso: {str(e)[:80]}'


def _contar(colecao):
    try:
        return len(listar(colecao, basico=True))
    except ee.EEException:
        return 0


def somar_subarvore(pastas: dict) -> dict:
    """Total de assets em cada pasta contando as subpastas (cada pasta conta 1)."""
    filhos = defaultdict(list)
    for q in pastas:
        filhos[q.rsplit('/', 1)[0]].append(q)
    total = {}
    for p in sorted(pastas, key=lambda x: x.count('/'), reverse=True):   # das folhas para a raiz
        st = pastas[p]
        proprio = st['tabelas'] + st['imagens'] + st['colecoes'] + st['imgs_em_colecoes']
        total[p] = proprio + sum(total[q] + 1 for q in filhos[p])
    return total


def consultar_tamanhos(tabelas: list[dict], pastas_alvo: set) -> None:
    alvo = [t for t in tabelas if t['pasta'] in pastas_alvo]
    print(f'Consultando o tamanho de {len(alvo)} tabelas…')

    def tamanho(t):
        try:
            return int(ee.data.getAsset(t['id']).get('sizeBytes', 0)) / 1024
        except ee.EEException:
            return None
    with ThreadPoolExecutor(TRABALHADORES) as ex:
        for i, (t, kb) in enumerate(zip(alvo, ex.map(tamanho, alvo)), 1):
            t['tamanho_kb'] = '' if kb is None else round(kb, 1)
            if i % 500 == 0:
                print(f'  {i}/{len(alvo)}')


# ═══════════════════════════════════════════════════════════════════════════
def perguntar(texto: str, padrao: str = '') -> str:
    r = input(f'{texto}' + (f' [{padrao}]' if padrao else '') + ': ').strip()
    return r or padrao


def main():
    print('══════ Pastas com muitas FeatureCollections (só leitura) ══════')
    for i, c in enumerate(CONTAS, 1):
        print(f'  {i}) {c}')
    conta = CONTAS[int(perguntar('Conta', str(CONTAS.index('superconta') + 1))) - 1]
    raiz = perguntar('Raiz a varrer', RAIZ_PADRAO).rstrip('/')
    ativar(conta)
    cota(raiz)

    inicio = datetime.now()
    pastas, tabelas = varrer(raiz)
    total = somar_subarvore(pastas)
    print(f'Varredura: {len(pastas)} pastas, {len(tabelas)} tabelas em '
          f'{(datetime.now() - inicio).seconds // 60} min')

    ordem = sorted(pastas, key=lambda p: pastas[p]['tabelas'], reverse=True)
    if perguntar(f'Consultar o tamanho das tabelas das {MOSTRAR} pastas com mais tabelas? '
                 '(1 chamada por tabela) (s/n)', 'n').lower().startswith('s'):
        consultar_tamanhos(tabelas, set(ordem[:MOSTRAR]))

    tam_por_pasta = defaultdict(list)
    for t in tabelas:
        if t['tamanho_kb'] != '':
            tam_por_pasta[t['pasta']].append(t['tamanho_kb'])

    linhas = []
    for p in ordem:
        st, kb = pastas[p], tam_por_pasta.get(p, [])
        datas = sorted(d for d in st['datas_tab'] if d)
        linhas.append({
            'pasta': p, 'tabelas': st['tabelas'], 'imagens': st['imagens'],
            'colecoes': st['colecoes'], 'imgs_em_colecoes': st['imgs_em_colecoes'],
            'subpastas': st['subpastas'], 'assets_subarvore': total[p],
            'tabela_mais_antiga': datas[0] if datas else '', 'tabela_mais_nova': datas[-1] if datas else '',
            'tabelas_consultadas': len(kb),
            'mediana_kb': round(statistics.median(kb), 1) if kb else '',
            'total_mb': round(sum(kb) / 1024, 1) if kb else '',
            f'pequenas_<{LIMITE_PEQUENA_KB}kb': sum(k < LIMITE_PEQUENA_KB for k in kb) if kb else '',
            'erro': st['erro'],
        })

    nome = raiz.replace('projects/', '').replace('/assets', '').replace('/', '_')[:60]
    carimbo = f'{datetime.now():%Y%m%d_%H%M}'
    arq_p = PASTA / f'pastas_{nome}_{carimbo}.csv'
    arq_t = PASTA / f'tabelas_{nome}_{carimbo}.csv'
    with open(arq_p, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0]))
        w.writeheader()
        w.writerows(linhas)
    with open(arq_t, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['pasta', 'id', 'atualizado', 'tamanho_kb'])
        w.writeheader()
        w.writerows(tabelas)

    print(f'\n{"tabelas":>8s} {"subárv.":>8s} {"mediana":>8s} {"pequenas":>8s}  {"período":21s}  pasta')
    for ln in linhas[:MOSTRAR]:
        if not ln['tabelas']:
            break
        print(f"{ln['tabelas']:8d} {ln['assets_subarvore']:8d} {str(ln['mediana_kb']):>8s} "
              f"{str(ln[f'pequenas_<{LIMITE_PEQUENA_KB}kb']):>8s}  "
              f"{ln['tabela_mais_antiga']}→{ln['tabela_mais_nova']:10s}  "
              f"{ln['pasta'].replace(raiz, '…')}")
    print(f'\nCSV: {arq_p}\n     {arq_t}')


if __name__ == '__main__':
    main()
