"""
Transfere as coleções ASTER de mapbiomas-brazil (CAATINGA) para mapbiomas-arida.

Por que copiar e não mover (ee.data.renameAsset):
  - renameAsset MOVE: apaga a origem. Exige escrita no mapbiomas-brazil, e a
    caatinga04 só tem LEITURA lá (as coleções estão com all_users_can_read).
  - Mover não tem volta e, se falhar no meio, as coleções ficam divididas.
  - copyAsset só precisa de leitura na origem e escrita no destino. A cópia é
    feita no servidor (sem task, sem baixar nada) e leva as propriedades da imagem.
  A remoção dos originais (para liberar espaço da superconta) é uma etapa à parte,
  só para imagens que já estão no destino, com a superconta e confirmação.

Etapas (perguntadas ao rodar):
  1) conferir  : lista origem × destino (o que falta copiar), não altera nada
  2) copiar    : com a caatinga04; cria a ImageCollection de destino se faltar,
                 copia as propriedades da coleção e cada imagem que ainda não
                 está no destino
  3) apagar originais : confere com a caatinga04 que a cópia existe no destino
                 com as mesmas bandas e, só para essas, apaga a origem com a
                 superconta (que não lê mapbiomas-arida). Pede confirmação
                 digitando APAGAR. NÃO TEM VOLTA.

Depois de copiar, os scripts que leem as pastas antigas precisam apontar para
as novas (ex.: build_mosaico_final.py, estimar_swir_rf.py, build_mosaicos_semestrais.py).

Uso: python transferir_assets_brazil_para_arida.py
"""
import csv
import sys
import time
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
ORIGEM  = 'projects/mapbiomas-brazil/assets/LAND-COVER/COLLECTION-11/GENERAL/SAMPLES/CAATINGA'
DESTINO = 'projects/mapbiomas-arida/mine'

# coleção na origem → coleção no destino (comente as que não quiser transferir)
COLECOES = {
    'ASTER':        'ASTER',
    'ASTER_0cloud': 'ASTER_0cloud',
    'ASTER_MOSAIC': 'ASTER_MOSAIC',
    'ASTER_v4':     'ASTER_v4',
    # 'ASTER_v3':            'ASTER_v3',              # vazia em 2026-10-12
    # 'aster-mosaico-final': 'aster-mosaico-final',   # vazia na origem; o destino já existe
}

CONTA_COPIA  = 'caatinga04'     # lê mapbiomas-brazil e escreve em mapbiomas-arida
CONTA_APAGAR = 'superconta'     # dona do mapbiomas-brazil
TENTATIVAS   = 4                # por imagem, em erro transitório

ARQ_CONTA_ATUAL = Path.home() / '.config' / 'earthengine' / 'current_user.txt'


# ═══════════════════════════════════════════════════════════════════════════
# UTILITÁRIOS
# ═══════════════════════════════════════════════════════════════════════════
def ativar(conta: str) -> None:
    if ARQ_CONTA_ATUAL.read_text().strip() != conta:
        switch_user(conta)
    ee.Initialize(project=get_project_from_account(conta))
    print(f'Earth Engine com {conta} ({get_project_from_account(conta)})')


def listar(pasta: str) -> dict[str, dict]:
    """{nome: asset} dos filhos de `pasta`; {} se não existir ou sem acesso."""
    try:
        ee.data.getAsset(pasta)
    except ee.EEException:
        return {}
    saida, token = {}, None
    while True:
        r = ee.data.listAssets({'parent': pasta, **({'pageToken': token} if token else {})})
        for a in r.get('assets', []):
            saida[a['id'].rsplit('/', 1)[-1]] = a
        token = r.get('nextPageToken')
        if not token:
            return saida


def com_tentativas(func, *args):
    for i in range(TENTATIVAS):
        try:
            return func(*args)
        except ee.EEException as e:
            transitorio = any(t in str(e).lower() for t in ('timeout', 'deadline', 'unavailable',
                                                            'internal', 'too many', 'quota'))
            if not transitorio or i == TENTATIVAS - 1:
                raise
            time.sleep(5 * 2 ** i)


def perguntar(texto: str, padrao: str = '') -> str:
    r = input(f'{texto}' + (f' [{padrao}]' if padrao else '') + ': ').strip()
    return r or padrao


def situacao() -> dict[str, tuple[dict, dict]]:
    """origem/destino de cada coleção: {col_origem: (imagens_origem, imagens_destino)}."""
    estado = {}
    print(f'\n{"coleção":15s} {"origem":>7s} {"destino":>8s} {"faltam":>7s}   destino')
    for col_o, col_d in COLECOES.items():
        orig = listar(f'{ORIGEM}/{col_o}')
        dest = listar(f'{DESTINO}/{col_d}')
        faltam = len(set(orig) - set(dest))
        print(f'{col_o:15s} {len(orig):7d} {len(dest):8d} {faltam:7d}   {DESTINO}/{col_d}')
        estado[col_o] = (orig, dest)
    return estado


def gravar_registro(linhas: list[dict], etapa: str) -> None:
    if not linhas:
        return
    arq = PASTA / f'registro_transferencia_{etapa}_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['colecao', 'imagem', 'origem', 'destino', 'status'])
        w.writeheader()
        w.writerows(linhas)
    print(f'Registro → {arq}')


# ═══════════════════════════════════════════════════════════════════════════
# ETAPAS
# ═══════════════════════════════════════════════════════════════════════════
def copiar(estado: dict) -> None:
    registro, total, erros = [], 0, 0
    for col_o, (orig, dest) in estado.items():
        col_d = COLECOES[col_o]
        destino_col = f'{DESTINO}/{col_d}'
        faltam = sorted(set(orig) - set(dest))
        if not faltam:
            print(f'{col_o}: nada a copiar')
            continue

        if not dest:
            try:
                ee.data.getAsset(destino_col)
            except ee.EEException:
                ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, destino_col)
                print(f'{destino_col}: ImageCollection criada')
        # propriedades da coleção (descrição, metadados), se houver
        props = ee.data.getAsset(f'{ORIGEM}/{col_o}').get('properties') or {}
        if props:
            ee.data.updateAsset(destino_col, {'properties': props}, ['properties'])

        print(f'{col_o}: copiando {len(faltam)} imagens → {destino_col}')
        for i, nome in enumerate(faltam, 1):
            src, dst = f'{ORIGEM}/{col_o}/{nome}', f'{destino_col}/{nome}'
            try:
                com_tentativas(ee.data.copyAsset, src, dst)
                status = 'copiado'
                total += 1
            except ee.EEException as e:
                status = f'erro: {str(e)[:120]}'
                erros += 1
                print(f'  ERRO {nome}: {str(e)[:120]}')
            registro.append({'colecao': col_o, 'imagem': nome, 'origem': src,
                             'destino': dst, 'status': status})
            if i % 25 == 0 or i == len(faltam):
                print(f'  {i}/{len(faltam)}')
    gravar_registro(registro, 'copia')
    print(f'Copiadas: {total} | erros: {erros}')


def apagar_originais(estado: dict) -> None:
    """Apaga da origem só as imagens cuja cópia existe no destino com as mesmas bandas.

    A conferência é feita com a CONTA_COPIA (a superconta não lê mapbiomas-arida);
    a superconta entra só para apagar."""
    conferidas, registro = [], []
    for col_o, (orig, dest) in estado.items():
        for nome in sorted(set(orig) & set(dest)):
            src = f'{ORIGEM}/{col_o}/{nome}'
            dst = f'{DESTINO}/{COLECOES[col_o]}/{nome}'
            linha = {'colecao': col_o, 'imagem': nome, 'origem': src, 'destino': dst}
            try:
                bandas_o = [b['id'] for b in ee.data.getAsset(src).get('bands', [])]
                bandas_d = [b['id'] for b in ee.data.getAsset(dst).get('bands', [])]
            except ee.EEException as e:
                registro.append({**linha, 'status': f'mantido (não conferido: {str(e)[:80]})'})
                continue
            if bandas_o and bandas_o == bandas_d:
                conferidas.append(linha)
            else:
                registro.append({**linha, 'status': 'mantido (bandas diferentes)'})

    sem_copia = sum(len(set(o) - set(d)) for o, d in estado.values())
    print(f'\n{len(conferidas)} imagens conferidas (cópia com as mesmas bandas) podem ser apagadas da origem.')
    if sem_copia:
        print(f'{sem_copia} imagens AINDA NÃO têm cópia e ficam na origem.')
    if len(registro):
        print(f'{len(registro)} imagens com cópia divergente ou não conferida ficam na origem.')
    if not conferidas:
        gravar_registro(registro, 'apagar')
        return

    print('ATENÇÃO: apagar da origem NÃO TEM VOLTA. Scripts que leem as pastas antigas param de funcionar.')
    if perguntar('Digite APAGAR para confirmar') != 'APAGAR':
        print('Confirmação não recebida — nada apagado.')
        return

    ativar(CONTA_APAGAR)
    apagadas = 0
    for linha in conferidas:
        try:
            com_tentativas(ee.data.deleteAsset, linha['origem'])
            apagadas += 1
            registro.append({**linha, 'status': 'apagado da origem'})
        except ee.EEException as e:
            registro.append({**linha, 'status': f'erro ao apagar: {str(e)[:100]}'})
    gravar_registro(registro, 'apagar')
    print(f'Apagadas da origem: {apagadas} de {len(conferidas)}')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    print('══════ Transferência mapbiomas-brazil → mapbiomas-arida ══════')
    print(f'  origem : {ORIGEM}\n  destino: {DESTINO}')
    print('Etapas:\n  1) conferir (não altera nada)\n  2) copiar (conta '
          f'{CONTA_COPIA})\n  3) apagar originais já copiados (conta {CONTA_APAGAR})')
    etapa = perguntar('Etapa', '1')

    ativar(CONTA_COPIA)
    estado = situacao()

    if etapa == '2':
        if perguntar('Copiar as imagens que faltam? (s/n)', 's').lower().startswith('s'):
            copiar(estado)
            print('\nSituação depois da cópia:')
            situacao()
    elif etapa == '3':
        apagar_originais(estado)


if __name__ == '__main__':
    main()
