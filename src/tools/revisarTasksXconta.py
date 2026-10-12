#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Revisão das tasks do GEE por conta.

Fluxo (tudo por perguntas):
  1. Lista as contas numeradas e pergunta se é uma conta só ou todas.
  2. Pergunta se quer filtrar pelo nome (ex.: "ASTER_zero_cloud_"): só entram
     as tasks cuja descrição contém o texto.
  3. Pergunta se mostra todas as tasks (qualquer estado) ou só as ativas
     (RUNNING + READY), e quantas mostrar.
  4. Para cada conta: resumo por estado (RUNNING / READY / COMPLETED / FAILED…)
     e a lista das tasks selecionadas, das mais recentes para as mais antigas.
  5. Opcional: cancelar as tasks ATIVAS selecionadas (pede confirmação digitando
     CANCELAR — não tem volta).

#Produzido por Geodatin - Dados e Geoinformacao
#DISTRIBUIDO COM GPLv2
'''
import os
import sys
import collections
import datetime as dt
from pathlib import Path

import ee

pathparent = str(Path(os.getcwd()).parents[0])
sys.path.append(pathparent)
from configure_account_projects_ee import get_project_from_account
from gee_tools import switch_user

contas = {
    '1': 'caatinga01',
    '2': 'caatinga02',
    '3': 'caatinga03',
    '4': 'caatinga04',
    '5': 'caatinga05',
    '6': 'solkan1201',
    '7': 'solkanGeodatin',
    '8': 'superconta',
    '9': 'solkanCengine',
}
ATIVAS = ('RUNNING', 'READY')


# ------- PERGUNTAS -------
def perguntar(texto, validos=None, padrao=None):
    while True:
        sufixo = f' [{padrao}]' if padrao is not None else ''
        resp = input(f'{texto}{sufixo}: ').strip()
        if not resp and padrao is not None:
            return padrao
        if validos is None or resp.upper() in validos:
            return resp
        print(f'  → opções válidas: {", ".join(validos)}')


def sim(texto, padrao='N'):
    return perguntar(f'{texto} (Y/N)', {'Y', 'N'}, padrao).upper() == 'Y'


def coletar_opcoes():
    print('\nContas disponíveis:')
    for k, nome in contas.items():
        print(f'  {k}) {nome}  ({get_project_from_account(nome)})')

    unica = sim('É única conta?', 'Y')
    escolhidas = ([contas[perguntar('Qual conta (número)', set(contas))]] if unica
                  else list(contas.values()))

    filtro = ''
    if sim('Filtrar pelo nome da task? (ex.: ASTER_zero_cloud_)'):
        filtro = perguntar('Texto que o nome deve conter').strip()

    so_ativas = perguntar('Mostrar: 1) só ativas (RUNNING + READY)  2) todas', {'1', '2'}, '1') == '1'
    n = int(perguntar('Quantas tasks mostrar no máximo por conta', None, '50') or 50)
    return escolhidas, filtro, so_ativas, n


# ------- TASKS -------
def usar_conta(nome):
    switch_user(nome)
    ee.Initialize(project=get_project_from_account(nome))


def hora(ms):
    return dt.datetime.fromtimestamp(ms / 1000).strftime('%m-%d %H:%M') if ms else '—'


def revisar(conta, filtro, so_ativas, n):
    """Mostra o resumo e a lista; devolve as tasks ativas selecionadas."""
    try:
        usar_conta(conta)
    except ee.EEException as e:
        print(f'\n--- {conta}: falhou ao inicializar ({e}) ---')
        return []

    todas = ee.data.getTaskList()                       # mais recente primeiro
    selecionadas = [t for t in todas if filtro in t.get('description', '')]
    resumo = collections.Counter(t['state'] for t in selecionadas)
    ativas = [t for t in selecionadas if t['state'] in ATIVAS]

    print(f'\n═══ Conta: {conta} ({get_project_from_account(conta)}) ═══')
    if filtro:
        print(f'Filtro: nome contém "{filtro}"')
    print('Por estado: ' + (', '.join(f'{s} {c}' for s, c in resumo.most_common()) or 'nenhuma task'))
    print(f'Ativas: {len(ativas)} (RUNNING {resumo.get("RUNNING", 0)} + READY {resumo.get("READY", 0)})')

    mostrar = ativas if so_ativas else selecionadas
    for t in mostrar[:n]:
        erro = f'  ← {" ".join(t.get("error_message", "").split())[:70]}' if t['state'] == 'FAILED' else ''
        print(f'  {t["state"]:10s} criada {hora(t.get("creation_timestamp_ms"))}  '
              f'início {hora(t.get("start_timestamp_ms"))}  {t["description"]}{erro}')
    if len(mostrar) > n:
        print(f'  … e mais {len(mostrar) - n}')
    return ativas


def cancelar(conta, ativas):
    usar_conta(conta)
    for t in ativas:
        ee.data.cancelTask(t['id'])
        print(f'  cancelada: {t["description"]}')


# ------- EXECUÇÃO -------
escolhidas, filtro, so_ativas, n = coletar_opcoes()
ativas_por_conta = {c: revisar(c, filtro, so_ativas, n) for c in escolhidas}

total = sum(len(a) for a in ativas_por_conta.values())
if total and sim(f'\nCancelar as {total} task(s) ATIVAS listadas'
                 + (f' com "{filtro}"' if filtro else ' (SEM filtro de nome)') + '?'):
    print('ATENÇÃO: cancelar não tem volta.')
    if input('Digite CANCELAR para confirmar: ').strip() == 'CANCELAR':
        for conta, ativas in ativas_por_conta.items():
            if ativas:
                cancelar(conta, ativas)
    else:
        print('Confirmação não recebida — nada cancelado.')
