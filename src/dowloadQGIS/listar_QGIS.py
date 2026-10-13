# -*- coding: utf-8 -*-
"""
Lista as camadas do plugin "Acervo SIG - CPRM" (acervo_cprm) que caem na área
de estudo (mesmo retângulo de select_black_list_save_mosaicSem.js).

O catálogo do plugin não tem geometria: a posição de cada camada é deduzida
pelo próprio plugin a partir do código de folha da CIM (ex.: SC.24-V-A) ou do
número do projeto aerogeofísico. Este script usa essas mesmas funções
(acervo_cprm.folhas) e separa as camadas em:
  - dentro : a caixa da folha/projeto cruza a área de estudo
  - fora   : a caixa não cruza
  - sem_caixa : o catálogo não diz onde a camada fica (LiDAR, ortofoto, parte
                da geoquímica) — conferir manualmente

Saídas:
  - CSV com todas as camadas do catálogo + situação + caixa
  - CSV só com as que cruzam a área, ordenadas por tipo e escala
  - No QGIS: camada de memória com as caixas das folhas que cruzam a área e
    o retângulo da área de estudo

Uso:
  - QGIS > Complementos > Terminal Python > Editor > abrir este arquivo > Executar
  - Fora do QGIS também roda (só gera os CSV): python listar_QGIS.py
"""

import csv
import sys
from collections import Counter
from pathlib import Path

# ------------------------------------------------------------------
# Configuração
# ------------------------------------------------------------------
# Área de estudo: (oeste, sul, leste, norte) em graus
AREA_ESTUDO = (-44.460483519818084, -14.790793773956535,
               -39.203525512005584, -9.416525285369607)

PASTA_PLUGINS = (Path.home() / '.local/share/QGIS/QGIS4/profiles/default'
                 / 'python/plugins')

# No editor do QGIS __file__ nem sempre existe; aí vale o caminho fixo
try:
    RAIZ = Path(__file__).resolve().parents[1]
except NameError:
    RAIZ = Path('/home/superuser/Dados/proj/prosprectra/modelagem_prospectra_mine/src')
PASTA_SAIDA = RAIZ / 'Dados' / 'acervo_cprm'

ADICIONAR_AO_PROJETO = True   # cria as camadas de memória no QGIS


# ------------------------------------------------------------------
# Carga do catálogo do plugin
# ------------------------------------------------------------------
if str(PASTA_PLUGINS) not in sys.path:
    sys.path.append(str(PASTA_PLUGINS))

from acervo_cprm import catalogo, folhas   # noqa: E402


def caminho_catalogo():
    """Cópia atualizada pelo botão 'Atualizar catálogo', se houver; senão o embutido."""
    try:
        from acervo_cprm import config
        atualizado = config.pasta_catalogo_atualizado()
        if atualizado.exists():
            return atualizado
    except Exception:
        pass
    return catalogo.EMBUTIDO


def escala_do_codigo(codigo):
    """Escala da folha pelo número de sufixos do código CIM."""
    if not codigo:
        return ''
    n = codigo.count('-')
    return {0: '1:1.000.000', 1: '1:500.000', 2: '1:250.000',
            3: '1:100.000', 4: '1:50.000', 5: '1:25.000'}.get(n, '')


def area_caixa(c):
    return max(c[2] - c[0], 0) * max(c[3] - c[1], 0)


def intersecao(a, b):
    return (max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3]))


def classificar(camadas):
    linhas = []
    for cam in camadas:
        cx = folhas.caixa_da_camada(cam)
        codigo = None
        for texto in (cam.titulo,) + tuple(cam.pastas):
            codigo = folhas.codigo_do_texto(texto)
            if codigo:
                break
        if cx is None:
            situacao = 'sem_caixa'
        elif folhas.cruza(cx, AREA_ESTUDO):
            situacao = 'dentro'
        else:
            situacao = 'fora'

        # fração da folha que cai na área e fração da área coberta pela folha
        frac_folha = frac_area = None
        if situacao == 'dentro':
            inter = area_caixa(intersecao(cx, AREA_ESTUDO))
            frac_folha = round(inter / area_caixa(cx), 3) if area_caixa(cx) else None
            frac_area = round(inter / area_caixa(AREA_ESTUDO), 3)

        linhas.append({
            'situacao': situacao,
            'id': cam.id,
            'tipo': cam.tipo,
            'titulo': cam.titulo,
            'pasta': ' / '.join(cam.pastas),
            'codigo_folha': codigo or '',
            'escala': escala_do_codigo(codigo),
            'formato': cam.formato,
            'tamanho_mb': round(cam.tamanho_bytes / 2**20, 2),
            'frac_folha_na_area': frac_folha,
            'frac_area_coberta': frac_area,
            'oeste': cx[0] if cx else None,
            'sul': cx[1] if cx else None,
            'leste': cx[2] if cx else None,
            'norte': cx[3] if cx else None,
            'link': cam.link,
            '_caixa': cx,
        })
    return linhas


def salvar_csv(linhas, arquivo):
    campos = [k for k in linhas[0] if not k.startswith('_')]
    with open(arquivo, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        w.writeheader()
        w.writerows(linhas)


def adicionar_camadas_qgis(dentro):
    """Caixas das camadas que cruzam a área + retângulo da área, em memória."""
    from qgis.core import (QgsFeature, QgsGeometry, QgsProject,
                           QgsRectangle, QgsVectorLayer)

    uri = ('Polygon?crs=EPSG:4326'
           '&field=id:string&field=tipo:string&field=titulo:string(254)'
           '&field=codigo:string&field=escala:string&field=tamanho_mb:double'
           '&field=link:string(254)')
    lyr = QgsVectorLayer(uri, 'Acervo CPRM - camadas na área', 'memory')
    feats = []
    for ln in dentro:
        f = QgsFeature(lyr.fields())
        f.setGeometry(QgsGeometry.fromRect(QgsRectangle(*ln['_caixa'])))
        f.setAttributes([ln['id'], ln['tipo'], ln['titulo'], ln['codigo_folha'],
                         ln['escala'], ln['tamanho_mb'], ln['link']])
        feats.append(f)
    lyr.dataProvider().addFeatures(feats)
    lyr.updateExtents()

    area = QgsVectorLayer('Polygon?crs=EPSG:4326&field=nome:string',
                          'Área de estudo', 'memory')
    fa = QgsFeature(area.fields())
    fa.setGeometry(QgsGeometry.fromRect(QgsRectangle(*AREA_ESTUDO)))
    fa.setAttributes(['area_estudo'])
    area.dataProvider().addFeatures([fa])
    area.updateExtents()

    QgsProject.instance().addMapLayers([lyr, area])


# ------------------------------------------------------------------
# Execução
# ------------------------------------------------------------------
arq_catalogo = caminho_catalogo()
camadas = [c for c in catalogo.carregar(arq_catalogo)]
linhas = classificar(camadas)

dentro = sorted((ln for ln in linhas if ln['situacao'] == 'dentro'),
                key=lambda ln: (ln['tipo'], ln['codigo_folha'].count('-'),
                                ln['titulo']))
cont = Counter(ln['situacao'] for ln in linhas)

PASTA_SAIDA.mkdir(parents=True, exist_ok=True)
salvar_csv(linhas, PASTA_SAIDA / 'acervo_cprm_catalogo_completo.csv')
if dentro:
    salvar_csv(dentro, PASTA_SAIDA / 'acervo_cprm_na_area_estudo.csv')

print('=' * 70)
print(f'Catálogo : {arq_catalogo}')
print(f'Camadas  : {len(linhas)}  | dentro: {cont["dentro"]}  '
      f'| fora: {cont["fora"]}  | sem caixa: {cont["sem_caixa"]}')
print('=' * 70)

print('\nCamadas que cruzam a área de estudo, por tipo:')
por_tipo = Counter(ln['tipo'] for ln in dentro)
tam_tipo = Counter()
for ln in dentro:
    tam_tipo[ln['tipo']] += ln['tamanho_mb']
for tipo, n in por_tipo.most_common():
    print(f'  {tipo:28s} {n:5d} camadas  {tam_tipo[tipo] / 1024:8.2f} GB')

tipo_atual = None
for ln in dentro:
    if ln['tipo'] != tipo_atual:
        tipo_atual = ln['tipo']
        print(f'\n--- {tipo_atual} ---')
    print(f"  [{ln['id']:>5}] {ln['escala']:12s} {ln['tamanho_mb']:9.1f} MB  "
          f"cobre {ln['frac_area_coberta'] * 100:5.1f}% da área  {ln['titulo']}")

print('\nSem localização no catálogo (conferir manualmente), por tipo:')
for tipo, n in Counter(ln['tipo'] for ln in linhas
                       if ln['situacao'] == 'sem_caixa').most_common():
    print(f'  {tipo:28s} {n:5d}')

print(f'\nCSV salvos em: {PASTA_SAIDA}')

def dentro_do_qgis():
    try:
        from qgis.utils import iface
    except ImportError:
        return False
    return iface is not None


if ADICIONAR_AO_PROJETO and dentro and dentro_do_qgis():
    adicionar_camadas_qgis(dentro)
    print('Camadas de memória adicionadas ao projeto.')
