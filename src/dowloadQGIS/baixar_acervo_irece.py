"""
Baixa do Acervo CPRM (GEOSGB) as camadas selecionadas para a área de estudo
de Irecê e grava em /mnt/Almacen/shpIRECE/<série>/<pacote>/.

Não baixa o .zip inteiro: lê o diretório central por HTTP Range e extrai só
os arquivos escolhidos (temas geológicos + arquivos auxiliares do shapefile).

Seleção (ver camadas_acervo_cprm_area_estudo.txt):
  milionesimo/   4 cartas geológicas ao milionésimo que cruzam a área
  250k/          Jacobina SC.24-Y-C (contém Irecê)
  100k/          todas as cartas 1:100.000 que cruzam a área (só a versão
                 mais nova quando a folha tem duas) — filtrar depois
  50k/           todas as cartas 1:50.000 que cruzam a área
  aerogeofisica/ 1042 Norte da Chapada Diamantina — grids Geotif
  geoquimica/    Morro do Chapéu PLGB, Agrominerais Irecê-Jaguarari,
                 Calcários da Bahia, Geoquímica do Bambuí I e II (BA)
  regional/      mapas sem código de folha que cobrem a área

Temas guardados dos vetores: lito, recmin, estrutura, afloramento e o limite
da folha (para filtrar pela área depois).

Entrada : src/Dados/acervo_cprm/acervo_cprm_catalogo_completo.csv (listar_QGIS.py)
Saída   : /mnt/Almacen/shpIRECE/ + manifesto.csv

Uso: python baixar_acervo_irece.py [--dry-run]
"""
import csv
import io
import re
import sys
import urllib.request
import zipfile
from pathlib import Path

from inspecionar_pacotes import HEADERS, ArquivoHTTP

RAIZ = Path(__file__).resolve().parents[1]
CATALOGO = RAIZ / 'Dados' / 'acervo_cprm' / 'acervo_cprm_catalogo_completo.csv'
DESTINO = Path('/mnt/Almacen/shpIRECE')

# ids escolhidos à mão (sem código de folha ou fora das séries por regra)
IDS_FIXOS = {
    'aerogeofisica': ['942'],                               # 1042-Geotif
    'geoquimica': ['4599', '7420', '6719', '4643', '4646'],
    'regional': ['5984', '6691', '6021', '6714', '3594', '3395', '3398', '3528'],
}

# Temas geológicos: casam com o NOME do shapefile (sem extensão)
RE_TEMA = re.compile(
    r'lito|recmin|rec_min|ocormin|estrut|aflora|limite_d?a?_?folha|limite_projeto'
    r'|limite_do_projeto|dolinas|dominios_tectonicos', re.I)
EXT_SHP = {'.shp', '.shx', '.dbf', '.prj', '.cpg', '.sbn', '.sbx', '.qix', '.xml'}
EXT_RASTER = {'.tif', '.ipj', '.gi'}          # aerogeofísica: grids + projeção
EXT_TABELA = {'.xlsx', '.xls', '.csv'}


def serie_por_regra(ln):
    if ln['situacao'] != 'dentro' or ln['tipo'] != 'SIG (Vetores)':
        return None
    pasta = ln['pasta']
    if 'Cartas Geológicas do Brasil ao Milionésimo' in pasta:
        return 'milionesimo'
    if 'Cartas Geológicas, Escala 1:250.000' in pasta and 'SC.24-Y-C' in ln['codigo_folha']:
        return '250k'
    if 'Cartas Geológicas, Escala 1:100.000' in pasta:
        return '100k'
    if 'Cartas Geológicas, Escala 1:50.000' in pasta:
        return '50k'
    return None


def ano(titulo):
    m = re.search(r'(19|20)\d{2}', titulo)
    return int(m.group(0)) if m else 0


def selecionar():
    cat = list(csv.DictReader(open(CATALOGO, encoding='utf-8')))
    por_id = {ln['id']: ln for ln in cat}
    sel = []
    for ln in cat:
        serie = serie_por_regra(ln)
        if serie:
            sel.append((serie, ln))

    # folha com duas versões (2014/2019): fica a mais nova
    melhor = {}
    for serie, ln in sel:
        chave = (serie, ln['codigo_folha'])
        atual = melhor.get(chave)
        if atual is None or (ano(ln['titulo']), int(ln['id'])) > (ano(atual['titulo']), int(atual['id'])):
            melhor[chave] = ln
    sel = [(s, ln) for s, ln in sel if melhor[(s, ln['codigo_folha'])] is ln]

    for serie, ids in IDS_FIXOS.items():
        sel += [(serie, por_id[i]) for i in ids]
    return sel


def nome_pasta(ln):
    nome = re.sub(r'^SIG \(Vetores\)\s*-\s*', '', ln['titulo'])
    nome = re.sub(r'[^\w\-.]+', '_', nome, flags=re.UNICODE).strip('_')
    return f"{ln['id']}_{nome[:70]}"


def escolher_membros(nomes, serie):
    """Arquivos do zip a extrair, conforme a série."""
    if serie == 'aerogeofisica':
        return [n for n in nomes if Path(n).suffix.lower() in EXT_RASTER]
    if serie == 'geoquimica':
        return [n for n in nomes if Path(n).suffix.lower() in EXT_TABELA | EXT_SHP]

    def base(n):          # "pasta/x.shp.xml" e "pasta/x.dbf" -> "pasta/x"
        n = n.lower().removesuffix('.xml')
        return n.rsplit('.', 1)[0]

    temas = {base(n) for n in nomes
             if n.lower().endswith('.shp') and RE_TEMA.search(Path(n).stem)}
    return [n for n in nomes
            if Path(n).suffix.lower() in EXT_SHP and base(n) in temas]


def baixar_arquivo(url, destino):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=600) as r:
        destino.write_bytes(r.read())


def processar(serie, ln, dry_run):
    pasta = DESTINO / serie / nome_pasta(ln)
    url = ln['link']
    if url.lower().endswith(tuple(EXT_TABELA)):
        alvo = pasta / url.rsplit('/', 1)[-1]
        if not dry_run and not alvo.exists():
            pasta.mkdir(parents=True, exist_ok=True)
            baixar_arquivo(url, alvo)
        return [alvo.name]

    remoto = io.BufferedReader(ArquivoHTTP(url), buffer_size=4 * 2**20)
    with zipfile.ZipFile(remoto) as z:
        nomes = [i.filename for i in z.infolist() if not i.is_dir()]
        membros = escolher_membros(nomes, serie)
        if dry_run:
            return membros
        for m in membros:
            alvo = pasta / Path(m).name          # achata as subpastas do zip
            if alvo.exists():
                continue
            alvo.parent.mkdir(parents=True, exist_ok=True)
            alvo.write_bytes(z.read(m))
    return membros


def main():
    dry_run = '--dry-run' in sys.argv
    sel = selecionar()
    print(f'{len(sel)} pacotes selecionados -> {DESTINO}' + ('  (dry-run)' if dry_run else ''))
    manifesto = []
    for k, (serie, ln) in enumerate(sel, 1):
        try:
            membros = processar(serie, ln, dry_run)
            n_shp = sum(m.lower().endswith('.shp') for m in membros)
            print(f"[{k:3d}/{len(sel)}] {serie:13s} [{ln['id']}] {ln['titulo'][:55]:55s} "
                  f"{len(membros):3d} arq ({n_shp} shp)")
            erro = ''
        except Exception as e:   # noqa: BLE001 — registra e segue
            membros, erro = [], f'{type(e).__name__}: {e}'
            print(f"[{k:3d}/{len(sel)}] {serie:13s} [{ln['id']}] ERRO {erro}")
        manifesto.append({'serie': serie, 'id': ln['id'], 'titulo': ln['titulo'],
                          'codigo_folha': ln['codigo_folha'], 'pasta': nome_pasta(ln),
                          'n_arquivos': len(membros),
                          'shapefiles': ';'.join(Path(m).name for m in membros
                                                 if m.lower().endswith(('.shp', '.tif', '.xlsx', '.xls'))),
                          'link': ln['link'], 'erro': erro})
    if not dry_run:
        DESTINO.mkdir(parents=True, exist_ok=True)
        with open(DESTINO / 'manifesto.csv', 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=list(manifesto[0]))
            w.writeheader()
            w.writerows(manifesto)
        print(f'Manifesto: {DESTINO / "manifesto.csv"}')


if __name__ == '__main__':
    main()
