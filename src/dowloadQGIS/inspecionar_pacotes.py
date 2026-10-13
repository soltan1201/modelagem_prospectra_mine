# -*- coding: utf-8 -*-
"""
Lista o conteúdo dos pacotes do Acervo CPRM que cruzam a área de estudo, sem
baixá-los: lê só o diretório central de cada .zip por HTTP Range.

Entrada : src/Dados/acervo_cprm/acervo_cprm_na_area_estudo.csv (gerado por listar_QGIS.py)
Saída   : src/dowloadQGIS/conteudo_pacotes_area_estudo.csv  (uma linha por arquivo interno)
          src/Dados/acervo_cprm/conteudo_pacotes.json       (cache, não relê o que já leu)

Mapas em PDF ficam de fora (tipo "Mapa").

Uso: python inspecionar_pacotes.py
"""

import csv
import io
import json
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[1]
ENTRADA = RAIZ / 'Dados' / 'acervo_cprm' / 'acervo_cprm_na_area_estudo.csv'
CACHE = RAIZ / 'Dados' / 'acervo_cprm' / 'conteudo_pacotes.json'
SAIDA = Path(__file__).resolve().parent / 'conteudo_pacotes_area_estudo.csv'

TIPOS_IGNORADOS = {'Mapa'}
NAO_ZIP = ('.xlsx', '.xls', '.csv', '.pdf')   # o resto é tratado como zip
HEADERS = {'User-Agent': 'Mozilla/5.0'}


class ArquivoHTTP(io.RawIOBase):
    """Arquivo remoto só-leitura com seek, para o zipfile ler o diretório central."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        req = urllib.request.Request(url, headers={**HEADERS, 'Range': 'bytes=0-0'})
        with urllib.request.urlopen(req, timeout=60) as r:
            self.tamanho = int(r.headers['Content-Range'].split('/')[-1])

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        base = {0: 0, 1: self.pos, 2: self.tamanho}[whence]
        self.pos = base + offset
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.tamanho - self.pos
        if n == 0 or self.pos >= self.tamanho:
            return b''
        fim = min(self.pos + n, self.tamanho) - 1
        req = urllib.request.Request(
            self.url, headers={**HEADERS, 'Range': f'bytes={self.pos}-{fim}'})
        with urllib.request.urlopen(req, timeout=120) as r:
            dados = r.read()
        self.pos += len(dados)
        return dados

    def readinto(self, buffer):        # necessário para io.BufferedReader
        dados = self.read(len(buffer))
        buffer[:len(dados)] = dados
        return len(dados)


def listar_zip(url):
    """[(caminho_interno, tamanho_bytes), ...] ou {'erro': ...}."""
    try:
        with zipfile.ZipFile(ArquivoHTTP(url)) as z:
            return [(i.filename, i.file_size) for i in z.infolist() if not i.is_dir()]
    except Exception as e:   # noqa: BLE001 — servidor fora, zip corrompido etc.
        return {'erro': f'{type(e).__name__}: {e}'}


def main():
    linhas = [ln for ln in csv.DictReader(open(ENTRADA, encoding='utf-8'))
              if ln['tipo'] not in TIPOS_IGNORADOS]
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

    pendentes = [ln for ln in linhas
                 if not ln['link'].lower().endswith(NAO_ZIP) and ln['id'] not in cache]
    print(f'{len(linhas)} pacotes, {len(pendentes)} zips para ler')
    with ThreadPoolExecutor(8) as ex:
        for ln, conteudo in zip(pendentes, ex.map(lambda ln: listar_zip(ln['link']), pendentes)):
            cache[ln['id']] = conteudo
            print(f"  [{ln['id']}] {ln['titulo'][:60]}: "
                  f"{conteudo['erro'] if isinstance(conteudo, dict) else len(conteudo)}")
    CACHE.write_text(json.dumps(cache, ensure_ascii=False, indent=1))

    with open(SAIDA, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['id', 'tipo', 'titulo', 'arquivo_interno', 'extensao', 'tamanho_mb'])
        for ln in linhas:
            conteudo = cache.get(ln['id'])
            if conteudo is None:          # não é zip: o próprio arquivo
                w.writerow([ln['id'], ln['tipo'], ln['titulo'], ln['link'].rsplit('/', 1)[-1],
                            Path(ln['link']).suffix.lower(), ln['tamanho_mb']])
            elif isinstance(conteudo, dict):
                w.writerow([ln['id'], ln['tipo'], ln['titulo'], conteudo['erro'], '', ''])
            else:
                for nome, tam in conteudo:
                    w.writerow([ln['id'], ln['tipo'], ln['titulo'], nome,
                                Path(nome).suffix.lower(), round(tam / 2**20, 3)])
    print(f'Salvo: {SAIDA}')


if __name__ == '__main__':
    main()
