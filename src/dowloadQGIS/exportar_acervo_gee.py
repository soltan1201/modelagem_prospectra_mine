"""
Exporta para assets do GEE os dados do Acervo CPRM selecionados em
ACERVO_CPRM_IRECE.md (seção 10).

Regras (decididas em 2026-10-12):
  - resolução NATIVA (gama do 1042 em ~490 m, MAG em ~220 m), sem reamostrar;
  - extensão do dado BAIXADO, sem recortar pela área de estudo;
  - vetores reprojetados para EPSG:4326 (os originais misturam 4674/4326/31984).

Etapas (perguntadas ao rodar):
  1) preparar : monta localmente, em PASTA_PREP, os arquivos a subir
                (GeoTIFF multibanda do gama, GeoTIFF da MAG, shapefiles .zip)
  2) enviar   : copia os arquivos para o bucket do Cloud Storage e dispara a
                ingestão de cada um no GEE (ee.data.startIngestion /
                startTableIngestion). A ingestão é uma task: acompanhar em
                tools/revisarTasksXconta.py ou no Code Editor.

O GEE só ingere arquivos que estão no Cloud Storage (gs://). Pré-requisitos do envio:
  - um bucket em que a conta escolhida possa escrever e o GEE possa ler
    (BUCKET abaixo, ou digitado ao rodar), num projeto com faturamento
    (PROJETO_STORAGE; ver GUIA_GCLOUD_BUCKET_GEE.md);
  - credenciais com acesso ao Storage: as da própria conta do EE, se ela foi
    autenticada com o escopo devstorage; senão as do gcloud
    (`gcloud auth application-default login`).

Assets criados em DESTINO (pasta criada se faltar):
  1042_gama            ee.Image, 7 bandas: eU_ppm, eTh_ppm, K_pct, CT, U_Th, U_K, Th_K
  1042_mag             ee.Image, 1 banda: MAG_nT
  ocorrencias_fosfato  pontos de fósforo de todos os recmin (sem duplicatas)
  recursos_minerais    todos os recmin juntos (sem duplicatas)
  geoquimica_p2o5      amostras de rocha de 7420 + 6719 com P2O5 numérico
  lito_jacobina_250k, estrut_jacobina_250k   folha SC.24-Y-C (contém Irecê)
  lito_bahia, estrut_bahia                   Geologia do Estado da Bahia (5984)
  dolinas_morro_chapeu                        Projeto Morro do Chapéu (6691)
  (+ lito_100k, estrut_100k se INCLUIR_100K = True)

Uso: python exportar_acervo_gee.py
"""
import csv
import math
import re
import shutil
import sys
import unicodedata
import zipfile
from datetime import datetime
from pathlib import Path

import openpyxl
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()

PASTA = Path(__file__).resolve().parent
sys.path.append(str(PASTA.parent))      # src/: configure_account_projects_ee, gee_tools

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
DADOS      = Path('/mnt/Almacen/shpIRECE')
PASTA_PREP = DADOS / 'gee_upload'
DESTINO    = 'projects/mapbiomas-arida/mine/acervo_cprm'
BUCKET     = 'gs://cengine17-acervo-cprm/acervo_cprm'   # vazio = perguntar ao rodar
# projeto do Google Cloud dono do bucket (o que tem faturamento). Vazio = o projeto da
# conta do EE. Ex.: o projeto do teste gratuito, deixando o projeto do EE sem faturamento
PROJETO_STORAGE = 'project-fadb9cba-032b-4fb9-882'   # teste gratuito de solkan.cengine17
INCLUIR_100K = False               # sobe também lito/estrutura das 56 folhas 1:100k

IRECE = (-41.86, -11.30)           # referência para dist_irece_km e grupo
AREA  = (-44.460483519818084, -14.790793773956535,
         -39.203525512005584, -9.416525285369607)   # só para o campo na_area
NODATA = -9999.0

# mesmas contas de tools/revisarTasksXconta.py; superconta não escreve em mapbiomas-arida
CONTAS = ['caatinga01', 'caatinga02', 'caatinga03', 'caatinga04', 'caatinga05',
          'solkan1201', 'solkanGeodatin', 'solkanCengine']

GAMA = ['eU_ppm', 'eTh_ppm', 'K_pct', 'CT', 'U_Th', 'U_K', 'Th_K']
GRIDS = DADOS / 'aerogeofisica' / '1042_grids'
PROPS_1042 = {
    'projeto': '1042 Norte da Chapada Diamantina',
    'fonte': 'SGB/CPRM Acervo GEOSGB, pacote 943 (1042-XYZ); grids gama do projeto BARMP',
    'conversao': 'geosoft_grd_para_tif.py (grid Geosoft .grd -> GeoTIFF float32, sem reamostragem)',
    'unidades': 'eU ppm, eTh ppm, K %, CT contagem total, razões adimensionais, MAG nT',
    'razoes_limiares': 'U_K e Th_K onde K >= 0.1 %; U_Th onde eTh >= 0.5 ppm',
    'crs': 'EPSG:4326',
}

# séries na ordem de preferência ao tirar duplicatas (a mais detalhada primeiro)
SERIES = ['50k', '100k', '250k', 'regional', 'milionesimo']


# ═══════════════════════════════════════════════════════════════════════════
# UTILITÁRIOS
# ═══════════════════════════════════════════════════════════════════════════
def sem_acento(s) -> str:
    return ''.join(c for c in unicodedata.normalize('NFD', str(s or ''))
                   if unicodedata.category(c) != 'Mn').lower()


def wgs84():
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def dist_irece(lon, lat) -> float:
    return round(math.hypot(lon - IRECE[0], lat - IRECE[1]) * 108, 1)


def na_area(lon, lat) -> int:
    return int(AREA[0] <= lon <= AREA[2] and AREA[1] <= lat <= AREA[3])


def achar(padrao: str) -> Path:
    achados = sorted(DADOS.glob(padrao))
    if not achados:
        raise FileNotFoundError(f'{DADOS}/{padrao}')
    return achados[0]


def zipar_shp(shp: Path) -> Path:
    """Zipa o shapefile e os arquivos auxiliares; apaga os soltos."""
    destino = shp.with_suffix('.zip')
    partes = [p for p in shp.parent.glob(shp.stem + '.*') if p.suffix != '.zip']
    with zipfile.ZipFile(destino, 'w', zipfile.ZIP_DEFLATED) as z:
        for p in partes:
            z.write(p, p.name)
    for p in partes:
        p.unlink()
    return destino


def reprojetar(origem: Path, nome: str) -> Path:
    """Shapefile em EPSG:4326 (atributos como estão) → nome.zip."""
    destino = PASTA_PREP / f'{nome}.shp'
    gdal.VectorTranslate(str(destino), str(origem), format='ESRI Shapefile',
                         dstSRS='EPSG:4326', layerCreationOptions=['ENCODING=UTF-8'])
    return zipar_shp(destino)


def juntar(origens: list[Path], nome: str) -> Path:
    """Vários shapefiles (mesmo tema) num só, em EPSG:4326, com campo 'pacote'."""
    destino = PASTA_PREP / f'{nome}.shp'
    for i, origem in enumerate(origens):
        tmp = PASTA_PREP / '_tmp.shp'
        gdal.VectorTranslate(str(tmp), str(origem), format='ESRI Shapefile', dstSRS='EPSG:4326',
                             SQLStatement=f"SELECT *, '{origem.parent.name[:80]}' AS pacote "
                                          f'FROM "{origem.stem}"', SQLDialect='SQLITE',
                             layerName=destino.stem, layerCreationOptions=['ENCODING=UTF-8'])
        gdal.VectorTranslate(str(destino), str(tmp), format='ESRI Shapefile',
                             accessMode='append' if i else None, addFields=True,
                             layerName=destino.stem, layerCreationOptions=['ENCODING=UTF-8'])
        for p in PASTA_PREP.glob('_tmp.*'):
            p.unlink()
    return zipar_shp(destino)


def novo_shp_pontos(nome: str, campos: list[tuple[str, int]]):
    """(datasource, layer) de um shapefile de pontos em EPSG:4326."""
    caminho = PASTA_PREP / f'{nome}.shp'
    ds = ogr.GetDriverByName('ESRI Shapefile').CreateDataSource(str(caminho))
    ly = ds.CreateLayer(nome, wgs84(), ogr.wkbPoint, options=['ENCODING=UTF-8'])
    for c, tipo in campos:
        fd = ogr.FieldDefn(c, tipo)
        if tipo == ogr.OFTString:
            fd.SetWidth(254)
        ly.CreateField(fd)
    return ds, ly, caminho


def gravar_ponto(ly, lon, lat, valores: dict):
    f = ogr.Feature(ly.GetLayerDefn())
    for k, v in valores.items():
        if v is not None:
            f.SetField(k, v)
    g = ogr.Geometry(ogr.wkbPoint)
    g.AddPoint_2D(lon, lat)
    f.SetGeometry(g)
    ly.CreateFeature(f)


# ═══════════════════════════════════════════════════════════════════════════
# 1) PREPARAR
# ═══════════════════════════════════════════════════════════════════════════
def preparar_rasters() -> None:
    vrt = PASTA_PREP / '_gama.vrt'
    gdal.BuildVRT(str(vrt), [str(GRIDS / f'1042_{b}.tif') for b in GAMA], separate=True)
    gdal.Translate(str(PASTA_PREP / '1042_gama.tif'), str(vrt), noData=NODATA,
                   creationOptions=['COMPRESS=DEFLATE', 'PREDICTOR=3', 'TILED=YES'])
    vrt.unlink()
    shutil.copy(GRIDS / '1042_MAG_nT.tif', PASTA_PREP / '1042_mag.tif')
    print('  1042_gama.tif (7 bandas), 1042_mag.tif')


# prefixos de nome de campo, comparados sem acento e sem caixa. 'DESCRI' pega o
# 'DESCRIÇÃ' do Projeto Morro do Chapéu (nome com acento, cortado em 10 bytes);
# 'RECURSO' é o campo da folha Cabeça no Tempo
CAMPOS_SUBST  = ('SUBS_PRINC', 'SUBST_PRIN', 'SUBSTANCIA', 'SUBSTANCI', 'RECURSO', 'DESCRI')
CAMPOS_STATUS = ('STATUS_ECO', 'TIPO')
CAMPOS_TOPON  = ('TOPONIMIA',)
# camadas que ficam de fora da leitura de ocorrências:
#   'amostragem'                → Mo_Chap_Amostragem_Rec_Min: amostras de geoquímica, não ocorrências
#   'p_mo_chapeu_ocormin_ponto' → cópia exata (mesmas 49 coordenadas) de Mo_Chap_OcorMin_Lito_RecMin,
#                                 mas com a descrição em branco (só a sigla em TEXTO)
CAMADAS_NAO_RECMIN = ('amostragem', 'p_mo_chapeu_ocormin_ponto')


def valor(f, defn, prefixos) -> str:
    """Valor do primeiro campo (na ordem de `prefixos`) cujo nome, sem acento e
    sem caixa, começa com o prefixo e está preenchido."""
    campos = [(sem_acento(defn.GetFieldDefn(i).GetName()).upper(), i)
              for i in range(defn.GetFieldCount())]
    for p in prefixos:
        for nome, i in campos:
            if nome.startswith(p) and f.IsFieldSetAndNotNull(i) and str(f.GetField(i)).strip():
                return str(f.GetField(i)).strip()
    return ''


def eh_fosforo(subst: str) -> bool:
    tokens = [t.strip() for t in re.split(r'[,;/()]', subst)]
    return 'P' in tokens or any(w in sem_acento(subst) for w in ('fosf', 'apatit'))


def ler_recmin() -> list[dict]:
    """Todos os pontos recmin/ocormin de todas as séries, em EPSG:4326."""
    pontos = []
    destino = wgs84()
    for serie in SERIES:
        for shp in sorted((DADOS / serie).rglob('*.shp')):
            s = shp.stem.lower()
            if (not any(t in s for t in ('recmin', 'ocormin', 'rec_min'))
                    or any(t in s for t in CAMADAS_NAO_RECMIN)):
                continue
            ds = ogr.Open(str(shp))
            ly = ds.GetLayer()
            srs = ly.GetSpatialRef()
            if srs:
                srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            ct = osr.CoordinateTransformation(srs, destino) if srs else None
            defn = ly.GetLayerDefn()
            for f in ly:
                g = f.GetGeometryRef()
                if g is None:
                    continue
                g = g.Clone()
                if ct:
                    g.Transform(ct)
                c = g.Centroid()
                pontos.append({
                    'lon': c.GetX(), 'lat': c.GetY(), 'serie': serie,
                    'pacote': shp.parent.name[:80], 'camada': shp.stem[:80],
                    'subst': valor(f, defn, CAMPOS_SUBST)[:254],
                    'status': valor(f, defn, CAMPOS_STATUS)[:254],
                    'toponimia': valor(f, defn, CAMPOS_TOPON)[:254],
                })
    return pontos


def sem_duplicatas(pontos: list[dict]) -> list[dict]:
    """Mesmo lugar (~100 m) e mesma substância = mesma ocorrência; fica a série
    mais detalhada (ordem de SERIES), que é a ordem de leitura."""
    vistos, saida = set(), []
    for p in pontos:
        chave = (round(p['lon'], 3), round(p['lat'], 3), sem_acento(p['subst']))
        if chave not in vistos:
            vistos.add(chave)
            saida.append(p)
    return saida


def preparar_recmin() -> None:
    todos = sem_duplicatas(ler_recmin())
    campos = [('serie', ogr.OFTString), ('pacote', ogr.OFTString), ('camada', ogr.OFTString),
              ('subst', ogr.OFTString), ('status', ogr.OFTString), ('toponimia', ogr.OFTString),
              ('na_area', ogr.OFTInteger), ('dist_irece', ogr.OFTReal)]

    ds, ly, shp = novo_shp_pontos('recursos_minerais', campos)
    for p in todos:
        gravar_ponto(ly, p['lon'], p['lat'], {**{k: p[k] for k in
                     ('serie', 'pacote', 'camada', 'subst', 'status', 'toponimia')},
                     'na_area': na_area(p['lon'], p['lat']), 'dist_irece': dist_irece(p['lon'], p['lat'])})
    ds = None
    zipar_shp(shp)

    # fósforo: a mesma ocorrência às vezes vem com 'P' numa série e 'Fósforo' noutra,
    # deslocada de ~100 m; a menos de ~300 m é a mesma (Lajedão I e II distam ~450 m)
    fosf = []
    for p in todos:
        if eh_fosforo(p['subst']) and not any(
                abs(p['lon'] - q['lon']) < 0.003 and abs(p['lat'] - q['lat']) < 0.003 for q in fosf):
            fosf.append(p)
    ds, ly, shp = novo_shp_pontos('ocorrencias_fosfato', campos + [('grupo', ogr.OFTString)])
    for p in fosf:
        d = dist_irece(p['lon'], p['lat'])
        grupo = 'bacia_irece' if d <= 100 else ('norte' if p['lat'] > -10.0 else 'leste_outros')
        gravar_ponto(ly, p['lon'], p['lat'], {**{k: p[k] for k in
                     ('serie', 'pacote', 'camada', 'subst', 'status', 'toponimia')},
                     'na_area': na_area(p['lon'], p['lat']), 'dist_irece': d, 'grupo': grupo})
    ds = None
    zipar_shp(shp)
    print(f'  recursos_minerais: {len(todos)} pontos | ocorrencias_fosfato: {len(fosf)} pontos '
          f'({sum(dist_irece(p["lon"], p["lat"]) <= 100 for p in fosf)} até 100 km de Irecê)')


def numero(v):
    """'7,3' → 7.3; '<0,01' → (0.01, True); vazio → (None, False)."""
    if v is None or str(v).strip() == '':
        return None, False
    s = str(v).strip().replace(',', '.')
    abaixo = s.startswith('<')
    m = re.search(r'-?\d+(\.\d+)?', s)
    return (float(m.group(0)), abaixo) if m else (None, False)


def preparar_geoquimica() -> None:
    campos = [('projeto', ogr.OFTString), ('amostra', ogr.OFTString), ('classe', ogr.OFTString),
              ('p2o5_pct', ogr.OFTReal), ('p2o5_lt', ogr.OFTInteger),
              ('na_area', ogr.OFTInteger), ('dist_irece', ogr.OFTReal)]
    ds, ly, shp = novo_shp_pontos('geoquimica_p2o5', campos)
    n = 0
    for pasta, projeto in [('geoquimica/7420_*', 'Agrominerais Irece-Jaguarari 2016'),
                           ('geoquimica/6719_*', 'Calcarios da Bahia 2021')]:
        xlsx = next(achar(pasta).glob('*.xlsx'))
        linhas = openpyxl.load_workbook(xlsx, read_only=True).worksheets[0].iter_rows(values_only=True)
        cab = [sem_acento(c) for c in next(linhas)]
        i_lat, i_lon = cab.index('latitude'), cab.index('longitude')
        i_p = cab.index('p2o5_pct')
        i_am = next((i for i, c in enumerate(cab) if c.startswith('numero_de_campo')), None)
        i_cl = next((i for i, c in enumerate(cab) if c == 'classe'), None)
        for r in linhas:
            (lat, _), (lon, _) = numero(r[i_lat]), numero(r[i_lon])
            p2o5, abaixo = numero(r[i_p])
            if lat is None or lon is None or p2o5 is None:
                continue
            gravar_ponto(ly, lon, lat, {
                'projeto': projeto, 'amostra': str(r[i_am]) if i_am is not None else None,
                'classe': str(r[i_cl]) if i_cl is not None else None,
                'p2o5_pct': p2o5, 'p2o5_lt': int(abaixo),
                'na_area': na_area(lon, lat), 'dist_irece': dist_irece(lon, lat)})
            n += 1
    ds = None
    zipar_shp(shp)
    print(f'  geoquimica_p2o5: {n} amostras com P2O5')


def preparar_vetores() -> None:
    reprojetar(achar('250k/2617_*/jacobina_lito.shp'), 'lito_jacobina_250k')
    reprojetar(achar('250k/2617_*/jacobina_estrut.shp'), 'estrut_jacobina_250k')
    reprojetar(achar('regional/5984_*/bahia_lito.shp'), 'lito_bahia')
    reprojetar(achar('regional/5984_*/bahia_estrut.shp'), 'estrut_bahia')
    reprojetar(achar('regional/6691_*/Dolinas_area.shp'), 'dolinas_morro_chapeu')
    print('  lito/estrut Jacobina 250k e Bahia, dolinas Morro do Chapéu')
    if INCLUIR_100K:
        juntar(sorted(DADOS.glob('100k/*/*_lito.shp')), 'lito_100k')
        juntar(sorted(p for p in DADOS.glob('100k/*/*_estrutura.shp')), 'estrut_100k')
        print('  lito/estrut 100k')


def preparar() -> None:
    if PASTA_PREP.exists():
        shutil.rmtree(PASTA_PREP)
    PASTA_PREP.mkdir(parents=True)
    print(f'Preparando em {PASTA_PREP}')
    preparar_rasters()
    preparar_recmin()
    preparar_geoquimica()
    preparar_vetores()


# ═══════════════════════════════════════════════════════════════════════════
# 2) ENVIAR (Cloud Storage → ingestão no GEE)
# ═══════════════════════════════════════════════════════════════════════════
def manifesto_imagem(asset: str, uri: str, bandas: list[str], pyramiding: str) -> dict:
    import ee
    return {
        'name': ee._cloud_api_utils.convert_asset_id_to_asset_name(asset),
        'tilesets': [{'id': 't0', 'sources': [{'uris': [uri]}]}],
        'bands': [{'id': b, 'tilesetId': 't0', 'tilesetBandIndex': i} for i, b in enumerate(bandas)],
        'missingData': {'values': [NODATA]},
        'pyramidingPolicy': pyramiding,
        'properties': PROPS_1042,
    }


def manifesto_tabela(asset: str, uri: str, descricao: str) -> dict:
    import ee
    return {
        'name': ee._cloud_api_utils.convert_asset_id_to_asset_name(asset),
        'sources': [{'uris': [uri], 'charset': 'UTF-8'}],
        'properties': {'fonte': 'SGB/CPRM Acervo GEOSGB', 'descricao': descricao,
                       'preparo': 'exportar_acervo_gee.py (EPSG:4326, sem recorte)'},
    }


DESCRICOES = {
    'ocorrencias_fosfato': 'ocorrencias de fosforo (P/Fosforo) de todos os recmin, sem duplicatas; '
                           'grupo bacia_irece = ate 100 km de (-41.86,-11.30)',
    'recursos_minerais': 'todos os recmin/ocormin baixados, sem duplicatas (~100 m + substancia)',
    'geoquimica_p2o5': 'rocha: Agrominerais Irece-Jaguarari (7420) + Calcarios da Bahia (6719); '
                       'p2o5_lt = 1 se abaixo do limite de deteccao',
    'lito_jacobina_250k': 'litologia folha Jacobina SC.24-Y-C 1:250.000 (pacote 2617)',
    'estrut_jacobina_250k': 'estruturas folha Jacobina SC.24-Y-C 1:250.000 (pacote 2617)',
    'lito_bahia': 'litologia Geologia do Estado da Bahia 1:1.000.000 (pacote 5984)',
    'estrut_bahia': 'estruturas Geologia do Estado da Bahia 1:1.000.000 (pacote 5984)',
    'dolinas_morro_chapeu': 'dolinas do Projeto Morro do Chapeu 1995 (pacote 6691)',
    'lito_100k': 'litologia das cartas 1:100.000 que cruzam a area (pacotes 100k)',
    'estrut_100k': 'estruturas das cartas 1:100.000 que cruzam a area (pacotes 100k)',
}


def perguntar(texto: str, padrao: str = '') -> str:
    resp = input(f'{texto}' + (f' [{padrao}]' if padrao else '') + ': ').strip()
    return resp or padrao


def escolher_conta() -> str:
    print('\nContas (todas escrevem em mapbiomas-arida/mine):')
    for i, c in enumerate(CONTAS, 1):
        print(f'  {i}) {c}')
    while True:
        r = perguntar('Conta', '4')
        if r.isdigit() and 1 <= int(r) <= len(CONTAS):
            return CONTAS[int(r) - 1]
        print('  → número da lista.')


def ativar(conta: str) -> None:
    import ee
    from configure_account_projects_ee import get_project_from_account
    from gee_tools import switch_user
    arq = Path.home() / '.config' / 'earthengine' / 'current_user.txt'
    if arq.read_text().strip() != conta:
        switch_user(conta)
    ee.Initialize(project=get_project_from_account(conta))


def cliente_storage(projeto: str):
    """Storage com as credenciais da conta do EE; se não der, as do gcloud (ADC).

    As requisições são cobradas do `projeto` (o dono do bucket), qualquer que seja
    o quota project gravado nas credenciais: o ADC costuma vir com o projeto do EE,
    que não tem faturamento, e aí o upload falha com 'billing account … absent'."""
    import ee
    import google.auth
    from google.cloud import storage

    def no_projeto(cred):
        return cred.with_quota_project(projeto) if hasattr(cred, 'with_quota_project') else cred

    try:
        cliente = storage.Client(project=projeto,
                                 credentials=no_projeto(ee.data.get_persistent_credentials()))
        list(cliente.list_buckets(max_results=1))
        return cliente
    except Exception as e:    # noqa: BLE001 — sem escopo de Storage nas credenciais do EE
        print(f'  credenciais do EE sem acesso ao Storage ({str(e)[:80]}); tentando gcloud ADC')
        cred, _ = google.auth.default(scopes=['https://www.googleapis.com/auth/devstorage.read_write'])
        return storage.Client(project=projeto, credentials=no_projeto(cred))


def enviar() -> None:
    import ee
    from configure_account_projects_ee import get_project_from_account

    arquivos = sorted(PASTA_PREP.glob('*.tif')) + sorted(PASTA_PREP.glob('*.zip'))
    if not arquivos:
        raise SystemExit(f'Nada em {PASTA_PREP}: rode a etapa "preparar" antes.')

    bucket_uri = perguntar('Bucket do Cloud Storage (gs://bucket/pasta)', BUCKET).rstrip('/')
    if not bucket_uri.startswith('gs://'):
        raise SystemExit('Bucket precisa começar com gs://')
    bucket_nome, _, prefixo = bucket_uri[5:].partition('/')

    conta = escolher_conta()
    ativar(conta)
    sobrescrever = perguntar('Se o asset já existir, sobrescrever? (s/n)', 'n').lower().startswith('s')

    try:
        ee.data.getAsset(DESTINO)
    except ee.EEException:
        ee.data.createFolder(DESTINO)
        print(f'Pasta criada: {DESTINO}')

    bucket = cliente_storage(PROJETO_STORAGE or get_project_from_account(conta)).bucket(bucket_nome)
    registro = []
    for arq in arquivos:
        nome = arq.stem
        asset = f'{DESTINO}/{nome}'
        try:
            ee.data.getAsset(asset)
            existe = True
        except ee.EEException:
            existe = False
        if existe and not sobrescrever:
            print(f'  {nome}: já existe — pulado')
            registro.append({'asset': asset, 'status': 'existe'})
            continue

        blob = f'{prefixo}/{arq.name}' if prefixo else arq.name
        bucket.blob(blob).upload_from_filename(str(arq))
        uri = f'gs://{bucket_nome}/{blob}'

        if nome == '1042_gama':
            manifesto = manifesto_imagem(asset, uri, GAMA, 'MEAN')
            inicio = ee.data.startIngestion
        elif nome == '1042_mag':
            manifesto = manifesto_imagem(asset, uri, ['MAG_nT'], 'MEAN')
            inicio = ee.data.startIngestion
        else:
            manifesto = manifesto_tabela(asset, uri, DESCRICOES.get(nome, nome))
            inicio = ee.data.startTableIngestion
        try:
            resp = inicio(ee.data.newTaskId()[0], manifesto, allow_overwrite=sobrescrever)
        except ee.EEException as e:
            print(f'  {nome}: ERRO ao iniciar a ingestão: {e}')
            registro.append({'asset': asset, 'uri': uri, 'status': f'erro: {str(e)[:100]}'})
            continue
        print(f'  {nome}: {uri} → {asset}  (task {resp.get("id")})')
        registro.append({'asset': asset, 'uri': uri, 'status': 'ingerindo', 'task_id': resp.get('id')})

    arq_reg = PASTA / f'registro_exportar_gee_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq_reg, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['asset', 'uri', 'status', 'task_id'], extrasaction='ignore')
        w.writeheader()
        w.writerows(registro)
    print(f'Conta {conta}. Registro → {arq_reg}')


# ═══════════════════════════════════════════════════════════════════════════
def main():
    print('══════ Acervo CPRM → assets GEE ══════')
    print(f'  origem : {DADOS}\n  preparo: {PASTA_PREP}\n  destino: {DESTINO}')
    print('Etapas:\n  1) só preparar os arquivos locais\n  2) só enviar (já preparados)\n  3) preparar e enviar')
    etapa = perguntar('Etapa', '1')
    if etapa in ('1', '3'):
        preparar()
    if etapa in ('2', '3'):
        enviar()


if __name__ == '__main__':
    main()
