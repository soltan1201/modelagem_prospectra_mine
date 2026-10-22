"""
Sobe os mapas de p da CNN (inferir_cnn.py --tiles) para o GEE: envia cada
p_<LxCy>.tif ao bucket do Cloud Storage e cria a ingestão como imagem numa
ImageCollection, no mesmo formato do favorabilidade_ml (banda 'probability' =
p × 10000, 0 = sem dado) para comparar com RF/GTB no visualizar_favorabilidade_ml.js.

Entrada   <run>/mapas_passo<k>/p_<LxCy>.tif   (UInt16, EPSG:4326)
Destino   COLECAO/<rodada>_p<k>_<LxCy>  com propriedades modelo (CNN_<backbone>),
          alvo, mosaico, celula, rodada, passo, canais (A/B1), escala_p
Bucket    BUCKET (o mesmo do acervo); cobrança no PROJETO_STORAGE — credenciais
          como no dowloadQGIS/exportar_acervo_gee.py (EE ou gcloud ADC; ver
          dowloadQGIS/GUIA_GCLOUD_BUCKET_GEE.md)

Pula imagens que já existem na coleção. Registro em <run>/registro_subida_gee_<data>.csv.

Uso
  python subir_predicao_gee.py --run /srv/almacen/db_images/runs_cnn/A_regnety_64_teste --passo 1
"""

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import ee

COLECAO         = 'projects/mapbiomas-arida/mine/cnn_favorabilidade'
BUCKET          = 'gs://cengine17-acervo-cprm/cnn_favorabilidade'
PROJETO_STORAGE = 'project-fadb9cba-032b-4fb9-882'   # mesmo projeto de cobrança do acervo
PROJETO_EE      = 'mapbiomas-brazil'
ALVO, MOSAICO   = 'fosfato', 'aster-mosaico-final'


def cliente_storage():
    """Storage com as credenciais do EE; se não der, as do gcloud (ADC). Cobrança sempre no
    PROJETO_STORAGE (sem isso: 'billing account … absent') — igual ao exportar_acervo_gee.py."""
    import google.auth
    from google.cloud import storage

    def no_projeto(cred):
        return cred.with_quota_project(PROJETO_STORAGE) if hasattr(cred, 'with_quota_project') else cred

    try:
        cliente = storage.Client(project=PROJETO_STORAGE,
                                 credentials=no_projeto(ee.data.get_persistent_credentials()))
        list(cliente.list_buckets(max_results=1))
        return cliente
    except Exception as e:    # noqa: BLE001 — credenciais do EE sem escopo de Storage
        print(f'  credenciais do EE sem acesso ao Storage ({str(e)[:80]}); tentando gcloud ADC')
        cred, _ = google.auth.default(scopes=['https://www.googleapis.com/auth/devstorage.read_write'])
        return storage.Client(project=PROJETO_STORAGE, credentials=no_projeto(cred))


def manifesto(asset, uri, props):
    return {
        'name': ee._cloud_api_utils.convert_asset_id_to_asset_name(asset),
        'tilesets': [{'id': 't0', 'sources': [{'uris': [uri]}]}],
        'bands': [{'id': 'probability', 'tilesetId': 't0', 'tilesetBandIndex': 0}],
        'missingData': {'values': [0]},
        'pyramidingPolicy': 'MEAN',
        'properties': props,
    }


def main(lista=None):
    a = argparse.ArgumentParser(description='Sobe os mapas de p da CNN para uma ImageCollection')
    a.add_argument('--run', required=True)
    a.add_argument('--passo', type=int, default=1)
    a.add_argument('--colecao', default=COLECAO)
    a.add_argument('--bucket', default=BUCKET)
    a.add_argument('--projeto', default=PROJETO_EE)
    a.add_argument('--sobrescrever', action='store_true')
    args = a.parse_args(lista)

    run = Path(args.run)
    pasta = run / f'mapas_passo{args.passo}'
    arquivos = sorted(pasta.glob('p_L*C*.tif'))
    if not arquivos:
        raise SystemExit(f'Nenhum p_<LxCy>.tif em {pasta} — rode inferir_cnn.py --tiles … --passo {args.passo}')
    cfg = json.loads((run / 'config.json').read_text())
    rodada = run.name

    ee.Initialize(project=args.projeto)
    try:
        ee.data.getAsset(args.colecao)
    except ee.EEException:
        ee.data.createAsset({'type': 'IMAGE_COLLECTION'}, args.colecao)
        print(f'Coleção criada: {args.colecao}')

    bucket_nome, _, prefixo = args.bucket[5:].partition('/')
    bucket = cliente_storage().bucket(bucket_nome)
    registro = []
    for arq in arquivos:
        cel = arq.stem.removeprefix('p_')
        nome = f'{rodada}_p{args.passo}_{cel}'.replace('.', '_')
        asset = f'{args.colecao}/{nome}'
        try:
            ee.data.getAsset(asset)
            if not args.sobrescrever:
                print(f'  {nome}: já existe — pulado')
                registro.append({'asset': asset, 'status': 'existe'})
                continue
        except ee.EEException:
            pass
        blob = f'{prefixo}/{rodada}/{arq.name}' if prefixo else f'{rodada}/{arq.name}'
        bucket.blob(blob).upload_from_filename(str(arq))
        uri = f'gs://{bucket_nome}/{blob}'
        props = {'modelo': f'CNN_{cfg.get("backbone")}', 'alvo': ALVO, 'mosaico': MOSAICO, 'celula': cel,
                 'rodada': rodada, 'passo': args.passo, 'canais': cfg.get('modelo'),
                 'patch': cfg.get('patch'), 'escala_p': 10000}
        try:
            resp = ee.data.startIngestion(ee.data.newTaskId()[0], manifesto(asset, uri, props),
                                          allow_overwrite=args.sobrescrever)
            print(f'  {nome}: {uri} → task {resp.get("id")}')
            registro.append({'asset': asset, 'uri': uri, 'status': 'ingerindo', 'task_id': resp.get('id')})
        except ee.EEException as e:
            print(f'  {nome}: ERRO {e}')
            registro.append({'asset': asset, 'uri': uri, 'status': f'erro: {str(e)[:100]}'})

    arq_reg = run / f'registro_subida_gee_{datetime.now():%Y%m%d_%H%M}.csv'
    with open(arq_reg, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['asset', 'uri', 'status', 'task_id'], extrasaction='ignore')
        w.writeheader()
        w.writerows(registro)
    print(f'Registro → {arq_reg}')


if __name__ == '__main__':
    main()
