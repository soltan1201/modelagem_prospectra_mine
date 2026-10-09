# Spec — sample-mosaico-dataset

## Objetivo
Rodar/editar `src/building_dataset/sample_quality_mosaic.py`, que amostra
pontos aleatórios dos mosaicos semestrais ASTER (harmonizados a 30m) e os
exporta para Drive (CSV) e/ou asset GEE (FeatureCollection).

## Pré-condições
- `ee` (earthengine-api) instalado e autenticado para o projeto
  `mapbiomas-caatinga-cloud02` (`ee.Initialize(project=...)` deve suceder).
- Asset `projects/mapbiomas-arida/mosaic_aster` existe e contém ao menos um
  mosaico semestral com as bandas `B01-B14` + `quality`.
- Se `EXPORTAR_PARA` incluir `'asset'`: usuário tem permissão de escrita em
  `projects/mapbiomas-arida/mine/points`.

## Entradas
| Campo | Origem | Obrigatório |
|---|---|---|
| `EXPORTAR_PARA` (`drive`/`asset`/`both`) | usuário/uso pretendido | sim |
| `N_PONTOS`, `SEED`, `ESCALA_AMOST` | usuário (default preservado se omitido) | não |

## Saídas
- Uma `Export.table` task por imagem semestral, iniciada
  (`task.start()`), destino Drive (`ASTER_Samples/<system:index>.csv`)
  e/ou asset (`projects/mapbiomas-arida/mine/points/samples_v2_<system:index>`).
- Resumo impresso no console: total de tasks iniciadas / total de imagens.

## Efeitos colaterais
- **Cria tasks assíncronas reais no GEE** (consome cota de export do
  projeto) — isto NÃO é uma operação local reversível. Confirme com o
  usuário antes de rodar para todas as imagens se a intenção era testar com
  uma única imagem (o script atual não tem flag de "dry run" ou "apenas N
  imagens" — se o usuário quiser testar, seria necessário fatiar
  `lista_imgs` manualmente, ex.: `range(1)` em vez de `range(n_total)`).

## Modos de falha e como reagir
- **`ee.Initialize` falha**: verificar autenticação
  (`earthengine authenticate`) antes de qualquer outra coisa — não
  prosseguir tentando contornar.
- **Asset de pontos já existe** (reexecução): `Export.table.toAsset` falha
  se o assetId já existir — avisar o usuário, sugerir sufixo de versão ou
  deletar o asset antigo (ação destrutiva — só com confirmação explícita).
- **Quantidade grande de imagens**: cada imagem gera uma task separada;
  para uma coleção grande isso pode ser dezenas de exports simultâneos —
  avisar sobre cota antes de rodar.

## Não-objetivos
- Não aplica filtro de qualidade nos pontos amostrados — isso é
  propositalmente adiado para pós-processamento (pandas/R local, ou o
  `QUALITY_MIN` do skill `estimate-swir-gee`).
- Não converte INT16 → unidade física — os CSVs saem em INT16 bruto; a
  conversão (`× 1e-4` reflectância, `× 0.1` TIR) é feita nos scripts de
  modelagem (`train_swir_models.py` / `estimate_swir_gee.py`).
