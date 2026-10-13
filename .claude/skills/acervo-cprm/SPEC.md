# Spec — acervo-cprm

## Objetivo
Manter o fluxo Acervo SIG CPRM → disco local → GCS → assets GEE para a
área de Irecê, e o documento de decisão `ACERVO_CPRM_IRECE.md`.

## Pré-condições
- Downloads: rede para o servidor do SGB; disco `/mnt/Almacen/shpIRECE/`.
- Export: `gcloud auth application-default login` feito com
  solkan.cengine17; bucket existe; conta EE escreve em
  `projects/mapbiomas-arida/mine`.

## Entradas
| Campo | Origem |
|---|---|
| seleção de camadas/cartas | usuário + regras de `baixar_acervo_irece.py` |
| conta EE da ingestão | usuário (pergunta ao rodar) |
| bucket | constante `BUCKET` ou pergunta |

## Saídas
- Arquivos em `/mnt/Almacen/shpIRECE/` por série (`50k`, `100k`, `250k`,
  `regional`, `milionesimo`), `aerogeofisica/1042_grids`, `manifesto.csv`.
- `gee_upload/` preparado; objetos no bucket; assets em
  `projects/mapbiomas-arida/mine/acervo_cprm`.
- `registro_exportar_gee_<data>.csv`.

## Efeitos colaterais
Uploads ao GCS (custo no crédito do teste) e tasks de ingestão GEE —
o usuário dispara.

## Modos de falha
- 403 billing no bucket/upload: projeto de cota errado → `PROJETO_STORAGE`.
- Campos com acento truncado no DBF: casar por prefixo.
- Camada duplicada de ocorrência: excluir em `CAMADAS_NAO_RECMIN`.

## Não-objetivos
- Mapas PDF não entram (sem uso no modelo).
- Não decide referências × covariáveis sem o usuário (seção 7 do md).
