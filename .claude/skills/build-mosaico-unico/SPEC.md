# Spec — build-mosaico-unico

## Objetivo
Editar/planejar `src/building_dataset/build_mosaico_unico.py`: mosaicos base
ASTER por célula separados por SWIR (etapa 1) e SWIR estimado por RF no
mosaico sem SWIR (etapa 2).

## Pré-condições
- `ee` inicializado com conta que lê ASTER L1T, `grid_irece_5x5_solape100`,
  `points_zeroCC_00_08` e escreve em `projects/mapbiomas-arida/mine`
  (superconta NÃO serve).
- Listas: `lista_imagens_Notaceites.txt`, lista zero-cloud ativa em
  `src/show_mosaics/visualizar_zero_cloud.js` (ou `list_zero_cloud.txt`).
- Etapa 2: `ASTER_base_sem_swir_<célula>` já existe;
  `models/zero_cloud/normalizacao_ols_cena_coef.csv` e CSV de pares legíveis.

## Entradas (interativas)
| Campo | Padrão |
|---|---|
| etapa | 1 |
| contas | lista de `b.perguntar_contas` |
| tipos | com_swir + sem_swir |
| meses | jul–out (`MESES_SECOS`) |
| células | todas |
| conferir / exportar | conferir |

## Saídas
- Etapa 1: `.../mine/ASTER_mosaico_unico_base/ASTER_base_{com,sem}_swir_<LxCy>`.
- Etapa 2: `.../mine/ASTER_mosaico_unico/ASTER_sem_swir_estimado_<LxCy>`.
- CSV `registro_mosaico_unico_etapa<N>_<AAAAMMDD_HHMM>.csv`.
- Propriedade `versao = 'unico_v2'`.

## Efeitos colaterais
Tasks `Export.image.toAsset` reais, rodízio entre contas. Só o usuário
dispara; Claude pode rodar testes só leitura (`getInfo`, contagens).

## Modos de falha
- Célula sem cenas após blacklist: pulada (`apos_blacklist == 0`).
- Conta sem acesso: descartada pelo script; não contornar.
- Filtro `ee` no nível de módulo quebra o import (antes de `ee.Initialize`).
- `ee.Filter.And()` vazio é inválido — usar `nos_meses()`.

## Não-objetivos
- Não combina com_swir e sem_swir_estimado (decisão em aberto).
- Não altera limiares de nuvem/quality de `build_mosaicos_semestrais_v2`.
