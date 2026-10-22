# Spec — build-mosaico-unico

## Objetivo
Editar/planejar `src/building_dataset/build_mosaico_unico.py`: mosaicos base
ASTER por célula separados por SWIR (etapa 1), SWIR estimado por RF no
mosaico sem SWIR (etapa 2) e mosaico único juntando os dois (etapa 3).

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
- Etapa 3: `.../mine/ASTER_mosaico_unico_final/ASTER_mosaico_unico_<LxCy>`
  (B01–B14, quality, data, camada, swir_estimado, normalizacao) e
  `retas_mosaico_unico_<AAAAMMDD_HHMM>.csv` (retas por célula e cena).
- CSV `registro_mosaico_unico_etapa<N>_<AAAAMMDD_HHMM>.csv`.
- Propriedade `versao = 'unico_v2'`.

## Efeitos colaterais
Tasks `Export.image.toAsset` reais, rodízio entre contas. Só o usuário
dispara; Claude pode rodar testes só leitura (`getInfo`, contagens).

## Modos de falha
- Célula sem cenas (`zero_cloud + preenchimento == 0`): status `sem_cenas`
  no registro, sem export.
- Etapa 2 sem `ASTER_base_sem_swir_<célula>`: status `sem_base`, pulada.
  Causa usual: etapa 1 rodada com "só com_swir" ou tasks não terminadas.
  Rodar etapa 1 com sem_swir antes; não contornar na etapa 2.
- Task duplicada de asset já salvo (enviada antes de outra conta
  terminar): falha com "already exists". O script só deduplica contra
  tasks READY/RUNNING das contas escolhidas na rodada; conferir as demais
  contas antes de reenviar.
- Conta sem acesso: descartada pelo script; não contornar.
- Filtro `ee` no nível de módulo quebra o import (antes de `ee.Initialize`).
- `ee.Filter.And()` vazio é inválido — usar `nos_meses()`.

- Etapa 3 sem um dos bases: `sem_base`; sem pares estáveis: `sem_retas`.
- Retas calculadas dentro da imagem: "Too many concurrent aggregations" —
  manter o cálculo por `getInfo` antes do export.

## Não-objetivos
- A etapa 3 não usa o `ASTER_sem_swir_estimado` da etapa 2 (re-estima).
- Não altera limiares de nuvem/quality de `build_mosaicos_semestrais_v2`.
