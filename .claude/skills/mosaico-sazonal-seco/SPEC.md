# Spec — mosaico-sazonal-seco

## Objetivo
Manter o fluxo do mosaico sazonal seco ASTER (períodos A e B) — construção, cópia,
visualização, harmonização com o MODIS e insumos para estimar o SWIR do B — e o
checklist de andamento.

## Pré-condições
- Grid `projects/mapbiomas-arida/mine/grid_irece_5x5_solape100` legível pela conta.
- `lista_imagens_Notaceites.txt` e `list_zero_cloud.txt` iguais aos do export: os
  mapas de cenas dependem disso (`confere` ≈ 1).
- Harmonização: mosaicos em `mine/ASTER_sazonal`; caatinga04 com uso do projeto
  mapbiomas-brazil; para a etapa 3, CSVs baixados do Drive em
  `src/Dados/ASTER_sazonal_harmonizacao/ASTER_sazonal_harmonizacao/`.

## Entradas
| Campo | Origem |
|---|---|
| etapa, períodos, células | menus interativos dos scripts |
| sobrescrever (build) | usuário, digita APAGAR |
| CSVs de pares | download manual do Drive pelo usuário |
| parâmetros de ajuste (`CV_MAX`, `FRAC_MIN`, `DIF_*`, `N_MIN`, `R2_MIN`, `GANHO_FAIXA`, faixas do sorteio) | topo do script |

## Saídas
| Etapa | Saída |
|---|---|
| build | `.../CAATINGA/ASTER_sazonal/*` + `registro_sazonal_seco_*.csv` |
| cópia | `mine/ASTER_sazonal/*` + `registro_transferencia_*.csv` |
| harm. 1 | `mine/ASTER_sazonal_cenas/*` (`cena_t`, `confere`) |
| harm. 2 | `Drive/ASTER_sazonal_harmonizacao/<P>_<ano>_<mês>.csv` |
| harm. 3 | `src/Dados/ASTER_sazonal_harmonizacao/coef_modis_<P>.csv` |
| harm. 4 | `mine/ASTER_sazonal_harmonizado/*` (+ `cena_t`, `harmonizado`) |
| checklist | `src/building_dataset/CHECKLIST_sazonal_seco.md` (editado por Claude e pelo usuário) |

## Efeitos colaterais
Tasks GEE (toAsset/toDrive) e criação de ImageCollections — só o usuário dispara.
Os scripts pulam o que já existe ou está READY/RUNNING; não apagam assets (exceto o
build com sobrescrever + APAGAR). Os originais nunca são alterados pela harmonização.

## Modos de falha e como reagir
- `confere` < 0,99 numa célula: listas mudaram depois do export → mapa inválido;
  refazer o mosaico ou aceitar a mediana nessa célula (perguntar).
- Task de pares com out of memory: dividir o mês (por cena) — não voltar ao getInfo
  com máscara de nuvem.
- Poucas cenas com reta válida ou ganho fora de 0,8–1,2: revisar `CV_MAX`, `DIF_*` e
  a faixa do sorteio (`FAIXA_CLASSE`; a Caatinga quase não passa de 4000 no vermelho).
- Emenda que sobra depois de harmonizar = fenologia: só remontar resolve; avisar.
- superconta não lê mapbiomas-arida: scripts que leem `mine/` usam caatinga04.

## Não-objetivos
- Não roda scripts nem dispara tasks.
- Não mistura A harmonizado com B original (escalas diferentes: NBAR × TOA).
- Não corrige o TIR (MCD43A4 não tem termal).
