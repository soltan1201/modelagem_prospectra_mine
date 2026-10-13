# Spec — gee-contas-assets

## Objetivo
Operar com segurança as várias contas GEE do projeto: escolher contas de
export, mover dados entre projetos, medir e liberar cota de assets.

## Pré-condições
- Credenciais de cada conta salvas (via `configure_account_projects_ee`).
- Para copiar: conta lê origem e escreve destino.

## Entradas
| Campo | Origem |
|---|---|
| contas | usuário (menus interativos) |
| raiz/pasta alvo | usuário |
| confirmação de apagar | usuário, explícita (digitar APAGAR / flag) |

## Saídas
- CSVs em `src/tools/` (`pastas_*`, `tabelas_*`, `registro_transferencia_*`).
- Assets copiados/apagados conforme a etapa.

## Efeitos colaterais
`copyAsset`, `deleteAsset`, `createAsset` são reais. Apagar não tem volta.

## Modos de falha
- `EEException` "not found"/permission no delete: registrar e seguir.
- Superconta em rodízio com mapbiomas-arida: falha de acesso — remover.
- Cota do projeto cheia: listar pastas e propor (não executar) limpeza.

## Não-objetivos
- Não apaga nada sem confirmação explícita do usuário para aquele caminho.
- Não decide o que apagar em pastas de outras equipes.
