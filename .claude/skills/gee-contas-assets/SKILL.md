---
name: gee-contas-assets
description: Use when the user needs multi-account Earth Engine handling for this project — choosing/rotating accounts for exports, checking read/write access to mapbiomas-arida vs mapbiomas-brazil, reviewing tasks per account, copying assets between projects, listing folders with many small tables (asset quota), or deleting assets. Triggers on "contas", "rodízio", "superconta lotada", "revisarTasks", "transferir assets", "copyAsset", "renameAsset", "cota de assets", "apagar assets", "deletfiles".
---

# Contas e assets GEE

Ver [SPEC.md](SPEC.md).

## Contas (`src/configure_account_projects_ee.py`, `src/gee_tools.switch_user`)

caatinga01–05, solkan1201, solkanGeodatin, solkanCengine
(`ee-solkancengine17`), superconta (mapbiomas-brazil).
- Todas exceto superconta leem e escrevem `projects/mapbiomas-arida/mine`.
- Superconta NÃO lê nem escreve mapbiomas-arida → nunca misturar com as
  outras num rodízio; scripts v2 a excluem.
- `mapbiomas-brazil/.../CAATINGA` está `all_users_can_read`.
- Conta ativa: `~/.config/earthengine/current_user.txt`.
- Projeto `mapbiomas-brazil` como projeto EE (2026-10-17): caatinga01–05,
  solkan1201 e superconta inicializam; solkanCengine falha (serviceusage);
  solkanGeodatin inicializa mas não lê arida. Usado no rodízio de
  `build_mosaico_unico.py`.
- Teste só leitura sem trocar conta: `Credentials(None, refresh_token=...,
  client_id=d.get('client_id') or ee.oauth.CLIENT_ID, client_secret=...,
  token_uri=...)` a partir do arquivo da conta.

- caatinga04 no projeto mapbiomas-brazil (cota/billing de lá) lê e escreve
  `mine/`: usada em `harmonizar_sazonal_modis.py` e na cópia de `ASTER_sazonal`.
- Exports de tabela: `toAsset` recusa feição com geometria nula ("Unable to export
  features with null geometry"); `toDrive` CSV aceita. Tabela com muitas cenas numa
  task só estoura memória: dividir (por mês).

## Ferramentas (`src/tools/`)

| Script | Uso |
|---|---|
| `revisarTasksXconta.py` | tasks por conta |
| `transferir_assets_brazil_para_arida.py` | conferir / copiar (`copyAsset`, caatinga04) / apagar originais (superconta, digita APAGAR); `COLECOES` ativa hoje: só `ASTER_sazonal` (as outras transferidas em 2026-10-12) |
| `listar_pastas_tabelas_gee.py` | cota + pastas com muitas tabelas, CSV; não apaga |
| `deletfilesFOLDERs.py` | apaga filhos de uma pasta (try/except por asset) |

## Regras

- Copiar, não mover: `renameAsset` apaga a origem e exige escrita lá.
- mapbiomas-brazil estoura pela QUANTIDADE de assets (200 000), não pelo
  tamanho; cada imagem de ImageCollection conta 1.
- Apagar: sempre confirmar com o usuário o caminho exato antes; mostrar o
  que será apagado. Pastas de outras equipes no mapbiomas-brazil.
- Assets "fantasmas" (listados mas `not found` para todas as contas, ex.
  `ALERTAS`): apagar a pasta pai ou abrir issue no Google; não insistir.
- Scripts de export multi-conta deduplicam contra assets salvos e tasks
  READY/RUNNING antes de enviar.
