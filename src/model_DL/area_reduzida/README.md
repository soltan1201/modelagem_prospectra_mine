# CNN — área reduzida (Irecê, ASTER + CPRM)

Pasta reservada para a CNN da área reduzida. Ainda não há scripts.

Diferenças em relação a `../regiao_completa/`:

- Pseudo-rótulos de `projects/mapbiomas-arida/mine/favorabilidade_AR`
  (`<RF|GTB>_fosfato_AR_irece_shm_selecao_fcl`), não de `favorabilidade_ml`.
- Canais: mosaico harmonizado + variáveis CPRM (gama, magnetometria, geologia),
  conforme `src/modeling/models/variaveis_AR_irece_shm.json`.
- Área de `src/modeling/area_reduzida/js/novas_areas_estudo.js`, não a grid 5 × 5.

Os `.md` de plano e rodadas ficam em `src/model_DL/`.
