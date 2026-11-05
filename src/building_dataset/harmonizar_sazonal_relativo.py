"""
Segunda passada da harmonização dos mosaicos sazonais secos: puxa células que
continuaram destoando depois de harmonizar_sazonal_modis.py (ex.: coluna 1, linhas 1
e 2 — cenas mais verdes que as vizinhas) para o padrão das células boas.

Por que a 1ª passada não resolveu: a reta de cada cena é ajustada contra o MCD43A4
do MESMO DIA. Isso tira atmosfera e calibração, mas não fenologia — se a cena é de
um junho ainda verde, o MODIS daquele dia também está verde e a reta não muda nada.

Ideia: usar as células boas (ex.: colunas 2 e 3) como referência, com o MODIS
mediano da estação seca (vários anos) como ponte entre as células, que não se
sobrepõem:
  1) nas células de REFERÊNCIA, a 500 m, ajusta-se por banda
       ASTER_harmonizado = a × MODIS_mediano + b
     ou seja, "como um mosaico bom se relaciona com o seco típico";
  2) em cada CENA das células ALVO, o valor esperado é z = a × MODIS_mediano + b
     (o que a referência teria ali) e ajusta-se z = ganho × ASTER_harmonizado + offset
     em pixels homogêneos da cena;
  3) aplica-se a reta da cena só nas células alvo.
Diferença real de cobertura (mata mais densa, lagoa) continua, porque o MODIS
mediano a enxerga; o que sai é o desvio da data da cena em relação ao seco típico.

Pares de bandas iguais aos da 1ª passada (B01↔4, B02↔1, B3N↔2, B04↔6, B05–B09 um
ganho só, ajustado em B05↔7). TIR não muda.

Etapas (perguntadas ao rodar):
  1) amostras → Drive: uma task por período, com pares das células de referência e
     das alvo → Drive/PASTA_DRIVE/rel_<P>.csv. Baixe para PASTA_PONTOS.
  2) retas (local): relação da referência e reta por cena alvo → coef_rel_<P>.csv
  3) aplicar: células alvo → COL_REL/<mesmo nome> (bandas do harmonizado + ajuste_rel:
     1 = reta da cena, 2 = mediana das cenas alvo). As demais células ficam só em
     COL_HARM; o visualizador junta as duas coleções.

Conta: caatinga04, projeto mapbiomas-brazil (igual a harmonizar_sazonal_modis.py).

Uso: python harmonizar_sazonal_relativo.py
"""

import logging
import statistics

import ee
import numpy as np
import pandas as pd

import build_mosaico_sazonal_seco as s
import harmonizar_sazonal_modis as h
from build_mosaicos_semestrais_v2 import (ler_grid, listar_destino, perguntar_opcao,
                                          perguntar_sim_nao)
from gee_tools import switch_user

log = logging.getLogger(__name__)

# ═══════════════════════════════════════════════════════════════════════════
# PARÂMETROS
# ═══════════════════════════════════════════════════════════════════════════
ALVOS       = ['L1C1', 'L2C1']                          # células a corrigir (padrão)
REFERENCIAS = ['L1C2', 'L1C3', 'L2C2', 'L2C3']          # células boas (padrão)
COL_REL     = f'{h.PASTA_MINE}/ASTER_sazonal_harm_rel'
PASTA_DRIVE = 'ASTER_sazonal_harmonizacao_rel'
PASTA_PONTOS = h.SAIDA / 'relativo'                     # CSVs baixados de Drive/PASTA_DRIVE

N_REF_CLASSE  = 150    # pixels por faixa do vermelho MODIS em cada célula de referência
N_ALVO_CENA   = 300    # pixels por cena em cada célula alvo
N_MIN_MODIS   = 3      # observações QA = 0 mínimas no MODIS mediano
N_MIN_REL     = 30     # pixels mínimos para aceitar a reta de uma cena alvo
R2_MIN_REL    = 0.5


def chaves(periodo: str) -> list[str]:
    return h.CHAVES[periodo]


def bandas_x(periodo: str) -> list[str]:
    return s.VNIR + (['B04', 'B05'] if periodo == 'A' else [])


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 1 — AMOSTRAS
# ═══════════════════════════════════════════════════════════════════════════
def modis_mediano(periodo: str) -> ee.Image:
    """MCD43A4 mediano da estação seca do período, só observações com QA = 0."""
    cfg = s.CONFIG[periodo]
    nums = sorted({h.AJUSTES[k][1] for k in chaves(periodo)} | {h.CLASSE_MODIS})
    col = (ee.ImageCollection(h.MODIS).filterDate(cfg['data_inicio'], cfg['data_fim'])
             .filter(s.filtro_meses(s.MESES_SECO))
             .map(lambda i: ee.Image.cat([
                 i.select(h.banda_modis(n)).updateMask(i.select(h.qa_modis(n)).eq(0))
                 for n in nums])))
    med = col.median().multiply(h.MODIS_ESC).rename([f'm{n}' for n in nums])
    n_obs = col.select(h.banda_modis(h.CLASSE_MODIS)).count()
    return med.updateMask(n_obs.gte(N_MIN_MODIS))


def para_500(harm: ee.Image, periodo: str, proj_modis) -> ee.Image:
    """Harmonizado (× 1000) → reflectância 0–1 a 500 m, com a cena e a máscara de
    pixel homogêneo de UMA cena só (15 m → 90 m → 500 m, como na 1ª passada)."""
    p90 = harm.select('B01').projection().atScale(h.ESCALA_INTER)
    red = lambda im, r, p: im.reduceResolution(r, maxPixels=64).reproject(p)
    m90 = red(harm.select(bandas_x(periodo)).divide(s.ESCALA_REFL), ee.Reducer.mean(), p90)
    f90 = red(harm.select('B3N').mask(), ee.Reducer.mean(), p90)
    t = harm.select('cena_t')
    tmin = red(red(t, ee.Reducer.min(), p90), ee.Reducer.min(), proj_modis)
    tmax = red(red(t, ee.Reducer.max(), p90), ee.Reducer.max(), proj_modis)
    media = red(m90, ee.Reducer.mean(), proj_modis)
    desvio = red(m90.select(['B02', 'B3N']), ee.Reducer.stdDev(), proj_modis)
    cv = desvio.divide(media.select(['B02', 'B3N'])).reduce(ee.Reducer.max())
    frac = red(f90, ee.Reducer.mean(), proj_modis)
    homog = cv.lt(h.CV_MAX).And(frac.gte(h.FRAC_MIN)).And(tmin.eq(tmax))
    x = media.select([h.AJUSTES[k][0] for k in chaves(periodo)],
                     [f'x_{k}' for k in chaves(periodo)])
    return x.addBands(tmin.toInt().rename('cena_t')).updateMask(homog)


def amostras_celula(periodo: str, nome: str, geom, tipo: str, med, proj_modis):
    harm = ee.Image(f'{h.COL_HARM}/{s.nome_mosaico(periodo, nome)}')
    m = med.select([f'm{h.AJUSTES[k][1]}' for k in chaves(periodo)],
                   [f'm_{k}' for k in chaves(periodo)])
    classe = (med.select(f'm{h.CLASSE_MODIS}').multiply(1 / h.MODIS_ESC)
                 .unitScale(*h.FAIXA_CLASSE).multiply(h.N_CLASSES).floor()
                 .clamp(0, h.N_CLASSES - 1).toInt().rename('classe'))
    img = para_500(harm, periodo, proj_modis).addBands(m).addBands(classe)
    # referência: espalhada pelo brilho; alvo: o mesmo número de pixels por cena
    banda, n = ('classe', N_REF_CLASSE) if tipo == 'ref' else ('cena_t', N_ALVO_CENA)
    return (img.stratifiedSample(numPoints=n, classBand=banda, region=geom,
                                 projection=proj_modis, seed=h.SEMENTE, dropNulls=True,
                                 tileScale=16, geometries=False)
               .map(lambda f: f.set({'tipo': tipo, 'celula': nome})))


def colunas_amostras(periodo: str) -> list[str]:
    return (['tipo', 'celula', 'cena_t', 'classe']
            + [f'x_{k}' for k in chaves(periodo)] + [f'm_{k}' for k in chaves(periodo)])


def etapa1(periodos, alvos, refs, por_nome):
    harm = listar_destino(h.COL_HARM, criar=False)
    proj_modis = ee.ImageCollection(h.MODIS).first().select(h.banda_modis(1)).projection()
    for p in periodos:
        faltam = [c for c in alvos + refs if s.nome_mosaico(p, c) not in harm]
        if faltam:
            log.warning(f'Período {p}: sem harmonizado em {", ".join(faltam)} — período pulado.')
            continue
        med = modis_mediano(p)
        fc = ee.FeatureCollection([
            amostras_celula(p, c, por_nome[c], tipo, med, proj_modis)
            for tipo, lista in (('ref', refs), ('alvo', alvos)) for c in lista]).flatten()
        task = ee.batch.Export.table.toDrive(
            collection=fc, description=f'rel_{p}', folder=PASTA_DRIVE,
            fileNamePrefix=f'rel_{p}', fileFormat='CSV', selectors=colunas_amostras(p))
        task.start()
        log.info(f'  task {task.id} → Drive/{PASTA_DRIVE}/rel_{p}.csv '
                 f'(ref {",".join(refs)}; alvo {",".join(alvos)})')
    print(f'\nQuando terminar, baixe Drive/{PASTA_DRIVE}/rel_<P>.csv para\n  {PASTA_PONTOS}\n'
          'e rode a etapa 2.')


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 2 — RETAS (local)
# ═══════════════════════════════════════════════════════════════════════════
def _ols_robusto(x, y, com_offset: bool) -> dict[str, float]:
    """OLS em duas passadas (corte K_RESIDUO × RMSE), igual à 1ª harmonização."""
    r = h._ols(x, y, com_offset)
    if r['n'] >= 3:
        manter = np.abs(y - (r['g'] * x + r['o'])) <= h.K_RESIDUO * r['rmse']
        r = h._ols(x[manter], y[manter], com_offset)
    return r


def arq_coef(periodo: str):
    return h.SAIDA / f'coef_rel_{periodo}.csv'


def etapa2(periodos):
    for p in periodos:
        arq = PASTA_PONTOS / f'rel_{p}.csv'
        if not arq.exists():
            log.warning(f'Período {p}: {arq} não existe.')
            continue
        pts = pd.read_csv(arq)
        ref, alvo = pts[pts['tipo'] == 'ref'], pts[pts['tipo'] == 'alvo']
        # 1) relação da referência: x = a × m + b (sempre com offset)
        rel = {}
        for k in chaves(p):
            d = ref[[f'm_{k}', f'x_{k}']].dropna().to_numpy(dtype=float)
            rel[k] = _ols_robusto(d[:, 0], d[:, 1], True)
            log.info(f'  {p} ref {k}: x = {rel[k]["g"]:.3f} × MODIS {rel[k]["o"]:+.4f} '
                     f'(n {rel[k]["n"]}, R² {rel[k]["r2"]:.2f})')
        # 2) reta de cada cena alvo: z = ganho × x + offset, z = a × m + b
        linhas = [{'cena_t': -1, 'celulas': 'REF', 'n_pontos': len(ref),
                   **{f'{c}_{k}': rel[k][c] for k in chaves(p) for c in ('n', 'g', 'o', 'r2', 'rmse')}}]
        for ct, g in alvo.groupby('cena_t'):
            linha = {'cena_t': int(ct), 'celulas': ','.join(sorted(g['celula'].unique())),
                     'n_pontos': len(g)}
            for k in chaves(p):
                d = g[[f'x_{k}', f'm_{k}']].dropna().to_numpy(dtype=float)
                z = rel[k]['g'] * d[:, 1] + rel[k]['o']
                r = _ols_robusto(d[:, 0], z, h.AJUSTES[k][3])
                linha.update({f'{c}_{k}': v for c, v in r.items()})
            linhas.append(linha)
        df = pd.DataFrame(linhas)
        df.to_csv(arq_coef(p), index=False)
        cenas, reserva = retas(p)
        log.info(f'Período {p}: {int(cenas["valida"].sum())}/{len(cenas)} cenas alvo com reta '
                 f'válida → {arq_coef(p)}')
        for k, (gm, om) in reserva.items():
            log.info(f'    {k}: ganho mediano {gm:.3f}, offset {om:+.4f}')


def retas(periodo: str):
    """Retas das cenas alvo (sem a linha da referência), coluna `valida` e mediana."""
    df = pd.read_csv(arq_coef(periodo))
    df = df[df['cena_t'] >= 0].copy()
    valida = pd.Series(True, index=df.index)
    for k in chaves(periodo):
        valida &= ((df[f'n_{k}'] >= N_MIN_REL) & (df[f'r2_{k}'] >= R2_MIN_REL)
                   & df[f'g_{k}'].between(*h.GANHO_FAIXA))
    df['valida'] = valida
    if not valida.any():
        raise SystemExit(f'Período {periodo}: nenhuma cena alvo com reta válida.')
    reserva = {k: (statistics.median(df.loc[valida, f'g_{k}']),
                   statistics.median(df.loc[valida, f'o_{k}'])) for k in chaves(periodo)}
    return df, reserva


# ═══════════════════════════════════════════════════════════════════════════
# ETAPA 3 — APLICAR
# ═══════════════════════════════════════════════════════════════════════════
def aplicar(periodo: str, harm: ee.Image, df: pd.DataFrame, reserva: dict) -> ee.Image:
    cfg = s.CONFIG[periodo]
    t = harm.select('cena_t')
    de = [int(v) for v in df['cena_t']]
    usa = t.remap(de, [1 if v else 2 for v in df['valida']], 2).unmask(2)
    esc = h.ESCALA_INT
    corrigidas = []
    for k in chaves(periodo):
        g_res, o_res = reserva[k]
        g = [r[f'g_{k}'] if r['valida'] else g_res for _, r in df.iterrows()]
        o = [r[f'o_{k}'] if r['valida'] else o_res for _, r in df.iterrows()]
        img_g = t.remap(de, [round(v * esc) for v in g], round(g_res * esc)).divide(esc)
        img_o = t.remap(de, [round(v * esc) for v in o], round(o_res * esc)).divide(esc)
        corrigidas.append(harm.select(h.AJUSTES[k][2]).multiply(img_g)
                              .add(img_o.multiply(s.ESCALA_REFL)).round().toInt32())
    saida = (ee.Image.cat(corrigidas).select(cfg['bandas_toa'])
               .addBands(harm.select(s.TIR + ['n_obs', 'cena_t', 'harmonizado']))
               .addBands(usa.updateMask(harm.select('B01').mask()).rename('ajuste_rel').toInt32()))
    return (ee.Image(saida.copyProperties(harm, harm.propertyNames()))
              .set({'ajuste_rel': 'harmonizar_sazonal_relativo.py: reta por cena contra a '
                                  'referência (MCD43A4 mediano seco como ponte)',
                    'ajuste_rel_banda': '1 = reta da cena, 2 = mediana das cenas alvo'}))


def etapa3(periodos, alvos, por_nome):
    harm = listar_destino(h.COL_HARM, criar=False)
    feitos = listar_destino(COL_REL, criar=True)
    fila = []
    for p in periodos:
        if not arq_coef(p).exists():
            log.warning(f'Período {p}: {arq_coef(p)} não existe — rode a etapa 2.')
            continue
        for c in alvos:
            nome = s.nome_mosaico(p, c)
            if nome not in harm:
                log.info(f'  {p} {c}: pulado (sem harmonizado)')
            elif nome in feitos:
                log.info(f'  {p} {c}: pulado (já existe em {COL_REL}; apague para refazer)')
            else:
                fila.append((p, c, nome))
    if not fila or not perguntar_sim_nao(f'Exportar {len(fila)} células ajustadas?', padrao=True):
        return
    tab = {p: retas(p) for p in {p for p, _, _ in fila}}
    for p, c, nome in fila:
        img = aplicar(p, ee.Image(f'{h.COL_HARM}/{nome}'), *tab[p])
        h.exportar_imagem(img, f'rel_{p}_{c}', f'{COL_REL}/{nome}', por_nome[c],
                          {'.default': 'mean', 'cena_t': 'mode', 'harmonizado': 'mode',
                           'ajuste_rel': 'mode', 'n_obs': 'mode'})


# ═══════════════════════════════════════════════════════════════════════════
def perguntar_nomes(texto: str, padrao: list[str], validos: set[str]) -> list[str]:
    while True:
        resp = input(f'{texto} [{",".join(padrao)}]: ').strip().upper()
        nomes = padrao if not resp else [n for n in resp.replace(' ', '').split(',') if n]
        erro = [n for n in nomes if n not in validos]
        if not erro:
            return nomes
        print(f'  → células inexistentes: {", ".join(erro)}')


def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    switch_user(h.CONTA)
    ee.Initialize(project=h.PROJETO)
    log.info(f'Earth Engine com {h.CONTA} (projeto {h.PROJETO})')
    por_nome = {c['nome']: c['geom'] for c in ler_grid().values()}

    print('\n══════ Harmonização relativa (2ª passada) das células que destoam ══════')
    etapa = perguntar_opcao('Etapa:', [
        ('amostras referência + alvo → Drive', 1),
        ('retas a partir do CSV baixado (local)', 2),
        ('aplicar nas células alvo', 3),
    ])
    periodos = perguntar_opcao('Períodos:', [('A e B', ['A', 'B']), ('só A', ['A']), ('só B', ['B'])])
    alvos = perguntar_nomes('Células alvo', ALVOS, set(por_nome)) if etapa != 2 else ALVOS
    if etapa == 1:
        refs = perguntar_nomes('Células de referência', REFERENCIAS, set(por_nome))
        etapa1(periodos, alvos, refs, por_nome)
    elif etapa == 2:
        etapa2(periodos)
    else:
        etapa3(periodos, alvos, por_nome)


if __name__ == '__main__':
    main()
