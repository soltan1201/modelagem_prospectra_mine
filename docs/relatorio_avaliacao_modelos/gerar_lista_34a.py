"""
Gera lista_34a.tex (lista explicada das variáveis, célula [34a] do
avaliar_maxent.ipynb) a partir do CSV ranking_features_shap_explicado_<data>.csv.
O relatório inclui o arquivo com \\input{lista_34a}.

Uso: python gerar_lista_34a.py [caminho_do_csv]
     (sem argumento: o CSV mais recente em src/modeling/models/)
"""

import sys
from pathlib import Path

import pandas as pd

AQUI = Path(__file__).resolve().parent
MODELS = AQUI.parents[1] / 'src' / 'modeling' / 'models'
SAIDA = AQUI / 'lista_34a.tex'

TROCAS = {'\\': r'\textbackslash{}', '_': r'\_', '%': r'\%', '&': r'\&', '#': r'\#',
          '✔': r'\ding{51}', '✕': r'\ding{55}', '↑': r'$\uparrow$', '↓': r'$\downarrow$',
          '×': r'$\times$', '→': r'$\rightarrow$', '³⁺': r'$^{3+}$', '~': r'$\sim$',
          'µ': r'\textmu{}', '—': '---', '–': '--'}


def tex(t: str) -> str:
    t = str(t)
    for a, b in TROCAS.items():
        t = t.replace(a, b)
    return t


def main():
    csv = Path(sys.argv[1]) if len(sys.argv) > 1 else max(MODELS.glob('ranking_features_shap_explicado_*.csv'))
    df = pd.read_csv(csv)
    linhas = [f'% gerado por gerar_lista_34a.py a partir de {csv.name} — não editar à mão',
              r'\begin{itemize}\small']
    for r in df.itertuples():
        marca = r'\ding{51}' if r.selecionada else r'\ding{55}'
        linhas.append(f'  \\item[{r.posto}.] {marca} \\texttt{{{tex(r.feature)}}} '
                      f'[{tex(r.tipo)}, posto médio {f"{r.posto_medio + 1e-9:.1f}".replace(".", ",")}]\\\\\n    {tex(r.porque)}')
    linhas.append(r'\end{itemize}')
    SAIDA.write_text('\n'.join(linhas) + '\n')
    print(f'{len(df)} variáveis → {SAIDA}')


if __name__ == '__main__':
    main()
