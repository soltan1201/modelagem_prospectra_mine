# Prospecção mineral sob incerteza — análise dos documentos e plano para Irecê

Base: os 10 textos distintos da pasta `docs/` (os demais são duplicatas):
*Prospeção Mineral*, *Hipóteses Geológicas*, *Inversão Geofísica* (.pdf e .docx),
*Modelos Substitutos na IA*, *Simulação Condicional*, *Risco Mineiro e Incerteza*
(= *Curvas de Incerteza Mineira*), *Gestão de Incerteza Mineira* (2 cópias),
*Inteligência Mineira*, o *Manual de Metodologia* (12 p.) e o *Guia de Estudo* (4 p.).
Todos descrevem a abordagem da KoBold Metals (casos: cobre sedimentar de Mingomba,
Zâmbia, por gravimetria; Ni-Cu-PGE magmático por eletromagnetismo).

---

## 1. Análise lógica do fluxo

O fluxo é um ciclo fechado de **dados → hipóteses → ensemble → decisão → novos
dados**, em que a incerteza é calculada e carregada em todas as etapas, em vez de
ser descartada na passagem de uma etapa para a outra.

```
 (1) Dados             (2) Hipóteses           (3) Modelagem direta     (4) Inversão
     furos  +              geológicas              + emulador IA            probabilística
     geofísica  ──────►    concorrentes  ──────►   (surrogate)   ──────►   (Bayes / MCMC)
     (pontual × amplo)     (priors P(m))           d = g(m) rápido          P(m|d) ∝ P(d|m)P(m)
         ▲                                                                        │
         │                                                                        ▼
 (8) Novo furo         (7) Decisão             (6) Mapas de              (5) Ensemble
     assimilado  ◄──────   sequencial    ◄──────   confiança      ◄──────   N > 1.000
     (atualiza P)          (alvo 50/50,            SNR = E/σ                realizações
                            MCTS/POMDP)                                     (+ simulação condicional)
                                                         │
                                                         ▼
                                              (9) Propagação a jusante:
                                                  curvas teor-tonelagem P10/P50/P90,
                                                  plano de lavra, RQD, água, densidade
```

| Etapa | O que faz | Lógica (por que existe) |
|---|---|---|
| 1. Dados | Furos (química, teor por ICP/XRF, densidade, condutividade, orientação) + levantamentos de superfície (EM, gravimetria) | Furo = muito preciso, mas pontual ("tocar o elefante no escuro"); geofísica = visão ampla, mas desfocada. Nenhum dos dois basta sozinho |
| 2. Hipóteses | Formaliza **todas** as interpretações plausíveis (ex.: condutor = sulfeto de Ni-Cu *ou* grafite; anomalia leve = aquífero espesso *ou* não) como priors | Evita a **ancoragem** cognitiva no primeiro modelo |
| 3. Modelagem direta + surrogate | Calcula o sinal de superfície de cada subsolo candidato; uma rede (CNN/FNO) treinada com milhões de simulações FEM substitui o solucionador (~10.000× mais rápido) | A inversão probabilística precisa de 10⁶–10⁸ avaliações; sem o emulador é inviável |
| 4. Inversão probabilística | Em vez de uma solução suavizada única (Tikhonov), amostra a posteriori P(m\|d); usa superfícies explícitas e covariância por unidade litológica, não voxels livres | O problema inverso não tem solução única; a suavização cria "alucinações" e alvos falsos (Mingomba: 4.000 m estimados × 300 m reais; 5 alvos falsos → 1 real) |
| 5. Ensemble | Milhares de realizações equiprováveis; simulação condicional (SGS / indicadores) com **vários variogramas** (parâmetros do variograma tratados como variáveis aleatórias) | As realizações coincidem nos furos e divergem longe deles — a divergência **é** a incerteza |
| 6. Confiança | SNR(x) = média / desvio-padrão do ensemble | Mostra onde o modelo sabe e onde só "chuta" |
| 7. Decisão | Fura onde P(minério) ≈ 0,5 (máximo ganho de informação ΔH); árvore de decisão sequencial (POMDP + MCTS), olhando N furos à frente | Furar onde já se sabe o resultado (minério certo — "furo de PR" — ou estéril certo) tem ganho de informação nulo |
| 8. Assimilação | O resultado do furo atualiza priors e ensemble | Fecha o ciclo |
| 9. A jusante | O ensemble é propagado (não a média) para curvas de teor e tonelagem com bandas, cortes P10/P50/P90 e plano de lavra; uma ferramenta interna reduz o ensemble a restrições para o software clássico (Surpac/Datamine) | Uma curva única gera a "ilusão de certeza"; o real sempre difere |

**Coerência interna:** o fluxo é consistente porque cada etapa resolve o gargalo
da anterior. A não-unicidade (4) exige hipóteses (2) e ensembles (5); ensembles
exigem velocidade (3); ensembles geram incerteza mensurável (6); a incerteza
mensurável permite decidir pelo ganho de informação (7); e só faz sentido
propagá-la (9) se ela não for jogada fora no meio.

**Pontos frágeis (para levar em conta ao adaptar):**

- A qualidade do resultado depende dos **priors** (hipóteses e variogramas). Se a
  hipótese certa não estiver no conjunto, o ensemble é confiante e errado.
- O surrogate tem erro próprio fora do domínio em que foi treinado; isso também é
  incerteza (epistêmica) e precisa ser medido.
- SNR = E/σ só é interpretável para grandezas positivas e afastadas de zero; para
  probabilidade, o mais útil é o par **(P média, σ)** ou a entropia.
- Os números de impacto (10.000×, ~3 furos por descoberta) são da própria empresa,
  sem validação independente nos textos.

---

## 2. Risco mineiro e incerteza — os métodos descritos

### 2.1 Taxonomia

| Tipo | Definição | Exemplos nos textos | Reduz com… |
|---|---|---|---|
| **Aleatória** (estocástica) | Ruído irredutível de medição e amostragem | XRF portátil × ICP-MS; desvio angular do furo a 1.000 m; ruído EM (correntes telúricas, linhas de alta tensão) | repetição e calibração — não some |
| **Epistêmica** (estrutural) | Falta de conhecimento; várias hipóteses explicam os mesmos dados | camada contínua × cortada por falha cega; condutor de Ni-Cu × grafite/pirrotita | mais dados no lugar certo + testar hipóteses |

### 2.2 Métodos

1. **Paradigma multimodelo** — todas as hipóteses plausíveis viram restrições
   matemáticas e entram juntas na amostragem; nada de "modelo ótimo" único.
2. **Geoestatística bayesiana** — patamar, alcance e pepita do variograma são
   distribuições P(θ), atualizadas com os dados: P(θ|D).
3. **Simulação condicional** (sequencial gaussiana, por indicadores) sobre cada θ
   amostrado → N > 1.000 realizações honrando os furos.
4. **Inversão probabilística com surrogate** — verossimilhança P(d|m) calculada
   pelo emulador; descartam-se as realizações incompatíveis com os dados.
5. **Discriminação por geometria** — o condutor é ajustado como um elipsoide;
   a distribuição da razão de aspecto AR = λmin/λmax é comparada à de depósitos
   conhecidos por divergência de Kullback-Leibler. Por exemplo, AR < 0,05 (laminar)
   indica grafite e o alvo é descartado; AR > 0,15 (compacto) indica sulfeto e o
   alvo é priorizado.
6. **Mapa de confiança SNR** — SNR ≫ 1: confiável; SNR ≤ 1: lacuna de
   informação ou possível "alucinação".
7. **Ganho de informação e alvo 50/50** — ΔH = H(antes) − E[H(depois)]. Quatro
   regras:
   - não furar onde SNR ≫ 3 com minério (furo de PR);
   - não furar onde o estéril já é certo;
   - furar na fronteira em que P ≈ 0,5;
   - escolher o furo que reduz a variância **global**, não só a local.
8. **Decisão sequencial** (POMDP + MCTS) — planeja a campanha N passos à frente.
9. **Curvas de teor e tonelagem com bandas** (±1σ, ±2σ; P10/P50/P90) e
   propagação para RQD, entrada de água, densidade e óxido × sulfeto → **risco
   do projeto** visível para o investimento.

**Risco mineiro**, nos textos, é a incerteza geológica propagada até o resultado
econômico (tonelagem, teor, custo). O argumento central: ao decidir com a média
de um modelo determinístico, a incerteza não deixa de existir; só deixa de ser
vista.

---

## 3. O que muda em Irecê

O fluxo da KoBold é de **subsuperfície** (geofísica e furos, alvos de centenas
de metros a quilômetros). O nosso é de **superfície** (ASTER multiespectral) e
está numa fase anterior: **gerar e priorizar áreas favoráveis** para fosfato.
A lógica de incerteza transfere-se bem; as ferramentas, só em parte.

| Conceito (KoBold) | Equivalente em Irecê |
|---|---|
| Furo (verdade pontual) | ocorrências e depósitos SGB (29 pontos, erro de 250–1.000 m) + visitas de campo, pXRF de solo e rocha |
| Geofísica (visão ampla) | mosaico ASTER (VNIR, SWIR, TIR) + índices espectrais e GLCM; **aerogeofísica do SGB** (gamaespectrometria: K, eTh, **eU**; magnetometria), se houver cobertura |
| Hipóteses geológicas | modelos de fosfato na Fm. Salitre (fosforito estromatolítico em litofácies específicas; enriquecimento supergênico/cárstico; solos residuais), cada um com uma assinatura esperada |
| Inversão / surrogate | o RF que estima o SWIR em 2009–2020 é um "emulador" estatístico; não há inversão física de subsuperfície nesta fase |
| Incerteza aleatória | ruído e saturação do sensor, nuvem residual, TOA × superfície, erro de posição dos pontos, erro do SWIR estimado |
| Incerteza epistêmica | qual mosaico, período e estação usar; normalização; conjunto de variáveis; escolha do fundo no Maxent; representatividade das 29 ocorrências; qual modelo de depósito |
| Ensemble | N rodadas do Maxent/RF variando presenças (bootstrap), fundo, variáveis, hiperparâmetros e realizações do SWIR |
| SNR | mapa de P média e σ do ensemble (e entropia) |
| Alvo 50/50 | locais de **visita de campo** com P ≈ 0,5 e σ alto (máximo aprendizado), não os depósitos já conhecidos |
| Curva teor-tonelagem | **curva área-favorabilidade**: área acima de cada limiar de P, com bandas P10/P50/P90 do ensemble |
| Decisão sequencial | campanhas de campo em rodadas; cada visita vira presença ou **ausência verdadeira** e o ensemble é refeito |

**Por que o eU importa:** fosforitos costumam concentrar urânio, então o canal
eU da gamaespectrometria é, para fosfato, o que a condutividade é para sulfetos.
É o candidato mais forte de dado complementar ao ASTER. É preciso verificar se
há levantamento aerogeofísico do SGB cobrindo o bloco.

---

## 4. Plano de implementação em Irecê

O que já existe e serve de base:

- `aster-mosaico-final`, com as bandas `origem` e `swir_estimado`;
- SWIR estimado por RF com normalização por cena (`estimar_swir_rf.py`);
- `maxent_fosfato.py`, com bandas, índices, GLCM e validação por célula;
- a grid 5 × 5 (processamento por célula).

### Fase 0 — Fechar a base (curto prazo)

- Exportar o `aster-mosaico-final` completo e conferir as camadas por célula.
- Levantar e incorporar dados complementares:
  - aerogeofísica do SGB (K, eTh, **eU**, razões eU/eTh e eU/K; campo magnético
    anômalo);
  - geologia 1:100k ou 1:250k (unidades da Fm. Salitre);
  - SRTM ou Copernicus DEM (relevo cárstico, declividade).
- Revisar os 29 pontos: classe (depósito × ocorrência), erro de posição, unidade
  geológica.

**Entrega:** a pilha de variáveis por célula e uma tabela de pontos revisada.

### Fase 1 — Hipóteses explícitas

- Escrever de 2 a 4 hipóteses de depósito, cada uma com as variáveis e a
  assinatura esperadas. Por exemplo:
  - H1: fosforito na rocha carbonática, com feição de carbonato (B08/B06,
    B13/B14) e eU alto;
  - H2: enriquecimento supergênico em solo e cárstico, com óxidos de ferro,
    argilas e relevo deprimido;
  - H3: controle estrutural, com lineamentos e magnetometria.
- Cada hipótese vira um conjunto de variáveis e uma restrição: por exemplo,
  P = 0 fora da Fm. Salitre na H1.

**Entrega:** uma tabela de hipóteses com as variáveis de cada uma.

### Fase 2 — Ensemble de favorabilidade

Para cada hipótese, N = 100 a 500 realizações do Maxent, variando:

- presenças por bootstrap e jitter de posição dentro do erro de 250–1.000 m;
- fundo, com sementes diferentes e opção de fundo restrito à unidade geológica
  (correção de viés);
- variáveis (sorteio entre as correlacionadas);
- `betaMultiplier` e tipos de função;
- realizações do SWIR, com K florestas de sementes diferentes no
  `estimar_swir_rf.py` → K versões do SWIR e dos índices.

Por célula, salvar a média, o desvio-padrão, os percentis P10/P50/P90 e a
entropia. Combinar as hipóteses pela média ponderada pela validação (AUC por
célula de fora) ou mantê-las separadas.

**Entrega:** a ImageCollection `fosfato_ensemble` (`p_media`, `p_dp`, `p_p10`,
`p_p90`, `entropia`, `hipotese`).

### Fase 3 — Confiança e risco

- Fazer o mapa de confiança: classes alta P com σ baixo (alvo robusto), P ≈ 0,5
  com σ alto (alvo de aprendizado) e baixa P com σ baixo (estéril provável).
- Fazer a curva área-favorabilidade com bandas P10–P90 por hipótese e para o
  bloco.
- Cruzar com `swir_estimado` e `origem`: medir quanto da favorabilidade alta vem
  de SWIR estimado; isso é risco epistêmico explícito.

**Entrega:** mapas e curvas, e uma tabela de risco por célula.

### Fase 4 — Campo em rodadas (decisão sequencial)

1. Selecionar de 10 a 20 locais por rodada pelo **ganho de informação**:
   P ≈ 0,5, σ alto, sem repetir vizinhança e cobrindo hipóteses diferentes. Os
   depósitos já conhecidos (furos de PR) e as áreas de estéril certo ficam de
   fora.
2. No campo: litologia, pXRF de P₂O₅ e U em solo e rocha, fotos e coordenadas
   GNSS. As coordenadas também reduzem o erro de posição dos pontos antigos.
3. Assimilar os resultados: positivos entram como presença e negativos como
   **ausência verdadeira**, o que permite passar do Maxent (presença-fundo) a
   modelos presença-ausência. Recalcular o ensemble e medir quanto o σ global
   caiu.
4. Repetir até o σ estabilizar ou o orçamento acabar.

**Entrega:** um registro de campo por rodada e a evolução do σ médio e da AUC.

### Fase 5 — Priorização final

- Listar as áreas com P10 alto (favoráveis mesmo no cenário conservador),
  apoiadas por SWIR medido e por pelo menos duas hipóteses.
- Fazer a ficha de cada alvo com P10/P50/P90, área, hipóteses que o sustentam,
  dados de campo e lacunas.

**Entrega:** um relatório de alvos para sondagem ou detalhamento.

### Mapeamento nos scripts

| Fase | Script / ação |
|---|---|
| 0 | `build_mosaico_final.py` (exportar); novo `stack_variaveis.py` (aerogeofísica, geologia, MDE) |
| 1 | `docs/hipoteses_fosfato.md` + `HIPOTESES` no Maxent |
| 2 | `maxent_fosfato.py` → modo ensemble (loop de realizações, export por célula); `estimar_swir_rf.py` → K sementes |
| 3 | novo `risco_favorabilidade.py` (curvas e tabela) |
| 4 | novo `selecionar_alvos_campo.py` (ganho de informação, sem vizinhança repetida) |

### Riscos do plano

- **29 pontos agrupados:** o ensemble mede a variabilidade, mas não cria
  informação nova. O campo (Fase 4) é o que reduz a incerteza epistêmica de
  verdade.
- **Mancha de carbonato:** o mapa pode refletir "carbonato da Fm. Salitre", e
  não "fosfato". Sem o eU ou um dado de campo, H1 e "calcário estéril" ficam
  indistinguíveis — o mesmo dilema de grafite × sulfeto dos textos.
- **SWIR estimado:** nas células sem SWIR medido, a favorabilidade herda o erro
  do RF; isso precisa aparecer no σ, não ficar escondido.
