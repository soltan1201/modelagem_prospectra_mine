# Revisão — Normalização de histograma dos mosaicos semestrais

Documentos revisados:

- [`NormalizacaoHistograma.pdf`](NormalizacaoHistograma.pdf) — fórmula, fluxo de trabalho e resultados por semestre.
- [`normalizacao_mosaicos.js`](normalizacao_mosaicos.js) — script de laboratório no GEE.

**Resumo:** a fórmula central está correta nos dois documentos. O script implementa uma
versão diferente (robusta, com mediana e percentis), também matematicamente válida, mas
que não está descrita no PDF. Os problemas sérios não estão na conta, e sim na escolha do
polígono âncora e nas "provas" usadas para validar a correção.

---

## 1. As equações

### 1.1 Fórmula do PDF — correta

$$
Pixel_{novo} = \frac{Pixel_{alvo} - \mu_{alvo}}{\sigma_{alvo}} \times \sigma_{ref} + \mu_{ref}
$$

É o ajuste clássico de média e desvio padrão (*mean/std matching*): o alvo passa a ter a
mesma média e o mesmo desvio padrão da referência.

### 1.2 Erros no texto do PDF (página 1)

| Trecho no PDF | Problema | Correção |
|---|---|---|
| "Ômega sendo o contraste geral da imagem" | O símbolo é **σ (sigma)**, o desvio padrão. ω (ômega) é outra letra. | "σ (sigma) é o desvio padrão, o contraste geral da imagem" |
| "Xnovo = z x ômegaNovo = miNova" | Há um "=" no lugar do "+". | **Xnovo = z × σ_ref + μ_ref** |
| "para que a variação máxima não ultrapasse 1" | Depois do Z-score, o que fica igual a 1 é o **desvio padrão**, não o valor máximo. Valores entre ±2 e ±3 continuam normais. | "para que o desvio padrão passe a ser 1" |

### 1.3 O script usa uma fórmula diferente da do PDF

Em [`normalizacao_mosaicos.js`](normalizacao_mosaicos.js) (linhas 59–81):

$$
Pixel_{novo} = \frac{Pixel_{alvo} - P50_{alvo}}{\max(\hat\sigma_{alvo},\ \varepsilon)} \times \hat\sigma_{ref} + P50_{ref},
\qquad \hat\sigma = \frac{P95 - P5}{3{,}29}
$$

- **Mediana (P50) no lugar da média** — correto, e menos sensível a nuvem residual dentro do polígono.
- **Desvio padrão estimado por (P95 − P5) / 3,29** — correto: numa distribuição normal,
  P95 − P5 = 2 × 1,645σ = 3,29σ.
- **Epsilon (1e-6) contra divisão por zero** — correto.

**Ação:** atualizar o PDF para descrever a versão com mediana e percentis, que é a que o
script realmente usa.

---

## 2. Problemas de método (os mais sérios)

### 2.1 O polígono âncora precisa ser de um alvo que não muda entre as datas

A normalização só corrige diferenças de atmosfera e iluminação se o chão dentro do
polígono for **o mesmo material nas duas datas** (*pseudo-invariant features*: rocha,
pedreira, água funda).

No semestre 2001.1 (páginas 3 e 4 do PDF), o polígono é descrito como *"solo mais exposto
com menos vegetação"* no mosaico geral e *"área com muita vegetação"* no semestre alvo. O
próprio script orienta: *"Desenhe a roi_ancora estritamente sobre a mancha vegetal (para
secar)"*.

Nesse caso a correção não remove efeito de atmosfera: ela transforma estatisticamente
**vegetação em solo** e aplica esse mesmo ganho na imagem inteira. Todos os outros alvos
(água, rocha, outros solos) ficam distorcidos. O "vermelho" de um semestre úmido, na
composição falsa-cor, é vegetação de verdade, não erro do sensor.

### 2.2 Um único ganho e deslocamento para ~340 mil km²

O mosaico semestral junta dezenas de cenas ASTER, cada uma com a sua data e iluminação. O
ajuste é calculado num polígono de poucos km² e por isso só representa a cena que cobre
esse polígono. Nas demais cenas o erro radiométrico é diferente, e a mesma reta não o
corrige.

### 2.3 As "provas" do PDF não provam a correção

- **Prova 1 (histograma):** o ajuste foi calculado nesse mesmo polígono, então o
  histograma normalizado (verde) ficar em cima do de referência (azul) é garantido por
  construção. Para validar, é preciso comparar em **outros polígonos** de alvos estáveis
  que não foram usados no ajuste.
- **Prova 2 (dispersão com r² = 1):** a saída é uma reta da entrada
  (normalizado = a × original + b), então r² = 1 sempre. Isso mostra apenas que a
  transformação é linear, não que o espectro ficou íntegro.

### 2.4 Os índices minerais não se preservam

Cada banda recebe o seu próprio ganho e o seu próprio deslocamento, então as **razões
entre bandas mudam de valor**. Os índices do projeto (Al-OH, carbonato, ferro férrico)
são razões de banda, então são afetados diretamente.

Os ganhos medidos nos próprios gráficos do PDF (banda B3N) mostram o tamanho do efeito:

| Semestre | Reta da Prova 2 | Ganho |
|---|---|---|
| 2001.1 | 0,776 × original + 0,0043 | 0,776 |
| 2004.1 | 0,383 × original + 0,05 | **0,383** (contraste comprimido 2,6×) |
| 2006.1 | 1,078 × original − 0,081 | 1,078 |
| 2007.1 | 0,989 × original − 0,049 | 0,989 |

Mudanças desse tamanho indicam diferença de cobertura do solo, não de atmosfera.

### 2.5 A referência é circular

- O `mosaico_mestre` do script é montado com `qualityMosaic` de **todos** os semestres de
  `mosaic_aster`, inclusive o próprio semestre alvo.
- O passo 5 do fluxo do PDF gera a "Versão 2.0" do mosaico geral a partir dos semestres
  corrigidos com base na versão 1.0.

---

## 3. Bugs no script

1. **Todas as bandas são normalizadas**, inclusive `quality`, `SAVI`, `NDWI` e o TIR
   (B10–B14), e `divide(10000)` é aplicado a todas. Após a correção do
   `converterPara16Bit` (TIR passa a ter escala ×10, não ×10000), dividir o TIR por 10000
   fica errado. A normalização deveria valer só para **B01–B09**.
2. **O filtro de "bandas válidas" provavelmente não filtra nada** (linhas 32–44). O
   `reduceRegion` com `Reducer.first()` devolve a chave da banda mesmo quando o valor é
   `null`. Se uma banda não tiver dado dentro do polígono, `ee.Number(null)` dá erro mais
   adiante. (Não testado.)
3. **Caminho errado no PDF:** a entrada do semestre 2007.2 aponta para
   `ASTER_QualityMosaic_2007_semestre_1_INT16`; deveria ser `..._2007_semestre_2_INT16`.

---

## 4. Recomendações

1. Usar **apenas alvos estáveis** como âncora: vários polígonos de rocha, pedreira ou água
   funda, espalhados pela área de estudo.
2. **Validar em polígonos independentes**, que não entraram no cálculo do ajuste.
3. Normalizar **somente B01–B09**, cada uma na sua escala correta.
4. Documentar no PDF a versão com mediana e percentis.
5. Se o objetivo é remover o "vermelho" (vegetação) dos semestres úmidos para a análise
   geológica, normalização de histograma não é a ferramenta adequada. O mais indicado é
   **preferir pixels de época seca no `qualityMosaic`**, por exemplo dando peso a NDVI
   baixo na banda de `quality`.
