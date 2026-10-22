# Ranking Atual: Top 10 Modelos CNN para Classificação de Imagens

O benchmark padrão-ouro para avaliação de modelos de classificação de imagens continua sendo o **ImageNet-1K** (ILSVRC 2012). Embora os *Vision Transformers* (ViTs) ocupem frequentemente o topo absoluto dos leaderboards gerais, as arquiteturas **CNN (Redes Neurais Convolucionais)** e híbridas evoluíram drasticamente. Elas mantêm desempenho de estado da arte (SOTA), muitas vezes com melhor custo-benefício computacional, menor latência e maior facilidade de implantação em produção.

Abaixo está o ranking dos 10 melhores modelos baseados em CNN (ou híbridos CNN-dominantes), ordenados por sua **acurácia Top-1** no ImageNet (considerando pré-treino em datasets amplos como ImageNet-21K ou JFT).

---

## 🏆 Resumo do Ranking

| # | Modelo | Tipo | Acurácia Top-1 (ImageNet) | Foco Principal |
|:--|:---|:---|:---|:---|
| **1** | **CoAtNet** | Híbrida (CNN + Attention) | ~88.5% – 90.8% | Escalabilidade e fusão de indutividade local com contexto global. |
| **2** | **ConvNeXt V2** | CNN Pura | ~88.9% | Modernização extrema com pré-treino auto-supervisionado (MAE). |
| **3** | **EfficientNetV2** | CNN Pura | ~87.3% | Equilíbrio ideal entre acurácia, parâmetros e velocidade de treino. |
| **4** | **NFNet** | CNN Pura | ~86.0% | Remoção total de BatchNorm para treinamento distribuído estável. |
| **5** | **ConvNeXt (v1)** | CNN Pura | ~85.5% – 86.0% | "ResNet modernizado" com design inspirado em Vision Transformers. |
| **6** | **MobileNetV4** | CNN Pura / Híbrida Leve | >86.0% (variantes L) | Otimização agressiva para latência em CPUs, GPUs e NPUs móveis. |
| **7** | **ResNet-RS** | CNN Pura | ~84.4% | Re-escalonamento sistemático de resolução e regularização do ResNet. |
| **8** | **RegNetY** | CNN Pura | ~84.4% | Espaço de busca de arquitetura neural (NAS) com largura quantizada. |
| **9** | **RepVGG** | CNN Pura | ~83.0% – 84.0% | Reparametrização estrutural: multi-ramos no treino, único ramo na inferência. |
| **10**| **ResNet-152/200**| CNN Pura | ~83.0% – 84.0% | Arquitetura residual clássica e profunda, padrão da indústria. |

---

## 🧠 Descrição Detalhada das Arquiteturas

### 1. CoAtNet (Convolution and Attention Network)
- **Arquitetura**: Híbrida. Combina blocos **MBConv** (convolução depthwise separável com expansão, herdada do EfficientNet) nos estágios iniciais da rede, e blocos de **Multi-Head Self-Attention (MHSA)** nos estágios finais.
- **Inovação**: Aproveita o *inductive bias* e a eficiência de memória das convoluções para extração de características locais em baixas resoluções, e a capacidade de modelagem de dependências de longo alcance dos Transformers em altas resoluções. Utiliza *relative attention* para melhor generalização.

### 2. ConvNeXt V2
- **Arquitetura**: CNN Pura. Evolução do ConvNeXt v1, mantendo a macro-arquitetura em estágios (stem de patchificação, bottlenecks invertidos, kernels 7x7 depthwise), mas introduzindo a **Global Response Normalization (GRN)**.
- **Inovação**: A GRN é aplicada após a convolução depthwise para evitar o "colapso de características" (feature collapse) e aumentar a diversidade de canais sem a sobrecarga computacional de normalizações complexas. É projetada para ser pré-treinada com **Masked Autoencoders (MAE)**, o que eleva drasticamente sua acurácia com menos dados supervisionados.

### 3. EfficientNetV2
- **Arquitetura**: CNN Pura. Introduz o bloco **Fused-MBConv** nos estágios iniciais da rede.
- **Inovação**: O Fused-MBConv funde a convolução de expansão 1x1 e a convolução depthwise 3x3 em uma única convolução padrão 3x3. Isso reduz a fragmentação de memória e acelera significativamente o treinamento em aceleradores modernos (GPUs/TPUs). Também emprega **Aprendizado Progressivo**, ajustando dinamicamente o tamanho da imagem e a força da regularização (como dropout e augmentation) ao longo do treinamento.

### 4. NFNet (Normalizer-Free Networks)
- **Arquitetura**: CNN Pura baseada em blocos residuais, mas **sem nenhuma camada de Batch Normalization (BatchNorm)**.
- **Inovação**: A BatchNorm é um gargalo em treinamento distribuído (devido à sincronização entre GPUs). A NFNet a substitui por **Scaled Weight Standardization (SWS)** nas convoluções e utiliza **Adaptive Gradient Clipping (AGC)** para estabilizar a norma dos gradientes. Isso permite o uso de tamanhos de lote (batch sizes) massivos com utilização de hardware próxima de 100%.

### 5. ConvNeXt (v1)
- **Arquitetura**: CNN Pura. Um "ResNet modernizado" que adota deliberadamente escolhas de design dos Vision Transformers.
- **Inovação**: Substitui o *stem* inicial por uma convolução 4x4 com stride 4 (patchify); usa bottlenecks invertidos (expansão 4x); substitui BatchNorm por **LayerNorm**; utiliza apenas uma ativação GELU por bloco; e emprega convoluções depthwise com **kernel grande (7x7)** para aumentar o campo receptivo, imitando a janela de atenção dos ViTs.

### 6. MobileNetV4
- **Arquitetura**: CNN Pura (com variantes híbridas opcionais). Introduz o bloco **Universal Inverted Bottleneck (UIB)**.
- **Inovação**: O UIB permite flexibilidade no tamanho do kernel (de 3x3 até 7x7 ou kernels assimétricos) e na expansão de canais, adaptando-se dinamicamente ao hardware alvo. A arquitetura foi co-projetada com ferramentas de NAS (Neural Architecture Search) específicas para backends móveis (CPU ARM, GPU móvel, DSP/NPU), otimizando operações de memória e evitando gargalos de *memory-bound*.

### 7. ResNet-RS (ResNet Re-Scaling)
- **Arquitetura**: CNN Pura baseada no bloco *bottleneck* residual clássico (1x1 → 3x3 → 1x1).
- **Inovação**: Não é uma nova topologia, mas um conjunto de regras de re-escalonamento descobertas empiricamente. As principais mudanças são: (1) aumentar a resolução da imagem de entrada mais rapidamente do que a profundidade da rede; (2) ajustar os *strides* para manter mapas de características de alta resolução por mais tempo; (3) aplicar regularização moderna (RandAugment, Dropout) de forma calibrada para o novo tamanho.

### 8. RegNetY
- **Arquitetura**: CNN Pura derivada de Neural Architecture Search (NAS).
- **Inovação**: Em vez de buscar uma rede única, os autores definiram um "espaço de design" parametrizado. O RegNetY segue uma regra simples: a largura da rede (número de canais) aumenta linearmente (ou via função quantizada) a cada estágio, e **todos os blocos dentro de um mesmo estágio têm exatamente a mesma largura e razão de bottleneck**. Utiliza convoluções em grupo (grouped convolutions) para eficiência.

### 9. RepVGG
- **Arquitetura**: CNN Pura com **Reparametrização Estrutural**.
- **Inovação**: Durante o **treinamento**, a rede possui uma topologia multi-ramo (semelhante ao ResNet), com três ramos paralelos: uma convolução 3x3, uma convolução 1x1 e uma conexão de identidade. Isso facilita o fluxo de gradiente e a capacidade de representação. Durante a **inferência**, esses três ramos são fundidos matematicamente em uma única camada de convolução 3x3 equivalente. O resultado é uma arquitetura em linha reta (estilo VGG) com velocidade de inferência máxima e zero sobrecarga de memória.

### 10. ResNet-152 / ResNet-200
- **Arquitetura**: CNN Pura profunda baseada em blocos residuais *bottleneck*.
- **Inovação**: A arquitetura clássica que resolveu o problema do vanishing gradient em redes muito profundas através de **skip connections** (adição da entrada à saída do bloco). O design *bottleneck* (1x1 para reduzir dimensão, 3x3 para processamento espacial, 1x1 para restaurar dimensão) permite que redes com 152 ou 200 camadas sejam treinadas de forma estável, servindo como *backbone* confiável para inúmeras tarefas de visão computacional.

---

## 📌 Observações Importantes

- **Top-1 vs. Top-5**: A "Acurácia Top-1" exige que a classe com a maior probabilidade prevista seja exatamente a correta. A "Acurácia Top-5" (onde a resposta correta está entre as 5 previsões mais prováveis) para esses modelos de ponta geralmente supera **98%**.
- **Impacto do Pré-treino**: As acurácias mais altas (>87%) quase sempre exigem pré-treino em datasets maiores (como ImageNet-21K ou JFT-300M), seguido de *fine-tuning* no ImageNet-1K. Treinar do zero apenas no ImageNet-1K resulta em números ligeiramente menores.
- **Híbridos vs. Puros**: Modelos como o CoAtNet são tecnicamente híbridos, mas são incluídos neste ranking porque a base convolucional é fundamental para seu desempenho e eficiência em comparação a ViTs puros.

---

## 💡 Recomendações Práticas de Escolha

| Cenário de Uso | Modelo Recomendado | Motivo |
|:---|:---|:---|
| **Máxima Acurácia** (sem restrições de hardware) | **CoAtNet** ou **ConvNeXt V2 Huge** | Oferecem o teto de performance mais alto para arquiteturas baseadas em convolução. |
| **Melhor Equilíbrio** (Acurácia vs. Custo Computacional) | **EfficientNetV2** (M ou L) ou **ConvNeXt Base** | Escolhas mais sensatas para a maioria dos projetos de produção, oferecendo excelente desempenho sem custo proibitivo. |
| **Dispositivos Móveis, IoT ou Tempo Real** | **MobileNetV4** ou **EfficientNetV2-S** | Otimizados para latência ultrabaixa e baixo consumo de memória, rodando eficientemente em CPUs ou NPUs de borda. |
| **Sistemas Legados ou Compatibilidade Total** | **ResNet-50** ou **ResNet-152** | A escolha mais segura, com documentação vasta, suporte nativo em ONNX, TensorRT, OpenVINO e todas as bibliotecas de ML. |

> **Nota:** Se você tiver um cenário específico (como restrição severa de memória, necessidade de rodar exclusivamente em CPU, ou um dataset de domínio específico como imagens médicas/satélite), a ordem de preferência pode mudar. 