
🌟 Mosaico Multi-Temporal ASTER - Quality Score
Script Google Earth Engine para criação de mosaico multi-temporal de alta qualidade utilizando imagens ASTER L1T, com seleção inteligente dos melhores pixels disponíveis em uma série temporal de mosaicos semestrais.
📋 Índice
    Visão Geral
    Funcionalidades
    Requisitos
    Como Usar
    Estrutura do Script
    Bandas de Saída
    Quality Score - Critérios
    Exportação
    Visualização
    Notas Técnicas

🌍 Visão Geral
Este script processa uma coleção de mosaicos ASTER semestrais (2002-2008) e gera um mosaico final otimizado selecionando pixel a pixel a melhor observação disponível, eliminando automaticamente:
    ☁️ Nuvens residuais
    🌑 Sombras de nuvens
    🔲 Bordas de nuvem (artefatos de alta textura)
    📊 Pixels saturados ou inválidos

O resultado é uma imagem contínua de alta qualidade para mapeamento geológico, prospecção mineral e análise ambiental.
⚡ Funcionalidades
🎯 Quality Score Multi-Critério

Sistema de pontuação que avalia cada pixel baseado em 5 critérios:
Critério	Peso	Descrição
Brilho (NIR)	30%	Penaliza pixels escuros (sombras)
Textura (GLCM)	25%	Penaliza bordas de nuvem (alto contraste)
NDVI	20%	Penaliza valores extremos (nuvens/vegetação densa)
SWIR Válido	15%	Garante valores dentro do range esperado
Variância Local	10%	Penaliza pixels isolados (ruído)

🔍 Detecção de Pixels Ruins
Máscara automática para remover:
    Sombras residuais (NIR < 0.06)
    Nuvens residuais (Verde > 0.35 + SWIR > 0.4)
    Bordas de nuvem (Contraste GLCM > 180)
    Pixels saturados

📊 Banda de Frequência
Cria uma banda adicional que indica quantos mosaicos semestrais contribuíram para cada pixel final:
    Valores altos: Área com muitas observações (alta confiabilidade)
    Valores baixos: Área com poucas observações (menor confiabilidade)

🗺️ Índices Geológicos
Calcula automaticamente índices para análise mineral:
Índice	Fórmula	Aplicação
AlOH_Clay	(B05 + B07) / (2 × B06)	Detecção de Argilas (Caulinita)
Carbonate	B08 / B06	Detecção de Carbonatos
Ferric_Fe	B02 / B01	Óxidos de Ferro (Hematita/Goethita)
NDVI	(B3N - B02) / (B3N + B02)	Índice de Vegetação
SAVI	1.5 × (B3N - B02) / (B3N + B02 + 0.5)	Vegetação (corrigido para solo)

💾 Exportação Otimizada (16-bit)
Conversão automática para INT16 reduzindo em 50% o espaço de armazenamento:
Tipo de Dado	Range	Escala	Formato
Reflectância	0-1.5	× 10000	INT16
Índices	-1 a 2	× 10000	INT16
Frequência	1-N	× 1	INT16
Quality Score	0-1	× 10000	INT16

📈 Análise de Contribuição
Gera relatório automático com:
    Estatísticas do Quality Score (média, desvio padrão)
    Contribuição de cada mosaico semestral
    Número de pixels válidos por período
    🔧 Preenchimento de Gaps

Interpolação espacial opcional para preencher pixels sem dados usando filtro Gaussiano.
📦 Requisitos
Assets Necessários
javascript
    // Coleção de mosaicos semestrais ASTER (2002-2008)
    var id_asset_mosaicos = 'projects/seu-projeto/assets/mosaic_aster';

Bandas Esperadas
Cada mosaico semestral deve conter as bandas:
    Banda	Descrição	Resolução
    B01	Verde (556 nm)	15m
    B02	Vermelho (661 nm)	15m
    B3N	NIR (807 nm)	15m
    B04	SWIR 1 (1656 nm)	30m
    B05	SWIR 2 (2167 nm)	30m
    B06	SWIR 3 (2209 nm)	30m
    B07	SWIR 4 (2263 nm)	30m
    B08	SWIR 5 (2336 nm)	30m
    B09	SWIR 6 (2395 nm)	30m
🚀 Como Usar
1. Configuração Inicial
javascript

    // Definir área de estudo (coordenadas do retângulo)
    var list_coord = [
        [-44.46, -14.79],  // NW
        [-39.20, -14.79],  // NE
        [-39.20, -9.41],   // SE
        [-44.46, -9.41]    // SW
    ];

    // ID da coleção de mosaicos
    var id_asset_mosaicos = 'projects/seu-projeto/assets/mosaic_aster';

2. Executar Processamento
javascript
    // O script automaticamente:
    // 1. Carrega todos os mosaicos semestrais
    // 2. Calcula Quality Score para cada pixel
    // 3. Seleciona os melhores pixels
    // 4. Gera índices geológicos
    // 5. Exporta resultado

3. Visualizar Resultados

O script adiciona automaticamente as seguintes camadas ao mapa:

    🌟 MOSAICO FINAL - RGB Natural: Composição colorida verdadeira

    🌟 MOSAICO FINAL - SWIR Geologia: Falsa-cor para geologia

    🌟 Al-OH (Argilas): Índice de argilominerais

    📊 Frequência: Número de observações por pixel

    📈 Quality Score Final: Pontuação de qualidade

4. Exportar Resultados
javascript

    // Asset (para uso no GEE)
    Export.image.toAsset({
        image: mosaico_16bit,
        description: 'MOSAICO_FINAL_ASTER_2002_2008',
        assetId: 'MOSAICO_FINAL_ASTER_2002_2008',
        region: area_estudo,
        scale: 15,
        maxPixels: 1e13
    });

    // Metadados (CSV)
    Export.table.toDrive({
        collection: metadados,
        description: 'METADADOS_MOSAICO_FINAL_2002_2008',
        folder: 'ASTER_Metadados'
    });

🏗️ Estrutura do Script
text

📁 Script Principal
├── 📐 1. Parâmetros e Configuração
│   ├── Área de estudo
│   └── Bandas para processamento
│
├── 📂 2. Carregar Mosaicos Semestrais
│   └── ImageCollection do Asset
│
├── 🧮 3. Funções de Processamento
│   ├── calcularQualityScore()      ← Score multi-critério
│   ├── detectarPixelsRuins()       ← Máscara de artefatos
│   └── adicionarIndicesGeo()       ← Índices minerais
│
├── 🔄 4. Processamento da Coleção
│   ├── Aplicar Quality Score
│   ├── Detectar pixels ruins
│   └── Adicionar timestamp
│
├── 🎯 5. Quality Mosaic Multi-Temporal
│   └── qualityMosaic('quality_score')
│
├── 📊 6. Pós-Processamento
│   ├── Preencher gaps (opcional)
│   ├── Adicionar índices
│   └── Banda de frequência
│
├── 👁️ 7. Visualização
│   ├── RGB Natural
│   ├── SWIR Geologia
│   ├── Índices minerais
│   └── Métricas de qualidade
│
└── 💾 8. Exportação
    ├── Converter para 16-bit
    ├── Exportar para Asset
    └── Exportar metadados

📊 Bandas de Saída
    Bandas Espectrais (Reflectância TOA)
    Banda	Nome	Range	Escala
    B01	AST_Green_556nm	0-1.5	×10000
    B02	AST_Red_661nm	0-1.5	×10000
    B3N	AST_NIR_807nm	0-1.5	×10000
    B04	AST_SWIR1_1656nm	0-1.5	×10000
    B05	AST_SWIR2_2167nm_AlOH	0-1.5	×10000
    B06	AST_SWIR3_2209nm_Ref	0-1.5	×10000
    B07	AST_SWIR4_2263nm_AlMgOH	0-1.5	×10000
    B08	AST_SWIR5_2336nm_CO3	0-1.5	×10000
    B09	AST_SWIR6_2395nm	0-1.5	×10000
    Índices Calculados
    Banda	Descrição	Range Típico
    NDVI	Vegetation Index	-1 a +1
    AlOH_Clay	Índice de Argilas	0.8-1.4
    Carbonate	Índice de Carbonatos	0.5-1.5
    Ferric_Fe	Óxidos de Ferro	0.5-2.0
    SAVI	Soil Adjusted VI	-1 a +1
    Bandas de Qualidade
    Banda	Descrição	Range
    quality_score	Pontuação do pixel (0-1)	0-10000
    frequency	Nº de observações	1-N
    timestamp	Data da observação (ms)	-

📐 Quality Score - Critérios
Fórmula Completa
text
Quality = (Brilho × 0.30) + (Textura × 0.25) + (NDVI × 0.20) + (SWIR × 0.15) + (Variância × 0.10)

Interpretação dos Valores
Score	Classificação	Significado
    0.8 - 1.0	⭐⭐⭐⭐⭐ Excelente	Sem nuvens/sombras, textura ideal
    0.5 - 0.8	⭐⭐⭐⭐ Bom	Pequenos artefatos, ainda utilizável
    0.2 - 0.5	⭐⭐⭐ Regular	Sombras leves ou bordas de nuvem
    0.0 - 0.2	⭐⭐ Ruim	Nuvens/sombras fortes (descartado)

💾 Exportação
Para Google Earth Engine Asset
javascript

    // Uso posterior no GEE
    var mosaico = ee.Image('projects/seu-projeto/assets/MOSAICO_FINAL_ASTER_2002_2008');

    // Converter de volta para float
    var reflectancia = mosaico.select(['B01','B02','B3N','B04','B05','B06','B07','B08','B09'])
        .divide(10000).float();

    var indices = mosaico.select(['NDVI','AlOH_Clay','Carbonate','Ferric_Fe','SAVI'])
        .divide(10000).float();

    var quality = mosaico.select('quality_score').divide(10000).float();

Para Google Drive
    Formato: GeoTIFF (Asset) ou CSV (Metadados)
    Resolução: 15 metros
    Projeção: Mantida do mosaico original
    Compressão: Automática (GEE)

👁️ Visualização
Parâmetros Recomendados
    Composição	Bandas	Min	Max	Gamma
    RGB Natural	B3N, B02, B01	0.05	0.35	1.4
    SWIR Geologia	B04, B06, B3N	0.10	0.50	1.2
    Al-OH Clay	AlOH_Clay	0.80	1.40	1.0
    Frequência	frequency	1	N	1.0
    Quality Score	quality_score	0	1	1.0
    Paletas de Cor
    javascript

    // Argilas
    palette: ['blue', 'cyan', 'green', 'yellow', 'red']

    // Frequência
    palette: ['red', 'orange', 'yellow', 'green', 'darkgreen']

    // Quality Score
    palette: ['red', 'yellow', 'green']

    // NDVI
    palette: ['brown', 'yellow', 'green']

🔬 Notas Técnicas
Performance
    Tempo de processamento: ~10-30 minutos (depende da área)
    Memória: Otimizado para áreas grandes (até 100.000 km²)
    Escala: 15m (VNIR) / 30m (SWIR upsampled)

Limitações Conhecidas
    SWIR pós-2008: ASTER SWIR parou de funcionar em 2008
    Geometria de iluminação: Sombras podem persistir em áreas montanhosas
    Nuvens persistentes: Áreas com cobertura de nuvens >90% podem ter gaps

Melhorias Futuras
    Integração com Landsat para preencher gaps SWIR
    Correção atmosférica DOS (Dark Object Subtraction)
    BRDF correction para normalizar ângulos de visada
    Machine Learning para detecção de nuvens

📝 Referências
    ASTER L1T: NASA LP DAAC
    Índice Al-OH: Rowan & Mars (2003) - Lithologic mapping in the Mountain Pass, California area using ASTER data
    Quality Mosaic: Google Earth Engine Documentação

👥 Contribuição
Script desenvolvido para processamento de dados ASTER com foco em mapeamento geológico e prospecção mineral.

Projeto: Prospectra 4.0
Autores: [Seu Nome]
Data: 2026
Versão: 1.0.0
📄 Licença
Este script é distribuído sob licença MIT. Consulte o arquivo LICENSE para mais detalhes.