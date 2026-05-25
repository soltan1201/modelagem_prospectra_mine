# ============================================================
# 1. Definir listas de palavras-chave por categoria
# ============================================================
lstwordMine = [
    'mineral', 'ore', 'deposit', 'prospect', 'exploration', 'mining', 'resource',
    'geochemical', 'geophysical', 'alteration', 'mineralization', 'prospecting'
]

lstwordDL = [
    'deep learning', 'neural network', 'cnn', 'rnn', 'lstm', 'transformer',
    'autoencoder', 'gan', 'generative adversarial', 'machine learning',
    'convolutional', 'deep neural', 'resnet', 'mamba', 'attention'
]

lstwordRS = [
    'remote sensing', 'hyperspectral', 'multispectral', 'satellite', 'uav',
    'drone', 'landsat', 'sentinel', 'aster', 'prisma', 'enmap', 'aviris',
    'imaging spectroscopy', 'sar', 'radar', 'lidar', 'spectral', 'earth observation'
]

# ============================================================
# 2. Função para contar ocorrências de uma lista em um texto
# ============================================================
def count_keywords(text, keyword_list):
    if not isinstance(text, str):
        return 0
    text_lower = text.lower()
    # Remove pontuação simples (opcional, mas evita 'remote-sensing' não ser encontrado)
    text_clean = re.sub(r'[^\w\s]', ' ', text_lower)
    words = set(text_clean.split())
    count = 0
    for kw in keyword_list:
        if kw in text_clean:  # substring matching, captura variações como 'mineralization'
            count += 1
    return count

# ============================================================
# 3. Aplicar contagem em cada linha do DataFrame
# ============================================================
import re

df_scopus['Mine'] = df_scopus.apply(
    lambda row: count_keywords(row['Title'], lstwordMine) +
                count_keywords(row['Abstract'], lstwordMine) +
                count_keywords(row['Author Keywords'], lstwordMine) +
                count_keywords(row['Index Keywords'], lstwordMine),
    axis=1
)

df_scopus['DeepLearning'] = df_scopus.apply(
    lambda row: count_keywords(row['Title'], lstwordDL) +
                count_keywords(row['Abstract'], lstwordDL) +
                count_keywords(row['Author Keywords'], lstwordDL) +
                count_keywords(row['Index Keywords'], lstwordDL),
    axis=1
)

df_scopus['RemoteSensing'] = df_scopus.apply(
    lambda row: count_keywords(row['Title'], lstwordRS) +
                count_keywords(row['Abstract'], lstwordRS) +
                count_keywords(row['Author Keywords'], lstwordRS) +
                count_keywords(row['Index Keywords'], lstwordRS),
    axis=1
)

# ============================================================
# 4. Filtrar artigos que tenham pelo menos 1 ocorrência em cada categoria
# ============================================================
df_filtered = df_scopus[
    (df_scopus['Mine'] > 0) &
    (df_scopus['DeepLearning'] > 0) &
    (df_scopus['RemoteSensing'] > 0)
].copy()

# ============================================================
# 5. Exibir resultados
# ============================================================
print(f"Total de artigos originais: {len(df_scopus)}")
print(f"Artigos filtrados (mineração + DL + sensoriamento remoto): {len(df_filtered)}")

# (Opcional) Mostrar primeiros títulos dos artigos filtrados
print("\nExemplos de artigos selecionados:")
for i, (idx, row) in enumerate(df_filtered.head(10).iterrows()):
    print(f"{i+1}. {row['Title']}")