#!/bin/python3
import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples
from pathlib import Path
from io import StringIO  # Import StringIO directly from the 'io' module

try:
    path = Path('./snp-dists.txt')
    path.resolve(strict=True)
except FileNotFoundError:
    path = Path('./snp_matrix.txt')

print("Using file path:", path)

# Helper function to perform replacements using regular expressions
def replace_with_regex(text):
    return re.sub(
        r"\.consensus_threshold_\d+\.\d+_quality_\d+",
        '',
        text.replace(",", "\t")
            .replace('snp-dists 0.8.2\t', '')
            .replace('_contigs', '')
            .replace('_genomic', '')
            .replace("Consensus_", '')
    )

# Read the entire content of the file and perform replacements
with open(path, 'r') as file:
    content = file.read()
    content = replace_with_regex(content)

# Create a DataFrame from the modified content
df = pd.read_csv(StringIO(content), sep='\t', index_col=0)  # Use StringIO directly

# Handle any potential non-finite values in the DataFrame
df = df.replace([np.inf, -np.inf], np.nan).dropna(inplace=False)

numSamples = len(df.columns)

print("Found", numSamples, "samples in snp_matrix.txt")

if numSamples <= 2:
    print("This matrix has too few samples or has been melted. Sorry!")
    exit(0)
else:
    
# Calculate the appropriate figure size based on the number of samples

    if numSamples <= 20:
        figureSize = (10, 8)
    else:
        figureSize = (round(numSamples / 2), round(numSamples / 2.5))
        print("\n\nNumber of samples:", numSamples, "\nFigure size:", figureSize)

# Now, let's proceed to compute clusters
clusters = sch.linkage(df.values, method='complete', metric='euclidean')

# Create clustermap and get the order of rows and columns based on clustering
clustergrid = sns.clustermap(
    df, xticklabels=True, yticklabels=True, vmin=0, vmax=50, center=20,
    annot=True, annot_kws={'size': 6}, cbar_kws={"orientation": "vertical", "pad": 0.5},
    cmap='Reds_r', linecolor="white", linewidths=.01, fmt='d', dendrogram_ratio=0.1,
    col_cluster=True, row_cluster=True, figsize=figureSize
)
plt.setp(clustergrid.ax_heatmap.get_xticklabels(), rotation=45, ha='right')

# Suppress printing of dendrogram along the y-axis
clustergrid.ax_row_dendrogram.set_visible(False)
clustergrid.ax_col_dendrogram.set_visible(True)
row_order = clustergrid.dendrogram_row.reordered_ind
col_order = row_order

# Sort the DataFrame based on the cluster order
sorted_df = df.iloc[row_order, col_order]

# Compute the number of SNPs within the cluster per row
within_cluster_snps = sorted_df.apply(lambda row: row[row < 500].sum(), axis=1)

# Add 'Within_Cluster_SNPs' column to the sorted DataFrame
sorted_df['Within_Cluster_SNPs'] = within_cluster_snps.values

# Calculate silhouette scores for different numbers of clusters
silhouette_scores = []
upper_range = min(numSamples, 11)
for n_clusters in range(2, upper_range):
    kmeans = KMeans(n_clusters=n_clusters, n_init=10)
    cluster_labels = kmeans.fit_predict(sorted_df.values)
    silhouette_avg = silhouette_samples(sorted_df.values, cluster_labels).mean()
    silhouette_scores.append(silhouette_avg)

# Find the optimal number of clusters with the highest silhouette score
optimal_num_clusters = silhouette_scores.index(max(silhouette_scores)) + 2

# Use the optimal number of clusters to assign cluster labels and sort the DataFrame
kmeans = KMeans(n_clusters=optimal_num_clusters, n_init=10)
cluster_labels = kmeans.fit_predict(sorted_df.values)
sorted_df['Cluster'] = cluster_labels

# Sort the DataFrame first by 'Cluster' and then by 'Within_Cluster_SNPs'
sorted_df = sorted_df.sort_values(by=['Cluster', 'Within_Cluster_SNPs'], ascending=[True, True], kind="mergesort")

# Drop 'Cluster' and 'Within_Cluster_SNPs' columns
sorted_df = sorted_df.drop(['Cluster', 'Within_Cluster_SNPs'], axis='columns')
sorted_df = sorted_df.reindex(columns=sorted_df.index)

# Save the finally sorted, tab-delimited SNP matrix
sorted_df.to_csv('Final_snp_matrix.tsv', sep='\t', header=True, index=True)

# Create the reordered heatmap with correct values
heatmap = sns.clustermap(
    sorted_df, xticklabels=True, yticklabels=True, vmin=0, vmax=50, center=20,
    annot=True, annot_kws={'size': 6}, cbar_kws={"orientation": "vertical", "pad": 0.5},
    cmap='Reds_r', linecolor="white", linewidths=.01, fmt='d', dendrogram_ratio=0.15,
    col_cluster=True, row_cluster=True, figsize=figureSize
)
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
heatmap.ax_row_dendrogram.set_visible(False)
heatmap.ax_col_dendrogram.set_visible(False)

# Save the reordered heatmap as a PDF and PNG
heatmap.savefig('SNP_matrix.pdf')
heatmap.savefig('SNP_matrix.png')
heatmap.savefig('SNP_matrix_mqc.png')
plt.show()
plt.close()

print("Saved heatmap as SNP_matrix.{pdf,png}")
