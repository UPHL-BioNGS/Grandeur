#!/bin/python3

##########################################
# written by Stephen Beckstrom-Sternberg #
# for creating SNP heatmaps for Grandeur #
##########################################

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from io import StringIO

# Read the SNP matrix file
with open("snp_matrix.txt", "r") as infile:
    lines = infile.readlines()
numSamples = len(lines) -1 #counts data lines

# Remove 'snp-dists 0.8.2', '_contigs' and '_genomic', & replace commas with tabs
cleaned_lines = [line.replace('snp-dists 0.8.2\t', '').replace('snp-dists 0.8.2,', '').
                 replace(",", "\t").replace('_contigs', '').replace('_genomic', '').replace("^\t", '')
                 for line in lines]

# Combine the cleaned lines into a single string instead of a file
snp_matrix_string = "\n".join(cleaned_lines)

# Read the tab-delimited string into a DataFrame
df = pd.read_csv(StringIO(snp_matrix_string), sep='\t')

#Define colormap for heatmap
cmap = 'Reds_r'

# Add Total_SNPs column for the sum of SNPs in each row
df['Total_SNPs'] = df.sum(numeric_only=True, axis="columns")

# Sort the DataFrame by Total SNPs in each row
TotalRow_snps = df.sort_values(by='Total_SNPs', kind="mergesort", axis="rows")

#Remove Total_SNPs column
sorted_cluster_matrix = TotalRow_snps.drop(['Total_SNPs'], axis="columns")

# Reorder the columns to mirror the row order
sorted_cluster_matrix=sorted_cluster_matrix.reindex(columns=sorted_cluster_matrix.index)

#Change output figure size tuple based on number of samples
if   (numSamples <= 20): 
    figureSize = (10, 8)
elif (numSamples <= 40): 
    figureSize = (20, 16)
elif (numSamples <= 60): 
    figureSize = (30, 24)
else: 
    figureSize = (40, 32)
print("\n\nNumber of samples: ", numSamples,"\nFigure size: ", figureSize)

# Compute clusters
clusters = sch.linkage(sorted_cluster_matrix.values, method='complete', metric='euclidean')

# Create clustermap and get the order of rows and columns based on clustering
clustergrid = sns.clustermap(
    sorted_cluster_matrix, 
    xticklabels=True, 
    vmin=0, 
    vmax=50, 
    center=20,
    annot=True, 
    annot_kws={'size': 6}, 
    cbar_kws={"pad": 0.5},
    cmap=cmap, 
    linecolor="white", 
    linewidths=.2, 
    fmt='d', 
    dendrogram_ratio=0.1,
    col_cluster=True, 
    row_cluster=True, 
    figsize=figureSize
    )
plt.setp(clustergrid.ax_heatmap.get_xticklabels(), rotation=45, ha='right')

# Suppress printing of dendrogram along the y-axis
clustergrid.ax_row_dendrogram.set_visible(False)
clustergrid.ax_col_dendrogram.set_visible(False)

row_order = clustergrid.dendrogram_row.reordered_ind
col_order = row_order
# Sort the DataFrame based on the cluster order
sorted_df = sorted_cluster_matrix.iloc[row_order, col_order]

# Compute the number of SNPs within the cluster per row
within_cluster_snps = sorted_df.apply(lambda row: row[row < 500].sum(), axis=1)

# Add 'Within_Cluster_SNPs' column to the sorted DataFrame
sorted_df['Within_Cluster_SNPs'] = within_cluster_snps.values

# Calculate silhouette scores for different numbers of clusters
silhouette_scores = []

if numSamples < 11:
    upper_range = numSamples
else:
    upper_range = 11

for n_clusters in range(2, upper_range):
    kmeans = KMeans(n_clusters=n_clusters, n_init=10)
    cluster_labels = kmeans.fit_predict(sorted_df.values)
    silhouette_scores.append(silhouette_score(sorted_df.values, cluster_labels))

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
    cmap=cmap, linecolor="white", linewidths=.1, fmt='d', dendrogram_ratio=0.1,
    col_cluster=False, row_cluster=False, figsize=figureSize
)
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right')

# Save the reordered heatmap as a PDF
heatmap.savefig('SNP_matrix.pdf')
heatmap.savefig('SNP_matrix.png')
heatmap.savefig('SNP_matrix_mqc.png')

plt.show()
plt.close()

print("Saved heatmap as Heatmap.{pdf,png}")