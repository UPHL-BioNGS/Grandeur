import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

# Read the SNP matrix file and write a tab-delimited version
with open("snp_matrix.txt", "r") as infile, open("snp_matrix.tsv", "w") as outfile:
    for line in infile:
        text = line.replace('snp-dists 0.8.2,', '').replace(",", "\t")
        outfile.write(text)

# Read the tab-delimited file into a DataFrame
df = pd.read_csv("snp_matrix.tsv", sep='\t')

# Add Total_SNPs column for sum of SNPs in each row
df['Total_SNPs'] = df.sum(numeric_only=True, axis=1)

# Sort dataframe by Total SNPs in each row
TotalRow_snps = df.sort_values(by='Total_SNPs', kind="mergesort")
print("totalRowSnps index: ", TotalRow_snps.index)
# Reorder the columns to mirror the row order
sorted_cluster_matrix = TotalRow_snps.reindex(columns=TotalRow_snps.index)
print("reindexed: ",sorted_cluster_matrix)
# Export sorted DataFrame to tab-delimited text file
sorted_cluster_matrix.to_csv('1_snp_matrix_initially_sorted.tsv', sep='\t', header=True, index=True)

###
# Prepare data for clustermap
mask_df = sorted_cluster_matrix
mySeqs = mask_df.columns.tolist()
print("mask_df", mask_df)
print("mySeqs: \n", mySeqs) 
cmap = 'Reds_r'

# Create clustermap
clustergrid = sns.clustermap(
    mask_df, xticklabels=True, yticklabels=True, vmin=0, vmax=100, center=30,
    annot=True, annot_kws={'size': 8}, cbar_kws={"pad": 0.5},
    cmap=cmap, linecolor="white", linewidths=.2, fmt='d', dendrogram_ratio=0.1, 
    col_cluster=True, row_cluster=True
) 

# Rotate x-labels 45 degrees to the right
plt.setp(clustergrid.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
plt.title('SNP matrix')
# Save the sorted clustermap
plt.savefig("SNPmatrix.pdf", dpi=120)
plt.savefig("SNPmatrix_mqc.png", dpi=120)
print("Initial Clustermap saved as 120 dpi in pdf and png format")
plt.show()
plt.close()
