import scanpy as sc
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
cell="Inh"
adata_combined=sc.read_h5ad("/data/projects/11003054/e0950183/brain_sQTL/h5ad/harmony_Inhibitory.h5ad")

import pandas as pd

##

gene_expression = adata_combined.X.toarray() if hasattr(adata_combined.X, "toarray") else adata_combined.X


for j in range(8):
    gene_pc1_correlations = []

# Compute Pearson correlation for each gene
    pc1_scores = adata_combined.obsm['X_pca_harmony'][:, j]
    adata_combined.obs['PC'+str(j+1)]=pc1_scores
    for i in range(gene_expression.shape[1]):
        corr, _ = pearsonr(gene_expression[:, i], pc1_scores)
        gene_pc1_correlations.append(corr)

# Convert to a numpy array for easier handling
    gene_pc1_correlations = np.array(gene_pc1_correlations)
    gene_pc1_assoc = pd.DataFrame({
     'Gene': adata_combined.var_names,
     'PC1_Correlation': gene_pc1_correlations})
    gene_pc1_assoc = gene_pc1_assoc.reindex(gene_pc1_assoc['PC1_Correlation'].abs().sort_values(ascending=False).index)
    gene_pc1_assoc.to_csv(cell+"_PC"+str(j+1)+".csv",index=False)
    top_genes = gene_pc1_assoc.head(5)
    gene_name=top_genes['Gene'].to_list()
    sc.pl.umap(adata_combined,color=["PC"+str(1+j),gene_name[0],gene_name[1]])
    plt.savefig(cell+"_PC"+str(j+1)+"_state.png",dpi=300)
    plt.close()
###save correlation
