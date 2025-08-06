cell="Inhibitory"
import pandas as pd

meta_info=pd.read_csv("/data/projects/11003054/e0950183/brain_sQTL/metacell/"+cell+"_meta.txt",sep="\t")
import scanpy as sc
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt

adata_combined=sc.read_h5ad("/data/projects/11003054/e0950183/brain_sQTL/h5ad/harmony_Inhibitory.h5ad")

for i in range(8):
    pc1_scores = adata_combined.obsm['X_pca_harmony'][:, i]
    adata_combined.obs['PC1']=pc1_scores

### PC1 for each metacell computing

    total_sample=list(set(meta_info['identify'].to_list()))

    total_PC=pd.DataFrame(np.nan,index=total_sample,columns=["PC"])
    for sample in total_sample:
      num=meta_info[meta_info['identify']==sample]['barcode'].to_list()
      PC=np.mean(adata_combined.obs['PC1'][num].tolist())
      total_PC.loc[sample,'PC']=PC

    total_PC.to_csv("meta"+cell+"PC"+str(i+1)+"_value.csv")
