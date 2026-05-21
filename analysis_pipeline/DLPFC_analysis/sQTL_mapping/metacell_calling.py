import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import pairwise_distances
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import networkx as nx
import anndata
import community as community_louvain
from sklearn.metrics import silhouette_score, calinski_harabasz_score

import sys

meta_size = int(sys.argv[1])

cell_type=sys.argv[2]
norm_choice=sys.argv[3]


if norm_choice=="log_normalize":
   adata = sc.read_h5ad("/data/projects/11003054/e0950183/AIDA_phaseIfreezeII/rds/"+cell_type+".h5ad")

if norm_choice=="SCTransform":
   adata = sc.read_h5ad("/data/projects/11003054/e0950183/AIDA_phaseIfreezeII/rds/"+cell_type+"_SCT.h5ad")

X = pd.DataFrame(adata.obsm['X_pca'])
X.index = adata.obs['DCP_ID'].index ##sample name

adata.obs["meta_cell"] = ""
meta_info = []


whole_label=pd.DataFrame(np.zeros((0,3)))
whole_label.columns=['meta','cell_id','ind_id']

inter_donor=np.unique(adata.obs['DCP_ID']).tolist()

performance = pd.DataFrame(np.nan, index=inter_donor, columns=['Silhouette','CH'])
for i in range(len(inter_donor)):
    sample_id = inter_donor[i]
    print(sample_id)
    num = adata.obs['DCP_ID'][adata.obs['DCP_ID']==sample_id]
    if len(num)<2*meta_size:
       newlabel = pd.DataFrame(np.zeros((len(num),3)))
       newlabel.columns=['meta','cell_id','ind_id']
       newlabel.loc[:,'meta']="0"
       newlabel.loc[:,'cell_id']=adata.obs['DCP_ID'].index[adata.obs['DCP_ID']==sample_id]
       newlabel.loc[:,'ind_id']=sample_id
       whole_label = pd.concat([whole_label,newlabel])
       continue
    print("pass1")
    pc_matrix = X.loc[num.index]
    tmp_adata = anndata.AnnData(X=np.array(pc_matrix))
    nbrs = NearestNeighbors(n_neighbors=meta_size).fit(pc_matrix)
    _,indices = nbrs.kneighbors(pc_matrix)
    G=nx.Graph()
    for m, neighbors in enumerate(indices):
      for n in neighbors[1:]:
        G.add_edge(m, n)
    partition = community_louvain.best_partition(G)
    labels = np.array([partition[i] for i in range(len(pc_matrix))])
    cluster = np.unique(labels)
    ### collapse clusters with less than meta size to larger clusters
    for unique_group in cluster:
       group_num = len(labels[labels==unique_group])
       if group_num<meta_size:
          labels[labels==unique_group]=10000
    newlabel = pd.DataFrame(labels)
    newlabel.loc[:,'cell_id']=adata.obs['DCP_ID'].index[adata.obs['DCP_ID']==sample_id]
    unassigned = newlabel.loc[newlabel.iloc[:,0]==10000,'cell_id']
    newcluster=np.unique(labels)
    if len(newcluster)==1:
       newlabel.columns=['meta','cell_id']
       newlabel.loc[:,'meta']="0"
       newlabel.loc[:,'ind_id']=sample_id
       whole_label = pd.concat([whole_label,newlabel])
       continue
    print("pass2")
    pc_centroid = pd.DataFrame(np.zeros(((len(newcluster)-1),50)))
    pc_centroid.index=newcluster[0:(len(newcluster)-1)]
    for pc_centroid_index in pc_centroid.index:
       cell_id_clu = newlabel.loc[newlabel.loc[:,0]==pc_centroid_index,'cell_id']
       cen_PC = pc_matrix.loc[cell_id_clu,:]
       pc_centroid.loc[pc_centroid_index,:]=np.mean(cen_PC,axis=0)
    
    for cell in unassigned:
       pc_cell=np.array(pc_matrix.loc[cell,:]).reshape(1,-1)
       D = pairwise_distances(pc_cell,pc_centroid,metric='euclidean')
       min_centroid = pc_centroid.index[D[0]==min(D[0])]
       newlabel.loc[newlabel.loc[:,'cell_id']==cell,0]=min_centroid
    newlabel.loc[:,'ind_id']=sample_id
    newlabel.loc[:,0]=np.array(newlabel.iloc[:,0]).astype(str)
    newlabel.columns=['meta','cell_id','ind_id']
    whole_label = pd.concat([whole_label,newlabel])
    ###compute silhouette_score
    if len(np.unique(newlabel.loc[:,'meta']))>1:
     performance.loc[sample_id,'Silhouette'] = silhouette_score(pc_matrix, np.array(newlabel.iloc[:,0]).astype(str), metric='euclidean')
    ###compute Calinski–Harabasz index
     performance.loc[sample_id,'CH'] = calinski_harabasz_score(pc_matrix, np.array(newlabel.iloc[:,0]).astype(str))


whole_label.loc[:,'combine_meta']=whole_label['ind_id']+':'+whole_label['meta']

len(np.unique(whole_label.loc[:,'combine_meta']))
       
len(np.unique(whole_label.loc[:,'ind_id']))

whole_label.to_csv("/data/projects/11003054/e0950183/ISSAC_revise/02_metacell_size/"+cell_type+norm_choice+str(meta_size)+"_meta.csv",index=False) ##metacell label

performance.to_csv("/data/projects/11003054/e0950183/ISSAC_revise/02_metacell_size/"+cell_type+norm_choice+str(meta_size)+"_performance.csv",index=True) ##performance of metacell construction

