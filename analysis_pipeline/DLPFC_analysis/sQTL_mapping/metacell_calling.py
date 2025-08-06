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

meta_size = 10

cell="Endothelial"
cohort="NG"
adata = sc.read_h5ad("harmony_excitatory_neurons.h5ad")


X = pd.DataFrame(adata.obsm['X_pca'])
X.index = adata.obs['individualID'].index

adata.obs["meta_cell"] = ""
meta_info = []

donor=pd.read_csv("/data/projects/11003054/e0950183/brain_sQTL/junc/total/Nature_sample",header=None)
donor_h5ad = np.unique(adata.obs['individualID'])
inter_donor = np.intersect1d(donor,donor_h5ad)

whole_label=pd.DataFrame(np.zeros((0,3)))
whole_label.columns=['meta','cell_id','ind_id']

for i in range(len(inter_donor)):
    sample_id = inter_donor[i]
    num = adata.obs['individualID'][adata.obs['individualID']==sample_id]
    if len(num)<meta_size:
       continue
    if len(num)<2*meta_size:
       newlabel = pd.DataFrame(np.zeros((len(num),3)))
       newlabel.columns=['meta','cell_id','ind_id']
       newlabel.loc[:,'meta']="0"
       newlabel.loc[:,'cell_id']=adata.obs['individualID'].index[adata.obs['individualID']==sample_id]
       newlabel.loc[:,'ind_id']=sample_id
       whole_label = pd.concat([whole_label,newlabel])
       continue
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
    newlabel.loc[:,'cell_id']=adata.obs['individualID'].index[adata.obs['individualID']==sample_id]
    unassigned = newlabel.loc[newlabel.iloc[:,0]==10000,'cell_id']
    newcluster=np.unique(labels)
    if len(newcluster)==1:
       newlabel.columns=['meta','cell_id']
       newlabel.loc[:,'meta']="0"
       newlabel.loc[:,'ind_id']=sample_id
       whole_label = pd.concat([whole_label,newlabel])
       continue
    
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


whole_label.loc[:,'combine_meta']=whole_label['ind_id']+':'+cohort+"_"+whole_label['meta']

len(np.unique(whole_label.loc[:,'combine_meta']))
       
len(np.unique(whole_label.loc[:,'ind_id']))

whole_label.to_csv("/home/users/nus/e0950183/scratch/brain_annot/metacell/"+cell+".csv",index=False)
