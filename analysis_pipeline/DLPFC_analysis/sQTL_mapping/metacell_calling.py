import sys
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import community as community_louvain
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import pairwise_distances, silhouette_score, calinski_harabasz_score

# Usage: python metacell_construction.py <meta_size> <cell_type_prefix> <log_normalize|SCTransform> [donor_column] [seed]
meta_size   = int(sys.argv[1])                                  # minimum number of cells per metacell
cell_type   = sys.argv[2]                                       # prefix of the h5ad file
norm_choice = sys.argv[3]                                       # log_normalize or SCTransform
donor_col   = sys.argv[4] if len(sys.argv) > 4 else "DCP_ID"    # donor/sample ID column in adata.obs
seed        = int(sys.argv[5]) if len(sys.argv) > 5 else 0      # random seed for Louvain

if norm_choice == "log_normalize":
    adata = sc.read_h5ad(cell_type + "_log.h5ad")
elif norm_choice == "SCTransform":
    adata = sc.read_h5ad(cell_type + "_SCT.h5ad")
else:
    sys.exit("norm_choice must be 'log_normalize' or 'SCTransform'")

if not adata.obs_names.is_unique:
    sys.exit("Cell barcodes (obs_names) must be unique.")

X = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
n_pc = X.shape[1]

label_list = []
donors = sorted(adata.obs[donor_col].unique().tolist())
performance = pd.DataFrame(np.nan, index=donors, columns=["Silhouette", "CH"])

for sample_id in donors:
    print(sample_id)
    cells = adata.obs_names[adata.obs[donor_col] == sample_id]
    n_cells = len(cells)

    # Donors with too few cells: all cells form a single metacell
    if n_cells < 2 * meta_size:
        label_list.append(pd.DataFrame({"meta": "0", "cell_id": cells, "ind_id": sample_id}))
        continue

    pc_matrix = X.loc[cells]

    # kNN graph (first neighbour is the cell itself and is skipped)
    nbrs = NearestNeighbors(n_neighbors=meta_size).fit(pc_matrix.values)
    _, indices = nbrs.kneighbors(pc_matrix.values)
    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    for m, neighbors in enumerate(indices):
        for n in neighbors[1:]:
            if n != m:
                G.add_edge(m, int(n))

    partition = community_louvain.best_partition(G, random_state=seed)
    labels = np.array([partition[k] for k in range(n_cells)])

    # Clusters smaller than meta_size are dissolved and their cells reassigned
    cluster_ids, cluster_sizes = np.unique(labels, return_counts=True)
    keep = cluster_ids[cluster_sizes >= meta_size]

    if len(keep) <= 1:
        # No (or only one) cluster large enough: the donor forms a single metacell
        label_list.append(pd.DataFrame({"meta": "0", "cell_id": cells, "ind_id": sample_id}))
        continue

    # Assign cells in dissolved clusters to the nearest retained centroid
    unassigned = ~np.isin(labels, keep)
    if unassigned.any():
        centroids = np.vstack([pc_matrix.values[labels == c].mean(axis=0) for c in keep])
        D = pairwise_distances(pc_matrix.values[unassigned], centroids, metric="euclidean")
        labels[unassigned] = keep[np.argmin(D, axis=1)]

    meta = labels.astype(str)
    label_list.append(pd.DataFrame({"meta": meta, "cell_id": cells, "ind_id": sample_id}))

    # Clustering quality within this donor
    performance.loc[sample_id, "Silhouette"] = silhouette_score(pc_matrix.values, meta, metric="euclidean")
    performance.loc[sample_id, "CH"] = calinski_harabasz_score(pc_matrix.values, meta)

whole_label = pd.concat(label_list, ignore_index=True)
whole_label["combine_meta"] = whole_label["ind_id"].astype(str) + ":" + whole_label["meta"]

print("Number of metacells:", whole_label["combine_meta"].nunique())
print("Number of donors:", whole_label["ind_id"].nunique())

prefix = cell_type + norm_choice + str(meta_size)
whole_label.to_csv(prefix + "_meta.csv", index=False)              # metacell labels
performance.to_csv(prefix + "_performance.csv", index=True)        # metacell construction performance

