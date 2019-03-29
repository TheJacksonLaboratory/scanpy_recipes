import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

classification.index = classification.index.str.slice(0, 16)
classification.index = classification.index + "-1"
intersect = classification.index.intersection(adata_raw.obs_names)

adata_raw.obs["qc_classification"] = "negative"
adata_raw.obs.loc[intersect, "qc_classification"] = classification.loc[intersect,
        "global_class"]
adata_raw.obs["qc_tag"] = "negative"
adata_raw.obs.loc[intersect, "qc_tag"] = classification.loc[intersect, "tag_class"]
adata_raw.obs.qc_tag = adata_raw.obs.qc_tag.astype("category")
adata_raw.obs.qc_tag.cat.reorder_categories(["tag_3", "tag_4", "tag_5", "negative"],
        inplace=True)
adata_raw.obs["qc_classification"] = adata_raw.obs["qc_classification"].astype("category")
adata_raw.obs["qc_classification"].cat.reorder_categories(["singlet", "multiplet",
    "negative"], inplace=True)

adata_raw = adata_raw[adata_raw.obs.qc_classification != "multiplet", :]


hto_raw = pd.read_csv("raw_data/AW18002-hto/AW18002_hto_counts.tsv", header=0,
        index_col=0).T

hto = hto_raw.sort_values("total_reads", ascending=False).head(24000).iloc[:,1:4]
hto.columns = ["tag_3", "tag_4", "tag_5"]



def CLR_norm(counts):
    inner = np.log1p(counts.values / np.exp(np.sum(np.log1p(counts[counts > 0]), axis=1) / counts.shape[1])[:,None])
    return pd.DataFrame(data=inner, index=counts.index, columns=counts.columns)


def load_tag_counts(counts_file, ):
    """
    Assumes tag counts matrix was generated using `CITE-seq-Count`:
    https://github.com/Hoohm/CITE-seq-Count
    """
    pass


def hto_demux(raw_htos, n_clusters=None):
    # 1
    htos = raw_htos.columns
    norm_htos = CLR_norm(raw_htos.copy())
    if n_clusters is None:
        n_clusters = len(htos) + 1

    # using Kmeans instead of kmediods to start
    classification = pd.DataFrame(index=raw_htos.index)
    k = KMeans(n_clusters=n_clusters).fit(raw_htos.values)
    classification["k_label"] = k.labels_

    # 2
    norms = pd.DataFrame(index=np.unique(k.labels_), columns=htos)
    for cluster, group in classification.groupby("k_label"):
        cells = group.index
        norms.loc[cluster, :] = norm_htos.loc[cells, :].mean(axis=0)

    # 3
    # Maybe use sm.ECDF instead of fitting a nbinom
    for hto in htos:
        min_cluster = norms[hto].astype(float).idxmin()
        print(min_cluster)
        min_cells = classification[classification.k_label == min_cluster].index
        background_expression = raw_htos.loc[min_cells, hto]
        threshold = np.percentile(background_expression.values, 95)
        indicator = raw_htos[hto] > threshold
        classification[hto] = indicator.astype(int)
        print(hto, threshold)

    hto_sums = classification.iloc[:, 1:].sum(axis=1)
    classification["global_class"] = "negative"
    classification.loc[hto_sums == 1, "global_class"] = "singlet"
    classification.loc[hto_sums > 1, "global_class"] = "multiplet"

    singlets = classification.loc[hto_sums == 1, classification.columns[1:len(htos)+1]]
    classification["tag_class"] = "negative"
    classification.loc[singlets.index, "tag_class"] = singlets.idxmax(axis=1)
    return classification
