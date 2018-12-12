import os
import numpy as np
import pandas as pd
from scanpy.tools.leiden import leiden
from scanpy.tools.louvain import louvain


def cluster(adata_filt, resolution=1.0, key_added="cluster", use_louvain=False):
    if use_louvain:
        louvain(adata_filt, resolution=resolution, key_added=key_added)
    else:
        leiden(adata_filt, resolution=resolution, key_added=key_added)

    shift_clusters(adata_filt, key_added)


def subcluster(adata_filt, cluster, resolution=0.4, cluster_key="cluster"):
    dtype = adata_filt.obs[cluster_key].dtype
    if dtype in ("object", str, pd.api.types.CategoricalDtype):
        cluster = str(cluster)
    #else:
    #    cluster = int(cluster)

    if cluster_key[-1].isdigit():
        base, version = cluster_key.split("_R")
        key_added = f"{base}_R{int(version)+1}"
    else:
        key_added = f"{cluster_key}_R1"
    louvain(
        adata_filt,
        resolution=resolution,
        key_added=key_added,
        restrict_to=(cluster_key, [cluster]),
    )

    # we want to preserve the ordering of old clusters and just add more.
    all_new_clusters = adata_filt.obs[key_added].cat.categories
    just_new_clusters = all_new_clusters[all_new_clusters.str.startswith(cluster)]
    if len(just_new_clusters) == 1:
        print(f"Wasn't able to subcluster with resolution `{resolution}`].")
        return

    old_max = adata_filt.obs[cluster_key].astype(int).max()

    tmp = pd.DataFrame()
    tmp["old"] = adata_filt.obs[cluster_key].astype("object")
    tmp["tmp"] = adata_filt.obs[key_added].astype("object")
    tmp["suff"] = adata_filt.obs[key_added].str.extract(",(\d+)", expand=False)
    tmp["new"] = tmp.old
    remapper = lambda row: row.old if row.suff == "0" else str(old_max + int(row.suff))
    tmp.loc[tmp.old != tmp.tmp, "new"] = tmp.loc[tmp.old != tmp.tmp].apply(
        remapper, axis=1
    )

    adata_filt.obs[key_added] = tmp["new"].astype(str).astype("category")
    order_clusters(adata_filt, key_added)
    print(f"Updated clusters under `adata_redux.obs['{key_added}']`.")


def shift_clusters(adata, key):
    if adata.obs[key].astype(int).min() == 0:
        adata.obs[key] = (
            (adata.obs[key].astype(int) + 1).astype("str").astype("category")
        )


def order_clusters(adata, key):
    adata.obs[key].cat.reorder_categories(
        adata.obs[key].cat.categories.astype(int).sort_values().astype(str),
        inplace=True
    )


__api_objects__ = {
    "cluster": cluster,
    "subcluster": subcluster,
}
