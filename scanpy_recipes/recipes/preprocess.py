import os
import numpy as np
import pandas as pd
from scanpy.preprocessing import simple as pp
from scanpy.preprocessing.bbknn import bbknn
from scanpy.neighbors import neighbors
from scanpy.tools.umap import umap
from scanpy.tools.leiden import leiden
from scanpy.tools.louvain import louvain


def preprocess(adata_raw, n_top_genes=1000, scale=False):
    adata = pp.normalize_per_cell(adata_raw, copy=True)
    adata.uns["raw_dtype"] = "normalized count"
    adata.raw = adata

    adata.uns["n_top_genes"] = n_top_genes
    adata_filt = pp.filter_genes_dispersion(adata, n_top_genes=n_top_genes,
                                               flavor="cell_ranger", log=True, copy=True)
    pp.normalize_per_cell(adata_filt)

    pp.log1p(adata)
    pp.log1p(adata_filt)
    if scale:
        pp.scale(adata_filt)

    return adata, adata_filt


def dimensionality_reduction(adata_filt, n_comps=50, n_neighbors=10, min_dist=0.5, metric="correlation",
                             is_aggregation=False):
    pp.pca(adata_filt, n_comps=n_comps, svd_solver="arpack")

    if is_aggregation:
        n_batches = len(adata.obs["batch"].unique())
        bbknn(
            adata_filt,
            metric="angular",
            neighbors_within_batch=n_neighbors//n_batches,
            approx=True
        )
    else:
        neighbors(adata_filt, n_neighbors=n_neighbors, metric=metric)

    adata_filt.uns["neighbors"]["params"]["metric"] = metric
    umap(adata_filt, min_dist=min_dist, n_components=3)
    adata_filt.obsm["X_umap_3d"] = adata_filt.obsm["X_umap"].copy()
    umap(adata_filt, min_dist=min_dist, n_components=2)


def cluster(adata_filt, key_added="cluster", use_louvain=False):
    if use_louvain:
        louvain(adata_filt, key_added=key_added)
    else:
        leiden(adata_filt, key_added=key_added)

    #sc.tl.paga(adata_filt)
    if adata_filt.obs[key_added].astype(int).min() == 0:
        adata_filt.obs[key_added] = (adata_filt.obs[key_added].astype(int) + 1).astype("category")

__api_objects__ = {
    "preprocess": preprocess,
    "dimensionality_reduction": dimensionality_reduction,
    "cluster": cluster
}


#def find_markers(adata_filt):
#    sc.tl.rank_genes_groups(adata_filt, "louvain", method="wilcoxon", use_raw=False,
#                            key_added="rank_genes_groups-wilcoxon")
#    sc.tl.rank_genes_groups(adata_filt, "louvain", method="t-test_overestim_var", use_raw=False,
#                            key_added="rank_genes_groups-ttest")
#    save_file = os.path.join(OUTPUT_DIRS[adata_filt.uns["sampleid"]], "top5_genes-wilcox.pdf")
#    fig = sc.pl.rank_genes_groups_matrixplot(adata_filt, groupby="louvain", n_genes=5,
#                                             key="rank_genes_groups-wilcoxon", use_raw=False,
#                                             cmap=cmo.cm.tempo,
#                                             figsize=(36, 6), save=True)
#    os.rename(os.path.join("figures", "matrixplot.pdf"), save_file)
#
#    save_file = os.path.join(OUTPUT_DIRS[adata_filt.uns["sampleid"]], "top5_genes-ttest.pdf")
#    fig = sc.pl.rank_genes_groups_matrixplot(adata_filt, groupby="louvain", n_genes=5,
#                                             key="rank_genes_groups-ttest", use_raw=False,
#                                             cmap=cmo.cm.tempo,
#                                             figsize=(36, 6), save=True)
#    os.rename(os.path.join("figures", "matrixplot.pdf"), save_file)
#
#
#def print_markers(adata_filt, method, outfile):
#    names = adata_filt.uns[f"rank_genes_groups-{method}"]["names"]
#    marker_df = pd.DataFrame(names)
#    marker_df.index.name = "Rank"
#    marker_df.columns.name = "Cluster ID"
#
#    marker_df.to_csv(os.path.join(OUTPUT_DIRS[adata_filt.uns["sampleid"]], outfile))
