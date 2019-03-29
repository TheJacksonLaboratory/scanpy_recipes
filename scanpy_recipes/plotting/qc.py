import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator
import seaborn as sns
import cmocean as cmo

from scanpy.api import pl
from scanpy.readwrite import read_10x_h5



def umi_rank_plot(adata_redux, return_fig=False):
    for h5_file in glob.glob(f"{adata_redux.uns['input_dir']}/raw_*matri*.h5"):
        if os.path.exists(h5_file):
            break
    else:
        raise IOError(f"Can't find raw matrix h5 file under [{input_dir}].")
    raw = read_10x_h5(h5_file, adata_redux.uns["genome"])

    min_umis = adata_redux.uns["10x_umi_cutoff"]

    umi_counts = raw.X.sum(axis=1).A1
    all_barcodes = np.sort(umi_counts)[::-1]
    cells = all_barcodes[all_barcodes > min_umis]
    n_cells = len(cells)

    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    #fig, ax = plt.subplots(figsize=(4, 3.5), dpi=300)
    ax.plot(cells, color="green", lw=3, label="Called cells", zorder=4)
    ax.plot(all_barcodes, color="0.7", lw=2, label="All barcodes")
    ax.plot([1, n_cells], [min_umis]*2, ls="-", color="0.7", lw=1)
    ax.plot([n_cells]*2, [0, min_umis], ls="-", color="0.7", lw=1)

    ax.set_xlabel("Barcodes")
    ax.set_ylabel("UMI Counts")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, ax.get_xlim()[1])
    ax.set_title(adata_redux.uns["sampleid"])

    ax.xaxis.set_major_locator(LogLocator())
    ax.yaxis.set_major_locator(LogLocator())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(LogLocator(subs=(2,5)))
    ax.yaxis.set_minor_formatter(FuncFormatter(lambda x, pos: f"{str(int(x))[0]}"))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x)}" if x < 1e4 else f"{int(x/1e3)}k"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{int(x)}" if x < 1e4 else f"{int(x/1e3)}k"))

    sns.despine(fig, ax)
    ax.legend(loc="upper right", frameon=False)
    ax.grid(which="both", axis="both", ls=":")

    fig.tight_layout()

    if return_fig:
        return fig


def cut_violin(ax, threshold=0., cut_above=False, color="r"):
    def cut_verts(polycol, threshold, cut_above, color):
        v = polycol.get_paths()[0].vertices
        if cut_above:
            ind = v[:,1] <= threshold
        else:
            ind = v[:,1] >= threshold
        polycol.set_verts([v[ind]])
        polycol.set(facecolor=color)
    # ax.collections[0] = first violins
    # ax.collections[1] = first jitter
    cut_verts(ax.collections[2], threshold, cut_above=cut_above, color=color)


def qc_violins(adata, return_fig=False):
    keys = adata.uns["obs_titles"].keys()
    #keys = sorted(filter(lambda s: not s.startswith("qc"), adata.obs_keys()))
    N = len(keys)

    thresholds = adata.uns.get("qc_cell_filter", None)
    use_thresholds = dict(zip(keys, [False]*N))
    if thresholds:
        use_thresholds = dict((key, thresholds.get("threshold_" + key, None)) for key in keys)

    flipped_keys = {"percent_mito", "hemoglobin_counts"}

    fig, axs = plt.subplots(1, N, figsize=(4*N, N), dpi=200)
    for ax, key in zip(axs.flat, keys):
        threshold = use_thresholds[key]
        print(key, threshold)
        params = dict(color="0.9", show=False, ax=ax, cut=0, gridsize=300,
                      linewidth=0.50)
        if threshold:
            ax = pl.violin(adata, key, **params)
            ax = pl.violin(adata, key, jitter=False, **params)
            cut_violin(ax, threshold=threshold,
                       cut_above=False if key in flipped_keys else True)
            ax.axhline(threshold, xmin=0.25, xmax=0.75, color="r")
            #return ax
        else:
            pl.violin(adata, key, **params)

        if adata.obs[key].sum() < 1:
            ax.set_ylim(0, 1)

        ax.set_title(adata.uns["obs_titles"][key])
        ax.set_xticks([])
        ax.set_ylabel("")
        sns.despine(fig, ax)

    fig.tight_layout()
    fig.text(
        x=0.02, y=0.5, s=f"Sample: {adata.uns['sampleid']}",
        rotation=90, va='center', ha='right', fontsize='xx-large'
    )
    fig.subplots_adjust(left=0.05)

    if return_fig:
        return fig


def _scat(fig, ax, adata_trial, key, cmap="Reds", cbar=True, sort_top=True):
    pdata = adata_trial.obs[["total_counts", "n_genes_by_counts", key]]
    if isinstance(pdata[key].dtype, pd.api.types.CategoricalDtype):
        pdata[key] = pdata[key].cat.codes
    pdata = pdata.sort_values(key, ascending=sort_top)

    p = ax.scatter(pdata["total_counts"],
                   pdata["n_genes_by_counts"],
                   c=pdata[key].values, linewidths=0.05, edgecolors="k",
                   s=12, alpha=0.6, cmap=cmap)
    if cbar:
        fig.colorbar(p, ax=ax)
    ax.set_title(adata_trial.uns["obs_titles"].get(key, "Overall QC"))
    ax.set_xlabel("UMIs")
    ax.set_ylabel("Genes")
    sns.despine(fig, ax)


def genes_umis_scatter(adata_trial, return_fig=False):
    keys = sorted(filter(lambda s: not s.startswith("qc") and not s.startswith("n_"), adata_trial.obs_keys()))
    L = len(keys) + 1

    fig, axs = plt.subplots(2, L//2, figsize=(L*4//2, 2*3))
    for ax, key in zip(axs.flatten(), keys):
        _scat(fig, ax, adata_trial, key)
    redblue = sns.blend_palette([sns.xkcd_rgb["light red"], "0.9"], 2, as_cmap=True)
    _scat(fig, axs.flatten()[-1], adata_trial, "qc_fail", cmap=redblue, cbar=False, sort_top=False)

    fig.tight_layout()

    if return_fig:
        return fig


def qc_pass_fail(adata_trial, return_fig=False):
    pdata = adata_trial.obs[["total_counts", "n_genes_by_counts", "qc_fail"]]
    passing = pdata.qc_fail == "pass"
    pdata_pass = pdata.loc[passing, :]
    pdata_fail = pdata.loc[~passing, :]
    params = dict(linewidths=0.05, edgecolors="k", s=12, alpha=0.6)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(pdata_pass["total_counts"], pdata_pass["n_genes_by_counts"],
               color="0.9", label="Pass", **params)
    ax.scatter(pdata_fail["total_counts"], pdata_fail["n_genes_by_counts"],
               color=sns.xkcd_rgb["light red"], label="Fail", **params)

    ax.set_xlabel("UMIs")
    ax.set_ylabel("Genes")
    ax.set_title(adata_trial.uns["sampleid"])
    sns.despine(fig, ax)

    ax.legend(bbox_to_anchor=(1.01, 1), title="QC", frameon=False)
    fig.subplots_adjust(right=0.75)

    if return_fig:
        return fig


__api_objects__ = {
    "umi_rank_plot": umi_rank_plot,
    "qc_pass_fail": qc_pass_fail,
    "qc_violins": qc_violins,
    "genes_umis_scatter": genes_umis_scatter,
}
