import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator
import seaborn as sns
import cmocean as cmo

from scanpy.api import pl
from scanpy.readwrite import read_10x_h5


#def plot_qc(adata, config):
#    fig, axarr = plt.subplots(1, 2, dpi=300, figsize=(10, 4))
#    fig.patch.set_facecolor('white')
#    kwargs = dict(s=16, lw=0.25, edgecolors='none')
#
#    obs = adata.obs.copy()
#    obs = obs.sort_values('sequencing_saturation', ascending=False)
#    p1 = axarr[0].scatter(obs['n_counts'], obs['n_genes'], c=obs['sequencing_saturation'],
#                          cmap=cmo.cm.matter, **kwargs)
#    axarr[0].set_xlabel('UMIs')
#    axarr[0].set_ylabel('Genes')
#    fig.colorbar(p1, ax=axarr[0], label='Sequencing Saturation')
#
#    obs = obs.sort_values('percent_mito', ascending=True)
#    p2 = axarr[1].scatter(obs['n_counts'], obs['n_genes'], c=obs['percent_mito'],
#                       cmap=cmo.cm.matter, **kwargs)
#    axarr[1].set_xlabel('UMIs')
#    axarr[1].set_ylabel('Genes')
#    fig.colorbar(p2, ax=axarr[1], label='mtRNA content')
#
#    [sns.despine(fig=fig, ax=ax) for ax in axarr.flat]
#    fig.tight_layout()
#    fig.savefig(os.path.join(config.output_dirs[adata.uns['sampleid']], 'qc.pdf'),
#                format='pdf', dpi=300)
#
#
#def plot_scatters(adata, basis='tsne', raw=True, components='1,2', **kwargs):
#    sampleid = adata.uns['sampleid']
#    pl.scatter(adata, basis=basis, color='louvain', projection='3d',
#                  size=4, title='Clusters', save="_clusters_3d.pdf", **kwargs)
#    os.rename('figures/umap_clusters_3d.pdf',
#              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_clusters_3d.pdf"))
#    pl.scatter(adata, basis=basis, color='louvain',
#                  size=4, title='Clusters', save="_clusters-labeled.pdf",
#                  components=adata.uns['umap_view_components'], legend_loc='on data', **kwargs)
#    os.rename('figures/umap_clusters-labeled.pdf',
#              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_clusters.pdf"))
#    pl.scatter(adata, basis=basis, color='n_counts', color_map=cmo.cm.dense,
#                  size=4, title='UMI counts', save="_umi_counts.pdf",
#                  components=adata.uns['umap_view_components'], **kwargs)
#    os.rename('figures/umap_umi_counts.pdf',
#              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_umi_counts.pdf"))
#    pl.scatter(adata, basis=basis, color='n_genes', color_map=cmo.cm.dense,
#                  size=4, title='Gene counts', save="_gene_counts.pdf",
#                  components=adata.uns['umap_view_components'], **kwargs)
#    os.rename('figures/umap_gene_counts.pdf',
#              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_gene_counts.pdf"))
#    pl.scatter(adata, basis=basis, color='percent_mito', color_map=cmo.cm.matter,
#                  size=4, title='mtRNA fraction', save="_mtrna_fraction.pdf",
#                  components=adata.uns['umap_view_components'], **kwargs)
#    os.rename('figures/umap_mtrna_fraction.pdf',
#              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_mtrna_fraction.pdf"))
#    pl.scatter(adata, basis=basis, color='sequencing_saturation', color_map=cmo.cm.matter,
#                  size=4, title='Sequencing Saturation', save="_sequencing_saturation.pdf",
#                  components=adata.uns['umap_view_components'], **kwargs)
#    os.rename('figures/umap_sequencing_saturation.pdf',
#              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_sequencing_saturation.pdf"))
#
#    cell_cycle_genes = pd.read_csv('/home/flynnb/projects/singlecell/cell_cycle_gene_sets.csv', header=0)
#    cell_cycle_genes = cell_cycle_genes.apply(lambda x: x.str.strip())
#    n_groups = cell_cycle_genes.shape[1]
#    for k in range(n_groups):
#        group = cell_cycle_genes.columns[k]
#        group_genes = adata.var_names[adata.var_names.isin(cell_cycle_genes.iloc[:, k])]
#        key = f'cell_cycle_{group}'.replace('/', '-')
#        adata.obs[key] = adata[:, group_genes].X.sum() / adata.obs['n_counts']
#        if not (adata.obs[key] > 0).any():
#            adata.obs.ix[-1, key] = 1e-6
#
#        pl.scatter(adata, basis=basis, color=key, color_map=cmo.cm.tempo, legend_fontsize=0,
#                      size=4, title=group, save=f"_{key}.pdf",
#                      components=adata.uns['umap_view_components'], **kwargs)
#        os.rename(f'figures/umap_{key}.pdf',
#                  os.path.join(OUTPUT_DIRS[sampleid], f'{sampleid}_{key}.pdf'))


def umi_rank_plot(adata_redux, return_fig=False):
    raw_h5 = os.path.join(os.path.dirname(adata_redux.uns["input_file"]),
                          "raw_gene_bc_matrices_h5.h5")
    raw = read_10x_h5(raw_h5, adata_redux.uns["genome"])

    min_umis = adata_redux.uns["10x_umi_cutoff"]

    umi_counts = raw.X.sum(axis=1).A1
    all_barcodes = np.sort(umi_counts)[::-1]
    cells = all_barcodes[all_barcodes > min_umis]
    n_cells = len(cells)

    #fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    fig, ax = plt.subplots(figsize=(4, 3.5), dpi=300)
    ax.plot(cells, color="green", lw=3, label="Called cells", zorder=4)
    ax.plot(all_barcodes, color="0.7", lw=2, label="All barcodes")
    ax.plot([1, n_cells], [min_umis]*2, ls="-", color="0.7", lw=1)
    ax.plot([n_cells]*2, [0, min_umis], ls="-", color="0.7", lw=1)

    ax.set_xlabel("Barcodes")
    ax.set_ylabel("UMI Counts")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, ax.get_xlim()[1])

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
    keys = sorted(filter(lambda s: not s.startswith("qc"), adata.obs_keys()))
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
        if threshold:
            params = dict(color="0.9", show=False, ax=ax, cut=0, gridsize=300
                          linewidth=0.25)
            ax = pl.violin(adata, key, **params)
            ax = pl.violin(adata, key, jitter=False, **params)
            cut_violin(ax, threshold=threshold,
                       cut_above=False if key in flipped_keys else True)
            ax.axhline(threshold, xmin=0.25, xmax=0.75, color="r")
            #return ax
        else:
            pl.violin(adata, key, show=False, ax=ax, cut=0, color="0.9")

        ax.set_title(adata.uns["obs_titles"][key])
        ax.set_xticks([])
        ax.set_ylabel("")
        sns.despine(fig, ax)

    fig.tight_layout()

    if return_fig:
        return fig


def _scat(fig, ax, adata_trial, key, cmap="Reds", cbar=True, sort_top=True):
    pdata = adata_trial.obs[["n_counts", "n_genes", key]]
    if isinstance(pdata[key].dtype, pd.api.types.CategoricalDtype):
        pdata[key] = pdata[key].cat.codes
    pdata = pdata.sort_values(key, ascending=sort_top)

    p = ax.scatter(pdata["n_counts"],
                   pdata["n_genes"],
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
    pdata = adata_trial.obs[["n_counts", "n_genes", "qc_fail"]]
    passing = pdata.qc_fail == "pass"
    pdata_pass = pdata.loc[passing, :]
    pdata_fail = pdata.loc[~passing, :]
    params = dict(linewidths=0.05, edgecolors="k", s=12, alpha=0.6)

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(pdata_pass["n_counts"], pdata_pass["n_genes"],
               color="0.9", label="Pass", **params)
    ax.scatter(pdata_fail["n_counts"], pdata_fail["n_genes"],
               color=sns.xkcd_rgb["light red"], label="Fail", **params)

    ax.set_xlabel("UMIs")
    ax.set_ylabel("Genes")
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
