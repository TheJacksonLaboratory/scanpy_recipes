import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo
import pandas as pd
from scanpy.api import pl


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


def cut_violin(ax, threshold=0., cut_above=[True, False], colors=["0.9", "r"]):
    def cut_verts(polycol, threshold, cut_above, color):
        v = polycol.get_paths()[0].vertices
        if cut_above:
            ind = v[:,1] <= threshold
        else:
            ind = v[:,1] > threshold
        polycol.set_verts([v[ind]])
        polycol.set(facecolor=color)
    cut_verts(ax.collections[0], threshold, cut_above=cut_above[0], color=colors[0])
    cut_verts(ax.collections[2], threshold, cut_above=cut_above[1], color=colors[1])


def qc_violins(adata, return_fig=False):
    keys = sorted(filter(lambda s: not s.startswith("qc"), adata.obs_keys()))
    N = len(keys)

    thresholds = adata.uns.get("qc_cell_filter", None)
    use_thresholds = dict(zip(keys, [False]*N))
    if thresholds:
        use_thresholds = dict((key, thresholds.get("threshold_" + key, None)) for key in keys)

    fig, axs = plt.subplots(1, N, figsize=(4*N, N), dpi=200)
    for ax, key in zip(axs.flat, keys):
        threshold = use_thresholds[key]
        print(key, threshold)
        if threshold:
            ax = pl.violin(adata, key, color="blue", show=False, ax=ax, cut=0)
            ax = pl.violin(adata, key, color="blue", show=False, ax=ax, cut=0, jitter=False)
            cut_violin(ax, threshold=threshold,
                       cut_above=[True, False] if key == "percent_mito" else [False, True])
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


def genes_umis_scatter(adata_trial, return_fig=False):
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


__api_objects__ = {
    "qc_violins": qc_violins,
    "genes_umis_scatter": genes_umis_scatter,
}
