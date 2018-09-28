import matplotlib.pyplot as plt
import seaborn as sns
import cmocean as cmo


def plot_qc(adata, config):
    fig, axarr = plt.subplots(1, 2, dpi=300, figsize=(10, 4))
    fig.patch.set_facecolor('white')
    kwargs = dict(s=16, lw=0.25, edgecolors='none')

    obs = adata.obs.copy()
    obs = obs.sort_values('sequencing_saturation', ascending=False)
    p1 = axarr[0].scatter(obs['n_counts'], obs['n_genes'], c=obs['sequencing_saturation'],
                          cmap=cmo.cm.matter, **kwargs)
    axarr[0].set_xlabel('UMIs')
    axarr[0].set_ylabel('Genes')
    fig.colorbar(p1, ax=axarr[0], label='Sequencing Saturation')

    obs = obs.sort_values('percent_mito', ascending=True)
    p2 = axarr[1].scatter(obs['n_counts'], obs['n_genes'], c=obs['percent_mito'],
                       cmap=cmo.cm.matter, **kwargs)
    axarr[1].set_xlabel('UMIs')
    axarr[1].set_ylabel('Genes')
    fig.colorbar(p2, ax=axarr[1], label='mtRNA content')

    [sns.despine(fig=fig, ax=ax) for ax in axarr.flat]
    fig.tight_layout()
    fig.savefig(os.path.join(config.output_dirs[adata.uns['sampleid']], 'qc.pdf'),
                format='pdf', dpi=300)


def plot_scatters(adata, basis='tsne', raw=True, components='1,2', **kwargs):
    sampleid = adata.uns['sampleid']
    sc.pl.scatter(adata, basis=basis, color='louvain', projection='3d',
                  size=4, title='Clusters', save="_clusters_3d.pdf", **kwargs)
    os.rename('figures/umap_clusters_3d.pdf',
              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_clusters_3d.pdf"))
    sc.pl.scatter(adata, basis=basis, color='louvain',
                  size=4, title='Clusters', save="_clusters-labeled.pdf",
                  components=adata.uns['umap_view_components'], legend_loc='on data', **kwargs)
    os.rename('figures/umap_clusters-labeled.pdf',
              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_clusters.pdf"))
    sc.pl.scatter(adata, basis=basis, color='n_counts', color_map=cmo.cm.dense,
                  size=4, title='UMI counts', save="_umi_counts.pdf",
                  components=adata.uns['umap_view_components'], **kwargs)
    os.rename('figures/umap_umi_counts.pdf',
              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_umi_counts.pdf"))
    sc.pl.scatter(adata, basis=basis, color='n_genes', color_map=cmo.cm.dense,
                  size=4, title='Gene counts', save="_gene_counts.pdf",
                  components=adata.uns['umap_view_components'], **kwargs)
    os.rename('figures/umap_gene_counts.pdf',
              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_gene_counts.pdf"))
    sc.pl.scatter(adata, basis=basis, color='percent_mito', color_map=cmo.cm.matter,
                  size=4, title='mtRNA fraction', save="_mtrna_fraction.pdf",
                  components=adata.uns['umap_view_components'], **kwargs)
    os.rename('figures/umap_mtrna_fraction.pdf',
              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_mtrna_fraction.pdf"))
    sc.pl.scatter(adata, basis=basis, color='sequencing_saturation', color_map=cmo.cm.matter,
                  size=4, title='Sequencing Saturation', save="_sequencing_saturation.pdf",
                  components=adata.uns['umap_view_components'], **kwargs)
    os.rename('figures/umap_sequencing_saturation.pdf',
              os.path.join(OUTPUT_DIRS[sampleid], f"{sampleid}_sequencing_saturation.pdf"))

    cell_cycle_genes = pd.read_csv('/home/flynnb/projects/singlecell/cell_cycle_gene_sets.csv', header=0)
    cell_cycle_genes = cell_cycle_genes.apply(lambda x: x.str.strip())
    n_groups = cell_cycle_genes.shape[1]
    for k in range(n_groups):
        group = cell_cycle_genes.columns[k]
        group_genes = adata.var_names[adata.var_names.isin(cell_cycle_genes.iloc[:, k])]
        key = f'cell_cycle_{group}'.replace('/', '-')
        adata.obs[key] = adata[:, group_genes].X.sum() / adata.obs['n_counts']
        if not (adata.obs[key] > 0).any():
            adata.obs.ix[-1, key] = 1e-6

        sc.pl.scatter(adata, basis=basis, color=key, color_map=cmo.cm.tempo, legend_fontsize=0,
                      size=4, title=group, save=f"_{key}.pdf",
                      components=adata.uns['umap_view_components'], **kwargs)
        os.rename(f'figures/umap_{key}.pdf',
                  os.path.join(OUTPUT_DIRS[sampleid], f'{sampleid}_{key}.pdf'))


