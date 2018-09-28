
def preprocess(adata, n_top_genes=1000):
    sc.pp.normalize_per_cell(adata)
    adata.uns['raw_dtype'] = 'normalized count'
    adata.raw = adata
    sc.pp.log1p(adata)

    adata.uns['n_top_genes'] = n_top_genes
    adata_filt = sc.pp.filter_genes_dispersion(adata, n_top_genes=n_top_genes, copy=True)

    return adata, adata_filt


def dimensionality_reduction(adata_filt, n_comps=50, n_neighbors=10, min_dist=0.5, metric="correlation"):
    sc.pp.pca(adata_filt, n_comps=n_comps, svd_solver='arpack')
    sc.pp.neighbors(adata_filt, n_neighbors=n_neighbors, metric=metric)
    adata_filt.uns['neighbors']['params']['metric'] = metric
    sc.tl.umap(adata_filt, min_dist=0.5, n_components=3)


def cluster(adata_filt):
    sc.tl.louvain(adata_filt)
    sc.tl.paga(adata_filt)
    if adata_filt.obs['louvain'].astype(int).min() == 0:
        adata_filt.obs['louvain'] = (adata_filt.obs['louvain'].astype(int) + 1).astype('category')
    _ = scatter(adata_filt, basis='umap', color='louvain', title='', components='1,2', legend_loc='on data')
    _ = scatter(adata_filt, basis='umap', color='louvain', title='', components='1,3', legend_loc='on data')
    _ = scatter(adata_filt, basis='umap', color='louvain', title='', components='2,3', legend_loc='on data')


def find_markers(adata_filt):
    sc.tl.rank_genes_groups(adata_filt, 'louvain', method='wilcoxon', use_raw=False,
                            key_added="rank_genes_groups-wilcoxon")
    sc.tl.rank_genes_groups(adata_filt, 'louvain', method='t-test_overestim_var', use_raw=False,
                            key_added="rank_genes_groups-ttest")
    save_file = os.path.join(OUTPUT_DIRS[adata_filt.uns['sampleid']], 'top5_genes-wilcox.pdf')
    fig = sc.pl.rank_genes_groups_matrixplot(adata_filt, groupby='louvain', n_genes=5,
                                             key='rank_genes_groups-wilcoxon', use_raw=False,
                                             cmap=cmo.cm.tempo,
                                             figsize=(36, 6), save=True)
    os.rename(os.path.join('figures', 'matrixplot.pdf'), save_file)

    save_file = os.path.join(OUTPUT_DIRS[adata_filt.uns['sampleid']], 'top5_genes-ttest.pdf')
    fig = sc.pl.rank_genes_groups_matrixplot(adata_filt, groupby='louvain', n_genes=5,
                                             key='rank_genes_groups-ttest', use_raw=False,
                                             cmap=cmo.cm.tempo,
                                             figsize=(36, 6), save=True)
    os.rename(os.path.join('figures', 'matrixplot.pdf'), save_file)


def print_markers(adata_filt, method, outfile):
    names = adata_filt.uns[f"rank_genes_groups-{method}"]["names"]
    marker_df = pd.DataFrame(names)
    marker_df.index.name = "Rank"
    marker_df.columns.name = "Cluster ID"

    marker_df.to_csv(os.path.join(OUTPUT_DIRS[adata_filt.uns['sampleid']], outfile))
