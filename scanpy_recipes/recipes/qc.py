import numpy as np
from collections import OrderedDict
from scanpy import queries
from scanpy.preprocessing.simple import filter_cells, filter_genes

def gen_qc(raw_adata):
    if raw_adata.uns['genome'].lower().startswith('m'):
        mt_query = queries.mitochondrial_genes('useast.ensembl.org', 'mmusculus')
    else:
        mt_query = queries.mitochondrial_genes('useast.ensembl.org', 'hsapiens')

    mito_genes = raw_adata.var_names.intersection(mt_query)
    raw_adata.obs['percent_mito'] = np.sum(
        raw_adata[:, mito_genes].X, axis=1).A1 / np.sum(raw_adata.X, axis=1).A1 * 100

    raw_adata.obs['n_counts'] = np.sum(raw_adata.X, axis=1).A1
    raw_adata.obs['n_genes'] = np.sum(raw_adata.X > 0, axis=1).A1
    raw_adata.var['n_counts'] = np.sum(raw_adata.X, axis=0).A1
    raw_adata.var['n_cells'] = np.sum(raw_adata.X > 0, axis=0).A1

    raw_cells, raw_genes = raw_adata.shape
    raw_adata.uns['raw_cells'] = raw_cells
    raw_adata.uns['raw_genes'] = raw_genes
    raw_adata.uns['empty_genes'] = np.sum(raw_adata.var['n_cells'] == 0).astype(int)


def run_qc(adata, min_cells=3, min_genes=200, min_counts_per_cell=500, min_counts_per_gene=3,
           sequencing_saturation=50.0, percent_mito=50.0):

    orig_shape = adata.shape
    filter_genes(adata, min_cells=min_cells)
    filter_genes(adata, min_counts=min_counts_per_gene)
    filter_cells(adata, min_counts=min_counts_per_cell)
    filter_cells(adata, min_genes=min_genes)
    #redux = adata.copy()
    redux = adata[adata.obs['sequencing_saturation'] > sequencing_saturation]
    redux = redux[redux.obs['percent_mito'] < percent_mito]
    qc_shape = redux.shape
    print(orig_shape, qc_shape)

    #redux.uns = adata.uns.copy()
    redux.uns['qc_gene_filter'] = {'min_cells': min_cells, 'min_counts': min_counts_per_gene}
    redux.uns['qc_cell_filter'] = {'min_genes': min_genes, 'min_counts': min_counts_per_cell,
                                   'min_sequencing_saturation': sequencing_saturation,
                                   'max_percent_mito': percent_mito}

    #return redux
