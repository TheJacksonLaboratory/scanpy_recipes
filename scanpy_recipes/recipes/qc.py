import numpy as np
import matplotlib.pyplot as plt
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

    redux = adata.copy()
    redux = adata[adata.obs['sequencing_saturation'] > sequencing_saturation]
    redux = redux[redux.obs['percent_mito'] < percent_mito]
    qc_shape = redux.shape
    print("Original dims: {}\nFiltered dims: {}".format(orig_shape, qc_shape))

    redux.uns = adata.uns.copy()
    redux.uns['qc_gene_filter'] = {'min_cells': min_cells, 'min_counts': min_counts_per_gene}
    redux.uns['qc_cell_filter'] = {'min_genes': min_genes, 'min_counts': min_counts_per_cell,
                                   'min_sequencing_saturation': sequencing_saturation,
                                   'max_percent_mito': percent_mito}

    return redux


def detect_rbc(adata, threshold=10):
    if adata.uns["species"] == 'mmusculus':
        key_gene = 'Hbb-bs'
    else:
        key_gene = 'HBB'

    cell_subset = adata[:, key_gene].X < threshold
    n_rbcs = adata.shape[0] - sum(cell_subset)

    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(8, 3))
    fig.patch.set_facecolor("white")
    p1 = ax1.scatter(adata.obs['n_counts'], adata.obs['n_genes'], c=adata[:, key_gene].X, cmap='Reds', s=12, alpha=0.8)
    ax1.set_xlabel('UMIs')
    ax1.set_ylabel('Genes')
    ax1.set_title('Before filtering')
    fig.colorbar(p1, ax=ax1, label=key_gene)


    adata._inplace_subset_obs(cell_subset)
    adata.uns['red_blood_cells_removed'] = n_rbcs
    p2 = ax2.scatter(adata.obs['n_counts'], adata.obs['n_genes'], c=adata[:, key_gene].X, cmap='Reds', s=12, alpha=0.8)
    ax2.set_xlabel('UMIs')
    ax2.set_ylabel('Genes')
    ax2.set_title('After filtering')
    fig.colorbar(p2, ax=ax2, label=key_gene)
    fig.tight_layout()

    print(f"Removed {n_rbcs} red blood cells.")


__api_objects__ = {
    "gen_qc": gen_qc,
    "run_qc": run_qc,
}
