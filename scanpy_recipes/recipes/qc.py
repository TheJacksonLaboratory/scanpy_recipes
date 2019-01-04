import pkg_resources
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from scanpy import queries
from scanpy.preprocessing.simple import filter_cells, filter_genes


def read_mito_file(genome):
    return pd.read_csv(pkg_resources.resource_filename("scanpy_recipes",
                                                       f"data/mito_genes.{genome}.csv"),
                       index_col=1, header=0).index

def gen_qc(raw_adata):
    """
    Appends new calculated metrics, sequencing_saturation, 
    percent_mito, hemoglobin_counts, n_counts and n_cells to the raw data. 
    
    Parameters
    ----------
    raw_adata
           Annotated data matrix.

    Returns
    -------
    None

    """
    genome = raw_adata.uns["genome"]
    mt_query = read_mito_file(genome)
    if genome.lower().startswith("m"):
        #mt_query = queries.mitochondrial_genes("useast.ensembl.org", "mmusculus")
        rbc_gene = "Hbb-bs"
    else:
        #mt_query = queries.mitochondrial_genes("useast.ensembl.org", "hsapiens")
        rbc_gene = "HBB"

    mito_genes = raw_adata.var_names.intersection(mt_query)
    raw_adata.obs["percent_mito"] = np.sum(
        raw_adata[:, mito_genes].X, axis=1).A1 / np.sum(raw_adata.X, axis=1).A1 * 100

    raw_adata.obs["hemoglobin_counts"] = raw_adata[:, rbc_gene].X

    raw_adata.obs["n_counts"] = np.sum(raw_adata.X, axis=1).A1
    raw_adata.obs["n_genes"] = np.sum(raw_adata.X > 0, axis=1).A1
    raw_adata.var["n_counts"] = np.sum(raw_adata.X, axis=0).A1
    raw_adata.var["n_cells"] = np.sum(raw_adata.X > 0, axis=0).A1
    raw_adata.uns["10x_umi_cutoff"] = np.min(raw_adata.obs["n_counts"].astype(int))

    raw_cells, raw_genes = raw_adata.shape
    raw_adata.uns["raw_cells"] = raw_cells
    raw_adata.uns["raw_genes"] = raw_genes
    raw_adata.uns["empty_genes"] = np.sum(raw_adata.var["n_cells"] == 0).astype(int)
    raw_adata.uns["10x_metrics"]["important"]["Median Sequencing Saturation per Cell"] = \
        f"{raw_adata.obs['sequencing_saturation'].median():.1f}%"

    raw_adata.uns["obs_titles"] = dict(
        n_counts="UMIs", n_genes="Genes", percent_mito="mtRNA content",
        sequencing_saturation="Sequencing saturation",
        hemoglobin_counts="Hemoglobin counts"
    )


def run_qc(adata_raw,
           min_cells_per_gene=3,
           min_counts_per_gene=3,
           min_genes_per_cell=200,
           min_counts_per_cell=500,
           sequencing_saturation=50.0,
           percent_mito=50.0,
           rbc_threshold=10,
           trial=False):
    """

    """
    orig_shape = adata_raw.shape
    adata = adata_raw.copy()
    adata_qc = adata_raw.copy()

    filter_genes(adata_qc, min_cells=min_cells_per_gene)
    filter_genes(adata_qc, min_counts=min_counts_per_gene)
    adata.var["qc_fail_counts"] = ~adata.var_names.isin(adata_qc.var_names)
    count_subset, n_counts = filter_cells(adata_qc.X, min_counts=min_counts_per_cell)
    gene_subset, n_genes = filter_cells(adata_qc.X, min_genes=min_genes_per_cell)
    adata.obs["qc_fail_counts"] = ~count_subset
    adata.obs["qc_fail_genes"] = ~gene_subset

    seqsat_subset, mito_subset, rbc_subset = True, True, True
    if sequencing_saturation:
        seqsat_subset = adata_qc.obs["sequencing_saturation"] > sequencing_saturation
    if percent_mito:
        mito_subset = adata_qc.obs["percent_mito"] < percent_mito
    if rbc_threshold:
        rbc_subset = adata_qc.obs["hemoglobin_counts"] < rbc_threshold
    keep_subset = seqsat_subset & mito_subset & rbc_subset & gene_subset & count_subset
    adata_qc._inplace_subset_obs(keep_subset)

    adata.obs["qc_fail_seqsat"] = ~seqsat_subset
    adata.obs["qc_fail_mito"] = ~mito_subset
    adata.obs["qc_fail_rbc"] = ~rbc_subset
    adata.obs["qc_fail"] = "fail"
    adata.obs.loc[keep_subset, "qc_fail"] = "pass"
    adata.obs["qc_fail"] = adata.obs["qc_fail"].astype("category")

    n_rbcs = int(sum(~rbc_subset))
    adata.uns["qc_gene_filter"] = {
        "threshold_n_cells": min_cells_per_gene,
        "threshold_n_counts": min_counts_per_gene
    }
    adata.uns["qc_cell_filter"] = {
        "threshold_n_genes": min_genes_per_cell,
        "threshold_n_counts": min_counts_per_cell,
        "threshold_sequencing_saturation": sequencing_saturation,
        "threshold_percent_mito": percent_mito,
        "threshold_hemoglobin_counts": rbc_threshold
    }

    qc_metrics = {
        "red_blood_cells_removed": n_rbcs,
        "cells_removed": int(sum(adata.obs.qc_fail == "fail")),
        "high_mtrna_cells_removed": int(sum(~mito_subset & rbc_subset)),
        "low_sequencing_saturation_cells_removed": \
                int(sum(~seqsat_subset & rbc_subset & mito_subset)),
        "low_count_cells_removed": int(sum(adata.obs.qc_fail_counts)),
        "low_count_genes_removed": int(sum(adata.var.qc_fail_counts)),
        "cells": adata_qc.shape[0],
        "genes": adata_qc.shape[1]
    }
    adata.uns["qc_metrics"] = qc_metrics

    qc_shape = adata_qc.shape
    print("Original dims: {}\nFiltered dims: {}".format(orig_shape, qc_shape))

    if trial:
        return adata

    adata_qc.uns = adata.uns.copy()
    return adata_qc


__api_objects__ = {
    "gen_qc": gen_qc,
    "run_qc": run_qc,
}
