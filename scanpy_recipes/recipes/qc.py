import pkg_resources
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from scanpy import queries
from scanpy.preprocessing import filter_cells, filter_genes
from scanpy.preprocessing import calculate_qc_metrics


# This is a stopgap while scanpy issues #242
# https://github.com/theislab/scanpy/issues/242
# gets resolved
def read_data_file(var_type, genome):
    return pd.read_csv(
        pkg_resources.resource_filename(
            "scanpy_recipes",
            f"data/{var_type}_genes.{genome}.csv"),
        index_col=1,
        header=0
    ).index


def read_hemoglobin_file(genome):
    return read_data_file("hemoglobin", genome)


def read_mito_file(genome):
    return read_data_file("mito", genome)


def read_ribo_file(genome):
    return read_data_file("ribo", genome)


def read_proliferation_file(genome):
    return read_data_file("proliferation", genome)


def gen_qc(raw_adata):
    """
    Appends calculated metrics, sequencing_saturation,
    percent_mito, hemoglobin_counts, n_counts and n_cells to the raw data.

    Parameters
    ----------
    raw_adata
        AnnData object which stores a data matrix (`adata.X`), dataframe-like annotation
        of observations (`adata.obs`) and variables (`adata.var`) and a dictionary of
        unstructured annotations (`adata.uns`).

    Returns
    -------
    Appends
    - the percentage of mitochondrial genes `percent_mito`,
    - the counts of hemoglobin genes `hemoglobin_counts`,
    - the UMI counts `n_counts`, and
    - the gene counts `n_genes`
    to the observation dataframe (`adata.obs`) of the AnnData object

    Appends
    - the UMI counts `n_counts` and
    - the cell counts `n_cells`
    to the variable dataframe (`adata.var`) of the AnnData object

    Appends plot titles (`obs_titles`) to `adata.uns` for the data in `adata.obs`.
    """

    genome = raw_adata.uns["genome"]
    mt_query = read_mito_file(genome)
    hemo_query = read_hemoglobin_file(genome)

    raw_adata.var["mitochondrial"] = raw_adata.var_names.isin(mt_query)
    raw_adata.var["hemoglobin"] = raw_adata.var_names.isin(hemo_query)

    calculate_qc_metrics(
        raw_adata,
        qc_vars=("mitochondrial", "hemoglobin"),
        inplace=True
    )

    raw_adata.uns["10x_umi_cutoff"] = np.min(raw_adata.obs["total_counts"].astype(int))

    raw_cells, raw_genes = raw_adata.shape
    raw_adata.uns["raw_cells"] = raw_cells
    raw_adata.uns["raw_genes"] = raw_genes
    raw_adata.uns["empty_genes"] = np.sum(raw_adata.var["n_cells_by_counts"] == 0).astype(int)
    raw_adata.uns["10x_metrics"]["important"]["Median Sequencing Saturation per Cell"] = \
        f"{raw_adata.obs['sequencing_saturation'].median():.1f}%"

    raw_adata.uns["obs_titles"] = dict(
        total_counts="UMIs", n_genes_by_counts="Genes",
        pct_counts_mitochondrial="mtRNA content",
        sequencing_saturation="Sequencing saturation",
        total_counts_hemoglobin="Hemoglobin counts"
    )


def run_qc(adata_raw,
           min_cells_per_gene=3,
           min_counts_per_gene=3,
           min_genes_per_cell=200,
           max_genes_per_cell=np.inf,
           min_counts_per_cell=500,
           max_counts_per_cell=np.inf,
           min_sequencing_saturation=-np.inf,
           max_sequencing_saturation=np.inf,
           min_pct_mitochondrial=-np.inf,
           max_pct_mitochondrial=np.inf,
           min_counts_hemoglobin=-np.inf,
           max_counts_hemoglobin=np.inf,
           trial=False):
    """
    Applies quality control thresholds to the raw data (`adata_raw`) using the metrics
    calculated by `gen_qc()`.

    Parameters
    ----------
    adata_raw
        AnnData object
    min_cells_per_gene
        Minimum number of cells expressing a gene in order for the gene to pass QC filter.
    min_counts_per_gene
        Minimum number of UMI counts required for a gene to pass QC filter.
    {min,max}_genes_per_cell
        Minimum number of genes expressed in a cell for the cell to pass QC filter.
    {min,max}_counts_per_cell
        Minimum number of UMI counts required for a cell to pass QC filter.
    {min,max}_sequencing_saturation
        Sets minimum per-cell sequencing saturation percentage for each cell to pass
        filter.  The sequencing saturation is defined as the fraction of confidently
        mapped, valid cell-barcode, valid UMI reads that are non-unique. It is 1-UMI
        counts/total reads.  A 50% sequencing saturation means, total reads = 2 X UMI
        counts; while, a 90% sequencing saturation means, total reads = 10 X UMI counts
    {min,max}_pct_mitochondrial
        Sets maximum fraction of mitochondrial RNA that should be present in a cell to
        pass QC filter.
    {min,max}_counts_hemoglobin
        Sets maximum hemoglobin gene counts for a cell to pass QC filter
    trial
        If `trial == False`, a subset of the raw data with all selected filters applied
        (adata_qc) is returned.  If `trial == True`, the input raw data is retuned with
        additional fields flagging the portion of the raw data that failed QC.

    Returns
    -------
    AnnData object (subject to the value of `trial`)
    """

    keywords = (
        min_cells_per_gene, min_counts_per_gene, min_genes_per_cell, max_genes_per_cell,
        min_counts_per_cell, max_counts_per_cell, min_sequencing_saturation,
        max_sequencing_saturation, min_pct_mitochondrial, max_pct_mitochondrial,
        min_counts_hemoglobin, max_counts_hemoglobin
    )

    if None in keywords:
        raise ValueError(
            "Cannot use `None` as a keyword.  Please use `-np.inf` to represent lower "
            "bounds and `np.inf` to represent upper bounds."
        )

    orig_shape = adata_raw.shape
    adata = adata_raw.copy()
    adata_qc = adata_raw.copy()

    count_subset_min = np.ones_like(adata_qc.obs_names, dtype=bool)
    count_subset_max = np.ones_like(adata_qc.obs_names, dtype=bool)
    gene_subset_min = np.ones_like(adata_qc.obs_names, dtype=bool)
    gene_subset_max = np.ones_like(adata_qc.obs_names, dtype=bool)
    if min_counts_per_cell:
        count_subset_min, n_counts_min = filter_cells(adata_qc.X, min_counts=min_counts_per_cell)
    if max_counts_per_cell:
        count_subset_max, n_counts_max = filter_cells(adata_qc.X, max_counts=max_counts_per_cell)
    if min_genes_per_cell:
        gene_subset_min, n_genes_min = filter_cells(adata_qc.X, min_genes=min_genes_per_cell)
    if max_genes_per_cell:
        gene_subset_max, n_genes_max = filter_cells(adata_qc.X, max_genes=max_genes_per_cell)
    adata.obs["qc_fail_counts"] = ~(count_subset_min & count_subset_max)
    adata.obs["qc_fail_genes"] = ~(gene_subset_min & gene_subset_max)

    seqsat_subset, mito_subset, rbc_subset = True, True, True
    if min_sequencing_saturation:
        seqsat_subset &= adata_qc.obs["sequencing_saturation"].fillna(0) > min_sequencing_saturation
    if max_sequencing_saturation:
        seqsat_subset &= adata_qc.obs["sequencing_saturation"].fillna(0) < max_sequencing_saturation
    if min_pct_mitochondrial:
        mito_subset &= adata_qc.obs["pct_counts_mitochondrial"].fillna(0) > min_pct_mitochondrial
    if max_pct_mitochondrial:
        mito_subset &= adata_qc.obs["pct_counts_mitochondrial"].fillna(0) < max_pct_mitochondrial
    if min_counts_hemoglobin:
        rbc_subset &= adata_qc.obs["total_counts_hemoglobin"].fillna(0) > min_counts_hemoglobin
    if max_counts_hemoglobin:
        rbc_subset &= adata_qc.obs["total_counts_hemoglobin"].fillna(0) < max_counts_hemoglobin
    keep_subset = count_subset_min & count_subset_max
    keep_subset &= gene_subset_min & gene_subset_max
    keep_subset &= seqsat_subset & mito_subset & rbc_subset
    adata_qc._inplace_subset_obs(keep_subset)

    adata.obs["qc_fail_seqsat"] = ~seqsat_subset
    adata.obs["qc_fail_mito"] = ~mito_subset
    adata.obs["qc_fail_rbc"] = ~rbc_subset
    adata.obs["qc_fail"] = "fail"
    adata.obs.loc[keep_subset, "qc_fail"] = "pass"
    adata.obs["qc_fail"] = adata.obs["qc_fail"].astype("category")

    filter_genes(adata_qc, min_cells=min_cells_per_gene)
    filter_genes(adata_qc, min_counts=min_counts_per_gene)
    adata.var["qc_fail_counts"] = ~adata.var_names.isin(adata_qc.var_names)

    n_rbcs = 0 if isinstance(rbc_subset, bool) else int(sum(~rbc_subset))
    adata.uns["qc_gene_filter"] = {
        "threshold_n_cells": min_cells_per_gene,
        "threshold_n_counts": min_counts_per_gene
    }
    adata.uns["qc_cell_filter"] = {
        "threshold_n_genes_by_counts": (min_genes_per_cell, max_genes_per_cell),
        "threshold_total_counts": (min_counts_per_cell, max_counts_per_cell),
        "threshold_sequencing_saturation": (min_sequencing_saturation,
                                            max_sequencing_saturation),
        "threshold_pct_counts_mitochondrial": (min_pct_mitochondrial,
                                               max_pct_mitochondrial),
        "threshold_total_counts_hemoglobin": (min_counts_hemoglobin,
                                              max_counts_hemoglobin),
    }

    qc_metrics = {
        "red_blood_cells_removed": n_rbcs,
        "cells_removed": int(sum(adata.obs.qc_fail == "fail")),
        "high_mtrna_cells_removed": int(sum(~mito_subset & rbc_subset)),
        "low_sequencing_saturation_cells_removed": \
                int(sum(~seqsat_subset & rbc_subset & mito_subset)),
        "low_count_cells_removed": int(sum(adata.obs.qc_fail_counts)),
        "low_gene_cells_removed": int(sum(adata.obs.qc_fail_genes)),
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
