from matplotlib.testing import setup
setup()
import pytest

from scanpy_recipes import sc

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#import pkg_resources

from .prep_test_dataset import TEST_DATASET_LOCATION as INPUT_DIR
from .prep_test_dataset import download_10x_data

SAMPLE_NAME = "XX12345"
#INPUT_DIR = pkg_resources.resource_filename("scanpy_recipes", f"tests/datasets/{SAMPLE_NAME}")
OUTPUT_DIR = "/tmp/scanpy-tests/test-outputs"

config_string = f"""
[names]
customer_name = Anonymous person
analyst_name = Analyst name
analysis_name = Test-analysis

[sample_names]
{SAMPLE_NAME} = 10X_1k_PMBCs

[genomes]
{SAMPLE_NAME} = GRCh38

[species]
hg19 = hsapiens
GRCh38 = hsapiens
mm10 = mmusculus

[input_dirs]
{SAMPLE_NAME} = {INPUT_DIR}

[output_dirs]
{SAMPLE_NAME} = {OUTPUT_DIR}

[sample_info]
target_cells = 1000
10x_chemistry = v3
"""


def setup():
    download_10x_data()


def test_base():
    ac = sc.AnalysisConfig()
    config = ac.read(config_string)

    adata_raw = sc.load_10x_data(SAMPLE_NAME, config)

    sc.qc.gen_qc(adata_raw)
    qc_params = dict(
        min_cells_per_gene=3,
        min_counts_per_gene=3,
        min_counts_per_cell=1000,
        min_genes_per_cell=500,
        sequencing_saturation=None,
        percent_mito=20.0,
        rbc_threshold=10
    )
    adata_qc = sc.qc.run_qc(adata_raw, trial=False, **qc_params)

    adata = sc.pp.preprocess(adata_qc, n_top_genes=1000, scale=True)
    sc.pp.dimensionality_reduction(adata, n_neighbors=10, min_dist=0.5)

    sc.tl.cluster(adata)

    markers = sc.tl.find_marker_genes(adata, "cluster", log_fold_change=1.0)
    sc.export_markers(adata, "cluster")

    sc.save_adata_to_rds(adata, cluster_key="cluster", submit=False)

    redux_save_file = sc.save_adata(adata, "test")


@pytest.mark.skip(reason="figure handling is broken")
def test_report():
    ac = sc.AnalysisConfig()
    config = ac.read(config_string)
    adata_raw = sc.load_10x_data(SAMPLE_NAME, config)

    sc.qc.gen_qc(adata_raw)
    qc_params = dict(
        min_cells_per_gene=3,
        min_counts_per_gene=3,
        min_counts_per_cell=1000,
        min_genes_per_cell=500,
        sequencing_saturation=None,
        percent_mito=20.0,
        rbc_threshold=10
    )
    trial = sc.qc.run_qc(adata_raw, trial=True, **qc_params)

    qc_fig1 = sc.pl.umi_rank_plot(adata_raw, return_fig=True)
    qc_fig1 = sc.pl.qc_violins(trial, return_fig=True)
    qc_fig3 = sc.pl.qc_pass_fail(trial, return_fig=True)

    adata_qc = sc.qc.run_qc(adata_raw, trial=False, **qc_params)

    adata = sc.pp.preprocess(adata_qc, n_top_genes=1000, scale=True)
    sc.pp.dimensionality_reduction(adata, n_neighbors=10, min_dist=0.5)

    sc.tl.cluster(adata)

    markers = sc.tl.find_marker_genes(adata, "cluster", log_fold_change=1.0)
    sc.export_markers(adata, "cluster")

    sc.save_adata_to_rds(adata, cluster_key="cluster", submit=False)

    redux_save_file = sc.save_adata(adata, "test")

    report = sc.SCBLReport()
    report.add_report_figures(
        adata_redux,
        cluster_key=cluster_key,
        ranks=[qc_fig1],
        violins=qc_fig2,
        scatters=qc_fig3,
    )
    report.generate_report(adata_redux)
