import os
import re
import typing
import inspect
import numpy as np
import pandas as pd
import pkg_resources
from configparser import ConfigParser, ExtendedInterpolation

import anndata
from scanpy.readwrite import read_10x_h5
from scanpy.readwrite import read as scread, write as scwrite
from scanpy.preprocessing.simple import normalize_per_cell, log1p
import scanpy.api.logging as logg

from .utils import datestamp, timestamp, shift_metadata, silence, quantile_limit
from .qsub import submit_rds_job


class AnalysisConfig(object):
    """
    """
    placeholder_prefix = "EXAMPLE_"
    default_template = pkg_resources.resource_string(__name__, "config_template.ini")
    default_template = default_template.decode("ascii")

    def __init__(self):
        #self.parser = ConfigParser(interpolation=ExtendedInterpolation())
        self.parser = ConfigParser(allow_no_value=True)
        self.parser.optionxform = lambda option: option
        self._load_default_template()

    def _load_default_template(self):
        resource_package = __name__
        self.parser.read_string(self.default_template)
        self.required_sections = self.parser.sections()

    def _validate_string_config(self, input_string):
        config = ConfigParser(allow_no_value=True)
        config.optionxform = lambda option: option
        config.read_string(input_string)
        for section in self.required_sections:
            if section not in config:
                raise Exception(f"Required section {section} not found.")
            if self.placeholder_prefix in section:
                raise Exception(
                    f"Placeholder {placeholder_prefix}* found under section {section}"
                )
            for key, value in config.items():
                assert self.placeholder_prefix not in key
                assert self.placeholder_prefix not in value
        sample_names = set(config["sample_names"].keys())
        for section in self.required_sections:
            if section in ("names", "species"): continue
            assert sample_names - set(config[section].keys()) == set()
        return config


    def print_template(self):
        print("config_string = \"\"\"")
        print(self.default_template)
        print("\"\"\"")
        print("The following sections are required:")
        print(self.required_sections)

    def read(self, input_string):
        self.config = self._validate_string_config(input_string)
        return self.config

    def read_file(self, input_file):
        with open(input_file, 'r') as fin:
            return self.read(fin.read())


def load_10x_data(sample_name: str, config: AnalysisConfig):
    """
    Parameters
    ----------
    sample_name: str
        A sample name present in `config.sample_names`

    config: AnalysisConfig object

    Returns
    -------
    - creates `config["output_dirs"][sample_name]`
    - reads h5 file and sequencing saturation csv (if available)
    - adds the following `uns` metadata:
      - 'sampleid', 'genome', 'species'
      - 'analyst', 'customer_name', 'analysis_version', 'date_created',
        'input_file', 'input_dir', 'output_dir'
    """
    # create output dir
    output_dir = config["output_dirs"][sample_name]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    input_dir = config["input_dirs"][sample_name]
    h5s = ["filtered_gene_bc_matrices_h5.h5",
           "filtered_gene_bc_matrices_mex_h5.h5"]
    h5_file = os.path.join(input_dir, h5s[0])
    if not os.path.exists(h5_file):
        h5_file = os.path.join(input_dir, h5s[1])

    genome = config["genomes"][sample_name]
    adata = read_10x_h5(h5_file, genome)
    adata.var_names_make_unique()
    logg.info("Ran `.var_names_make_unique()` for you.")

    adata.obs['sequencing_saturation'] = np.nan
    seqsat_file = os.path.join(input_dir, "sequencing_saturation.csv")
    if os.path.exists(seqsat_file):
        seqsat = pd.read_csv(seqsat_file, index_col=0, header=0)
        adata.obs.loc[seqsat.index, 'sequencing_saturation'] = seqsat['saturation'].values

    metrics_file = os.path.join(input_dir, "metrics_summary.csv")
    if os.path.exists(metrics_file):
        metrics = pd.read_csv(metrics_file, header=0)
        adata.uns["sequencing_metrics"] = metrics.to_dict(orient="index")[0]

    adata.uns['sampleid'] = sample_name
    adata.uns['genome'] = genome
    adata.uns['species'] = config["species"][genome]

    adata.uns['analyst'] = config["names"]["analyst_name"]
    adata.uns['customer_name'] = config["names"]["customer_name"]
    adata.uns['analysis_version'] = 1
    adata.uns['date_created'] = datestamp()
    adata.uns['input_file'] = os.path.abspath(h5_file)
    adata.uns['input_dir'] = os.path.abspath(input_dir)
    adata.uns['output_dir'] = os.path.abspath(output_dir)

    return adata


def save_anndata(adata, version=None):
    """
    Wrapper around `scanpy.readwrite.write`

    Parameters
    ----------
    adata : AnnData

    version : (default: `None`)
    Adds a new version number to `uns` and uses in the save file name.

    Returns
    -------
    Writes a `h5ad` file and stores the new save file under `uns`.
    Old versions and save files are suffixed under `uns` with `_previous`.
    """
    if version is not None:
        shift_metadata(adata, "version")
        adata.uns["version"] = version

    shift_metadata(adata, "save_file")
    adata.uns["save_file"] = os.path.join(adata.uns["output_dir"],
                                          f"analysis_{adata.uns['version']}_{datestamp()}.h5ad")
    scwrite(adata.uns["save_file"], adata)


def load_anndata(infile):
    """
    Wrapper around `scanpy.readwrite.read`

    Parameters
    ----------
    infile : str

    Returns
    -------
    An `AnnData` object and stores `infile` under `adata.uns["save_file_previous"]`.
    """
    adata = scread(infile)
    adata.uns["save_file_previous"] = infile
    return adata


def save_adata_to_rds(adata, cluster_key="cluster", aggr=False, n_dims=3):
    #tmpd = pd.DataFrame(np.asarray(adata.raw.X.todense()),
    #                    index=adata.obs_names,
    #                    columns=adata.raw.var_names)
    tmpd = anndata.AnnData(np.asarray(adata.raw.X.todense()),
                           var=adata.raw.var, obs=adata.obs)

    # if `.raw` counts aren't normalized, normalize them.
    # if they are, then normalization will return the same values (they're already normalized)
    normalize_per_cell(tmpd)
    log1p(tmpd)
    counts = pd.DataFrame(tmpd.X, columns=tmpd.var.gene_ids.values, index=tmpd.obs_names).T

    features = pd.DataFrame({"Associated.Gene.Name": tmpd.var.index, "Chromosome.Name": 1}, index=counts.index)

    counts.loc["ENSGGENES", :] = quantile_limit(adata.obs, "n_genes")
    counts.loc["ENSGUMI", :] = quantile_limit(adata.obs, "n_counts_total")
    counts.loc["ENSGMITO", :] = quantile_limit(adata.obs, "percent_mito")
    counts.loc["ENSGSEQSAT", :] = quantile_limit(adata.obs, "sequencing_saturation")
    if aggr:
        counts.loc["ENSGSAMP", :] = adata.obs.batch.cat.codes
        features.loc["ENSGSAMP", :] = ["Sample", 1]

    features.loc["ENSGGENES", :] = ["Genes", 1]
    features.loc["ENSGUMI", :] = ["Umi", 1]
    features.loc["ENSGMITO", :] = ["PercentMito", 1]
    features.loc["ENSGSEQSAT", :] = ["Saturation", 1]

    if "qc_tag" in adata.obs_keys():
        counts.loc["ENSGHASHTAG", :] = adata.obs["qc_tag"].str.split("_", expand=True).iloc[:,1].astype(float)
        counts.loc["ENSGHASHTAG", :].fillna(0, inplace=True)
        features.loc["ENSGHASHTAG", :] = ["Tag", 1]

    if n_dims == 3:
        umap = adata.obsm["X_umap_3d"]
    else:
        umap = adata.obsm["X_umap"]
        umap = np.column_stack((umap, np.zeros(umap.shape[0])))
    tsne = pd.DataFrame(umap, columns=["V1", "V2", "V3"], index=adata.obs_names)
    tsne["dbCluster"] = adata.obs[cluster_key].astype(int)

    sampleid = adata.uns["sampleid"]
    outdir = adata.uns["output_dir"]
    for data, out_type in zip((counts, features, tsne),
                              ("counts", "features", "umap3d")):
        outname = f"{sampleid}_{out_type}.csv"
        outfile = os.path.join(outdir, outname)
        data.to_csv(outfile)
        logg.info(f"Saved {out_type} to {outfile}.")

    submit_rds_job(sampleid, outdir, f"{sampleid}_{datestamp()}.Rds")


def export_markers(adata, cluster_key):
    excelname = f"{adata.uns['sampleid']}_markers_{datestamp()}.xlsx"
    excelfile = os.path.join(adata.uns["output_dir"], excelname)

    csvname = f"{adata.uns['sampleid']}_markers_{datestamp()}.csv"
    csvfile = os.path.join(adata.uns["output_dir"], csvname)
    markers = pd.DataFrame.from_records(adata.uns["auroc_markers"])
    markers.to_csv(csvfile)
    logg.info(f"CSV file saved to [{csvfile}].")

    with pd.ExcelWriter(excelfile) as writer:
        for cluster in markers[cluster_key].unique():
            inds = markers[cluster_key] == cluster
            name = f"Cluster {cluster}"
            markers.loc[inds].to_excel(writer, name)
        writer.save()
    logg.info(f"Excel file saved to [{excelfile}].")


def save_adata(obj, suffix):
    outname = f"{obj.uns['sampleid']}-{suffix}_{datestamp()}.h5ad"
    outfile = os.path.join(obj.uns["output_dir"], outname)

    logg.info(f"Saving {outname} to {obj.uns['output_dir']}.")
    scwrite(outfile, obj)

    return outfile


def save_all_adata():
    frame = inspect.currentframe()
    try:
        objects = frame.f_back.f_locals
        adata_objects = filter(lambda x: x.startswith("adata_"), objects.keys())
        for object_name in adata_objects:
            suffix = object_name.split("_")[1]
            obj = objects[object_name]
            save_adata(obj, suffix)
    finally:
        del frame


def save_loom(adata):
    outname = f"{obj.uns['sampleid']}_{datestamp()}.loom"
    outfile = os.path.join(obj.uns["output_dir"], outname)

    logg.info(f"Saving {outname} to {obj.uns['output_dir']}.")
    adata.write_loom(outfile)

    return outfile


__api_objects__ = {
    "AnalysisConfig": AnalysisConfig,
    #"load_anndata": load_anndata,
    #"save_anndata": save_anndata,
    "save_adata": save_adata,
    "save_all_adata": save_all_adata,
    "save_adata_to_rds": save_adata_to_rds,
    "save_loom": save_loom,
    "load_10x_data": load_10x_data,
    "export_markers": export_markers,
}
