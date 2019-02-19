import os
import re
import glob
import typing
import pathlib
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

from .versions import REPORT_SCHEMA_VERSION, ANALYSIS_PIPELINE_VERSION
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
        sample_ids = set(config["sample_names"].keys())
        for section in self.required_sections:
            if section in ("names", "species"): continue
            assert sample_ids - set(config[section].keys()) == set()
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


def parse_10x_metrics(metrics_file):
    if not os.path.exists(metrics_file):
        return dict()

    metrics = pd.read_csv(metrics_file, header=0).T
    alignment = metrics.index.str.startswith("Reads Mapped")
    sequencing = metrics.index.str.startswith("Q30")
    important = metrics.index.str.contains("per Cell") | \
                metrics.index.str.contains("of Cells")
    sample_level = ~alignment & ~sequencing & ~important

    return {
        "alignment": metrics[alignment].to_dict()[0],
        "sequencing": metrics[sequencing].to_dict()[0],
        "important": metrics[important].to_dict()[0],
        "sample": metrics[sample_level].to_dict()[0],
    }


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
    for h5_file in glob.glob(f"{input_dir}/filtered_*matri*.h5"):
        if os.path.exists(h5_file):
            break
    else:
        raise IOError(f"Can't find filtered matrix h5 file under [{input_dir}].")

    genome = config["genomes"][sample_name]
    adata = read_10x_h5(h5_file, genome)
    adata.var_names_make_unique()
    logg.info("Ran `.var_names_make_unique()` for you.")

    adata.obs["sequencing_saturation"] = np.nan
    seqsat_file = os.path.join(input_dir, "sequencing_saturation.csv")
    if os.path.exists(seqsat_file):
        seqsat = pd.read_csv(seqsat_file, index_col=0, header=0)
        adata.obs.loc[seqsat.index, "sequencing_saturation"] = seqsat["saturation"].values

    metrics_file = os.path.join(input_dir, "metrics_summary.csv")
    adata.uns["10x_metrics"] = parse_10x_metrics(metrics_file)
    adata.uns["10x_metrics"]["target_cells"] = config["target_cells"].get(sample_name, None)
    adata.uns["10x_chemistry"] = config["10x_chemistry"].get(sample_name, None)
    adata.uns["cellranger_version"] = config["cellranger_version"].get(sample_name, None)
    adata.uns["cellranger_reference_version"] = config["reference_version"].get(sample_name, None)

    adata.uns["sampleid"] = sample_name
    customer_sample_name = config["sample_names"].get(sample_name, None)
    if customer_sample_name:
        adata.uns["sample_name"] = customer_sample_name
    adata.uns["genome"] = genome
    adata.uns["species"] = config["species"][genome]

    adata.uns["analyst"] = config["names"]["analyst_name"]
    adata.uns["customer_name"] = config["names"]["customer_name"]
    adata.uns["principal_investigator_name"] = config["names"].get(
        "pi_name", adata.uns["customer_name"]
    )
    adata.uns["date_created"] = datestamp()

    adata.uns["input_file"] = os.path.abspath(h5_file)
    adata.uns["input_dir"] = os.path.abspath(input_dir)
    adata.uns["output_dir"] = os.path.abspath(output_dir)

    adata.uns["report_schema_version"] = REPORT_SCHEMA_VERSION
    adata.uns["analysis_pipeline_version"] = ANALYSIS_PIPELINE_VERSION

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


def get_rds_embedding(adata, cluster_key="cluster", n_dims=3, embedding_type="umap"):
    if n_dims == 3:
        embedding = adata.obsm[f"X_{embedding_type}_3d"]
    else:
        embedding = adata.obsm[f"X_{embedding_type}"]
        embedding = np.column_stack((umap, np.zeros(umap.shape[0])))
    tsne_df = pd.DataFrame(embedding, columns=["V1", "V2", "V3"], index=adata.obs_names)
    tsne_df["dbCluster"] = adata.obs[cluster_key].astype(int)

    return tsne_df


def save_rds_embedding_info(adata, outpath, cluster_key="cluster", n_dims=3):
    df = get_rds_embedding(adata, cluster_key=cluster_key, n_dims=n_dims)
    df.to_csv(outpath, sep=",")


def write_csvs(adata, outdir, uns_keys=None):
    """
    More controlled version of `anndata.AnnData.write_csvs`

    Parameters
    ----------
    adata : AnnData object
    outdir : custom outdir to write to
    uns_keys : iterable containing `.uns` keys to export.  If None, it'll take a good
        guess.

    Returns
    -------
    """
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if uns_keys is None:
        uns_keys = [
            "sampleid", "sample_name", "principal_investigator_name", "year",
            "submitter_name", "analyst", "customer_name", "gt_project_id",
            "genome", "species", "strain", "tissue_type", "cell_line",
            "chemistry", "target_cells",
            "raw_genes", "raw_cells", "n_top_genes"
        ]

    uns_df = {}
    for key in uns_keys:
        val = adata.uns.get(key, None)
        if val is not None:
            uns_df[key] = val

    metrics_10x = adata.uns.get("10x_metrics", None)
    if metrics_10x is not None:
        metrics = {}
        for k, v in metrics_10x.items():
            metrics.update(**v)
        uns_df.update(**metrics)
    uns_df = pd.DataFrame(uns_df, index=["value"]).T

    to_write = {
        "X": pd.DataFrame(adata._X.toarray() if issparse(adata._X) else adata._X),
        "var": adata._var.to_df(),
        "obs": adata._obs.to_df(),
        "varm": adata._varm.to_df(),
        "obsm": adata._obsm.to_df(),
        "uns": uns_df
    }
    if adata.raw is not None:
        to_write["X_raw"] = pd.DataFrame(adata.raw.X.toarray()
                                         if issparse(adata.raw.X) else adata.raw.X)

    for key, value in to_write.items():
        filename = os.path.join(outdir, f"{key}.csv")
        df.to_csv(
            filename, sep=",",
            header=key in {"obs", "var", "obsm", "varm"},
            index=key in {"obs", "var"},
        )


def save_adata_to_rds(adata, cluster_key="cluster", n_dims=3):
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
    if adata.uns.get("is_aggregation", None):
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

    tsne = get_rds_embedding(adata, cluster_key, n_dims=n_dims)

    sampleid = adata.uns["sampleid"]
    outdir = adata.uns["output_dir"]
    for data, out_type in zip((counts, features, tsne),
                              ("counts", "features", "umap3d")):
        outname = f"{sampleid}_{out_type}.csv"
        outfile = os.path.join(outdir, outname)
        data.to_csv(outfile)
        logg.info(f"Saved {out_type} to {outfile}.")

    submit_rds_job(sampleid, outdir, f"{sampleid}_{datestamp()}.Rds")


def get_marker_dataframe(adata, marker_key="auroc_markers"):
    return pd.DataFrame.from_dict(adata.uns[marker_key])


def export_markers(adata, cluster_key, marker_key="auroc_markers"):
    excelname = f"{adata.uns['sampleid']}_markers_{datestamp()}.xlsx"
    excelfile = os.path.join(adata.uns["output_dir"], excelname)

    csvname = f"{adata.uns['sampleid']}_markers_{datestamp()}.csv"
    csvfile = os.path.join(adata.uns["output_dir"], csvname)

    #markers = pd.DataFrame.from_records(adata.uns["auroc_markers"])
    markers = get_marker_dataframe(adata, marker_key)
    markers[cluster_key] = markers[cluster_key].astype(int)

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
    "get_marker_dataframe": get_marker_dataframe,
    "write_csvs": write_csvs,
    "save_rds_embedding_info": save_rds_embedding_info,
}
