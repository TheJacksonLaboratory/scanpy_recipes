import os
import re
import anndata
import typing
from scanpy.readwrite import read_10x_h5
from scanpy.readwrite import read as scread, write as scwrite
from datetime import datetime
from argparse import Namespace
from configparser import ConfigParser, ExtendedInterpolation
import pkg_resources



class AnalysisConfig(object):
    """
    """
    placeholder_prefix = "EXAMPLE_"
    default_template = pkg_resources.resource_string(__name__, "config_template.ini")
    default_template = default_template.decode("ascii")

    def __init__(self):
        #self.parser = ConfigParser(interpolation=ExtendedInterpolation())
        self.parser = ConfigParser(allow_no_value=True)
        self._load_default_template()

    def _load_default_template(self):
        resource_package = __name__
        self.parser.read_string(self.default_template)
        self.required_sections = self.parser.sections()

    def _validate_string_config(self, input_string):
        config = ConfigParser(allow_no_value=True)
        config.parse_string(input_string)
        for section in self.required_sections:
            if section not in config:
                raise Exception(f"Required section {section} not found.")
            if self.placeholder_prefix in section:
                raise Exception(
                    f"Placeholder {placeholder_prefix}* found under section {section}"
                )
            for key, value in config.items("section"):
                assert self.placeholder_prefix not in key
                assert self.placeholder_prefix not in value
        sample_names = set(config["sample_names"].keys())
        for section in self.required_sections:
            if section in ("names", "species"): continue
            assert sample_names - set(config[section].keys()) == set()
        return config


    def print_template(self):
        print("config = \"\"\"")
        print(self.default_template)
        print("\"\"\"")
        print("The following sections are required:")
        print(self.required_sections)

    def read(self, input_string):
        self.config = self._validate_string_template(input_string)
        return self.config

    def read_file(self, input_file):
        with open(input_file, 'r') as fin:
            return self.read(fin.read())


def timestamp():
    """
    Returns current time in format: YYYY-MM-DDThh-mm-ss
    """
    return datetime.now().strftime('%Y-%m-%dT%H-%M-%S')


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

    adata.obs['sequencing_saturation'] = np.nan
    seqsat = pd.read_csv(os.path.join(input_dir,
                                      "sequencing_saturation.csv"),
                         index_col=0, header=0)
    adata.obs.loc[seqsat.index, 'sequencing_saturation'] = seqsat['saturation'].values

    adata.uns['sampleid'] = sample_name
    adata.uns['genome'] = genome
    adata.uns['species'] = config["species"][genome]

    adata.uns['analyst'] = config["names"]["analyst_name"]
    adata.uns['customer_name'] = config["names"]["customer_name"]
    adata.uns['analysis_version'] = 1
    adata.uns['date_created'] = timestamp()
    adata.uns['input_file'] = os.path.abspath(h5_file)
    adata.uns['input_dir'] = os.path.abspath(input_dir)
    adata.uns['output_dir'] = os.path.abspath(output_dir)

    return adata


def shift_metadata(adata, uns_key):
    prev_key = adata.get(uns_key, None)
    if prev_key is not None:
        adata.uns[f"{uns_key}_previous"] = prev_key


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
                                          f"analysis_{adata.uns['version']}_{timestamp()}.h5ad")
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


def update_h5_metadata(fin: typing.Union[anndata.AnnData, str],
                       metadata: dict):
    adata_in = fin
    if isinstance(fin, str):
        adata_in = load_anndata(fin)

    for key, value in dict.items():
        adata_in.uns[key] = value

    save_anndata(adata_in)


def save_rds_file(adata):
    """
    """
    counts = pd.DataFrame(adata.raw.X.todense(), columns=adata.raw.var_names, index=adata.obs_names).T
    counts = np.log1p(counts)

    features = pd.DataFrame({'Associated.Gene.Name': counts.index, 'Chromosome.Name': 1}, index=counts.index)

    genes = adata.obs['n_genes']
    genes[genes > genes.quantile(0.99)] = genes.quantile(0.99)
    counts.ix['ENSGGENES'] = genes
    umis = adata.obs['n_counts']
    umis[umis > umis.quantile(0.99)] = umis.quantile(0.99)
    counts.ix['ENSGUMI'] = umis

    features.ix['ENSGGENES'] = ['Genes', 1]
    features.ix['ENSGUMI'] = ['Umi', 1]

    tsne = pd.DataFrame(adata.obsm['X_umap'], columns=['V1', 'V2', 'V3'], index=adata.obs_names)
    tsne['dbCluster'] = adata.obs['louvain'].astype(int)

    spl.read_write._make_rds(counts, features, tsne,
                             os.path.join(adata.uns["output_dir"],
                                          f"{adata.uns['sampleid']}_{timestamp()}.Rds"))
