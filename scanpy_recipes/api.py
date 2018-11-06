import scanpy.api as sc
from .read_write import __api_objects__ as read_write_api
from .recipes import preprocess, qc
from . import plotting as pl
from cmocean import cm as cmo

for object_key, object in read_write_api.items():
    sc.__dict__[object_key] = object
# sc.AnalysisConfig = AnalysisConfig
# sc.save_rds_file = save_rds_file

sc.processing = preprocess
sc.qc = qc

del read_write_api
del preprocess, qc, pl


sc.api.set_figure_params(
    scanpy=True,
    dpi=150,
    dpi_save=450,
    frameon=False,
    vector_friendly=True,
    color_map=cmo.matter,
    format="pdf",
    transparent=False,
    ipython_format="png2x",
)
