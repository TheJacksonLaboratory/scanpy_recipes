import scanpy.api as sc
from .read_write import __api_objects__ as read_write_api
from .recipes import preprocess, qc
from . import plotting as pl

for object_key, object in read_write_api.items():
    sc.__dict__[object_key] = object
#sc.AnalysisConfig = AnalysisConfig
#sc.save_rds_file = save_rds_file

sc.processing = preprocess
sc.qc = qc

del read_write_api
del preprocess, qc, pl
