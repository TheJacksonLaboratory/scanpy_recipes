import scanpy.api as sc
from .read_write import AnalysisConfig, save_rds_file
from .recipes import preprocess, qc
from . import plotting as pl

sc.AnalysisConfig = AnalysisConfig
sc.save_rds_file = save_rds_file

sc.processing = preprocess
sc.qc = qc

del AnalysisConfig, save_rds_file
del preprocess, qc, pl
