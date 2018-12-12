import scanpy.api as sc
from .read_write import __api_objects__ as read_write_api
from .recipes import preprocess, qc, tools
from . import plotting as pl
from .plotting.rcmod import update_figure_params

for object_key, object_ in read_write_api.items():
    sc.__dict__[object_key] = object_
# sc.AnalysisConfig = AnalysisConfig
# sc.save_rds_file = save_rds_file

sc.processing = preprocess
sc.qc = qc

for object_key, object_ in qc.__api_objects__.items():
    sc.qc.__dict__[object_key] = object_

for object_key, object_ in preprocess.__api_objects__.items():
    sc.pp.__dict__[object_key] = object_

for object_key, object_ in tools.__api_objects__.items():
    sc.tl.__dict__[object_key] = object_

update_figure_params()

sc.settings.verbosity = 3

del read_write_api
del preprocess, qc, tools, pl
