
#from ._version import get_versions
#__version__ = get_versions()['version']
#del get_versions
__version__ = "0.1"

import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
from .read_write import __api_objects__ as read_write_api
from .report import __api_objects__ as report_api
from .recipes import aggregate, preprocess, qc, tools
from .plotting import qc as plqc
from .plotting.rcmod import update_figure_params
from .utils import __api_objects__ as utils_api

for object_key, object_ in utils_api.items():
    sc.__dict__[object_key] = object_

for object_key, object_ in read_write_api.items():
    sc.__dict__[object_key] = object_

for object_key, object_ in report_api.items():
    sc.__dict__[object_key] = object_

for object_key, object_ in aggregate.__api_objects__.items():
    sc.__dict__[object_key] = object_

#sc.processing = preprocess
sc.qc = qc

for object_key, object_ in qc.__api_objects__.items():
    sc.qc.__dict__[object_key] = object_

for object_key, object_ in preprocess.__api_objects__.items():
    sc.pp.__dict__[object_key] = object_

for object_key, object_ in tools.__api_objects__.items():
    sc.tl.__dict__[object_key] = object_

for object_key, object_ in plqc.__api_objects__.items():
    sc.pl.__dict__[object_key] = object_

update_figure_params()

sc.settings.verbosity = 3

del read_write_api, report_api
del preprocess, qc, tools, plqc

sc.logging.print_version_and_date()
