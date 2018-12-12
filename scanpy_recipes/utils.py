import sys
import contextlib
from datetime import datetime
from anndata import logging


def timestamp():
    """
    Returns current time in format: YYYY-MM-DDThh-mm-ss
    """
    return datetime.now().strftime("%Y-%m-%dT%H-%M-%S")


def datestamp():
    return datetime.now().strftime("%Y%m%d")


def quantile_limit(obj, key, q=0.99):
    obj = obj[key].copy()
    obj[obj > obj.quantile(0.99)] = obj.quantile(0.99)
    return obj


def shift_metadata(adata, uns_key):
    prev_key = adata.get(uns_key, None)
    if prev_key is not None:
        adata.uns[f"{uns_key}_previous"] = prev_key


class _DummyFile(object):
    def write(self, x): pass
@contextlib.contextmanager
def silence():
    save_stdout = sys.stdout
    #logging.logging.setLevel('ERROR')
    sys.stdout = _DummyFile()
    yield
    sys.stdout = save_stdout
    #logging.logging.setLevel('INFO')
