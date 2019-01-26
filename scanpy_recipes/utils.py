import sys
import contextlib
import functools
import numpy as np
from datetime import datetime
from anndata import logging

import scanpy.api.logging as logg

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


def play_nice_with_categories(has_int=False):
    """
    Any function that modifies a dataframe categorical dataframe column, e.g.
    `adata.obs['cluster']`, will first have the column transformed to 'object' dtype (or
    optionally 'int' dtype) for it to operate on and return.  The returning column will
    then be transformed into an ordered categorical column again.
    """
    def decorator_play_nice_with_categories(func):
        @functools.wraps(func)
        def wrapper(df, *args, **kwargs):
            """
            expects just a single adata.obs[key] column
            """
            df = df.astype("str")
            if has_int:
                df = df.astype("int")
            new_df = func(df, *args, **kwargs)
            new_df = new_df.astype("str").astype("category")
            if has_int:
                new_df.cat.reorder_categories(
                    new_df.cat.categories.astype(int).sort_values().astype(str),
                    inplace=True
                )
            return new_df
        return wrapper
    return decorator_play_nice_with_categories


@play_nice_with_categories(has_int=True)
def shift_clusters(obs):
    if obs.min() == 0:
        obs += 1
    return obs


@play_nice_with_categories(has_int=True)
def order_clusters(obs):
    return obs


@play_nice_with_categories(has_int=True)
def reset_int_category(obs):
    uniq = sorted(obs.unique())
    new_ints = np.arange(1, len(uniq) + 1, dtype=int)
    mapping = dict(zip(uniq, new_ints))
    return obs.map(mapping)


__api_objects__ = {
    "datestamp": datestamp,
    "timestamp": timestamp,
}
