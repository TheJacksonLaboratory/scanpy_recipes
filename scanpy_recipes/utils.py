import sys
import contextlib
import functools
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


def play_nice_with_categories(func):
    """
    Any function that modifies a dataframe categorical dataframe column, e.g.
    `adata.obs['cluster']`, will first have the column transformed to 'object' dtype (or
    optionally 'int' dtype) for it to operate on and return.  The returning column will
    then be transformed into an ordered categorical column again.
    """
    @functools.wraps(func)
    def wrapper(df, *args, has_int=False, **kwargs):
        """
        expects just a single adata.obs[key] column
        """
        df = df.astype("str")
        if has_int:
            df = df.astype("int")
        new_df = func(df, *args, **kwargs)
        if has_int:
            new_df.cat.reorder_categories(
                new_df.cat.categories.astype(int).sort_values().astype(str),
                ordered=True, inplace=True
            )
        new_df = df.astype("str").astype("category")
        return new_df
    return wrapper


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
    mapping = dict(zip(uniq, new_ints)
    return obs.map(mapping)


@play_nice_with_categories(has_int=True)
def reorder_clusters(obs, subset):
    """
    Specify new ordering.  Can only be a subset
    """
    subset_index = pd.Index(subset)
    subset_in_df = subset_index.isin(obs)
    if not subset_in_df.all():
        logg.error(f"Cluster(s) [{subset_index[~subset_in_df]}] not found in obs!")
        raise ValueError

    new_ints = np.arange(min(subset), max(subset)+1, dtype=int)
    mapping = dict(zip(subset, new_ints))
    return obs.map(mapping)


__api_objects__ = {
    "reorder_clusters": reorder_clusters,
    "datestamp": datestamp,
    "timestamp": timestamp,
}
