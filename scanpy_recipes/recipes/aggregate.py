import os
import collections

def _combine_uns_data(adatas, key):
    comb_uns = dict()
    for adata in adatas:
        sampleid = adata.uns["sampleid"]
        obj = adata.uns.get(key, None)
        if hasattr(obj, "__copy__"):
            obj = obj.copy()
        comb_uns[sampleid] = obj
    return comb_uns

def aggregate(*adatas, combined_output_dir=None,
              combined_sampleid=None, combined_sample_name=None,
              batch_key="batch", del_batch_var=False,
              make_output_dir=True, del_batch_uns=False):
    """

    """

    if len(adatas) <= 1:
        if isinstance(adatas, collections.Sequence):
            adatas = adatas[0]
        if len(adatas) <= 1:
            return adatas
    combined = adatas[0].concatenate(*adatas[1:], join="outer", batch_key=batch_key)

    sampleids = [adata.uns["sampleid"] for adata in adatas]
    combined.obs["sampleid"] = combined.obs.batch.map(
        dict(zip(combined.obs.batch.cat.categories, sampleids))
    )

    sample_names = [adata.uns.get("sample_name", None) for adata in adatas]
    if any(filter(None, sample_names)):
        combined.obs["sample_name"] = combined.obs.batch.map(
            dict(zip(combined.obs.batch.cat.categories, sample_names))
        )

    if not combined_sampleid:
        combined_sampleid = "-".join(sampleids)
    if not combined_sample_name:
        combined_sample_name = f"{combined_sampleid} aggregation"
    if not combined_output_dir:
        combined_output_dir = os.path.join(
            os.path.dirname(adatas[0].uns["output_dir"]),
            f"{combined_sampleid}_outputs"
        )
    if make_output_dir and not os.path.exists(combined_output_dir):
        os.mkdir(combined_output_dir)

    combined.var["gene_ids"] = combined.var["gene_ids-0"]
    n_counts_cols = [k for k in combined.var_keys() if k.startswith("n_counts")]
    combined.var["n_counts"] = combined.var[n_counts_cols].sum(axis=1)
    n_cells_cols = [k for k in combined.var_keys() if k.startswith("n_cells")]
    combined.var["n_cells"] = combined.var[n_cells_cols].sum(axis=1)

    # I guess you might be interested in knowing, for example, how many counts you have
    # for each gene coming from each sample
    if del_batch_var:
        combined.var.drop(n_counts_cols + n_cells_cols, inplace=True)

    # now we have to reconcile what's going on from the `uns` perspective
    combined.uns["sampleid"] = combined_sampleid
    combined.uns["sampleids"] = _combine_uns_data(adatas, "sampleid")
    combined.uns["sample_name"] = combined_sample_name
    combined.uns["sample_names"] = _combine_uns_data(adatas, "sample_name")

    combined.uns["output_dir"] = combined_output_dir
    # don't care about individual samples outputdirs so don't capture that

    if del_batch_uns:
        return combined

    copy_keys = ["analysis_version", "analyst", "customer_name", "date_created", "obs_titles"]
    for key in copy_keys:
        combined.uns[key] = adatas[0].uns[key]
    merge_keys = set(adatas[0].uns_keys()) - set(combined.uns_keys())
    for key in merge_keys:
        combined.uns[key] = _combine_uns_data(adatas, key)

    return combined


__api_objects__ = {
    "aggregate": aggregate,
}
