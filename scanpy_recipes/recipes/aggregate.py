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


def _fix_var_nans(obj, final_key, key_prefix, del_col=False):
    affected_cols = [key for key in obj.var.columns if key.startswith(key_prefix)]

    is_numeric = True
    if isinstance(obj.var[affected_cols[0]][0], str):
        is_numeric = False

    if is_numeric:
        obj.var[final_key] = obj.var[affected_cols].fillna(0).sum(axis=1)
    else:
        # this is primarily to get `gene_ids` to work
        # this should be guaranteed to work since the index is the union of all indexes,
        # so there should be at least one column where each row is defined.
        obj.var[final_key] = obj.var[affected_cols[0]].astype(str)
        for key in affected_cols[1:]:
            nan_locs = obj.var[final_key] == "nan"
            obj.var.loc[nan_locs, final_key] = obj.var.loc[nan_locs, key].astype(str)

    if del_col:
        obj.var.drop(affected_cols, inplace=True)


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

    # This is needed to fix NaNs cropping up due to `join="outer"` used in
    # `AnnData.concatenate`
    # as documented in the warning here:
    # https://anndata.readthedocs.io/en/latest/anndata.AnnData.concatenate.html#anndata.AnnData.concatenate
    _fix_var_nans(combined, "gene_ids", "gene_ids-", del_col=del_batch_var)
    _fix_var_nans(combined, "n_counts", "n_counts-", del_col=del_batch_var)
    _fix_var_nans(combined, "n_cells", "n_cells-", del_col=del_batch_var)

    # now we have to reconcile what's going on from the `uns` perspective
    combined.uns["sampleid"] = combined_sampleid
    combined.uns["sampleids"] = _combine_uns_data(adatas, "sampleid")
    combined.uns["sample_name"] = combined_sample_name
    combined.uns["sample_names"] = _combine_uns_data(adatas, "sample_name")

    combined.uns["output_dir"] = combined_output_dir
    combined.uns["is_aggregation"] = True
    # don't care about individual samples outputdirs so don't capture that

    if del_batch_uns:
        return combined

    copy_keys = [
        "analysis_pipeline_version", "analyst", "customer_name",
        "principal_investigator_name", "date_created", "obs_titles"
    ]
    for key in copy_keys:
        combined.uns[key] = adatas[0].uns[key]
    merge_keys = set(adatas[0].uns_keys()) - set(combined.uns_keys())
    for key in merge_keys:
        combined.uns[key] = _combine_uns_data(adatas, key)

    return combined


__api_objects__ = {
    "aggregate": aggregate,
}
