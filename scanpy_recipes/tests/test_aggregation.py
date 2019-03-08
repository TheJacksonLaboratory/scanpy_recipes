#
## coding: utf-8
#
## # Example aggregation analysis notebook
## ### Author: Bill Flynn (bill.flynn@jax.org)
##
## Updates:
## - 2018-12-31: creation
#
## This notebook assumes you have some familiarity with the `scanpy_recipes` extensions.  If not, check out `examples/Analysis-Template.ipynb` for more details.  This will be an abbreviated notebook just showing how things work with an aggregated sample.
#
## In[1]:
#
#
#from scanpy_recipes.api import sc
#
#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
#
#
## ## Loading data
#
## In[2]:
#
#
#ac = sc.AnalysisConfig()
#
#
## In[3]:
#
#
#config_string = """
#[names]
#customer_name = Anonymous person
#analyst_name = Bill Flynn
#analysis_name = Test-analysis
#
#[sample_names]
#PR18016 = CD45+
#PR18017 = CD45-
#
#[genomes]
#PR18016 = GRCh38
#PR18017 = GRCh38
#
#[species]
#hg19 = hsapiens
#GRCh38 = hsapiens
#mm10 = mmusculus
#
#[input_dirs]
#PR18016 = /projects/flynnb/singlecell/IBC/PR18016/
#PR18017 = /projects/flynnb/singlecell/IBC/PR18017/
#
#[output_dirs]
#PR18016 = /fastscratch/flynnb/test-outputs
#PR18017 = /fastscratch/flynnb/test-outputs
#"""
#config = ac.read(config_string)
#
#
## Below, I'll be storing the adata objects in an iterable so that you can get the sense for how this might scale with many samples.
#
## In[4]:
#
#
#combined_sampleid = "PR18016-PR18017"
#combined_sample_name = "CD45+/-"
#combined_output_dir = "/fastscratch/flynnb/test-aggr-outputs"
#
#
## In[5]:
#
#
#sampleids = config["sample_names"]
#
#
## In[6]:
#
#
#adata_raws = dict((sampleid, sc.load_10x_data(sampleid, config)) for sampleid in sampleids)
#
#
## ## Quality control and filtering
##
## Analysis first starts by generating **per-cell** and **per-gene** metrics, that we can then use to filter the data.  We compute these before we do any aggregation
#
## In[7]:
#
#
#for adata_raw in adata_raws.values():
#    sc.qc.gen_qc(adata_raw)
#
#
## In[8]:
#
#
#for adata_raw in adata_raws.values():
#    sc.pl.qc_violins(adata_raw)
#
#
## Either use one set of QC params that can be used all samples or use individual params. Below is one way to do this.
#
## In[9]:
#
#
#qc_params = dict(
#    min_cells_per_gene=3,
#    min_counts_per_gene=3,
#    min_counts_per_cell=500,
#    min_genes_per_cell=200,
#    sequencing_saturation=None,
#    percent_mito=10.0,
#    rbc_threshold=10
#)
#qc_params = dict((sampleid, qc_params.copy()) for sampleid in adata_raws.keys())
#
#
## In[10]:
#
#
#qc_params["PR18017"].update(percent_mito=5.0)
#
#
## In[11]:
#
#
#trials = [sc.qc.run_qc(adata_raw, trial=True, **qc_params[sampleid])
#          for sampleid, adata_raw in adata_raws.items()]
#
#
## In[12]:
#
#
#qc_figs = [sc.pl.qc_violins(trial, return_fig=True) for trial in trials]
#
#
## In[13]:
#
#
#adata_qcs = [sc.qc.run_qc(adata_raw, trial=False, **qc_params[sampleid])
#             for sampleid, adata_raw in adata_raws.items()]
#
#
## In[14]:
#
#
#qc_figs2 = [sc.pl.genes_umis_scatter(trial, return_fig=True) for trial in trials]
#
#
## ## Aggregation time
## Now let's aggregate them together.
#
## In[15]:
#
#
#combined_qc = sc.aggregate(
#    adata_qcs,
#    combined_output_dir=combined_output_dir,
#    combined_sample_name=combined_sample_name,
#    combined_sampleid=combined_sampleid,
#    make_output_dir=True,
#    del_batch_var=False
#)
#
#
## You can see that the following keys are added to `.obs`:
## * `batch`
## * `sampleid`
## * `sample_name`
##
## The keys added to `.var` are pretty redundant:
## * `gene_ids == gene_ids-0 == gene_ids-1 == ...`
## * `n_counts = n_counts-0 + n_counts-1 + ...`
## * `n_genes = n_genes-0 + n_genes-1 + ...`
##
## I can envision how having the counts/genes on a per-sample basis could be useful, so these are left in by default, but all these extra `.var` keys can be removed by passing `del_batch_var=True`.
#
## In[16]:
#
#
#print("New `.obs` keys")
#print(f"old: {adata_raws['PR18016'].obs_keys()}\nnew: {combined_qc.obs_keys()}")
#print("\nNew `.var` keys")
#print(f"old: {adata_raws['PR18016'].var_keys()}\nnew: {combined_qc.var_keys()}")
#
#
## The `.uns` data also gets an overhaul.  Entries that describe this analysis remain the same, but others that describe per-sample data have now become nested dictionaries, accessible via sampleID, e.g.
## ```
## In : combined_qc.uns["10x_metrics"].keys()
## Out: dict_keys(['PR18016', 'PR18017'])
## In : combined_qc.uns["qc_cell_filter"].keys()
## Out: dict_keys(['PR18016', 'PR18017'])
## ```
## Important entries like `sample_name` and `sampleid` have been replaced by the `combined_` versions defined during `sc.aggregate` and the previous per-sample versions are moved to `sample_names` and `sampleids`, respectively.
#
## In[17]:
#
#
#combined_qc.uns_keys()
#
#
## Everything from here continues **mostly** as usual.  Importantly, there is one additional key added to `.uns`:
## ```
## In : combined_qc.uns["is_aggregation"]
## Out: True
## ```
## This signals to a few other methods to handle the sample differently.  For example, dimensionality reduction uses `bbknn` to construct the neighborhood graph and RDS saving will add an extra pseudogene denoting `batch`.
#
## In[18]:
#
#
#qc_save_file = sc.save_adata(combined_qc, "qc")
#
#
## ## HVG selection, dimensionality reduction, and clustering
#
## In[19]:
#
#
#adata_qc = sc.read_h5ad(qc_save_file)
#
#
## In[20]:
#
#
#adata_full, adata_redux = sc.pp.preprocess(adata_qc, n_top_genes=1000, scale=True)
#
#
## In[21]:
#
#
#sc.pp.dimensionality_reduction(adata_redux,
#                               n_neighbors=10, min_dist=0.5)
#
#
## We can visualize the aggregation by `"sample_name"`, `"sampleid"`, or `"batch"`!
#
## In[54]:
#
#
#sc.pl.scatter(adata_redux, basis="umap", color="sample_name", show=False)
#
#
## In[23]:
#
#
#sc.tl.cluster(adata_redux,)
#
#
## In[24]:
#
#
#sc.pl.scatter(adata_redux, basis="umap", color="cluster", legend_loc="on data")
#
#
## ## Saving outputs
#
## In[27]:
#
#
#sc.save_adata_to_rds(adata_redux, cluster_key="cluster")
#
#
## In[28]:
#
#
#full_save_file = sc.save_adata(adata_full, "full")
#redux_save_file = sc.save_adata(adata_redux, "redux")
#
#
## ## Report generation
#
## Report generation currently under construction for aggregated datasets.
#
## In[25]:
#
#
#report = sc.SCBLReport()
#
#
## In[27]:
#
#
#report.add_report_figures(
#    adata_redux,
#    violins=qc_figs,
#    scatters=qc_figs2,
#    cluster_key="cluster",
#    batch_key="batch"
#)
#
#
## In[59]:
#
#
#report.generate_report(adata_redux)
#
