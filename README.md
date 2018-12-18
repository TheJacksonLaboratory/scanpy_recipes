# scanpy_recipes
Generalized and aggregated analysis functions for working with [`Anndata`][anndata] and
[`Scanpy`][scanpy].

## Introduction

This package is an extension to the [`scanpy`][scanpy] API tentatively called
`scanpy_recipes`, inspired by the way `scanpy` references their implementations of
previous workflows as recipes.

### Aims

There are three main aims of the additional functions of `scanpy_recipes`:

1. simplify scRNA-seq analysis in JAX's Single Cell Biology Lab.
2. use the `AnnData` object to capture and retain as much metadata as possible.  In
   theory, I'd like the `AnnData` object to be a self-contained description of the
   experiment, analysis, and results.
3. integrate seamlessly with the normal scanpy API and as such only extend functionality.

Specifically, `scanpy` is normally loaded using `import scanpy.api as sc`. Instead, this
package provides this api object as
```{python}
from scanpy_recipes.api import sc
```
This is the same `scanpy.api` object, does not change any `scanpy` functionality, and only
adds additional functionality.

### Workflow

The workflow has 3 components with all three generating their own `AnnData` object.  By
default, `scanpy` will modify `AnnData` objects inplace.  Because we're working in a
collaborative environment where multiple users of varying skill levels will be using this
framework, I wanted this inplace modification to be crystal clear, especially since
finding highly variable genes will remove all genes except the highly variable ones from
the default expression matrix (they can still be found in the `AnnData.raw`).

These three components are:

-   Specify sample metadata, input/output directories up from and loading `adata_raw`
-   Iteratively filter and QC data, ultimately yielding `adata_qc`
-   Generate `adata_redux` through HVG selection, dimensionality reduction, and clustering.

Each of these objects will be saved, and you can restart any step by loading the object
output by the previous step.

See [examples/Analysis-Template.ipynb] for an example analysis.

## Getting Started

### Prerequisites

This library doesn't have too many requirements besides `scanpy` and the basic Python
scientific stack (`numpy`, `pandas`, `matplotlib`).

### Installing

Until we hit a major release, installation should be done through pip egg installs:
```{bash}
pip install -e git+https://github.com/TheJacksonLaboratory/scanpy_recipes.git#egg=scanpy_recipes
```
or
```{bash}
git clone https://github.com/TheJacksonLaboratory/scanpy_recipes.git
cd scanpy_recipes
pip install -e .
```

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the
process for submitting pull requests to us.

## Authors

-   Bill Flynn [@wflynny][wflynny]

## License

Since this package extends [`scanpy`][scanpy], this work is under the same BSD-3 license.
See [LICENSE](LICENSE) for more details.


[anndata]: github.com/theislab/anndata
[scanpy]: github.com/theislab/scanpy
[wflynny]: github.com/wflynny
