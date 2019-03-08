import os
from collections import namedtuple
from urllib.request import urlretrieve

from sklearn.datasets.base import _fetch_remote


RemoteFileMetadata = namedtuple(
    "RemoteFileMetadata",
    ["filename", "url", "checksum"]
)
module_dir = os.path.dirname(__file__)
TEST_DATASET_LOCATION = os.path.join(module_dir, "datasets/XX12345")


def download_10x_data():
    base_url = "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3"
    dataset_files = [
        RemoteFileMetadata(
            filename="filtered_feature_bc_matrix.h5",
            url=f"{base_url}/pbmc_1k_v3_filtered_feature_bc_matrix.h5",
            checksum="8191f576550c1b449d03441b9eb098ee9f73fa82513d171ba87d31d551e3ffda"
        ),
        RemoteFileMetadata(
            filename="raw_feature_bc_matrix.h5",
            url=f"{base_url}/pbmc_1k_v3_raw_feature_bc_matrix.h5",
            checksum="ed17bec5896a8208acfd6202e035fc36c65da14b3605aa577e93496b5d1ebd81"
        ),
        RemoteFileMetadata(
            filename="metrics_summary.csv",
            url=f"{base_url}/pbmc_1k_v3_metrics_summary.csv",
            checksum="992fa10cd3f7e1829c4a666f082f4a0e7e6d2c43c646f5c8cbf3991fda8fa172"
        )
    ]

    if not os.path.exists(TEST_DATASET_LOCATION):
        os.makedirs(TEST_DATASET_LOCATION)

    if not all([os.path.exists(os.path.join(TEST_DATASET_LOCATION, metadata.filename))
        for metadata in dataset_files]):
        print(f"downloading 10x dataset from {base_url} to {TEST_DATASET_LOCATION}.")
        for metadata in dataset_files:
            _fetch_remote(metadata, dirname=TEST_DATASET_LOCATION)
    else:
        print(f"10x dataset already in {TEST_DATASET_LOCATION}.")
