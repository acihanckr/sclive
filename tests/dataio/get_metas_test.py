from sclive.dataio import get_metas
from scanpy import read_h5ad
import pytest

def test_get_metas(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    get_metas(adata, meta_ids = ["leiden"])
    get_metas(adata, meta_ids = ["leiden"], cat=True)
    get_metas(adata, meta_ids = ["leiden", "fake_cats"])
    get_metas(adata, meta_ids = ["leiden", "fake_cats"], cat=True)
    get_metas(adata, meta_ids = ["total_counts"])
    get_metas(adata, meta_ids = ["total_counts"], cat=False)
    get_metas(adata, meta_ids = ["total_counts", "total_counts_mt"])
    get_metas(adata, meta_ids = ["total_counts", "total_counts_mt"], cat=False)
    with pytest.raises(ValueError):
        get_metas(adata, meta_ids = ["total_counts", "fake_cats"])
print("done")