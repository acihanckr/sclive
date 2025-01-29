from sclive.adata import get_metas
from scanpy import read_h5ad

def test_get_metas():
    adata = read_h5ad("test_data/pbmc.h5ad")
    adata_meta = get_metas(adata, meta_ids = ["leiden"])
    adata_meta = get_metas(adata, meta_ids = ["leiden"], cat=True)
    adata_meta = get_metas(adata, meta_ids = ["leiden", "fake_cats"])
    adata_meta = get_metas(adata, meta_ids = ["leiden", "fake_cats"], cat=True)
    adata_meta = get_metas(adata, meta_ids = ["total_counts"])
    adata_meta = get_metas(adata, meta_ids = ["total_counts"], cat=False)
    adata_meta = get_metas(adata, meta_ids = ["total_counts", "total_counts_mt"])
    adata_meta = get_metas(adata, meta_ids = ["total_counts", "total_counts_mt"], cat=False)
    adata_meta = get_metas(adata, meta_ids = ["total_counts", "fake_cats"])
print("done")