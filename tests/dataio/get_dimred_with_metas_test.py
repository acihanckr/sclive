from sclive.dataio.get_dimred_with_metas_func import get_dimred_with_metas
from scanpy import read_h5ad


def test_get_dimred_with_metas(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    get_dimred_with_metas(adata, "X_umap", ["leiden"])
    get_dimred_with_metas(adata, dimred_id = "umap", comps=[0,1], dimred_id_suffix="X_", meta_ids=["total_counts", "n_genes_by_counts"])
