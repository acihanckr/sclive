from sclive.dataio.get_dimred_with_exprs_func import get_dimred_with_exprs
from scanpy import read_h5ad

def test_get_dimred_with_exprs(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    adata_meta = get_dimred_with_exprs(adata, dimred_id = "umap", comps=[0,1], dimred_id_suffix="X_", genes=["SRM", "SIK1"])
    adata_meta = get_dimred_with_exprs(adata, dimred_id = "umap", comps=[0,1], dimred_id_suffix="X_", genes=["SRM", "FAKE1"])
    adata_meta = get_dimred_with_exprs(adata, dimred_id = "umap", comps=[0,1], dimred_id_suffix="X_", genes=["FAKE2", "FAKE1"])