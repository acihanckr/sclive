from anndata import read_h5ad
from sclive.plotting.dimred_plt_2d_func import dimred_plt_2d


def test_dimred_plt_2d(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    fig = dimred_plt_2d(adata, "X_pca", "leiden")
