import random
import pandas as pd
from anndata import read_h5ad
from sclive.plotting.dimred_coexprs_2d_func import dimred_coexprs_2d


def test_dimred_plt_2d(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3").write_html("fig.html")
    

fig = test_dimred_plt_2d(None)
print(fig)