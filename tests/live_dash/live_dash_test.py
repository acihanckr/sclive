from anndata import read_h5ad
from sclive.live_dash import create_dash_app
from sclive.live_dash import ScLiveDash


def test_box_plt(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    sclive_dash = ScLiveDash(adata)
    create_dash_app(sclive_dash)
