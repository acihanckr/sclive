import json
import polars as pl
import polars.selectors as cs
from anndata import read_h5ad
from sclive.live_dash import create_dash_app
from sclive.live_dash import ScLiveDash

adata = read_h5ad("test_data/pbmc.h5ad")
sclive_dash = ScLiveDash(adata)
app = create_dash_app(sclive_dash)

