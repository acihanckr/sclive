# ScLive
ScLive provides interactive and highly customizable plots commonly used in single cell analysis.
It also has live_dash module to create single cell analysis dashboards from anndata objects.
The resulting dashboard allows to customize various aspects of the plot using an user interface.

ScLive can be installed using pip

```bash 
pip install sclive
```

Once installed and annotated data is loaded, the following code is enough to create a single cell analysis dashboard:
```python
from sclive.live_dash import create_dash_app, ScLiveDash

def test_box_plt(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    sclive_dash = ScLiveDash(adata)
    create_dash_app(sclive_dash)

```