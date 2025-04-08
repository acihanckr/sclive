from random import sample
from anndata import read_h5ad
from sclive.plotting.dot_plt_func import dot_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_dot_plt(scanpy_std_pipeline):
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"]).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], scale_features=True).write_html("fig.html")
    dot_plt(adata, "fake_dogs", ["CST3", "SLC39A3"], scale_features=True).write_html("fig.html")
    dot_plt(adata, "fake_dogs", ["CST3", "SLC39A3"]).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], legend_size = 20).write_html("fig.html")
    dot_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), legend_size = 20).write_html("fig.html")
    dot_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10, )).write_html("fig.html")
    dot_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cont_color="magma").write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], legend_size = 40).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], cont_color="magma").write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], use_raw=True).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], layer="counts").write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], ticks_font_size=20).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], width = 400).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], width = 400, height= 80).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], title_size=12).write_html("fig.html")
    dot_plt(adata, "leiden", ["CST3", "SLC39A3"], title_size=12, title="Dotplot").write_html("fig.html")
