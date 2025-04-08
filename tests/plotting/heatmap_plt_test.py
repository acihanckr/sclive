from anndata import read_h5ad
from sclive.plotting.heatmap_plt_func import heatmap_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_heatmap_plt(scanpy_std_pipeline):
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cont_color="reds").write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], use_raw=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], layer="counts").write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], ticks_font_size=20).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_columns=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True, width = 400).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True, width = 400, height= 80).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True, title_size=12).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True, title_size=12, title="Heatmap").write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True, scale_features=True).write_html("fig.html")
