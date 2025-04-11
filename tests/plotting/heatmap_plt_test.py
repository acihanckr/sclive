from random import sample
from anndata import read_h5ad
from sclive.plotting.heatmap_plt_func import heatmap_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_heatmap_plt(scanpy_std_pipeline):
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], ticks_font_size=40).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], legend_font_size=40,cluster_columns=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10)).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cont_color="magma").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), use_raw=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), layer="counts").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), ticks_font_size=40).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_columns=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10),cluster_columns=True, cluster_rows=True, legend_font_size=12).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True, cluster_columns=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True, cluster_columns=True, width = 400).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True, cluster_columns=True, width = 400, height= 80).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True, cluster_columns=True, title_size=12).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True, cluster_columns=True, title_size=12, title="Heatmap").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cluster_rows=True, cluster_columns=True, scale_features=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], cluster_rows=True, cluster_columns=True, scale_features=True).write_html("fig.html")
