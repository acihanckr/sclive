from random import sample
from anndata import read_h5ad
from sclive.plotting.heatmap_plt_func import heatmap_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_heatmap_plt(scanpy_std_pipeline):
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], meta_order=["0", "1", "2", "3", "4", "5", "6"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], meta_order=["0", "1", "2", "3", "4", "5", "6"], gene_order=["SLC39A3", "CST3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], meta_order=["0", "1", "2", "3", "4", "5", "6"], gene_order=["CST3", "SLC39A3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], gene_order=["CST3", "SLC39A3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], gene_order=["SLC39A3", "CST3"]).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], scale_features=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], ticks_font_size=40).write_html("fig.html")
    heatmap_plt(adata, "leiden", ["CST3", "SLC39A3"], legend_font_size=40,legend_title="Gene Expression").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10)).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), scale_features=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10),legend_font_size=12).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), ticks_font_size=40).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), cont_color="magma").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), use_raw=True).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), layer="counts").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), width = 400).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), width = 400, height= 80).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), title_size=12).write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), title_size=40, title="Heatmap").write_html("fig.html")
    heatmap_plt(adata, "leiden", sample(adata.var_names.to_list(), k=10), scale_features=True).write_html("fig.html")
