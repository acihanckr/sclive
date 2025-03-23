from anndata import read_h5ad
from sclive.plotting.box_plt_func import box_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_box_plt(scanpy_std_pipeline):
    box_plt(adata, "leiden", "n_genes").write_html("fig.html")
    box_plt(adata, "leiden", "n_genes", "fake_cats", box_type="grouped").write_html("fig.html")
    box_plt(adata, "leiden", "n_genes", pts="all").write_html("fig.html")
    box_plt(adata, "leiden", "n_genes", "fake_cats", box_type="grouped", pts="all").write_html("fig.html")
    box_plt(adata, "leiden", "CST3", axis_font_size=20).write_html("fig.html")
    box_plt(adata, "leiden", "CST3", axis_font_size=None).write_html("fig.html")
    box_plt(adata, "leiden", "CST3", ticks_font_size=20).write_html("fig.html")
    box_plt(adata, "leiden", "n_genes", pts="all", pt_size=12).write_html("fig.html")
    box_plt(adata, "leiden", "n_genes", legend_size=20).write_html("fig.html")
    box_plt(adata, "leiden", "n_genes", "fake_cats", box_type="grouped", pts="all", legend_size=20).write_html("fig.html")
    box_plt(adata, "leiden", "CST3", width = 400).write_html("fig.html")
    box_plt(adata, "leiden", "CST3", width = 400, height= 80).write_html("fig.html")
    box_plt(adata, "leiden", "CST3", axis_font_size=12, title_size=12).write_html("fig.html")
    box_plt(adata, "leiden", "CST3", axis_font_size=12, title_size=12, title="Box Plot").write_html("fig.html")

test_box_plt(None)