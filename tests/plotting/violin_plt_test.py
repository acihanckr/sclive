from anndata import read_h5ad
from sclive.plotting.violin_plt_func import violin_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_violin_plt(scanpy_std_pipeline):
    violin_plt(adata, "leiden", "n_genes").write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", "fake_cats", vln_type="grouped").write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", "fake_dogs", vln_type="split", legend_font_size = 12).write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", "fake_dogs", pts="all", vln_type="split").write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", pts="all", vln_type="split", jitter=0.5).write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", "fake_cats", vln_type="grouped", pts="all").write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", axis_font_size=20).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", ticks_font_size=40).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", legend_font_size=20).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", legend_font_size=40).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", legend_font_size=20, legend_title="CST3 Expression").write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", legend_font_size=40, legend_title="CST3 Expression").write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", pts="all", pt_size=12).write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", legend_font_size=20).write_html("fig.html")
    violin_plt(adata, "leiden", "n_genes", "fake_cats", vln_type="grouped", pts="all", legend_font_size=20).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", width = 400).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", width = 400, height= 80).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", axis_font_size=12, title_size=12).write_html("fig.html")
    violin_plt(adata, "leiden", "CST3", axis_font_size=12, title_size=40, title="Violin Plot").write_html("fig.html")
