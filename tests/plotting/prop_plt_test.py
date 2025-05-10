from anndata import read_h5ad
from sclive.plotting.prop_plt_func import prop_plt

adata = read_h5ad("test_data/pbmc.h5ad")

def test_prop_plt(scanpy_std_pipeline):
    prop_plt(adata, "leiden", "leiden", plt_type="pct").write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="pct").write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="pct", x_order=["0", "1", "2", "3", "4", "5", "6"]).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="pct", x_order=["0", "1", "2", "3", "4", "5", "6"], group_order=["C","A", "B"]).write_html("fig.html")
    prop_plt(adata, "fake_dogs", "fake_dogs", plt_type="count").write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count").write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", stacked=False).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", coord_flip=True, axis_font_size=None).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", width = 400).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", width = 400, height= 80).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", axis_font_size=40, axis_labels=["Group1", "Group2"]).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", ticks_font_size=40).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", ticks_font_size=None).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", axis_font_size=12, title_size=12).write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", axis_font_size=12, title_size=12, title="Fake Cats").write_html("fig.html")
    prop_plt(adata, "leiden", "fake_cats", plt_type="count", axis_font_size=12, title_size=40, title="Fake Cats").write_html("fig.html")
