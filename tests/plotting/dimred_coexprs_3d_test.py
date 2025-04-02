import random
from anndata import read_h5ad
from sclive.plotting.dimred_coexprs_3d_func import dimred_coexprs_3d


def test_dimred_coexprs_3d(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3").write_html("fig.html")
    
    #different components
    dimred_coexprs_3d(adata, "X_pca", "CST3", "SLC39A3", comps=[1,2, 3]).write_html("fig.html")
    
    #selected barcodes
    dimred_coexprs_3d(adata, "X_pca", "CST3", "SLC39A3", selected_barcodes=random.choices(adata.obs_names, k=1000)).write_html("fig.html")
    
    #different layers and use_raw
    dimred_coexprs_3d(adata, "X_pca", "CST3", "SLC39A3", layer="counts").write_html("fig.html")
    dimred_coexprs_3d(adata, "X_pca", "CST3", "SLC39A3", layer="counts", use_raw=True).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_pca", "CST3", "SLC39A3", use_raw=True).write_html("fig.html")
    
    #various aesthetics
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube").write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", plt_size=80).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="data",plt_size=400).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="auto",plt_size=800).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="manual",plt_size=800, pt_size=5).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5, title_size=12).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5, title_size=12, title="Test").write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5, title_size=12, title="Test", axis_font_size=12).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5, title_size=12, title="Test", axis_font_size=12, dimred_labels="PCA").write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5, title_size=12, title="Test", axis_font_size=12, dimred_labels=["PC1", "PC2", "PC3"], ticks_font_size=12).write_html("fig.html")
    dimred_coexprs_3d(adata, "X_umap", "CST3", "SLC39A3", aspectmode="cube", plt_size=800, pt_size=5, title_size=32, title="Test", axis_font_size=12, dimred_labels=["PC1", "PC2", "PC3"], ticks_font_size=12).write_html("fig.html")
