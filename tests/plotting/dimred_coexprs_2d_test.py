import random
from anndata import read_h5ad
from sclive.plotting.dimred_coexprs_2d_func import dimred_coexprs_2d

def test_dimred_plt_2d(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3").write_html("fig.html")
    
    #different components
    dimred_coexprs_2d(adata, "X_pca", "CST3", "SLC39A3", comps=[1,2]).write_html("fig.html")
    
    #selected barcodes
    dimred_coexprs_2d(adata, "X_pca", "CST3", "SLC39A3", selected_barcodes=random.choices(adata.obs_names, k=1000)).write_html("fig.html")
    
    #different layers and use_raw
    dimred_coexprs_2d(adata, "X_pca", "CST3", "SLC39A3", layer="counts").write_html("fig.html")
    dimred_coexprs_2d(adata, "X_pca", "CST3", "SLC39A3", layer="counts", use_raw=True).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_pca", "CST3", "SLC39A3", use_raw=True).write_html("fig.html")
    
    #various aesthetics
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", height=80).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=400).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height="true_asp_ratio").write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=800, pt_size=5).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=800, pt_size=5, title_size=12).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=800, pt_size=5, title_size=12, title="Test").write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=800, pt_size=5, title_size=12, title="Test", axis_font_size=12).write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=800, pt_size=5, title_size=12, title="Test", axis_font_size=12, dimred_labels="PCA").write_html("fig.html")
    dimred_coexprs_2d(adata, "X_umap", "CST3", "SLC39A3", width=800,height=800, pt_size=5, title_size=12, title="Test", axis_font_size=12, dimred_labels=["PC1", "PC2"], ticks_font_size=12).write_html("fig.html")
