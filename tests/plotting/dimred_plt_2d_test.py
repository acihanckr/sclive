import random
import pandas as pd
import scanpy as sc
from anndata import read_h5ad
from sclive.plotting.dimred_plt_2d_func import dimred_plt_2d


def test_dimred_plt_2d(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")

    #standard use
    dimred_plt_2d(adata, "X_umap", "leiden").write_html("fig.html")
    dimred_plt_2d(adata, "X_pca", "leiden").write_html("fig.html")
    dimred_plt_2d(adata, "umap", "leiden", dimred_id_suffix="X_").write_html("fig.html")
    
    #numeric vs categorical metadata
    adata.obs["leiden_numeric"] = pd.to_numeric(adata.obs["leiden"])
    dimred_plt_2d(adata, "X_pca", "leiden_numeric").write_html("fig.html")
    dimred_plt_2d(adata, "X_pca", "leiden_numeric", cat=True).write_html("fig.html")
    
    #meta in obs and gene expressions
    adata.obs["CST3"] = adata.obs["leiden"]
    dimred_plt_2d(adata, "X_umap", "CST3", title_size=12).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "CST3", title_size=12, is_gene_exp=True).write_html("fig.html")
    
    #different components
    dimred_plt_2d(adata, "X_pca", "leiden", comps=[1,2]).write_html("fig.html")
    
    #selected barcodes
    dimred_plt_2d(adata, "X_pca", "leiden", selected_barcodes=random.choices(adata.obs_names, k=1000)).write_html("fig.html")
    
    #different layers and use_raw
    dimred_plt_2d(adata, "X_pca", "CST3", layer="counts").write_html("fig.html")
    dimred_plt_2d(adata, "X_pca", "CST3", layer="counts", use_raw=True).write_html("fig.html")
    
    #meta order and colors
    dimred_plt_2d(adata, "X_pca", "leiden", meta_order=random.choices(adata.obs["leiden"].cat.categories, k=10)).write_html("fig.html")
    dimred_plt_2d(adata, "X_pca", "leiden", meta_colors=["#FF0000"]*len(adata.obs["leiden"].cat.categories)).write_html("fig.html")
    
    #various aesthetics
    dimred_plt_2d(adata, "X_umap", "leiden", width=800).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", height=80).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=400).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height="true_asp_ratio").write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12, labels_size=12).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12, labels_size=12, title_size=12).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12, labels_size=12, title_size=12, title="Test").write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12, labels_size=12, title_size=12, title="Test", axis_font_size=12).write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12, labels_size=12, title_size=12, title="Test", axis_font_size=12, dimred_labels="PCA").write_html("fig.html")
    dimred_plt_2d(adata, "X_umap", "leiden", width=800,height=800, pt_size=5, legend_size=12, labels_size=12, title_size=12, title="Test", axis_font_size=12, dimred_labels=["PC1", "PC2"], ticks_font_size=12).write_html("fig.html")

test_dimred_plt_2d(None)
print("dimred_plt_2d test passed")