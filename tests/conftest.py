import os
import numpy as np
import scanpy as sc
import pytest

@pytest.fixture(scope="session")
def scanpy_std_pipeline():
    if not os.path.exists("test_data/pbmc.h5ad"):         
        adata = sc.read_10x_mtx("test_data/pbmc")
        adata.var_names_make_unique() 
        
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )

        fake_cats = ["A", "B", "C"]
        fake_dogs = ["ALp", "Bet"]
        adata.obs["fake_cats"] = np.random.choice(fake_cats, adata.n_obs)
        adata.obs["fake_dogs"] = np.random.choice(fake_dogs, adata.n_obs)

        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)

        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]

        sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
        adata.write("test_data/pbmc.h5ad")
