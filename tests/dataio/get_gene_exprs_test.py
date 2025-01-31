from sclive.dataio.get_gene_exprs_func import get_gene_exprs
from scanpy import read_h5ad

def test_get_gene_exprs(scanpy_std_pipeline):
    adata = read_h5ad("test_data/pbmc.h5ad")
    adata_meta = get_gene_exprs(adata, genes=["SRM", "SIK1"])
    adata_meta = get_gene_exprs(adata, genes=["SRM", "SIK1"], use_raw=True)
    adata_meta = get_gene_exprs(adata, genes=["SRM", "SIK1"], use_raw=False)
    adata_meta = get_gene_exprs(adata, genes=["SRM", "SIK1"], layer="counts", use_raw=True)
    adata_meta = get_gene_exprs(adata, genes=["SRM", "SIK1"], layer="counts")
