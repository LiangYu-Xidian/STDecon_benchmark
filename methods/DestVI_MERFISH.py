import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scvi.model import CondSCVI,DestVI
import time
start =time.time()


sc_adata = sc.read("sc.h5ad")
st_adata = sc.read("ST_100.h5ad")
'''
G = 2000
sc.pp.filter_genes(sc_adata, min_counts=1)


sc_adata.layers["counts"] = sc_adata.X.copy()
sc.pp.highly_variable_genes(
    sc_adata,
    n_top_genes=G,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)
'''
sc_adata.layers["counts"] = sc_adata.X.copy()
sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata
st_adata.layers["counts"] = st_adata.X.copy()
st_adata.obsm['spatial'] = st_adata.obsm['location']

sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)

CondSCVI.setup_anndata(sc_adata, layer="counts", labels_key="cell_types")

sc_model = CondSCVI(sc_adata, weight_obs=False)
sc_model.view_anndata_setup()

sc_model.train(max_epochs=1000,lr=0.001)

DestVI.setup_anndata(st_adata, layer="counts")
st_model = DestVI.from_rna_model(st_adata, sc_model)
st_model.view_anndata_setup()
st_model.train(max_epochs=2000)
st_model.save("D:/MERFISH/DestVI/st_100",overwrite=True)
st_adata.obsm["proportions"] = st_model.get_proportions()
output_matrix = st_adata.obsm["proportions"]
st_adata.write_h5ad("st_final_100.h5ad")

output_matrix.to_csv("DestVI_MERFISH_100.csv",index=True)
end =time.time()
print(end-start)
