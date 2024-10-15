from CellDART.pred_cellf_celldart import pred_cellf_celldart
import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import time
start =time.time()

adata_sc = sc.read("sc.h5ad")
adata_sp = sc.read("ST_100.h5ad")

adata_sc.layers["counts"] = adata_sc.X.copy()
adata_sc.raw = adata_sc
adata_sp.layers["counts"] = adata_sp.X.copy()
adata_sp.obsm['spatial'] = adata_sp.obsm['location']
'''
G = 2000
sc.pp.filter_genes(adata_sc, min_counts=10)

adata_sc.layers["counts"] = adata_sc.X.copy()

sc.pp.highly_variable_genes(
    adata_sc,
    n_top_genes=G,
    subset=True,
    layer="counts",
    flavor="seurat_v3"
)


adata_sc.raw = adata_sc
adata_sp.layers["counts"] = adata_sp.X.copy()
adata_sp.obsm['spatial'] = adata_sp.obsm['location']
'''
adata_sp.raw = adata_sp
intersect = np.intersect1d(adata_sc.var_names, adata_sp.var_names)
adata_sp = adata_sp[:, intersect].copy()
adata_sc = adata_sc[:, intersect].copy()

adata_sp = pred_cellf_celldart(adata_sp=adata_sp, adata_sc=adata_sc, count_from_raw = False, gpu=True, celltype='cell_types', num_markers=15, nmix=8, npseudo=10000, alpha=0.6, alpha_lr=5, batch_size=64,
emb_dim=64, n_iterations=10000, init_train_epoch=30,  outdir='./CellDART_output_100', return_anndata=True)

end=time.time()
print(end-start)