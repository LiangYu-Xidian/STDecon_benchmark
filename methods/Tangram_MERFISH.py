import pandas as pd
import numpy as np
from copy import deepcopy
import os
import anndata
import scanpy as sc
import logging
import time
start =time.time()

sc_adata = sc.read_h5ad("sc.h5ad")
st_adata = sc.read_h5ad("ST_100.h5ad")


celltype_key = 'cell_types'

sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

st_adata.layers["counts"] = st_adata.X.copy()
st_adata.obsm['spatial'] = st_adata.obsm['location']
sc.pp.normalize_total(st_adata, target_sum=10e4)
sc.pp.log1p(st_adata)
st_adata.raw = st_adata


'''
G = len(intersect)
print(G)
'''
import tangram as tg
tg.pp_adatas(sc_adata, st_adata, genes=None)
ad_map = tg.map_cells_to_space(
                   sc_adata,
                   st_adata,
                   mode='clusters',
                   cluster_label=celltype_key)
tg.project_cell_annotations(ad_map, st_adata, annotation=celltype_key)
celltype_density = st_adata.obsm['tangram_ct_pred']
celltype_density = (celltype_density.T/celltype_density.sum(axis=1)).T
celltype_density.to_csv('D:\MERFISH\Tangram\Tangram_MERFISH_100.tsv',sep='\t')
end = time.time()
print(end - start)