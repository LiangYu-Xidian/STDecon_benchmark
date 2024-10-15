import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata

import scvi
from scvi.external import RNAStereoscope, SpatialStereoscope
import time
start = time.time()
scvi.settings.reset_logging_handler()

sc_adata = sc.read_h5ad("sc.h5ad")
st_adata = sc.read_h5ad("ST_100.h5ad")
#print(st_adata)

sc_adata.layers["counts"] = sc_adata.X.copy()


st_adata.layers["counts"] = st_adata.X.copy()
st_adata.obsm['spatial'] = st_adata.obsm['location']


intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()

RNAStereoscope.setup_anndata(sc_adata, labels_key="cell_types")

sc_stereo = RNAStereoscope(sc_adata, )
sc_stereo.train(lr=0.01, max_epochs=30000,batch_size=1000)
SpatialStereoscope.setup_anndata(st_adata)
st_stereo = SpatialStereoscope.from_rna_model(st_adata, sc_stereo)

st_stereo.train(lr=0.01,max_epochs=30000,batch_size=1000)
proportions =st_stereo.get_proportions()
end = time.time()
proportions.to_csv("D:\MERFISH\Stereoscope\Stereoscope_MERFISH_100.csv",index=True)
print(end - start)