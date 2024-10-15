import tempfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvi
import seaborn as sns
import torch
from cell2location.models import Cell2location, RegressionModel
from cell2location.plt import plot_spatial
from cell2location.utils import select_slide
from cell2location.utils.filtering import filter_genes
import time
start =time.time()
results_folder = "D:\MERFISH\Cell2location/"
sc_adata = sc.read_h5ad( "sc.h5ad")
st_adata = sc.read_h5ad( "ST_100.h5ad")
'''
from cell2location.utils.filtering import filter_genes
selected = filter_genes(sc_adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
sc_adata = sc_adata[:, selected].copy()
'''
RegressionModel.setup_anndata(adata=sc_adata,
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_types',
                       )
mod = RegressionModel(sc_adata)
mod.train(
    max_epochs=250,
    batch_size=2500,
    train_size=1,
    lr=0.002,
)

# mod.plot_history(20)
sc_adata = mod.export_posterior(
    sc_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
# Results saving folder
ref_run_name = f'./reference_signatures_100'
run_name = f"{results_folder}/cell2location_map_100"
# Save model
mod.save(f"{ref_run_name}", overwrite=True)

adata_file = f"{ref_run_name}/sc.h5ad"
sc_adata.write(adata_file)

if "means_per_cluster_mu_fg" in sc_adata.varm.keys():
    inf_aver = sc_adata.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in sc_adata.uns["mod"]["factor_names"]]
    ].copy()
else:
    inf_aver = sc_adata.var[
        [f"means_per_cluster_mu_fg_{i}" for i in sc_adata.uns["mod"]["factor_names"]]
    ].copy()

inf_aver.columns = sc_adata.uns["mod"]["factor_names"]

intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
inf_aver= inf_aver.loc[intersect, :].copy()



Cell2location.setup_anndata(adata=st_adata)
mod = Cell2location(
    st_adata,
    cell_state_df= inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=6,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200,
)
mod.train(
    max_epochs=10000,
    # train using full data (batch_size=None)
    batch_size=None,
    # use all data points in training because
    # we need to estimate cell abundance at all locations
    train_size=1,
)

st_adata = mod.export_posterior(
    st_adata, sample_kwargs={'num_samples': 500, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

mod.save(f"{run_name}", overwrite=True)
end = time.time()
adata_file = f"{run_name}/sp.h5ad"
st_adata.write(adata_file)

st_adata.obs[st_adata.uns['mod']['factor_names']] = st_adata.obsm['q05_cell_abundance_w_sf']
result1 = st_adata.obsm['q05_cell_abundance_w_sf']
result2 = st_adata.obsm['q95_cell_abundance_w_sf']
result3 = st_adata.obsm['means_cell_abundance_w_sf']
sum_result_3 = result3.sum(axis=1)
result3_percent = result3.div(result3.assign(total=sum_result_3)['total'], axis='index')
result3_percent.to_csv("D:\MERFISH\Cell2location\Cell2location_MERFISH_100.csv")
print(end - start)
