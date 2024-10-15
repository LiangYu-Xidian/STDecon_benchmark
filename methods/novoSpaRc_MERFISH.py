import pandas as pd
import numpy as np
from copy import deepcopy
import os

import novosparc as nc
from scipy.spatial.distance import cdist
import time
start = time.time()

st_coords = pd.read_csv('D:\MERFISH\locations_100.csv')
st_coords.columns = ['barcode','X','Y']
st_counts = pd.read_csv('D:\MERFISH\spatialcount_100.csv')
st_counts = st_counts.rename({'Unnamed: 0':'barcode'}, axis=1)

sc_counts = pd.read_csv('D:\MERFISH/raw_somatosensory_sc_exp.txt',sep='\t')
sc_counts = sc_counts.set_index('cell_id')
sc_counts = sc_counts.T
sc_labels = pd.read_csv('D:\MERFISH\somatosensory_sc_labels.txt',header=None)
sc_labels.columns = ['celltype']
celltype = list(set(sc_labels.celltype))
celltype_dict = dict(zip([x+1 for x in range(len(celltype))],celltype))
metacell_dict = dict(zip([str(x+1) for x in range(len(celltype))],celltype))
sc_labels['cluster'] = [celltype.index(x)+1 for x in sc_labels.celltype]
sc_labels['barcode'] = sc_counts.index


RNA_data = sc_counts
Spatial_data = st_counts
locations = st_coords.loc[:, ['X', 'Y']]



RNA_data = RNA_data.T
Spatial_data = Spatial_data.drop('barcode', axis=1)
common_genes = set(Spatial_data.columns.values).intersection(RNA_data.index.values)
RNA_data = RNA_data.loc[common_genes, :]
Spatial_data = Spatial_data.loc[:, common_genes]

gene_names = np.array(RNA_data.index.values)
dge = RNA_data.values
dge = dge.T
num_cells = dge.shape[0]
print('number of cells and genes in the matrix:', dge.shape)

hvg = np.argsort(np.divide(np.var(dge, axis=0), np.mean(dge, axis=0) + 0.0001))
dge_hvg = dge[:, hvg[-2000:]]

num_locations = locations.shape[0]

p_location, p_expression = nc.rc.create_space_distributions(num_locations, num_cells)
cost_expression, cost_locations = nc.rc.setup_for_OT_reconstruction(dge_hvg, locations, num_neighbors_source=5,
                                                                    num_neighbors_target=5)

insitu_matrix = np.array(Spatial_data)
insitu_genes = np.array(Spatial_data.columns)

markers_in_sc = np.array([], dtype='int')
for marker in insitu_genes:
    marker_index = np.where(gene_names == marker)[0]
    if len(marker_index) > 0:
        markers_in_sc = np.append(markers_in_sc, marker_index[0])

cost_marker_genes = cdist(dge[:, markers_in_sc] / np.amax(dge[:, markers_in_sc]),
                          insitu_matrix / np.amax(insitu_matrix))
alpha_linear = 0.45
gw = nc.rc._GWadjusted.gromov_wasserstein_adjusted_norm(cost_marker_genes, cost_expression, cost_locations,
                                                        alpha_linear, p_expression, p_location, 'square_loss',
                                                        epsilon=5e-3, verbose=True)

df = pd.DataFrame(gw)
df.columns = st_coords.loc[:, 'barcode']
df = df / df.sum(axis=0)
df['class'] = sc_labels.loc[:, 'celltype']
df = df.groupby('class').sum().T
df.to_csv('D:\MERFISH/novoSpaRc/novoSpaRc_MERFISH_100.csv')

end = time.time()
print(end - start)
