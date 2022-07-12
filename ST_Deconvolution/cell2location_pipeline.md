import sys
import datetime
IN_COLAB = "google.colab" in sys.modules
if IN_COLAB and branch == "stable":
    !pip install --quiet git+https://github.com/BayraktarLab/cell2location#egg=cell2location[tutorials]

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import psutil
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

import cell2location
import scvi

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

results_folder = './results/merfish_analysis/'
#read spatial data
#import time
start_time = datetime.datetime.now()
#start = time.clock()
adata_vis = sc.read_csv("/home/yll/cell2loction/merfish_spatial_countdata.csv")
print(adata_vis)

adata_vis.obs['sample'] = adata_vis.obs_names
adata_vis.var['SYMBOL'] = adata_vis.var_names#.var_names views gene names
print(adata_vis)

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'
print(ref_run_name)
print(run_name)
os.getcwd() #check the current path
output_file_path = os.getcwd()
#read singlecell data
adata_ref = sc.read_csv("/home/yll/cell2loction/mefish_cell_countdata.csv")
print(adata_ref)

ref_cell_metadata = pd.read_csv("/home/yll/cell2loction/mefish_cell_metadata.csv")
celltype = ref_cell_metadata['Cell_class'].astype('category')
celltype.index = adata_ref.obs.index
#print(celltype)

adata_ref.obs['CellType'] = celltype
adata_ref.obs['sample'] = adata_ref.obs_names
adata_ref.var['SYMBOL'] = adata_ref.var.index
print(adata_ref)

# prepare anndata for the regression model
scvi.data.setup_anndata(adata=adata_ref,
                        # cell type, covariate used for constructing signatures
                        labels_key = 'CellType',
                       )
scvi.data.view_anndata_setup(adata_ref)

# create and train the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=False)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
scvi.data.setup_anndata(adata=adata_vis)
scvi.data.view_anndata_setup(adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection (using default here):
    detection_alpha=200
)

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=False)

pid = os.getpid()
p = psutil.Process(pid)
info = p.memory_full_info()
sys_memory1 = info.uss/ 1024 / 1024 / 1024
print("memory used: %f GB" % sys_memory1)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

# Save model
mod.save(f"{run_name}", overwrite=True)
# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)
end_time = datetime.datetime.now()
print("start time is",start_time)
print("end time is",end_time)

import os 
pid = os.getpid()
p = psutil.Process(pid)
info = p.memory_full_info()
sys_memory2 = info.uss/ 1024 / 1024 / 1024
print("memory used: %f GB" % sys_memory2)
# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)

print(adata_vis)
print(adata_ref)
means_cell_abundance = adata_vis.obsm['means_cell_abundance_w_sf']
q05_cell_abundance = adata_vis.obsm['q05_cell_abundance_w_sf']

#means_cell_abundance
means_cell_abundance.to_csv(output_file_path+'/cell2location_means_result.csv',sep=',')
q05_cell_abundance.to_csv(output_file_path + '/cell2location_q05_result.csv',sep=',')
