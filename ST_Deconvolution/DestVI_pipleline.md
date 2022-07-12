import scanpy as sc
import matplotlib.pyplot as plt
from scvi.model import CondSCVI, DestVI
import numpy as np
import pandas as pd
import scvi
import psutil
import time
import os
os.getcwd()
output_file_path = os.getcwd()
def show_info(start):
    pid = os.getpid()
    #模块名比较容易理解：获得当前进程的pid
    p = psutil.Process(pid)
    #根据pid找到进程，进而找到占用的内存值
    info = p.memory_full_info()
    memory = info.uss/1024/1024/1024
    return memory

first = time.time()
start = show_info('strat')

#load single cell dataset
sc_adata = sc.read_csv("/home/yll/cell2loction/mefish_cell_countdata.csv")

sc_adata.obs['sample'] = sc_adata.obs_names
sc_adata.var['SYMBOL'] = sc_adata.var_names
print(sc_adata)

ref_cell_metadata = pd.read_csv("/home/yll/cell2loction/mefish_cell_metadata.csv")
#print(ref_cell_metadata)

#embryo:final_cluster_label
celltype = ref_cell_metadata['Cell_class'].astype('category')
celltype.index = sc_adata.obs.index
print(celltype)
sc_adata.obs['CellType'] = celltype

sc_adata.layers["counts"] = sc_adata.X.copy()

sc.pp.normalize_total(sc_adata, target_sum=10e4)
sc.pp.log1p(sc_adata)
sc_adata.raw = sc_adata

st_adata = sc.read_csv("/home/yll/cell2loction/merfish_spatial_countdata.csv")
#print(adata_vis)

st_adata.obs['sample'] = st_adata.obs_names
st_adata.var['SYMBOL'] = st_adata.var_names#.var_names views gene names
print(st_adata)

# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
emb_sc_adata = emb_sc_adata[:, intersect].copy()
G = len(intersect)

scvi.data.setup_anndata(sc_adata, layer="counts", labels_key="CellType")

scvi.model.CondSCVI.setup_anndata(sc_adata,labels_key="CellType")
sc_model = scvi.model.CondSCVI(sc_adata)

#deconvolution 
st_adata.layers["counts"] = st_adata.X.copy()

# get dataset ready
scvi.data.setup_anndata(st_adata, layer="counts")

DestVI.setup_anndata(st_adata)
spatial_model = DestVI.from_rna_model(st_adata, sc_model)
spatial_model.train(max_epochs=2000)
st_adata.obsm["proportions"] = spatial_model.get_proportions(st_adata)

end = show_info('end')
final = time.time()
print('Running time '+str(round((final-first)/60,4)) + 'min')
print('total memory used '+str(end-start) + 'GB')

decon_mtrx = st_adata.obsm["proportions"]
decon_mtrx.to_csv(output_file_path +'/DestVI_result.csv',sep=',')