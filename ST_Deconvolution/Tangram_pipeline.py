### Tangram usage guide
### system requirment: python >= 3.8.5
### the package installation refer to https://github.com/broadinstitute/Tangram

### Step-1: prepare conda environment for Tangram
conda env create -f environment.yml

### Step-2: install STRIDE package
conda activate tangram-env
pip install tangram-sc

### Step-3: deconvolution using Tangram
import pandas as pd
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import time
import psutil
import os

output_file_path = os.getcwd()
# add `Tangram` to path
import tangram as tg

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
adata_sc = sc.read_csv("/home/yll/cell2loction/mefish_singlecell_new_data.csv")
adata_sc.obs['sample'] = adata_sc.obs_names
adata_sc.var['SYMBOL'] = adata_sc.var_names
print(adata_sc)

ref_cell_metadata = pd.read_csv("/home/yll/cell2loction/mefish_singlecell_new_metadata.csv")
#embryo:final_cluster_label, brain:cell_type, merfish:Cell_class
celltype = ref_cell_metadata['Cell_class'].astype('category')
celltype.index = adata_sc.obs.index
print(celltype)
adata_sc.obs['CellType'] = celltype
sc.pp.normalize_total(adata_sc)

celltype_counts = adata_sc.obs['CellType'].value_counts()

sc.tl.rank_genes_groups(adata_sc, groupby="CellType", use_raw=False)
markers_df = pd.DataFrame(adata_sc.uns["rank_genes_groups"]["names"]).iloc[0:200, :]
#print(markers_df)

genes_sc = np.unique(markers_df.melt().value.values)
print(genes_sc)

#load spatial data
adata_st = sc.read_csv("/home/yll/cell2loction/merfish_spatial_data.csv")
#print(adata_vis)

adata_st.obs['sample'] = adata_st.obs_names
adata_st.var['SYMBOL'] = adata_st.var_names#.var_names views gene names
print(adata_st)

genes_st = adata_st.var_names.values
genes = list(set(genes_sc).intersection(set(genes_st)))

tg.pp_adatas(adata_sc, adata_st, genes=genes)

ad_map = tg.map_cells_to_space(
                   adata_sc,
                   adata_st,
                   mode='clusters',
                   cluster_label='CellType')

tg.project_cell_annotations(ad_map, adata_st, annotation='CellType')

end = show_info('end')
final = time.time()
print('Running time '+str(round((final-first)/60,4)) + 'min')
print('total memory used '+str(end-start) + 'GB')

celltype_density = adata_st.obsm['tangram_ct_pred']
celltype_density = (celltype_density.T/celltype_density.sum(axis=1)).T
celltype_density.to_csv(output_file_path+'/Tangram_result.csv',sep=',')

