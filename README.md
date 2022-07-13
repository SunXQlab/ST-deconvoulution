Benchmarking and integration of methods for deconvoluting spatial transcriptomic data
===========================
Implementation descripttion
--------------------------
we collected three synthetic ST datasets with known single-cell compositions and a human heart ST dataset with known regional cell type information
to benchmark 14 different deconvlution methods. Furthermore, we investigate the robustness of different methods to sequencing depth, spot size, and 
the choice of normalization. Moreover, we propose a new ensemble learning-based deconvolution method (EnDecon) by integrating multiple individual 
methods for more accurate deconvolution.
The major new findings include: 
1) cell2loction, RCTD, and spatialDWLS are more accurate than other ST deconvolution methods, based on the evaluation of three metrics: RMSE, PCC, and JSD;<br> 
2) cell2location and spatialDWLS are more robust to the variation of sequencing depth than RCTD;<br>
3) the accuracy of the existing methods tends to decrease as the spot size becomes smaller;<br> 
4) most decon-volution methods perform best when they normalize ST data using the method described in their origi-nal papers;<br> 
5) the integrative method, EnDecon, could achieve more accurate ST deconvolution.<br>
--------------------------
Datasets
All datasets can be downloaded from their respective sources:<br>
* mouse embryo was downloaded from the National Center for Biotechnology Information (NCBI) under GSE166692<br>
* MPOA was downloaded from https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248/<br>
* mouse brain ST data was download from https://www.dropbox.com/s/azjysbt7lbpmbew/brain_st_cortex.rds?dl=0 <br>
* mouse brain scRNA-seq data was downloaded from https://www.dropbox.com/s/ruseq3necn176c7/brain_sc.rds?dl=0<br>
* human developing heart ST data was downloaded from https://data.mendeley.com/datasets/mbvhhf8m62/2<br>
---------------------------------
Dependencies and requirements for disentangle discrete cell types from spatial mixtures<br>
**R**<br>
RCTD(by spacexr of version 2.0.0)<br>
Giotto/PAGE/RANK/Hypergeometric(Version 1.1.0)<br>
spatialDWLS(by Giotto of Version 1.0.4)<br>
SPOTlight(Version 0.1.7)<br>
Seurat(Version 4.0.5)<br>
STdeconvolve(Version 0.1.0)<br>
EnDecon<br>
**Python**<br>
cell2location(Version 0.7a0)<br>
DestVI(scvi-tools Version 0.11.0)<br>
stereoscope(Version 0.3.1)<br>
STRIDE(Version 0.0.1b0)<br>
Tangram(Version 1.0.0)<br>





