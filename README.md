# neonat-inflam-microglia


This repository houses the codes for the single-cell transcriptome analyses presented in the study **Neonatal inflammation impairs developmentally associated subsets of microglia and promotes a highly reactive microglial subset**. The raw and processed data are available on GEO with accession number [GSE165113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165113) and will be made publically available at the paper acceptance. 

The repository **neonat-inflam-microglia** is divided into five folders: 

1- Seurat

2- Model predicting clusters 

3- DEG

4- SCPA

5- SCENIC


## Abstract

The brain-resident myeloid cells, the microglia, and border-associated macrophages (BAMs) play crucial roles in both immunity and neurodevelopment. Microglia regulate neurogenesis, synapse formation and myelination and their heterogeneity is tightly regulated across developmental stages. The disruption of microglial development trajectories by neonatal inflammatory challenges is an important issue in research on neurodevelopmental disorders (NDDs), as models have suggested a strong association between inflammation and the onset of cognitive defects. We explored the impact of neonatal inflammation induced by interleukin-1 beta injections in a mouse NDD model on the heterogeneity of brain myeloid cell subsets. Sparc-expressing microglia were hyperreactive and had high levels of the complement receptor C5AR1, whereas other microglial subsets, including the proliferative and Spp1-enriched microglia, but also BAMs, had lower levels of this receptor. These findings suggest major changes in microglial development trajectories potentially contributing to NDDs induced by neonatal inflammation, and indicate possible treatment strategies targeting microglia.

## 1- Seurat
The Jupyter notebooks used for the analyses of the single-cell RNA sequencing data of the study :
- [QC filtering and data merging](https://github.com/Goultard59/microglia_paper_2024/blob/main/Seurat/merge_data.ipynb)
- [preprocessing, integration and UMAP visualisation of data](https://github.com/Goultard59/microglia_paper_2024/blob/main/Seurat/microglia_preprocess.ipynb)
- [building the heatmap of expression levels for the top markergenes of CD11B+ cell clusters in PBS condition](https://github.com/Goultard59/microglia_paper_2024/blob/main/Seurat/microglia_HeatMap.ipynb)
- [evaluating cell cluster proportions and UMAP visualisation with cell cycle scores and SingleR cell annotation predictions](https://github.com/Goultard59/microglia_paper_2024/blob/main/Seurat/main_plot.ipynb)
- UMAP visualisation with SingleR cell annotation predictions using the Li et al dataset
- subclustering of the C0 cluster

  ## 2-

  ## 3- DEG
  The Jupyter notebooks used to build the supplementary table 4:
  - IL1-beta induced differentially expressed genes in each cluster at each time points
  - gene cell-cluster markers in the IL1-beta condition
  - gene cell-cluster markers in the PBS condition

  ## 4- SCPA
  The Jupyter notebook Fig3
  For the subclusters

  ## 5- SCENIC

