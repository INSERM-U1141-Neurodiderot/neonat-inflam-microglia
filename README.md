# neonat-inflam-microglia


This repository houses the codes for the single-cell transcriptome analyses presented in the study **Neonatal inflammation impairs developmentally associated subsets of microglia and promotes a highly reactive microglial subset**. The raw and processed data are available on GEO with accession number [GSE165113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165113) and will be made publically available at the paper acceptance. 

The repository **neonat-inflam-microglia** is divided into five folders: 

  - Seurat
  - Model predicting clusters 
  - DEG
  - SCPA
  - SCENIC


## Abstract

The brain-resident myeloid cells, the microglia, and border-associated macrophages (BAMs) play crucial roles in both immunity and neurodevelopment. Microglia regulate neurogenesis, synapse formation and myelination and their heterogeneity is tightly regulated across developmental stages. The disruption of microglial development trajectories by neonatal inflammatory challenges is an important issue in research on neurodevelopmental disorders (NDDs), as models have suggested a strong association between inflammation and the onset of cognitive defects. We explored the impact of neonatal inflammation induced by interleukin-1 beta injections in a mouse NDD model on the heterogeneity of brain myeloid cell subsets. Sparc-expressing microglia were hyperreactive and had high levels of the complement receptor C5AR1, whereas other microglial subsets, including the proliferative and Spp1-enriched microglia, but also BAMs, had lower levels of this receptor. These findings suggest major changes in microglial development trajectories potentially contributing to NDDs induced by neonatal inflammation, and indicate possible treatment strategies targeting microglia.

## Seurat
The package Seurat was used for the analyses of the single-cell RNA sequencing data of the study:
- 1-[merge_data.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/1-merge_data.ipynb) for QC filtering and data merging
- 2-[microglia_preprocess.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/2-microglia_preprocess.ipynb) for preprocessing, integration and UMAP visualisation of data
- 3-[microglia_HeatMap.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/3-microglia_HeatMap.ipynb) to build the heatmap of expression levels for the top markergenes of CD11B+ cell clusters in PBS condition
- 4-[main_plot.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/4-main_plot.ipynb) to evaluate cell cluster proportions and for UMAP visualisation with cell cycle scores and SingleR cell annotation predictions
- 5-[microglia_align_ref_Li.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/5-microglia_align_ref_Li.ipynb) for UMAP visualisation with SingleR cell annotation predictions using the dataset of [Li et al. 2019](https://doi.org/10.1016/j.neuron.2018.12.006)
- 6-[microglia_C0subclustering.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/6-microglia_C0subclustering.ipynb) for the subclustering of the C0 cluster

## Model predicting clusters 

## DEG
Differential expression analyses were performed with the <i>FindAllMarkers<i> presented in the supplementary table 3:
  - notebook to obtain the Supp Table [3a](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/3-DEG/supp_table_3a.ipynb) for IL1-beta induced differentially expressed genes in each cluster at each time points
  - notebook to obtain the Supp Table [3b](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/3-DEG/supp_table_3b_IL1.ipynb) for gene markers of clusters in the IL1-beta condition
  - notebook to obtain the Supp Table [3c](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/3-DEG/supp_table_3c_PBS.ipynb) for gene markers of clusters in the PBS condition

## SCPA
  The Jupyter notebook Fig3
  For the subclusters

## SCENIC


