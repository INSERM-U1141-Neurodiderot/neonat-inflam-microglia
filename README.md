# neonat-inflam-microglia

<!-- [![DOI](https://zenodo.org/badge/696752316.svg)](https://zenodo.org/doi/10.5281/zenodo.TO CREATE) -->

This repository houses the codes for the single-cell transcriptome analyses presented in the study **Neonatal inflammation impairs developmentally associated subsets of microglia and promotes a highly reactive microglial subset**. The raw and processed data are available on GEO with accession number [GSE165113](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165113) and will be made publically available at the paper acceptance. 

The brain-resident myeloid cells, the microglia, and border-associated macrophages (BAMs) play crucial roles in both immunity and neurodevelopment. The disruption of microglial development trajectories by neonatal inflammatory challenges is an important issue in research on neurodevelopmental disorders (NDDs). We explored the impact of neonatal inflammation induced by interleukin-1 beta injections in a mouse NDD model on the heterogeneity of brain myeloid cell subsets. 

The repository **neonat-inflam-microglia** is divided into five folders: 

  - Seurat
  - Model predicting clusters 
  - DEG
  - SCPA
  - SCENIC

## Seurat
The package [Seurat v3.0](https://doi.org/10.1016/j.cell.2019.05.031) was used for the analyses of the single-cell RNA sequencing data of the study (47,211 cells):
- 1-[merge_data.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/1-merge_data.ipynb) for QC filtering and data merging
- 2-[microglia_preprocess.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/2-microglia_preprocess.ipynb) for preprocessing, integration and UMAP visualisation of data
- 3-[microglia_HeatMap.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/3-microglia_HeatMap.ipynb) to build the heatmap of expression levels for the top markergenes of CD11B+ cell clusters in PBS condition
- 4-[main_plot.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/4-main_plot.ipynb) to evaluate cell cluster proportions and for UMAP visualisation with cell cycle scores and SingleR cell annotation predictions
- 5-[microglia_align_ref_Li.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/5-microglia_align_ref_Li.ipynb) for UMAP visualisation with SingleR cell annotation predictions using the dataset of [Li et al. 2019](https://doi.org/10.1016/j.neuron.2018.12.006)
- 6-[microglia_C0subclustering.ipynb](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/1-Seurat/6-microglia_C0subclustering.ipynb) for the subclustering of the C0 cluster

## Model predicting clusters 
We used a in-house scRNA-seq pipeline targeting 46 genes to replicate the cell cluster proportions obtained by whole-transcriptome RNA sequencing on 19,061 independently obtained cells. A multinomial logistic model predicting the cluster membership for each cell was fitted to a training dataset containing data restricted to the gene panel for 95% of the cells from the whole transcriptome expression data set. For the remaining 5% of cells used for the whole-trancriptome data set, the fitted model had a good cluster membership prediction performance (AUC = 0.923 and ARI = 0.452).The code to analyse the fitted model is presented in a pdf [document](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/2-Model/Analyse_modele.pdf).

## DEG
Differential expression analyses were performed with the <i>FindAllMarkers</i> of the Seurat package. The results are presented in the supplementary table 3:
  - notebook to obtain the Supp Table [3a](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/3-DEG/supp_table_3a.ipynb) for IL1-beta induced differentially expressed genes in each cluster at each time points
  - notebook to obtain the Supp Table [3b](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/3-DEG/supp_table_3b_IL1.ipynb) for gene markers of clusters in the IL1-beta condition
  - notebook to obtain the Supp Table [3c](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/3-DEG/supp_table_3c_PBS.ipynb) for gene markers of clusters in the PBS condition

## SCPA
Functional enrichment was analyzed with the [SCPA](https://doi.org/10.1016/j.celrep.2022.111697) package, with Hallmark functional annotations from the [msigdb](https://doi.org/10.1016/j.cels.2015.12.004) database. Seurat objects were split between PBS and IL-1 and Seurat objects. SCPA for Hallmark functional annotations were extracted at each postnatal day for:
- each cluster ([notebook](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/4-SCPA/SCPA_fig3.ipynb) to obtain the heatmap of Figure 3A)
- each subcluster of the cluster C0 ([notebook](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/4-SCPA/SCPA_fig4.ipynb) to obtain the heatmap of Figure 4F)

## SCENIC
We investigated the gene regulatory networks of the subcluster C0-2 by using [SCENIC](https://doi.org/10.1038/nmeth.4463) tools to identify regulons. The results are presented in the Figure 4G. The [configuration file](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/5-SCENIC/microglia.vsn-pipelines.complete.config) to execute the pipeline in Nextflow and the [notebook](https://github.com/INSERM-U1141-Neurodiderot/neonat-inflam-microglia/blob/main/5-SCENIC/microglia_C0_2subcluster_scenic.ipynb) of the code used to obtain the heatmap are available.

## Cite

+ To cite this work please use the reference of the following paper:

[Neonatal inflammation impairs developmentally associated subsets of microglia and promotes a highly reactive microglial subset](preprintweblink to add)
<!-- + BibTeX citation: -->
