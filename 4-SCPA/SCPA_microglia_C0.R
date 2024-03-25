library(ggplot2)
library(cowplot)
library(dplyr)
library(SCPA)
library(tidyverse)
library(msigdbr)
library(tibble)
library(AnnotationHub)
library(patchwork)
library(Seurat)
library(rliger)
library(SeuratData)
library(magrittr)
library(SeuratWrappers)
library(SeuratDisk)
library(viridis)
library(tidyr)

immune.integrated <- readRDS("/home/adufour/work/rds_storage/microglia/microglia_IEG.rds")

tb1 <- msigdbr("Mus musculus", category = "C2", subcategory = "CP:KEGG")
tb2 <- msigdbr("Mus musculus", category = "C5", subcategory = "GO:BP")
tb3 <- msigdbr("Mus musculus", category = "H")

pathways <- rbind(tb1, tb2, tb3) %>% format_pathways()

cell_types <- unique(immune.integrated@meta.data$seurat_clusters)

scpa_out <- list()
bl_bm <- list(); bl_ln <- list()
for (i in cell_types) {
  cluster_0 <- seurat_extract(immune.integrated, meta1 = "seurat_clusters", value_meta1 = 0)
  cluster_1 <- seurat_extract(immune.integrated, meta1 = "seurat_clusters", value_meta1 = 1)
  cluster_2 <- seurat_extract(immune.integrated, meta1 = "seurat_clusters", value_meta1 = 2)

  bl_bm[[i]] <- compare_pathways(list(cluster_0, cluster_1), pathways, parallel = TRUE, core = 7)
  bl_ln[[i]] <- compare_pathways(list(cluster_0, cluster_2), pathways, parallel = TRUE, core = 7)
}

# For faster analysis with parallel processing, use 'parallel = TRUE' and 'cores = x' arguments

saveRDS(scpa_out, file = "/home/adufour/work/rds_storage/microglia_scpa_ieg_0.rds")
saveRDS(scpa_out, file = "/home/adufour/work/rds_storage/microglia_scpa_ieg_1.rds")