# neonat-inflam-microglia

This repository houses the codes for the single-cell transcriptome analyses presented in the study **Neonatal inflammation impairs developmentally associated subsets of microglia and promotes a highly reactive microglial subset**. We used a in-house scRNA-seq pipeline targeting 46 genes to replicate the cell cluster proportions obtained by whole-transcriptome RNA sequencing on 19,061 independently obtained cells. A multinomial logistic model predicting the cluster membership for each cell was fitted to a training dataset containing data restricted to the gene panel for 95% of the cells from the whole transcriptome expression data set. For the remaining 5% of cells used for the whole-trancriptome data set, the fitted model had a good cluster membership prediction performance (AUC = 0.923 and ARI = 0.452). The code to fit the model, to test its performance and infer cluster annotation using datta restricted to the gene panel are presented in this directory.

## Dependencies

OS: Debian Linux.

```bash
python==3.6
numpy==1.19.5
pandas==1.1.5
scipy==1.5.4
scikit-learn==0.24.2
seaborn==0.11.2
mygene==3.2.2
umap-learn==0.5.3
umap-learn[plot]
Levenshtein==0.18.1
```

## Files

+ *params.py* (parameter file)
+ *cluster_alignment.py* (fits and saves the model to infer whole clusters from samples restricted to whole/targeted/facs panel genes)
+ *cluster_prediction.py* (predicts using the fitted model clusters on whole/targeted/facs data)
+ *cluster_barplot.py* (plots barplot of cell proportions in some clusters in whole vs whole/targeted/facs data)
+ *cluster_curves.py* (plots average and cluster-wise ROC and Precision-Recall curves)

## Installation 

Run (requires Conda and Pip)

```bash
conda create -n microglia python=3.6 -y
conda activate microglia
python3 -m pip install -r requirements.txt
```

## Reproduce the results shown in the paper

Download the following files from GEO

+ *info_experience.csv* (metadata on whole data)
+ *raw_counts.csv* (whole gene expression matrix on the targeted gene panel)
+ *raw_counts_variable_feature.csv* (whole gene expression matrix, restricted to the 2,000 most variable genes)
+ *metadata_targeted.csv* (metadata on targeted data)
+ *raw_counts_targeted.csv* (targeted gene expression matrix on the targeted gene panel)

and put them in a folder called DATA/ in this folder. Then run

```bash
bash dataprocess_commands.sh ## data processing from files
bash analysis_commands.sh    ## obtain the results
```

To obtain the figures shown in the paper, run the Jupyter notebook

```bash
jupyter notebook model_analysis.ipynb
```
