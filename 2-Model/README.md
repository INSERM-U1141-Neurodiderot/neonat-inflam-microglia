Python 3.6

Necessary code:
+ *params.py* (parameter file)
+ *cluster_alignment.py* (fits and saves the model to infer whole clusters from samples restricted to whole/targeted/facs panel genes)
+ *cluster_prediction.py* (predicts using the fitted model clusters on whole/targeted/facs data)
+ *cluster_statistics.py* (tests whether distributions of cluster sizes is the same in whole vs whole/targeted/facs data using predicted clusters)
+ *cluster_barplot.py* (plots barplot of cell proportions in some clusters in whole vs whole/targeted/facs data)
+ *cluster_curves.py* (plots average and cluster-wise ROC and Precision-Recall curves)

Necessary files: (in DATA/)
+ *info_experience.csv* (metadata on whole data)
+ *raw_counts.csv* (whole gene expression matrix on the targeted gene panel)
+ *raw_counts_variable_feature.csv* (whole gene expression matrix, restricted to the 2,000 most variable genes)
+ *metadata_targeted.csv* (metadata on targeted data)
+ *raw_counts_targeted.csv* (targeted gene expression matrix on the targeted gene panel)
+ *export_PBS_CD11B_MATRIX_concat_CD11B_PBS_v3_1.csv* (FACS data on P5-PBS cells)
+ *export_IL1_CD11B_MATRIX_concat_IL1_CD11B_v3_1.csv* (FACS data on P5-IL1 cells)

Execute

```bash
conda create -n microglia python=3.6 -y
conda activate microglia
python3 -m pip install -r requirements.txt
bash dataprocess_commands.sh ## data processing from raw files
bash analysis_commands.sh ## obtain the results
```