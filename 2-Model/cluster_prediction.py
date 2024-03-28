#coding:utf-8

import subprocess as sb

import random
random.seed(0)
import numpy as np
import pandas as pd
import scipy
import pickle as pkl
import os
import seaborn as sns

###############################
## PARAMETERS & METHODS      ##
###############################

from params import *
from cluster_alignment import data_normalisation

with open(rfolder+"gene_list.pck", "rb") as f:
    di = pkl.load(f)
targeted_panel, whole_panel = di[pred_data_type+"_panel"], di["whole_panel"]

if (pred_data_type=="whole"):
    if (gene_set == "2000mostvariable"): ## restrict to the 2,000 most variable
        info_experience = pd.read_csv(dfolder+dfolder_files["metadata_whole"], index_col=0)
        whole_fname = dfolder+dfolder_files["data_whole2000"]
    elif (gene_set == "full_genes"): ## do not restrict the input set of genes
        info_experience = pd.read_csv(dfolder+dfolder_files["metadata_whole"], index_col=0)
        whole_fname = dfolder+dfolder_files["data_whole"]
    else:
        raise ValueError
    counts = pd.read_csv(rfolder+"raw_counts_subset_feature.csv", index_col=0).loc[whole_panel].dropna().T
elif (pred_data_type=="targeted"):
    metadata = pd.read_csv(dfolder+dfolder_files["metadata_targeted"], index_col=0)
    counts = pd.read_csv(dfolder+dfolder_files["data_targeted"], index_col=0)
    counts = counts[["X"+(".".join(m.split("#"))) for m in metadata.index]].loc[targeted_panel].fillna(0.).T
else:
    raise ValueError

whole_matrix = pd.read_csv(rfolder+"raw_counts_subset_feature.csv", index_col=0).loc[whole_panel].dropna().T
_, X_pred_fit = data_normalisation(whole_matrix, normalisation_type, counts)

model_fname = rfolder+method+"_model.pck"
with open(model_fname, "rb") as f:
    max_auc_test_di = pkl.load(f)
clf = max_auc_test_di["estimator"]
print(clf) ## fitted estimator

Y_targeted_pred = clf.predict(X_pred_fit).ravel()
pd.DataFrame([Y_targeted_pred], columns=counts.index, index=["Y_targeted_pred"]).T.to_csv(rfolder+"Y_targeted_pred.csv")
