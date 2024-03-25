#coding:utf-8

import subprocess as sb

#first_time=False

## Modules
#if (first_time):
#    pyversion = sb.check_output("python3 --version", shell=True).decode("utf-8") 
#    assert ".".join(pyversion.split(".")[:2]) == "Python 3.6"
#    sb.call("python -m pip install -U pip", shell=True) #newest version of pip
#    pipversion = sb.check_output("echo $(python3 -m pip --version) |  cut -d" " -f2", shell=True).decode("utf-8")
#    assert int(pipversion.split(".")) >= 19
#    sb.call("python3 -m pip install setuptools==51.0.0 numpy==1.16.6 pandas==0.24.2 scipy==1.2.3 scikit-learn==0.24.0 seaborn==0.9.1 mygene==3.1.0 umap-learn==0.4.6 umap-learn[plot] Levenshtein==0.18.1", shell=True)
#    #sb.call("python3 -m pip install --no-binary :all: nmslib", shell=True)

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

#print(metadata)
if (pred_data_type=="whole"):
    if (gene_set == "2000mostvariable"): ## restrict to the 2,000 most variable
        info_experience = pd.read_csv(dfolder+"info_experience.csv", index_col=0)
        whole_fname = dfolder+"raw_counts_variable_feature.csv"
    elif (gene_set == "full_genes"): ## do not restrict the input set of genes
        info_experience = pd.read_csv(dfolder+"metadata_microglia.csv", index_col=0)
        whole_fname = dfolder+"counts_microglia.csv"
    else:
        raise ValueError
    counts = pd.read_csv(rfolder+"raw_counts_subset_feature.csv", index_col=0)
    counts = counts[info_experience.index].loc[whole_panel].dropna().T
elif (pred_data_type=="targeted"):
    #metadata = pd.read_csv(dfolder+"metadata_targeted.csv", index_col=0)
    #counts = pd.read_csv(dfolder+"raw_counts_targeted.csv", index_col=0)
    metadata = pd.read_csv(dfolder+"final_metadata_targeted.csv", index_col=0)
    counts = pd.read_csv(dfolder+"final_targeted.csv", index_col=0).T
    #counts.columns = ["#".join(c[1:].split(".")) for c in counts.columns]
    counts = counts[metadata.index].loc[targeted_panel].fillna(0.).T
elif (pred_data_type=="facs"):
    panel = [matchings_WHOLE2FACS[g].upper() if (matchings_WHOLE2FACS[g]!="Ki67") else matchings_WHOLE2FACS[g] for g in whole_panel]
    targeted_matrix_IL1 = pd.read_csv(dfolder+"export_IL1_CD11B_MATRIX_concat_IL1_CD11B_v3_1.csv").T.loc[panel].dropna().T
    targeted_matrix_PBS = pd.read_csv(dfolder+"export_PBS_CD11B_MATRIX_concat_CD11B_PBS_v3_1.csv").T.loc[panel].dropna().T
    targeted_matrix_IL1.index = ["P5-IL1_%d" % (i+1) for i in range(targeted_matrix_IL1.shape[0])]
    targeted_matrix_PBS.index = ["P5-PBS_%d" % (i+1) for i in range(targeted_matrix_PBS.shape[0])]
    counts = pd.concat((targeted_matrix_IL1,targeted_matrix_PBS), axis=0)
    #targeted_matrix = pd.read_csv("raw_counts.csv").loc[targeted_panel].dropna().T
    #Nsamples_targeted, N_features = np.shape(targeted_matrix)
else:
    raise ValueError
#targeted_matrix = counts.loc[targeted_panel].dropna().T

#whole_matrix = pd.read_csv(rfolder+"raw_counts_variable_feature.csv", index_col=0).loc[whole_panel].dropna().T
whole_matrix = pd.read_csv(rfolder+"raw_counts_subset_feature.csv", index_col=0).loc[whole_panel].dropna().T
_, X_pred_fit = data_normalisation(whole_matrix, normalisation_type, counts)

model_fname = rfolder+method+"_model.pck"
with open(model_fname, "rb") as f:
    #clf = pkl.load(f)
    max_auc_test_di = pkl.load(f)
clf = max_auc_test_di["estimator"]
print(clf) ## fitted estimator

Y_targeted_pred = clf.predict(X_pred_fit).ravel()
pd.DataFrame([Y_targeted_pred], columns=counts.index, index=["Y_targeted_pred"]).T.to_csv(rfolder+"Y_targeted_pred.csv")
