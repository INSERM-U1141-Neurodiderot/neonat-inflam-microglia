#coding:utf-8

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import seaborn as sns
import numpy as np

from params import *

import pickle as pkl
import subprocess as sb
import pandas as pd

np.random.seed(0)
import random
random.seed(0)

with open(rfolder+method+"_model.pck", "rb") as f:
	model = pkl.load(f)
print("%s: %f" % ("AUC (test)", model["AUC (test)"]))
print("%s: %f" % ("ARI (test)", model["ARI (test)"]))
fpr, tpr, pre, rec = model["fpr"], model["tpr"], model["pre"], model["rec"]
fprs = pd.DataFrame({i: np.array(list(fpr[i])+[1]*(3-len(fpr[i]))) for i in fpr}).values.T
tprs = pd.DataFrame({i: np.array(list(tpr[i])+[1]*(3-len(tpr[i]))) for i in tpr}).values.T
pres = pd.DataFrame({i: np.array(list(pre[i])+[1]*(3-len(pre[i]))) for i in pre}).values.T
recs = pd.DataFrame({i: np.array(list(rec[i])+[0]*(3-len(rec[i]))) for i in rec}).values.T

## Compute average values across users
mean_fprs = fprs.mean(axis=0)
mean_pres = pres.mean(axis=0)
mean_recs = recs.mean(axis=0)
mean_tprs = tprs.mean(axis=0)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
## ROC curve and Precision-Recall curve
for i in all_clusters:
	axes[0].plot(fpr[i], tpr[i], 'r-', alpha = 0.2)
	axes[1].plot(pre[i], rec[i], 'b-', alpha = 0.2)
axes[0].plot(mean_fprs, mean_tprs, 'r-', alpha = 0.8, label="Model")
axes[1].plot(mean_pres, mean_recs, 'b-', alpha = 0.8, label="Model")
axes[0].plot([0, 1], [0, 1], linestyle = '--', lw = 2, color = 'r', alpha= 0.8, label="Constant")
axes[1].plot([0, 1], [1, 0], linestyle = '--', lw = 2, color = 'b', alpha= 0.8, label="Constant")
axes[0].set_ylabel('True Positive Rate')
axes[0].set_xlabel('False Positive Rate')
axes[0].legend(loc="lower right")
axes[0].set_title('Avg. cluster ROC curve (ARI test = %f)' % model["ARI (test)"])
axes[1].set_ylabel('Precision')
axes[1].set_xlabel('Recall')
axes[1].legend(loc="lower left")
axes[1].set_title('Avg. cluster PR curve (AUC test = %f)' % model["AUC (test)"])
plt.savefig(rfolder+"curves.png", bbox_inches="tight")
