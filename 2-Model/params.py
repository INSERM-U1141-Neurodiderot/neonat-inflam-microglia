#coding: utf-8

import numpy
from scipy import stats
from multiprocessing import cpu_count
import numpy as np
import os

def rank_magnitude(x, y):
    def asym_rm(x,y):
        ids = np.argsort(y.flatten().tolist())
        ids_x = np.argsort(x.flatten().tolist()).tolist()
        y = y.flatten()[ids]
        min_rank = np.sum([v*(len(y)-vi+1) for vi, v in enumerate(y)])
        max_rank = np.sum([v*vi for vi, v in enumerate(y)])
        coef = 2.*np.sum([y[sidx]*ids_x.index(idx) for sidx, idx in enumerate(ids)])
        return (coef-max_rank-min_rank)/float(max_rank-min_rank)
    return (asym_rm(x,y)+asym_rm(y,x))*0.5

method = "LogisticRegression"
metric = "correlation"
Niter_CV=1 # number of times the cross-validation is run
Nfold_CV=5 # number of parts in cross-validation
Ntest_size=0.05 # test size proportion compared to training dataset
pred_data_type = ["whole","targeted"][1]
normalisation_type = ["standard", "standard-batch", "standard-training", "centered", "log-shifted", "censored"][0]
npermutations = 100
gene_set = ["full_genes", "2000mostvariable"][0] ## potentially consider all genes in whole data, or restrict to the 2,000 most variable markers
clusters_annot = "12clusters"
all_clusters = list(map(str,list(range(8)))) 
assert (np.max(list(map(int,all_clusters)))<=11 and np.min(list(map(int,all_clusters)))>=0)

np.random.seed(npermutations)

dist = metric if (metric not in ["rank_magnitude", "spearman"]) else (rank_magnitude if (metric=="rank_magnitude") else lambda x,y : stats.spearmanr(x,y)[0])
params = {"dist": dist, "shrink_threshold": None, "gamma": 20., "l1_ratio": 1, "metric": metric, "max_iter": 2000, "n_jobs": cpu_count()-2, "tol": 1e-4,"alpha": 1.0, "class_weight": None, "multiclass": "multinomial"}

conditions = ["PBS", "IL1"]
condition = conditions[1]
factors = ["P1", "P3", "P5"]
age = factors[-1]
batch = ["a", "b", "c", "d", "e", "f"][0]
sample_type = [age.upper()+"-"+condition.upper()+batch.lower(),age.upper()+"-"+condition.upper(),""][-1]

dfolder="data/" ## data folder
dfolder_files = {
	"metadata_whole": "microglia_info_experience.csv",
	"data_whole": "microglia_raw_counts.csv",
	"data_whole2000": "microglia_raw_counts_variable_feature.csv",
	"metadata_targeted": "microglia_metadata_targeted.csv",
	"data_targeted": "microglia_raw_counts_targeted.csv",
}
assert all([os.path.exists(dfolder+"/"+dfolder_files[f]) for f in dfolder_files])
rfolder="results_"+("_".join([pred_data_type,normalisation_type,sample_type,method]))+"/"

palette = {"0": ["#16a085", "C0"], #0 ieg
"1": ["#2980b9", "C1"], #1 macs
"2": ["#4834d4", "C2"], #2 hom
"3": ["#f39c12", "C3"], #3 div
"4": ["#c0392b", "C4"], #4 monocyte
"5": ["#130f40", "C5"], #5 cd11c
"6": ["#f78fb3", "C6"], #6 granulocyte
"7": ["#65d6ce", "C7"], #7 macs dividing
"8": ["#d1c145", "C8"], #8
"9": ["#2ecc71", "C9"], #9
"10": ["#d35400", "C10"], #10
"11": ["#c44569", "C11"], #11
"12": ['#ff7f50', "12"],
"13": ['#706fd3', "13"],
"14": ['#f9ca24', "14"],
"15": ['#34ace0', "15"],
"16": ['#33d9b2', "16"],
"17": ['#2c2c54', "17"],
"18": ['#be2edd', "18"]
}

