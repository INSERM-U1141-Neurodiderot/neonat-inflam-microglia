#coding:utf-8

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import seaborn as sns
import numpy as np
from scipy.stats import chisquare, ks_2samp#, kstest
from params import *

import pickle as pkl
import subprocess as sb
import pandas as pd

np.random.seed(0)
import random
random.seed(0)

metadata_fname = dfolder+"info_experience.csv"
info_experience = pd.read_csv(metadata_fname, index_col=0)
metadata_fname_ = dfolder+"metadata_microglia.csv"
info_experience_ = pd.read_csv(metadata_fname_, index_col=0)
if (gene_set == "2000mostvariable"): ## restrict to the 2,000 most variable
	## Clustering using Seurat on whole transcriptome data
	clustering_whole = info_experience["seurat_clusters"].astype(str)
elif (gene_set == "full_genes"): ## do not restrict the input set of genes
	## Clustering using Seurat on whole transcriptome data
	clustering_whole = info_experience.loc[info_experience_.index]["seurat_clusters"].astype(str)
else:
	raise ValueError

clustering_targeted = pd.read_csv(rfolder+"Y_targeted_pred.csv", index_col=0)[["Y_targeted_pred"]].astype(str)
if (clusters_annot=="new_clusters"):
	assert len(sample_type)==0
	subclusters = pd.read_csv(dfolder+"metadata_microglia.csv", index_col=0)["new_clusters"]
	clustering_whole[subclusters.index] = list(subclusters)
	clustering_whole = clustering_whole[subclusters.index]
	metadata_whole = pd.read_csv(dfolder+"metadata_microglia.csv", index_col=0) 
else:
	assert clusters_annot=="12clusters"
	metadata_whole = info_experience.loc[info_experience_.index]
	metadata_whole = metadata_whole[~metadata_whole.index.duplicated()]
samples_whole = np.array([x[:-1] for x in metadata_whole["SAMP"]]) 
clustering_targeted = clustering_targeted.loc[~clustering_targeted.index.duplicated()]
clustering_targeted = clustering_targeted["Y_targeted_pred"]

## <AGE>-<COND> samples
if (pred_data_type=="whole"):
	samples_targeted = samples_whole
	#samples_targeted = [metadata_whole.loc[x]["PND"].upper()+"-"+metadata_whole.loc[x]["STIM"].upper() for x in metadata_whole.index] ## TODO
elif (pred_data_type=="targeted"):
	#metadata = pd.read_csv(dfolder+"metadata_targeted.csv", index_col=0)
	metadata = pd.read_csv(dfolder+"final_metadata_targeted.csv", index_col=0)
	metadata = metadata.loc[~metadata.index.duplicated()]
	samples_targeted = [metadata.loc[x]["PND"].upper()+"-"+metadata.loc[x]["STIM"].upper() for x in metadata.index]
elif (pred_data_type=="facs"):
	samples_targeted = [i.split("_")[0] for i in clustering_targeted.index]
else:
	raise ValueError
samples_targeted = np.array(samples_targeted)

assert samples_targeted.shape==clustering_targeted.shape
assert samples_whole.shape==clustering_whole.shape

in_clusters = np.vectorize(lambda x : x in all_clusters)
ids_whole = in_clusters(clustering_whole)
ids_targeted = in_clusters(clustering_targeted)
clustering_whole = clustering_whole.loc[ids_whole]
clustering_targeted = clustering_targeted.loc[ids_targeted]
samples_targeted = np.array(samples_targeted)[ids_targeted.flatten()]
samples_whole = samples_whole[ids_whole.flatten()]

ages = list(sorted(list(set([i.split("-")[0] for i in samples_targeted]))))
conditions = list(sorted(list(set([i.split("-")[1] for i in samples_targeted]))))
get_subset = lambda cls, sp, s : cls.loc[sp==s].astype(str).values.flatten()
get_proportion = lambda cls, sp, s, c : 100*np.mean(get_subset(cls, sp, s)==c)
df = {factor: {cond: pd.DataFrame(
	{
	"predicted": {"Cluster %s" % c : get_proportion(clustering_targeted,samples_targeted,age+"-"+cond,c) for c in all_clusters},
	"whole": {"Cluster %s" % c : get_proportion(clustering_whole,samples_whole,age+"-"+cond,c) for c in all_clusters},
	})
	for cond in conditions} for factor in ages}

## Prints out
#for factor in list(sorted(df.keys())):
#	for cond in list(sorted(df[list(df.keys())[0]].keys())):
#		print({factor+"-"+cond: df[factor][cond]})

## Statistical tests
res_di = {}
for factor in list(sorted(df.keys())):
	for cond in list(sorted(df[list(df.keys())[0]].keys())):
		## Kolmogorov-Smirnov test
		r1 = get_subset(clustering_targeted, samples_targeted, factor+"-"+cond)
		r2 = get_subset(clustering_whole, samples_whole, factor+"-"+cond)
		res_ks = ks_2samp(r1, r2, alternative="two-sided")
		#p_ = []
		#for _ in range(npermutations):
		#	np.random.seed(int(1e6))
		#	np.random.shuffle(r1)
		#	np.random.shuffle(r2)
		#	res_ks_ = ks_2samp(r1, r2, alternative="two-sided")
		#	p_.append(res_ks_.pvalue)
		#print((np.mean(p_),np.var(p_)))
		## Chi2 test
		n_obs = get_subset(clustering_targeted, samples_targeted, factor+"-"+cond).shape[0]
		ids = ["Cluster %s" % c for c in all_clusters]
		f_obs = np.array(df[factor][cond].loc[ids]["predicted"])
		f_obs /= np.round(np.sum(f_obs), 7)/100
		f_exp = np.array(df[factor][cond].loc[ids]["whole"])
		f_exp /= np.round(np.sum(f_exp), 7)/100
		f_obs = f_obs*n_obs/100
		f_exp = f_exp*n_obs/100
		if not ((np.min(f_obs)>=5) and (np.min(f_exp)>=5) and (n_obs>=13)): # filter clusters
			ids = [i for ii, i in enumerate(all_clusters) if (f_obs[ii]<5)]
			f_obs = np.array([f_obs[ii]/n_obs for ii, i in enumerate(all_clusters) if (i not in ids)])
			f_obs = np.array([x/np.sum(f_obs)*n_obs for x in f_obs])
			f_exp = np.array([f_exp[ii]/n_obs for ii, i in enumerate(all_clusters) if (i not in ids)])
			f_exp = np.array([x/np.sum(f_exp)*n_obs for x in f_exp])
			print("Missing clusters %s %s %s" % (factor, cond, str(ids)))
		assert np.isclose(np.sum(f_obs),np.sum(f_exp)) ## ensure same sum of freq.
		res = chisquare(f_obs=f_obs, f_exp=f_exp)
		#res_ks2 = kstest(r1, r2) ## same result as in ks_2samp
		res_di.setdefault(factor+"-"+cond, {"chisq": res.statistic, "pvalue(chisq)": res.pvalue, "KSstat": res_ks.statistic, "pvalue(KSstat)": res_ks.pvalue})
print("The null hypothesis is that the proportions of clusters in %s data are equal to the proportions of clusters in whole sequencing data" % pred_data_type)
print(pd.DataFrame(res_di))
pd.DataFrame(res_di).to_csv(rfolder+"tests.csv")