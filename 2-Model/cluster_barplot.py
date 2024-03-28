#coding: utf-8

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
import pandas as pd
import numpy as np
import seaborn as sns

from params import *

categories = all_clusters

def plot_stacked(data, title="percent_stacked_barchart", folder=""):
    barWidth = 0.85
    conditions = list(data.keys())
    categories = list(sorted([c for c in list(data[conditions[0]].columns) if (c != "Ncells")]))
    r = range(data[conditions[0]].shape[0])
    fig, axes = plt.subplots(figsize=(10*len(conditions), 8), nrows=1, ncols=len(conditions))
    if (len(conditions) == 1):
        axes = [axes]
    for ic, cond in enumerate(conditions):
        for ict, cat in enumerate(categories):
            color = palette[str(cat)][0]
            if (ict<1):
                axes[ic].bar(r, list(data[cond][cat]), color=color, edgecolor='white', width=barWidth)
            else:
                axes[ic].bar(r, list(data[cond][cat]), bottom=[sum(ls) for ls in zip(*[data[cond][categories[ct]] for ct in range(ict)])], color=color, edgecolor='white', width=barWidth)
        h = [Line2D([0], [0], color=palette[str(cat)][0], linewidth=5) for ict, cat in enumerate(categories)]
        axes[ic].legend(h, [palette[str(category)][1] for category in categories])
        if (pred_data_type=="facs"):
                axes[ic].set_xticklabels([0]+[y for if_, factor in enumerate(data[cond].index) for y in [""]*2+[factor+"\n(N="+str(data[cond].loc[factor]["Ncells"])+")"]+[""]*1 ])
        else:
                axes[ic].set_xticklabels([0]+[factor+"\n(N="+str(data[cond].loc[factor]["Ncells"])+")" for if_, factor in enumerate(data[cond].index)])
        axes[ic].set_xlabel("Age")
        axes[ic].set_title("Proportions of cells according to "+cond)
    plt.savefig(folder+title+".png", bbox_inches="tight")

#########################################################
## Cumulated histogram of proportions in each category ##
#########################################################

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
elif (pred_data_type=="targeted"):
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

factors = list(sorted(set(list([i.split("-")[0] for i in samples_targeted]))))
conditions = list(sorted(list(set([i.split("-")[1] for i in samples_targeted]))))

total_nb = lambda x : sum([np.sum(np.array(x)==cat) for cat in categories])
percentage = lambda x, cat : 100*sum(np.array(x)==cat)/total_nb(x)
get_subset = lambda cls, sp, s : cls.loc[sp==s].astype(str).values.flatten()

## Build data
data = {}
for cond in conditions:
    data_ = pd.DataFrame([], columns=categories+["Ncells"], index=[factor+"-"+type_ for factor in factors for type_ in ["whole", pred_data_type+" (model)"]])
    for factor in factors:
        cl_whole = get_subset(clustering_whole, samples_whole, factor+"-"+cond)
        cl_target = get_subset(clustering_targeted, samples_targeted, factor+"-"+cond)
        Nt, Nw = [total_nb(x) for x in [cl_target, cl_whole]]
        data_.loc[factor+"-"+pred_data_type+" (model)"] = [percentage(cl_target,category) for category in categories]+[Nt]
        data_.loc[factor+"-whole"] = [percentage(cl_whole,category) for category in categories]+[Nw]
    data.setdefault(cond, data_)

plot_stacked(data, title="percent_stacked_barchart", folder=rfolder)

for d in data:
    print(d)
    data[d].columns = [palette[str(cat)][1] for cat in categories]+["Ncells"]
    print(data[d])
    data[d].to_csv(rfolder+"percentages_cond="+d+".csv")
