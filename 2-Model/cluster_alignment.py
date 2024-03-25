# coding: utf-8

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

from scipy.optimize import minimize
from sklearn.metrics import pairwise_distances as pdist

from scipy import stats
from sklearn import neighbors
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
from multiprocessing import cpu_count
#https://scikit-learn.org/stable/modules/multiclass
#Can't use any method relying on pairwise distances, because #samples is too large (kNN, SVM, ...)
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import cross_validate
from sklearn.metrics import roc_curve, adjusted_rand_score, precision_recall_curve

## Train to get the whole data clustering on panel-restricted whole count data

def data_normalisation(wm, type_, tm=None):#tm, type_):
    assert type_ in ["standard", "standard-batch", "standard-training", "centered", "log-shifted", "censored"]
    X = np.matrix(wm.values, dtype=np.float64)
    if (tm is not None):
        X_pred = np.matrix(tm.values, dtype=np.float64)
    if (type_=="standard"):
        X_fit = StandardScaler().fit_transform(X)
        if (tm is None):
            X_pred_fit = None
        else:
            X_pred_fit = StandardScaler().fit_transform(X_pred)
    elif (type_ == "standard-batch"):
        X_fit = StandardScaler().fit_transform(X)
        if (tm is None):
            X_pred_fit = None
        else:
            assert pred_data_type=="targeted"
            metadata_stdb = pd.read_csv(dfolder+"metadata_targeted.csv", index_col=0) 
            X_pred_fit = np.empty(X_pred.shape) ##
            for s in list(set(list(metadata_stdb["SAMP"]))):
                ids = np.argwhere(metadata_stdb[["SAMP"]].values.flatten()==s).flatten().tolist() 
                X_pred_fit[ids,:] =  StandardScaler().fit_transform(X_pred[ids,:]) 
    elif (type_ == "standard-training"):
        scaler = StandardScaler()
        scaler.fit(X)
        X_fit = scaler.transform(X)
        if (tm is None):
            X_pred_fit = None
        else:
            X_pred_fit = scaler.transform(X_pred)
    ## Ahlmann-Eltze, C., Huber, W. Comparison of transformations for single-cell RNA-seq data. Nat Methods 20, 665–672 (2023). https://doi.org/10.1038/s41592-023-01814-1
    elif (type_ == "log-shifted"):
        s_train = np.sum(X, axis=0) ## size factor for each cell
        L = float(np.sum(s_train)/np.prod(s_train.shape))
        s_train /= L
        y_0 = 1
        s_test = np.sum(X_pred, axis=0)/L ## size factor for each cell
        X_fit = np.log(X/s_train+y_0)
        if (tm is None):
            X_pred_fit = None
        else:
            X_pred_fit = np.log(X_pred/s_test+y_0)
    elif (type_ == "centered"):
        nb_UMI_per_cell = np.mean(X, axis=1)/(1+np.sum(X, axis=1)) 
        X_fit = np.multiply(nb_UMI_per_cell, X)
        if (tm is None):
            X_pred_fit = None
        else:
            order_in_tm = [list(tm.index).index(c) for c in list(wm.index)]
            X_pred_fit = np.multiply(nb_UMI_per_cell[order_in_tm], X_pred)
    elif (type_ == "censored"):
        thres=0.90 # don't normalize expr. values of genes which expr. is larger than the 100*thres th percentile of expr in a given cell
        percs=np.percentile(X, int(100*thres), axis=1)
        means=np.mean(X, axis=1)
        X_fit = np.empty(X.shape) 
        for cell in range(X.shape[0]):
            p = float(percs[cell])
            m = float(means[cell])
            ids = np.argwhere(X[cell,:] < p)
            X_fit[cell, ids] = X[cell, ids]/m
        if (tm is None):
            X_pred_fit = None
        else:
            X_pred_fit = np.empty(X_pred.shape)                
            order_in_tm = [list(tm.index).index(c) for c in list(wm.index)]
            for cell in range(X.shape[0]):
                p = float(percs[cell])
                m = float(means[cell])
                ids_pred = np.argwhere(X[order_in_tm[cell],:] < p)
                X_pred_fit[order_in_tm[cell], ids_pred] = X_pred[order_in_tm[cell], ids_pred]/m
    return X_fit, X_pred_fit

if __name__ == "__main__":

    from params import *

    #####################################
    ## DATA IMPORT                     ##
    #####################################

    if (gene_set == "2000mostvariable"): ## restrict to the 2,000 most variable
        whole_fname = dfolder+"raw_counts_variable_feature.csv"
    elif (gene_set == "full_genes"): ## do not restrict the input set of genes
        whole_fname = dfolder+"counts_microglia.csv"
    else:
        raise ValueError
    metadata_fname = dfolder+"info_experience.csv"
    panel_genes = {
        "whole": list(pd.read_csv(dfolder+"final_targeted.csv", index_col=0).T.index), 
        ## list(pd.read_csv(dfolder+"raw_counts.csv", index_col=0).index), 
        "targeted": list(pd.read_csv(dfolder+"final_targeted.csv", index_col=0).T.index), 
        ## list(pd.read_csv(dfolder+"raw_counts_targeted.csv", index_col=0).index),
        "facs": list(matchings_WHOLE2FACS.keys())
        }[pred_data_type]

    sb.call("mkdir -p "+rfolder,shell=True)

    print("Selection of samples..."),

    ###### GET SAMPLES
    ## Select samples according to condition and age
    info_experience = pd.read_csv(metadata_fname, index_col=0)
    if (len(sample_type) > 0):
        metadata = info_experience.loc[[idx for idx in info_experience.index if (sample_type in info_experience.loc[idx]["SAMP"])]]
        if (condition.upper() in sample_type):
            assert np.all(metadata["STIM"] == condition.lower())
        if (age.upper() in sample_type):
            assert np.all(metadata["PND"] == age.upper())
        if (batch.lower() in sample_type):
            assert np.all([samp[-1] == batch.lower() for samp in metadata["SAMP"]])
    else:
        metadata = info_experience.copy()
    samples = list(metadata.index)
    ## Clustering using Seurat on whole transcriptome data
    whole_clustering = info_experience.loc[samples]["seurat_clusters"].astype(str)
    if (clusters_annot=="new_clusters"):
        assert len(sample_type)==0
        subclusters = pd.read_csv(dfolder+"metadata_microglia.csv", index_col=0)["new_clusters"].to_dict()
        whole_clustering = whole_clustering.to_dict()
        whole_clustering.update(subclusters)
        whole_clustering = pd.Series(whole_clustering).loc[samples]
    else:
        assert clusters_annot=="12clusters"
    #whole_genes = list(pd.read_csv(whole_fname, index_col=0).index)
    whole_genes = sb.check_output("cut -d',' -f1 %s" % whole_fname, shell=True).decode("utf-8").split("\n")[:-1]
    whole_genes = [x[1:-1] for x in whole_genes if (len(x[1:-1])>0)]
    #targeted_clustering = info_experience["Seurat_clusters_targeted"]

    assert ((gene_set=="2000mostvariable") and (len(whole_genes)==2000)) or ((gene_set=="full_genes") and (len(whole_genes)==26693))
    assert (len(samples)==47211) and (whole_clustering.shape[0]==len(samples))
    assert ((clusters_annot=="12clusters") and (whole_clustering.unique().shape[0]==12)) or ((clusters_annot=="new_clusters") and (whole_clustering.unique().shape[0]==16))

    print("... done")

    #get_genes = lambda x : sb.check_output("cat "+x+".csv | cut -d',' -f1", shell=True).decode("utf-8").split("\n")[1:-1]
    #get_samples = lambda x : sb.check_output("head -n1 "+x+".csv", shell=True).decode("utf-8").split("\n")[0].split(",")

    #print("Building targeted_matrix..."),

    #targeted_samples, targeted_genes = get_samples("raw_counts"), get_genes("raw_counts")
    #whole_fname = "raw_counts_variable_feature" if (os.path.exists("raw_counts_variable_feature.csv")) else "raw_counts_subset_feature"
    ## Restriction to panel genes
    #whole_samples, whole_genes = get_samples(whole_fname), get_genes(whole_fname)

    ###### GET GENES
    if (not os.path.exists(rfolder+"gene_list.pck")):
        if (pred_data_type=="facs"):
            missing = list(set(panel_genes).difference(whole_genes))
            #print(missing)
            whole_panel = list(set(panel_genes).intersection(whole_genes)) #+ list(filter(lambda x:x,missing_whole_genes))
            targeted_panel = None
        else:
            ## Ensure that all genes in targeted data are present in the whole dataset
            missing = list(set(panel_genes).difference(whole_genes))
            from Levenshtein import distance
            def get_closest_id(missing_gene, whole_genes):
                dists = [distance(missing_gene, gene) for gene in whole_genes]
                max_idx = np.argmin(dists)
                return [whole_genes[max_idx], dists[max_idx]]
            missing_ids = [[missing_gene]+get_closest_id(missing_gene, whole_genes) for missing_gene in missing]
    
            import mygene
            mg = mygene.MyGeneInfo()
            res = mg.querymany(missing+whole_genes, scopes='symbol', species='mouse', as_dataframe=True, returnall=True)
            miss = res["missing"]
            missing_in_index = list(set(missing).intersection(list(miss["query"])))
            if (len(missing_in_index) > 0):
                missing = list(set(missing).difference(missing_in_index))
            ## manually checked
            missing = [mis for mis in missing]# if (mis not in ['Jak3', 'IFNGR1', 'Map3k8', 'Nfkb1', 'Stat3'])]
            out = res["out"][["symbol", "_score", "entrezgene"]]
            missing_entrezgene = list(out.loc[missing]["entrezgene"])
            out = out.loc[[idx for idx in out.index if (idx not in missing)]]
            missing_whole_genes = [None]*len(missing)
            for msid, missing_id in enumerate(missing_entrezgene):
                match = out.loc[out["entrezgene"] == missing_id]
                if (len(match.index) > 0):
                    match = match["symbol"] if ("str" in str(type(match["symbol"]))) else list(match["symbol"])[0]
                else:
                    match = [ms_id[1] for ms_id in missing_ids if (ms_id[0]==missing[msid])][0]
                missing_whole_genes[msid] = match if ("str" in str(type(match))) else list(match)[0]
            targeted_panel = list(set(panel_genes).intersection(whole_genes)) + [mis for mis in missing if (str(missing_whole_genes[missing.index(mis)]) != "None")]
            whole_panel = list(set(panel_genes).intersection(whole_genes)) + list(filter(lambda x:x,missing_whole_genes))

        #missing_samples = list(set(targeted_samples).difference(whole_samples))
        #assert len(missing_samples) == 0
        #missing_samples = list(set(whole_samples).difference(targeted_samples))
        #assert len(missing_samples) == 0

        with open(rfolder+"gene_list.pck", "wb+") as f:
            pkl.dump({pred_data_type+"_panel": targeted_panel, "whole_panel": whole_panel}, f)

    else:
        with open(rfolder+"gene_list.pck", "rb") as f:
            di = pkl.load(f)
        targeted_panel, whole_panel = di[pred_data_type+"_panel"], di["whole_panel"]

    print("#genes = %d (target), %d (source)" % (len(targeted_panel) if (targeted_panel is not None) else 10, len(whole_panel)))

################################
## Counts for targeted RNAseq ##
################################

#targeted_matrix = pd.read_csv("raw_counts.csv").loc[targeted_panel][samples].dropna().T
#targeted_matrix = pd.read_csv("raw_counts.csv").loc[targeted_panel].dropna().T
#Nsamples_targeted, N_features = np.shape(targeted_matrix)

#print(set(targeted_genes).difference(set(targeted_panel)))

#print("... done")

    print("Building whole_matrix..."),

    #######################################################################################
    ## Counts for the targeted genes among the 2,000 most variable genes in whole RNAseq ##
    #######################################################################################

    if (not os.path.exists(rfolder+"raw_counts_subset_feature.csv")):
        sb.call("if [ -d 'ids.csv' ]; then rm -f ids.csv; fi", shell=True)
        with open("ids.csv", "w+") as f:
            f.write("\n".join(whole_panel))
        sb.call("head -n1 "+whole_fname+" > "+rfolder+"raw_counts_subset_feature.csv", shell=True)
        sb.call("grep -Fwf ids.csv "+whole_fname+" >> "+rfolder+"raw_counts_subset_feature.csv", shell=True)
        sb.call("rm -f ids.csv", shell=True)
    #whole_matrix = pd.read_csv("raw_counts_subset_feature.csv").loc[whole_panel][samples].dropna().T
    whole_matrix = pd.read_csv(rfolder+"raw_counts_subset_feature.csv", index_col=0).loc[whole_panel]
    in_samples = [x for x in list(whole_matrix.columns) if (x in samples)]
    whole_matrix = whole_matrix[in_samples].dropna().T
    whole_clustering = whole_clustering.loc[in_samples]
    #Nsamples_whole, N_features_ = np.shape(targeted_matrix)

    #assert Nsamples_whole == Nsamples_targeted
    #assert N_features == N_features_
    assert whole_clustering.shape[0]==whole_matrix.shape[0]

    print("... done")

    print("Processing input data..."),

    ####################################
    ## Normalizing data & Input data  ##
    ####################################

    X_fit, _ = data_normalisation(whole_matrix, normalisation_type, None)
    Y = np.ravel(list(whole_clustering))
    #Y_true = np.ravel(list(targeted_clustering))

    ##############################
    ## CROSS VALIDATION ROUTINE ##
    ##############################

    print("Running training..."),
    ## Number of neighbors
    N = int(0.5*np.sqrt(X_fit.shape[0]))
    params["n_neighbors"] = N

    def testing_routine(MODEL, YY, XX, verbose=False, extreme_verbose=False):
        seeds0 = [random.sample(range(int(1e8)),k=1)[0] for k in range(Niter_CV)]
        mean_auc_test, mean_auc_train = [None]*Niter_CV, [None]*Niter_CV
        mean_ari_test, mean_ari_train = [None]*Niter_CV, [None]*Niter_CV
        max_auc_test_di = {"fpr": None, "tpr": None}
        max_auc_test = -float("inf")
        for ns, seed0 in enumerate(seeds0):
            if (verbose):
                print("Iteration %d/%d" % (ns+1, len(seeds0)))
            random.seed(seed0)
            np.random.seed(seed0)
            cv = StratifiedShuffleSplit(n_splits=Nfold_CV, test_size=Ntest_size, random_state=seed0)
    
            results_cv = cross_validate(MODEL(0), X_fit, Y, cv=cv, scoring="roc_auc_"+(params["multiclass"] if (params["multiclass"]!="multinomial") else "ovr")+"_weighted", n_jobs=params["n_jobs"], return_train_score=True, return_estimator=True)
            mean_auc_test[ns] = np.mean(results_cv['test_score'])
            mean_auc_train[ns] = np.mean(results_cv['train_score'])
            ari_train, ari_test = [], []
            for i, (train_ids, test_ids) in enumerate(cv.split(X_fit, Y)):
                cls = results_cv["estimator"][i]
                ari_train.append(adjusted_rand_score(Y[train_ids], cls.predict(X_fit[train_ids,:])))
                ari_test.append(adjusted_rand_score(Y[test_ids], cls.predict(X_fit[test_ids,:])))
            mean_ari_test[ns] = np.mean(ari_train)
            mean_ari_train[ns] = np.mean(ari_test)
    
            max_auc_id = np.argmax(results_cv['test_score'])
            if (results_cv['test_score'][max_auc_id]>max_auc_test):
                max_estimator = results_cv['estimator'][max_auc_id]
                ids_test = [test_ids for (_, test_ids) in cv.split(X_fit, Y)][max_auc_id]
                Y_pred = max_estimator.predict(X_fit[test_ids,:])
                fpr, tpr, pre, rec = {}, {}, {}, {}
                for nclass in np.unique(Y).flatten().tolist():
                    class_true, class_pred = np.array(Y[test_ids]==nclass, dtype=int), np.array(Y_pred==nclass, dtype=int)
                    fpr_class, tpr_class, _ = roc_curve(class_true, class_pred)
                    pre_class, rec_class, _ = precision_recall_curve(class_true, class_pred)
                    fpr.setdefault(str(nclass), fpr_class)
                    tpr.setdefault(str(nclass), tpr_class)
                    pre.setdefault(str(nclass), pre_class)
                    rec.setdefault(str(nclass), rec_class)
                max_auc_test_di.update({"AUC (test)": results_cv["test_score"][max_auc_id], "ARI (test)": ari_test[max_auc_id], "pre": pre, "rec": rec, "fpr": fpr, "tpr": tpr, "estimator": max_estimator})
                max_auc_test = results_cv['test_score'][max_auc_id]
            if (verbose):
                print("mean AUC on test set: %.2f\tmax AUC on test set so far: %.2f" % (mean_auc_test[ns],max_auc_test))

        print("mean AUC Train: %.2f (+- %.2f)\tTest: %.2f (+- %.2f)" % (np.mean(mean_auc_train), np.std(mean_auc_train),np.mean(mean_auc_test), np.std(mean_auc_test)))
        print("mean ARI Train: %.2f (+- %.2f)\tTest: %.2f (+- %.2f)" % (np.mean(mean_ari_train), np.std(mean_ari_train),np.mean(mean_ari_test), np.std(mean_ari_test)))
        return max_auc_test_di, max_auc_test

#if (method == "centroids"):
#    centroids_fname = "centroids_metric="+params["metric"]+".csv"
#    x, y, s = X_fit_train, Y_train, samples_train
#    if (not os.path.exists(centroids_fname)):
#        all_clusters = list(set(y.tolist()))
#        centroids = pd.DataFrame([], index=[cluster for cluster in all_clusters], columns=whole_matrix.columns)
#        for sc, cluster in enumerate(all_clusters):
#            print("Computing centroid #"+str(sc+1)+" out of "+str(len(all_clusters)))
#            ids_cluster = [sidx for sidx,_ in enumerate(s) if (y[sidx] == cluster)]
#            barycenter = np.mean(x[ids_cluster,:], axis=0)
#            if (metric == "euclidean"):
#                centroids.loc[cluster] = barycenter
#            elif (metric == "cityblock"):
#                centroids.loc[cluster] = np.median(x[ids_cluster,:], axis=0)
#            else:
#                f = lambda xx : np.sum(pdist(xx.reshape(1,-1), x[ids_cluster, :], metric=params["dist"]))
#                res = minimize(f, barycenter, method='Nelder-Mead', tol=1e-6)
#                centroids.loc[cluster] = res.x
#        centroids.to_csv(centroids_fname)
#    else:
#        centroids = pd.read_csv(centroids_fname)
#        centroids.index = centroids[centroids.columns[0]]
#        centroids = centroids.drop(columns=centroids.columns[:1])
#    def predict(X):
#        ## rows = samples in X, columns = clusters in centroids
#        dist_mat = pdist(X, centroids.values, metric=params["dist"])
#        clusters = [centroids.index[idx] for idx in np.argmin(dist_mat, axis=1).tolist()]
#        return np.ravel(clusters)
#else:

#from sklearn.multiclass import OneVsRestClassifier
    if (True):
        model_fname = rfolder+method+"_model.pck"
        if (not os.path.exists(method+"_model.pck")):
            if (method == "NearestCentroid"):
               clf = lambda _ : neighbors.NearestCentroid(metric=params["metric"], shrink_threshold=params["shrink_threshold"])
            elif (method == "LogisticRegression"):
                clf = lambda _ : LogisticRegression(penalty="elasticnet", l1_ratio=params["l1_ratio"], tol=params["tol"], random_state=0, max_iter=params["max_iter"], multi_class=params["multiclass"], solver="saga",class_weight=params["class_weight"])
            elif (method == "RidgeClassifier"):
               clf = lambda _ : RidgeClassifier(alpha=params["alpha"], tol=params["tol"], random_state=0,class_weight="balanced")
            elif (method == "DecisionTree"):
                clf = lambda _ : DecisionTreeClassifier(criterion="gini", splitter="best", random_state=0)
            elif (method == "QuadraticDiscriminantAnalysis"):
                clf = lambda _ : QDA(tol=params["tol"])
            else:
                raise ValueError
            #clf.fit(X_fit,Y)
            max_auc_test_di, max_auc_test = testing_routine(clf, Y, X_fit, verbose=True)
            with open(model_fname, "wb+") as f:
                #pkl.dump(clf,f)
                pkl.dump(max_auc_test_di, f)
        else:
            with open(model_fname, "rb") as f:
                #clf = pkl.load(f)
                max_auc_test_di = pkl.load(f)

    clf = max_auc_test_di["estimator"]
    print(clf)

    print("... done!")

##############################
## COMPARE TO TARGETED DATA ##
##############################

#samples_targeted = [list(targeted_matrix.index).index(s) for s in samples]
#samples_whole = [list(whole_matrix.index).index(s) for s in samples]

#print("Compute scores..."),

## Compare clusters and clusterings
#from sklearn.metrics import adjusted_rand_score as ARI
#from sklearn.metrics import roc_auc_score as AUC
#from sklearn.preprocessing import OneHotEncoder
#from sklearn.model_selection import RepeatedStratifiedKFold
#onehot = OneHotEncoder(handle_unknown='error')
#onehot.fit(np.concatenate((Y,Y_true),axis=0).reshape(-1, 1))
#print(onehot.categories_)
## same clusters between (predicted) targeted data and (true) whole data?
#Y_whole = Y[samples_whole]
#Y_targeted_pred = clf.predict(X_pred_fit)[samples_targeted]
## same clusters between (predicted) targeted data and (predicted) whole data?
#Y_whole_pred = clf.predict(X_fit)[samples_whole]
## same clusters between (predicted) and (true) targeted data?
#Y_targeted = Y_true[samples_targeted]
## same clusters between (predicted) and (true) whole data?
#values_cols = ["Seurat(target)","Predicted(target)","Seurat(whole)","Predicted(whole)"]
#values = [Y_targeted,Y_targeted_pred,Y_whole,Y_whole_pred]
#values_onehot = [onehot.transform(y.reshape(-1, 1)).toarray() for y in values]
#ari_df, auc_df = [pd.DataFrame([], index=values_cols) for _ in range(2)]
#for iv, v in enumerate(values_cols):
#    print(v)
#    ari_df[v] = [ARI(values[iv], y) for y in values]
#    auc_list = []
#    for y in values_onehot:
#        keep_ids = np.argwhere(np.sum(values_onehot[iv],axis=0)).flatten().tolist()
#        true_rest,pred_rest = [x[:,keep_ids] for x in [values_onehot[iv], y]]
#        auc = AUC(true_rest,pred_rest,average='weighted',multi_class=params["multiclass"])
#        auc_list.append(auc)
#    auc_df[v] = auc_list

#print("ARI values")
#print(ari_df)
#print("AUC values")
#print(auc_df)

#print("... done!")

############################
## SAVE RESULTS           ##
############################

#print("Saving results..."),

#suffix = "_method=%s_metric=%s_%s" % (method,metric,sample_type)

## Save ARI/AUC results 
#ari_df.to_csv(rfolder+"aridf"+suffix+".csv")
#auc_df.to_csv(rfolder+"aucdf"+suffix+".csv")

## Save clusterings
#cols=["Seurat_clusters(whole)", "Learnt_clusters(whole)", "Seurat_clusters(targeted)", "Learnt_clusters(targeted)"]
#clusterings = pd.DataFrame([], index=samples, columns=cols)
#clusterings["Seurat_clusters(whole)"] = Y_whole
#clusterings["Learnt_clusters(whole)"] = Y_whole_pred
#clusterings["Seurat_clusters(targeted)"] = Y_targeted
#clusterings["Learnt_clusters(targeted)"] = Y_targeted_pred
#clusterings.to_csv(rfolder+"results"+suffix+".csv")

#print("... done!")

## Save matrices
#np.savetxt(rfolder+"X_fit"+suffix+".txt", X_fit)
#np.savetxt(rfolder+"X_pred_fit"+suffix+".txt", X_pred_fit)
