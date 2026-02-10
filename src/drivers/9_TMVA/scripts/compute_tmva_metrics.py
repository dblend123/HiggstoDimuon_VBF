#!/usr/bin/env python3
import ROOT, numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import ks_2samp

infile = "Outputs/2024/scored_train.root"   # change if needed
treename = "mva"
vars_list = ["bdt_score","fold","weight","isSignal"]  # replace isSignal with your truth branch name

rdf = ROOT.RDataFrame(treename, infile)
arrs = rdf.AsNumpy(vars_list)
scores = np.asarray(arrs["bdt_score"])
folds  = np.asarray(arrs["fold"], dtype=int)
weights= np.asarray(arrs["weight"], dtype=float)
y      = np.asarray(arrs["isSignal"], dtype=int)

def weighted_auc(y_true, y_score, w):
    try:
        return roc_auc_score(y_true, y_score, sample_weight=w)
    except Exception as e:
        return None

print("Overall AUC:", weighted_auc(y, scores, weights))

for f in sorted(set(folds)):
    idx = (folds == f)
    if idx.sum() > 0:
        auc = weighted_auc(y[idx], scores[idx], weights[idx])
        print(f"Fold {f}: entries={idx.sum():d}, AUC={auc:.5f}")

# KS: compare train (folds 0,1) vs val (fold 2)
train_idx = (folds==0)|(folds==1)
val_idx   = (folds==2)
if val_idx.sum()>0:
    stat,p = ks_2samp(scores[train_idx], scores[val_idx])
    print("KS train vs val: stat=%.4f p=%.4g" % (stat,p))
else:
    print("No fold 2 entries to run KS.")

# KS train vs test (fold 3)
test_idx = (folds==3)
if test_idx.sum()>0:
    stat,p = ks_2samp(scores[train_idx], scores[test_idx])
    print("KS train vs test: stat=%.4f p=%.4g" % (stat,p))
else:
    print("No fold 3 entries to run KS.")


