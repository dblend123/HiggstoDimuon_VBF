# XGBoost must be imported before ROOT
# source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
# python3 tmva101_Training.py

import os
import ROOT
import numpy as np
from xgboost import XGBClassifier

from tmva100_DataPreparation import VARIABLES

print(ROOT.gROOT.GetVersion())
print(ROOT.gSystem.Getenv("ROOTSYS"))

TREENAME = "mva"
WEIGHT_BRANCH = "weight"   # branch in your mva tree

def _sanitize_weights(w):
    w = np.asarray(w, dtype=np.float64)
    # replace NaN/±Inf with 0
    w[~np.isfinite(w)] = 0.0
    # use absolute value (NLO negatives -> positive for training)
    w = np.abs(w)
    # clip to strictly positive (xgboost disallows 0)
    eps = 1e-6
    w = np.where(w > eps, w, eps)
    return w

def load_data(signal_filename, background_filename):
    ds = ROOT.RDataFrame(TREENAME, signal_filename).AsNumpy()
    db = ROOT.RDataFrame(TREENAME, background_filename).AsNumpy()

    x_sig = np.vstack([ds[v] for v in VARIABLES]).T
    x_bkg = np.vstack([db[v] for v in VARIABLES]).T
    x = np.vstack([x_sig, x_bkg])

    ns, nb = x_sig.shape[0], x_bkg.shape[0]
    y = np.hstack([np.ones(ns), np.zeros(nb)])

    # event weights (sanitized) × class balancing
    if WEIGHT_BRANCH in ds and WEIGHT_BRANCH in db:
        ws_raw, wb_raw = ds[WEIGHT_BRANCH], db[WEIGHT_BRANCH]
    else:
        ws_raw, wb_raw = np.ones(ns), np.ones(nb)

    ws = _sanitize_weights(ws_raw)
    wb = _sanitize_weights(wb_raw)

    scale_s = (ns + nb) / max(ns, 1)
    scale_b = (ns + nb) / max(nb, 1)
    w = np.hstack([ws * scale_s, wb * scale_b])

    # (optional) quick diagnostics
    neg_s = np.sum(np.asarray(ws_raw) < 0)
    neg_b = np.sum(np.asarray(wb_raw) < 0)
    print(f"[info] weights: ns={ns} (neg={neg_s}), nb={nb} (neg={neg_b})")

    return x, y, w

if __name__ == "__main__":
    # list of (signal, background) files for all eras
    DIR = "DataPreparation_Files"
    eras = [
        (os.path.join(DIR, "12022_train_signal.root"), os.path.join(DIR, "12022_train_background.root")),
        (os.path.join(DIR, "22022_train_signal.root"), os.path.join(DIR, "22022_train_background.root")),
        (os.path.join(DIR, "12023_train_signal.root"), os.path.join(DIR, "12023_train_background.root")),
        (os.path.join(DIR, "22023_train_signal.root"), os.path.join(DIR, "22023_train_background.root")),
    ]

    xs, ys, ws = [], [], []
    for sig, bkg in eras:
        x, y, w = load_data(sig, bkg)
        xs.append(x)
        ys.append(y)
        ws.append(w)

    # concatenate across eras
    x_all = np.vstack(xs)
    y_all = np.hstack(ys)
    w_all = np.hstack(ws)

    bdt = XGBClassifier(
        n_estimators=500,
        max_depth=3,
        learning_rate=0.05,
        subsample=0.5,
        colsample_bytree=0.9,
        tree_method="hist",
        random_state=42,
        n_jobs=4,
    )

    bdt.fit(x_all, y_all, sample_weight=w_all)

    ROOT.gSystem.Load("libTMVA")
    print("Training done on", x_all.shape[0], "events. Saving model in tmva101.root")
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, "myBDT", "tmva101.root", num_inputs=x_all.shape[1])

