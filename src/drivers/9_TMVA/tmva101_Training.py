# XGBoost must be imported before ROOT
# source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
from xgboost import XGBClassifier

import ROOT
import numpy as np

print(ROOT.gROOT.GetVersion())
print(ROOT.gSystem.Getenv("ROOTSYS"))


variables = [
    "z","pt_balance","mu1_pt_over_m","mu2_pt_over_m",
    "mu1_eta","mu2_eta","jet1_eta","m_mumu","deta_jj","dphi_jj",
    "y_mumu","pt_j1","pt_j2",
]

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

    x_sig = np.vstack([ds[v] for v in variables]).T
    x_bkg = np.vstack([db[v] for v in variables]).T
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
    x, y, w = load_data("train_signal.root", "train_background.root")

    bdt = XGBClassifier(
        n_estimators=500,
        max_depth=3,
        learning_rate=0.05,
        subsample=0.7,
        colsample_bytree=0.9,
        tree_method="hist",
        random_state=42,
        n_jobs=4,
    )

    bdt.fit(x, y, sample_weight=w)

    # Ensure TMVA is loaded; no header includes needed
    ROOT.gSystem.Load("libTMVA")

    print("Training done on", x.shape[0], "events. Saving model in tmva101.root")
    print("SaveXGBoost available? ", hasattr(ROOT.TMVA.Experimental, "SaveXGBoost"))
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, "myBDT", "tmva101.root", num_inputs=x.shape[1])

