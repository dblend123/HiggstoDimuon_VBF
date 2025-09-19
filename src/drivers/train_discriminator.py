# train_discriminator.py
# Usage: python train_discriminator.py train_signal.root train_background.root
# Requires: pip install xgboost uproot numpy; and ROOT >= 6.24 (for TMVA SaveXGBoost)
# Tree name is assumed to be "Events".

import sys
import numpy as np
import uproot
import ROOT
from xgboost import XGBClassifier

# --- Edit branch names here to match your ntuples ---
BRANCHES = dict(
    y_mumu      = "y_mumu",   
    y_j1        = "y_j1",     
    y_j2        = "y_j2",     
    pt_mumu     = "pt_mumu", phi_mumu="phi_mumu",
    pt_j1  ="pt_j1",   phi_j1  ="phi_j1",
    pt_j2  ="pt_j2",   phi_j2  ="phi_j2",
    # Optional per-event weight branch; set to None to use 1.0
    weight ="eventWeight"
)

def load_arrays(filename):
    with uproot.open(filename) as f:
        tree = f["Events"]
        arr = {k: tree[BRANCHES[k]].array(library="np") for k in BRANCHES if BRANCHES[k] is not None}
    # helpers: px, py from pt, phi
    def px(pt, phi): return pt * np.cos(phi)
    def py(pt, phi): return pt * np.sin(phi)

    # Zeppenfeld variable: z = |y_mumu - 0.5*(y_j1 + y_j2)|
    z = np.abs(arr["y_mumu"] - 0.5*(arr["y_j1"] + arr["y_j2"]))

    # pT-balance ratio: |vector sum| / (scalar sum)
    px_sum = px(arr["pt_mumu"], arr["phi_mumu"]) + px(arr["pt_j1"], arr["phi_j1"]) + px(arr["pt_j2"], arr["phi_j2"])
    py_sum = py(arr["pt_mumu"], arr["phi_mumu"]) + py(arr["pt_j1"], arr["phi_j1"]) + py(arr["pt_j2"], arr["phi_j2"])
    vec_sum = np.sqrt(px_sum**2 + py_sum**2)
    scalar_sum = arr["pt_mumu"] + arr["pt_j1"] + arr["pt_j2"]
    pt_balance = np.divide(vec_sum, scalar_sum, out=np.zeros_like(vec_sum), where=scalar_sum>0)

    X = np.vstack([z, pt_balance]).T
    w = arr["eventWeight"] if "eventWeight" in arr else np.ones_like(z)
    return X, w

def main(sig_file, bkg_file, out_model="tmva_discriminator.root"):
    Xs, ws = load_arrays(sig_file)
    Xb, wb = load_arrays(bkg_file)

    X = np.vstack([Xs, Xb])
    y = np.hstack([np.ones(len(Xs)), np.zeros(len(Xb))])

    # balance classes via sample weights (optional but helpful)
    w_sig_scale = (len(X) / len(Xs)) if len(Xs) > 0 else 1.0
    w_bkg_scale = (len(X) / len(Xb)) if len(Xb) > 0 else 1.0
    w = np.hstack([ws * w_sig_scale, wb * w_bkg_scale])

    # Simple, strong baseline BDT
    clf = XGBClassifier(
        n_estimators=600,
        max_depth=3,
        subsample=0.7,
        colsample_bytree=0.9,
        learning_rate=0.05,
        reg_lambda=1.0,
        n_jobs=4,
        random_state=42,
        tree_method="hist",
    )
    clf.fit(X, y, sample_weight=w)

    # Save model to TMVA format for fast inference inside ROOT
    # Note: variable order must match training order: ["z", "pt_balance"]
    ROOT.TMVA.Experimental.SaveXGBoost(clf, "VBF_Discr", out_model, num_inputs=X.shape[1])
    print(f"Training done on {len(X)} events. Wrote TMVA model to {out_model}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python train_discriminator.py signal.root background.root [out.root]")
        sys.exit(1)
    sig, bkg = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else "tmva_discriminator.root"
    main(sig, bkg, out)

