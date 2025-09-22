#!/usr/bin/env python3

import numpy as np
import ROOT

# Reuse the same variable order + loader you trained with
from tmva101_Training import load_data, TREENAME, WEIGHT_BRANCH

# --- 1) Load test data ---
x, y_true, w = load_data("test_signal.root", "test_background.root")

# xgboost isn't required here (we use TMVA reader), but TMVA must be loaded
ROOT.gSystem.Load("libTMVA")

# --- 2) Load the trained model (TMVA format) ---
MODEL_FILE = "tmva101.root"
MODEL_NAME = "myBDT"

if ROOT.gSystem.AccessPathName(MODEL_FILE, ROOT.kReadPermission):
    ROOT.Error("tmva102_Testing.py", f"{MODEL_FILE} does not exist or is not readable")
    raise SystemExit(1)

# RBDT is the TMVA runtime model
bdt = ROOT.TMVA.Experimental.RBDT[""](MODEL_NAME, MODEL_FILE)

# --- 3) Predict scores on the test set ---
# Ensure float32 for fastest path
x = np.asarray(x, dtype=np.float32)
y_pred = bdt.Compute(x)  # returns numpy array of scores

# --- 4) Weighted ROC / AUC (scikit-learn) ---
from sklearn.metrics import roc_curve, auc

fpr, tpr, _ = roc_curve(y_true, y_pred, sample_weight=w)
auc_score = auc(fpr, tpr)
print(f"[tmva102] Weighted AUC: {auc_score:.5f}")

# Optional: a couple of operating points
# (TPR at 1% FPR and Youden’s J)
def interp_tpr_at_fpr(target_fpr):
    if target_fpr <= fpr[0]: return tpr[0]
    if target_fpr >= fpr[-1]: return tpr[-1]
    return np.interp(target_fpr, fpr, tpr)

tpr_at_1pct = interp_tpr_at_fpr(0.01)
youden_idx = np.argmax(tpr - fpr)
print(f"[tmva102] TPR @ 1% FPR: {tpr_at_1pct:.3f}")
print(f"[tmva102] Best threshold (Youden J): idx={youden_idx}, TPR={tpr[youden_idx]:.3f}, FPR={fpr[youden_idx]:.3f}")

# --- 5) Plot and save ROC curve ---
c = ROOT.TCanvas("roc", "", 700, 650)
gr = ROOT.TGraph(len(fpr), np.array(fpr, dtype="f8"), np.array(tpr, dtype="f8"))
gr.SetTitle(f"ROC;False-positive rate;True-positive rate")
gr.SetLineWidth(3)
gr.SetLineColor(ROOT.kRed)
gr.Draw("AL")
gr.GetXaxis().SetRangeUser(0.0, 1.0)
gr.GetYaxis().SetRangeUser(0.0, 1.0)

# AUC label
latex = ROOT.TLatex()
latex.SetNDC(True)
latex.SetTextSize(0.04)
latex.DrawLatex(0.15, 0.92, f"AUC = {auc_score:.4f}")

c.Update()
c.SaveAs("tmva102_roc.pdf")
c.SaveAs("tmva102_roc.png")
print("[tmva102] Wrote tmva102_roc.pdf and tmva102_roc.png")

