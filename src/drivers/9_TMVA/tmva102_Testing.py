#!/usr/bin/env python3
# python3 tmva102_Testing.py --model 12022_tmva101.root --era 12022

import os
import argparse
import numpy as np
import ROOT
from sklearn.metrics import roc_curve, auc

# Reuse the same tree name and loader you trained with
from tmva101_Training import load_folds, TREENAME  # no load_data here

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="DataPreparation_Files", help="Directory with prepared ROOT files")
    ap.add_argument("--era", required=True, help="Era, e.g. 12022")
    ap.add_argument("--model", default="tmva101.root", help="TMVA model file")
    ap.add_argument("--name",  default="myBDT", help="Model name inside TMVA file")
    args = ap.parse_args()

    DIR = args.dir
    ERA = args.era
    OUTDIR = "Outputs"
    
    sig_all = os.path.join(DIR, f"{ERA}_all_signal.root")
    bkg_all = os.path.join(DIR, f"{ERA}_all_background.root")
    sig_test = os.path.join(DIR, f"{ERA}_test_signal.root")
    bkg_test = os.path.join(DIR, f"{ERA}_test_background.root")

    # Prefer the unified _all_ files (contain 'fold'); fallback to _test_ files
    sig_files = [sig_all] if os.path.exists(sig_all) else [sig_test]
    bkg_files = [bkg_all] if os.path.exists(bkg_all) else [bkg_test]

    print("[tmva102] using files:")
    print("  signal:", sig_files)
    print("  bkg   :", bkg_files)

    # We only need the test split. No extra signal reweighting at test time.
    (_,_,_), (_,_,_), (x_te, y_te, w_te) = load_folds(
        sig_files, bkg_files,
        folds_train=(), folds_valid=(), folds_test=(3,),
        sig_extra_weight_branch=None
    )

    if x_te.shape[0] == 0:
        print("[tmva102][WARN] test split is empty. Does your input have a non-constant 'fold' column?")
        return

    # --- Load TMVA and the runtime model ---
    # Load libs explicitly to avoid unresolved symbols on some builds
    ROOT.gSystem.Load("libRIO")
    ROOT.gSystem.Load("libTree")
    ROOT.gSystem.Load("libTreePlayer")
    ROOT.gSystem.Load("libTMVA")

    model_path = args.model
    if ROOT.gSystem.AccessPathName(model_path, ROOT.kReadPermission):
        ROOT.Error("tmva102_Testing.py", f"{model_path} does not exist or is not readable")
        raise SystemExit(1)

    MODEL_NAME = args.name
    bdt = ROOT.TMVA.Experimental.RBDT[""](MODEL_NAME, model_path)

    # --- Predict and evaluate ---
    x = np.asarray(x_te, dtype=np.float32)
    y_pred = bdt.Compute(x)

    print("[tmva102] y_pred mean/std/min/max:",
          float(np.mean(y_pred)), float(np.std(y_pred)),
          float(np.min(y_pred)), float(np.max(y_pred)))

    fpr, tpr, _ = roc_curve(y_te, y_pred, sample_weight=w_te)
    auc_score = auc(fpr, tpr)
    print(f"[tmva102] Weighted AUC: {auc_score:.5f}")

    # Operating points
    def interp_tpr_at_fpr(target_fpr):
        if target_fpr <= fpr[0]: return tpr[0]
        if target_fpr >= fpr[-1]: return tpr[-1]
        return np.interp(target_fpr, fpr, tpr)

    tpr_at_1pct = interp_tpr_at_fpr(0.01)
    youden_idx = int(np.argmax(tpr - fpr))
    print(f"[tmva102] TPR @ 1% FPR: {tpr_at_1pct:.3f}")
    print(f"[tmva102] Best threshold (Youden J): idx={youden_idx}, "
          f"TPR={tpr[youden_idx]:.3f}, FPR={fpr[youden_idx]:.3f}")

    # --- Plot ROC ---
    OUTDIR = "Outputs"
    os.makedirs(OUTDIR, exist_ok=True)
    
    c = ROOT.TCanvas("roc", "", 700, 650)
    gr = ROOT.TGraph(len(fpr), np.asarray(fpr, dtype="f8"), np.asarray(tpr, dtype="f8"))
    gr.SetTitle("ROC;False-positive rate;True-positive rate")
    gr.SetLineWidth(3)
    gr.SetLineColor(ROOT.kRed)
    gr.Draw("AL")
    gr.GetXaxis().SetRangeUser(0.0, 1.0)
    gr.GetYaxis().SetRangeUser(0.0, 1.0)

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.15, 0.92, f"AUC = {auc_score:.4f}")

    c.Update()
    c.SaveAs(os.path.join(OUTDIR, f"{ERA}_tmva102_roc.pdf"))
    c.SaveAs(os.path.join(OUTDIR, f"{ERA}_tmva102_roc.png"))

    with open(os.path.join(OUTDIR, "tmva102_metrics.txt"), "w") as f:
        f.write(f"AUC {auc_score:.6f}\n")
        f.write(f"TPR_at_FPR1pct {tpr_at_1pct:.6f}\n")
        f.write(f"Youden_index {youden_idx}\n")
        f.write(f"TPR_at_Youden {tpr[youden_idx]:.6f}\n")
        f.write(f"FPR_at_Youden {fpr[youden_idx]:.6f}\n")
    with open(os.path.join(OUTDIR, "tmva102_roc_points.csv"), "w") as f:
        f.write("fpr,tpr\n")
        for a,b in zip(fpr, tpr):
            f.write(f"{a:.8f},{b:.8f}\n")

    print(f"[tmva102] Wrote Outputs/{ERA}_tmva102_roc.(pdf|png), {ERA}_tmva102_metrics.txt, {ERA}_tmva102_roc_points.csv")

if __name__ == "__main__":
    main()

