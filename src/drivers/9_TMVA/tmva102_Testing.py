#!/usr/bin/env python3
"""
tmva102_Testing.py

Evaluate a trained TMVA RBDT on the held-out fold (fold==3) of the prepared
signal/background snapshots. Writes weighted ROC, AUC, TPR@1%FPR and Youden J.

Example:
    python3 tmva102_Testing.py --era 12022 \
        --dir DataPreparation_Files/12022 \
        --model Outputs/Training/12022/12022_tmva101_background.root
"""

import os, argparse
import numpy as np
import ROOT
from sklearn.metrics import roc_curve, auc

from tmva101_Training import load_folds, TREENAME  # reuse loader / weighting


def _prefer_all_or_test(directory, era, kind):
    allf  = os.path.join(directory, f"{era}_all_{kind}.root")
    testf = os.path.join(directory, f"{era}_test_{kind}.root")
    return [allf] if os.path.exists(allf) else [testf]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir",   default="DataPreparation_Files")
    ap.add_argument("--era",   required=True)
    ap.add_argument("--model", default="tmva101.root")
    ap.add_argument("--name",  default="myBDT")
    args = ap.parse_args()

    sig_files = _prefer_all_or_test(args.dir, args.era, "signal")
    bkg_files = _prefer_all_or_test(args.dir, args.era, "background")
    print("[tmva102] signal:", sig_files, "\n[tmva102] bkg   :", bkg_files)

    # Only need fold==3 (test). The unused training/validation tuples are dropped.
    _, _, (x_te, y_te, w_te) = load_folds(
        sig_files, bkg_files,
        folds_train=(), folds_valid=(), folds_test=(3,),
    )
    if x_te.shape[0] == 0:
        print("[tmva102][WARN] empty test split — does the input contain a 'fold' column?")
        return

    for lib in ("libRIO", "libTree", "libTreePlayer", "libTMVA"):
        ROOT.gSystem.Load(lib)
    if ROOT.gSystem.AccessPathName(args.model, ROOT.kReadPermission):
        raise SystemExit(f"[tmva102] model not readable: {args.model}")

    bdt = ROOT.TMVA.Experimental.RBDT[""](args.name, args.model)
    y_pred = bdt.Compute(np.asarray(x_te, dtype=np.float32))
    print(f"[tmva102] score range: mean={float(np.mean(y_pred)):.3f} "
          f"std={float(np.std(y_pred)):.3f} "
          f"min={float(np.min(y_pred)):.3f} max={float(np.max(y_pred)):.3f}")

    fpr, tpr, _ = roc_curve(y_te, y_pred, sample_weight=w_te)
    auc_score   = auc(fpr, tpr)
    tpr_at_1pct = np.interp(0.01, fpr, tpr)
    youden_idx  = int(np.argmax(tpr - fpr))
    print(f"[tmva102] Weighted AUC: {auc_score:.5f}")
    print(f"[tmva102] TPR @ 1% FPR: {tpr_at_1pct:.3f}")
    print(f"[tmva102] Youden J: idx={youden_idx} TPR={tpr[youden_idx]:.3f} FPR={fpr[youden_idx]:.3f}")

    OUT = "Outputs"
    os.makedirs(OUT, exist_ok=True)

    c  = ROOT.TCanvas("roc", "", 700, 650)
    gr = ROOT.TGraph(len(fpr), np.asarray(fpr, dtype="f8"), np.asarray(tpr, dtype="f8"))
    gr.SetTitle("ROC;False-positive rate;True-positive rate")
    gr.SetLineWidth(3); gr.SetLineColor(ROOT.kRed)
    gr.Draw("AL")
    gr.GetXaxis().SetRangeUser(0.0, 1.0)
    gr.GetYaxis().SetRangeUser(0.0, 1.0)
    lat = ROOT.TLatex(); lat.SetNDC(True); lat.SetTextSize(0.04)
    lat.DrawLatex(0.15, 0.92, f"AUC = {auc_score:.4f}")
    c.Update()
    for ext in ("pdf", "png"):
        c.SaveAs(os.path.join(OUT, f"{args.era}_tmva102_roc.{ext}"))

    with open(os.path.join(OUT, "tmva102_metrics.txt"), "w") as f:
        f.write(f"AUC {auc_score:.6f}\n")
        f.write(f"TPR_at_FPR1pct {tpr_at_1pct:.6f}\n")
        f.write(f"Youden_index {youden_idx}\n")
        f.write(f"TPR_at_Youden {tpr[youden_idx]:.6f}\n")
        f.write(f"FPR_at_Youden {fpr[youden_idx]:.6f}\n")
    with open(os.path.join(OUT, "tmva102_roc_points.csv"), "w") as f:
        f.write("fpr,tpr\n")
        for a, b in zip(fpr, tpr):
            f.write(f"{a:.8f},{b:.8f}\n")
    print(f"[tmva102] Wrote ROC + metrics for era {args.era} into {OUT}/")


if __name__ == "__main__":
    main()
