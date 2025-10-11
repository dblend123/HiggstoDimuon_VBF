#!/usr/bin/env python3
# python3 tmva_KS_perfold.py --era 12022 --alpha 1e-3 --dthr-pct 0.20
#
# 4-fold KS test between training-eval (signal) and validation (signal).
# Produces: per-fold CDF plots, a 2x2 grid, and a PASS/FAIL table.

import os, argparse
import numpy as np
import ROOT 

# IMPORTANT: import xgboost BEFORE ROOT to avoid symbol clashes
from xgboost import XGBClassifier

# Reuse your fold-aware loader and tree name from training
from tmva101_Training import load_folds, TREENAME  # :contentReference[oaicite:3]{index=3}

OUTDIR = "Outputs"

# ---------- small internal helpers ----------

def ks_2samp(x, y):
    """
    Two-sample Kolmogorov–Smirnov statistic (D) and asymptotic p-value.
    No SciPy required.
    """
    x = np.sort(np.asarray(x, dtype=np.float64))
    y = np.sort(np.asarray(y, dtype=np.float64))
    nx, ny = x.size, y.size

    # Walk both CDFs to compute sup |Fn - Gm|
    i = j = 0
    d = 0.0
    while i < nx and j < ny:
        if x[i] <= y[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i/nx - j/ny))

    # Smirnov asymptotic p-value (two-sided)
    en = np.sqrt(nx * ny / (nx + ny))
    lam = (en + 0.12 + 0.11/en) * d
    # p = 2 * sum_{k=1..∞} (-1)^{k-1} exp(-2 k^2 lam^2)
    p = 0.0
    for k in range(1, 101):
        p += 2.0 * ((-1.0) ** (k-1)) * np.exp(-2.0 * (k*k) * (lam*lam))
    p = float(np.clip(p, 0.0, 1.0))
    return d, p

def ecdf(a):
    a = np.sort(np.asarray(a, dtype=np.float64))
    n = a.size
    return a, np.arange(1, n+1)/n

def train_xgb(x, y, w):
    # Ensure numeric & finite
    x = np.asarray(x, dtype=np.float32)
    y = np.asarray(y, dtype=np.int32)
    w = np.asarray(w, dtype=np.float32)
    w[~np.isfinite(w)] = 0.0
    
    bdt = XGBClassifier(
        objective="binary:logistic",      # explicit
        eval_metric="logloss",            # explicit
        tree_method="hist",
        missing=np.nan, 
        max_bin=30, 
        max_depth=4,
        learning_rate=0.1,
        n_estimators=1000,
        booster="gbtree",
        n_jobs=30,
        nthread=None, 
        gamma=0, 
        min_child_weight=0.03, 
        max_delta_step=0, 
        subsample=0.5,
        colsample_bytree=1,
        colsample_bylevel=1, 
        reg_alpha=0,
        reg_lambda=1, 
        scale_pos_weight=0.7, 
        base_score=0.5, 
        random_state=0, 
        seed=None, 
        verbosity=1                      # show training progress
    )
    bdt.fit(x, y, sample_weight=w)
    return bdt

# ---------- plotting helpers ----------

def draw_cdf_pair(x_trainS, x_validS, title, out_png, out_pdf):
    # Compute ECDFs
    x1, F1 = ecdf(x_trainS)
    x2, F2 = ecdf(x_validS)

    c = ROOT.TCanvas("c", "", 700, 650)
    g1 = ROOT.TGraph(len(x1), x1.astype("f8"), F1.astype("f8"))
    g2 = ROOT.TGraph(len(x2), x2.astype("f8"), F2.astype("f8"))

    g1.SetTitle(f"{title};BDT score;CDF")
    g1.SetLineWidth(3); g1.SetLineColor(ROOT.kOrange+7)
    g2.SetLineWidth(3); g2.SetLineColor(ROOT.kAzure+2)
    g1.Draw("AL"); g2.Draw("L SAME")

    leg = ROOT.TLegend(0.15, 0.75, 0.45, 0.88)
    leg.AddEntry(g2, "Validation CDF", "l")
    leg.AddEntry(g1, "Train-eval CDF", "l")
    leg.Draw()

    c.SaveAs(out_pdf)
    c.SaveAs(out_png)

def draw_grid(figs, out_png, out_pdf):
    # figs: list of (x_trainS, x_validS, title)
    c = ROOT.TCanvas("grid", "", 1100, 900)
    c.Divide(2, 2)
    for i, (xtr, xva, title) in enumerate(figs, start=1):
        c.cd(i)
        x1, F1 = ecdf(xtr); x2, F2 = ecdf(xva)
        g1 = ROOT.TGraph(len(x1), x1.astype("f8"), F1.astype("f8"))
        g2 = ROOT.TGraph(len(x2), x2.astype("f8"), F2.astype("f8"))
        g1.SetTitle(f"{title};BDT score;CDF")
        g1.SetLineWidth(3); g1.SetLineColor(ROOT.kOrange+7)
        g2.SetLineWidth(3); g2.SetLineColor(ROOT.kAzure+2)
        g1.Draw("AL"); g2.Draw("L SAME")
        leg = ROOT.TLegend(0.15, 0.75, 0.45, 0.88)
        leg.AddEntry(g2, "Validation CDF", "l")
        leg.AddEntry(g1, "Train-eval CDF", "l")
        leg.Draw()
    c.SaveAs(out_pdf)
    c.SaveAs(out_png)

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Per-fold KS test: train-eval S vs validation S")
    ap.add_argument("--dir", default="DataPreparation_Files", help="Dir containing *_all_*.root")
    ap.add_argument("--era", required=True, help="Era like 12022")
    ap.add_argument("--alpha", type=float, default=1e-3, help="significance level for p-value test")
    ap.add_argument("--dthr-pct", type=float, default=0.20,
                    help="Optional absolute D(%) threshold for PASS (e.g. 0.20 = 0.2%)")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)

    sig_all = os.path.join(args.dir, f"{args.era}_all_signal.root")
    bkg_all = os.path.join(args.dir, f"{args.era}_all_background.root")
    if not (os.path.exists(sig_all) and os.path.exists(bkg_all)):
        raise SystemExit("Missing *_all_signal.root or *_all_background.root with 'fold' column.")

    # For each k in 0..3:
    #   Train on folds {(k+1)%4, (k+2)%4}, Validate on {k}, (Optionally Test on {(k+3)%4})
    rows = []
    grid_figs = []

    for k in range(4):
        folds_tr = ( (k+1) % 4, (k+2) % 4 )
        folds_va = ( k, )
        folds_te = ( (k+3) % 4, )  # not used here, but we show the intent

        # Load splits using your loader (annihilation, class-sum balancing built in):contentReference[oaicite:5]{index=5}
        (x_tr, y_tr, w_tr), (x_va, y_va, w_va), _ = load_folds(
            [sig_all], [bkg_all],
            folds_train=folds_tr, folds_valid=folds_va, folds_test=folds_te,
            sig_extra_weight_branch=None
        )

        # Mask signal only
        m_trS = (y_tr == 1)
        m_vaS = (y_va == 1)
        n_trS = int(m_trS.sum())
        n_vaS = int(m_vaS.sum())

        # Train & score
        bdt = train_xgb(x_tr, y_tr, w_tr)
        yhat_trS = bdt.predict_proba(np.asarray(x_tr[m_trS], dtype=np.float32))[:, 1]
        yhat_vaS = bdt.predict_proba(np.asarray(x_va[m_vaS], dtype=np.float32))[:, 1]
        
        std_tr = float(np.std(yhat_trS)); std_va = float(np.std(yhat_vaS))
        print(f"[debug] fold {k}: preds min/max trainS=({yhat_trS.min():.3f},{yhat_trS.max():.3f}) "
            f"validS=({yhat_vaS.min():.3f},{yhat_vaS.max():.3f}) stds=({std_tr:.4f},{std_va:.4f})")
        if std_tr == 0.0 or std_va == 0.0:
            print(f"[ERROR] fold {k}: constant predictions; check training inputs/labels/weights")
            continue

        # KS statistic & p-value
        D, p = ks_2samp(yhat_vaS, yhat_trS)
        D_pct = 100.0 * D

        # Pass/Fail criteria:
        #  (1) p-value >= alpha  AND/OR  (2) D <= threshold
        pass_p   = (p >= args.alpha)
        pass_d   = (D_pct <= args.dthr_pct)
        passed   = pass_p and pass_d

        rows.append((k, D_pct, n_vaS, args.alpha, args.dthr_pct, passed))

        png = os.path.join(OUTDIR, f"{args.era}_tmva_KS_fold{k}.png")
        pdf = os.path.join(OUTDIR, f"{args.era}_tmva_KS_fold{k}.pdf")
        draw_cdf_pair(yhat_trS, yhat_vaS,
                f"BDT distribution of signal sample {args.era}_(fold {k})",png, pdf)
       
        grid_figs.append((yhat_trS, yhat_vaS, f"{args.era}_Fold {k}"))

        # Debug Checks 
        print(f"[debug] fold {k}: train S preds mean={yhat_trS.mean():.3f}, "
             f"valid S preds mean={yhat_vaS.mean():.3f}, "
             f"len(trainS)={len(yhat_trS)}, len(validS)={len(yhat_vaS)}")
        print(f"[debug] fold {k}: preds min={yhat_vaS.min():.3f}, max={yhat_vaS.max():.3f}, "
             f"unique≈{np.unique(yhat_vaS.round(3)).size}")

    # 2x2 grid
    draw_grid(grid_figs,
              os.path.join(OUTDIR, f"{args.era}_tmva_KS_grid.png"),
              os.path.join(OUTDIR, f"{args.era}_tmva_KS_grid.pdf"))

    # Print table like your example
    print("\nFold | KS statistic | N events | alpha | threshold | PASS")
    for (k, D_pct, n, alpha, thr, passed) in rows:
        print(f"{k:>4} | {D_pct:10.2f}% | {n:8d} | {alpha:.3g} | {thr:.2f}%   | {'TRUE' if passed else 'FALSE'}")

if __name__ == "__main__":
    main()

