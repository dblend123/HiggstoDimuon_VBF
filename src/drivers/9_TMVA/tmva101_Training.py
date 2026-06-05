#!/usr/bin/env python3
"""
tmva101_Training.py

Train one XGBoost BDT per (era, background-category) against the VBF H->mumu
signal, using the prepared snapshots from tmva100_DataPreparation.

Folds: train=(0,1), valid=(2), test=(3) → ~50/25/25.
Weights: |w| with negative-MC-weight annihilation, then sqrt(SwS * SwB)
class re-balancing so logloss is unbiased on weighted MC.
Model is saved both as XGBoost JSON and as a TMVA RBDT ROOT file.

Example:
    python3 tmva101_Training.py --era 12022 \
        --base-dir DataPreparation_Files/12022 --bkg-cats background
"""

# XGBoost must be imported before ROOT to avoid symbol clashes.
from xgboost import XGBClassifier
import argparse, os, numpy as np, ROOT

from tmva100_DataPreparation import VARIABLES

TREENAME = "mva"
WEIGHT_BRANCH = "weight"


# ----------------------------------------------------------------------
# Weight helpers
# ----------------------------------------------------------------------
def _sanitize(w):
    w = np.asarray(w, dtype=np.float64)
    w[~np.isfinite(w)] = 0.0
    return np.clip(np.abs(w), 1e-9, None)


def _annihilate(w_signed):
    """Cancel equal sums of +/- MC weights; keep the surplus on the majority sign."""
    w = np.asarray(w_signed, dtype=np.float64)
    pos = np.clip(w, 0.0, None)
    neg = -np.clip(w, None, 0.0)
    sp, sn = pos.sum(), neg.sum()
    if sp <= 0 or sn <= 0:
        return np.clip(w, 0.0, None)
    return pos * (1.0 - min(sp, sn) / sp)


def _balance(wS, wB):
    sS, sB = wS.sum(), wB.sum()
    if sS <= 0 or sB <= 0:
        return wS, wB
    g = np.sqrt(sS * sB)
    return wS * (g / sS), wB * (g / sB)


# ----------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------
def _load_numpy(files, vars_list):
    vs = ROOT.std.vector("string")()
    for f in files:
        vs.push_back(f)
    rdf = ROOT.RDataFrame(TREENAME, vs)
    cols = set(rdf.GetColumnNames())
    fetch = vars_list + [c for c in (WEIGHT_BRANCH, "fold") if c in cols]
    arr = rdf.AsNumpy(fetch)
    x = np.vstack([arr[v] for v in vars_list]).T
    w = np.asarray(arr[WEIGHT_BRANCH], dtype=np.float64) if WEIGHT_BRANCH in arr else np.ones(x.shape[0])
    f = np.asarray(arr["fold"], dtype=np.int32) if "fold" in arr else np.zeros(x.shape[0], np.int32)
    return x, w, f


def _prep(xS, wS, xB, wB):
    wS = _sanitize(_annihilate(_sanitize(wS) * np.sign(wS)))
    wB = _sanitize(_annihilate(_sanitize(wB) * np.sign(wB)))
    wS, wB = _balance(wS, wB)
    y = np.hstack([np.ones(xS.shape[0]), np.zeros(xB.shape[0])])
    x = np.vstack([xS, xB])
    w = np.hstack([wS, wB])
    return x, y, w


def load_folds(sig_files, bkg_files,
               folds_train=(0, 1), folds_valid=(2,), folds_test=(3,)):
    xS, wS, fS = _load_numpy(sig_files, VARIABLES)
    xB, wB, fB = _load_numpy(bkg_files, VARIABLES)

    def slc(x, w, f, allowed):
        m = np.isin(f, np.asarray(allowed, dtype=np.int32))
        return x[m], w[m]

    tr = _prep(*slc(xS, wS, fS, folds_train),  *slc(xB, wB, fB, folds_train))
    va = _prep(*slc(xS, wS, fS, folds_valid),  *slc(xB, wB, fB, folds_valid))
    te = _prep(*slc(xS, wS, fS, folds_test),   *slc(xB, wB, fB, folds_test))

    print(f"[info] train: S={int((tr[1]==1).sum())} B={int((tr[1]==0).sum())}  "
          f"valid: S={int((va[1]==1).sum())} B={int((va[1]==0).sum())}  "
          f"test:  S={int((te[1]==1).sum())} B={int((te[1]==0).sum())}")
    return tr, va, te


def find_prep_files(base_dir, era, kind):
    allf = os.path.join(base_dir, f"{era}_all_{kind}.root")
    if os.path.exists(allf):
        return [allf]
    return [os.path.join(base_dir, f"{era}_train_{kind}.root"),
            os.path.join(base_dir, f"{era}_test_{kind}.root")]


# ----------------------------------------------------------------------
# Training
# ----------------------------------------------------------------------
def _xgb_classifier():
    return XGBClassifier(
        tree_method="hist", max_bin=30,
        max_depth=4, learning_rate=0.1, n_estimators=1000,
        eval_metric="logloss", booster="gbtree", n_jobs=30,
        gamma=0, min_child_weight=0.03,
        subsample=0.5, colsample_bytree=1, colsample_bylevel=1,
        reg_alpha=0, reg_lambda=1, scale_pos_weight=0.7,
        base_score=0.5, random_state=0,
    )


def _save_feature_importance(bdt, txt_path, cat):
    raw = bdt.get_booster().get_score(importance_type="gain")
    gains = {VARIABLES[i]: raw.get(f"f{i}", 0.0) for i in range(len(VARIABLES))}
    tot = sum(gains.values()) or 1.0
    pct = sorted(((k, 100.0 * v / tot) for k, v in gains.items()), key=lambda x: -x[1])
    print(f"\n[info][{cat}] Feature importance (gain, % of total):")
    for k, p in pct:
        print(f"  {k:20s} : {p:6.2f}%")
    with open(txt_path, "w") as f:
        f.write(f"Feature importance (gain, normalized to 100%) — category {cat}\n")
        for k, p in pct:
            f.write(f"{k:20s} : {p:6.2f}%\n")
    print(f"[info] Wrote {txt_path}")


def _save_tmva(bdt, out_path, n_inputs):
    ROOT.gSystem.Load("libTMVA")
    print(f"[info] Saving TMVA model to {out_path}")
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, "myBDT", out_path, num_inputs=n_inputs)
    sz = os.path.getsize(out_path) if os.path.exists(out_path) else 0
    if sz < 1000:
        if os.path.exists(out_path):
            os.remove(out_path)
        raise RuntimeError(f"SaveXGBoost produced an incomplete file ({sz} bytes)")
    print(f"[info] Model saved ({sz} bytes)")


def train_one(era, base_dir, bkg_cat, outdir_root="Outputs/Training"):
    sig_files = find_prep_files(base_dir, era, "signal")
    bkg_files = find_prep_files(base_dir, era, bkg_cat)
    for p in sig_files + bkg_files:
        if not os.path.exists(p):
            raise FileNotFoundError(p)
    print(f"[debug] category={bkg_cat}\n        sig={sig_files}\n        bkg={bkg_files}")

    (x_tr, y_tr, w_tr), (x_va, y_va, w_va), _ = load_folds(sig_files, bkg_files)
    if x_va.shape[0] == 0:
        print("[WARN] empty validation set — check that 'fold' is present in the inputs")

    bdt = _xgb_classifier()
    bdt.fit(x_tr, y_tr, sample_weight=w_tr,
            eval_set=[(x_va, y_va)], sample_weight_eval_set=[w_va],
            verbose=False)

    outdir = os.path.join(outdir_root, str(era))
    os.makedirs(outdir, exist_ok=True)
    bdt.save_model(os.path.join(outdir, f"{era}_xgb_{bkg_cat}.json"))
    _save_feature_importance(bdt, os.path.join(outdir, f"{era}_featureimportance_{bkg_cat}.txt"), bkg_cat)
    _save_tmva(bdt, os.path.join(outdir, f"{era}_tmva101_{bkg_cat}.root"), x_tr.shape[1])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--era", default="12022")
    p.add_argument("--base-dir", default=None,
                   help="Default: DataPreparation_Files/<era>/")
    p.add_argument("--bkg-cats", nargs="*",
                   default=["DY", "EWK", "TT", "VV", "VVV"],
                   help='Use ["background"] to train against the merged background file.')
    p.add_argument("--outdir", default="Outputs/Training")
    return p.parse_args()


def main():
    args = parse_args()
    base_dir = args.base_dir or os.path.join("DataPreparation_Files", str(args.era))
    if not os.path.isdir(base_dir):
        raise NotADirectoryError(base_dir)
    print(f"[setup] era={args.era}  base_dir={base_dir}  cats={args.bkg_cats}")
    for cat in args.bkg_cats:
        try:
            train_one(args.era, base_dir, cat, outdir_root=args.outdir)
        except FileNotFoundError as e:
            print(f"[WARN] skipping {cat}: missing {e}")


if __name__ == "__main__":
    main()
