#!/usr/bin/env python3
"""
tmva101_Training.py

Train one XGBoost BDT per (era, background-category) against the VBF H->mumu
signal, using the prepared snapshots from tmva100_DataPreparation.

Five-fold out-of-fold training: each model uses three folds for
training, one for validation, and one held-out fold for testing.

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
def _training_weights(w):
    """XGBoost requires nonnegative sample weights."""
    w = np.asarray(w, dtype=np.float64).copy()
    w[~np.isfinite(w)] = 0.0
    return np.abs(w)

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
    wS = _training_weights(wS)
    wB = _training_weights(wB)
    wS, wB = _balance(wS, wB)

    y = np.hstack([
        np.ones(xS.shape[0], dtype=np.int32),
        np.zeros(xB.shape[0], dtype=np.int32),
    ])
    x = np.vstack([xS, xB])
    w = np.hstack([wS, wB])

    finite = np.all(np.isfinite(x), axis=1)
    finite &= np.isfinite(w)
    finite &= w >= 0.0

    return x[finite], y[finite], w[finite]

def load_folds(
    sig_files,
    bkg_files,
    folds_train,
    folds_valid,
    folds_test,
):
    xS, wS, fS = _load_numpy(sig_files, VARIABLES)
    xB, wB, fB = _load_numpy(bkg_files, VARIABLES)

    def slc(x, w, f, allowed):
        allowed = np.asarray(allowed, dtype=np.int32)
        mask = np.isin(f, allowed)
        return x[mask], w[mask]

    tr = _prep(
        *slc(xS, wS, fS, folds_train),
        *slc(xB, wB, fB, folds_train),
    )
    va = _prep(
        *slc(xS, wS, fS, folds_valid),
        *slc(xB, wB, fB, folds_valid),
    )
    te = _prep(
        *slc(xS, wS, fS, folds_test),
        *slc(xB, wB, fB, folds_test),
    )

    print(
        f"[info] folds: train={tuple(folds_train)} "
        f"valid={tuple(folds_valid)} test={tuple(folds_test)}"
    )
    print(
        f"[info] train: S={np.count_nonzero(tr[1] == 1)} "
        f"B={np.count_nonzero(tr[1] == 0)} | "
        f"valid: S={np.count_nonzero(va[1] == 1)} "
        f"B={np.count_nonzero(va[1] == 0)} | "
        f"test: S={np.count_nonzero(te[1] == 1)} "
        f"B={np.count_nonzero(te[1] == 0)}"
    )

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

    for path in sig_files + bkg_files:
        if not os.path.exists(path):
            raise FileNotFoundError(path)

    print(
        f"[debug] category={bkg_cat}\n"
        f"        signal={sig_files}\n"
        f"        background={bkg_files}"
    )

    outdir = os.path.join(outdir_root, str(era))
    os.makedirs(outdir, exist_ok=True)

    all_folds = set(range(5))

    for test_fold in range(5):
        valid_fold = (test_fold + 1) % 5
        train_folds = tuple(sorted(all_folds - {test_fold, valid_fold}))

        (x_tr, y_tr, w_tr), (x_va, y_va, w_va), (x_te, y_te, w_te) = load_folds(
            sig_files,
            bkg_files,
            folds_train=train_folds,
            folds_valid=(valid_fold,),
            folds_test=(test_fold,),
        )

        if x_tr.shape[0] == 0:
            raise RuntimeError(
                f"Empty training sample for held-out fold {test_fold}"
            )
        if x_va.shape[0] == 0:
            raise RuntimeError(
                f"Empty validation sample for held-out fold {test_fold}"
            )
        if x_te.shape[0] == 0:
            raise RuntimeError(
                f"Empty test sample for held-out fold {test_fold}"
            )

        print(
            f"[train] model fold={test_fold}: "
            f"train={train_folds}, validation={valid_fold}, "
            f"test={test_fold}"
        )

        bdt = _xgb_classifier()
        bdt.fit(
            x_tr,
            y_tr,
            sample_weight=w_tr,
            eval_set=[(x_va, y_va)],
            sample_weight_eval_set=[w_va],
            verbose=False,
        )

        stem = f"{era}_tmva101_{bkg_cat}_fold{test_fold}"

        bdt.save_model(
            os.path.join(
                outdir,
                f"{era}_xgb_{bkg_cat}_fold{test_fold}.json",
            )
        )

        _save_feature_importance(
            bdt,
            os.path.join(
                outdir,
                f"{era}_featureimportance_{bkg_cat}_fold{test_fold}.txt",
            ),
            f"{bkg_cat}_fold{test_fold}",
        )

        _save_tmva(
            bdt,
            os.path.join(outdir, f"{stem}.root"),
            x_tr.shape[1],
        )

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
