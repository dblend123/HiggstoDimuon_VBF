# XGBoost must be imported before ROOT
# source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh && python3 tmva101_Training.py --era 22023

from xgboost import XGBClassifier
import argparse, os, numpy as np, ROOT
import matplotlib.pyplot as plt

from tmva100_DataPreparation import VARIABLES

TREENAME = "mva"
WEIGHT_BRANCH = "weight"
OUTDIR = "Outputs"

def _sanitize_weights(w):
    w = np.asarray(w, dtype=np.float64)
    w[~np.isfinite(w)] = 0.0
    w = np.abs(w)
    w = np.where(w > 1e-9, w, 1e-9)
    return w

def _annihilate_signed_weights(w_signed):
    """
    Pair negative with positive weights within a class and annihilate them.
    Returns non-negative effective weights (float64) of same shape.
    Strategy: cancel total negative mass against positives proportionally.
    """
    w_signed = np.asarray(w_signed, dtype=np.float64)
    pos = np.clip(w_signed, 0.0, None)
    neg = -np.clip(w_signed, None, 0.0)

    sum_pos = pos.sum()
    sum_neg = neg.sum()
    if sum_neg <= 0 or sum_pos <= 0:
        # nothing to do; also handles all-positive or all-negative edge cases
        return np.clip(w_signed, 0.0, None)

    # cancel min(sum_pos, sum_neg) proportionally from the positive side
    cancel_from_pos = min(sum_pos, sum_neg)
    if sum_pos > 0:
        pos_scaled = pos * (1.0 - cancel_from_pos / sum_pos)
    else:
        pos_scaled = pos

    # If negatives exceed positives, residual negatives can't be used -> set to 0
    # (this is standard "annihilation": we do not keep negative effective weights)
    eff = pos_scaled  # residual positives after annihilation
    return eff

def _normalize_class_weights(w_sig, w_bkg):
    """Scale so that sum(w_sig) == sum(w_bkg)."""
    sw_s = w_sig.sum()
    sw_b = w_bkg.sum()
    if sw_s <= 0 or sw_b <= 0:
        return w_sig, w_bkg
    scale = np.sqrt(sw_s * sw_b)
    return w_sig * (scale / sw_s), w_bkg * (scale / sw_b)

def _load_numpy(filelist, treename, vars_list, folds_allowed=None, extra_weight_branch=None, is_signal=False):
    """
    Load arrays from ROOT, select folds, build X matrix and weight vector.
    For signal, multiply by extra weight branch if provided.
    """
    # Support multiple files
    if isinstance(filelist, (list, tuple)):
        rdf = ROOT.RDataFrame(treename, ROOT.std.vector("string")(list(filelist)))
    else:
        rdf = ROOT.RDataFrame(treename, filelist)

    cols = set(rdf.GetColumnNames())
    if folds_allowed is not None and "fold" in cols:
        # Filter on allowed folds
        fold_expr = " || ".join([f"fold=={int(f)}" for f in sorted(folds_allowed)])
        rdf = rdf.Filter(f"({fold_expr})", "select folds")

    # pull variables
    arrs = rdf.AsNumpy(vars_list + ([WEIGHT_BRANCH] if WEIGHT_BRANCH in cols else []) + (["fold"] if "fold" in cols else []) + ([extra_weight_branch] if (is_signal and extra_weight_branch and extra_weight_branch in cols) else []))

    x = np.vstack([arrs[v] for v in vars_list]).T
    # base weights
    if WEIGHT_BRANCH in arrs:
        w_signed = np.asarray(arrs[WEIGHT_BRANCH], dtype=np.float64)
    else:
        w_signed = np.ones(x.shape[0], dtype=np.float64)

    # For signal: multiply by mass-resolution (or any extra) branch if requested
    if is_signal and extra_weight_branch and extra_weight_branch in arrs:
        sig_extra = np.asarray(arrs[extra_weight_branch], dtype=np.float64)
        # sanitize multiplicative factor (non-negative, finite)
        sig_extra[~np.isfinite(sig_extra)] = 0.0
        sig_extra = np.clip(sig_extra, 0.0, None)
        w_signed = w_signed * sig_extra

    return x, w_signed

def load_folds(sig_files, bkg_files, folds_train=(0,1), folds_valid=(2,), folds_test=(3,), sig_extra_weight_branch=None):
    # Load signal/background for train+valid+test, then slice by folds
    xS, wS_signed = _load_numpy(sig_files, TREENAME, VARIABLES, folds_allowed=None, extra_weight_branch=sig_extra_weight_branch, is_signal=True)
    xB, wB_signed = _load_numpy(bkg_files, TREENAME, VARIABLES, folds_allowed=None, extra_weight_branch=None, is_signal=False)

    # Also fetch folds to split (re-load only fold column to avoid big memory hit)
    def get_folds(files):
        rdf = ROOT.RDataFrame(TREENAME, ROOT.std.vector("string")(list(files)) if isinstance(files, (list,tuple)) else files)
        cols = set(rdf.GetColumnNames())
        if "fold" in cols:
            return np.asarray(rdf.AsNumpy(["fold"])["fold"], dtype=np.int32)
        else:
            # fallback: 0 for all -> all go to 'train' if folds not present
            return np.zeros(len(xS if files is sig_files else xB), dtype=np.int32)
    fS = get_folds(sig_files)
    fB = get_folds(bkg_files)

    def _slice(x, w_signed, f, allowed):
        m = np.isin(f, np.array(allowed, dtype=np.int32))
        return x[m], w_signed[m]

    # sanitize & annihilate within each split separately to keep consistency
    def _prep_split(xS, wS_signed, xB, wB_signed):
        # sanitize (finite, abs, eps)
        wS_san = _sanitize_weights(wS_signed) * np.sign(wS_signed)  # keep sign for annihilation
        wB_san = _sanitize_weights(wB_signed) * np.sign(wB_signed)

        # annihilate (per class)
        wS_eff = _annihilate_signed_weights(wS_san)
        wB_eff = _annihilate_signed_weights(wB_san)

        # after annihilation, sanitize to strictly positive eps
        wS_eff = _sanitize_weights(wS_eff)
        wB_eff = _sanitize_weights(wB_eff)

        # class-balance by weight sum (not counts)
        wS_bal, wB_bal = _normalize_class_weights(wS_eff, wB_eff)

        # labels
        y = np.hstack([np.ones(xS.shape[0]), np.zeros(xB.shape[0])])
        x = np.vstack([xS, xB])
        w = np.hstack([wS_bal, wB_bal])
        return x, y, w

    # Train
    xS_tr, wS_tr = _slice(xS, wS_signed, fS, folds_train)
    xB_tr, wB_tr = _slice(xB, wB_signed, fB, folds_train)
    x_tr, y_tr, w_tr = _prep_split(xS_tr, wS_tr, xB_tr, wB_tr)

    # Valid
    xS_va, wS_va = _slice(xS, wS_signed, fS, folds_valid)
    xB_va, wB_va = _slice(xB, wB_signed, fB, folds_valid)
    x_va, y_va, w_va = _prep_split(xS_va, wS_va, xB_va, wB_va)

    # Test
    xS_te, wS_te = _slice(xS, wS_signed, fS, folds_test)
    xB_te, wB_te = _slice(xB, wB_signed, fB, folds_test)
    x_te, y_te, w_te = _prep_split(xS_te, wS_te, xB_te, wB_te)

    # quick diagnostics
    print(f"[info] train sizes: S={xS_tr.shape[0]} B={xB_tr.shape[0]} | sums(w): S={w_tr[y_tr==1].sum():.3e} B={w_tr[y_tr==0].sum():.3e}")
    print(f"[info] valid sizes: S={xS_va.shape[0]} B={xB_va.shape[0]}")
    print(f"[info] test  sizes: S={xS_te.shape[0]} B={xB_te.shape[0]}")
    return (x_tr, y_tr, w_tr), (x_va, y_va, w_va), (x_te, y_te, w_te)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="DataPreparation_Files")
    ap.add_argument("--era", default="12022")
    ap.add_argument("--sig-extra-weight-branch", default=None, help="Extra multiplicative weight for SIGNAL only (e.g. mass resolution)")
    args = ap.parse_args()

    sig_files = [os.path.join(args.dir, f"{args.era}_all_signal.root")] if os.path.exists(os.path.join(args.dir, f"{args.era}_all_signal.root")) else [os.path.join(args.dir, f"{args.era}_train_signal.root"), os.path.join(args.dir, f"{args.era}_test_signal.root")]
    bkg_files = [os.path.join(args.dir, f"{args.era}_all_background.root")] if os.path.exists(os.path.join(args.dir, f"{args.era}_all_background.root")) else [os.path.join(args.dir, f"{args.era}_train_background.root"), os.path.join(args.dir, f"{args.era}_test_background.root")]

    print("[debug] sig_files:", sig_files)
    print("[debug] bkg_files:", bkg_files)
    for p in sig_files + bkg_files:
        rdf = ROOT.RDataFrame(TREENAME, p)
        print(f"[debug] {p} has fold? {'fold' in set(rdf.GetColumnNames())}")


    # Load folds (train {0,1}, valid {2}, test {3}); apply annihilation and class-sum normalization
    (x_tr, y_tr, w_tr), (x_va, y_va, w_va), (x_te, y_te, w_te) = load_folds(
        sig_files, bkg_files,
        folds_train=(0,1), folds_valid=(2,), folds_test=(3,),
        sig_extra_weight_branch=args.sig_extra_weight_branch
    )
    
    # loud guard if folds are missing (everything went to train)
    if x_va.shape[0] == 0 and x_te.shape[0] == 0:
        print("[WARN] validation/test are empty. Did your snapshots include a 'fold' column?")

    # Train with validation monitoring (optional early stopping if you like)
    bdt = XGBClassifier(
        tree_method="hist",
        max_bin=30, 
        max_depth=4,
        learning_rate=0.1,
        n_estimators=1000,
        eval_metric="logloss", 
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
        missing=None
    )
    bdt.fit(x_tr, y_tr, sample_weight=w_tr,
            eval_set=[(x_va, y_va)],
            sample_weight_eval_set=[w_va],
            verbose=False)

    # --- Feature importance (gain) printout (sorted, highest→lowest) ---

    # raw (unitless) gain per feature
    raw = bdt.get_booster().get_score(importance_type="gain")  # keys like 'f0','f1',...
    feat_gain = {VARIABLES[i]: raw.get(f"f{i}", 0.0) for i in range(len(VARIABLES))}

    tot = sum(feat_gain.values()) or 1.0  # avoid div-by-zero
    feat_pct = {k: 100.0 * v / tot for k, v in feat_gain.items()}  # percent

    # sort high→low and print
    sorted_pct = sorted(feat_pct.items(), key=lambda x: x[1], reverse=True)
    print("\n[info] Feature importance (gain, % of total):")
    for k, p in sorted_pct:
        print(f"  {k:20s} : {p:6.2f}%")

    # write to file
    out_txt = os.path.join(OUTDIR, f"{args.era}_featureimportance.txt")
    with open(out_txt, "w") as f:
        f.write("Feature importance (gain, normalized to 100%)\n")
        for k, p in sorted_pct:
            f.write(f"{k:20s} : {p:6.2f}%\n")
    print(f"[info] Wrote {out_txt}")

    ROOT.gSystem.Load("libTMVA")
    print("Training done on", x_tr.shape[0], f"train events. Saving model in {args.era}_tmva101.root")
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, "myBDT", f"{args.era}_tmva101.root", num_inputs=x_tr.shape[1])

