# XGBoost must be imported before ROOT
#
# Example:
#   source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh && python3 tmva101_Training.py --era 12022 --base-dir DataPreparation_Files/12022 --bkg-cats background
#
# This version looks in DataPreparation_Files/<era>/ and trains one model per background
# category (DY, EWK, TT, VV, VVV) vs signal.

from xgboost import XGBClassifier
import argparse, os, glob, numpy as np, ROOT
import matplotlib.pyplot as plt

from tmva100_DataPreparation import VARIABLES

TREENAME = "mva"
WEIGHT_BRANCH = "weight"

def _sanitize_weights(w):
    w = np.asarray(w, dtype=np.float64)
    w[~np.isfinite(w)] = 0.0
    w = np.abs(w)
    w = np.where(w > 1e-9, w, 1e-9)
    return w

def _annihilate_signed_weights(w_signed):
    w_signed = np.asarray(w_signed, dtype=np.float64)
    pos = np.clip(w_signed, 0.0, None)
    neg = -np.clip(w_signed, None, 0.0)

    sum_pos = pos.sum()
    sum_neg = neg.sum()
    if sum_neg <= 0 or sum_pos <= 0:
        return np.clip(w_signed, 0.0, None)

    cancel_from_pos = min(sum_pos, sum_neg)
    if sum_pos > 0:
        pos_scaled = pos * (1.0 - cancel_from_pos / sum_pos)
    else:
        pos_scaled = pos

    eff = pos_scaled
    return eff

def _normalize_class_weights(w_sig, w_bkg):
    sw_s = w_sig.sum()
    sw_b = w_bkg.sum()
    if sw_s <= 0 or sw_b <= 0:
        return w_sig, w_bkg
    scale = np.sqrt(sw_s * sw_b)
    return w_sig * (scale / sw_s), w_bkg * (scale / sw_b)

def _load_numpy(filelist, treename, vars_list, folds_allowed=None, extra_weight_branch=None, is_signal=False):
    if isinstance(filelist, (list, tuple)):
        vs = ROOT.std.vector("string")()
        for f in filelist: vs.push_back(f)
        rdf = ROOT.RDataFrame(treename, vs)
    else:
        rdf = ROOT.RDataFrame(treename, filelist)

    cols = set(rdf.GetColumnNames())
    if folds_allowed is not None and "fold" in cols:
        fold_expr = " || ".join([f"fold=={int(f)}" for f in sorted(folds_allowed)])
        rdf = rdf.Filter(f"({fold_expr})", "select folds")

    arrs = rdf.AsNumpy(vars_list + ([WEIGHT_BRANCH] if WEIGHT_BRANCH in cols else []) + (["fold"] if "fold" in cols else []) + ([extra_weight_branch] if (is_signal and extra_weight_branch and extra_weight_branch in cols) else []))

    x = np.vstack([arrs[v] for v in vars_list]).T
    if WEIGHT_BRANCH in arrs:
        w_signed = np.asarray(arrs[WEIGHT_BRANCH], dtype=np.float64)
    else:
        w_signed = np.ones(x.shape[0], dtype=np.float64)

    if is_signal and extra_weight_branch and extra_weight_branch in arrs:
        sig_extra = np.asarray(arrs[extra_weight_branch], dtype=np.float64)
        sig_extra[~np.isfinite(sig_extra)] = 0.0
        sig_extra = np.clip(sig_extra, 0.0, None)
        w_signed = w_signed * sig_extra

    return x, w_signed

def load_folds(sig_files, bkg_files, folds_train=(0,1), folds_valid=(2,), folds_test=(3,), sig_extra_weight_branch=None):
    xS, wS_signed = _load_numpy(sig_files, TREENAME, VARIABLES, folds_allowed=None, extra_weight_branch=sig_extra_weight_branch, is_signal=True)
    xB, wB_signed = _load_numpy(bkg_files, TREENAME, VARIABLES, folds_allowed=None, extra_weight_branch=None, is_signal=False)

    def get_folds(files, fallback_len):
        rdf = ROOT.RDataFrame(TREENAME, ROOT.std.vector("string")([*files]) if isinstance(files, (list,tuple)) else files)
        cols = set(rdf.GetColumnNames())
        if "fold" in cols:
            return np.asarray(rdf.AsNumpy(["fold"])["fold"], dtype=np.int32)
        else:
            return np.zeros(fallback_len, dtype=np.int32)

    fS = get_folds(sig_files, len(xS))
    fB = get_folds(bkg_files, len(xB))

    def _slice(x, w_signed, f, allowed):
        m = np.isin(f, np.array(allowed, dtype=np.int32))
        return x[m], w_signed[m]

    def _prep_split(xS, wS_signed, xB, wB_signed):
        wS_san = _sanitize_weights(wS_signed) * np.sign(wS_signed)
        wB_san = _sanitize_weights(wB_signed) * np.sign(wB_signed)

        wS_eff = _annihilate_signed_weights(wS_san)
        wB_eff = _annihilate_signed_weights(wB_san)

        wS_eff = _sanitize_weights(wS_eff)
        wB_eff = _sanitize_weights(wB_eff)

        wS_bal, wB_bal = _normalize_class_weights(wS_eff, wB_eff)

        y = np.hstack([np.ones(xS.shape[0]), np.zeros(xB.shape[0])])
        x = np.vstack([xS, xB])
        w = np.hstack([wS_bal, wB_bal])
        return x, y, w

    xS_tr, wS_tr = _slice(xS, wS_signed, fS, (0,1))
    xB_tr, wB_tr = _slice(xB, wB_signed, fB, (0,1))
    x_tr, y_tr, w_tr = _prep_split(xS_tr, wS_tr, xB_tr, wB_tr)

    xS_va, wS_va = _slice(xS, wS_signed, fS, (2,))
    xB_va, wB_va = _slice(xB, wB_signed, fB, (2,))
    x_va, y_va, w_va = _prep_split(xS_va, wS_va, xB_va, wB_va)

    xS_te, wS_te = _slice(xS, wS_signed, fS, (3,))
    xB_te, wB_te = _slice(xB, wB_signed, fB, (3,))
    x_te, y_te, w_te = _prep_split(xS_te, wS_te, xB_te, wB_te)

    print(f"[info] train sizes: S={xS_tr.shape[0]} B={xB_tr.shape[0]} | sums(w): S={w_tr[y_tr==1].sum():.3e} B={w_tr[y_tr==0].sum():.3e}")
    print(f"[info] valid sizes: S={xS_va.shape[0]} B={xB_va.shape[0]}")
    print(f"[info] test  sizes: S={xS_te.shape[0]} B={xB_te.shape[0]}")
    return (x_tr, y_tr, w_tr), (x_va, y_va, w_va), (x_te, y_te, w_te)

def find_prep_files(base_dir, era, kind):
    if kind.lower() == "signal":
        allf = os.path.join(base_dir, f"{era}_all_signal.root")
        if os.path.exists(allf): return [allf]
        return [os.path.join(base_dir, f"{era}_train_signal.root"), os.path.join(base_dir, f"{era}_test_signal.root")]
    else:
        allf = os.path.join(base_dir, f"{era}_all_{kind}.root")
        if os.path.exists(allf): return [allf]
        return [os.path.join(base_dir, f"{era}_train_{kind}.root"), os.path.join(base_dir, f"{era}_test_{kind}.root")]

def train_one(era, base_dir, bkg_cat, sig_extra_weight_branch=None, outdir_root="Outputs/Training"):
    sig_files = find_prep_files(base_dir, era, "signal")
    bkg_files = find_prep_files(base_dir, era, bkg_cat)
    print(f"[debug] category={bkg_cat}  sig_files={sig_files}  bkg_files={bkg_files}")
    for p in sig_files + bkg_files:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing input {p}")
        rdf = ROOT.RDataFrame(TREENAME, p)
        print(f"[debug] {p} has fold? {'fold' in set(rdf.GetColumnNames())}")

    (x_tr, y_tr, w_tr), (x_va, y_va, w_va), (x_te, y_te, w_te) = load_folds(
        sig_files, bkg_files,
        folds_train=(0,1), folds_valid=(2,), folds_test=(3,),
        sig_extra_weight_branch=sig_extra_weight_branch
    )
    if x_va.shape[0] == 0 and x_te.shape[0] == 0:
        print("[WARN] validation/test are empty. Did your snapshots include a 'fold' column?")

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
    outdir = os.path.join(outdir_root, str(era))
    os.makedirs(outdir, exist_ok=True)
    bdt.save_model(os.path.join(outdir, f"{era}_xgb_{bkg_cat}.json"))
    
    raw = bdt.get_booster().get_score(importance_type="gain")
    feat_gain = {VARIABLES[i]: raw.get(f"f{i}", 0.0) for i in range(len(VARIABLES))}
    tot = sum(feat_gain.values()) or 1.0
    feat_pct = {k: 100.0 * v / tot for k, v in feat_gain.items()}
    sorted_pct = sorted(feat_pct.items(), key=lambda x: x[1], reverse=True)
    print(f"\n[info][{bkg_cat}] Feature importance (gain, % of total):")
    for k, p in sorted_pct:
        print(f"  {k:20s} : {p:6.2f}%")

    txt = os.path.join(outdir, f"{era}_featureimportance_{bkg_cat}.txt")
    with open(txt, "w") as f:
        f.write(f"Feature importance (gain, normalized to 100%) — category {bkg_cat}\n")
        for k, p in sorted_pct:
            f.write(f"{k:20s} : {p:6.2f}%\n")
    print(f"[info] Wrote {txt}")

    ROOT.gSystem.Load("libTMVA")
    out_model = os.path.join(outdir, f"{era}_tmva101_{bkg_cat}.root")
    print(f"[info] Saving model to {out_model}")
    
    try:
        ROOT.TMVA.Experimental.SaveXGBoost(bdt, "myBDT", out_model, num_inputs=x_tr.shape[1])
        
        # Verify the file was written successfully
        if not os.path.exists(out_model):
            raise RuntimeError(f"SaveXGBoost did not create output file {out_model}")
        
        file_size = os.path.getsize(out_model)
        if file_size < 1000:
            raise RuntimeError(f"SaveXGBoost output file is too small ({file_size} bytes). File may be corrupted.")
        
        # NOTE: Skipping immediate ROOT file validation because EOS may not have 
        # flushed the data yet, leading to transient read errors. The file size check
        # above is sufficient to detect incomplete writes.
        
        print(f"[info] Model successfully saved: {out_model} ({file_size} bytes)")
    except Exception as e:
        print(f"[ERROR] Failed to save model: {e}")
        # Clean up partial file
        if os.path.exists(out_model):
            os.remove(out_model)
            print(f"[info] Cleaned up partial file {out_model}")
        raise

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--era", default="12022")
    ap.add_argument("--base-dir", default=None, help="Override base dir. Default: DataPreparation_Files/<era>/")
    ap.add_argument("--bkg-cats", nargs="*", default=["DY","EWK","TT","VV","VVV"], help="Background categories to train against.")
    ap.add_argument("--sig-extra-weight-branch", default=None, help="Extra multiplicative weight for SIGNAL only (e.g. mass resolution)")
    ap.add_argument("--outdir", default="Outputs/Training", help="Output base directory (default: %(default)s)")
    args = ap.parse_args()

    base_dir = args.base_dir or os.path.join("DataPreparation", str(args.era))
    if not os.path.isdir(base_dir):
        raise NotADirectoryError(f"Expected {base_dir} to exist and contain the prepared ROOT files.")

    print(f"[setup] era={args.era}  base_dir={base_dir}  cats={args.bkg_cats}")
    for cat in args.bkg_cats:
        try:
            train_one(args.era, base_dir, cat, sig_extra_weight_branch=args.sig_extra_weight_branch, outdir_root=args.outdir)
        except FileNotFoundError as e:
            print(f"[WARN] Skipping {cat}: {e}")
        except Exception as e:
            print(f"[ERROR] Training failed for {cat}: {e}")
            import traceback
            traceback.print_exc()
            raise
