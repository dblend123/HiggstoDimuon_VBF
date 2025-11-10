#!/usr/bin/env python3
# tmva100_DataPreparation.py  (split backgrounds)

# Example:
# python3 tmva100_DataPreparation.py  --era 12022 --treename-in mva --signal-subdir VBF --weight-branch weight --make-all --keep-branches m_mumu

import os, glob, argparse
import ROOT

# --------------------------
# Config / CLI
# --------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Prepare train/test ROOT files (mva tree) for TMVA/XGBoost. Backgrounds are written per-process: DY, TTbar, VV, VVV.")
    p.add_argument("--era", required=True, help="Era string, e.g. 12022")
    p.add_argument("--base", default="/eos/user/d/dblend/Analyzed_SkimmedFiles/{era}/Rochester/",
                   help="Base path pattern with {era} placeholder (default: %(default)s)")
    p.add_argument("--signal-name", default="VBF_all.root",
                   help="Signal ROOT file name pattern (default: %(default)s)")
    p.add_argument("--signal-subdir", default="VBF",
                   help="Subdirectory under base where signal file lives (default: %(default)s)")
    # We keep bg-subdirs for compatibility; we will partition them into DY/TTbar/VV/VVV.
    p.add_argument("--bg-subdirs", nargs="*", default=[
        "DY",
        "EWK",
        "TT",
        "WWW","WWZ","WZZ","ZZZ",
        "WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"
    ], help="Background subdirectory names under the base path (will be split into DY/TTbar/VV/VVV).")
    p.add_argument("--treename-in", default="Events", help="Input tree name (default: %(default)s)")
    p.add_argument("--treename-out", default="mva", help="Output tiny tree name (default: %(default)s)")
    p.add_argument("--weight-branch", default="weight", help="Column name to write as per-event weight. Use 'None' to skip. (default: %(default)s)")
    p.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    p.add_argument("--force-define", action="store_true", help="Force recomputation of VARIABLES even if present")
    p.add_argument("--make-all", action="store_true", help="Also snapshot *_all_*.root that include all folds")
    p.add_argument("--keep-branches", default="", help="Comma-separated extra branches to copy into snapshots (not used for training)")
    return p.parse_args()

# --------------------------
# Variables for the BDT
# --------------------------
VARIABLES = [
    "z", "pt_balance", "mu1_pt_over_m", "mu2_pt_over_m", "mu1_eta", "mu2_eta", "y_mumu", "pt_mumu",
    "jet1_eta", "deta_jj", "dphi_jj", "njets", "pt_j1", "pt_j2", "dphi_mumu_jj", "deta_mumu_jj", "phi_CS", "costheta_CS",
]

# --------------------------
# Selection & feature defs
# --------------------------
def filter_events(df):
    # Put any baseline selection here if needed; otherwise leave as no-op.
    return df

ROOT.gInterpreter.Declare('#include "TVector2.h"')

def define_variables(df):
    # Zeppenfeld: |y_mumu - 0.5*(y_j1 + y_j2)|
    df = df.Define("z", "fabs(y_mumu - 0.5f*(y_j1 + y_j2))")

    # pT-balance = | vector sum pT(mumu+j1+j2) | / (scalar sum pT)
    df = (df.Define("px_sum", "pt_mumu*cos(phi_mumu) + pt_j1*cos(phi_j1) + pt_j2*cos(phi_j2)")
            .Define("py_sum", "pt_mumu*sin(phi_mumu) + pt_j1*sin(phi_j1) + pt_j2*sin(phi_j2)")
            .Define("pt_vecsum", "sqrt(px_sum*px_sum + py_sum*py_sum)")
            .Define("pt_scalar", "pt_mumu + pt_j1 + pt_j2")
            .Define("pt_balance", "pt_scalar>0 ? float(pt_vecsum/pt_scalar) : 0.f"))

    # Extra inputs
    df = (df.Define("mu1_pt_over_m", "m_mumu>0 ? float(mu1_pt/m_mumu) : 0.f")
            .Define("mu2_pt_over_m", "m_mumu>0 ? float(mu2_pt/m_mumu) : 0.f")
            .Define("mu1_eta", "float(mu1_eta)")
            .Define("mu2_eta", "float(mu2_eta)")
            .Define("jet1_eta", "float(eta_j1)")
            .Define("m_mumu",  "float(m_mumu)")
            .Define("deta_jj", "fabs(eta_j1 - eta_j2)")
            .Define("dphi_jj", "fabs(TVector2::Phi_mpi_pi(phi_j1 - phi_j2))")
            .Define("y_mumu",  "float(y_mumu)")
            .Define("pt_j1",   "float(pt_j1)")
            .Define("pt_j2",   "float(pt_j2)"))

    # Weight (optional): expose eventWeight as 'weight'; else unit weights
    cols = set(df.GetColumnNames())
    if "eventWeight" in cols:
        df = df.Define("weight", "float(eventWeight)")
    else:
        df = df.Define("weight", "1.f")
    return df

# --------------------------
# Helpers
# --------------------------
def add_fold(df):
    cols = set(df.GetColumnNames())
    if {"run","luminosityBlock","event"}.issubset(cols):
        return df.Define(
            "fold",
            "( ( ( (int)run * 1315423911 + (int)luminosityBlock * 2654435761u + (int)event ) & 0x7fffffff ) % 4 )"
        )
    elif "event" in cols:
        return df.Define("fold", "((int)event % 4)")
    else:
        return df.Define("fold", "int(rdfentry_ % 4)")

def collect_root_files(base_dir, subdirs):
    files = []
    for sd in subdirs:
        d = os.path.join(base_dir, sd)
        files.extend(sorted(glob.glob(os.path.join(d, "*_all.root"))))
    return files

def build_df(files, treename_in, force_define=False):
    flist = ROOT.std.vector("string")()
    for f in files: flist.push_back(f)
    df = ROOT.RDataFrame(treename_in, flist)
    cols = set(df.GetColumnNames())
    if (not force_define) and set(VARIABLES).issubset(cols):
        return df
    df = filter_events(df)
    df = define_variables(df)
    return df

def add_split_flags(df):
    cols = set(df.GetColumnNames())
    if {"run","luminosityBlock","event"}.issubset(cols):
        df = df.Define(
            "ev_hash",
            "((ULong64_t)run * 1315423911ull) ^ ((ULong64_t)luminosityBlock * 2654435761ull) ^ (ULong64_t)event"
        ).Define("is_train", "(ev_hash & 1ull) == 0ull")
    else:
        df = df.Define("is_train", "(rdfentry_ % 2) == 0")
    return df

def snapshot(df, outname, treename_out, weight_branch, keep_extra):
    cols_out = ROOT.std.vector("string")()
    for v in VARIABLES:
        cols_out.push_back(v)
    if weight_branch and (weight_branch in df.GetColumnNames()):
        cols_out.push_back(weight_branch)
    elif weight_branch and weight_branch != "None":
        print(f"[WARN] requested weight branch '{weight_branch}' not present; writing without it.")
    if "fold" in df.GetColumnNames():
        cols_out.push_back("fold")
    for b in keep_extra:
        if b in df.GetColumnNames():
            cols_out.push_back(b)
    df.Snapshot(treename_out, outname, cols_out)

# ----------------------------------
# Combine backgrounds into one file 
# ----------------------------------
def chain_and_snapshot(file_list, out_path, treename, weight_branch, keep_extra):
    if not file_list: 
        print(f"[WARN] no files to combine for {out_path}")
        return
    vs = ROOT.std.vector("string")()
    for f in file_list:
        if os.path.exists(f): vs.push_back(f)
    if vs.empty():
        print(f"[WARN] none of the inputs exist for {out_path}")
        return
    rdf = ROOT.RDataFrame(treename, vs)
    snapshot(rdf, out_path, treename, weight_branch, KEEP)


# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    args = parse_args()
    KEEP = [b.strip() for b in args.keep_branches.split(",") if b.strip()] 

    base_dir = args.base.format(era=args.era)
    outdir = f"DataPreparation_Files/{args.era}"
    os.makedirs(outdir, exist_ok=True)

    # ----- SIGNAL -----
    signal_file = os.path.join(base_dir, args.signal_subdir, args.signal_name.format(era=args.era))
    print(">>> Preparing SIGNAL from:", signal_file)
    if not os.path.isfile(signal_file):
        raise FileNotFoundError("Signal file not found: {}".format(signal_file))

    df_sig = build_df([signal_file], args.treename_in, force_define=args.force_define)
    df_sig = add_split_flags(df_sig)
    rep_sig = df_sig.Report()

    # Write signal train/test (unchanged)
    train_sig_out = os.path.join(outdir, f"{args.era}_train_signal.root")
    test_sig_out  = os.path.join(outdir, f"{args.era}_test_signal.root")
    snapshot(df_sig.Filter("is_train"), train_sig_out, args.treename_out, args.weight_branch, KEEP)
    snapshot(df_sig.Filter("!is_train"), test_sig_out,  args.treename_out, args.weight_branch, KEEP)
    rep_sig.Print()

    # Partition background subdirs into DY / TTbar / VV / VVV
    all_bgs = list(dict.fromkeys(args.bg_subdirs))  # unique, preserve order

    def match_any(s, keys):
        sl = s.lower()
        return any(k.lower() in sl for k in keys)

    dy_subdirs   = [sd for sd in all_bgs if match_any(sd, ["dy"])]
    ewk_subdirs  = [sd for sd in all_bgs if match_any(sd, ["ewk"])]
    tt_subdirs   = [sd for sd in all_bgs if match_any(sd, ["tt"])]
    # diboson: WW/WZ/ZZ channels, excluding explicit tribosons
    vv_subdirs   = [sd for sd in all_bgs if match_any(sd, ["ww","wz","zz"]) and sd not in ["WWW","WWZ","WZZ","ZZZ"]]
    vvv_subdirs  = [sd for sd in all_bgs if sd in ["WWW","WWZ","WZZ","ZZZ"]]

    CATS = [
        ("DY",      dy_subdirs),
        ("EWK",     ewk_subdirs),
        ("TT",      tt_subdirs),
        ("VV",      vv_subdirs),
        ("VVV",     vvv_subdirs),
    ]

    print(">>> Background category mapping:")
    for name, subs in CATS:
        print(f"    {name}: {subs}")

    # ----- BACKGROUNDS (per-category) -----
    for cat_name, subdirs in CATS:
        files = collect_root_files(base_dir, subdirs)
        print(f">>> Preparing BACKGROUND[{cat_name}] from {len(files)} files under {base_dir}")
        if not files:
            print(f"[WARN] No files found for {cat_name}; skipping this category.")
            continue

        df_cat = build_df(files, args.treename_in, force_define=args.force_define)
        df_cat = add_split_flags(df_cat)
        df_cat = add_fold(df_cat)
        rep_cat = df_cat.Report()

        # Optionally print a few feature stats to confirm content
        try:
            for v in ["pt_mumu","njets","z"]:
                if v in df_cat.GetColumnNames():
                    m = float(df_cat.Mean(v).GetValue()); s = float(df_cat.StdDev(v).GetValue())
                    print(f"[{cat_name}] {v}: mean={m:.3f} std={s:.3f}")
        except Exception as e:
            print(f"[{cat_name}] stats skipped: {e}")

        train_out = os.path.join(outdir, f"{args.era}_train_{cat_name}.root")
        test_out  = os.path.join(outdir, f"{args.era}_test_{cat_name}.root")
        snapshot(df_cat.Filter("is_train"), train_out, args.treename_out, args.weight_branch, KEEP)
        snapshot(df_cat.Filter("!is_train"), test_out,  args.treename_out, args.weight_branch, KEEP)

        if args.make_all:
            all_out = os.path.join(outdir, f"{args.era}_all_{cat_name}.root")
            snapshot(df_cat, all_out, args.treename_out, args.weight_branch, KEEP)
            print("      wrote:", all_out)

        rep_cat.Print()
    
    # --- gather per-category outputs we just wrote ---
    cats = [name for name, _ in CATS]  # DY, EWK, TT, VV, VVV, etc.
    train_bkgs = [os.path.join(outdir, f"{args.era}_train_{c}.root") for c in cats]
    test_bkgs  = [os.path.join(outdir, f"{args.era}_test_{c}.root")  for c in cats]
    all_bkgs   = [os.path.join(outdir, f"{args.era}_all_{c}.root")   for c in cats]

    # --- make overall background snapshots ---
    chain_and_snapshot(train_bkgs,
                       os.path.join(outdir, f"{args.era}_train_background.root"),
                       args.treename_out, args.weight_branch, KEEP)

    chain_and_snapshot(test_bkgs,
                       os.path.join(outdir, f"{args.era}_test_background.root"),
                       args.treename_out, args.weight_branch, KEEP)

    if args.make_all:
        chain_and_snapshot(all_bkgs,
                           os.path.join(outdir, f"{args.era}_all_background.root"),
                           args.treename_out, args.weight_branch, KEEP)

    print("Done. Wrote signal and per-category background files into", outdir)
