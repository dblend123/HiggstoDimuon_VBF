#!/usr/bin/env python3
"""
tmva100_DataPreparation.py

Read NanoAOD-derived analysis ntuples, define the 18 BDT input features, split
each event deterministically into one of 4 folds, and snapshot per-era ROOT
files used by training/testing/scoring (tmva101/102/103).

Backgrounds are written per process category (DY, EWK, TT, VV, VVV) and also
chained into a single "background" file.

Example:
    python3 tmva100_DataPreparation.py --era 12022 \
        --signal-subdir VBF --treename-in mva --make-all --keep-branches m_mumu
"""

import os, glob, argparse
import ROOT

# -----------------------------------------------------------------------------
# BDT input features (order must match tmva101/103)
# -----------------------------------------------------------------------------
VARIABLES = [
    "z", "pt_balance",
    "mu1_pt_over_m", "mu2_pt_over_m", "mu1_eta", "mu2_eta",
    "y_mumu", "pt_mumu",
    "jet1_eta", "deta_jj", "dphi_jj", "njets", "pt_j1", "pt_j2",
    "dphi_mumu_jj", "deta_mumu_jj",
    "phi_CS", "costheta_CS",
]

ROOT.gInterpreter.Declare('#include "TVector2.h"')


def define_variables(df):
    """Define the 18 BDT inputs + an event weight column."""
    df = (df
        .Define("z",          "fabs(y_mumu - 0.5f*(y_j1 + y_j2))")
        .Define("px_sum",     "pt_mumu*cos(phi_mumu) + pt_j1*cos(phi_j1) + pt_j2*cos(phi_j2)")
        .Define("py_sum",     "pt_mumu*sin(phi_mumu) + pt_j1*sin(phi_j1) + pt_j2*sin(phi_j2)")
        .Define("pt_vecsum",  "sqrt(px_sum*px_sum + py_sum*py_sum)")
        .Define("pt_scalar",  "pt_mumu + pt_j1 + pt_j2")
        .Define("pt_balance", "pt_scalar>0 ? float(pt_vecsum/pt_scalar) : 0.f")
        .Define("mu1_pt_over_m", "m_mumu>0 ? float(mu1_pt/m_mumu) : 0.f")
        .Define("mu2_pt_over_m", "m_mumu>0 ? float(mu2_pt/m_mumu) : 0.f")
        .Define("mu1_eta",  "float(mu1_eta)")
        .Define("mu2_eta",  "float(mu2_eta)")
        .Define("jet1_eta", "float(eta_j1)")
        .Define("m_mumu",   "float(m_mumu)")
        .Define("deta_jj",  "fabs(eta_j1 - eta_j2)")
        .Define("dphi_jj",  "fabs(TVector2::Phi_mpi_pi(phi_j1 - phi_j2))")
        .Define("y_mumu",   "float(y_mumu)")
        .Define("pt_j1",    "float(pt_j1)")
        .Define("pt_j2",    "float(pt_j2)"))
    cols = set(df.GetColumnNames())
    df = df.Define("weight", "float(eventWeight)" if "eventWeight" in cols else "1.f")
    return df


def add_fold_and_split(df):
    """Add a deterministic 4-fold hash plus a train/test flag."""
    cols = set(df.GetColumnNames())
    if {"run", "luminosityBlock", "event"}.issubset(cols):
        return (df
            .Define("fold",
                    "( ( ( (int)run * 1315423911 + (int)luminosityBlock * 2654435761u + (int)event ) & 0x7fffffff ) % 4 )")
            .Define("ev_hash",
                    "((ULong64_t)run * 1315423911ull) ^ ((ULong64_t)luminosityBlock * 2654435761ull) ^ (ULong64_t)event")
            .Define("is_train", "(ev_hash & 1ull) == 0ull"))
    if "event" in cols:
        return df.Define("fold", "((int)event % 4)").Define("is_train", "((int)event % 2) == 0")
    return df.Define("fold", "int(rdfentry_ % 4)").Define("is_train", "(rdfentry_ % 2) == 0")


def build_df(files, treename, force_define=False):
    flist = ROOT.std.vector("string")()
    for f in files:
        flist.push_back(f)
    df = ROOT.RDataFrame(treename, flist)
    if (not force_define) and set(VARIABLES).issubset(set(df.GetColumnNames())):
        return df
    return define_variables(df)


def collect_files(base_dir, subdirs):
    out = []
    for sd in subdirs:
        out.extend(sorted(glob.glob(os.path.join(base_dir, sd, "*_all.root"))))
    return out


def snapshot(df, outname, treename, weight_branch, keep_extra):
    present = set(df.GetColumnNames())
    cols = ROOT.std.vector("string")()
    for v in VARIABLES:
        cols.push_back(v)
    if weight_branch and weight_branch in present:
        cols.push_back(weight_branch)
    if "fold" in present:
        cols.push_back("fold")
    for b in keep_extra:
        if b in present:
            cols.push_back(b)
    os.makedirs(os.path.dirname(outname) or ".", exist_ok=True)
    df.Snapshot(treename, outname, cols)
    if not os.path.exists(outname) or os.path.getsize(outname) < 500:
        raise RuntimeError(f"Snapshot incomplete: {outname}")
    print(f"      wrote: {outname} ({os.path.getsize(outname)} bytes)")


def chain_snapshot(files, outname, treename, weight_branch, keep_extra):
    existing = [f for f in files if os.path.exists(f)]
    if not existing:
        print(f"[WARN] no inputs for {outname}; skipping")
        return
    vs = ROOT.std.vector("string")()
    for f in existing:
        vs.push_back(f)
    snapshot(ROOT.RDataFrame(treename, vs), outname, treename, weight_branch, keep_extra)


def categorize(subdirs):
    """Map raw subdirectory names to background categories."""
    has = lambda s, k: k.lower() in s.lower()
    vvv = {"WWW", "WWZ", "WZZ", "ZZZ"}
    return [
        ("DY",  [s for s in subdirs if has(s, "dy")]),
        ("EWK", [s for s in subdirs if has(s, "ewk")]),
        ("TT",  [s for s in subdirs if has(s, "tt")]),
        ("VV",  [s for s in subdirs if (has(s, "ww") or has(s, "wz") or has(s, "zz")) and s not in vvv]),
        ("VVV", [s for s in subdirs if s in vvv]),
    ]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--era", required=True)
    p.add_argument("--base", default="/eos/user/d/dblend/Analyzed_SkimmedFiles/{era}/Rochester/")
    p.add_argument("--signal-name", default="VBF_all.root")
    p.add_argument("--signal-subdir", default="VBF")
    p.add_argument("--bg-subdirs", nargs="*", default=[
        "DY", "EWK", "TT",
        "WWW", "WWZ", "WZZ", "ZZZ",
        "WZto2L2Q", "WZto3LNu", "ZZto2L2Nu", "ZZto2L2Q", "ZZto4L",
    ])
    p.add_argument("--treename-in",  default="Events")
    p.add_argument("--treename-out", default="mva")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--force-define", action="store_true")
    p.add_argument("--make-all", action="store_true",
                   help="Also write {era}_all_*.root that contain all four folds.")
    p.add_argument("--keep-branches", default="",
                   help="Comma-separated extra branches to copy through.")
    return p.parse_args()


def main():
    args = parse_args()
    keep = [b.strip() for b in args.keep_branches.split(",") if b.strip()]
    base_dir = args.base.format(era=args.era)
    outdir = f"DataPreparation_Files/{args.era}"
    os.makedirs(outdir, exist_ok=True)

    # ---- SIGNAL ----
    sig_file = os.path.join(base_dir, args.signal_subdir, args.signal_name.format(era=args.era))
    print(">>> SIGNAL:", sig_file)
    if not os.path.isfile(sig_file):
        raise FileNotFoundError(sig_file)
    df_sig = add_fold_and_split(build_df([sig_file], args.treename_in, args.force_define))
    snapshot(df_sig.Filter("is_train"),  f"{outdir}/{args.era}_train_signal.root", args.treename_out, args.weight_branch, keep)
    snapshot(df_sig.Filter("!is_train"), f"{outdir}/{args.era}_test_signal.root",  args.treename_out, args.weight_branch, keep)
    if args.make_all:
        snapshot(df_sig, f"{outdir}/{args.era}_all_signal.root", args.treename_out, args.weight_branch, keep)

    # ---- BACKGROUNDS (per category) ----
    cats = categorize(list(dict.fromkeys(args.bg_subdirs)))
    print(">>> Background category mapping:")
    for name, subs in cats:
        print(f"    {name}: {subs}")

    for cat, subs in cats:
        files = collect_files(base_dir, subs)
        print(f">>> BKG[{cat}]: {len(files)} files")
        if not files:
            continue
        df = add_fold_and_split(build_df(files, args.treename_in, args.force_define))
        snapshot(df.Filter("is_train"),  f"{outdir}/{args.era}_train_{cat}.root", args.treename_out, args.weight_branch, keep)
        snapshot(df.Filter("!is_train"), f"{outdir}/{args.era}_test_{cat}.root",  args.treename_out, args.weight_branch, keep)
        if args.make_all:
            snapshot(df, f"{outdir}/{args.era}_all_{cat}.root", args.treename_out, args.weight_branch, keep)

    # ---- Combined background snapshots (chain across categories) ----
    names = [n for n, _ in cats]
    chain_snapshot([f"{outdir}/{args.era}_train_{n}.root" for n in names],
                   f"{outdir}/{args.era}_train_background.root", args.treename_out, args.weight_branch, keep)
    chain_snapshot([f"{outdir}/{args.era}_test_{n}.root" for n in names],
                   f"{outdir}/{args.era}_test_background.root", args.treename_out, args.weight_branch, keep)
    if args.make_all:
        chain_snapshot([f"{outdir}/{args.era}_all_{n}.root" for n in names],
                       f"{outdir}/{args.era}_all_background.root", args.treename_out, args.weight_branch, keep)

    print("Done. Outputs in", outdir)


if __name__ == "__main__":
    main()
