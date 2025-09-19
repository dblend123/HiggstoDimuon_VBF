#!/usr/bin/env python3
# tmva100_DataPreparation.py
# python3 tmva100_DataPreparation.py   --era 12023   --treename-in mva   --signal-subdir VBF   --weight-branch weight

import os, glob, argparse
import ROOT

# --------------------------
# Config / CLI
# --------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Prepare train/test ROOT files (mva tree) for TMVA/XGBoost.")
    p.add_argument("--era", required=True, help="Era string, e.g. 12022")
    p.add_argument("--base", default="/eos/user/d/dblend/Analyzed_SkimmedFiles/{era}/MC/",
                   help="Base path pattern with {era} placeholder (default: %(default)s)")
    p.add_argument("--signal-name", default="VBF_MCRun{era}.root",
                   help="Signal ROOT file name pattern (default: %(default)s)")
    p.add_argument("--signal-subdir", default="VBF",
                   help="Subdirectory under base where signal file lives (default: %(default)s)")
    p.add_argument("--bg-subdirs", nargs="*", default=[
        "DY","EWK","TT","WWW","WWZ","WZto2L2Q","WZto3LNu",
        "WZZ","ZZto2L2Nu","ZZto2L2Q","ZZto4L","ZZZ"
    ], help="Background subdirectory names under the base path.")
    p.add_argument("--treename-in", default="Events", help="Input tree name (default: %(default)s)")
    p.add_argument("--treename-out", default="mva", help="Output tiny tree name (default: %(default)s)")
    p.add_argument("--weight-branch", default="weight",
                   help="Column name to write the per-event weight as. Use 'None' to skip. (default: %(default)s)")
    p.add_argument("--outdir", default=".", help="Output directory (default: current dir)")
    return p.parse_args()

# --------------------------
# Variables for the BDT
# --------------------------
VARIABLES = [
    "z",
    "pt_balance",
    "mu1_pt_over_m", "mu2_pt_over_m",
    "mu1_eta", "mu2_eta", "jet1_eta",
    "m_mumu", "deta_jj", "dphi_jj",
    "y_mumu", "pt_j1", "pt_j2",
]

# --------------------------
# Selection & feature defs
# --------------------------
def filter_events(df):
    # Put any baseline selection here if needed; otherwise leave as no-op.
    return df

def define_variables(df):
    # If TVector2 header isn't autoloaded on your setup, uncomment:
    # ROOT.gInterpreter.Declare('#include "TVector2.h"')

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
def collect_root_files(base_dir, subdirs):
    files = []
    for sd in subdirs:
        d = os.path.join(base_dir, sd)
        files.extend(sorted(glob.glob(os.path.join(d, "*.root"))))
    return files

def build_df(files, treename_in):
    flist = ROOT.std.vector("string")()
    for f in files: flist.push_back(f)
    df = ROOT.RDataFrame(treename_in, flist)

    # If the input tree already has all VARIABLES, don’t redefine them.
    cols = set(df.GetColumnNames())
    if set(VARIABLES).issubset(cols):
        return df

    # Otherwise compute them (original behavior)
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

def snapshot(df, outname, treename_out, weight_branch):
    cols_out = ROOT.std.vector("string")()
    for v in VARIABLES:
        cols_out.push_back(v)
    if weight_branch and (weight_branch in df.GetColumnNames()):
        cols_out.push_back(weight_branch)
    elif weight_branch and weight_branch != "None":
        print("[WARN] requested weight branch '{}' not present; writing without it.".format(weight_branch))
    df.Snapshot(treename_out, outname, cols_out)

# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    args = parse_args()

    base_dir = args.base.format(era=args.era)
    os.makedirs(args.outdir, exist_ok=True)

    # ----- SIGNAL -----
    signal_file = os.path.join(base_dir, args.signal_subdir, args.signal_name.format(era=args.era))
    print(">>> Preparing SIGNAL from:", signal_file)
    if not os.path.isfile(signal_file):
        raise FileNotFoundError("Signal file not found: {}".format(signal_file))

    df_sig = build_df([signal_file], args.treename_in)
    df_sig = add_split_flags(df_sig)
    rep_sig = df_sig.Report()

    train_sig_out = os.path.join(args.outdir, "train_signal.root")
    test_sig_out  = os.path.join(args.outdir, "test_signal.root")

    snapshot(df_sig.Filter("is_train"), train_sig_out, args.treename_out, args.weight_branch)
    snapshot(df_sig.Filter("!is_train"), test_sig_out,  args.treename_out, args.weight_branch)

    rep_sig.Print()

    # ----- BACKGROUNDS -----
    bg_files = collect_root_files(base_dir, args.bg_subdirs)
    print(">>> Preparing BACKGROUND from {} files under {}".format(len(bg_files), base_dir))
    if not bg_files:
        raise RuntimeError("No background ROOT files found under {} subdirs {}".format(base_dir, args.bg_subdirs))

    df_bkg = build_df(bg_files, args.treename_in)
    df_bkg = add_split_flags(df_bkg)
    rep_bkg = df_bkg.Report()

    train_bkg_out = os.path.join(args.outdir, "train_background.root")
    test_bkg_out  = os.path.join(args.outdir, "test_background.root")

    snapshot(df_bkg.Filter("is_train"), train_bkg_out, args.treename_out, args.weight_branch)
    snapshot(df_bkg.Filter("!is_train"), test_bkg_out,  args.treename_out, args.weight_branch)

    rep_bkg.Print()

    print("Done. Wrote:")
    print("  ", train_sig_out)
    print("  ", test_sig_out)
    print("  ", train_bkg_out)
    print("  ", test_bkg_out)

