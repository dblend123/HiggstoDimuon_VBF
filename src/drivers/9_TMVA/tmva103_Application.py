#!/usr/bin/env python3
# tmva103_Application.py
# Usage examples:
#   python3 tmva103_Application.py --input test_signal.root --output scored_signal.root
#   python3 tmva103_Application.py --input test_background.root --output scored_bkg.root
#   python3 tmva103_Application.py --input train_signal.root train_background.root --output scored_train.root

import argparse, array
import numpy as np
import ROOT

# --- keep variable order IDENTICAL to training ---
from tmva100_DataPreparation import VARIABLES

def parse_args():
    p = argparse.ArgumentParser(description="Apply TMVA RBDT model to ROOT file(s) with an 'mva' tree.")
    p.add_argument("--input", nargs="+", required=True, help="Input ROOT file(s) that contain the 'mva' tree.")
    p.add_argument("--output", required=True, help="Output ROOT file to write.")
    p.add_argument("--treename", default="mva", help="Tree name (default: mva).")
    p.add_argument("--model-file", default="tmva101.root", help="TMVA model file produced in training.")
    p.add_argument("--model-name", default="myBDT", help="Model name used in SaveXGBoost.")
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Load TMVA + model
    ROOT.gSystem.Load("libTMVA")
    if ROOT.gSystem.AccessPathName(args.model_file, ROOT.kReadPermission):
        ROOT.Error("tmva103", f"Model file not found or unreadable: {args.model_file}")
        raise SystemExit(1)
    model = ROOT.TMVA.Experimental.RBDT[""](args.model_name, args.model_file)

    # Build an RDataFrame to pull NumPy arrays efficiently
    flist = ROOT.std.vector("string")()
    for f in args.input: flist.push_back(f)
    df = ROOT.RDataFrame(args.treename, flist)

    cols = list(df.GetColumnNames())
    missing = [v for v in VARIABLES if v not in cols]
    if missing:
        ROOT.Error("tmva103", "Missing VARIABLES in input tree: " + ", ".join(missing))
        raise SystemExit(2)

    # Fetch all columns and compute scores in one go
    arrays = df.AsNumpy(VARIABLES)  # only VARIABLES needed for scoring
    X = np.vstack([arrays[v] for v in VARIABLES]).T.astype(np.float32, copy=False)
    scores = model.Compute(X).astype(np.float32, copy=False)
    n_entries = scores.shape[0]

    # Now clone the original tree and append the new branch
    chain = ROOT.TChain(args.treename)
    for f in args.input: chain.Add(f)
    if chain.GetEntries() != n_entries:
        ROOT.Error("tmva103", f"Entry mismatch: chain={chain.GetEntries()} vs scores={n_entries}")
        raise SystemExit(3)

    out = ROOT.TFile(args.output, "RECREATE")
    outtree = chain.CloneTree(0)  # clone structure only
    bdt_buf = array.array('f', [0.0])
    outtree.Branch("bdt_score", bdt_buf, "bdt_score/F")

    # Fill
    for i in range(n_entries):
        chain.GetEntry(i)
        bdt_buf[0] = float(scores[i])
        outtree.Fill()

    outtree.Write()
    out.Close()
    print(f"[tmva103] Wrote {args.output} with branch 'bdt_score' in tree '{args.treename}'")
