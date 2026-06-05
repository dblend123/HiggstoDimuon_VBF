#!/usr/bin/env python3
"""
tmva103_Application.py

Apply a trained TMVA RBDT model to any ROOT file that contains the 'mva' tree
(signal, background, or data) and write a clone with one new branch:
    bdt_score / F

Example (MC):
    python3 tmva103_Application.py \
        --input DataPreparation_Files/12022/12022_all_signal.root \
        --output Outputs/12022/12022_scored_signal.root \
        --model-file Outputs/Training/12022/12022_tmva101_background.root

Example (Data):
    python3 tmva103_Application.py \
        --input /eos/.../12022/MC/DATA/DATA_all.root \
        --output Outputs/12022/12022_data_scored_all.root \
        --model-file Outputs/Training/12022/12022_tmva101_background.root
"""

import argparse, array, os
import numpy as np
import ROOT

from tmva100_DataPreparation import VARIABLES  # MUST match training order


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--treename", default="mva")
    p.add_argument("--model-file", default="tmva101.root")
    p.add_argument("--model-name", default="myBDT")
    return p.parse_args()


def main():
    args = parse_args()

    ROOT.gSystem.Load("libTMVA")
    if not os.path.exists(args.model_file) or os.path.getsize(args.model_file) < 1000:
        raise SystemExit(f"[tmva103] Model missing or incomplete: {args.model_file}")
    model = ROOT.TMVA.Experimental.RBDT[""](args.model_name, args.model_file)

    flist = ROOT.std.vector("string")()
    for f in args.input:
        flist.push_back(f)
    df = ROOT.RDataFrame(args.treename, flist)

    missing = [v for v in VARIABLES if v not in list(df.GetColumnNames())]
    if missing:
        raise SystemExit(f"[tmva103] Missing variables in input tree: {missing}")

    arrays = df.AsNumpy(VARIABLES)
    X = np.vstack([arrays[v] for v in VARIABLES]).T.astype(np.float32, copy=False)
    scores = model.Compute(X).astype(np.float32, copy=False)
    n = scores.shape[0]

    chain = ROOT.TChain(args.treename)
    for f in args.input:
        chain.Add(f)
    if chain.GetEntries() != n:
        raise SystemExit(f"[tmva103] Entry mismatch: chain={chain.GetEntries()} vs scores={n}")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    out = ROOT.TFile(args.output, "RECREATE")
    outtree = chain.CloneTree(0)
    buf = array.array('f', [0.0])
    outtree.Branch("bdt_score", buf, "bdt_score/F")
    for i in range(n):
        chain.GetEntry(i)
        buf[0] = float(scores[i])
        outtree.Fill()
    outtree.Write()
    out.Close()
    print(f"[tmva103] Wrote {args.output} (+branch 'bdt_score' in tree '{args.treename}')")


if __name__ == "__main__":
    main()
