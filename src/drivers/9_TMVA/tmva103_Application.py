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

import argparse
import array
import os

import numpy as np
import ROOT

from tmva100_DataPreparation import VARIABLES


def parse_args():
    parser = argparse.ArgumentParser(
        description="Apply five held-out-fold TMVA RBDT models."
    )
    parser.add_argument("--input", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--treename", default="mva")
    parser.add_argument(
        "--model-pattern",
        required=True,
        help=(
            "Model path containing {fold}, for example "
            "Outputs/Training/2024/"
            "2024_tmva101_background_fold{fold}.root"
        ),
    )
    parser.add_argument("--model-name", default="myBDT")
    return parser.parse_args()


def main():
    args = parse_args()

    ROOT.gSystem.Load("libTMVA")

    models = {}
    for fold in range(5):
        model_path = args.model_pattern.format(fold=fold)

        if (
            not os.path.isfile(model_path)
            or os.path.getsize(model_path) < 1000
        ):
            raise SystemExit(
                f"[tmva103] Missing or incomplete model: {model_path}"
            )

        models[fold] = ROOT.TMVA.Experimental.RBDT[""](
            args.model_name,
            model_path,
        )
        print(f"[tmva103] Loaded fold {fold}: {model_path}")

    input_files = ROOT.std.vector("string")()
    for path in args.input:
        if not os.path.isfile(path):
            raise SystemExit(f"[tmva103] Missing input: {path}")
        input_files.push_back(path)

    df = ROOT.RDataFrame(args.treename, input_files)
    columns = {str(c) for c in df.GetColumnNames()}

    missing = [v for v in VARIABLES if v not in columns]
    if missing:
        raise SystemExit(
            f"[tmva103] Missing BDT input variables: {missing}"
        )

    if "fold" not in columns:
        raise SystemExit("[tmva103] Input tree has no fold branch")

    requested = list(VARIABLES) + ["fold"]
    arrays = df.AsNumpy(requested)

    X = np.vstack([
        np.asarray(arrays[v], dtype=np.float32)
        for v in VARIABLES
    ]).T

    folds = np.asarray(arrays["fold"], dtype=np.int32)

    if not np.all(np.isin(folds, np.arange(5))):
        bad = np.unique(folds[~np.isin(folds, np.arange(5))])
        raise SystemExit(
            f"[tmva103] Invalid fold values encountered: {bad}"
        )

    if not np.all(np.isfinite(X)):
        nbad = np.count_nonzero(~np.all(np.isfinite(X), axis=1))
        raise SystemExit(
            f"[tmva103] Found {nbad} events with nonfinite BDT inputs"
        )

    scores = np.full(X.shape[0], np.nan, dtype=np.float32)

    for fold in range(5):
        mask = folds == fold
        count = int(np.count_nonzero(mask))

        if count == 0:
            print(f"[tmva103] Fold {fold}: no events")
            continue

        scores[mask] = np.ravel(
            models[fold].Compute(X[mask]).astype(
                np.float32,
                copy=False,
            )
        )

        print(f"[tmva103] Fold {fold}: scored {count} events")

    if not np.all(np.isfinite(scores)):
        raise SystemExit("[tmva103] Some events were not scored")

    chain = ROOT.TChain(args.treename)
    for path in args.input:
        chain.Add(path)

    if chain.GetEntries() != scores.shape[0]:
        raise SystemExit(
            f"[tmva103] Entry mismatch: "
            f"chain={chain.GetEntries()} scores={scores.shape[0]}"
        )

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    output_file = ROOT.TFile(args.output, "RECREATE")
    if not output_file or output_file.IsZombie():
        raise SystemExit(
            f"[tmva103] Cannot create output: {args.output}"
        )

    output_tree = chain.CloneTree(0)

    score_buffer = array.array("f", [0.0])
    output_tree.Branch(
        "bdt_score",
        score_buffer,
        "bdt_score/F",
    )

    for entry in range(scores.shape[0]):
        chain.GetEntry(entry)
        score_buffer[0] = float(scores[entry])
        output_tree.Fill()

    output_file.cd()
    output_tree.Write()
    output_file.Write()
    output_file.Close()

    print(
        f"[tmva103] Wrote {args.output}: "
        f"{scores.shape[0]} out-of-fold scores"
    )


if __name__ == "__main__":
    main()
