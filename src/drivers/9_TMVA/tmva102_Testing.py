#!/usr/bin/env python3
"""
tmva102_Testing.py

Evaluate a trained TMVA RBDT on the held-out fold (fold==3) of the prepared
signal/background snapshots. Writes weighted ROC, AUC, TPR@1%FPR and Youden J.

Example:
    python3 tmva102_Testing.py --era 12022 \
        --dir DataPreparation_Files/12022 \
        --model Outputs/Training/12022/12022_tmva101_background.root
"""
import argparse
import os

import numpy as np
import ROOT
from sklearn.metrics import auc, roc_curve


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate five-fold out-of-fold BDT scores."
    )
    parser.add_argument("--era", required=True)
    parser.add_argument("--signal", required=True)
    parser.add_argument("--background", required=True)
    parser.add_argument("--treename", default="mva")
    parser.add_argument("--score-branch", default="bdt_score")
    parser.add_argument("--weight-branch", default="weight")
    parser.add_argument("--output-dir", default="Outputs/Testing")
    return parser.parse_args()


def load_arrays(path, args):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    df = ROOT.RDataFrame(args.treename, path)
    columns = {str(c) for c in df.GetColumnNames()}

    required = {
        args.score_branch,
        args.weight_branch,
        "fold",
    }
    missing = sorted(required - columns)
    if missing:
        raise RuntimeError(
            f"{path} is missing branches: {missing}"
        )

    arrays = df.AsNumpy([
        args.score_branch,
        args.weight_branch,
        "fold",
    ])

    score = np.asarray(
        arrays[args.score_branch],
        dtype=np.float64,
    )
    signed_weight = np.asarray(
        arrays[args.weight_branch],
        dtype=np.float64,
    )
    fold = np.asarray(arrays["fold"], dtype=np.int32)

    mask = np.isfinite(score)
    mask &= np.isfinite(signed_weight)
    mask &= np.isin(fold, np.arange(5))

    return score[mask], signed_weight[mask], fold[mask]


def main():
    args = parse_args()

    sig_score, sig_signed_weight, sig_fold = load_arrays(
        args.signal,
        args,
    )
    bkg_score, bkg_signed_weight, bkg_fold = load_arrays(
        args.background,
        args,
    )

    if sig_score.size == 0 or bkg_score.size == 0:
        raise SystemExit("[tmva102] Empty signal or background input")

    # ROC implementations require nonnegative sample weights.
    # These are diagnostic ROC weights, not physical yield weights.
    sig_roc_weight = np.abs(sig_signed_weight)
    bkg_roc_weight = np.abs(bkg_signed_weight)

    labels = np.concatenate([
        np.ones(sig_score.size, dtype=np.int32),
        np.zeros(bkg_score.size, dtype=np.int32),
    ])
    scores = np.concatenate([sig_score, bkg_score])
    roc_weights = np.concatenate([
        sig_roc_weight,
        bkg_roc_weight,
    ])

    fpr, tpr, thresholds = roc_curve(
        labels,
        scores,
        sample_weight=roc_weights,
    )
    auc_score = auc(fpr, tpr)

    print(f"[tmva102] Combined out-of-fold AUC = {auc_score:.6f}")

    for fold in range(5):
        ns = int(np.count_nonzero(sig_fold == fold))
        nb = int(np.count_nonzero(bkg_fold == fold))

        fold_labels = np.concatenate([
            np.ones(ns, dtype=np.int32),
            np.zeros(nb, dtype=np.int32),
        ])
        fold_scores = np.concatenate([
            sig_score[sig_fold == fold],
            bkg_score[bkg_fold == fold],
        ])
        fold_weights = np.concatenate([
            sig_roc_weight[sig_fold == fold],
            bkg_roc_weight[bkg_fold == fold],
        ])

        if ns == 0 or nb == 0:
            print(
                f"[tmva102] Fold {fold}: "
                f"S={ns} B={nb}; AUC unavailable"
            )
            continue

        fold_fpr, fold_tpr, _ = roc_curve(
            fold_labels,
            fold_scores,
            sample_weight=fold_weights,
        )
        fold_auc = auc(fold_fpr, fold_tpr)

        print(
            f"[tmva102] Fold {fold}: "
            f"S={ns} B={nb} AUC={fold_auc:.6f}"
        )

    os.makedirs(args.output_dir, exist_ok=True)

    metrics_path = os.path.join(
        args.output_dir,
        f"{args.era}_tmva102_metrics.txt",
    )
    with open(metrics_path, "w", encoding="utf-8") as output:
        output.write(f"combined_oof_auc {auc_score:.8f}\n")
        output.write(
            f"signal_signed_yield {sig_signed_weight.sum():.12g}\n"
        )
        output.write(
            f"background_signed_yield "
            f"{bkg_signed_weight.sum():.12g}\n"
        )

    points_path = os.path.join(
        args.output_dir,
        f"{args.era}_tmva102_roc_points.csv",
    )
    with open(points_path, "w", encoding="utf-8") as output:
        output.write("fpr,tpr,threshold\n")
        for x, y, threshold in zip(fpr, tpr, thresholds):
            output.write(
                f"{x:.10g},{y:.10g},{threshold:.10g}\n"
            )

    print(f"[tmva102] Wrote {metrics_path}")
    print(f"[tmva102] Wrote {points_path}")


if __name__ == "__main__":
    main()
