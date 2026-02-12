#!/usr/bin/env python3

import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os

# ============================================================
# Argument parsing
# ============================================================

def parse_args():
    p = argparse.ArgumentParser(description="Comprehensive BDT diagnostics")
    p.add_argument("--signal", required=True, help="Signal ROOT file")
    p.add_argument("--backgrounds", nargs="+", required=True,
                   help="Background ROOT files")
    p.add_argument("--data", help="Optional data ROOT file")
    p.add_argument("--treename", default="mva")
    p.add_argument("--output", default="bdt_comprehensive_analysis.pdf")
    return p.parse_args()

# ============================================================
# ROOT loading
# ============================================================

def load_process(filename, treename, label):

    print(f"Loading {label}: {filename}")

    df = ROOT.RDataFrame(treename, filename)
    cols = list(df.GetColumnNames())

    if "bdt_score" not in cols:
        raise ValueError(f"No bdt_score in {filename}")

    arrays = df.AsNumpy(cols)

    data = {c: np.asarray(arrays[c], dtype=float) for c in cols}

    if "weight" not in data:
        data["weight"] = np.ones(len(data["bdt_score"]))

    data["label"] = label
    return data

# ============================================================
# Significance scan
# ============================================================

def compute_significance_scan(signal, backgrounds):

    thresholds = np.linspace(0, 1, 200)

    results = []

    for t in thresholds:

        S = np.sum(signal["weight"][signal["bdt_score"] >= t])

        B = 0
        for b in backgrounds:
            B += np.sum(b["weight"][b["bdt_score"] >= t])

        Z = S / np.sqrt(B) if B > 0 else 0
        SB = S / B if B > 0 else 0

        results.append((t, S, B, SB, Z))

    return np.array(results)

# ============================================================
# Feature ranking
# ============================================================

def rank_features(signal, backgrounds):

    print("\nFeature–Score Correlation Ranking\n")

    combined_bkg = {}
    for key in signal:
        if key in ["bdt_score", "weight", "label"]:
            continue
        vals = []
        for b in backgrounds:
            if key in b:
                vals.append(b[key])
        if vals:
            combined_bkg[key] = np.concatenate(vals)

    ranking = []

    for feat in combined_bkg:

        if feat not in signal:
            continue

        sig_corr = np.corrcoef(signal[feat], signal["bdt_score"])[0, 1]
        bkg_corr = np.corrcoef(combined_bkg[feat],
                               np.concatenate([b["bdt_score"] for b in backgrounds]))[0, 1]

        ranking.append((feat, abs(sig_corr - bkg_corr)))

    ranking.sort(key=lambda x: x[1], reverse=True)

    for i, (f, val) in enumerate(ranking[:10], 1):
        print(f"{i:2d}. {str(f):20}  separation metric = {float(val):.4f}")


    return ranking

# ============================================================
# Main plotting
# ============================================================

def create_pdf(signal, backgrounds, data, output):

    with PdfPages(output) as pdf:

        # ====================================================
        # Page 1 — What the BDT means
        # ====================================================

        fig = plt.figure(figsize=(11, 8.5))
        plt.axis("off")

        text = """
BDT SCORE INTERPRETATION

The Boosted Decision Tree (BDT) score ranges from 0 to 1.

Low score (≈0):
    Event is background-like.

High score (≈1):
    Event is signal-like (VBF Higgs candidate).

The classifier was trained using kinematic and angular
variables including:
    • Δη_jj
    • pT balance
    • cosθ_CS
    • jet kinematics
    • dimuon variables

Interpretation strategy:
    1. Examine score distributions per process
    2. Determine optimal working point via S/√B
    3. Validate stability across folds
    4. Check mass independence
    5. Inspect high-score background composition

This document provides a complete validation.
"""
        plt.text(0.05, 0.95, text, va="top", fontsize=12)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 2 — Score distributions per process
        # ====================================================

        fig = plt.figure(figsize=(11, 8))

        for b in backgrounds:
            plt.hist(b["bdt_score"],
                     bins=60,
                     weights=b["weight"],
                     histtype="step",
                     density=True,
                     label=b["label"])

        plt.hist(signal["bdt_score"],
                 bins=60,
                 weights=signal["weight"],
                 histtype="step",
                 linewidth=2,
                 density=True,
                 label="Signal")

        plt.xlabel("BDT Score")
        plt.ylabel("Normalized Density")
        plt.title("BDT Score by Process")
        plt.legend()
        plt.grid(alpha=0.3)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 3 — Log scale
        # ====================================================

        fig = plt.figure(figsize=(11, 8))

        for b in backgrounds:
            plt.hist(b["bdt_score"],
                     bins=60,
                     weights=b["weight"],
                     histtype="step",
                     density=True,
                     label=b["label"])

        plt.hist(signal["bdt_score"],
                 bins=60,
                 weights=signal["weight"],
                 histtype="step",
                 linewidth=2,
                 density=True,
                 label="Signal")

        plt.yscale("log")
        plt.xlabel("BDT Score")
        plt.ylabel("Normalized Density (log)")
        plt.title("BDT Score by Process (Log Scale)")
        plt.legend()
        plt.grid(alpha=0.3)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 4 — Significance Scan
        # ====================================================

        results = compute_significance_scan(signal, backgrounds)

        thresholds = results[:, 0]
        S = results[:, 1]
        B = results[:, 2]
        SB = results[:, 3]
        Z = results[:, 4]

        optimal_idx = np.argmax(Z)
        optimal_cut = thresholds[optimal_idx]

        print(f"\nOptimal cut = {optimal_cut:.3f}")
        print(f"Max S/sqrt(B) = {Z[optimal_idx]:.3f}")

        fig = plt.figure(figsize=(11, 8))
        plt.plot(thresholds, Z, label="S/√B")
        plt.plot(thresholds, SB, label="S/B")
        plt.axvline(optimal_cut, linestyle="--",
                    label=f"Optimal cut = {optimal_cut:.3f}")
        plt.xlabel("BDT Threshold")
        plt.ylabel("Metric")
        plt.title("Significance Scan")
        plt.legend()
        plt.grid(alpha=0.3)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 5 — Yields after optimal cut
        # ====================================================

        fig = plt.figure(figsize=(11, 8))

        labels = []
        yields = []

        S_cut = np.sum(signal["weight"][signal["bdt_score"] >= optimal_cut])
        labels.append("Signal")
        yields.append(S_cut)

        for b in backgrounds:
            y = np.sum(b["weight"][b["bdt_score"] >= optimal_cut])
            labels.append(b["label"])
            yields.append(y)

        plt.bar(labels, yields)
        plt.yscale("log")
        plt.ylabel("Weighted Yield")
        plt.title(f"Yields after BDT > {optimal_cut:.3f}")
        plt.xticks(rotation=45)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 6 — Fold Stability
        # ====================================================

        if "fold" in signal:

            fig = plt.figure(figsize=(11, 8))

            for f in np.unique(signal["fold"]):
                mask = signal["fold"] == f
                plt.hist(signal["bdt_score"][mask],
                         bins=50,
                         histtype="step",
                         density=True,
                         label=f"Fold {int(f)}")

            plt.title("Signal Score by Training Fold")
            plt.xlabel("BDT Score")
            plt.ylabel("Normalized Density")
            plt.legend()
            plt.grid(alpha=0.3)
            pdf.savefig()
            plt.close()

        # ====================================================
        # Page 7 — Score vs Mass (if available)
        # ====================================================

        if "z" in signal:

            fig = plt.figure(figsize=(11, 8))
            plt.hist2d(signal["z"], signal["bdt_score"],
                       bins=[100, 50], cmap="viridis")
            plt.xlabel("Dimuon Mass")
            plt.ylabel("BDT Score")
            plt.title("Signal: Score vs Mass")
            plt.colorbar(label="Events")
            pdf.savefig()
            plt.close()

    print(f"\nPDF saved to {output}")

# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    signal = load_process(args.signal, args.treename, "Signal")

    backgrounds = []
    for bfile in args.backgrounds:
        label = os.path.basename(bfile).replace(".root", "")
        backgrounds.append(load_process(bfile, args.treename, label))

    data = None
    if args.data:
        data = load_process(args.data, args.treename, "Data")

    rank_features(signal, backgrounds)

    create_pdf(signal, backgrounds, data, args.output)


if __name__ == "__main__":
    main()
