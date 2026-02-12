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
    p = argparse.ArgumentParser(description="Physics-grade BDT diagnostics")
    p.add_argument("--signal", required=True)
    p.add_argument("--backgrounds", nargs="+", required=True)
    p.add_argument("--data")
    p.add_argument("--treename", default="mva")
    p.add_argument("--output", default="bdt_validation.pdf")
    return p.parse_args()


# ============================================================
# ROOT loading
# ============================================================

def load_process(filename, treename, label):

    print(f"Loading {label}: {filename}")

    df = ROOT.RDataFrame(treename, filename)
    cols = list(df.GetColumnNames())

    if "bdt_score" not in cols:
        raise ValueError(f"bdt_score missing in {filename}")

    arrays = df.AsNumpy(cols)
    data = {c: np.asarray(arrays[c], dtype=float) for c in cols}

    if "weight" not in data:
        data["weight"] = np.ones(len(data["bdt_score"]))

    data["label"] = label
    return data


# ============================================================
# Efficiency & Significance
# ============================================================

def compute_scan(signal, backgrounds):

    thresholds = np.linspace(0, 1, 200)

    results = []

    for t in thresholds:

        S = np.sum(signal["weight"][signal["bdt_score"] >= t])

        B = 0
        for b in backgrounds:
            B += np.sum(b["weight"][b["bdt_score"] >= t])

        eff_sig = S / np.sum(signal["weight"])
        eff_bkg = B / np.sum([np.sum(b["weight"]) for b in backgrounds])

        Z = S / np.sqrt(B) if B > 0 else 0

        results.append((t, S, B, eff_sig, eff_bkg, Z))

    return np.array(results)


# ============================================================
# Feature ranking (clean + robust)
# ============================================================

def rank_features(signal, backgrounds):

    print("\nTop Feature–Score Separation\n")

    ranking = []

    for feat in signal:
        if feat in ["bdt_score", "weight", "label"]:
            continue

        if feat not in signal:
            continue

        sig_corr = np.corrcoef(signal[feat], signal["bdt_score"])[0, 1]

        bkg_feat = []
        bkg_score = []

        for b in backgrounds:
            if feat in b:
                bkg_feat.append(b[feat])
                bkg_score.append(b["bdt_score"])

        if not bkg_feat:
            continue

        bkg_feat = np.concatenate(bkg_feat)
        bkg_score = np.concatenate(bkg_score)

        if len(bkg_feat) < 2:
            continue

        bkg_corr = np.corrcoef(bkg_feat, bkg_score)[0, 1]

        metric = abs(sig_corr - bkg_corr)
        ranking.append((feat, metric))

    ranking.sort(key=lambda x: x[1], reverse=True)

    for i, (feat, val) in enumerate(ranking[:10], 1):
        print(f"{i:2d}. {feat:20}  Δcorr = {val:.4f}")

    return ranking


# ============================================================
# PDF creation
# ============================================================

def create_pdf(signal, backgrounds, data, output):

    with PdfPages(output) as pdf:

        # ====================================================
        # Page 1 — Physics Overview
        # ====================================================

        fig = plt.figure(figsize=(11, 8.5))
        plt.axis("off")

        text = """
BDT VALIDATION OVERVIEW

Goal:
  Separate VBF H→μμ signal from dominant backgrounds.

Score interpretation:
  Score ~ 0 → background-like
  Score ~ 1 → signal-like

Validation strategy:
  • Score distribution per process
  • Background composition at high score
  • Signal efficiency vs background rejection
  • Significance scan (S/√B)
  • Mass sculpting check
  • Fold stability
  • Yields at chosen working point
"""

        plt.text(0.05, 0.95, text, va="top", fontsize=13)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 2 — Score by process
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
        plt.title("BDT Score Distribution by Process")
        plt.legend()
        plt.grid(alpha=0.3)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 3 — High-score background composition
        # ====================================================

        fig = plt.figure(figsize=(11, 8))

        cut = 0.9

        labels = []
        yields = []

        for b in backgrounds:
            y = np.sum(b["weight"][b["bdt_score"] >= cut])
            labels.append(b["label"])
            yields.append(y)

        plt.bar(labels, yields)
        plt.yscale("log")
        plt.ylabel("Weighted Yield")
        plt.title(f"Background Composition for BDT > {cut}")
        plt.xticks(rotation=45)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 4 — Efficiency + Significance
        # ====================================================

        scan = compute_scan(signal, backgrounds)

        thresholds = scan[:, 0]
        eff_sig = scan[:, 3]
        eff_bkg = scan[:, 4]
        Z = scan[:, 5]

        opt_idx = np.argmax(Z)
        opt_cut = thresholds[opt_idx]

        print(f"\nOptimal cut = {opt_cut:.3f}")
        print(f"Max S/√B = {Z[opt_idx]:.3f}")

        fig = plt.figure(figsize=(11, 8))
        plt.plot(thresholds, eff_sig, label="Signal efficiency")
        plt.plot(thresholds, eff_bkg, label="Background efficiency")
        plt.plot(thresholds, Z/np.max(Z), label="Normalized S/√B")
        plt.axvline(opt_cut, linestyle="--",
                    label=f"Optimal cut = {opt_cut:.3f}")
        plt.xlabel("BDT Threshold")
        plt.legend()
        plt.grid(alpha=0.3)
        pdf.savefig()
        plt.close()

        # ====================================================
        # Page 5 — Mass sculpting check
        # ====================================================

        if "z" in signal:

            fig = plt.figure(figsize=(11, 8))

            plt.hist2d(signal["z"], signal["bdt_score"],
                       bins=[100, 50], cmap="viridis")

            plt.xlabel("Dimuon Mass [GeV]")
            plt.ylabel("BDT Score")
            plt.title("Signal: Mass vs BDT Score")
            plt.colorbar(label="Events")
            pdf.savefig()
            plt.close()

        # ====================================================
        # Page 6 — Yields at optimal cut
        # ====================================================

        fig = plt.figure(figsize=(11, 8))

        labels = ["Signal"]
        yields = [
            np.sum(signal["weight"][signal["bdt_score"] >= opt_cut])
        ]

        for b in backgrounds:
            labels.append(b["label"])
            yields.append(
                np.sum(b["weight"][b["bdt_score"] >= opt_cut])
            )

        plt.bar(labels, yields)
        plt.yscale("log")
        plt.ylabel("Weighted Yield")
        plt.title(f"Yields for BDT > {opt_cut:.3f}")
        plt.xticks(rotation=45)
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
