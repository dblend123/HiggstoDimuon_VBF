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
    p = argparse.ArgumentParser(description="BDT validation (analysis-grade)")
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
# Significance scan
# ============================================================

def compute_scan(signal, backgrounds):

    thresholds = np.linspace(0, 1, 200)
    results = []

    total_sig = np.sum(signal["weight"])
    total_bkg = sum(np.sum(b["weight"]) for b in backgrounds)

    for t in thresholds:

        S = np.sum(signal["weight"][signal["bdt_score"] >= t])
        B = sum(np.sum(b["weight"][b["bdt_score"] >= t])
                for b in backgrounds)

        eff_sig = S / total_sig if total_sig > 0 else 0
        eff_bkg = B / total_bkg if total_bkg > 0 else 0
        Z = S / np.sqrt(B) if B > 0 else 0

        results.append((t, S, B, eff_sig, eff_bkg, Z))

    return np.array(results)


# ============================================================
# PDF creation
# ============================================================

def create_pdf(signal, backgrounds, data, output):

    with PdfPages(output) as pdf:

        # ====================================================
        # 1) Stacked MC + Signal + Data
        # ====================================================

        fig, (ax, ratio_ax) = plt.subplots(
            2, 1, figsize=(11, 10),
            gridspec_kw={"height_ratios": [3, 1]},
            sharex=True
        )

        bins = np.linspace(0, 1, 50)

        # stack backgrounds
        bkg_hist = []
        bkg_labels = []

        for b in backgrounds:
            hist, _ = np.histogram(
                b["bdt_score"], bins=bins,
                weights=b["weight"]
            )
            bkg_hist.append(hist)
            bkg_labels.append(b["label"])

        ax.hist([bins[:-1]] * len(bkg_hist),
                bins=bins,
                weights=bkg_hist,
                stacked=True,
                label=bkg_labels)

        # signal overlay
        ax.hist(signal["bdt_score"],
                bins=bins,
                weights=signal["weight"],
                histtype="step",
                linewidth=2,
                label="Signal")

        # data overlay
        if data:
            data_hist, _ = np.histogram(data["bdt_score"], bins=bins)
            centers = 0.5 * (bins[1:] + bins[:-1])
            ax.errorbar(centers, data_hist,
                        yerr=np.sqrt(data_hist),
                        fmt="o", label="Data")

            # ratio
            total_mc = np.sum(bkg_hist, axis=0)
            ratio = data_hist / total_mc
            ratio_ax.plot(centers, ratio, "o")
            ratio_ax.axhline(1, linestyle="--")
            ratio_ax.set_ylabel("Data/MC")
            ratio_ax.set_ylim(0, 2)

        ax.set_ylabel("Events")
        ax.set_title("BDT Score: Stacked MC + Data")
        ax.legend()
        ax.grid(alpha=0.3)

        pdf.savefig()
        plt.close()

        # ====================================================
        # 2) High-score background composition
        # ====================================================

        cut = 0.9
        fig = plt.figure(figsize=(11, 8))

        labels = []
        yields = []

        for b in backgrounds:
            y = np.sum(b["weight"][b["bdt_score"] >= cut])
            labels.append(b["label"])
            yields.append(y)

        plt.bar(labels, yields)
        plt.yscale("log")
        plt.xticks(rotation=45)
        plt.title(f"Background Composition (BDT > {cut})")
        plt.ylabel("Weighted Yield")
        pdf.savefig()
        plt.close()

        # ====================================================
        # 3) Efficiency + Significance
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
        plt.plot(thresholds, 1 - eff_bkg,
                 label="Background rejection")
        plt.plot(thresholds, Z/np.max(Z),
                 label="Normalized S/√B")
        plt.axvline(opt_cut, linestyle="--",
                    label=f"Optimal cut = {opt_cut:.3f}")
        plt.xlabel("BDT Threshold")
        plt.legend()
        plt.grid(alpha=0.3)
        pdf.savefig()
        plt.close()

        # ====================================================
        # 4) Background mass sculpting check
        # ====================================================

        if "z" in signal:

            fig = plt.figure(figsize=(11, 8))

            for t in [0.2, 0.5, 0.8]:
                mass_vals = []
                weights = []

                for b in backgrounds:
                    mask = b["bdt_score"] > t
                    if "z" in b:
                        mass_vals.append(b["z"][mask])
                        weights.append(b["weight"][mask])

                if mass_vals:
                    mass_vals = np.concatenate(mass_vals)
                    weights = np.concatenate(weights)

                    plt.hist(mass_vals,
                             bins=60,
                             weights=weights,
                             histtype="step",
                             label=f"BDT > {t}")

            plt.xlabel("Dimuon Mass [GeV]")
            plt.ylabel("Events")
            plt.title("Background Mass Shape vs BDT Cut")
            plt.legend()
            plt.grid(alpha=0.3)
            pdf.savefig()
            plt.close()

        # ====================================================
        # 5) Yields at optimal cut
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
        plt.xticks(rotation=45)
        plt.title(f"Yields for BDT > {opt_cut:.3f}")
        plt.ylabel("Weighted Yield")
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

    create_pdf(signal, backgrounds, data, args.output)


if __name__ == "__main__":
    main()
