#!/usr/bin/env python3

import argparse
import glob
import os

import numpy as np
import ROOT
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from tmva100_DataPreparation import VARIABLES

ROOT.gROOT.SetBatch(True)

DEFAULT_BDT_BINS = [0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 1.00]

def parse_args():
    p = argparse.ArgumentParser(description="Comprehensive BDT analysis PDF (scores from Outputs).")
    p.add_argument("--outputs-dir", default="Outputs", help="Base Outputs directory")
    p.add_argument("--eras", nargs="*", default=None, help="Limit to these eras (e.g. 12022 12023)")
    p.add_argument("--treename", default="mva")
    p.add_argument("--bdt-branch", default="bdt_score")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--mass-branch", default="m")
    p.add_argument("--mass-min", type=float, default=110)
    p.add_argument("--mass-max", type=float, default=150)
    p.add_argument("--bdt-bins", default=",".join([f"{x:.2f}" for x in DEFAULT_BDT_BINS]))
    p.add_argument("--alpha", type=float, default=1e-3, help="KS p-value threshold")
    p.add_argument("--dthr-pct", type=float, default=0.20, help="KS D(%) threshold")
    p.add_argument("--abs-weights", action="store_true", help="Use |weight| for yields/significance")
    p.add_argument("--output", default="bdt_analysis.pdf")
    return p.parse_args()


def list_branches(path, treename):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return []
    t = f.Get(treename)
    if not t:
        f.Close()
        return []
    names = [br.GetName() for br in t.GetListOfBranches()]
    f.Close()
    return names

def discover_eras(outputs_dir, requested):
    if requested:
        return requested
    patterns = [
        os.path.join(outputs_dir, "*_scored_signal.root"),
        os.path.join(outputs_dir, "*_scored_background.root"),
        os.path.join(outputs_dir, "*_scored_all.root"),
        os.path.join(outputs_dir, "*", "*_scored_signal.root"),
        os.path.join(outputs_dir, "*", "*_scored_background.root"),
        os.path.join(outputs_dir, "*", "*_scored_all.root"),
    ]
    eras = set()
    for pat in patterns:
        for p in glob.glob(pat):
            base = os.path.basename(p)
            if "_scored_" in base:
                eras.add(base.split("_scored_")[0])
    return sorted(eras)


def find_scored_files(outputs_dir, era):
    candidates = [
        os.path.join(outputs_dir, f"{era}_scored_signal.root"),
        os.path.join(outputs_dir, f"{era}_scored_background.root"),
        os.path.join(outputs_dir, f"{era}_scored_all.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_signal.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_background.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_all.root"),
    ]
    out = {"signal": None, "background": None, "all": None}
    for p in candidates:
        if os.path.isfile(p):
            if p.endswith("_scored_signal.root"):
                out["signal"] = p
            elif p.endswith("_scored_background.root"):
                out["background"] = p
            elif p.endswith("_scored_all.root"):
                out["all"] = p
    return out


def sanitize_weights(w, use_abs=False):
    w = np.asarray(w, dtype=np.float64)
    w[~np.isfinite(w)] = 0.0
    if use_abs:
        w = np.abs(w)
    return w


def load_arrays(files, treename, branches, selection):
    if isinstance(files, str):
        files = [files]
    flist = ROOT.std.vector("string")()
    for f in files:
        flist.push_back(f)
    df = ROOT.RDataFrame(treename, flist)
    if selection:
        df = df.Filter(selection)
    arrs = df.AsNumpy(branches)
    return {b: np.asarray(arrs[b], dtype=np.float64) for b in branches}


def ks_2samp(x, y):
    x = np.sort(np.asarray(x, dtype=np.float64))
    y = np.sort(np.asarray(y, dtype=np.float64))
    nx, ny = x.size, y.size
    if nx == 0 or ny == 0:
        return 0.0, 0.0
    i = j = 0
    d = 0.0
    while i < nx and j < ny:
        if x[i] <= y[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / nx - j / ny))
    en = np.sqrt(nx * ny / (nx + ny))
    lam = (en + 0.12 + 0.11 / en) * d
    p = 0.0
    for k in range(1, 101):
        p += 2.0 * ((-1.0) ** (k - 1)) * np.exp(-2.0 * (k * k) * (lam * lam))
    p = float(np.clip(p, 0.0, 1.0))
    return d, p


def ecdf(a):
    a = np.sort(np.asarray(a, dtype=np.float64))
    n = a.size
    if n == 0:
        return np.array([]), np.array([])
    return a, np.arange(1, n + 1) / n


def compute_roc(sig_scores, sig_w, bkg_scores, bkg_w):
    scores = np.concatenate([sig_scores, bkg_scores])
    labels = np.concatenate([np.ones_like(sig_scores), np.zeros_like(bkg_scores)])
    weights = np.concatenate([sig_w, bkg_w])
    order = np.argsort(scores)[::-1]
    scores = scores[order]
    labels = labels[order]
    weights = weights[order]
    tps = np.cumsum(weights * (labels == 1))
    fps = np.cumsum(weights * (labels == 0))
    tpr = tps / max(tps[-1], 1e-12)
    fpr = fps / max(fps[-1], 1e-12)
    return fpr, tpr


def weighted_quantile(values, weights, q):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.size == 0:
        return np.nan
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    csum = np.cumsum(w)
    if csum[-1] <= 0:
        return v[-1]
    target = q * csum[-1]
    return np.interp(target, csum, v)


def ams(S, B):
    if S <= 0 or B <= 0:
        return 0.0
    return float(np.sqrt(2.0 * ((S + B) * np.log(1.0 + S / max(B, 1e-12)) - S)))


def parse_bins(spec):
    return [float(x) for x in spec.split(",")]


def plot_table(pdf, title, columns, rows):
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axis("off")
    ax.set_title(title, fontsize=14)
    table = ax.table(cellText=rows, colLabels=columns, loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.3)
    pdf.savefig(fig)
    plt.close(fig)


def plot_score_hist(pdf, title, sig, bkg, bins):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.hist(bkg["score"], bins=bins, weights=bkg["weight"],
            histtype="stepfilled", alpha=0.4, label="Background")
    ax.hist(sig["score"], bins=bins, weights=sig["weight"],
            histtype="step", linewidth=2, label="Signal")
    ax.set_xlabel("BDT score")
    ax.set_ylabel("Weighted events")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    pdf.savefig(fig)
    plt.close(fig)


def plot_roc(pdf, title, fpr, tpr):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.plot(fpr, tpr, linewidth=2)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    pdf.savefig(fig)
    plt.close(fig)


def plot_cdf_grid(pdf, title, fold_pairs):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, (fold, tr_scores, va_scores) in zip(axes, fold_pairs):
        x1, F1 = ecdf(tr_scores)
        x2, F2 = ecdf(va_scores)
        ax.plot(x1, F1, label="Train-eval", linewidth=2)
        ax.plot(x2, F2, label="Validation", linewidth=2)
        ax.set_title(f"Fold {fold}")
        ax.set_xlabel("BDT score")
        ax.set_ylabel("CDF")
        ax.grid(alpha=0.3)
    axes[0].legend(loc="lower right")
    fig.suptitle(title, fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)


def plot_mass_sculpt(pdf, title, mass, score, weight, bdt_bins):
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(bdt_bins) - 1))
    for i, (lo, hi) in enumerate(zip(bdt_bins[:-1], bdt_bins[1:])):
        mask = (score >= lo) & (score < hi)
        if not np.any(mask):
            continue
        ax.hist(mass[mask], bins=40, weights=weight[mask], histtype="step",
                linewidth=2, color=colors[i], label=f"{lo:.2f} <= BDT < {hi:.2f}")
    ax.set_xlabel("m_mumu [TeV]")
    ax.set_ylabel("Weighted events")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)


def plot_ams(pdf, title, effs, ams_vals):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(effs, ams_vals, linewidth=2)
    ax.set_xlabel("Signal efficiency")
    ax.set_ylabel("Approx median significance (AMS)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    pdf.savefig(fig)
    plt.close(fig)


def plot_feature_grid(pdf, title, features, data_map, bdt_bins, mass_branch=None):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(bdt_bins) - 1))
    for ax, feat in zip(axes, features):
        vals = data_map[feat]
        score = data_map["score"]
        weight = data_map["weight"]
        for i, (lo, hi) in enumerate(zip(bdt_bins[:-1], bdt_bins[1:])):
            mask = (score >= lo) & (score < hi)
            if not np.any(mask):
                continue
            ax.hist(vals[mask], bins=40, weights=weight[mask], histtype="step",
                    linewidth=2, color=colors[i], label=f"{lo:.2f}-{hi:.2f}")
        ax.set_title(feat)
        if mass_branch and feat == mass_branch:
            ax.set_xlabel(f"{feat} [TeV]")
        else:
            ax.set_xlabel(feat)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.suptitle(title, fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)


def read_feature_importance(outputs_dir, era):
    pat = os.path.join(outputs_dir, "Training", era, f"{era}_featureimportance_*.txt")
    files = sorted(glob.glob(pat))
    if not files:
        return None, []
    path = files[0]
    items = []
    with open(path, "r") as f:
        for line in f:
            if ":" not in line:
                continue
            name, pct = line.split(":", 1)
            name = name.strip()
            pct = pct.strip().replace("%", "")
            try:
                items.append((name, float(pct)))
            except ValueError:
                continue
    items = sorted(items, key=lambda x: x[1], reverse=True)
    return path, items


def main():
    args = parse_args()
    outputs_dir = os.path.abspath(args.outputs_dir)
    eras = discover_eras(outputs_dir, args.eras)
    if not eras:
        raise SystemExit("No eras discovered. Check --outputs-dir or --eras.")

    bdt_bins = parse_bins(args.bdt_bins)

    with PdfPages(args.output) as pdf:
        for era in eras:
            files = find_scored_files(outputs_dir, era)
            sig_file = files["signal"]
            bkg_file = files["background"]
            all_file = files["all"]

            if not sig_file or not bkg_file:
                print(f"[warn] {era}: missing signal/background scored files; skipping signal-vs-background plots.")

            main_ref = sig_file or bkg_file or all_file
            if not main_ref:
                print(f"[warn] {era}: no scored files found.")
                continue

            mass_branch = infer_mass_branch(main_ref, args.treename, args.mass_branch)
            if mass_branch is None:
                print(f"[warn] {era}: could not infer mass branch; mass-based plots will be skipped.")

            branches = [args.bdt_branch, args.weight_branch]
            if mass_branch:
                branches.append(mass_branch)
            if sig_file:
                sig_branches = list(branches)
                sig_branches += ["fold"]
            else:
                sig_branches = list(branches)

            sel_mass = ""
            if mass_branch:
                sel_mass = (
                    f"({mass_branch}>={args.mass_min} && {mass_branch}<={args.mass_max} "
                    f"&& TMath::Finite({mass_branch}))"
                )
            sel_finite = f"TMath::Finite({args.bdt_branch})"
            selection = sel_finite if not sel_mass else f"({sel_finite}) && ({sel_mass})"

            if sig_file and bkg_file:
                sig = load_arrays(sig_file, args.treename, sig_branches, selection)
                bkg = load_arrays(bkg_file, args.treename, branches, selection)
                sig_scores = sig[args.bdt_branch]
                bkg_scores = bkg[args.bdt_branch]
                sig_w = sanitize_weights(sig[args.weight_branch], args.abs_weights)
                bkg_w = sanitize_weights(bkg[args.weight_branch], args.abs_weights)

                plot_score_hist(
                    pdf,
                    f"{era}: BDT scores in {args.mass_min:.0f}-{args.mass_max:.0f} GeV",
                    {"score": sig_scores, "weight": sig_w},
                    {"score": bkg_scores, "weight": bkg_w},
                    bins=np.linspace(0, 1, 60),
                )

                fpr, tpr = compute_roc(sig_scores, sig_w, bkg_scores, bkg_w)
                plot_roc(pdf, f"{era}: ROC (weighted)", fpr, tpr)

                # AMS vs signal efficiency
                effs = np.arange(0.01, 1.01, 0.01)
                ams_vals = []
                sig_w_q = sanitize_weights(sig[args.weight_branch], use_abs=True)
                for eff in effs:
                    thr = weighted_quantile(sig_scores, sig_w_q, 1.0 - eff)
                    s = float(sig_w[(sig_scores >= thr)].sum())
                    b = float(bkg_w[(bkg_scores >= thr)].sum())
                    ams_vals.append(ams(s, b))
                plot_ams(pdf, f"{era}: AMS vs signal efficiency", effs, ams_vals)

                # Expected signal yield per BDT bin (mass window)
                rows = []
                for lo, hi in zip(bdt_bins[:-1], bdt_bins[1:]):
                    mask = (sig_scores >= lo) & (sig_scores < hi)
                    rows.append([f"{lo:.2f}-{hi:.2f}", f"{sig_w[mask].sum():.3f}"])
                plot_table(
                    pdf,
                    f"{era}: Expected signal yield per BDT bin ({args.mass_min:.0f}-{args.mass_max:.0f} GeV)",
                    ["BDT bin", "Yield"],
                    rows,
                )

                # Mass sculpting (TeV)
                if mass_branch:
                    sig_mass = sig[mass_branch] / 1000.0
                    bkg_mass = bkg[mass_branch] / 1000.0
                    plot_mass_sculpt(
                        pdf,
                        f"{era}: Signal mass sculpting (TeV)",
                        sig_mass,
                        sig_scores,
                        sig_w,
                        bdt_bins,
                    )
                    plot_mass_sculpt(
                        pdf,
                        f"{era}: Background mass sculpting (TeV)",
                        bkg_mass,
                        bkg_scores,
                        bkg_w,
                        bdt_bins,
                    )

                # Top-3 feature distributions from signal, by BDT bin
                _, feat_items = read_feature_importance(outputs_dir, era)
                top_feats = [x[0] for x in feat_items[:3]] if feat_items else VARIABLES[:3]
                feat_branches = [args.bdt_branch, args.weight_branch] + top_feats
                if mass_branch and mass_branch not in feat_branches:
                    feat_branches.append(mass_branch)
                sig_feat = load_arrays(sig_file, args.treename, feat_branches, selection)
                sig_feat_map = {
                    "score": sig_feat[args.bdt_branch],
                    "weight": sanitize_weights(sig_feat[args.weight_branch], args.abs_weights),
                }
                for feat in top_feats:
                    sig_feat_map[feat] = sig_feat[feat]
                if mass_branch:
                    sig_feat_map[mass_branch] = sig_feat[mass_branch] / 1000.0
                feature_grid = [mass_branch] + top_feats if mass_branch else top_feats
                plot_feature_grid(
                    pdf,
                    f"{era}: Signal distributions by BDT bin",
                    feature_grid[:4],
                    sig_feat_map,
                    bdt_bins,
                    mass_branch=mass_branch,
                )

            # KS/CDF per fold using signal scores (if fold exists)
            if sig_file:
                base_cols = list_branches(sig_file, args.treename)
                if "fold" in base_cols:
                    sel_bdt = f"TMath::Finite({args.bdt_branch})"
                    sig_all = load_arrays(sig_file, args.treename, [args.bdt_branch, "fold"], sel_bdt)
                    scores = sig_all[args.bdt_branch]
                    folds = sig_all["fold"].astype(int)
                    rows = []
                    fold_pairs = []
                    for k in range(4):
                        va = scores[folds == k]
                        tr = scores[folds != k]
                        D, p = ks_2samp(va, tr)
                        D_pct = 100.0 * D
                        passed = (p >= args.alpha) and (D_pct <= args.dthr_pct)
                        rows.append([
                            str(k),
                            f"{D_pct:.2f}%",
                            f"{len(va):d}",
                            f"{args.alpha:.3g}",
                            f"{args.dthr_pct:.2f}%",
                            "TRUE" if passed else "FALSE",
                        ])
                        fold_pairs.append((k, tr, va))

                    plot_cdf_grid(pdf, f"{era}: CDF (train-eval vs validation)", fold_pairs)
                    plot_table(
                        pdf,
                        f"{era}: KS test summary",
                        ["Fold", "KS statistic", "N events", "alpha", "Pass threshold", "KS pass"],
                        rows,
                    )
                else:
                    print(f"[warn] {era}: fold branch missing; KS/CDF skipped.")

            # Feature importance table
            path, items = read_feature_importance(outputs_dir, era)
            if items:
                rows = [[name, f"{pct:.2f}%"] for name, pct in items[:10]]
                plot_table(
                    pdf,
                    f"{era}: Top-10 feature importance (gain)\n{os.path.basename(path)}",
                    ["Feature", "Gain %"],
                    rows,
                )

    print(f"[done] Wrote {args.output}")


if __name__ == "__main__":
    main()
