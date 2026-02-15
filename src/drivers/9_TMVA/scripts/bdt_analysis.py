#!/usr/bin/env python3

import argparse
import glob
import os
import sys

import numpy as np
import ROOT
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Add parent directory to path so we can import tmva100_DataPreparation
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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


def infer_mass_branch(path, treename, preferred):
    if preferred != "auto":
        return preferred
    cols = set(list_branches(path, treename))
    mass_candidates = ["m_mumu", "m", "mass", "mll", "mlljj", "mllg"]
    for c in mass_candidates:
        if c in cols:
            return c
    return None


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


def compute_auc(fpr, tpr):
    if len(fpr) == 0:
        return 0.0
    return float(np.trapz(tpr, fpr))


def weighted_mean_std(values, weights):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.size == 0:
        return np.nan, np.nan
    wsum = weights.sum()
    if wsum <= 0:
        return np.nan, np.nan
    mean = float(np.sum(weights * values) / wsum)
    var = float(np.sum(weights * (values - mean) ** 2) / wsum)
    return mean, np.sqrt(max(var, 0.0))


def weighted_corr(x, y, w):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    if x.size == 0:
        return np.nan
    wsum = w.sum()
    if wsum <= 0:
        return np.nan
    mx = np.sum(w * x) / wsum
    my = np.sum(w * y) / wsum
    cov = np.sum(w * (x - mx) * (y - my)) / wsum
    vx = np.sum(w * (x - mx) ** 2) / wsum
    vy = np.sum(w * (y - my) ** 2) / wsum
    if vx <= 0 or vy <= 0:
        return np.nan
    return float(cov / np.sqrt(vx * vy))


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


def format_ratio(val):
    if not np.isfinite(val):
        return "inf"
    return f"{val:.3g}"


def parse_bins(spec):
    return [float(x) for x in spec.split(",")]


def refine_bins(base_bins, focus=None):
    bins = [float(b) for b in base_bins]
    if focus is not None and np.isfinite(focus):
        extras = [
            focus - 0.10, focus - 0.05, focus - 0.03, focus - 0.02, focus - 0.01,
            focus, focus + 0.01, focus + 0.02, focus + 0.03, focus + 0.05, focus + 0.10,
        ]
        bins.extend([x for x in extras if 0.0 <= x <= 1.0])
    bins = sorted(set([round(x, 4) for x in bins if 0.0 <= x <= 1.0]))
    if bins and bins[0] != 0.0:
        bins.insert(0, 0.0)
    if bins and bins[-1] != 1.0:
        bins.append(1.0)
    return bins


def pick_top_cuts(thresholds, ams_vals, max_cuts=3, min_sep=0.02):
    order = np.argsort(ams_vals)[::-1]
    picked = []
    for idx in order:
        thr = float(thresholds[idx])
        if not np.isfinite(thr):
            continue
        if all(abs(thr - p) >= min_sep for p in picked):
            picked.append(thr)
        if len(picked) >= max_cuts:
            break
    return sorted(picked)


def plot_table(pdf, title, columns, rows, col_align=None):
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.axis("off")
    ax.set_title(title, fontsize=14)
    table = ax.table(cellText=rows, colLabels=columns, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.4)
    header_color = "#e6e6e6"
    stripe_color = "#f7f7f7"
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor(header_color)
            cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor(stripe_color)
        if col_align and col < len(col_align):
            cell._loc = col_align[col]
    pdf.savefig(fig)
    plt.close(fig)


def plot_score_hist(pdf, title, sig, bkg, bins):
    sig_total = float(np.sum(sig["weight"]))
    bkg_total = float(np.sum(bkg["weight"]))
    scale = bkg_total / max(sig_total, 1e-12)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True)

    ax = axes[0]
    ax.hist(bkg["score"], bins=bins, weights=bkg["weight"],
        histtype="stepfilled", alpha=0.35, label=f"Background (sum={bkg_total:.2f})")
    ax.hist(sig["score"], bins=bins, weights=sig["weight"] * scale,
        histtype="step", linewidth=2.2, label=f"Signal x{scale:.1e} (sum={sig_total:.2f})")
    ax.set_xlabel("BDT score")
    ax.set_ylabel("Weighted events")
    ax.set_title("Absolute yield (signal scaled)")
    ax.set_yscale("log")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1]
    ax.hist(bkg["score"], bins=bins, weights=bkg["weight"],
        histtype="step", linewidth=2, density=True, label="Background (norm)")
    ax.hist(sig["score"], bins=bins, weights=sig["weight"],
        histtype="step", linewidth=2, density=True, label="Signal (norm)")
    ax.set_xlabel("BDT score")
    ax.set_ylabel("Density")
    ax.set_title("Unit-normalized shapes")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)

    fig.suptitle(title, fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)


def plot_roc(pdf, title, fpr, tpr, auc):
    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.plot(fpr, tpr, linewidth=2, label=f"AUC = {auc:.4f}")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    pdf.savefig(fig)
    plt.close(fig)


def plot_cdf_grid(pdf, title, fold_pairs):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, item in zip(axes, fold_pairs):
        fold = item["fold"]
        tr_scores = item["tr_scores"]
        va_scores = item["va_scores"]
        D = item["D"]
        p = item["p"]

        x1, F1 = ecdf(tr_scores)
        x2, F2 = ecdf(va_scores)
        ax.plot(x1, F1, label="Train-eval", linewidth=2, linestyle="--", color="tab:blue")
        ax.plot(x2, F2, label="Validation", linewidth=2, linestyle="-", color="tab:orange")

        qs = [0.1, 0.5, 0.9]
        qtr = [np.quantile(tr_scores, q) for q in qs if tr_scores.size > 0]
        qva = [np.quantile(va_scores, q) for q in qs if va_scores.size > 0]
        if qtr:
            ax.scatter(qtr, qs, s=18, color="tab:blue", alpha=0.7)
        if qva:
            ax.scatter(qva, qs, s=18, color="tab:orange", alpha=0.7)

        ax.set_title(f"Fold {fold}")
        ax.set_xlabel("BDT score")
        ax.set_ylabel("CDF")
        ax.grid(alpha=0.3)
        ax.text(
            0.02, 0.08,
            f"KS D = {100.0 * D:.2f}%\n"
            f"p = {p:.3g}\n"
            f"N(tr)={tr_scores.size:d}, N(va)={va_scores.size:d}",
            transform=ax.transAxes,
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
    axes[0].legend(loc="lower right")
    fig.suptitle(title, fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)


def plot_mass_sculpt(pdf, title, mass, score, weight, bdt_bins):
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(bdt_bins) - 1))
    rows = []
    for i, (lo, hi) in enumerate(zip(bdt_bins[:-1], bdt_bins[1:])):
        mask = (score >= lo) & (score < hi)
        if not np.any(mask):
            continue
        m = mass[mask]
        w = weight[mask]
        mean, rms = weighted_mean_std(m, w)
        rows.append([f"{lo:.2f}-{hi:.2f}", f"{w.sum():.2f}", f"{mean:.2f}", f"{rms:.2f}"])
        ax.hist(m, bins=40, weights=w, histtype="step",
                linewidth=2, color=colors[i], label=f"{lo:.2f} <= BDT < {hi:.2f}")
    ax.set_xlabel("Mass [GeV]")
    ax.set_ylabel("Weighted events")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    if rows:
        ax.text(
            0.02, 0.98,
            "Bin  Yield  Mean  RMS\n" + "\n".join([" ".join(r) for r in rows[:6]]),
            transform=ax.transAxes,
            fontsize=8,
            va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
    pdf.savefig(fig)
    plt.close(fig)


def plot_mass_sculpt_cuts(pdf, title, mass, score, weight, cuts):
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.plasma(np.linspace(0.2, 0.9, max(len(cuts), 1)))
    rows = []
    for i, cut in enumerate(cuts):
        mask = score >= cut
        if not np.any(mask):
            continue
        m = mass[mask]
        w = weight[mask]
        mean, rms = weighted_mean_std(m, w)
        rows.append([f">= {cut:.3f}", f"{w.sum():.2f}", f"{mean:.2f}", f"{rms:.2f}"])
        ax.hist(m, bins=40, weights=w, histtype="step",
                linewidth=2, color=colors[i], label=f"BDT >= {cut:.3f}")
    ax.set_xlabel("Mass [GeV]")
    ax.set_ylabel("Weighted events")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    if rows:
        ax.text(
            0.02, 0.98,
            "Cut  Yield  Mean  RMS\n" + "\n".join([" ".join(r) for r in rows]),
            transform=ax.transAxes,
            fontsize=9,
            va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
    pdf.savefig(fig)
    plt.close(fig)


def plot_ams(pdf, title, effs, ams_vals, best_summary):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(effs, ams_vals, linewidth=2)
    ax.scatter([best_summary["eff"]], [best_summary["ams"]], color="tab:red", zorder=3)
    ax.set_xlabel("Signal efficiency")
    ax.set_ylabel("Approx median significance (AMS)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.text(
        0.02, 0.95,
        "Best AMS summary\n"
        f"AMS = {best_summary['ams']:.3f}\n"
        f"Eff = {best_summary['eff']:.2f}\n"
        f"BDT cut = {best_summary['thr']:.3f}\n"
        f"S = {best_summary['s']:.2f}, B = {best_summary['b']:.2f}\n"
        f"S/B = {best_summary['s_over_b']:.3g}",
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    pdf.savefig(fig)
    plt.close(fig)


def plot_feature_grid(pdf, title, features, data_map, bdt_bins, mass_branch=None, bin_stats=None):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(bdt_bins) - 1))
    for ax, feat in zip(axes, features):
        vals = data_map[feat]
        score = data_map["score"]
        weight = data_map["weight"]
        corr = weighted_corr(vals, score, weight)
        for i, (lo, hi) in enumerate(zip(bdt_bins[:-1], bdt_bins[1:])):
            mask = (score >= lo) & (score < hi)
            if not np.any(mask):
                continue
            ax.hist(vals[mask], bins=40, weights=weight[mask], histtype="step",
                    linewidth=2, color=colors[i], label=f"{lo:.2f}-{hi:.2f}")
        if np.isfinite(corr):
            ax.set_title(f"{feat} (corr={corr:.2f})")
        else:
            ax.set_title(feat)
        if mass_branch and feat == mass_branch:
            ax.set_xlabel(f"{feat} [GeV]")
        else:
            ax.set_xlabel(feat)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7)
    fig.suptitle(title, fontsize=14)
    if bin_stats:
        text_lines = ["BDT bin yield (signal)" ]
        text_lines += [f"{b}: {frac:.1f}%" for b, frac in bin_stats]
        fig.text(
            0.02, 0.02, "\n".join(text_lines), fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )
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
                auc = compute_auc(fpr, tpr)
                plot_roc(pdf, f"{era}: ROC (weighted)", fpr, tpr, auc)

                # AMS vs signal efficiency
                effs = np.arange(0.01, 1.01, 0.01)
                ams_vals = []
                thr_vals = []
                s_vals = []
                b_vals = []
                sig_w_q = sanitize_weights(sig[args.weight_branch], use_abs=True)
                for eff in effs:
                    thr = weighted_quantile(sig_scores, sig_w_q, 1.0 - eff)
                    s = float(sig_w[(sig_scores >= thr)].sum())
                    b = float(bkg_w[(bkg_scores >= thr)].sum())
                    thr_vals.append(thr)
                    s_vals.append(s)
                    b_vals.append(b)
                    ams_vals.append(ams(s, b))
                best_idx = int(np.argmax(ams_vals)) if ams_vals else 0
                best_summary = {
                    "eff": float(effs[best_idx]) if effs.size else 0.0,
                    "thr": float(thr_vals[best_idx]) if thr_vals else 0.0,
                    "s": float(s_vals[best_idx]) if s_vals else 0.0,
                    "b": float(b_vals[best_idx]) if b_vals else 0.0,
                    "ams": float(ams_vals[best_idx]) if ams_vals else 0.0,
                    "s_over_b": float(s_vals[best_idx] / b_vals[best_idx]) if b_vals and b_vals[best_idx] > 0 else np.inf,
                }
                plot_ams(pdf, f"{era}: AMS vs signal efficiency", effs, ams_vals, best_summary)

                suggested_cuts = pick_top_cuts(thr_vals, ams_vals, max_cuts=3, min_sep=0.02)

                # Expected signal yield per BDT bin (mass window)
                bins_yield = refine_bins(bdt_bins, best_summary["thr"])
                s_total = float(sig_w.sum())
                b_total = float(bkg_w.sum())
                rows = []
                cum_s = 0.0
                cum_b = 0.0
                for lo, hi in reversed(list(zip(bins_yield[:-1], bins_yield[1:]))):
                    mask_s = (sig_scores >= lo) & (sig_scores < hi)
                    mask_b = (bkg_scores >= lo) & (bkg_scores < hi)
                    s_bin = float(sig_w[mask_s].sum())
                    b_bin = float(bkg_w[mask_b].sum())
                    cum_s += s_bin
                    cum_b += b_bin
                    s_frac = 100.0 * s_bin / s_total if s_total > 0 else 0.0
                    s_cum_frac = 100.0 * cum_s / s_total if s_total > 0 else 0.0
                    sb_cum = cum_s / cum_b if cum_b > 0 else np.inf
                    rows.append([
                        f"{lo:.3f}-{hi:.3f}",
                        f"{s_bin:.3f}",
                        f"{s_frac:.1f}%",
                        f"{s_cum_frac:.1f}%",
                        f"{b_bin:.3f}",
                        format_ratio(sb_cum),
                    ])
                rows.reverse()
                plot_table(
                    pdf,
                    f"{era}: Expected yield per BDT bin ({args.mass_min:.0f}-{args.mass_max:.0f} GeV)\n"
                    f"S total={s_total:.2f}, B total={b_total:.2f}, best cut={best_summary['thr']:.3f}",
                    ["BDT bin", "S (bin)", "S frac", "S>=lo", "B (bin)", "S/B (>=lo)"],
                    rows,
                    col_align=["center", "right", "right", "right", "right", "right"],
                )

                # Mass sculpting (TeV)
                if mass_branch:
                    sig_mass = sig[mass_branch]
                    bkg_mass = bkg[mass_branch]
                    plot_mass_sculpt(
                        pdf,
                        f"{era}: Signal mass sculpting (GeV)",
                        sig_mass,
                        sig_scores,
                        sig_w,
                        bdt_bins,
                    )
                    plot_mass_sculpt(
                        pdf,
                        f"{era}: Background mass sculpting (GeV)",
                        bkg_mass,
                        bkg_scores,
                        bkg_w,
                        bdt_bins,
                    )
                    if suggested_cuts:
                        plot_mass_sculpt_cuts(
                            pdf,
                            f"{era}: Signal mass sculpting for high BDT cuts",
                            sig_mass,
                            sig_scores,
                            sig_w,
                            suggested_cuts,
                        )
                        plot_mass_sculpt_cuts(
                            pdf,
                            f"{era}: Background mass sculpting for high BDT cuts",
                            bkg_mass,
                            bkg_scores,
                            bkg_w,
                            suggested_cuts,
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
                    sig_feat_map[mass_branch] = sig_feat[mass_branch]
                feature_grid = [mass_branch] + top_feats if mass_branch else top_feats
                bin_stats = []
                if bdt_bins:
                    s_total_bins = float(sig_feat_map["weight"].sum())
                    for lo, hi in zip(bdt_bins[:-1], bdt_bins[1:]):
                        mask = (sig_feat_map["score"] >= lo) & (sig_feat_map["score"] < hi)
                        s_bin = float(sig_feat_map["weight"][mask].sum())
                        frac = 100.0 * s_bin / s_total_bins if s_total_bins > 0 else 0.0
                        bin_stats.append((f"{lo:.2f}-{hi:.2f}", frac))
                plot_feature_grid(
                    pdf,
                    f"{era}: Signal distributions by BDT bin",
                    feature_grid[:4],
                    sig_feat_map,
                    bdt_bins,
                    mass_branch=mass_branch,
                    bin_stats=bin_stats,
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
                        fold_pairs.append({
                            "fold": k,
                            "tr_scores": tr,
                            "va_scores": va,
                            "D": D,
                            "p": p,
                        })

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
