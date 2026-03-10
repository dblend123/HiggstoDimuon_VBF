#!/usr/bin/env python3
"""
Standalone BDT sideband validation workflow.

This script performs three analysis tasks:
1) Build a Run-3 sideband-only Data/MC BDT-score histogram and write a bin-count table.
2) Use multi-bin category optimization.
3) Produce per-category dimuon-mass Data/MC plots in sidebands only.

Mass sidebands default to:
- left sideband:  110 < m_mumu < 115 GeV
- right sideband: 135 < m_mumu < 150 GeV
"""

import argparse
import csv
import glob
import logging
import os
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
try:
    import ROOT  # type: ignore[import-not-found]
except Exception as exc:
    raise SystemExit(
        "[error] PyROOT is not available in this Python environment. "
        "Run this script in a ROOT-enabled environment (e.g. CMSSW/lxplus)."
    ) from exc

from bdt_analysis import (
    canonical_proc,
    find_scored_files,
    load_arrays,
    optimize_multibin_boundaries,
    parse_group_map,
    proc_color,
    proc_sort_key,
    resolve_mass_branch,
    sanitize_weights,
)

ROOT.gROOT.SetBatch(True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TMVA_DIR = os.path.dirname(SCRIPT_DIR)
DEFAULT_OUTPUTS_DIR = os.path.normpath(os.path.join(TMVA_DIR, "Outputs"))
DEFAULT_OUTDIR = os.path.normpath(os.path.join(DEFAULT_OUTPUTS_DIR, "SidebandValidation"))


def resolve_input_dir(path: str, reference_dir: str) -> str:
    """Resolve an input directory robustly for different launch directories."""
    p = os.path.expanduser(str(path).strip())
    if os.path.isabs(p):
        return os.path.normpath(p)

    cwd_candidate = os.path.normpath(os.path.abspath(p))
    script_candidate = os.path.normpath(os.path.join(reference_dir, p))
    for cand in (cwd_candidate, script_candidate):
        if os.path.isdir(cand):
            return cand
    return cwd_candidate


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run-3 BDT sideband validation and blinded category plots")

    p.add_argument("--outputs-dir", default=DEFAULT_OUTPUTS_DIR)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--output-prefix", default="run3_bdt_sideband_validation")

    p.add_argument("--eras", nargs="*", default=None,
                   help="Optional era list (e.g. 12022 22022 12023 22023 2024)")

    p.add_argument("--treename", default="mva")
    p.add_argument("--bdt-branch", default="bdt_score")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--mass-branch", default="auto")

    p.add_argument("--data-era", action="append", default=[],
                   help="ERA:/path/data1.root[,/path/data2.root]")
    p.add_argument("--mc-era", action="append", default=[],
                   help="ERA:PROC:/path/mc1.root[,/path/mc2.root] (repeat)")

    p.add_argument("--score-nbins", type=int, default=20)
    p.add_argument("--score-min", type=float, default=0.0)
    p.add_argument("--score-max", type=float, default=1.0)

    p.add_argument("--mass-min", type=float, default=110.0)
    p.add_argument("--mass-max", type=float, default=150.0)
    p.add_argument("--mass-nbins", type=int, default=40)

    p.add_argument("--sb1-min", type=float, default=110.0)
    p.add_argument("--sb1-max", type=float, default=115.0)
    p.add_argument("--sb2-min", type=float, default=135.0)
    p.add_argument("--sb2-max", type=float, default=150.0)

    p.add_argument("--blind-min", type=float, default=115.0)
    p.add_argument("--blind-max", type=float, default=135.0)

    p.add_argument("--abs-weights", dest="abs_weights", action="store_true", default=True,
                   help="Use |weight| for MC yields (default: enabled)")
    p.add_argument("--signed-weights", dest="abs_weights", action="store_false",
                   help="Disable abs(weight) transform")

    p.add_argument("--score-logy", action="store_true",
                   help="Use log scale for sideband BDT score panel")

    # Category boundaries
    p.add_argument("--category-edges", default="",
                   help="Comma-separated edges or a file path containing category edges")
    p.add_argument("--category-eras", nargs="*", default=None,
                   help="Eras to use for optimization from scored_signal/scored_background files")
    p.add_argument("--multibin-min-n", type=int, default=4)
    p.add_argument("--multibin-max-n", type=int, default=6)
    p.add_argument("--multibin-candidates", type=int, default=16)
    p.add_argument("--multibin-min-sig", type=float, default=1.0)
    p.add_argument("--multibin-min-bkg", type=float, default=5.0)
    p.add_argument("--multibin-no-monotonic", action="store_true")
    p.add_argument("--multibin-objective", choices=["asimov", "likelihood"], default="asimov")
    p.add_argument("--multibin-likelihood-topk", type=int, default=20)

    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    return p.parse_args()


def setup_logger(level: str) -> logging.Logger:
    logger = logging.getLogger("bdt_sideband_validation")
    logger.setLevel(getattr(logging, level.upper()))
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.handlers.clear()
    logger.addHandler(handler)
    return logger


def parse_csv_paths(text: str) -> List[str]:
    return [x.strip() for x in str(text).split(",") if x.strip()]


def infer_eras_from_outputs(outputs_dir: str, requested: Optional[Sequence[str]]) -> List[str]:
    if requested:
        return sorted(set([str(x) for x in requested if str(x).strip()]))

    eras = set()
    patterns = [
        os.path.join(outputs_dir, "*_data_scored.root"),
        os.path.join(outputs_dir, "*_scored_all.root"),
        os.path.join(outputs_dir, "*", "*_data_scored.root"),
        os.path.join(outputs_dir, "*", "*_scored_all.root"),
        os.path.join(outputs_dir, "*", "*_scored_signal.root"),
    ]
    for pat in patterns:
        for path in glob.glob(pat):
            base = os.path.basename(path)
            if base.endswith("_data_scored.root"):
                eras.add(base[: -len("_data_scored.root")])
            elif "_scored_" in base:
                eras.add(base.split("_scored_", 1)[0])

    return sorted(eras)


def discover_datamc_for_era(outputs_dir: str, era: str) -> Tuple[List[str], List[Tuple[str, List[str]]]]:
    data_files = []
    for p in [
        os.path.join(outputs_dir, f"{era}_data_scored.root"),
        os.path.join(outputs_dir, era, f"{era}_data_scored.root"),
    ]:
        if os.path.isfile(p):
            data_files.append(p)

    proc_files: Dict[str, List[str]] = {}
    for pat in [
        os.path.join(outputs_dir, f"{era}_scored_*.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_*.root"),
    ]:
        for path in glob.glob(pat):
            base = os.path.basename(path)
            suffix = base.replace(f"{era}_scored_", "").replace(".root", "")
            if suffix.lower() in {"signal", "background", "all", "test", "data"}:
                continue
            proc_files.setdefault(suffix, []).append(path)

    mc_entries = []
    if proc_files:
        for proc, files in proc_files.items():
            uniq = sorted(set(files))
            if uniq:
                mc_entries.append((proc, uniq))
    else:
        fallback = [
            os.path.join(outputs_dir, f"{era}_scored_all.root"),
            os.path.join(outputs_dir, era, f"{era}_scored_all.root"),
            os.path.join(outputs_dir, f"{era}_scored_background.root"),
            os.path.join(outputs_dir, era, f"{era}_scored_background.root"),
        ]
        found = [p for p in fallback if os.path.isfile(p)]
        if found:
            mc_entries.append(("Background", found))

    mc_entries = sorted(mc_entries, key=lambda x: proc_sort_key(x[0]))
    return data_files, mc_entries


def parse_user_inputs(args: argparse.Namespace, logger: logging.Logger) -> Tuple[Dict[str, List[str]], Dict[str, List[Tuple[str, List[str]]]]]:
    data_map: Dict[str, List[str]] = {}
    mc_map: Dict[str, List[Tuple[str, List[str]]]] = {}

    logger.info("Using outputs directory: %s", args.outputs_dir)
    if args.eras:
        logger.info("Requested eras: %s", ", ".join(sorted(set(args.eras))))

    if args.data_era:
        raw_data = parse_group_map(args.data_era, "--data-era")
        for era, values in raw_data.items():
            files: List[str] = []
            for v in values:
                files.extend(parse_csv_paths(v))
            files = [f for f in files if os.path.isfile(f)]
            if files:
                data_map[era] = files

    if args.mc_era:
        raw_mc = parse_group_map(args.mc_era, "--mc-era")
        for era, values in raw_mc.items():
            groups: List[Tuple[str, List[str]]] = []
            for entry in values:
                if ":" not in entry:
                    raise SystemExit(f"[error] --mc-era needs ERA:PROC:files, got {era}:{entry}")
                proc, files_txt = entry.split(":", 1)
                files = [f for f in parse_csv_paths(files_txt) if os.path.isfile(f)]
                if files:
                    groups.append((proc, files))
            if groups:
                mc_map[era] = groups

    if data_map or mc_map:
        return data_map, mc_map

    eras = infer_eras_from_outputs(args.outputs_dir, args.eras)
    for era in eras:
        dfiles, mcentries = discover_datamc_for_era(args.outputs_dir, era)
        if dfiles:
            data_map[era] = dfiles
        if mcentries:
            mc_map[era] = mcentries
        if not dfiles and not mcentries:
            logger.warning("%s: no Data/MC inputs found under %s", era, args.outputs_dir)

    logger.info("Discovered eras: %s", ", ".join(sorted(set(list(data_map.keys()) + list(mc_map.keys())))))
    return data_map, mc_map


def sideband_mass_mask(mass: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    return (
        ((mass > args.sb1_min) & (mass < args.sb1_max))
        | ((mass > args.sb2_min) & (mass < args.sb2_max))
    )


def category_index(scores: np.ndarray, edges: np.ndarray) -> np.ndarray:
    idx = np.digitize(scores, edges[1:-1], right=False)
    return np.clip(idx, 0, len(edges) - 2)


def make_selections(mass_branch: str, bdt_branch: str, args: argparse.Namespace) -> Tuple[str, str]:
    mass_sb = (
        f"(({mass_branch}>{args.sb1_min} && {mass_branch}<{args.sb1_max})"
        f" || ({mass_branch}>{args.sb2_min} && {mass_branch}<{args.sb2_max}))"
    )
    mass_full = f"({mass_branch}>={args.mass_min} && {mass_branch}<={args.mass_max})"

    finite = f"TMath::Finite({bdt_branch}) && TMath::Finite({mass_branch})"
    sel_sb = f"({finite}) && ({mass_sb})"
    sel_full = f"({finite}) && ({mass_full})"
    return sel_sb, sel_full


def resolve_mass_branch_for_era(
    data_files: List[str],
    mc_entries: List[Tuple[str, List[str]]],
    treename: str,
    preferred: str,
) -> Optional[str]:
    probes = []
    probes.extend(data_files)
    for _, files in mc_entries:
        probes.extend(files)

    for p in probes:
        mb = resolve_mass_branch(p, treename, preferred)
        if mb is not None:
            return mb
    return None


def concat_chunks(chunks: List[np.ndarray]) -> np.ndarray:
    if not chunks:
        return np.array([], dtype=np.float64)
    return np.concatenate(chunks).astype(np.float64, copy=False)


def load_era_payload(
    era: str,
    data_files: List[str],
    mc_entries: List[Tuple[str, List[str]]],
    args: argparse.Namespace,
    logger: logging.Logger,
) -> Optional[Dict[str, object]]:
    mass_branch = resolve_mass_branch_for_era(data_files, mc_entries, args.treename, args.mass_branch)
    if mass_branch is None:
        logger.warning("%s: could not resolve mass branch, skipping era", era)
        return None

    sel_sb, sel_full = make_selections(mass_branch, args.bdt_branch, args)

    data_side_scores = np.array([], dtype=np.float64)
    data_full_scores = np.array([], dtype=np.float64)
    data_full_mass = np.array([], dtype=np.float64)

    if data_files:
        try:
            arr_side = load_arrays(data_files, args.treename, [args.bdt_branch], sel_sb)
            data_side_scores = arr_side[args.bdt_branch]

            arr_full = load_arrays(data_files, args.treename, [args.bdt_branch, mass_branch], sel_full)
            data_full_scores = arr_full[args.bdt_branch]
            data_full_mass = arr_full[mass_branch]
        except Exception as exc:
            logger.warning("%s: failed loading data arrays: %s", era, exc)

    mc_side_groups: List[Tuple[str, np.ndarray, np.ndarray]] = []
    mc_full_groups: List[Tuple[str, np.ndarray, np.ndarray, np.ndarray]] = []

    for proc, files in mc_entries:
        try:
            arr_side = load_arrays(files, args.treename, [args.bdt_branch, args.weight_branch], sel_sb)
            s_side = arr_side[args.bdt_branch]
            w_side = sanitize_weights(arr_side[args.weight_branch], args.abs_weights)
            if s_side.size > 0 and np.sum(w_side) > 0:
                mc_side_groups.append((proc, s_side, w_side))

            arr_full = load_arrays(files, args.treename, [args.bdt_branch, args.weight_branch, mass_branch], sel_full)
            s_full = arr_full[args.bdt_branch]
            m_full = arr_full[mass_branch]
            w_full = sanitize_weights(arr_full[args.weight_branch], args.abs_weights)
            if s_full.size > 0 and np.sum(w_full) > 0:
                mc_full_groups.append((proc, s_full, m_full, w_full))

        except Exception as exc:
            logger.warning("%s: failed loading MC '%s': %s", era, proc, exc)

    if not mc_side_groups and data_side_scores.size == 0:
        logger.warning("%s: empty sideband payload", era)
        return None

    logger.info(
        "%s: data(side)=%d, mc_procs(side)=%d, mass_branch=%s",
        era,
        data_side_scores.size,
        len(mc_side_groups),
        mass_branch,
    )

    return {
        "era": era,
        "mass_branch": mass_branch,
        "data_side_scores": data_side_scores,
        "data_full_scores": data_full_scores,
        "data_full_mass": data_full_mass,
        "mc_side_groups": mc_side_groups,
        "mc_full_groups": mc_full_groups,
    }


def combine_sideband_inputs(
    era_payloads: Dict[str, Dict[str, object]]
) -> Tuple[np.ndarray, List[Tuple[str, np.ndarray, np.ndarray]]]:
    data_chunks: List[np.ndarray] = []
    proc_scores: Dict[str, List[np.ndarray]] = {}
    proc_weights: Dict[str, List[np.ndarray]] = {}

    for payload in era_payloads.values():
        data_chunks.append(payload["data_side_scores"])  # type: ignore[index]
        for proc, scores, weights in payload["mc_side_groups"]:  # type: ignore[index]
            cproc = canonical_proc(proc)
            proc_scores.setdefault(cproc, []).append(scores)
            proc_weights.setdefault(cproc, []).append(weights)

    data_scores = concat_chunks(data_chunks)
    groups: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for proc in sorted(proc_scores.keys(), key=proc_sort_key):
        s = concat_chunks(proc_scores[proc])
        w = concat_chunks(proc_weights[proc])
        if s.size > 0 and np.sum(w) > 0:
            groups.append((proc, s, w))

    return data_scores, groups


def compute_group_hists(
    groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    bins: np.ndarray,
) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    total = np.zeros(len(bins) - 1, dtype=np.float64)
    for proc, vals, w in groups:
        h, _ = np.histogram(vals, bins=bins, weights=w)
        h = h.astype(np.float64, copy=False)
        out[proc] = h
        total += h
    return out, total


def make_stacked_ratio_figure(
    title: str,
    xlabel: str,
    bins: np.ndarray,
    data_values: np.ndarray,
    mc_groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    logy: bool = False,
    blind_band: Optional[Tuple[float, float]] = None,
) -> Tuple[plt.Figure, np.ndarray, Dict[str, np.ndarray], np.ndarray]:
    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(11, 9),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    mc_hists, mc_sum = compute_group_hists(mc_groups, bins)

    hist_vals = []
    hist_weights = []
    colors = []
    labels = []
    for proc, vals, weights in mc_groups:
        h = mc_hists.get(proc, None)
        if h is None or np.sum(h) <= 0:
            continue
        hist_vals.append(vals)
        hist_weights.append(weights)
        colors.append(proc_color(proc))
        labels.append(proc)

    if hist_vals:
        ax_top.hist(
            hist_vals,
            bins=bins,
            weights=hist_weights,
            stacked=True,
            color=colors,
            label=labels,
            histtype="stepfilled",
            alpha=0.85,
        )

    data_hist, _ = np.histogram(data_values, bins=bins)
    data_hist = data_hist.astype(np.float64, copy=False)
    centers = 0.5 * (bins[:-1] + bins[1:])
    yerr = np.sqrt(np.maximum(data_hist, 1.0))

    if data_values.size > 0:
        ax_top.errorbar(centers, data_hist, yerr=yerr, fmt="o", color="black", label="Data")

    ax_top.set_title(title)
    ax_top.set_ylabel("Events / bin")
    if logy:
        ax_top.set_yscale("log")
    ax_top.grid(alpha=0.3)
    ax_top.legend(fontsize=8, ncol=2)

    ratio = np.divide(data_hist, mc_sum, out=np.zeros_like(data_hist), where=mc_sum > 0)
    ratio_err = np.divide(yerr, mc_sum, out=np.zeros_like(yerr), where=mc_sum > 0)
    ax_bot.axhline(1.0, color="gray", linestyle="--", linewidth=1.0)
    if data_values.size > 0:
        ax_bot.errorbar(centers, ratio, yerr=ratio_err, fmt="o", color="black")
    ax_bot.set_xlabel(xlabel)
    ax_bot.set_ylabel("Data/MC")
    ax_bot.set_ylim(0.4, 1.6)
    ax_bot.grid(alpha=0.3)

    if blind_band is not None:
        lo, hi = blind_band
        ax_top.axvspan(lo, hi, color="gray", alpha=0.15)
        ax_bot.axvspan(lo, hi, color="gray", alpha=0.15)

    fig.tight_layout()
    return fig, data_hist, mc_hists, mc_sum


def write_score_bin_table(
    out_csv: str,
    bins: np.ndarray,
    data_hist: np.ndarray,
    mc_hists: Dict[str, np.ndarray],
    mc_total: np.ndarray,
) -> None:
    procs = sorted(mc_hists.keys(), key=proc_sort_key)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["bin_low", "bin_high", "data_count", "mc_total_weight", "data_over_mc"] + procs)
        for i in range(len(bins) - 1):
            ratio = float(data_hist[i] / mc_total[i]) if mc_total[i] > 0 else 0.0
            row = [
                f"{bins[i]:.6f}",
                f"{bins[i+1]:.6f}",
                f"{int(round(data_hist[i]))}",
                f"{mc_total[i]:.6f}",
                f"{ratio:.6f}",
            ]
            for proc in procs:
                row.append(f"{mc_hists[proc][i]:.6f}")
            writer.writerow(row)


def parse_category_edges_arg(text: str) -> Optional[np.ndarray]:
    s = str(text).strip()
    if not s:
        return None

    edges: List[float] = []
    if os.path.isfile(s):
        with open(s, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if not lines:
            return None

        header = lines[0].lower()
        if "bin_low" in header and "bin_high" in header:
            lows = []
            highs = []
            for ln in lines[1:]:
                toks = [t.strip() for t in ln.split(",")]
                if len(toks) < 2:
                    continue
                try:
                    lows.append(float(toks[0]))
                    highs.append(float(toks[1]))
                except ValueError:
                    continue
            edges = sorted(set(lows + highs))
        else:
            for ln in lines:
                for tok in ln.replace(",", " ").split():
                    try:
                        edges.append(float(tok))
                    except ValueError:
                        continue
    else:
        for tok in s.split(","):
            tok = tok.strip()
            if not tok:
                continue
            edges.append(float(tok))

    if not edges:
        return None

    out = np.asarray(sorted(set(edges)), dtype=np.float64)
    if out.size < 2:
        return None
    return out


def optimize_category_edges(
    args: argparse.Namespace,
    fallback_eras: Sequence[str],
    logger: logging.Logger,
) -> np.ndarray:
    manual = parse_category_edges_arg(args.category_edges)
    if manual is not None:
        logger.info("Using user-provided category edges: %s", ", ".join([f"{x:.6f}" for x in manual]))
        return manual

    opt_eras = args.category_eras if args.category_eras else list(fallback_eras)
    if not opt_eras:
        opt_eras = infer_eras_from_outputs(args.outputs_dir, None)

    sig_scores_chunks: List[np.ndarray] = []
    sig_weight_chunks: List[np.ndarray] = []
    sig_mass_chunks: List[np.ndarray] = []

    bkg_scores_chunks: List[np.ndarray] = []
    bkg_weight_chunks: List[np.ndarray] = []
    bkg_mass_chunks: List[np.ndarray] = []

    used_eras: List[str] = []
    for era in opt_eras:
        files = find_scored_files(args.outputs_dir, era)
        sig_file = files.get("signal", None)
        bkg_file = files.get("background", None)
        if not sig_file or not bkg_file:
            continue

        mass_branch = resolve_mass_branch(sig_file, args.treename, args.mass_branch)
        if mass_branch is None:
            mass_branch = resolve_mass_branch(bkg_file, args.treename, args.mass_branch)
        if mass_branch is None:
            logger.warning("%s: could not resolve mass branch for optimization; skipping", era)
            continue

        sel = (
            f"TMath::Finite({args.bdt_branch}) && TMath::Finite({mass_branch})"
            f" && ({mass_branch}>={args.mass_min} && {mass_branch}<={args.mass_max})"
        )

        try:
            s_arr = load_arrays(sig_file, args.treename, [args.bdt_branch, args.weight_branch, mass_branch], sel)
            b_arr = load_arrays(bkg_file, args.treename, [args.bdt_branch, args.weight_branch, mass_branch], sel)
        except Exception as exc:
            logger.warning("%s: failed loading optimization arrays: %s", era, exc)
            continue

        s_sc = s_arr[args.bdt_branch]
        s_w = sanitize_weights(s_arr[args.weight_branch], args.abs_weights)
        s_m = s_arr[mass_branch]

        b_sc = b_arr[args.bdt_branch]
        b_w = sanitize_weights(b_arr[args.weight_branch], args.abs_weights)
        b_m = b_arr[mass_branch]

        if s_sc.size == 0 or b_sc.size == 0:
            continue

        sig_scores_chunks.append(s_sc)
        sig_weight_chunks.append(s_w)
        sig_mass_chunks.append(s_m)

        bkg_scores_chunks.append(b_sc)
        bkg_weight_chunks.append(b_w)
        bkg_mass_chunks.append(b_m)
        used_eras.append(era)

    if not used_eras:
        raise SystemExit(
            "[error] No scored signal/background inputs found for category optimization. "
            "Provide --category-edges or make sure *_scored_signal.root and *_scored_background.root exist."
        )

    sig_scores = concat_chunks(sig_scores_chunks)
    sig_weights = concat_chunks(sig_weight_chunks)
    sig_mass = concat_chunks(sig_mass_chunks)

    bkg_scores = concat_chunks(bkg_scores_chunks)
    bkg_weights = concat_chunks(bkg_weight_chunks)
    bkg_mass = concat_chunks(bkg_mass_chunks)

    logger.info(
        "Category optimization eras: %s | S events=%d, B events=%d",
        ", ".join(used_eras),
        sig_scores.size,
        bkg_scores.size,
    )

    scan = optimize_multibin_boundaries(
        sig_scores,
        sig_weights,
        bkg_scores,
        bkg_weights,
        min_n=args.multibin_min_n,
        max_n=args.multibin_max_n,
        n_candidates=args.multibin_candidates,
        min_sig=args.multibin_min_sig,
        min_bkg=args.multibin_min_bkg,
        enforce_monotonic=(not args.multibin_no_monotonic),
        objective=args.multibin_objective,
        sig_mass=sig_mass if args.multibin_objective == "likelihood" else None,
        bkg_mass=bkg_mass if args.multibin_objective == "likelihood" else None,
        mass_min=args.mass_min,
        mass_max=args.mass_max,
        mass_nbins=args.mass_nbins,
        likelihood_topk=args.multibin_likelihood_topk,
    )

    best = scan.get("best", None)
    if best is None:
        raise SystemExit("[error] Category optimization failed to find a valid boundary set.")

    edges = np.asarray(best["edges"], dtype=np.float64)
    metric_label = best.get("metric_label", "ZA(total)")
    metric_val = float(best.get("opt_metric", best.get("za_total", 0.0)))

    logger.info(
        "Best category boundaries (n=%d): %s | %s=%.4f",
        len(edges) - 1,
        ", ".join([f"{x:.6f}" for x in edges]),
        metric_label,
        metric_val,
    )
    return edges


def write_category_edges_csv(path: str, edges: np.ndarray) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["category", "bdt_low", "bdt_high"])
        for i in range(len(edges) - 1):
            writer.writerow([i, f"{edges[i]:.6f}", f"{edges[i+1]:.6f}"])


def build_category_mass_groups(
    era_payloads: Dict[str, Dict[str, object]],
    edges: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[List[np.ndarray], List[List[Tuple[str, np.ndarray, np.ndarray]]], Dict[int, Dict[str, float]]]:
    ncat = len(edges) - 1
    data_by_cat: List[List[np.ndarray]] = [[] for _ in range(ncat)]

    proc_mass_by_cat: List[Dict[str, List[np.ndarray]]] = [dict() for _ in range(ncat)]
    proc_w_by_cat: List[Dict[str, List[np.ndarray]]] = [dict() for _ in range(ncat)]

    # Yield summary per category
    summary: Dict[int, Dict[str, float]] = {i: {} for i in range(ncat)}

    for payload in era_payloads.values():
        d_sc = payload["data_full_scores"]  # type: ignore[index]
        d_m = payload["data_full_mass"]  # type: ignore[index]
        if d_sc.size > 0 and d_m.size > 0:
            d_idx = category_index(d_sc, edges)
            d_sb = sideband_mass_mask(d_m, args)
            for i in range(ncat):
                mask = (d_idx == i) & d_sb
                if np.any(mask):
                    data_by_cat[i].append(d_m[mask])
                    summary[i]["data_count"] = summary[i].get("data_count", 0.0) + float(np.count_nonzero(mask))

        for proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:  # type: ignore[index]
            cproc = canonical_proc(proc)
            idx = category_index(m_sc, edges)
            sb = sideband_mass_mask(m_mass, args)
            for i in range(ncat):
                mask = (idx == i) & sb
                if not np.any(mask):
                    continue
                proc_mass_by_cat[i].setdefault(cproc, []).append(m_mass[mask])
                proc_w_by_cat[i].setdefault(cproc, []).append(m_w[mask])
                key = f"mc_{cproc}"
                summary[i][key] = summary[i].get(key, 0.0) + float(np.sum(m_w[mask]))

    out_data = [concat_chunks(chs) for chs in data_by_cat]
    out_mc: List[List[Tuple[str, np.ndarray, np.ndarray]]] = []
    for i in range(ncat):
        groups: List[Tuple[str, np.ndarray, np.ndarray]] = []
        for proc in sorted(proc_mass_by_cat[i].keys(), key=proc_sort_key):
            masses = concat_chunks(proc_mass_by_cat[i][proc])
            weights = concat_chunks(proc_w_by_cat[i][proc])
            if masses.size > 0 and np.sum(weights) > 0:
                groups.append((proc, masses, weights))
        out_mc.append(groups)

        mc_total = 0.0
        for _, _, w in groups:
            mc_total += float(np.sum(w))
        summary[i]["mc_total"] = mc_total

    return out_data, out_mc, summary


def write_category_yield_table(
    out_csv: str,
    edges: np.ndarray,
    summary: Dict[int, Dict[str, float]],
) -> None:
    all_proc_keys = set()
    for vals in summary.values():
        for k in vals.keys():
            if k.startswith("mc_"):
                all_proc_keys.add(k)

    proc_cols = sorted(all_proc_keys, key=lambda x: proc_sort_key(x.replace("mc_", "")))

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "category",
            "bdt_low",
            "bdt_high",
            "data_sideband_count",
            "mc_sideband_weight",
            "data_over_mc",
        ] + proc_cols)

        for i in range(len(edges) - 1):
            vals = summary.get(i, {})
            d = float(vals.get("data_count", 0.0))
            m = float(vals.get("mc_total", 0.0))
            ratio = d / m if m > 0 else 0.0
            row = [
                i,
                f"{edges[i]:.6f}",
                f"{edges[i+1]:.6f}",
                f"{d:.0f}",
                f"{m:.6f}",
                f"{ratio:.6f}",
            ]
            for col in proc_cols:
                row.append(f"{float(vals.get(col, 0.0)):.6f}")
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    logger = setup_logger(args.log_level)

    args.outputs_dir = resolve_input_dir(args.outputs_dir, TMVA_DIR)
    if not os.path.isabs(args.outdir):
        args.outdir = os.path.normpath(os.path.abspath(args.outdir))

    os.makedirs(args.outdir, exist_ok=True)

    data_map, mc_map = parse_user_inputs(args, logger)
    eras = sorted(set(list(data_map.keys()) + list(mc_map.keys())))
    if not eras:
        raise SystemExit(
            "[error] No eras found from user input or auto-discovery. "
            f"Checked outputs-dir: {args.outputs_dir}"
        )

    era_payloads: Dict[str, Dict[str, object]] = {}
    for era in eras:
        payload = load_era_payload(
            era,
            data_map.get(era, []),
            mc_map.get(era, []),
            args,
            logger,
        )
        if payload is not None:
            era_payloads[era] = payload

    if not era_payloads:
        raise SystemExit("[error] No valid era payloads were loaded.")

    score_bins = np.linspace(args.score_min, args.score_max, args.score_nbins + 1)
    run3_data_scores, run3_mc_groups = combine_sideband_inputs(era_payloads)

    if run3_data_scores.size == 0 and not run3_mc_groups:
        raise SystemExit("[error] Run-3 sideband inputs are empty.")

    score_title = (
        "Run-3 BDT score validation in mass sidebands "
        f"[{args.sb1_min:.0f},{args.sb1_max:.0f}) U [{args.sb2_min:.0f},{args.sb2_max:.0f}] GeV"
    )
    fig_score, data_hist, mc_hists, mc_total = make_stacked_ratio_figure(
        title=score_title,
        xlabel="BDT score",
        bins=score_bins,
        data_values=run3_data_scores,
        mc_groups=run3_mc_groups,
        logy=args.score_logy,
    )

    score_pdf = os.path.join(args.outdir, f"{args.output_prefix}_score_sidebands.pdf")
    score_png = os.path.join(args.outdir, f"{args.output_prefix}_score_sidebands.png")
    fig_score.savefig(score_pdf)
    fig_score.savefig(score_png, dpi=160)
    plt.close(fig_score)

    score_csv = os.path.join(args.outdir, f"{args.output_prefix}_score_sideband_bin_table.csv")
    write_score_bin_table(score_csv, score_bins, data_hist, mc_hists, mc_total)

    edges = optimize_category_edges(args, eras, logger)
    edges_csv = os.path.join(args.outdir, f"{args.output_prefix}_category_edges.csv")
    write_category_edges_csv(edges_csv, edges)

    cat_data, cat_mc, cat_summary = build_category_mass_groups(era_payloads, edges, args)
    cat_yield_csv = os.path.join(args.outdir, f"{args.output_prefix}_category_sideband_yields.csv")
    write_category_yield_table(cat_yield_csv, edges, cat_summary)

    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    mass_pdf = os.path.join(args.outdir, f"{args.output_prefix}_mass_sidebands_by_category.pdf")

    with PdfPages(mass_pdf) as pdf:
        for i in range(len(edges) - 1):
            title = (
                f"Category {i}: BDT in [{edges[i]:.3f}, {edges[i+1]:.3f}) | "
                f"blinded [{args.blind_min:.0f}, {args.blind_max:.0f}] GeV"
            )
            fig_cat, _, _, _ = make_stacked_ratio_figure(
                title=title,
                xlabel="m_mumu [GeV]",
                bins=mass_bins,
                data_values=cat_data[i],
                mc_groups=cat_mc[i],
                logy=False,
                blind_band=(args.blind_min, args.blind_max),
            )
            pdf.savefig(fig_cat)
            plt.close(fig_cat)

    logger.info("Wrote sideband score plots: %s | %s", score_pdf, score_png)
    logger.info("Wrote score bin-count table: %s", score_csv)
    logger.info("Wrote category edges: %s", edges_csv)
    logger.info("Wrote category sideband yields: %s", cat_yield_csv)
    logger.info("Wrote blinded category mass plots: %s", mass_pdf)
    print(f"[done] Outputs written to {args.outdir}")


if __name__ == "__main__":
    main()
