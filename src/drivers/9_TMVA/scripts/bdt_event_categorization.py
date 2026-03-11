#!/usr/bin/env python3
"""
Standalone BDT sideband validation workflow.

This script performs three analysis tasks:
1) Build a Run-3 blinded Data/MC BDT-score validation plot with ratio and MC stat uncertainty.
2) Use multi-bin category optimization.
3) Produce a single professional PDF report with blinded mass control plot and category tables.

Mass sidebands default to:
- left sideband:  110 < m_mumu < 115 GeV
- right sideband: 135 < m_mumu < 150 GeV
"""

import argparse
import glob
import logging
import math
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


def configure_plot_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "savefig.bbox": "tight",
    })


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run-3 BDT sideband validation and blinded category plots")

    p.add_argument("--outputs-dir", default=DEFAULT_OUTPUTS_DIR)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--output-prefix", default="bdt_sideband_validation")
    p.add_argument("--report-name", default="",
                   help="Optional PDF file name (overrides output-prefix report name)")

    p.add_argument("--eras", nargs="*", default=None,
                   help="Optional era list (e.g. 12022 22022 12023 22023 2024)")

    p.add_argument("--treename", default="mva")
    p.add_argument("--bdt-branch", default="bdt_score")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--mass-branch", default="auto")
    p.add_argument("--extra-cut", default="",
                   help="Optional extra ROOT filter expression (e.g. 'njets>=2 && mjj>400')")

    p.add_argument("--data-era", action="append", default=[],
                   help="ERA:/path/data1.root[,/path/data2.root]")
    p.add_argument("--mc-era", action="append", default=[],
                   help="ERA:PROC:/path/mc1.root[,/path/mc2.root] (repeat)")

    p.add_argument("--score-nbins", type=int, default=40)
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

    p.add_argument("--mc-normalization", choices=["none", "sidebands"], default="sidebands",
                   help="Normalize total MC to observed data sidebands for visual comparison")

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
    if args.extra_cut:
        sel_sb = f"({sel_sb}) && ({args.extra_cut})"
        sel_full = f"({sel_full}) && ({args.extra_cut})"
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

    sig_full_scores = np.array([], dtype=np.float64)
    sig_full_mass = np.array([], dtype=np.float64)
    sig_full_weight = np.array([], dtype=np.float64)
    sig_file = find_scored_files(args.outputs_dir, era).get("signal", None)
    if sig_file and os.path.isfile(sig_file):
        try:
            sig_arr = load_arrays(
                sig_file,
                args.treename,
                [args.bdt_branch, args.weight_branch, mass_branch],
                sel_full,
            )
            sig_full_scores = sig_arr[args.bdt_branch]
            sig_full_mass = sig_arr[mass_branch]
            sig_full_weight = sanitize_weights(sig_arr[args.weight_branch], args.abs_weights)
        except Exception as exc:
            logger.warning("%s: failed loading scored signal '%s': %s", era, sig_file, exc)

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
        "sig_full_scores": sig_full_scores,
        "sig_full_mass": sig_full_mass,
        "sig_full_weight": sig_full_weight,
    }


def combine_data_side_scores(era_payloads: Dict[str, Dict[str, object]]) -> np.ndarray:
    chunks: List[np.ndarray] = []
    for payload in era_payloads.values():
        chunks.append(np.asarray(payload["data_side_scores"], dtype=np.float64))  # type: ignore[index]
    return concat_chunks(chunks)


def combine_data_side_masses(era_payloads: Dict[str, Dict[str, object]], args: argparse.Namespace) -> np.ndarray:
    chunks: List[np.ndarray] = []
    for payload in era_payloads.values():
        mass = np.asarray(payload["data_full_mass"], dtype=np.float64)  # type: ignore[index]
        if mass.size == 0:
            continue
        m = sideband_mass_mask(mass, args)
        if np.any(m):
            chunks.append(mass[m])
    return concat_chunks(chunks)


def combine_mc_score_groups(
    era_payloads: Dict[str, Dict[str, object]],
    region: str = "full",
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    proc_scores: Dict[str, List[np.ndarray]] = {}
    proc_weights: Dict[str, List[np.ndarray]] = {}

    for payload in era_payloads.values():
        if region == "side":
            for proc, scores, weights in payload["mc_side_groups"]:  # type: ignore[index]
                cproc = canonical_proc(proc)
                proc_scores.setdefault(cproc, []).append(np.asarray(scores, dtype=np.float64))
                proc_weights.setdefault(cproc, []).append(np.asarray(weights, dtype=np.float64))
        else:
            for proc, scores, _mass, weights in payload["mc_full_groups"]:  # type: ignore[index]
                cproc = canonical_proc(proc)
                proc_scores.setdefault(cproc, []).append(np.asarray(scores, dtype=np.float64))
                proc_weights.setdefault(cproc, []).append(np.asarray(weights, dtype=np.float64))

    groups: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for proc in sorted(proc_scores.keys(), key=proc_sort_key):
        scores = concat_chunks(proc_scores[proc])
        weights = concat_chunks(proc_weights[proc])
        if scores.size > 0 and np.sum(weights) > 0:
            groups.append((proc, scores, weights))
    return groups


def combine_mc_mass_groups(
    era_payloads: Dict[str, Dict[str, object]],
    args: argparse.Namespace,
    region: str = "full",
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    proc_mass: Dict[str, List[np.ndarray]] = {}
    proc_weights: Dict[str, List[np.ndarray]] = {}

    for payload in era_payloads.values():
        for proc, _scores, mass, weights in payload["mc_full_groups"]:  # type: ignore[index]
            cproc = canonical_proc(proc)
            mass = np.asarray(mass, dtype=np.float64)
            weights = np.asarray(weights, dtype=np.float64)
            if region == "side":
                m = sideband_mass_mask(mass, args)
                mass = mass[m]
                weights = weights[m]
            if mass.size == 0 or np.sum(weights) <= 0:
                continue
            proc_mass.setdefault(cproc, []).append(mass)
            proc_weights.setdefault(cproc, []).append(weights)

    groups: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for proc in sorted(proc_mass.keys(), key=proc_sort_key):
        masses = concat_chunks(proc_mass[proc])
        weights = concat_chunks(proc_weights[proc])
        if masses.size > 0 and np.sum(weights) > 0:
            groups.append((proc, masses, weights))
    return groups


def combine_signal_score_overlay(
    era_payloads: Dict[str, Dict[str, object]],
) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    sc_chunks: List[np.ndarray] = []
    w_chunks: List[np.ndarray] = []
    for payload in era_payloads.values():
        sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)  # type: ignore[index]
        ww = np.asarray(payload["sig_full_weight"], dtype=np.float64)  # type: ignore[index]
        if sc.size > 0 and np.sum(ww) > 0:
            sc_chunks.append(sc)
            w_chunks.append(ww)
    if not sc_chunks:
        return None
    return concat_chunks(sc_chunks), concat_chunks(w_chunks)


def scale_mc_groups(
    groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    scale: float,
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    if abs(scale - 1.0) < 1e-12:
        return [(p, s, w) for p, s, w in groups]
    out: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for proc, vals, weights in groups:
        out.append((proc, vals, weights * float(scale)))
    return out


def group_backgrounds_for_display(
    groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    major_labels: Optional[Sequence[str]] = None,
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    majors = set(["DY", "Top", "Diboson", "EWK"] if major_labels is None else [str(x) for x in major_labels])
    score_map: Dict[str, List[np.ndarray]] = {}
    weight_map: Dict[str, List[np.ndarray]] = {}

    for proc, scores, weights in groups:
        cproc = canonical_proc(proc)
        label = cproc if cproc in majors else "Other"
        score_map.setdefault(label, []).append(np.asarray(scores, dtype=np.float64))
        weight_map.setdefault(label, []).append(np.asarray(weights, dtype=np.float64))

    out: List[Tuple[str, np.ndarray, np.ndarray]] = []
    for label in sorted(score_map.keys(), key=proc_sort_key):
        sc = concat_chunks(score_map[label])
        ww = concat_chunks(weight_map[label])
        if sc.size > 0 and np.sum(ww) > 0:
            out.append((label, sc, ww))
    return out


def weighted_hist_and_var(
    values: np.ndarray,
    bins: np.ndarray,
    weights: Optional[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    if weights is None:
        h, _ = np.histogram(values, bins=bins)
        h = h.astype(np.float64, copy=False)
        return h, h.copy()

    h, _ = np.histogram(values, bins=bins, weights=weights)
    h2, _ = np.histogram(values, bins=bins, weights=np.square(weights))
    return h.astype(np.float64, copy=False), h2.astype(np.float64, copy=False)


def compute_group_hists(
    groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    bins: np.ndarray,
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    out_var: Dict[str, np.ndarray] = {}
    total = np.zeros(len(bins) - 1, dtype=np.float64)
    total_var = np.zeros(len(bins) - 1, dtype=np.float64)
    for proc, vals, w in groups:
        h, v = weighted_hist_and_var(vals, bins, w)
        out[proc] = h
        out_var[proc] = v
        total += h
        total_var += v
    return out, out_var, total, total_var


def compute_sideband_mc_scale(
    data_side_scores: np.ndarray,
    mc_side_groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
) -> float:
    d = float(data_side_scores.size)
    m = 0.0
    for _proc, _vals, w in mc_side_groups:
        m += float(np.sum(w))
    if d <= 0 or m <= 0:
        return 1.0
    return d / m


def make_stacked_ratio_figure(
    title: str,
    xlabel: str,
    bins: np.ndarray,
    data_values: np.ndarray,
    mc_groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    logy: bool = False,
    blind_band: Optional[Tuple[float, float]] = None,
    category_band: Optional[Tuple[float, float]] = None,
    signal_overlay: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    signal_label: str = "VBF expected",
    annotation_lines: Optional[List[str]] = None,
    show_mc_unc_band: bool = True,
    shade_category_band: bool = True,
) -> Tuple[plt.Figure, np.ndarray, np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(11.5, 9.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1]},
    )

    mc_hists, mc_vars, mc_sum, mc_var_sum = compute_group_hists(mc_groups, bins)
    centers = 0.5 * (bins[:-1] + bins[1:])

    hist_vals: List[np.ndarray] = []
    hist_weights: List[np.ndarray] = []
    colors: List[str] = []
    labels: List[str] = []
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
            alpha=0.9,
            linewidth=0.4,
        )

    mc_err = np.sqrt(np.maximum(mc_var_sum, 0.0))
    if show_mc_unc_band:
        ax_top.fill_between(
            centers,
            np.clip(mc_sum - mc_err, 0.0, None),
            mc_sum + mc_err,
            step="mid",
            color="#9e9e9e",
            alpha=0.35,
            label="MC stat. unc.",
            zorder=3,
        )

    if signal_overlay is not None:
        s_vals, s_w = signal_overlay
        h_sig, _ = weighted_hist_and_var(s_vals, bins, s_w)
        ax_top.step(
            centers,
            h_sig,
            where="mid",
            color="#1b9e77",
            linewidth=2.0,
            linestyle="--",
            label=signal_label,
            zorder=4,
        )

    data_hist, data_var = weighted_hist_and_var(data_values, bins, None)
    data_err = np.sqrt(np.maximum(data_var, 0.0))
    data_mask = data_hist > 0

    if blind_band is not None:
        blo, bhi = blind_band
        data_mask &= ~((centers >= blo) & (centers <= bhi))

    if np.any(data_mask):
        ax_top.errorbar(
            centers[data_mask],
            data_hist[data_mask],
            yerr=data_err[data_mask],
            fmt="o",
            color="black",
            markersize=4.0,
            linewidth=1.0,
            label="Data",
            zorder=5,
        )

    if blind_band is not None:
        blo, bhi = blind_band
        ax_top.axvspan(blo, bhi, color="#bdbdbd", alpha=0.18, zorder=1)
        ax_bot.axvspan(blo, bhi, color="#bdbdbd", alpha=0.18, zorder=1)

    if category_band is not None:
        clo, chi = category_band
        if shade_category_band:
            ax_top.axvspan(clo, chi, color="#9e9e9e", alpha=0.15, zorder=1)
            ax_bot.axvspan(clo, chi, color="#9e9e9e", alpha=0.15, zorder=1)
        ax_top.axvline(clo, color="#616161", linestyle="--", linewidth=1.0)
        ax_top.axvline(chi, color="#616161", linestyle="--", linewidth=1.0)

    ax_top.set_title(title)
    ax_top.set_ylabel("Events / bin")
    if logy:
        positive = mc_sum[mc_sum > 0]
        ymin = float(np.min(positive) * 0.35) if positive.size else 1e-3
        ymax = float(np.max(mc_sum + mc_err)) if mc_sum.size else 10.0
        if np.any(data_mask):
            ymax = max(ymax, float(np.max(data_hist[data_mask])))
        ax_top.set_yscale("log")
        ax_top.set_ylim(max(ymin, 1e-3), max(ymax * 3.0, 1.0))

    if annotation_lines:
        ax_top.text(
            0.012,
            0.975,
            "\n".join(annotation_lines),
            transform=ax_top.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#b0b0b0"),
        )

    ax_top.legend(loc="upper right", ncol=2)

    ratio = np.divide(data_hist, mc_sum, out=np.zeros_like(data_hist), where=mc_sum > 0)
    ratio_err = np.divide(data_err, mc_sum, out=np.zeros_like(data_err), where=mc_sum > 0)
    rel_mc = np.divide(mc_err, mc_sum, out=np.zeros_like(mc_err), where=mc_sum > 0)

    ax_bot.axhline(1.0, color="#4d4d4d", linestyle="--", linewidth=1.0)
    if show_mc_unc_band:
        ax_bot.fill_between(
            centers,
            1.0 - rel_mc,
            1.0 + rel_mc,
            step="mid",
            color="#9e9e9e",
            alpha=0.35,
            zorder=2,
        )
    ratio_mask = data_mask & (mc_sum > 0)
    if np.any(ratio_mask):
        ax_bot.errorbar(
            centers[ratio_mask],
            ratio[ratio_mask],
            yerr=ratio_err[ratio_mask],
            fmt="o",
            color="black",
            markersize=4.0,
            linewidth=1.0,
            zorder=3,
        )
    ax_bot.set_xlabel(xlabel)
    ax_bot.set_ylabel("Data/MC")
    ax_bot.set_ylim(0.4, 1.6)

    fig.tight_layout()
    return fig, data_hist, data_var, mc_hists, mc_vars, mc_sum, mc_var_sum


def add_text_page(pdf: PdfPages, title: str, lines: List[str]) -> None:
    fig, ax = plt.subplots(figsize=(11.5, 8.2))
    ax.axis("off")
    ax.set_title(title, loc="left", fontsize=16, fontweight="bold", pad=14)
    ax.text(0.02, 0.95, "\n".join(lines), va="top", ha="left", fontsize=11, family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def add_table_pages(
    pdf: PdfPages,
    title: str,
    columns: List[str],
    rows: List[List[str]],
    max_rows: int = 24,
) -> None:
    if not rows:
        add_text_page(pdf, title, ["No rows available."])
        return

    n_pages = int(math.ceil(len(rows) / float(max_rows)))
    for page in range(n_pages):
        chunk = rows[page * max_rows : (page + 1) * max_rows]
        fig, ax = plt.subplots(figsize=(14.0, 8.4))
        ax.axis("off")
        page_title = title if n_pages == 1 else f"{title} (page {page+1}/{n_pages})"
        ax.set_title(page_title, loc="left", fontsize=15, fontweight="bold", pad=14)
        table = ax.table(
            cellText=chunk,
            colLabels=columns,
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(8.5)
        table.scale(1.0, 1.25)
        for (r, _c), cell in table.get_celld().items():
            if r == 0:
                cell.set_facecolor("#e6e6e6")
                cell.set_text_props(weight="bold")
            elif r % 2 == 0:
                cell.set_facecolor("#f7f7f7")
        pdf.savefig(fig)
        plt.close(fig)


def build_bdt_bin_table(
    bins: np.ndarray,
    data_hist: np.ndarray,
    mc_hists: Dict[str, np.ndarray],
    mc_vars: Dict[str, np.ndarray],
    mc_total: np.ndarray,
    mc_total_var: np.ndarray,
) -> Tuple[List[str], List[List[str]]]:
    procs = sorted(mc_hists.keys(), key=proc_sort_key)
    columns = ["Bin", "BDT range", "Data", "MC total ± stat", "Data/MC"] + procs
    rows: List[List[str]] = []
    for i in range(len(bins) - 1):
        mc = float(mc_total[i])
        mc_err = math.sqrt(max(float(mc_total_var[i]), 0.0))
        ratio = (float(data_hist[i]) / mc) if mc > 0 else 0.0
        row = [
            str(i + 1),
            f"[{bins[i]:.3f}, {bins[i+1]:.3f})",
            f"{int(data_hist[i])}",
            f"{mc:.2f} ± {mc_err:.2f}",
            f"{ratio:.3f}",
        ]
        for proc in procs:
            row.append(f"{float(mc_hists[proc][i]):.2f}")
        rows.append(row)
    return columns, rows


def summarize_categories(
    era_payloads: Dict[str, Dict[str, object]],
    edges: np.ndarray,
    args: argparse.Namespace,
    mc_weight_scale: float = 1.0,
) -> Dict[int, Dict[str, float]]:
    ncat = len(edges) - 1
    summary: Dict[int, Dict[str, float]] = {i: {} for i in range(ncat)}

    for payload in era_payloads.values():
        d_sc = np.asarray(payload["data_full_scores"], dtype=np.float64)  # type: ignore[index]
        d_m = np.asarray(payload["data_full_mass"], dtype=np.float64)  # type: ignore[index]
        if d_sc.size > 0 and d_m.size > 0:
            d_idx = category_index(d_sc, edges)
            d_sb = sideband_mass_mask(d_m, args)
            for i in range(ncat):
                m_full = (d_idx == i)
                if np.any(m_full):
                    summary[i]["data_full_count"] = summary[i].get("data_full_count", 0.0) + float(np.count_nonzero(m_full))
                m_sb = m_full & d_sb
                if np.any(m_sb):
                    summary[i]["data_side_count"] = summary[i].get("data_side_count", 0.0) + float(np.count_nonzero(m_sb))

        for proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:  # type: ignore[index]
            cproc = canonical_proc(proc)
            idx = category_index(np.asarray(m_sc, dtype=np.float64), edges)
            mass = np.asarray(m_mass, dtype=np.float64)
            w = np.asarray(m_w, dtype=np.float64) * float(mc_weight_scale)
            sb = sideband_mass_mask(mass, args)
            for i in range(ncat):
                m_full = (idx == i)
                if np.any(m_full):
                    y_full = float(np.sum(w[m_full]))
                    e2_full = float(np.sum(np.square(w[m_full])))
                    pkey_full = f"mc_full_{cproc}"
                    pkey_full_e2 = f"mcerr2_full_{cproc}"
                    summary[i][pkey_full] = summary[i].get(pkey_full, 0.0) + y_full
                    summary[i][pkey_full_e2] = summary[i].get(pkey_full_e2, 0.0) + e2_full
                    summary[i]["mc_full_total"] = summary[i].get("mc_full_total", 0.0) + y_full
                    summary[i]["mc_full_total_err2"] = summary[i].get("mc_full_total_err2", 0.0) + e2_full

                m_sb = m_full & sb
                if np.any(m_sb):
                    y_sb = float(np.sum(w[m_sb]))
                    e2_sb = float(np.sum(np.square(w[m_sb])))
                    pkey_sb = f"mc_side_{cproc}"
                    pkey_sb_e2 = f"mcerr2_side_{cproc}"
                    summary[i][pkey_sb] = summary[i].get(pkey_sb, 0.0) + y_sb
                    summary[i][pkey_sb_e2] = summary[i].get(pkey_sb_e2, 0.0) + e2_sb
                    summary[i]["mc_side_total"] = summary[i].get("mc_side_total", 0.0) + y_sb
                    summary[i]["mc_side_total_err2"] = summary[i].get("mc_side_total_err2", 0.0) + e2_sb

    return summary


def build_category_summary_table(
    edges: np.ndarray,
    summary: Dict[int, Dict[str, float]],
) -> Tuple[List[str], List[List[str]]]:
    columns = [
        "Category",
        "BDT range",
        "Data sidebands",
        "Data full",
        "MC sidebands ± stat",
        "Data/MC (sidebands)",
        "MC full (110-150) ± stat",
    ]
    rows: List[List[str]] = []
    for i in range(len(edges) - 1):
        vals = summary.get(i, {})
        d_side = float(vals.get("data_side_count", 0.0))
        d_full = float(vals.get("data_full_count", 0.0))
        m_side = float(vals.get("mc_side_total", 0.0))
        me_side = math.sqrt(max(float(vals.get("mc_side_total_err2", 0.0)), 0.0))
        m_full = float(vals.get("mc_full_total", 0.0))
        me_full = math.sqrt(max(float(vals.get("mc_full_total_err2", 0.0)), 0.0))
        ratio = d_side / m_side if m_side > 0 else 0.0
        rows.append([
            str(i),
            f"[{edges[i]:.3f}, {edges[i+1]:.3f})",
            f"{d_side:.0f}",
            f"{d_full:.0f}",
            f"{m_side:.2f} ± {me_side:.2f}",
            f"{ratio:.3f}",
            f"{m_full:.2f} ± {me_full:.2f}",
        ])
    return columns, rows


def build_category_process_table(
    edges: np.ndarray,
    summary: Dict[int, Dict[str, float]],
) -> Tuple[List[str], List[List[str]]]:
    proc_keys = set()
    for vals in summary.values():
        for k in vals.keys():
            if k.startswith("mc_full_") and (not k.startswith("mcerr2_")):
                proc_keys.add(k)
    procs = sorted([k.replace("mc_full_", "") for k in proc_keys], key=proc_sort_key)

    columns = ["Category", "BDT range"] + procs
    rows: List[List[str]] = []
    for i in range(len(edges) - 1):
        vals = summary.get(i, {})
        row = [
            str(i),
            f"[{edges[i]:.3f}, {edges[i+1]:.3f})",
        ]
        for p in procs:
            row.append(f"{float(vals.get(f'mc_full_{p}', 0.0)):.2f}")
        rows.append(row)
    return columns, rows


def build_category_mass_grid(
    era_payloads: Dict[str, Dict[str, object]],
    edges: np.ndarray,
    args: argparse.Namespace,
    mass_bins: np.ndarray,
    mc_weight_scale: float = 1.0,
) -> plt.Figure:
    ncat = len(edges) - 1
    ncols = 2 if ncat > 1 else 1
    nrows = int(math.ceil(ncat / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 4.2 * nrows), sharex=True)
    axes_flat = np.atleast_1d(axes).ravel()

    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])

    for i in range(ncat):
        ax = axes_flat[i]

        data_mass_chunks: List[np.ndarray] = []
        bkg_mass_chunks: List[np.ndarray] = []
        bkg_weight_chunks: List[np.ndarray] = []
        sig_mass_chunks: List[np.ndarray] = []
        sig_weight_chunks: List[np.ndarray] = []

        lo, hi = float(edges[i]), float(edges[i + 1])

        for payload in era_payloads.values():
            d_sc = np.asarray(payload["data_full_scores"], dtype=np.float64)  # type: ignore[index]
            d_m = np.asarray(payload["data_full_mass"], dtype=np.float64)  # type: ignore[index]
            if d_sc.size > 0 and d_m.size > 0:
                d_cat = (d_sc >= lo) & (d_sc < hi)
                d_sb = sideband_mass_mask(d_m, args)
                dsel = d_cat & d_sb
                if np.any(dsel):
                    data_mass_chunks.append(d_m[dsel])

            for _proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:  # type: ignore[index]
                m_sc = np.asarray(m_sc, dtype=np.float64)
                m_mass = np.asarray(m_mass, dtype=np.float64)
                m_w = np.asarray(m_w, dtype=np.float64) * float(mc_weight_scale)
                msel = (m_sc >= lo) & (m_sc < hi)
                if np.any(msel):
                    bkg_mass_chunks.append(m_mass[msel])
                    bkg_weight_chunks.append(m_w[msel])

            s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)  # type: ignore[index]
            s_mass = np.asarray(payload["sig_full_mass"], dtype=np.float64)  # type: ignore[index]
            s_w = np.asarray(payload["sig_full_weight"], dtype=np.float64)  # type: ignore[index]
            if s_sc.size > 0 and s_mass.size > 0 and s_w.size > 0:
                ssel = (s_sc >= lo) & (s_sc < hi)
                if np.any(ssel):
                    sig_mass_chunks.append(s_mass[ssel])
                    sig_weight_chunks.append(s_w[ssel] * float(mc_weight_scale))

        data_mass = concat_chunks(data_mass_chunks)
        bkg_mass = concat_chunks(bkg_mass_chunks)
        bkg_w = concat_chunks(bkg_weight_chunks)
        sig_mass = concat_chunks(sig_mass_chunks)
        sig_w = concat_chunks(sig_weight_chunks)

        h_data, v_data = weighted_hist_and_var(data_mass, mass_bins, None)
        h_bkg, _ = weighted_hist_and_var(bkg_mass, mass_bins, bkg_w if bkg_w.size > 0 else None)
        h_sig, _ = weighted_hist_and_var(sig_mass, mass_bins, sig_w if sig_w.size > 0 else None)
        h_sb = h_bkg + h_sig

        # Sideband-only polynomial fit to data for the background curve.
        fit_mask = (centers < args.blind_min) | (centers > args.blind_max)
        fit_x = centers[fit_mask]
        fit_y = h_data[fit_mask]
        fit_w = 1.0 / np.sqrt(np.maximum(v_data[fit_mask], 1.0))
        if fit_x.size >= 4 and np.sum(fit_y) > 0:
            try:
                coeffs = np.polyfit(fit_x, fit_y, deg=2, w=fit_w)
                h_bkg_fit = np.polyval(coeffs, centers)
                h_bkg_fit = np.clip(h_bkg_fit, 0.0, None)
            except Exception:
                h_bkg_fit = h_bkg.copy()
        else:
            h_bkg_fit = h_bkg.copy()

        # S+B model: background fit + expected signal template
        h_sb_fit = h_bkg_fit + h_sig

        width = (mass_bins[1] - mass_bins[0]) if len(mass_bins) > 1 else 1.0
        ax.bar(
            centers,
            h_bkg,
            width=width,
            color="#4f79a7",
            alpha=0.35,
            edgecolor="#2f4f6f",
            linewidth=0.4,
            label="Bkg component" if i == 0 else None,
            zorder=1,
        )
        ax.step(
            centers,
            h_bkg_fit,
            where="mid",
            color="#1f77b4",
            linewidth=2.0,
            label="Bkg fit" if i == 0 else None,
            zorder=3,
        )
        ax.step(
            centers,
            h_sb_fit,
            where="mid",
            color="#d62728",
            linewidth=2.0,
            label="S+B fit" if i == 0 else None,
            zorder=4,
        )

        dmask = (centers < args.blind_min) | (centers > args.blind_max)
        if np.any(dmask):
            ax.errorbar(
                centers[dmask],
                h_data[dmask],
                yerr=np.sqrt(np.maximum(v_data[dmask], 1.0)),
                fmt="o",
                color="black",
                markersize=3.5,
                linewidth=1.0,
                label="Data" if i == 0 else None,
                zorder=5,
            )

        ax.axvspan(args.blind_min, args.blind_max, color="#bdbdbd", alpha=0.18, zorder=0)
        ax.set_title(f"Cat {i}: [{lo:.3f}, {hi:.3f})")
        ax.set_ylabel("Events / bin")
        ax.grid(alpha=0.25)

    for j in range(ncat, len(axes_flat)):
        axes_flat[j].axis("off")

    for ax in axes_flat[:max(ncat, 1)]:
        ax.set_xlabel("m_mumu [GeV]")

    if ncat > 0:
        axes_flat[0].legend(loc="upper right", ncol=2)

    fig.suptitle(
        "Per-category dimuon mass: Data (blinded), Bkg fit, and S+B fit",
        fontsize=15,
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.975])
    return fig


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


def main() -> None:
    args = parse_args()
    logger = setup_logger(args.log_level)
    configure_plot_style()

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

    edges = optimize_category_edges(args, eras, logger)
    interior = np.asarray(edges[1:-1], dtype=np.float64)
    category_band = None
    if interior.size >= 1:
        category_band = (float(np.min(interior)), float(np.max(interior)))

    run3_data_side_scores = combine_data_side_scores(era_payloads)
    run3_data_side_mass = combine_data_side_masses(era_payloads, args)

    run3_mc_score_side = combine_mc_score_groups(era_payloads, region="side")
    run3_mc_score_full = combine_mc_score_groups(era_payloads, region="full")
    run3_mc_mass_full = combine_mc_mass_groups(era_payloads, args, region="full")

    if run3_data_side_scores.size == 0:
        raise SystemExit("[error] No data events found in sidebands for the selected eras.")
    if not run3_mc_score_full:
        raise SystemExit("[error] No MC events found for BDT plotting.")

    mc_scale = 1.0
    if args.mc_normalization == "sidebands":
        mc_scale = compute_sideband_mc_scale(run3_data_side_scores, run3_mc_score_side)
        logger.info("MC normalization (sidebands): scale factor = %.6f", mc_scale)

    run3_mc_score_full = scale_mc_groups(run3_mc_score_full, mc_scale)
    run3_mc_mass_full = scale_mc_groups(run3_mc_mass_full, mc_scale)
    run3_mc_score_display = group_backgrounds_for_display(run3_mc_score_full)

    score_bins = np.linspace(args.score_min, args.score_max, args.score_nbins + 1)
    fig_score, data_hist, _data_var, mc_hists, mc_vars, mc_total, mc_total_var = make_stacked_ratio_figure(
        title="BDT output distribution",
        xlabel="BDT score",
        bins=score_bins,
        data_values=run3_data_side_scores,
        mc_groups=run3_mc_score_display,
        logy=args.score_logy,
        category_band=category_band,
        show_mc_unc_band=False,
        shade_category_band=False,
        annotation_lines=[
            f"Eras: {', '.join(sorted(era_payloads.keys()))}",
            f"Boundary range: {category_band[0]:.3f} to {category_band[1]:.3f}" if category_band else "Boundary range: N/A",
        ],
    )

    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    fig_mass, _m_dhist, _m_dvar, _m_hists, _m_vars, _m_tot, _m_tot_var = make_stacked_ratio_figure(
        title="Blinded dimuon mass: data sidebands with MC expectation",
        xlabel="m_mumu [GeV]",
        bins=mass_bins,
        data_values=run3_data_side_mass,
        mc_groups=run3_mc_mass_full,
        logy=False,
        blind_band=(args.blind_min, args.blind_max),
        show_mc_unc_band=False,
    )

    fig_cat_mass = build_category_mass_grid(
        era_payloads=era_payloads,
        edges=edges,
        args=args,
        mass_bins=mass_bins,
        mc_weight_scale=mc_scale,
    )

    bdt_cols, bdt_rows = build_bdt_bin_table(score_bins, data_hist, mc_hists, mc_vars, mc_total, mc_total_var)

    cat_summary = summarize_categories(era_payloads, edges, args, mc_weight_scale=mc_scale)
    cat_cols, cat_rows = build_category_summary_table(edges, cat_summary)
    cat_proc_cols, cat_proc_rows = build_category_process_table(edges, cat_summary)

    report_name = args.report_name.strip() if args.report_name else f"{args.output_prefix}_report.pdf"
    report_path = report_name if os.path.isabs(report_name) else os.path.join(args.outdir, report_name)

    with PdfPages(report_path) as pdf:
        add_text_page(
            pdf,
            "Run-3 BDT Sideband Validation Report",
            [
                f"Outputs directory: {args.outputs_dir}",
                f"Eras used: {', '.join(sorted(era_payloads.keys()))}",
                f"Extra cut: {args.extra_cut if args.extra_cut else '(none)'}",
                f"Blinding window: {args.blind_min:.0f} < m_mumu < {args.blind_max:.0f} GeV",
                f"Sidebands: ({args.sb1_min:.0f}, {args.sb1_max:.0f}) U ({args.sb2_min:.0f}, {args.sb2_max:.0f}) GeV",
                f"Optimized category edges: {', '.join([f'{x:.4f}' for x in edges])}",
            ],
        )
        pdf.savefig(fig_score)
        pdf.savefig(fig_mass)
        pdf.savefig(fig_cat_mass)
        add_table_pages(pdf, "BDT bin-by-bin comparison", bdt_cols, bdt_rows, max_rows=20)
        add_table_pages(pdf, "Optimized category summary", cat_cols, cat_rows, max_rows=20)
        add_table_pages(pdf, "Category composition by process (full mass window)", cat_proc_cols, cat_proc_rows, max_rows=20)

    plt.close(fig_score)
    plt.close(fig_mass)
    plt.close(fig_cat_mass)

    logger.info("Wrote report PDF: %s", report_path)
    print(f"[done] Wrote report: {report_path}")


if __name__ == "__main__":
    main()
