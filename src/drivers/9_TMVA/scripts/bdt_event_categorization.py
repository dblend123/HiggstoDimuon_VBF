#!/usr/bin/env python3
"""
Run-3 VBF BDT categorisation & blinded sideband validation — MAIN script.

PURPOSE
-------
This is the primary analysis script.  It ingests scored ROOT ntuples produced
by the TMVA training step, runs multibin category optimisation, performs a
blinded sideband Data/MC comparison, fits background models in each category,
and writes a 9-page PDF report.

REPORT PAGES
------------
  1.  Cover page  — paths, eras, blinding window, optimised edges
  2.  BDT score   — full 110-150 GeV window, Data vs stacked MC + ratio
  3.  BDT bin table — per-bin signal/background/sensitivity summary
  4.  Blinded mass — sideband data points overlaid on normalised MC stack
  5.  Per-category mass grid — background model fits + ratio-to-best
  6.  Fit comparison table  — chi2/ndf for all models per category
  7.  Category summary table — data/MC yields + expected significance
  8.  Process composition table — per-process MC yields per category

MASS SIDEBANDS (defaults)
--------------------------
  Left  sideband : 110 < m_µµ < 115 GeV
  Right sideband : 135 < m_µµ < 150 GeV
  Blind window   : 115 < m_µµ < 135 GeV

RELATED SCRIPTS
---------------
  bdt_analysis.py             — shared utility library (import only)
  vbf_bdt_validation.py       — single-era diagnostic PDF (KS, ROC, features)
  bdt_publishable_figures.py  — thin wrapper emitting a figure-only PDF
  significance_inference.py   — post-categorisation significance scan
"""

import argparse
from array import array
import glob
import logging
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
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
        "axes.facecolor": "#fafafa",
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.titlepad": 6,
        "axes.labelsize": 11,
        "axes.grid": True,
        "grid.alpha": 0.30,
        "grid.linestyle": ":",
        "grid.color": "#a0a0a0",
        # Ticks drawn inward on all four sides (CMS style)
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        # Legend with subtle frame
        "legend.frameon": True,
        "legend.framealpha": 0.88,
        "legend.edgecolor": "#c8c8c8",
        "legend.fontsize": 9,
        "savefig.bbox": "tight",
        "savefig.dpi": 150,
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
    p.add_argument("--multibin-likelihood-topk", type=int, default=13)

    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    return p.parse_args()


def setup_logger(level: str) -> logging.Logger:
    logger = logging.getLogger("bdt_validation")
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


def combine_data_full_scores(era_payloads: Dict[str, Dict[str, object]]) -> np.ndarray:
    chunks: List[np.ndarray] = []
    for payload in era_payloads.values():
        chunks.append(np.asarray(payload["data_full_scores"], dtype=np.float64))  # type: ignore[index]
    return concat_chunks(chunks)


def combine_data_side_scores(era_payloads: Dict[str, Dict[str, object]]) -> np.ndarray:
    """Return BDT scores for data events in the mass sidebands across all eras."""
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


def _weighted_linear_prediction(
    A_fit: np.ndarray,
    y_fit: np.ndarray,
    sigma_fit: np.ndarray,
    A_eval: np.ndarray,
) -> Optional[Tuple[np.ndarray, int]]:
    try:
        if A_fit.ndim != 2 or A_eval.ndim != 2:
            return None
        npar = int(A_fit.shape[1])
        if npar <= 0 or int(A_fit.shape[0]) <= npar:
            return None
        w = 1.0 / np.clip(sigma_fit, 1.0, None)
        Aw = A_fit * w[:, None]
        yw = y_fit * w
        beta, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
        pred = A_eval @ beta
        return pred.astype(np.float64, copy=False), npar
    except Exception:
        return None


def fit_background_models_full(
    centers: np.ndarray,
    h_data: np.ndarray,
    v_data: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, str, float, Dict[str, Tuple[float, int]], Dict[str, np.ndarray]]:
    fit_mask = np.ones_like(centers, dtype=bool)
    x_fit = centers[fit_mask]
    y_fit = h_data[fit_mask]
    sigma_fit = np.sqrt(np.maximum(v_data[fit_mask], 1.0))

    stats: Dict[str, Tuple[float, int]] = {}
    preds: Dict[str, np.ndarray] = {}
    if x_fit.size < 4 or np.sum(y_fit) <= 0:
        return h_data.copy(), "fallback_data", float("inf"), stats, preds

    x_min = float(args.mass_min)
    x_max = float(args.mass_max)
    x_span = max(x_max - x_min, 1.0)
    t_fit = np.clip((x_fit - x_min) / x_span, 0.0, 1.0)
    t_eval = np.clip((centers - x_min) / x_span, 0.0, 1.0)
    mz = 91.1876
    gamma = 2.4952
    bw_fit = 1.0 / (((x_fit - mz) ** 2) + (0.5 * gamma) ** 2)
    bw_eval = 1.0 / (((centers - mz) ** 2) + (0.5 * gamma) ** 2)

    fit_results: Dict[str, Tuple[np.ndarray, int]] = {}

    def _register_linear(name: str, basis_fit, basis_eval) -> None:
        out = _weighted_linear_prediction(
            basis_fit,
            y_fit,
            sigma_fit,
            basis_eval,
        )
        if out is not None:
            fit_results[name] = out

    _register_linear(
        "BWZ",
        np.column_stack([bw_fit, np.ones_like(x_fit), t_fit]),
        np.column_stack([bw_eval, np.ones_like(centers), t_eval]),
    )
    _register_linear(
        "BWZRedux",
        np.column_stack([bw_fit, np.ones_like(x_fit), t_fit, t_fit * t_fit]),
        np.column_stack([bw_eval, np.ones_like(centers), t_eval, t_eval * t_eval]),
    )
    _register_linear(
        "BWZGamma",
        np.column_stack([bw_fit, 1.0 / np.clip(x_fit, 1.0, None), np.ones_like(x_fit)]),
        np.column_stack([bw_eval, 1.0 / np.clip(centers, 1.0, None), np.ones_like(centers)]),
    )
    _register_linear(
        "Bernstein",
        np.column_stack([
            (1.0 - t_fit) ** 3,
            3.0 * t_fit * (1.0 - t_fit) ** 2,
            3.0 * (t_fit ** 2) * (1.0 - t_fit),
            t_fit ** 3,
        ]),
        np.column_stack([
            (1.0 - t_eval) ** 3,
            3.0 * t_eval * (1.0 - t_eval) ** 2,
            3.0 * (t_eval ** 2) * (1.0 - t_eval),
            t_eval ** 3,
        ]),
    )
    _register_linear(
        "Power Series",
        np.column_stack([np.ones_like(t_fit), t_fit, t_fit * t_fit, t_fit * t_fit * t_fit]),
        np.column_stack([np.ones_like(t_eval), t_eval, t_eval * t_eval, t_eval * t_eval * t_eval]),
    )

    y_pos = np.clip(y_fit, 1e-6, None)
    sigma_log = np.clip(sigma_fit / y_pos, 1e-3, None)
    exp_out = _weighted_linear_prediction(
        np.column_stack([
            np.ones_like(x_fit),
            np.clip((x_fit - x_min) / x_span, 0.0, 1.0),
            np.clip((x_fit - x_min) / x_span, 0.0, 1.0) ** 2,
        ]),
        np.log(y_pos),
        sigma_log,
        np.column_stack([
            np.ones_like(centers),
            np.clip((centers - x_min) / x_span, 0.0, 1.0),
            np.clip((centers - x_min) / x_span, 0.0, 1.0) ** 2,
        ]),
    )
    if exp_out is not None:
        pred_log, npar_exp = exp_out
        fit_results["Exponential Series"] = (np.exp(np.clip(pred_log, -50.0, 50.0)), npar_exp)

    best_name = "fallback_data"
    best_chi2ndf = float("inf")
    best_pred = h_data.copy()

    for name, (pred_raw, npar) in fit_results.items():
        pred = np.clip(pred_raw.astype(np.float64, copy=False), 0.0, None)
        preds[name] = pred
        resid = y_fit - pred[fit_mask]
        chi2 = float(np.sum(np.square(resid) / np.maximum(v_data[fit_mask], 1.0)))
        ndf = int(max(x_fit.size - npar, 1))
        chi2ndf = chi2 / float(ndf)
        stats[name] = (chi2, ndf)
        if np.isfinite(chi2ndf) and chi2ndf < best_chi2ndf:
            best_chi2ndf = chi2ndf
            best_name = name
            best_pred = pred

    return best_pred, best_name, best_chi2ndf, stats, preds


def fit_signal_dcb_template(
    centers: np.ndarray,
    mass_bins: np.ndarray,
    h_signal: np.ndarray,
    args: argparse.Namespace,
) -> np.ndarray:
    if h_signal.size == 0 or float(np.sum(h_signal)) <= 0.0:
        return np.zeros_like(centers, dtype=np.float64)

    try:
        edges = array("d", [float(x) for x in mass_bins])
        htmp = ROOT.TH1D("htmp_sig_dcb", "htmp_sig_dcb", len(mass_bins) - 1, edges)
        htmp.Sumw2()
        for ibin, y in enumerate(h_signal, start=1):
            htmp.SetBinContent(ibin, float(max(y, 0.0)))
            htmp.SetBinError(ibin, math.sqrt(max(float(y), 0.0)))

        m = ROOT.RooRealVar("m_sig_dcb", "m_sig_dcb", float(args.mass_min), float(args.mass_max))
        data = ROOT.RooDataHist("data_sig_dcb", "data_sig_dcb", ROOT.RooArgList(m), htmp)

        mean = ROOT.RooRealVar("mean_sig_dcb", "mean_sig_dcb", 125.0, 123.0, 127.0)
        sigma = ROOT.RooRealVar("sigma_sig_dcb", "sigma_sig_dcb", 1.6, 0.2, 6.0)
        alpha_l = ROOT.RooRealVar("alphaL_sig_dcb", "alphaL_sig_dcb", 1.5, 0.2, 6.0)
        n_l = ROOT.RooRealVar("nL_sig_dcb", "nL_sig_dcb", 5.0, 0.5, 80.0)
        alpha_r = ROOT.RooRealVar("alphaR_sig_dcb", "alphaR_sig_dcb", 2.0, 0.2, 6.0)
        n_r = ROOT.RooRealVar("nR_sig_dcb", "nR_sig_dcb", 8.0, 0.5, 80.0)
        dcb = ROOT.RooCrystalBall(
            "dcb_sig",
            "dcb_sig",
            m,
            mean,
            sigma,
            alpha_l,
            n_l,
            alpha_r,
            n_r,
        )
        dcb.fitTo(
            data,
            ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.Minimizer("Minuit2", "migrad"),
            ROOT.RooFit.Strategy(1),
        )

        norm_set = ROOT.RooArgSet(m)
        pred = np.zeros_like(centers, dtype=np.float64)
        for i, x in enumerate(centers):
            m.setVal(float(x))
            pred[i] = max(float(dcb.getVal(norm_set)), 0.0)

        total_sig = float(np.sum(h_signal))
        if np.sum(pred) > 0.0:
            pred *= total_sig / float(np.sum(pred))
        return pred
    except Exception:
        return h_signal.astype(np.float64, copy=True)


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
    show_mc_fraction_inset: bool = False,
    ratio_ylabel: str = "Data/MC",
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

    def _stack_sort_key(item: Tuple[str, np.ndarray, np.ndarray]) -> Tuple[int, float, int]:
        proc = item[0]
        cproc = canonical_proc(proc)
        yld = float(np.sum(mc_hists.get(proc, np.zeros(len(bins) - 1, dtype=np.float64))))
        # Keep DY at the bottom so sub-dominant backgrounds are visible on top.
        if cproc == "DY":
            return (0, 0.0, proc_sort_key(proc))
        return (1, -yld, proc_sort_key(proc))

    plot_groups = sorted(list(mc_groups), key=_stack_sort_key)

    # Per-process hatch patterns so small contributions (EWK, VV, VVV) remain
    # visually distinguishable even when their contribution is narrow on the
    # y-axis scale relative to DY.
    _HATCH_MAP: Dict[str, str] = {
        "EWK": "///",
        "VVV": "xxx",
        "VV": "...",
        "DYVBF": "\\\\\\",
    }
    bin_widths = bins[1:] - bins[:-1]
    base = np.zeros(len(bins) - 1, dtype=np.float64)
    for proc, _vals, _weights in plot_groups:
        hh = mc_hists.get(proc, None)
        if hh is None or np.sum(hh) <= 0:
            continue
        color = proc_color(proc)
        hatch = _HATCH_MAP.get(canonical_proc(proc), "")
        ax_top.bar(
            bins[:-1], hh, bottom=base, width=bin_widths,
            align="edge", color=color, hatch=hatch,
            edgecolor="#222222" if hatch else "#2e2e2e",
            linewidth=0.55 if not hatch else 0.8,
            alpha=0.87, label=proc, zorder=1,
        )
        base = base + hh

    if np.any(mc_sum > 0):
        ax_top.step(
            centers,
            mc_sum,
            where="mid",
            color="#2b2b2b",
            linewidth=1.3,
            label="MC total",
            zorder=3,
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

    if show_mc_fraction_inset and len(plot_groups) > 1 and np.any(mc_sum > 0):
        ax_frac = ax_top.inset_axes([0.60, 0.08, 0.36, 0.30])
        frac_base = np.zeros_like(mc_sum)
        for proc, _vals, _weights in plot_groups:
            hh = mc_hists.get(proc, None)
            if hh is None:
                continue
            frac = np.divide(100.0 * hh, mc_sum, out=np.zeros_like(hh), where=mc_sum > 0)
            frac_top = frac_base + frac
            ax_frac.fill_between(
                centers,
                frac_base,
                frac_top,
                step="mid",
                color=proc_color(proc),
                alpha=0.85,
                linewidth=0.0,
            )
            frac_base = frac_top
        ax_frac.set_ylim(0.0, 100.0)
        ax_frac.set_yticks([0.0, 50.0, 100.0])
        ax_frac.set_xticks([])
        ax_frac.set_title("MC composition [%]", fontsize=8)
        ax_frac.tick_params(axis="y", labelsize=7)
        ax_frac.grid(alpha=0.2)

    handles, leg_labels = ax_top.get_legend_handles_labels()
    if handles:
        ax_top.legend(
            handles, leg_labels,
            loc="upper right", ncol=2,
            frameon=True, framealpha=0.90,
            edgecolor="#c0c0c0", fontsize=8.5,
        )

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
    ax_bot.set_ylabel(ratio_ylabel)
    ax_bot.set_ylim(0.4, 1.6)

    fig.tight_layout()
    return fig, data_hist, data_var, mc_hists, mc_vars, mc_sum, mc_var_sum


def add_text_page(pdf: PdfPages, title: str, lines: List[str]) -> None:
    fig = plt.figure(figsize=(11.5, 8.2))
    fig.patch.set_facecolor("white")
    # Title
    fig.text(0.06, 0.93, title, ha="left", va="top",
             fontsize=15, fontweight="bold", color="#1a1a1a")
    # Separator rule under the title
    sep = plt.Line2D([0.05, 0.95], [0.895, 0.895],
                     transform=fig.transFigure, color="#4a4a4a", linewidth=1.2)
    fig.add_artist(sep)
    # Content in monospace for alignment
    fig.text(0.06, 0.87, "\n".join(lines), ha="left", va="top",
             fontsize=10.5, family="monospace", color="#2a2a2a")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def add_table_pages(
    pdf: PdfPages,
    title: str,
    columns: List[str],
    rows: List[List[str]],
    max_rows: int = 18,
) -> None:
    """Render *rows* as a styled table, paginating when necessary.

    Each cell is auto-sized to its content.  Wide tables (many columns)
    get a proportionally wider figure so text does not overflow.
    """
    if not rows:
        add_text_page(pdf, title, ["No rows available."])
        return

    ncols = len(columns)
    # Scale figure width with column count; cap at ~24 inches for readability.
    fig_w = max(13.0, min(ncols * 1.65, 24.0))
    # Reduce font to fit very wide tables without creating tiny text.
    fontsize = 8.5 if ncols <= 9 else (7.5 if ncols <= 12 else 7.0)

    n_pages = int(math.ceil(len(rows) / float(max_rows)))
    for page in range(n_pages):
        chunk = rows[page * max_rows : (page + 1) * max_rows]
        fig, ax = plt.subplots(figsize=(fig_w, 8.4))
        ax.axis("off")
        page_title = title if n_pages == 1 else f"{title}  (page {page + 1} / {n_pages})"
        ax.set_title(page_title, loc="left", fontsize=13, fontweight="bold", pad=12)
        table = ax.table(
            cellText=chunk,
            colLabels=columns,
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
        table.scale(1.0, 1.35)
        # Auto-size column widths to content so text does not overflow cells.
        table.auto_set_column_width(col=list(range(ncols)))
        for (r, c), cell in table.get_celld().items():
            if r == 0:
                # Distinct header: steel-blue background, white bold text.
                cell.set_facecolor("#2c5f8a")
                cell.set_text_props(weight="bold", color="white")
                cell.set_edgecolor("#1a3f60")
            elif r % 2 == 1:
                cell.set_facecolor("#ffffff")
                cell.set_edgecolor("#d0d0d0")
            else:
                cell.set_facecolor("#eef3f8")
                cell.set_edgecolor("#d0d0d0")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def build_bdt_signal_background_table(
    bins: np.ndarray,
    signal_hist: np.ndarray,
    bkg_hist: np.ndarray,
    bkg_var: np.ndarray,
    data_hist: np.ndarray,
) -> Tuple[List[str], List[List[str]]]:
    columns = [
        "BDT bin",
        "Total signal",
        "Bkg. ± DeltaB",
        "Data",
        "S/(S+B) (%)",
        "S/sqrt(B)",
    ]
    rows: List[List[str]] = []

    nbins = len(bins) - 1
    for i in range(nbins):
        s = float(signal_hist[i]) if i < signal_hist.size else 0.0
        b = float(bkg_hist[i]) if i < bkg_hist.size else 0.0
        db = math.sqrt(max(float(bkg_var[i]), 0.0)) if i < bkg_var.size else 0.0
        d = float(data_hist[i]) if i < data_hist.size else 0.0

        purity = (100.0 * s / (s + b)) if (s + b) > 0 else 0.0
        signif = (s / math.sqrt(b)) if b > 0 else 0.0

        rows.append([
            f"[{bins[i]:.3f}, {bins[i+1]:.3f})",
            f"{s:.3f}",
            f"{b:.2f} ± {db:.2f}",
            f"{int(round(d))}",
            f"{purity:.2f}",
            f"{signif:.3f}",
        ])

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

        s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)  # type: ignore[index]
        s_w = np.asarray(payload["sig_full_weight"], dtype=np.float64)  # type: ignore[index]
        if s_sc.size > 0 and s_w.size > 0:
            s_idx = category_index(s_sc, edges)
            s_w_scaled = s_w * float(mc_weight_scale)
            for i in range(ncat):
                m_sig = (s_idx == i)
                if np.any(m_sig):
                    s_full = float(np.sum(s_w_scaled[m_sig]))
                    s_e2_full = float(np.sum(np.square(s_w_scaled[m_sig])))
                    summary[i]["sig_full_total"] = summary[i].get("sig_full_total", 0.0) + s_full
                    summary[i]["sig_full_total_err2"] = summary[i].get("sig_full_total_err2", 0.0) + s_e2_full

    return summary


def build_category_summary_table(
    edges: np.ndarray,
    summary: Dict[int, Dict[str, float]],
) -> Tuple[List[str], List[List[str]]]:
    columns = [
        "Category",
        "BDT range",
        "Data sidebands",
        "MC sidebands ± stat",
        "Data/MC (sidebands)",
        "MC full (110-150) ± stat",
        "Expected signif. (S/sqrt(B))",
    ]
    rows: List[List[str]] = []
    for i in range(len(edges) - 1):
        vals = summary.get(i, {})
        d_side = float(vals.get("data_side_count", 0.0))
        m_side = float(vals.get("mc_side_total", 0.0))
        me_side = math.sqrt(max(float(vals.get("mc_side_total_err2", 0.0)), 0.0))
        m_full = float(vals.get("mc_full_total", 0.0))
        me_full = math.sqrt(max(float(vals.get("mc_full_total_err2", 0.0)), 0.0))
        s_full = float(vals.get("sig_full_total", 0.0))
        ratio = d_side / m_side if m_side > 0 else 0.0
        signif = (s_full / math.sqrt(m_full)) if m_full > 0.0 else 0.0
        rows.append([
            str(i),
            f"[{edges[i]:.3f}, {edges[i+1]:.3f})",
            f"{d_side:.0f}",
            f"{m_side:.2f} ± {me_side:.2f}",
            f"{ratio:.3f}",
            f"{m_full:.2f} ± {me_full:.2f}",
            f"{signif:.3f}",
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
) -> Tuple[plt.Figure, List[str], List[List[str]]]:
    ncat = len(edges) - 1
    ncols = 2 if ncat > 1 else 1
    nrows = int(math.ceil(ncat / float(ncols)))
    # Extra vertical space per row (6.2 in) + 0.5 in padding for suptitle
    fig = plt.figure(figsize=(14.0, 6.2 * nrows + 0.5))
    gs = fig.add_gridspec(
        nrows * 2,
        ncols,
        height_ratios=[3.0, 1.0] * nrows,
        hspace=0.35,
        wspace=0.28,
    )

    top_axes = []
    residual_axes = []

    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    fit_columns = [
        "Category",
        "BDT range",
        "Best model",
        "Best chi2/ndf",
        "S/sqrt(B)",
        "BWZ",
        "BWZRedux",
        "BWZGamma",
        "Bernstein",
        "Exponential Series",
        "Power Series",
    ]
    fit_rows: List[List[str]] = []

    for i in range(ncat):
        row = i // ncols
        col = i % ncols

        ax = fig.add_subplot(gs[2 * row, col])
        ax_res = fig.add_subplot(gs[2 * row + 1, col], sharex=ax)
        top_axes.append(ax)
        residual_axes.append(ax_res)

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
                dsel = d_cat
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

        h_bkg_fit, best_model, best_chi2ndf, model_stats, model_preds = fit_background_models_full(
            centers,
            h_data,
            v_data,
            args,
        )

        h_sig_dcb = fit_signal_dcb_template(centers, mass_bins, h_sig, args)

        # S+B model: best background fit + DSCB signal fit.
        h_sb_fit = h_bkg_fit + h_sig_dcb

        dmask = (centers < args.blind_min) | (centers > args.blind_max)
        signif = (float(np.sum(h_sig_dcb)) / math.sqrt(float(np.sum(h_bkg_fit)))) if float(np.sum(h_bkg_fit)) > 0.0 else 0.0

        model_styles = {
            "BWZ": ("#1f77b4", "-"),
            "BWZRedux": ("#d62728", "-"),
            "BWZGamma": ("#2ca02c", "-"),
            "Bernstein": ("#9467bd", "--"),
            "Exponential Series": ("#ff7f0e", "--"),
            "Power Series": ("#8c564b", ":"),
        }

        # Draw all background models with strong styling and white stroke to avoid overlap loss.
        for model_name in ["BWZ", "BWZRedux", "BWZGamma", "Bernstein", "Exponential Series", "Power Series"]:
            if model_name not in model_preds:
                continue
            pred = model_preds[model_name]
            color, ls = model_styles.get(model_name, ("#4d4d4d", "-"))
            lw = 2.8 if model_name == best_model else 1.9
            alpha = 1.0 if model_name == best_model else 0.92
            label = f"{model_name} ({'best' if model_name == best_model else 'alt'})"
            line = ax.step(
                centers,
                pred,
                where="mid",
                color=color,
                linewidth=lw,
                linestyle=ls,
                alpha=alpha,
                label=label,
                zorder=3 if model_name == best_model else 2,
            )[0]
            line.set_path_effects([pe.Stroke(linewidth=lw + 1.2, foreground="white"), pe.Normal()])

        ax.step(
            centers,
            h_sb_fit,
            where="mid",
            color="#111111",
            linewidth=2.1,
            linestyle="-.",
            label="Best B + DSCB(S)",
            zorder=4,
        )

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
        ax.set_title(f"Cat {i}:  BDT [{lo:.3f}, {hi:.3f})", fontsize=11, pad=4)
        ax.set_ylabel("Events / bin")
        ax.grid(alpha=0.25)
        ax.tick_params(labelbottom=False)
        # Stats box in lower-left to avoid overlapping with legend (upper-right)
        # and subplot title (top edge).
        ax.text(
            0.02,
            0.38,
            f"Best: {best_model}\n$\\chi^2$/ndf = {best_chi2ndf:.3f}\nS/$\\sqrt{{B}}$ = {signif:.3f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8.0,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85, edgecolor="#b0b0b0"),
        )

        ax.legend(
            loc="upper right", fontsize=7.0, ncol=2,
            frameon=True, framealpha=0.88, edgecolor="#c0c0c0",
        )

        model_order = ["BWZ", "BWZRedux", "BWZGamma", "Bernstein", "Exponential Series", "Power Series"]
        row = [str(i), f"[{lo:.3f}, {hi:.3f})", best_model, f"{best_chi2ndf:.3f}", f"{signif:.3f}"]
        for model_name in model_order:
            if model_name in model_stats:
                chi2, ndf = model_stats[model_name]
                row.append(f"{chi2/float(max(ndf,1)):.3f} ({int(ndf)})")
            else:
                row.append("n/a")
        fit_rows.append(row)

        ax_res.axhline(1.0, color="#4d4d4d", linestyle="--", linewidth=1.0)
        for model_name in model_order:
            if model_name not in model_preds:
                continue
            pred = model_preds[model_name]
            ratio = np.divide(pred, h_bkg_fit, out=np.ones_like(pred), where=h_bkg_fit > 0)
            color, ls = model_styles.get(model_name, ("#4d4d4d", "-"))
            lw = 1.9 if model_name == best_model else 1.4
            line = ax_res.step(
                centers,
                ratio,
                where="mid",
                color=color,
                linewidth=lw,
                linestyle=ls,
                alpha=0.95,
                zorder=3 if model_name == best_model else 2,
            )[0]
            line.set_path_effects([pe.Stroke(linewidth=lw + 1.0, foreground="white"), pe.Normal()])

        ax_res.axvspan(args.blind_min, args.blind_max, color="#bdbdbd", alpha=0.18, zorder=0)
        ax_res.set_ylabel("PDF / Best")
        ax_res.grid(alpha=0.25)
        model_ratios = []
        for model_name in model_order:
            if model_name in model_preds:
                rr = np.divide(model_preds[model_name], h_bkg_fit, out=np.ones_like(h_bkg_fit), where=h_bkg_fit > 0)
                model_ratios.append(rr)
        if model_ratios:
            all_ratio = np.concatenate(model_ratios)
            finite_ratio = all_ratio[np.isfinite(all_ratio)]
            if finite_ratio.size > 0:
                rmin = float(np.min(finite_ratio))
                rmax = float(np.max(finite_ratio))
                pad = max(0.01, 0.18 * max(abs(rmin - 1.0), abs(rmax - 1.0)))
                ax_res.set_ylim(max(0.94, rmin - pad), min(1.06, rmax + pad))
            else:
                ax_res.set_ylim(0.97, 1.03)
        else:
            ax_res.set_ylim(0.97, 1.03)

    for j in range(ncat, nrows * ncols):
        row = j // ncols
        col = j % ncols
        ax_empty_top = fig.add_subplot(gs[2 * row, col])
        ax_empty_res = fig.add_subplot(gs[2 * row + 1, col], sharex=ax_empty_top)
        ax_empty_top.axis("off")
        ax_empty_res.axis("off")

    for ax in residual_axes:
        ax.set_xlabel("$m_{\\mu\\mu}$ [GeV]", fontsize=10)

    # Each subplot already has its own legend drawn in the loop above.
    # Do NOT call top_axes[0].legend() again — it would replace that
    # subplot's legend with a differently-styled one.

    # Reserve ~7% at the top for the suptitle so it does not overlap
    # with the first row of subplots.
    fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.932])
    fig.suptitle(
        "Per-category dimuon mass: background model comparison and ratio-to-best",
        fontsize=13,
        fontweight="bold",
        y=0.968,
    )
    return fig, fit_columns, fit_rows


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
        # Fallback for low-stat / non-monotonic eras: derive robust quantile edges
        # so the report can still be produced per era.
        target_n = int(max(args.multibin_min_n, min(args.multibin_max_n, 6)))
        q = np.linspace(0.0, 1.0, target_n + 1)
        try:
            q_edges = np.quantile(sig_scores, q)
            q_edges = np.asarray(q_edges, dtype=np.float64)
            q_edges[0] = 0.0
            q_edges[-1] = 1.0
            q_edges = np.clip(q_edges, 0.0, 1.0)
            q_edges = np.unique(q_edges)
            if q_edges.size < 3:
                raise ValueError("degenerate quantile edges")
            logger.warning(
                "Category optimization found no valid set; using quantile fallback edges: %s",
                ", ".join([f"{x:.6f}" for x in q_edges]),
            )
            return q_edges
        except Exception:
            uni_edges = np.linspace(0.0, 1.0, target_n + 1, dtype=np.float64)
            logger.warning(
                "Category optimization found no valid set; using uniform fallback edges: %s",
                ", ".join([f"{x:.6f}" for x in uni_edges]),
            )
            return uni_edges

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

    run3_data_full_scores = combine_data_full_scores(era_payloads)
    run3_data_side_mass = combine_data_side_masses(era_payloads, args)

    run3_mc_score_side = combine_mc_score_groups(era_payloads, region="side")
    run3_mc_score_full = combine_mc_score_groups(era_payloads, region="full")
    run3_mc_mass_full = combine_mc_mass_groups(era_payloads, args, region="full")

    if run3_data_full_scores.size == 0:
        raise SystemExit("[error] No data events found in full mass window for the selected eras.")
    if not run3_mc_score_full:
        raise SystemExit("[error] No full-window MC events found for BDT plotting.")

    mc_scale = 1.0
    if args.mc_normalization == "sidebands":
        if not run3_mc_score_side:
            raise SystemExit("[error] No sideband MC events found for sideband normalization.")
        mc_scale = compute_sideband_mc_scale(run3_data_side_mass, run3_mc_score_side)
        logger.info("MC normalization (sidebands): scale factor = %.6f", mc_scale)

    run3_mc_score_full = scale_mc_groups(run3_mc_score_full, mc_scale)
    run3_mc_mass_full = scale_mc_groups(run3_mc_mass_full, mc_scale)
    run3_mc_score_display = group_backgrounds_for_display(run3_mc_score_full)

    score_bins = np.linspace(args.score_min, args.score_max, args.score_nbins + 1)
    fig_score, data_hist, _data_var, _score_hists, _score_vars, mc_total, mc_total_var = make_stacked_ratio_figure(
        title="BDT score in full 110-150 GeV window: Data vs stacked MC",
        xlabel="BDT score",
        bins=score_bins,
        data_values=run3_data_full_scores,
        mc_groups=run3_mc_score_display,
        logy=args.score_logy,
        category_band=category_band,
        show_mc_unc_band=False,
        shade_category_band=False,
        ratio_ylabel="Data/Bkg",
        annotation_lines=[
            f"Eras: {', '.join(sorted(era_payloads.keys()))}",
            f"BDT uses full mass window: {args.mass_min:.0f} to {args.mass_max:.0f} GeV",
            f"Boundary range: {category_band[0]:.3f} to {category_band[1]:.3f}" if category_band else "Boundary range: N/A",
        ],
    )

    signal_hist = np.zeros(len(score_bins) - 1, dtype=np.float64)
    score_signal_overlay = combine_signal_score_overlay(era_payloads)
    if score_signal_overlay is not None:
        sig_scores, sig_weights = score_signal_overlay
        scaled_sig_weights = np.asarray(sig_weights, dtype=np.float64) * float(mc_scale)
        if sig_scores.size > 0 and np.sum(scaled_sig_weights) > 0:
            signal_hist, _ = weighted_hist_and_var(
                np.asarray(sig_scores, dtype=np.float64),
                score_bins,
                scaled_sig_weights,
            )

    bdt_cols, bdt_rows = build_bdt_signal_background_table(
        score_bins,
        signal_hist,
        mc_total,
        mc_total_var,
        data_hist,
    )

    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    fig_mass, _m_dhist, _m_dvar, _m_hists, _m_vars, _m_tot, _m_tot_var = make_stacked_ratio_figure(
        title="Blinded m_mumu spectrum: sideband data points with MC expectation",
        xlabel="m_mumu [GeV]",
        bins=mass_bins,
        data_values=run3_data_side_mass,
        mc_groups=run3_mc_mass_full,
        logy=False,
        blind_band=(args.blind_min, args.blind_max),
        show_mc_unc_band=False,
    )

    fig_cat_mass, fit_cols, fit_rows = build_category_mass_grid(
        era_payloads=era_payloads,
        edges=edges,
        args=args,
        mass_bins=mass_bins,
        mc_weight_scale=mc_scale,
    )

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
        add_table_pages(pdf, "BDT bin summary: signal, background, and sensitivity", bdt_cols, bdt_rows, max_rows=20)
        pdf.savefig(fig_mass)
        pdf.savefig(fig_cat_mass)
        add_table_pages(pdf, "Per-category full-range fit comparison (chi2/ndf, ndf)", fit_cols, fit_rows, max_rows=20)
        add_table_pages(pdf, "Optimized category summary", cat_cols, cat_rows, max_rows=20)
        add_table_pages(pdf, "Category composition by process (full mass window)", cat_proc_cols, cat_proc_rows, max_rows=20)

    plt.close(fig_score)
    plt.close(fig_mass)
    plt.close(fig_cat_mass)

    logger.info("Wrote report PDF: %s", report_path)
    print(f"[done] Wrote report: {report_path}")


if __name__ == "__main__":
    main()
