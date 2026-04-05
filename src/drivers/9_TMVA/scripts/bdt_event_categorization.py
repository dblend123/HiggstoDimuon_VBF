#!/usr/bin/env python3
"""
Run-3 VBF BDT categorisation & blinded sideband validation — MAIN script.

PURPOSE
-------
This is the primary analysis script.  It ingests scored ROOT ntuples produced
by the TMVA training step, runs multibin category optimisation, performs a
blinded sideband Data/MC comparison, fits background models in each category,
and writes a multi-page PDF report.

REPORT PAGES
------------
  1.  Cover page        — paths, eras, blinding window, optimised edges
  2.  BDT score         — full 110-150 GeV window, Data vs stacked MC + ratio
  3.  BDT bin table     — per-bin signal/background/sensitivity summary
  4.  Significance scan — evolution of combined Z as categories are added (AN Fig. 33)
  5.  Atanh BDT + S/√B  — transformed score with category lines (AN Fig. 34)
  6.  Blinded mass      — sideband data points overlaid on normalised MC stack
  7.  Per-category mass — background model fits + ratio-to-best (AN Fig. 36)
  8.  Fit comparison    — chi2/ndf for all models per category
  9.  VBF DCB grid      — signal DCB extraction per category (AN Fig. 38)
  10. DCB parameter table
  11. Category summary  — data/MC yields + expected significance
  12. Process composition — per-process MC yields per category

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
from scipy.optimize import minimize as _sp_minimize
from scipy.stats import chi2 as _sp_chi2, moyal as _sp_moyal
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
DEFAULT_OUTDIR = os.path.normpath(os.path.join(DEFAULT_OUTPUTS_DIR, "Validation"))


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
    p.add_argument("--output-prefix", default="bdt_validation")
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
    p.add_argument("--bkg-exclude", default="ggh",
                   help="Comma-separated canonical process names to drop from background MC "
                        "(matched via canonical_proc; default: 'ggh' to exclude ggH->mumu). "
                        "Pass '' to disable.")

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
    p.add_argument("--multibin-improvement", type=float, default=0.01,
                   help="Minimum relative gain in combined significance to add a new category (default 0.01 = 1%%)"
                   )

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
    _MIN_ROOT_BYTES = 1024  # files smaller than this are stub / corrupt outputs

    def _valid_size(path: str) -> bool:
        try:
            return os.path.getsize(path) >= _MIN_ROOT_BYTES
        except OSError:
            return False

    data_files = []
    for p in [
        os.path.join(outputs_dir, f"{era}_data_scored.root"),
        os.path.join(outputs_dir, era, f"{era}_data_scored.root"),
    ]:
        if os.path.isfile(p) and _valid_size(p):
            data_files.append(p)

    proc_files: Dict[str, List[str]] = {}
    for pat in [
        os.path.join(outputs_dir, f"{era}_scored_*.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_*.root"),
    ]:
        for path in glob.glob(pat):
            if not _valid_size(path):
                continue
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
        found = [p for p in fallback if os.path.isfile(p) and _valid_size(p)]
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

    # Build set of canonical process names to exclude from the background.
    bkg_exclude: set = set()
    raw_exclude = getattr(args, "bkg_exclude", "ggh") or ""
    for _tok in raw_exclude.split(","):
        _tok = _tok.strip().lower()
        if _tok:
            bkg_exclude.add(_tok)

    for proc, files in mc_entries:
        cproc_check = canonical_proc(proc).lower()
        if cproc_check in bkg_exclude:
            logger.warning(
                "%s: skipping MC process '%s' (canonical='%s') — excluded from "
                "background by --bkg-exclude ('%s'). "
                "This process likely peaks near 125 GeV and must not be treated as background.",
                era, proc, canonical_proc(proc), raw_exclude,
            )
            continue
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

    proc_names = ", ".join(canonical_proc(p) for p, _, _, _ in mc_full_groups)
    logger.info(
        "%s: data(side)=%d, mc_procs(side)=%d [%s], mass_branch=%s",
        era,
        data_side_scores.size,
        len(mc_side_groups),
        proc_names,
        mass_branch,
    )

    # Per-era signal and per-process sideband yield diagnostic
    sig_w_sum = float(np.sum(sig_full_weight)) if sig_full_weight.size > 0 else 0.0
    mc_side_w_sum = sum(float(np.sum(w)) for _p, _v, w in mc_side_groups)
    logger.info(
        "%s: sum(sig_full_weight)=%.3f events, sum(MC_sideband_weights)=%.1f, "
        "N_data_sideband=%d",
        era, sig_w_sum, mc_side_w_sum, int(data_side_scores.size),
    )
    for sp, _sv, sw in mc_side_groups:
        logger.info("%s:   sideband MC proc %-10s: %.1f events (raw)", era, sp, float(np.sum(sw)))

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


def _bwz_core(m: np.ndarray) -> np.ndarray:
    """Unnormalized Breit-Wigner Z-peak shape: 1 / [(m - mZ)^2 + (ΓZ/2)^2]."""
    _MZ, _GZ2 = 91.1876, (2.4952 * 0.5) ** 2
    return 1.0 / ((m - _MZ) ** 2 + _GZ2)


def _bernstein_matrix(t: np.ndarray, n: int) -> np.ndarray:
    """Return (len(t), n+1) Bernstein polynomial design matrix.  t must be in [0, 1]."""
    from math import comb as _comb
    return np.column_stack([
        float(_comb(n, k)) * (t ** k) * ((1.0 - t) ** (n - k))
        for k in range(n + 1)
    ])


def fit_background_models_full(
    centers: np.ndarray,
    h_data: np.ndarray,
    v_data: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, str, float, Dict[str, Tuple[float, int]], Dict[str, np.ndarray]]:
    """Fit CMS background model families to the dimuon mass histogram.

    Physics-inspired models (from the CMS H->mumu analysis, equations 7-9, 14):
      BWZRedux       exp(a1*t + a2*t^2) x BW(m)                        [Eq. 9]
      BWZGamma       exp(a1*t) x BW(m) x (1 + f/m^2)                   [Eq. 8]
      BWZxBern_n     exp(a*t) x BW(m) x Bern_n(t), LRT over n=1,2,3   [Eq. 7]
      LandauxBern_n  Landau(m;mg,sg) x Bern_n(t),  LRT over n=0,1,2,3 [Eq. 14]

    Agnostic models (equations 11-13):
      Bernstein_n    sum_i ci*Bi,n(t),            LRT over n=3,4,5,6   [Eq. 11]
      S-PowerLaw     sum_i ai*m^bi,               LRT over 1->2 terms  [Eq. 12]
      S-Exponential  sum_i ai*exp(bi*m),           LRT over 1->2 terms  [Eq. 13]

    Fitting: binned Poisson NLL via scipy.optimize.minimize (L-BFGS-B).
    Normalization is profiled analytically at each optimizer step.
    DOF selection: CMS likelihood ratio test at 5%: 2*DeltalogL ~ chi2(1).
    Best model: lowest chi2/ndf across all LRT-selected models.
    """
    data = h_data.astype(np.float64, copy=True)
    n_data = float(np.sum(data))
    x_min = float(args.mass_min)
    x_max = float(args.mass_max)
    x_span = max(x_max - x_min, 1.0)
    t = np.clip((centers - x_min) / x_span, 0.0, 1.0)
    bw = _bwz_core(centers)

    stats: Dict[str, Tuple[float, int]] = {}
    preds: Dict[str, np.ndarray] = {}

    if data.size < 6 or n_data <= 0.0:
        return data.copy(), "fallback_data", float("inf"), stats, preds

    rng = np.random.default_rng(seed=42)

    def _fit_nll(f_shape, p0, bounds, n_restarts: int = 2):
        """Minimize Poisson NLL with analytically-profiled normalization.

        Returns (ok, logL, pred, npar) where npar = len(p0) + 1 (shape + N).
        """
        p0_arr = np.asarray(p0, dtype=np.float64)
        best_fun = np.inf
        best_x = None

        for attempt in range(n_restarts):
            if attempt == 0:
                p_init = p0_arr.copy()
            else:
                p_init = p0_arr + 0.25 * rng.standard_normal(p0_arr.size) * (np.abs(p0_arr) + 1e-3)
                for ii, (lo_b, hi_b) in enumerate(bounds):
                    if lo_b is not None:
                        p_init[ii] = max(p_init[ii], float(lo_b))
                    if hi_b is not None:
                        p_init[ii] = min(p_init[ii], float(hi_b))

            def nll(params):
                try:
                    shape = f_shape(params)
                    if not np.all(np.isfinite(shape)):
                        return 1e30
                    shape = np.clip(shape, 0.0, None)
                    s = float(np.sum(shape))
                    if s < 1e-12:
                        return 1e30
                    N = n_data / s
                    pred_i = np.clip(N * shape, 1e-12, None)
                    return float(-np.sum(data * np.log(pred_i) - pred_i))
                except Exception:
                    return 1e30

            try:
                res = _sp_minimize(
                    nll, p_init, method="L-BFGS-B", bounds=bounds,
                    options={"maxiter": 2000, "ftol": 1e-15, "gtol": 1e-9},
                )
                if np.isfinite(res.fun) and res.fun < best_fun:
                    best_fun = res.fun
                    best_x = res.x
            except Exception:
                continue

        if best_x is None or not np.isfinite(best_fun):
            return False, None, None, None

        shape = np.clip(f_shape(best_x), 0.0, None)
        s = float(np.sum(shape))
        if s < 1e-12:
            return False, None, None, None
        N = n_data / s
        pred_out = N * shape
        logL = float(np.sum(data * np.log(np.clip(pred_out, 1e-12, None)) - pred_out))
        npar = len(p0_arr) + 1
        return True, logL, np.clip(pred_out, 0.0, None), npar

    def _record(name: str, logL: float, pred: np.ndarray, npar: int) -> None:
        preds[name] = np.clip(pred, 0.0, None)
        resid = data - preds[name]
        chi2 = float(np.sum(resid ** 2 / np.maximum(v_data, 1.0)))
        ndf = max(int(data.size) - npar, 1)
        stats[name] = (chi2, ndf)

    def _lrt_scan(fit_func, orders):
        """Sequential LRT: scan orders, accept n when 2·ΔlogL ~ χ²(1) p≥5%."""
        prev_logL = prev_pred = prev_npar = None
        for order in orders:
            ok, logL, pred, npar = fit_func(order)
            if not ok:
                continue
            if prev_logL is None:
                prev_logL, prev_pred, prev_npar = logL, pred, npar
                continue
            stat = max(2.0 * (logL - prev_logL), 0.0)
            pval = float(_sp_chi2.sf(stat, df=1))
            if pval < 0.05:
                prev_logL, prev_pred, prev_npar = logL, pred, npar
            else:
                break
        return prev_logL, prev_pred, prev_npar

    # ── 1. BWZRedux: exp(a₁·t + a₂·t²) × BW(m) ─────────────────────────────
    try:
        def _bwzredux(params):
            a1, a2 = params
            return np.exp(np.clip(a1 * t + a2 * t * t, -30, 30)) * bw

        ok, logL, pred, npar = _fit_nll(_bwzredux, [0.0, 0.0], [(-20.0, 20.0), (-20.0, 20.0)])
        if ok:
            _record("BWZRedux", logL, pred, npar)
    except Exception:
        pass

    # ── 2. BWZGamma: exp(a₁·t) × BW(m) × (1 + f/m²) ────────────────────────
    try:
        def _bwzgamma(params):
            a1, ff = params
            factor = np.clip(1.0 + ff / np.maximum(centers ** 2, 1.0), 1e-12, None)
            return np.exp(np.clip(a1 * t, -30, 30)) * bw * factor

        ok, logL, pred, npar = _fit_nll(_bwzgamma, [0.0, 0.0], [(-20.0, 20.0), (-1e4, 1e4)])
        if ok:
            _record("BWZGamma", logL, pred, npar)
    except Exception:
        pass

    # ── 3. BWZ×Bern_n: exp(a·t) × BW(m) × Bern_n(t),  LRT n = 1,2,3 ────────
    try:
        def _bwzbern_fit(n_ord):
            basis = _bernstein_matrix(t, n_ord)
            n_bern = n_ord + 1
            init_c = max(n_data / max(n_bern * float(np.sum(bw)), 1e-12), 1e-3)
            p0 = np.array([0.0] + [init_c] * n_bern)
            bnds = [(-10.0, 10.0)] + [(0.0, None)] * n_bern

            def _shape(params):
                a = params[0]
                c = np.abs(params[1:])
                return np.exp(np.clip(a * t, -30, 30)) * bw * (basis @ c)

            return _fit_nll(_shape, p0, bnds)

        logL, pred, npar = _lrt_scan(_bwzbern_fit, [1, 2, 3])
        if pred is not None:
            _record("BWZ\u00d7Bern", logL, pred, npar)
    except Exception:
        pass

    # ── 4. Landau×Bern_n: Landau(m;mγ,σ) × Bern_n(t),  LRT n = 0,1,2,3 ────
    try:
        def _landaubern_fit(n_ord):
            basis = _bernstein_matrix(t, n_ord)
            n_bern = n_ord + 1
            mg0 = x_min + 0.25 * x_span
            sg0 = 0.15 * x_span
            if n_ord == 0:
                p0 = np.array([mg0, sg0])
                bnds = [(x_min - x_span, x_max + x_span), (0.5, x_span)]

                def _shape(params):
                    mg, sg = params[0], max(params[1], 0.5)
                    return np.clip(_sp_moyal.pdf(centers, loc=mg, scale=sg), 1e-30, None)
            else:
                init_c = max(n_data / max(n_bern, 1), 1e-3)
                p0 = np.array([mg0, sg0] + [init_c] * n_bern)
                bnds = [(x_min - x_span, x_max + x_span), (0.5, x_span)] + [(0.0, None)] * n_bern

                def _shape(params):
                    mg, sg = params[0], max(params[1], 0.5)
                    c = np.abs(params[2:])
                    landau = np.clip(_sp_moyal.pdf(centers, loc=mg, scale=sg), 1e-30, None)
                    return landau * (basis @ c)

            return _fit_nll(_shape, p0, bnds, n_restarts=3)

        logL, pred, npar = _lrt_scan(_landaubern_fit, [0, 1, 2, 3])
        if pred is not None:
            _record("Landau\u00d7Bern", logL, pred, npar)
    except Exception:
        pass

    # ── 5. Bernstein_n: Σ cᵢ·Bᵢ,ₙ(t),  LRT n = 3,4,5,6 ───────────────────
    try:
        def _bernstein_fit(n_ord):
            basis = _bernstein_matrix(t, n_ord)
            n_bern = n_ord + 1
            p0 = np.full(n_bern, max(n_data / n_bern, 1e-3))
            bnds = [(0.0, None)] * n_bern

            def _shape(params):
                return basis @ np.abs(params)

            return _fit_nll(_shape, p0, bnds)

        logL, pred, npar = _lrt_scan(_bernstein_fit, [3, 4, 5, 6])
        if pred is not None:
            _record("Bernstein", logL, pred, npar)
    except Exception:
        pass

    # ── 6. S-PowerLaw: Σᵢ aᵢ·m^bᵢ,  LRT 1→2 terms ─────────────────────────
    try:
        def _powlaw_fit(n_terms):
            p0 = []
            for i in range(n_terms):
                p0 += [max(n_data / n_terms, 1e-3), -4.0 - 0.5 * i]
            p0 = np.array(p0)
            bnds = [(0.0, None), (-30.0, 5.0)] * n_terms
            m_safe = np.maximum(centers, 1.0)

            def _shape(params):
                result = np.zeros_like(centers)
                params_clamped = params.copy()
                for i in range(n_terms):
                    a = max(params_clamped[2 * i], 0.0)
                    b = params_clamped[2 * i + 1]
                    result += a * np.power(m_safe, b)
                return result

            return _fit_nll(_shape, p0, bnds)

        logL, pred, npar = _lrt_scan(_powlaw_fit, [1, 2])
        if pred is not None:
            _record("S-PowerLaw", logL, pred, npar)
    except Exception:
        pass

    # ── 7. S-Exponential: Σᵢ aᵢ·exp(bᵢ·m),  LRT 1→2 terms ─────────────────
    try:
        def _sexp_fit(n_terms):
            p0 = []
            for i in range(n_terms):
                p0 += [max(n_data / n_terms, 1e-3), -0.06 * (i + 1)]
            p0 = np.array(p0)
            bnds = [(0.0, None), (-2.0, 0.0)] * n_terms

            def _shape(params):
                result = np.zeros_like(centers)
                for i in range(n_terms):
                    a = max(params[2 * i], 0.0)
                    b = params[2 * i + 1]
                    result += a * np.exp(np.clip(b * centers, -30, 30))
                return result

            return _fit_nll(_shape, p0, bnds)

        logL, pred, npar = _lrt_scan(_sexp_fit, [1, 2])
        if pred is not None:
            _record("S-Exponential", logL, pred, npar)
    except Exception:
        pass

    # ── Select best model by χ²/ndf ──────────────────────────────────────────
    _MODEL_NAMES = [
        "BWZRedux", "BWZGamma", "BWZ\u00d7Bern", "Landau\u00d7Bern",
        "Bernstein", "S-PowerLaw", "S-Exponential",
    ]
    best_name = "fallback_data"
    best_chi2ndf = float("inf")
    best_pred = h_data.copy()

    for name in _MODEL_NAMES:
        if name not in stats:
            continue
        chi2, ndf = stats[name]
        chi2ndf = chi2 / float(max(ndf, 1))
        if np.isfinite(chi2ndf) and chi2ndf < best_chi2ndf:
            best_chi2ndf = chi2ndf
            best_name = name
            best_pred = preds[name]

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


def _fit_signal_dcb_full(
    centers: np.ndarray,
    mass_bins: np.ndarray,
    h_signal: np.ndarray,
    args: argparse.Namespace,
    cat_idx: int = 0,
) -> Tuple[np.ndarray, Dict[str, float], float]:
    """Fit a double-sided Crystal Ball to a binned signal histogram.

    Returns (pred_curve_at_centers, params_dict, chi2ndf).
    params_dict keys: 'mean', 'sigma', 'alpha_l', 'n_l', 'alpha_r', 'n_r'.

    chi2/ndf is computed following the AN2019_090 convention (Sec. 7.1):
      - chi2 = sum_i (obs_i - exp_i)^2 / obs_i  over bins where obs_i > 0
      - ndf  = n_nonzero_bins - 6  (6 free DCB shape parameters)
    This gives chi2/ndf ~ 1 for a well-fitting DCB, consistent with the AN.

    The previous implementation used pred_i in the denominator and masked on
    pred > eps, which gave artificially low chi2 when the predicted curve had
    many nonzero bins but the sparse weighted histogram had few filled bins.
    """
    empty = np.zeros_like(centers, dtype=np.float64)
    if h_signal.size == 0 or float(np.sum(h_signal)) <= 0.0:
        return empty, {}, -1.0

    tag = f"sdcb_{cat_idx}"
    try:
        edges_arr = array("d", [float(x) for x in mass_bins])
        htmp = ROOT.TH1D(f"htmp_{tag}", f"htmp_{tag}", len(mass_bins) - 1, edges_arr)
        htmp.Sumw2()
        for ibin, y in enumerate(h_signal, start=1):
            htmp.SetBinContent(ibin, float(max(y, 0.0)))
            htmp.SetBinError(ibin, math.sqrt(max(float(y), 0.0)))

        m = ROOT.RooRealVar(f"m_{tag}", f"m_{tag}", float(args.mass_min), float(args.mass_max))
        data = ROOT.RooDataHist(f"data_{tag}", f"data_{tag}", ROOT.RooArgList(m), htmp)

        mean  = ROOT.RooRealVar(f"mean_{tag}",   f"mean_{tag}",  125.0, 123.0, 127.0)
        sigma = ROOT.RooRealVar(f"sigma_{tag}",  f"sigma_{tag}",   1.6,   0.2,   6.0)
        aL    = ROOT.RooRealVar(f"alphaL_{tag}", f"alphaL_{tag}",  1.5,   0.2,   6.0)
        nL    = ROOT.RooRealVar(f"nL_{tag}",     f"nL_{tag}",      5.0,   0.5,  80.0)
        aR    = ROOT.RooRealVar(f"alphaR_{tag}", f"alphaR_{tag}",  2.0,   0.2,   6.0)
        nR    = ROOT.RooRealVar(f"nR_{tag}",     f"nR_{tag}",      8.0,   0.5,  80.0)

        dcb = ROOT.RooCrystalBall(
            f"dcb_{tag}", f"dcb_{tag}", m, mean, sigma, aL, nL, aR, nR,
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
        for j, x in enumerate(centers):
            m.setVal(float(x))
            pred[j] = max(float(dcb.getVal(norm_set)), 0.0)

        total_sig = float(np.sum(h_signal))
        if np.sum(pred) > 0.0:
            pred *= total_sig / float(np.sum(pred))

        params: Dict[str, float] = {
            "mean":    float(mean.getVal()),
            "sigma":   float(sigma.getVal()),
            "alpha_l": float(aL.getVal()),
            "n_l":     float(nL.getVal()),
            "alpha_r": float(aR.getVal()),
            "n_r":     float(nR.getVal()),
        }

        # ----------------------------------------------------------------
        # chi2/ndf following AN2019_090 Sec. 7.1 convention:
        #   chi2 = sum_{i: obs_i > 0} (obs_i - exp_i)^2 / obs_i
        #   ndf  = n_nonzero_bins - 6   (6 free DCB shape params)
        # This matches "chi2/ndf ~ 1 for good DCB fits" stated in the AN.
        # Using obs_i (not exp_i) in the denominator is the standard
        # Neyman chi2 for histograms with non-negligible per-bin counts.
        # ----------------------------------------------------------------
        nonzero_mask = h_signal > 0.0
        n_nonzero = int(np.count_nonzero(nonzero_mask))
        ndf = max(n_nonzero - 6, 1)   # 6 free DCB shape parameters
        if n_nonzero > 0:
            chi2ndf = float(
                np.sum(
                    np.where(
                        nonzero_mask,
                        (h_signal - pred) ** 2 / np.maximum(h_signal, 1e-12),
                        0.0,
                    )
                )
            ) / float(ndf)
        else:
            chi2ndf = -1.0

        return pred, params, chi2ndf

    except Exception:
        return empty, {}, -1.0


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
    """Build per-category yield summary.

    Full-window yields (mc_full_total, sig_full_total) are retained for
    data/MC comparisons.  HWHM-windowed yields (sig_hwhm_total,
    mc_hwhm_total) follow the analysis note convention: S and B are
    integrated within the HWHM of the per-category DCB-fitted signal peak
    (AN2019_090 Sec. 6.1, Fig. 34 caption).  These are the quantities used
    for the expected significance S/sqrt(B) and Asimov Z in the summary table.
    """
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
            # Signal normalization (xsec x lumi x BR x eff) is independent
            # of the data-driven sideband MC scale; do NOT apply mc_weight_scale
            s_w_scaled = s_w
            for i in range(ncat):
                m_sig = (s_idx == i)
                if np.any(m_sig):
                    s_full = float(np.sum(s_w_scaled[m_sig]))
                    s_e2_full = float(np.sum(np.square(s_w_scaled[m_sig])))
                    summary[i]["sig_full_total"] = summary[i].get("sig_full_total", 0.0) + s_full
                    summary[i]["sig_full_total_err2"] = summary[i].get("sig_full_total_err2", 0.0) + s_e2_full

    # ------------------------------------------------------------------
    # HWHM-windowed S and B yields (AN2019_090 Sec. 6.1 convention).
    # For each category fit a DCB to the signal mass histogram, extract
    # sigma -> HWHM = 1.1774 * sigma, then integrate S and B within
    # [mean - HWHM, mean + HWHM].  This matches how S/sqrt(B) and the
    # combined Asimov Z are defined in the analysis note.
    # ------------------------------------------------------------------
    mass_bins_sum = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    centers_sum = 0.5 * (mass_bins_sum[:-1] + mass_bins_sum[1:])

    for i in range(ncat):
        lo_cat = float(edges[i])
        hi_cat = float(edges[i + 1])

        # Build signal mass histogram for this category (no mc_scale on signal)
        sig_m_chunks: List[np.ndarray] = []
        sig_w_chunks_hwhm: List[np.ndarray] = []
        for payload in era_payloads.values():
            s_sc   = np.asarray(payload["sig_full_scores"],  dtype=np.float64)  # type: ignore[index]
            s_mass = np.asarray(payload["sig_full_mass"],    dtype=np.float64)  # type: ignore[index]
            s_w    = np.asarray(payload["sig_full_weight"],  dtype=np.float64)  # type: ignore[index]
            if s_sc.size == 0 or s_mass.size == 0 or s_w.size == 0:
                continue
            sel = (s_sc >= lo_cat) & (s_sc < hi_cat)
            if np.any(sel):
                sig_m_chunks.append(s_mass[sel])
                sig_w_chunks_hwhm.append(s_w[sel])

        sig_m_cat  = concat_chunks(sig_m_chunks)
        sig_w_cat  = concat_chunks(sig_w_chunks_hwhm)

        h_sig_cat, _ = weighted_hist_and_var(
            sig_m_cat, mass_bins_sum,
            sig_w_cat if sig_w_cat.size > 0 else None,
        )

        # Fit DCB to extract mean and sigma; fall back to global defaults
        _pred_dcb, _params_dcb, _chi2_dcb = _fit_signal_dcb_full(
            centers_sum, mass_bins_sum, h_sig_cat, args, cat_idx=i,
        )
        if _params_dcb:
            sig_mean_i = _params_dcb["mean"]
            sig_sigma_i = _params_dcb["sigma"]
        else:
            sig_mean_i = 125.09
            sig_sigma_i = 2.0  # GeV fallback (wider window on fit failure)

        mass_lo_i, mass_hi_i = _hwhm_window_from_sigma(sig_mean_i, sig_sigma_i)
        summary[i]["hwhm_mass_lo"] = mass_lo_i
        summary[i]["hwhm_mass_hi"] = mass_hi_i

        s_hwhm, b_hwhm = _sig_bkg_in_hwhm(
            era_payloads, lo_cat, hi_cat, mass_lo_i, mass_hi_i, mc_weight_scale,
        )
        summary[i]["sig_hwhm_total"] = s_hwhm
        summary[i]["mc_hwhm_total"]  = b_hwhm

    return summary


def build_category_summary_table(
    edges: np.ndarray,
    summary: Dict[int, Dict[str, float]],
) -> Tuple[List[str], List[List[str]]]:
    """Build the per-category summary table.

    S/sqrt(B) and Asimov Z are computed using signal and background yields
    integrated within the HWHM of the per-category signal peak, following
    the analysis note convention (AN2019_090 Sec. 6.1, Fig. 34 caption).
    The full-window MC yield is still shown for data/MC comparison.
    """
    columns = [
        "Category",
        "BDT range",
        "Data sidebands",
        "MC sidebands ± stat",
        "Data/MC (sidebands)",
        "MC full (110-150) ± stat",
        "S/sqrt(B) [HWHM]",
        "Asimov Z (cat) [HWHM]",
    ]
    rows: List[List[str]] = []
    # Accumulate terms for combined Z = sqrt(sum_i 2[(s+b)log(1+s/b) - s])
    # using HWHM-windowed yields per category (AN2019_090 convention)
    combined_z_sum = 0.0
    for i in range(len(edges) - 1):
        vals = summary.get(i, {})
        d_side = float(vals.get("data_side_count", 0.0))
        m_side = float(vals.get("mc_side_total", 0.0))
        me_side = math.sqrt(max(float(vals.get("mc_side_total_err2", 0.0)), 0.0))
        m_full = float(vals.get("mc_full_total", 0.0))
        me_full = math.sqrt(max(float(vals.get("mc_full_total_err2", 0.0)), 0.0))
        ratio = d_side / m_side if m_side > 0 else 0.0

        # Use HWHM-windowed yields for sensitivity metrics (AN note convention)
        s_hwhm = float(vals.get("sig_hwhm_total", 0.0))
        b_hwhm = float(vals.get("mc_hwhm_total", 0.0))

        signif = (s_hwhm / math.sqrt(b_hwhm)) if b_hwhm > 0.0 else 0.0

        # Per-category Asimov Z term (HWHM window)
        if b_hwhm > 1e-12 and s_hwhm > 0.0:
            term = 2.0 * ((s_hwhm + b_hwhm) * math.log1p(s_hwhm / b_hwhm) - s_hwhm)
            combined_z_sum += max(term, 0.0)
            cat_z = math.sqrt(max(term, 0.0))
        else:
            cat_z = 0.0
        rows.append([
            str(i),
            f"[{edges[i]:.3f}, {edges[i+1]:.3f})",
            f"{d_side:.0f}",
            f"{m_side:.2f} ± {me_side:.2f}",
            f"{ratio:.3f}",
            f"{m_full:.2f} ± {me_full:.2f}",
            f"{signif:.4f}",
            f"{cat_z:.4f}",
        ])
    # Footer: combined Z across all categories
    combined_z = math.sqrt(max(combined_z_sum, 0.0))
    rows.append([
        "Combined", "all",
        "", "", "", "",
        "", f"{combined_z:.4f}",
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
        "BWZRedux",
        "BWZGamma",
        "BWZ\u00d7Bern",
        "Landau\u00d7Bern",
        "Bernstein",
        "S-PowerLaw",
        "S-Exponential",
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
                    sig_weight_chunks.append(s_w[ssel])  # no mc_scale for signal

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
            "BWZRedux":      ("#d62728", "-"),
            "BWZGamma":      ("#2ca02c", "-"),
            "BWZ\u00d7Bern":    ("#1f77b4", "-"),
            "Landau\u00d7Bern": ("#ff7f0e", "-"),
            "Bernstein":     ("#9467bd", "--"),
            "S-PowerLaw":    ("#8c564b", ":"),
            "S-Exponential": ("#e377c2", "-."),
        }

        # Draw all background models
        for model_name in ["BWZRedux", "BWZGamma", "BWZ\u00d7Bern", "Landau\u00d7Bern", "Bernstein", "S-PowerLaw", "S-Exponential"]:
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

        model_order = ["BWZRedux", "BWZGamma", "BWZ\u00d7Bern", "Landau\u00d7Bern", "Bernstein", "S-PowerLaw", "S-Exponential"]
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

    fig.subplots_adjust(left=0.07, right=0.97, top=0.930, bottom=0.055,
                        hspace=0.35, wspace=0.28)
    fig.suptitle(
        "Per-category dimuon mass: background model comparison and ratio-to-best  (AN Fig. 36 analog)",
        fontsize=13,
        fontweight="bold",
        y=0.968,
    )
    return fig, fit_columns, fit_rows


def build_vbf_signal_dcb_grid(
    era_payloads: Dict[str, Dict[str, object]],
    edges: np.ndarray,
    args: argparse.Namespace,
    mass_bins: np.ndarray,
    mc_weight_scale: float = 1.0,
    cat_summary: Optional[Dict[int, Dict[str, float]]] = None,
) -> Tuple[plt.Figure, List[str], List[List[str]]]:
    """Per-BDT-category VBF H(125)->mumu signal DCB extraction plots.

    Reproduces AN2019_090 Fig. 38 for the VBF channel.

    Histograms the VBF signal MC mass in each category summed across all eras
    (Events / GeV on the y-axis, matching the AN figures), fits a double-sided
    Crystal Ball (DCB), and draws the result.

    chi2/ndf is computed as sum_{i:obs>0}(obs-exp)^2/obs / (n_nonzero-6),
    giving chi2/ndf ~ 1 for good fits per AN2019_090 Sec. 7.1.

    Returns (figure, table_columns, table_rows).
    """
    ncat = len(edges) - 1
    ncols = 2 if ncat > 1 else 1
    nrows = int(math.ceil(ncat / float(ncols)))

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(14.0, 5.0 * nrows + 0.5),
        squeeze=False,
    )

    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    bin_width = float(mass_bins[1] - mass_bins[0])  # GeV per bin

    dcb_cols: List[str] = [
        "Category",
        "BDT range",
        "N_sig",
        "Mean [GeV]",
        "Sigma [GeV]",
        "alpha_L",
        "n_L",
        "alpha_R",
        "n_R",
        "FWHM [GeV]",
        "S/sqrt(B)",
        "Asimov Z",
    ]
    dcb_rows: List[List[str]] = []

    for i in range(ncat):
        row_idx = i // ncols
        col_idx = i % ncols
        ax = axes[row_idx, col_idx]

        lo, hi = float(edges[i]), float(edges[i + 1])

        sig_mass_chunks: List[np.ndarray] = []
        sig_weight_chunks: List[np.ndarray] = []

        for payload in era_payloads.values():
            s_sc   = np.asarray(payload["sig_full_scores"],  dtype=np.float64)
            s_mass = np.asarray(payload["sig_full_mass"],    dtype=np.float64)
            s_w    = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
            if s_sc.size == 0 or s_mass.size == 0 or s_w.size == 0:
                continue
            ssel = (s_sc >= lo) & (s_sc < hi)
            if np.any(ssel):
                sig_mass_chunks.append(s_mass[ssel])
                sig_weight_chunks.append(s_w[ssel])  # no mc_scale for signal

        sig_mass = concat_chunks(sig_mass_chunks)
        sig_w    = concat_chunks(sig_weight_chunks)

        h_sig, _ = weighted_hist_and_var(
            sig_mass, mass_bins, sig_w if sig_w.size > 0 else None,
        )

        pred_dcb, params, chi2ndf = _fit_signal_dcb_full(
            centers, mass_bins, h_sig, args, cat_idx=i,
        )

        n_sig_total = float(np.sum(h_sig))

        # Retrieve per-category S and B from pre-computed summary.
        # Use HWHM-windowed yields (AN2019_090 Sec. 6.1 convention).
        _summary_i = (cat_summary or {}).get(i, {})
        s_yield = float(_summary_i.get("sig_hwhm_total", 0.0))
        b_yield = float(_summary_i.get("mc_hwhm_total", 0.0))
        # Fall back to full-window sig yield if HWHM yields are unavailable
        if s_yield == 0.0:
            s_yield = float(_summary_i.get("sig_full_total", n_sig_total))
        sbsqrtb = s_yield / math.sqrt(b_yield) if b_yield > 0.0 else 0.0
        if b_yield > 0.0 and s_yield > 0.0:
            asimov_z = math.sqrt(
                max(0.0, 2.0 * ((s_yield + b_yield) * math.log1p(s_yield / b_yield) - s_yield))
            )
        else:
            asimov_z = 0.0

        # Normalize to Events / GeV (matching AN Fig. 38 y-axis)
        h_sig_gev    = h_sig / bin_width
        pred_dcb_gev = pred_dcb / bin_width

        ax.step(
            centers, h_sig_gev, where="mid",
            color="#4477CC", linewidth=1.4,
            label=f"VBF H(125) MC  (N = {n_sig_total:.1f})",
            zorder=2,
        )
        ax.fill_between(
            centers, h_sig_gev,
            step="mid", color="#4477CC", alpha=0.20, zorder=1,
        )

        if np.sum(pred_dcb_gev) > 0.0:
            ax.step(
                centers, pred_dcb_gev, where="mid",
                color="#CC4444", linewidth=2.2,
                label="DCB fit", zorder=3,
            )

        ax.set_title(f"Cat {i}:  BDT [{lo:.3f}, {hi:.3f})", fontsize=11, pad=4)
        ax.set_xlabel("$m_{\\mu\\mu}$ [GeV]", fontsize=10)
        ax.set_ylabel("Events / GeV")
        ax.grid(alpha=0.25)

        if params:
            fwhm = 2.355 * params["sigma"]
            stats_text = (
                f"$\\mu$ = {params['mean']:.3f} GeV\n"
                f"$\\sigma$ = {params['sigma']:.3f} GeV\n"
                f"$\\alpha_L$ = {params['alpha_l']:.3f},  $n_L$ = {params['n_l']:.2f}\n"
                f"$\\alpha_R$ = {params['alpha_r']:.3f},  $n_R$ = {params['n_r']:.2f}\n"
                f"FWHM $\\approx$ {fwhm:.3f} GeV"
            )
        else:
            fwhm = float("nan")
            stats_text = f"N = {n_sig_total:.1f}\n(DCB fit failed)"

        ax.text(
            0.97, 0.95,
            stats_text,
            transform=ax.transAxes,
            va="top", ha="right",
            fontsize=8.5,
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor="white", alpha=0.88, edgecolor="#b0b0b0",
            ),
        )
        ax.legend(loc="upper left", fontsize=9)

        chi2_str_row = f"{chi2ndf:.3f}" if chi2ndf >= 0.0 else "n/a"  # kept for internal use
        dcb_rows.append([
            str(i),
            f"[{lo:.3f}, {hi:.3f})",
            f"{n_sig_total:.2f}",
            f"{params['mean']:.4f}" if params else "n/a",
            f"{params['sigma']:.4f}" if params else "n/a",
            f"{params['alpha_l']:.4f}" if params else "n/a",
            f"{params['n_l']:.3f}" if params else "n/a",
            f"{params['alpha_r']:.4f}" if params else "n/a",
            f"{params['n_r']:.3f}" if params else "n/a",
            f"{fwhm:.3f}" if math.isfinite(fwhm) else "n/a",
            f"{sbsqrtb:.3f}",
            f"{asimov_z:.3f}",
        ])

    for j in range(ncat, nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    fig.subplots_adjust(
        left=0.07, right=0.97, top=0.925, bottom=0.08,
        hspace=0.42, wspace=0.28,
    )
    fig.suptitle(
        "VBF H(125) signal DCB extraction per BDT category  (AN Fig. 38 analog)",
        fontsize=13, fontweight="bold", y=0.965,
    )
    return fig, dcb_cols, dcb_rows


def build_vbf_bdt_transformed_figure(
    era_payloads: Dict[str, Dict[str, object]],
    edges: np.ndarray,
    args: argparse.Namespace,
    mass_bins: np.ndarray,
    mc_weight_scale: float = 1.0,
    n_score_bins: int = 25,
) -> plt.Figure:
    """Reproduce AN2019_090 Figure 34 (right panel) for the VBF channel.

    Top pad: signal (blue) and background (red) normalised histograms in
    atanh-transformed BDT space with category boundary lines.
    Bottom pad: S/sqrt(B) per BDT bin, where S and B are integrated
    within the HWHM of the expected signal peak (AN2019_090 convention).
    """
    # ---- collect signal (unscaled) and background (scaled) ----
    sig_score_chunks: List[np.ndarray] = []
    sig_mass_chunks:  List[np.ndarray] = []
    sig_w_chunks:     List[np.ndarray] = []
    bkg_score_chunks: List[np.ndarray] = []
    bkg_mass_chunks:  List[np.ndarray] = []
    bkg_w_chunks:     List[np.ndarray] = []

    for payload in era_payloads.values():
        s_sc  = np.asarray(payload["sig_full_scores"],  dtype=np.float64)
        s_m   = np.asarray(payload["sig_full_mass"],    dtype=np.float64)
        s_w   = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
        if s_sc.size > 0 and s_m.size > 0 and s_w.size > 0:
            sig_score_chunks.append(s_sc)
            sig_mass_chunks.append(s_m)
            sig_w_chunks.append(s_w)  # physics-normalised, no mc_scale

        for _proc, b_sc, b_m, b_w in payload["mc_full_groups"]:
            b_sc = np.asarray(b_sc, dtype=np.float64)
            b_m  = np.asarray(b_m,  dtype=np.float64)
            b_w  = np.asarray(b_w,  dtype=np.float64) * float(mc_weight_scale)
            if b_sc.size > 0:
                bkg_score_chunks.append(b_sc)
                bkg_mass_chunks.append(b_m)
                bkg_w_chunks.append(b_w)

    sig_scores = concat_chunks(sig_score_chunks)
    sig_mass   = concat_chunks(sig_mass_chunks)
    sig_w      = concat_chunks(sig_w_chunks)
    bkg_scores = concat_chunks(bkg_score_chunks)
    bkg_mass   = concat_chunks(bkg_mass_chunks)
    bkg_w      = concat_chunks(bkg_w_chunks)

    if sig_scores.size == 0 or bkg_scores.size == 0:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.text(0.5, 0.5, "No signal/background data", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    # ---- transform to atanh space ----
    t_sig = _bdt_transform(sig_scores)
    t_bkg = _bdt_transform(bkg_scores)
    t_edges = _bdt_transform(np.clip(edges[1:-1], 1e-6, 1.0 - 1e-6))

    t_min = float(np.percentile(t_sig, 1.0))
    t_max = float(np.percentile(t_sig, 99.5))
    t_min = min(t_min, float(np.percentile(t_bkg, 1.0)))
    t_max = max(t_max, float(np.percentile(t_bkg, 99.5)))
    score_edges = np.linspace(t_min, t_max, n_score_bins + 1)

    # ---- fit inclusive signal DCB to get mean and sigma -> HWHM ----
    centers_mass = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    h_sig_incl, _ = weighted_hist_and_var(sig_mass, mass_bins, sig_w if sig_w.size > 0 else None)
    _pred_incl, _par_incl, _chi2_incl = _fit_signal_dcb_full(
        centers_mass, mass_bins, h_sig_incl, args, cat_idx=999,
    )
    if _par_incl:
        sig_mean = _par_incl["mean"]
        hwhm = 1.177 * _par_incl["sigma"]
    else:
        sig_mean = 125.09
        hwhm = 2.0  # GeV fallback

    # ---- S/sqrt(B) per score bin, integrated within HWHM mass window ----
    mass_lo = sig_mean - hwhm
    mass_hi = sig_mean + hwhm

    sig_in_window = (sig_mass >= mass_lo) & (sig_mass <= mass_hi)
    bkg_in_window = (bkg_mass >= mass_lo) & (bkg_mass <= mass_hi)

    s_per_bin = np.zeros(n_score_bins, dtype=np.float64)
    b_per_bin = np.zeros(n_score_bins, dtype=np.float64)

    if np.any(sig_in_window):
        s_idx = np.clip(np.digitize(t_sig[sig_in_window], score_edges) - 1, 0, n_score_bins - 1)
        np.add.at(s_per_bin, s_idx, sig_w[sig_in_window])
    if np.any(bkg_in_window):
        b_idx = np.clip(np.digitize(t_bkg[bkg_in_window], score_edges) - 1, 0, n_score_bins - 1)
        np.add.at(b_per_bin, b_idx, bkg_w[bkg_in_window])

    sqrtb = np.sqrt(np.maximum(b_per_bin, 0.0))
    s_over_sqrtb = np.where(sqrtb > 0.0, s_per_bin / sqrtb, 0.0)

    # ---- normalised score histograms ----
    bin_widths = score_edges[1:] - score_edges[:-1]
    h_sig_score, _ = weighted_hist_and_var(t_sig, score_edges, sig_w if sig_w.size > 0 else None)
    h_bkg_score, _ = weighted_hist_and_var(t_bkg, score_edges, bkg_w if bkg_w.size > 0 else None)
    tot_s = float(np.sum(h_sig_score))
    tot_b = float(np.sum(h_bkg_score))
    h_sig_norm = h_sig_score / (tot_s * bin_widths) if tot_s > 0 else h_sig_score
    h_bkg_norm = h_bkg_score / (tot_b * bin_widths) if tot_b > 0 else h_bkg_score

    bin_centers = 0.5 * (score_edges[:-1] + score_edges[1:])

    # ---- figure ----
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(7, 6),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.04},
    )

    ax_top.step(bin_centers, h_sig_norm, where="mid",
                color="#1f77b4", linewidth=1.8, label="Signal $m_H=125$ GeV")
    ax_top.step(bin_centers, h_bkg_norm, where="mid",
                color="#d62728", linewidth=1.8, label="Background")

    for t_edge in t_edges:
        if score_edges[0] < t_edge < score_edges[-1]:
            ax_top.axvline(t_edge, color="black", linestyle="--", linewidth=1.1, alpha=0.7)
            ax_bot.axvline(t_edge, color="black", linestyle="--", linewidth=1.1, alpha=0.7)

    ax_top.set_ylabel("a.u.", fontsize=11)
    ax_top.tick_params(labelbottom=False)
    ax_top.legend(fontsize=10, loc="upper left")
    ax_top.grid(alpha=0.2)
    ax_top.set_xlim(score_edges[0], score_edges[-1])

    # CMS label
    ax_top.text(0.98, 0.96, "(13.6 TeV)", transform=ax_top.transAxes,
                ha="right", va="top", fontsize=9.5)
    ax_top.text(0.02, 0.96, "CMS", transform=ax_top.transAxes,
                ha="left", va="top", fontsize=11, fontweight="bold")

    ax_bot.step(bin_centers, s_over_sqrtb, where="mid",
                color="black", linewidth=1.6)
    ax_bot.fill_between(bin_centers, s_over_sqrtb, step="mid",
                        color="black", alpha=0.12)
    ax_bot.set_xlabel("atanh((BDT output+1)/2)", fontsize=11)
    ax_bot.set_ylabel("S/$\\sqrt{B}$", fontsize=11)
    ax_bot.grid(alpha=0.2)
    ax_bot.set_xlim(score_edges[0], score_edges[-1])
    ymax_bot = max(float(np.max(s_over_sqrtb)) * 1.25, 0.01)
    ax_bot.set_ylim(0.0, ymax_bot)

    hwhm_note = (
        f"HWHM window: [{mass_lo:.1f}, {mass_hi:.1f}] GeV\n"
        f"($\\sigma$={_par_incl['sigma']:.2f} GeV)"
        if _par_incl else "HWHM window (fallback)"
    )
    ax_bot.text(
        0.98, 0.95,
        hwhm_note,
        transform=ax_bot.transAxes,
        ha="right", va="top", fontsize=8,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
    )

    fig.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.10)
    fig.suptitle(
        "VBF BDT score and per-bin S/$\\sqrt{B}$ in signal HWHM window  (AN Fig. 34 analog)",
        fontsize=12, fontweight="bold", y=0.99,
    )
    return fig


def build_significance_evolution_figure(
    scan_history: List[Tuple[int, float]],
) -> plt.Figure:
    """Plot the evolution of combined Asimov Z as categories are added.

    Reproduces AN2019_090 Figure 33: each panel shows the significance as
    a function of the boundary position (expressed as 1 - signal efficiency
    accumulated from the high end of the BDT score), with the optimal
    boundary marked.  Here we plot a simpler but equivalent view: combined
    Z vs. number of categories, since the per-boundary scan points are not
    individually cached.

    Parameters
    ----------
    scan_history : list of (n_categories, combined_Z) pairs recorded by
        iterative_greedy_category_scan when logger is provided.
    """
    if not scan_history:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No scan history available\n(manual edges or single-era run)",
                ha="center", va="center", transform=ax.transAxes, fontsize=11)
        ax.set_title("Significance evolution  (AN Fig. 33 analog)", fontsize=12)
        return fig

    ns = [h[0] for h in scan_history]
    zs = [h[1] for h in scan_history]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ns, zs, "o-", color="#1f77b4", linewidth=2.0, markersize=7,
            markerfacecolor="white", markeredgewidth=2.0, zorder=3)

    # Mark the optimal (maximum) point
    best_idx = int(np.argmax(zs))
    ax.axvline(ns[best_idx], color="#d62728", linestyle="--", linewidth=1.3,
               label=f"Optimal: n={ns[best_idx]}, Z={zs[best_idx]:.4f}")
    ax.scatter([ns[best_idx]], [zs[best_idx]], color="#d62728", s=80, zorder=5)

    for n_val, z_val in zip(ns, zs):
        ax.annotate(
            f"{z_val:.3f}",
            xy=(n_val, z_val),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
            fontsize=9,
            color="#333333",
        )

    ax.set_xlabel("Number of BDT categories", fontsize=11)
    ax.set_ylabel("Combined Asimov $Z$ [HWHM window]", fontsize=11)
    ax.set_title(
        "Evolution of combined significance vs. number of VBF BDT categories\n"
        "(AN Fig. 33 analog — greedy iterative scan)",
        fontsize=11,
    )
    ax.set_xticks(ns)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # CMS label
    ax.text(0.98, 0.04, "(13.6 TeV)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9)
    ax.text(0.02, 0.96, "CMS", transform=ax.transAxes,
            ha="left", va="top", fontsize=11, fontweight="bold")

    fig.tight_layout()
    return fig


def build_fwhm_vs_bdt_figure(
    era_payloads: Dict[str, Dict[str, object]],
    args: argparse.Namespace,
    mass_bins: np.ndarray,
    n_score_bins: int = 10,
) -> plt.Figure:
    """Plot FWHM of the dimuon mass resolution vs. BDT score bin.

    Reproduces AN2019_090 Figure 22 (right panel) for the VBF channel.
    For each BDT score slice the signal mass distribution is fit with a DCB
    and the FWHM = 2.355 * sigma is extracted.  This demonstrates that the
    BDT score is correlated with dimuon mass resolution, motivating the
    categorization strategy.
    """
    sig_score_chunks: List[np.ndarray] = []
    sig_mass_chunks:  List[np.ndarray] = []
    sig_w_chunks:     List[np.ndarray] = []

    for payload in era_payloads.values():
        s_sc  = np.asarray(payload["sig_full_scores"],  dtype=np.float64)
        s_m   = np.asarray(payload["sig_full_mass"],    dtype=np.float64)
        s_w   = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
        if s_sc.size > 0 and s_m.size > 0 and s_w.size > 0:
            sig_score_chunks.append(s_sc)
            sig_mass_chunks.append(s_m)
            sig_w_chunks.append(s_w)

    sig_scores = concat_chunks(sig_score_chunks)
    sig_mass   = concat_chunks(sig_mass_chunks)
    sig_w      = concat_chunks(sig_w_chunks)

    if sig_scores.size == 0:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.text(0.5, 0.5, "No signal data available", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title("FWHM vs. VBF BDT score  (AN Fig. 22 analog)")
        return fig

    # Bin signal events by BDT score using equal-statistics quantile bins
    # (matches AN approach of plotting FWHM per BDT score slice)
    q_edges = np.quantile(sig_scores, np.linspace(0.0, 1.0, n_score_bins + 1))
    q_edges = np.unique(q_edges)
    n_actual = len(q_edges) - 1

    centers_mass = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    bin_centers_score = 0.5 * (q_edges[:-1] + q_edges[1:])
    # Use atanh-transformed bin centers for x-axis (matching AN Fig. 22)
    x_plot = _bdt_transform(np.clip(bin_centers_score, 1e-6, 1.0 - 1e-6))

    fwhm_vals: List[float] = []
    fwhm_errs: List[float] = []  # use sigma uncertainty proxy ~ sigma/sqrt(2N)

    for b in range(n_actual):
        lo_s = float(q_edges[b])
        hi_s = float(q_edges[b + 1])
        sel = (sig_scores >= lo_s) & (sig_scores < hi_s)
        if not np.any(sel):
            fwhm_vals.append(float("nan"))
            fwhm_errs.append(float("nan"))
            continue

        h_slice, _ = weighted_hist_and_var(
            sig_mass[sel], mass_bins,
            sig_w[sel] if sig_w.size > 0 else None,
        )

        _pred, _params, _chi2 = _fit_signal_dcb_full(
            centers_mass, mass_bins, h_slice, args, cat_idx=1000 + b,
        )

        if _params:
            sigma_val = _params["sigma"]
            fwhm_vals.append(2.355 * sigma_val)
            # Approximate sigma uncertainty from fit; use 10% as conservative proxy
            fwhm_errs.append(2.355 * sigma_val * 0.10)
        else:
            fwhm_vals.append(float("nan"))
            fwhm_errs.append(float("nan"))

    fwhm_arr = np.array(fwhm_vals, dtype=np.float64)
    fwhm_err_arr = np.array(fwhm_errs, dtype=np.float64)
    valid = np.isfinite(fwhm_arr)

    fig, ax = plt.subplots(figsize=(7, 4.5))

    if np.any(valid):
        ax.errorbar(
            x_plot[valid], fwhm_arr[valid],
            yerr=fwhm_err_arr[valid],
            fmt="o-", color="#1f77b4",
            linewidth=1.8, markersize=6,
            capsize=3, label="FWHM (2.355 × σ from DCB fit)",
        )

    ax.set_xlabel("atanh((VBF BDT output + 1) / 2)", fontsize=11)
    ax.set_ylabel("FWHM [GeV]", fontsize=11)
    ax.set_title(
        "Dimuon mass resolution (FWHM) vs. VBF BDT score\n"
        "(AN Fig. 22 analog — equal-statistics BDT slices)",
        fontsize=11,
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    # CMS label
    ax.text(0.98, 0.96, "(13.6 TeV)", transform=ax.transAxes,
            ha="right", va="top", fontsize=9)
    ax.text(0.02, 0.96, "CMS", transform=ax.transAxes,
            ha="left", va="top", fontsize=11, fontweight="bold")

    fig.tight_layout()
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


# ---------------------------------------------------------------------------
# Greedy iterative category boundary scan
# ---------------------------------------------------------------------------

_ATANH_EPS = 1e-7


def _bdt_transform(scores: np.ndarray) -> np.ndarray:
    """BDT' = atanh((BDT + 1) / 2), clipped to avoid ±inf.

    Works for BDT scores in either [-1, 1] or [0, 1].  Scores close to the
    boundary of [-1,1] map to atanh values near ±inf; clipping keeps things
    numerically stable.
    """
    x = np.clip(
        (np.asarray(scores, dtype=np.float64) + 1.0) / 2.0,
        _ATANH_EPS,
        1.0 - _ATANH_EPS,
    )
    return np.arctanh(x)


def _bdt_inverse(t: np.ndarray) -> np.ndarray:
    """Inverse of _bdt_transform: s = 2 * tanh(t) - 1."""
    return 2.0 * np.tanh(np.asarray(t, dtype=np.float64)) - 1.0


def _bin_yields_t(
    t_scores: np.ndarray,
    weights: np.ndarray,
    edges_t: np.ndarray,
) -> np.ndarray:
    """Weighted yields per bin, defined by *edges_t* in transformed BDT space."""
    bin_idx = np.clip(
        np.digitize(t_scores, edges_t) - 1,
        0,
        len(edges_t) - 2,
    )
    return np.bincount(bin_idx, weights=weights, minlength=len(edges_t) - 1).astype(np.float64)


def _bin_yields_t_hwhm(
    t_scores: np.ndarray,
    weights: np.ndarray,
    mass: np.ndarray,
    mass_lo: float,
    mass_hi: float,
    edges_t: np.ndarray,
) -> np.ndarray:
    """Weighted yields per bin within the HWHM mass window.

    Like _bin_yields_t but restricts to events with mass in [mass_lo, mass_hi].
    This implements the AN2019_090 Sec. 6.1 convention where the figure of merit
    for category optimization uses S and B integrated within the HWHM of the
    expected signal peak.
    """
    in_window = (mass >= mass_lo) & (mass <= mass_hi)
    if not np.any(in_window):
        return np.zeros(len(edges_t) - 1, dtype=np.float64)
    bin_idx = np.clip(
        np.digitize(t_scores[in_window], edges_t) - 1,
        0,
        len(edges_t) - 2,
    )
    return np.bincount(
        bin_idx, weights=weights[in_window], minlength=len(edges_t) - 1
    ).astype(np.float64)


def _combined_asimov_z(s_yields: np.ndarray, b_yields: np.ndarray) -> float:
    """Combined Asimov significance Z = sqrt(sum_i 2[(s+b)log(1+s/b) - s])."""
    total = 0.0
    for s, b in zip(s_yields, b_yields):
        if b < 1e-12 or s <= 0.0:
            continue
        total += 2.0 * ((s + b) * float(np.log1p(s / b)) - s)
    return float(np.sqrt(max(total, 0.0)))


def iterative_greedy_category_scan(
    sig_scores: np.ndarray,
    sig_weights: np.ndarray,
    bkg_scores: np.ndarray,
    bkg_weights: np.ndarray,
    max_n: int = 6,
    min_improvement: float = 0.01,
    min_sig: float = 1e-6,
    min_bkg: float = 1e-6,
    sig_mass: Optional[np.ndarray] = None,
    bkg_mass: Optional[np.ndarray] = None,
    hwhm_mass_lo: float = 0.0,
    hwhm_mass_hi: float = float("inf"),
    logger: logging.Logger = None,
) -> Tuple[np.ndarray, List[Tuple[int, float]]]:
    """Greedy iterative category boundary optimisation on atanh-transformed BDT.

    Algorithm
    ---------
    1. Transform BDT -> BDT' = atanh((BDT + 1) / 2).
    2. Generate 99 candidate boundary positions at every 1 % signal-efficiency
       quantile of the transformed signal distribution.
    3. Begin with n = 1 category (no interior boundaries).
    4. At each step try inserting each unused candidate; retain the one that
       maximises the combined Asimov significance
           Z = sqrt(sum_i 2 [(s_i+b_i) log(1 + s_i/b_i) - s_i]).
       where s_i and b_i are the yields integrated within the HWHM of the
       expected signal peak (AN2019_090 Sec. 6.1 convention).  If mass arrays
       are not provided, falls back to the full-window yields.
    5. A candidate is only considered if every resulting bin has at least
       *min_sig* weighted signal yield AND *min_bkg* weighted background yield.
    6. Add that boundary only when the relative improvement in Z is >=
       *min_improvement* (default 1 %).  Stop otherwise or when n == max_n.

    Returns
    -------
    edges : np.ndarray
        Interior edges back-transformed to the original BDT score space,
        bookended by 0.0 and 1.0.
    scan_history : list of (n_categories, combined_Z) pairs, one per accepted
        boundary step (plus the baseline with n=1).  Used by
        build_significance_evolution_figure to reproduce AN Fig. 33.
    """
    sig_scores = np.asarray(sig_scores, dtype=np.float64)
    bkg_scores = np.asarray(bkg_scores, dtype=np.float64)
    sig_weights = np.asarray(sig_weights, dtype=np.float64)
    bkg_weights = np.asarray(bkg_weights, dtype=np.float64)

    use_hwhm = (
        sig_mass is not None
        and bkg_mass is not None
        and np.isfinite(hwhm_mass_lo)
        and np.isfinite(hwhm_mass_hi)
        and hwhm_mass_hi > hwhm_mass_lo
    )
    if use_hwhm:
        sig_mass_arr = np.asarray(sig_mass, dtype=np.float64)
        bkg_mass_arr = np.asarray(bkg_mass, dtype=np.float64)
    else:
        sig_mass_arr = None
        bkg_mass_arr = None

    t_sig = _bdt_transform(sig_scores)
    t_bkg = _bdt_transform(bkg_scores)

    t_lo = float(np.min(t_sig))
    t_hi = float(np.max(t_sig))

    # Candidate boundaries: every 1 % signal-efficiency quantile
    q_levels = np.linspace(0.01, 0.99, 99)
    cands = np.unique(np.quantile(t_sig, q_levels))
    cands = cands[(cands > t_lo) & (cands < t_hi)]

    if logger:
        if use_hwhm:
            logger.info(
                "Greedy scan (atanh BDT, HWHM window [%.2f, %.2f] GeV): "
                "%d S events (yield=%.3f), %d B events (yield=%.3f), "
                "max_n=%d, min_improvement=%.1f%%",
                hwhm_mass_lo, hwhm_mass_hi,
                sig_scores.size, float(np.sum(sig_weights)),
                bkg_scores.size, float(np.sum(bkg_weights)),
                max_n, min_improvement * 100.0,
            )
        else:
            logger.info(
                "Greedy scan (atanh BDT, full mass window): "
                "%d S events (yield=%.3f), %d B events (yield=%.3f), "
                "max_n=%d, min_improvement=%.1f%%",
                sig_scores.size, float(np.sum(sig_weights)),
                bkg_scores.size, float(np.sum(bkg_weights)),
                max_n, min_improvement * 100.0,
            )

    def _get_yields(trial_edges):
        if use_hwhm:
            s_yields = _bin_yields_t_hwhm(
                t_sig, sig_weights, sig_mass_arr,
                hwhm_mass_lo, hwhm_mass_hi, trial_edges,
            )
            b_yields = _bin_yields_t_hwhm(
                t_bkg, bkg_weights, bkg_mass_arr,
                hwhm_mass_lo, hwhm_mass_hi, trial_edges,
            )
        else:
            s_yields = _bin_yields_t(t_sig, sig_weights, trial_edges)
            b_yields = _bin_yields_t(t_bkg, bkg_weights, trial_edges)
        return s_yields, b_yields

    # Baseline significance with 1 category
    baseline_edges = np.array([t_lo, t_hi])
    prev_z = _combined_asimov_z(*_get_yields(baseline_edges))

    # Record scan history: (n_categories, Z) — start with baseline
    scan_history: List[Tuple[int, float]] = [(1, prev_z)]

    current_interior: List[float] = []

    for step in range(max_n - 1):
        remaining = set(cands.tolist()) - set(current_interior)
        if not remaining:
            break

        best_cand_t: float | None = None
        best_z = -1.0
        for ct in remaining:
            trial_interior = sorted(current_interior + [ct])
            trial_edges = np.array([t_lo] + trial_interior + [t_hi])
            s_yields, b_yields = _get_yields(trial_edges)
            # Min-yield guard uses full-window yields
            s_full = _bin_yields_t(t_sig, sig_weights, trial_edges)
            b_full = _bin_yields_t(t_bkg, bkg_weights, trial_edges)
            if np.any(s_full < min_sig) or np.any(b_full < min_bkg):
                continue
            z = _combined_asimov_z(s_yields, b_yields)
            if z > best_z:
                best_z = z
                best_cand_t = ct

        if best_cand_t is None:
            break

        rel_gain = (best_z - prev_z) / (max(prev_z, 1e-12))
        if rel_gain < min_improvement:
            if logger:
                logger.info(
                    "  Step %d: best gain %.2f%% < threshold %.1f%%; "
                    "stopping at n=%d categories.",
                    step + 1, rel_gain * 100.0, min_improvement * 100.0,
                    len(current_interior) + 1,
                )
            break

        current_interior.append(best_cand_t)
        prev_z = best_z
        scan_history.append((len(current_interior) + 1, best_z))

        if logger:
            logger.info(
                "  Step %d: added t'=%.4f  -> n=%d cats, Z=%.4f (gain +%.2f%%)",
                step + 1, best_cand_t,
                len(current_interior) + 1, best_z, rel_gain * 100.0,
            )

    if not current_interior:
        return np.array([0.0, 1.0], dtype=np.float64), scan_history

    interior_orig = _bdt_inverse(np.array(current_interior))
    edges_orig = np.concatenate([[0.0], np.sort(interior_orig), [1.0]])
    edges_orig = np.unique(np.clip(edges_orig, 0.0, 1.0))
    return edges_orig, scan_history


def _gather_category_data(
    era_payloads: Dict[str, Dict[str, object]],
    lo: float,
    hi: float,
    mass_bins: np.ndarray,
    mc_weight_scale: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (h_data, v_data, h_bkg, h_sig) histograms for one BDT category."""
    data_chunks: List[np.ndarray] = []
    bkg_chunks: List[np.ndarray] = []
    bkg_w_chunks: List[np.ndarray] = []
    sig_chunks: List[np.ndarray] = []
    sig_w_chunks: List[np.ndarray] = []

    for payload in era_payloads.values():
        d_sc = np.asarray(payload["data_full_scores"], dtype=np.float64)  # type: ignore[index]
        d_m = np.asarray(payload["data_full_mass"], dtype=np.float64)  # type: ignore[index]
        if d_sc.size > 0 and d_m.size > 0:
            sel = (d_sc >= lo) & (d_sc < hi)
            if np.any(sel):
                data_chunks.append(d_m[sel])

        for _proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:  # type: ignore[index]
            m_sc = np.asarray(m_sc, dtype=np.float64)
            m_mass = np.asarray(m_mass, dtype=np.float64)
            m_w = np.asarray(m_w, dtype=np.float64) * mc_weight_scale
            sel = (m_sc >= lo) & (m_sc < hi)
            if np.any(sel):
                bkg_chunks.append(m_mass[sel])
                bkg_w_chunks.append(m_w[sel])

        s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)  # type: ignore[index]
        s_mass = np.asarray(payload["sig_full_mass"], dtype=np.float64)  # type: ignore[index]
        s_w = np.asarray(payload["sig_full_weight"], dtype=np.float64)  # type: ignore[index]
        if s_sc.size > 0 and s_mass.size > 0:
            sel = (s_sc >= lo) & (s_sc < hi)
            if np.any(sel):
                sig_chunks.append(s_mass[sel])
                sig_w_chunks.append(s_w[sel])  # no mc_scale for signal

    h_data, v_data = weighted_hist_and_var(
        concat_chunks(data_chunks), mass_bins, None,
    )
    h_bkg, _ = weighted_hist_and_var(
        concat_chunks(bkg_chunks), mass_bins,
        concat_chunks(bkg_w_chunks) if bkg_w_chunks else None,
    )
    h_sig, _ = weighted_hist_and_var(
        concat_chunks(sig_chunks), mass_bins,
        concat_chunks(sig_w_chunks) if sig_w_chunks else None,
    )
    return h_data, v_data, h_bkg, h_sig


def _hwhm_window_from_sigma(mean: float, sigma: float) -> Tuple[float, float]:
    """Return (mass_lo, mass_hi) of the HWHM window for a Gaussian-core signal.

    HWHM = 1.1774 * sigma  (exact half-width at half-maximum of a Gaussian).
    The analysis note (Fig. 34 caption) defines S and B as the integrals of
    signal and background within this window around the fitted signal mean.
    """
    hwhm = 1.1774 * sigma
    return mean - hwhm, mean + hwhm


def _sig_bkg_in_hwhm(
    era_payloads: Dict[str, Dict[str, object]],
    lo: float,
    hi: float,
    mass_lo: float,
    mass_hi: float,
    mc_weight_scale: float,
) -> Tuple[float, float]:
    """Return (s_yield, b_yield) for one BDT category integrated within a mass window."""
    s_total = 0.0
    b_total = 0.0
    for payload in era_payloads.values():
        # Signal
        s_sc   = np.asarray(payload["sig_full_scores"],  dtype=np.float64)  # type: ignore[index]
        s_mass = np.asarray(payload["sig_full_mass"],    dtype=np.float64)  # type: ignore[index]
        s_w    = np.asarray(payload["sig_full_weight"],  dtype=np.float64)  # type: ignore[index]
        if s_sc.size > 0 and s_mass.size > 0 and s_w.size > 0:
            sel = (s_sc >= lo) & (s_sc < hi) & (s_mass >= mass_lo) & (s_mass <= mass_hi)
            s_total += float(np.sum(s_w[sel]))

        # Background (all MC processes)
        for _proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:  # type: ignore[index]
            m_sc   = np.asarray(m_sc,   dtype=np.float64)
            m_mass = np.asarray(m_mass, dtype=np.float64)
            m_w    = np.asarray(m_w,    dtype=np.float64) * float(mc_weight_scale)
            sel = (m_sc >= lo) & (m_sc < hi) & (m_mass >= mass_lo) & (m_mass <= mass_hi)
            b_total += float(np.sum(m_w[sel]))

    return s_total, b_total


def trim_degenerate_top_categories(
    edges: np.ndarray,
    era_payloads: Dict[str, Dict[str, object]],
    args: argparse.Namespace,
    mc_weight_scale: float,
    logger: logging.Logger,
) -> np.ndarray:
    """Merge the top category into the one below if its fit is degenerate."""
    edges = np.asarray(edges, dtype=np.float64).copy()
    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])

    while edges.size > 2:
        lo = float(edges[-2])
        hi = float(edges[-1])

        h_data, v_data, h_bkg, h_sig = _gather_category_data(
            era_payloads, lo, hi, mass_bins, mc_weight_scale,
        )

        if float(np.sum(h_data)) == 0.0 and float(np.sum(h_bkg)) == 0.0:
            logger.info(
                "Top category [%.3f, %.3f) is empty; merging with category below.",
                lo, hi,
            )
            edges = np.concatenate([edges[:-2], [edges[-1]]])
            continue

        h_bkg_fit, best_model, best_chi2ndf, _, _ = fit_background_models_full(
            centers, h_data, v_data, args,
        )
        h_sig_dcb = fit_signal_dcb_template(centers, mass_bins, h_sig, args)

        b_total = float(np.sum(h_bkg_fit))
        s_total = float(np.sum(h_sig_dcb))
        signif = (s_total / math.sqrt(b_total)) if b_total > 0.0 else 0.0

        is_degenerate = (
            best_model == "fallback_data"
            or not np.isfinite(best_chi2ndf)
            or signif == 0.0
        )

        if is_degenerate:
            logger.info(
                "Top category [%.3f, %.3f) is degenerate "
                "(model=%s, chi2/ndf=%s, S/sqrtB=%.4f); merging with category below.",
                lo, hi, best_model,
                f"{best_chi2ndf:.3f}" if np.isfinite(best_chi2ndf) else "inf",
                signif,
            )
            edges = np.concatenate([edges[:-2], [edges[-1]]])
        else:
            break

    return edges


def optimize_category_edges(
    args: argparse.Namespace,
    fallback_eras: Sequence[str],
    logger: logging.Logger,
) -> Tuple[np.ndarray, List[Tuple[int, float]]]:
    """Run greedy category boundary optimisation.

    Returns
    -------
    edges : np.ndarray  — optimised (or user-provided) category edges
    scan_history : list of (n_categories, Z) pairs from the greedy scan.
        Empty when edges are provided manually.
    """
    manual = parse_category_edges_arg(args.category_edges)
    if manual is not None:
        logger.info("Using user-provided category edges: %s", ", ".join([f"{x:.6f}" for x in manual]))
        return manual, []

    opt_eras = args.category_eras if args.category_eras else list(fallback_eras)
    if not opt_eras:
        opt_eras = infer_eras_from_outputs(args.outputs_dir, None)

    sig_scores_chunks: List[np.ndarray] = []
    sig_weight_chunks: List[np.ndarray] = []
    sig_mass_chunks:   List[np.ndarray] = []

    bkg_scores_chunks: List[np.ndarray] = []
    bkg_weight_chunks: List[np.ndarray] = []
    bkg_mass_chunks:   List[np.ndarray] = []

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

        has_mass = mass_branch is not None

        sel_bdt = f"TMath::Finite({args.bdt_branch})"
        if has_mass:
            sel_full = (
                f"({sel_bdt}) && ({mass_branch}>={args.mass_min}"
                f" && {mass_branch}<={args.mass_max})"
            )
            if args.extra_cut:
                sel_full = f"({sel_full}) && ({args.extra_cut})"
            load_branches_s = [args.bdt_branch, args.weight_branch, mass_branch]
            load_branches_b = [args.bdt_branch, args.weight_branch, mass_branch]
        else:
            sel_full = sel_bdt
            if args.extra_cut:
                sel_full = f"({sel_full}) && ({args.extra_cut})"
            load_branches_s = [args.bdt_branch, args.weight_branch]
            load_branches_b = [args.bdt_branch, args.weight_branch]

        try:
            s_arr = load_arrays(sig_file, args.treename, load_branches_s, sel_full)
            b_arr = load_arrays(bkg_file, args.treename, load_branches_b, sel_full)
        except Exception as exc:
            logger.warning("%s: failed loading optimization arrays: %s", era, exc)
            continue

        s_sc = s_arr[args.bdt_branch]
        s_w = sanitize_weights(s_arr[args.weight_branch], args.abs_weights)

        b_sc = b_arr[args.bdt_branch]
        b_w = sanitize_weights(b_arr[args.weight_branch], args.abs_weights)

        if s_sc.size == 0 or b_sc.size == 0:
            continue

        sig_scores_chunks.append(s_sc)
        sig_weight_chunks.append(s_w)

        bkg_scores_chunks.append(b_sc)
        bkg_weight_chunks.append(b_w)

        if has_mass:
            sig_mass_chunks.append(np.asarray(s_arr[mass_branch], dtype=np.float64))
            bkg_mass_chunks.append(np.asarray(b_arr[mass_branch], dtype=np.float64))

        used_eras.append(era)

    if not used_eras:
        raise SystemExit(
            "[error] No scored signal/background inputs found for category optimization. "
            "Provide --category-edges or make sure *_scored_signal.root and *_scored_background.root exist."
        )

    sig_scores = concat_chunks(sig_scores_chunks)
    sig_weights = concat_chunks(sig_weight_chunks)

    bkg_scores = concat_chunks(bkg_scores_chunks)
    bkg_weights = concat_chunks(bkg_weight_chunks)

    total_sig_yield = float(np.sum(sig_weights))
    total_bkg_yield = float(np.sum(bkg_weights))

    # Determine HWHM window from inclusive DCB fit
    hwhm_mass_lo = 0.0
    hwhm_mass_hi = float("inf")
    sig_mass_all: Optional[np.ndarray] = None
    bkg_mass_all: Optional[np.ndarray] = None

    if sig_mass_chunks and bkg_mass_chunks:
        sig_mass_all = concat_chunks(sig_mass_chunks)
        bkg_mass_all = concat_chunks(bkg_mass_chunks)

        mass_bins_opt = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
        centers_opt = 0.5 * (mass_bins_opt[:-1] + mass_bins_opt[1:])
        h_sig_incl, _ = weighted_hist_and_var(
            sig_mass_all, mass_bins_opt,
            sig_weights if sig_weights.size > 0 else None,
        )
        _pred_opt, _par_opt, _chi2_opt = _fit_signal_dcb_full(
            centers_opt, mass_bins_opt, h_sig_incl, args, cat_idx=9999,
        )
        if _par_opt:
            hwhm_mass_lo, hwhm_mass_hi = _hwhm_window_from_sigma(
                _par_opt["mean"], _par_opt["sigma"]
            )
            logger.info(
                "Optimization HWHM window from inclusive DCB fit: "
                "mean=%.3f GeV, sigma=%.3f GeV -> [%.3f, %.3f] GeV",
                _par_opt["mean"], _par_opt["sigma"],
                hwhm_mass_lo, hwhm_mass_hi,
            )
        else:
            hwhm_mass_lo, hwhm_mass_hi = 123.09, 127.09
            logger.warning(
                "Inclusive DCB fit failed; using fallback HWHM window [%.2f, %.2f] GeV",
                hwhm_mass_lo, hwhm_mass_hi,
            )
    else:
        logger.warning(
            "Mass branch not available in scored files; category optimization "
            "will use full-window yields instead of HWHM-windowed yields."
        )

    adaptive_min_sig = min(args.multibin_min_sig,
                           max(1e-6, total_sig_yield * 0.05))
    adaptive_min_bkg = min(args.multibin_min_bkg,
                           max(1e-6, total_bkg_yield * 0.05))

    logger.info(
        "Category optimization eras: %s | S events=%d (yield=%.3f), B events=%d (yield=%.3f)",
        ", ".join(used_eras),
        sig_scores.size,
        total_sig_yield,
        bkg_scores.size,
        total_bkg_yield,
    )
    logger.info(
        "Per-bin minimums: min_sig=%.4f (user=%.4f), min_bkg=%.4f (user=%.4f)",
        adaptive_min_sig, args.multibin_min_sig,
        adaptive_min_bkg, args.multibin_min_bkg,
    )

    edges, scan_history = iterative_greedy_category_scan(
        sig_scores,
        sig_weights,
        bkg_scores,
        bkg_weights,
        max_n=args.multibin_max_n,
        min_improvement=args.multibin_improvement,
        min_sig=adaptive_min_sig,
        min_bkg=adaptive_min_bkg,
        sig_mass=sig_mass_all,
        bkg_mass=bkg_mass_all,
        hwhm_mass_lo=hwhm_mass_lo,
        hwhm_mass_hi=hwhm_mass_hi,
        logger=logger,
    )

    logger.info(
        "Greedy scan result: n=%d categories, edges: %s",
        len(edges) - 1,
        ", ".join([f"{x:.6f}" for x in edges]),
    )
    return edges, scan_history


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

    edges, scan_history = optimize_category_edges(args, eras, logger)

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
        n_data_side = int(run3_data_side_mass.size)
        sum_mc_side = sum(float(np.sum(w)) for _p, _v, w in run3_mc_score_side)
        logger.info(
            "MC normalization (sidebands): scale factor = %.6f  "
            "(N_data_sideband=%d, sum_MC_sideband_weights=%.1f)",
            mc_scale, n_data_side, sum_mc_side,
        )
        if mc_scale < 0.1 or mc_scale > 10.0:
            logger.warning(
                "mc_scale=%.4f is far from 1.  This means the raw MC sideband weights "
                "sum to %.1f events vs %d data sideband events.  The MC weights may not "
                "be normalized to the correct luminosity; S/sqrt(B) values in the "
                "summary table will be unreliable until the normalization is fixed.",
                mc_scale, sum_mc_side, n_data_side,
            )
            # Per-process breakdown: compare sideband vs blinded signal window density
            # to identify any process peaking near 125 GeV (e.g. ggH included as background).
            sb_width = (args.sb1_max - args.sb1_min) + (args.sb2_max - args.sb2_min)
            sigwin_width = args.blind_max - args.blind_min
            proc_b_side: Dict[str, float] = {}
            proc_b_sigwin: Dict[str, float] = {}
            for payload in era_payloads.values():
                for proc, _sc, m_mass, m_w in payload["mc_full_groups"]:  # type: ignore[index]
                    cproc = canonical_proc(proc)
                    m_mass = np.asarray(m_mass, dtype=np.float64)
                    m_w_sc = np.asarray(m_w, dtype=np.float64) * mc_scale
                    sb_mask = sideband_mass_mask(m_mass, args)
                    sigwin_mask = (m_mass >= args.blind_min) & (m_mass <= args.blind_max)
                    proc_b_side[cproc] = proc_b_side.get(cproc, 0.0) + float(np.sum(m_w_sc[sb_mask]))
                    proc_b_sigwin[cproc] = proc_b_sigwin.get(cproc, 0.0) + float(np.sum(m_w_sc[sigwin_mask]))
            logger.warning(
                "Per-process B breakdown (mc_scale corrected, scaled to data):  "
                "%-12s %10s %12s %12s",
                "Process", "B_sideband", "B_signalwin", "Ratio/GeV",
            )
            all_procs = sorted(set(list(proc_b_side) + list(proc_b_sigwin)))
            for cproc in all_procs:
                b_sb = proc_b_side.get(cproc, 0.0)
                b_sw = proc_b_sigwin.get(cproc, 0.0)
                dens_sb = b_sb / sb_width if sb_width > 0 else 0.0
                dens_sw = b_sw / sigwin_width if sigwin_width > 0 else 0.0
                ratio = dens_sw / dens_sb if dens_sb > 1e-9 else float("inf")
                logger.warning("  %-12s %10.2f %12.2f %12.3f", cproc, b_sb, b_sw, ratio)
            logger.warning(
                "  Ratio/GeV > 1 signals a process that peaks in the signal window "
                "(e.g. ggH->mumu); such processes should not be in the background MC."
            )
        # Log total signal yield across all eras for cross-check
        total_sig_w = sum(
            float(np.sum(np.asarray(p["sig_full_weight"], dtype=np.float64)))
            for p in era_payloads.values()
        )
        logger.info(
            "Signal yield check: sum(sig_full_weight) over all eras = %.3f events "
            "(cross-check: VBF H->mumu at full Run3 lumi expect O(10-100) events depending on selection)",
            total_sig_w,
        )

    # Trim degenerate top categories
    edges = trim_degenerate_top_categories(edges, era_payloads, args, mc_scale, logger)

    interior = np.asarray(edges[1:-1], dtype=np.float64)
    category_band = None
    if interior.size >= 1:
        category_band = (float(np.min(interior)), float(np.max(interior)))

    run3_mc_score_full = scale_mc_groups(run3_mc_score_full, mc_scale)
    run3_mc_mass_full = scale_mc_groups(run3_mc_mass_full, mc_scale)
    run3_mc_score_display = group_backgrounds_for_display(run3_mc_score_full)

    # ── Page 2: BDT score figure ────────────────────────────────────────────
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
        sig_scores_all, sig_weights_all = score_signal_overlay
        # Signal normalization (xsec × lumi × BR × eff) is physically independent
        # of the data-driven sideband MC scale.  Do NOT apply mc_scale here.
        sig_weights_overlay = np.asarray(sig_weights_all, dtype=np.float64)
        if sig_scores_all.size > 0 and np.sum(sig_weights_overlay) > 0:
            signal_hist, _ = weighted_hist_and_var(
                np.asarray(sig_scores_all, dtype=np.float64),
                score_bins,
                sig_weights_overlay,
            )

    bdt_cols, bdt_rows = build_bdt_signal_background_table(
        score_bins,
        signal_hist,
        mc_total,
        mc_total_var,
        data_hist,
    )

    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)

    # ── Page 4: Significance evolution (AN Fig. 33 analog) ──────────────────
    fig_sig_evolution = build_significance_evolution_figure(scan_history)

    # ── Page 5: Atanh BDT + S/√B (AN Fig. 34 analog) ───────────────────────
    fig_vbf_transformed = build_vbf_bdt_transformed_figure(
        era_payloads=era_payloads,
        edges=edges,
        args=args,
        mass_bins=mass_bins,
        mc_weight_scale=mc_scale,
    )

    # ── Page 6: Blinded mass sideband ───────────────────────────────────────
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

    # ── Page 7: Per-category mass + background model fits (AN Fig. 36) ──────
    fig_cat_mass, fit_cols, fit_rows = build_category_mass_grid(
        era_payloads=era_payloads,
        edges=edges,
        args=args,
        mass_bins=mass_bins,
        mc_weight_scale=mc_scale,
    )

    # ── Category summary (needs DCB fits per category) ───────────────────────
    cat_summary = summarize_categories(era_payloads, edges, args, mc_weight_scale=mc_scale)

    # ── Page 9: VBF DCB signal extraction (AN Fig. 38) ──────────────────────
    fig_vbf_bdt, dcb_cols, dcb_rows = build_vbf_signal_dcb_grid(
        era_payloads=era_payloads,
        edges=edges,
        args=args,
        mass_bins=mass_bins,
        mc_weight_scale=mc_scale,
        cat_summary=cat_summary,
    )

    # ── Additional: FWHM vs BDT score (AN Fig. 22 analog) ───────────────────
    fig_fwhm_vs_bdt = build_fwhm_vs_bdt_figure(
        era_payloads=era_payloads,
        args=args,
        mass_bins=mass_bins,
    )

    cat_cols, cat_rows = build_category_summary_table(edges, cat_summary)
    cat_proc_cols, cat_proc_rows = build_category_process_table(edges, cat_summary)

    # ── Assemble PDF ─────────────────────────────────────────────────────────
    report_name = args.report_name.strip() if args.report_name else None
    if not report_name:
        era_label = "_".join(sorted(era_payloads.keys())) if era_payloads else "run3"
        report_name = f"{args.output_prefix}_report_{era_label}.pdf"
    report_path = report_name if os.path.isabs(report_name) else os.path.join(args.outdir, report_name)

    with PdfPages(report_path) as pdf:
        # ── 1. Cover ──────────────────────────────────────────────────────────
        add_text_page(
            pdf,
            "Run-3 VBF BDT Sideband Validation Report",
            [
                f"Outputs directory : {args.outputs_dir}",
                f"Eras used         : {', '.join(sorted(era_payloads.keys()))}",
                f"Extra cut         : {args.extra_cut if args.extra_cut else '(none)'}",
                f"Blinding window   : {args.blind_min:.0f} < m_mumu < {args.blind_max:.0f} GeV",
                f"Sidebands         : ({args.sb1_min:.0f}, {args.sb1_max:.0f}) U "
                f"({args.sb2_min:.0f}, {args.sb2_max:.0f}) GeV",
                f"Optimized edges   : {', '.join([f'{x:.4f}' for x in edges])}",
                "",
                "Report contents:",
                "  p.2  BDT score: Data vs stacked MC (full mass window)",
                "  p.3  BDT bin table: signal / background / sensitivity",
                "  p.4  Significance evolution vs. n categories  [AN Fig. 33]",
                "  p.5  Atanh-transformed BDT + S/sqrtB per bin  [AN Fig. 34]",
                "  p.6  FWHM of signal peak vs. VBF BDT score    [AN Fig. 22]",
                "  p.7  Blinded m_mumu: sideband data vs. MC",
                "  p.8  Per-category mass: background model fits  [AN Fig. 36]",
                "  p.9  Background model fit comparison table",
                "  p.10 VBF signal DCB extraction per category   [AN Fig. 38]",
                "  p.11 DCB fit parameter table",
                "  p.12 Optimised category summary",
                "  p.13 Category process composition",
            ],
        )

        # ── 2. BDT score ──────────────────────────────────────────────────────
        pdf.savefig(fig_score)

        # ── 3. BDT bin table ──────────────────────────────────────────────────
        add_table_pages(pdf, "BDT bin summary: signal, background, and sensitivity",
                        bdt_cols, bdt_rows, max_rows=20)

        # ── 4. Significance evolution (AN Fig. 33) ───────────────────────────
        pdf.savefig(fig_sig_evolution)

        # ── 5. Atanh BDT + S/√B (AN Fig. 34) ────────────────────────────────
        pdf.savefig(fig_vbf_transformed)

        # ── 6. FWHM vs BDT score (AN Fig. 22) ────────────────────────────────
        pdf.savefig(fig_fwhm_vs_bdt)

        # ── 7. Blinded mass sideband ──────────────────────────────────────────
        pdf.savefig(fig_mass)

        # ── 8. Per-category mass + bkg fits (AN Fig. 36) ─────────────────────
        pdf.savefig(fig_cat_mass)

        # ── 9. Background model fit comparison table ──────────────────────────
        add_table_pages(pdf, "Per-category background model fit comparison (chi2/ndf, ndf)",
                        fit_cols, fit_rows, max_rows=20)

        # ── 10. VBF DCB signal extraction (AN Fig. 38) ───────────────────────
        pdf.savefig(fig_vbf_bdt)

        # ── 11. DCB parameter table ───────────────────────────────────────────
        add_table_pages(pdf, "VBF signal DCB fit parameters per category",
                        dcb_cols, dcb_rows, max_rows=20)

        # ── 12. Category summary ──────────────────────────────────────────────
        add_table_pages(pdf, "Optimised category summary (HWHM-windowed S and B)",
                        cat_cols, cat_rows, max_rows=20)

        # ── 13. Process composition ────────────────────────────────────────────
        add_table_pages(pdf, "Category composition by process (full mass window)",
                        cat_proc_cols, cat_proc_rows, max_rows=20)

    plt.close(fig_score)
    plt.close(fig_mass)
    plt.close(fig_cat_mass)
    plt.close(fig_vbf_bdt)
    plt.close(fig_vbf_transformed)
    plt.close(fig_sig_evolution)
    plt.close(fig_fwhm_vs_bdt)

    logger.info("Wrote report PDF: %s", report_path)
    print(f"[done] Wrote report: {report_path}")


if __name__ == "__main__":
    main()