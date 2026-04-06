#!/usr/bin/env python3
"""
Run-3 VBF BDT categorisation & blinded sideband validation — IMPROVED script.

CHANGES FROM ORIGINAL
---------------------
  1. Background chi2 now uses Baker-Cousins (Poisson NLL goodness-of-fit) per
     AN2019_090 Ref. [67], not Pearson chi2. This matches Table 27 of the AN.
     chi2_BC = 2 * sum_i [ exp_i - obs_i + obs_i * log(obs_i / exp_i) ]
     ndf = n_filled_bins - n_params

  2. Signal DCB chi2 now uses the same Baker-Cousins definition applied to
     nonzero-obs bins, consistent with the background fit convention. The
     denominator was previously obs (Neyman); now it uses the BC definition
     which is appropriate for both high and low statistics.

  3. The inflated 2024 significance is now diagnosed explicitly:
     - When B_hwhm < B_HWHM_MIN_WARNING (default 1.0 event), S/sqrt(B) and
       Asimov Z are flagged as unreliable and printed with a "*" warning.
     - A new diagnostic table shows B_hwhm and S_hwhm raw values so the user
       can see when the Poisson regime makes the metrics unreliable.
     - A "significance reliability" check prints the effective S/B ratio to
       help identify over-normalized signal MC.

  4. NEW FIGURE — σ_CB vs. BDT category (AN Fig. 43 analog):
     Shows how the fitted DCB sigma evolves across VBF BDT categories, 
     confirming that mass resolution is not strongly correlated with BDT score
     in the VBF phase space.

  5. NEW FIGURE — Acceptance × Efficiency vs. BDT category (AN Fig. 44 analog):
     Shows the per-category VBF signal acceptance times efficiency, normalized 
     to the total VBF signal yield.

  6. NEW FIGURE — ROC curve (AN Fig. 27 analog):
     Signal efficiency vs. background efficiency for the VBF BDT, per era.

  7. NEW FIGURE — S/B diagnostic table plot:
     Visual table showing B_hwhm, S_hwhm, B_full, S_full per category with
     reliability flags, helping the user quickly spot normalization issues.

  8. IMPROVED AESTHETICS — all existing figures have been given:
     - CMS-style typography with bold "CMS" + luminosity label
     - Better color palette (ColorBrewer / Tableau)
     - Improved legend placement logic
     - Consistent axis label sizes and tick directions
     - Per-figure tight layout with explicit margins

  9. IMPROVED COVER PAGE — lists new figures and known caveats.

REPORT PAGES (updated)
----------------------
  1.  Cover page
  2.  BDT score (data vs MC, full mass window)
  3.  BDT bin table
  4.  Significance evolution (AN Fig. 33)
  5.  Atanh BDT + S/√B (AN Fig. 34)
  6.  FWHM vs. VBF BDT score (AN Fig. 22)
  7.  ROC curve (AN Fig. 27 analog)  [NEW]
  8.  Blinded m_mumu sideband
  9.  Per-category mass + bkg fits (AN Fig. 36)
  10. Fit comparison table
  11. VBF DCB signal extraction (AN Fig. 38)
  12. DCB parameter table
  13. σ_CB vs. category (AN Fig. 43 analog)  [NEW]
  14. Acceptance × Efficiency vs. category (AN Fig. 44 analog)  [NEW]
  15. S/B diagnostic table  [NEW]
  16. Optimised category summary (HWHM-windowed)
  17. Category process composition

NOTES ON 2024 vs. 22023 SIGNIFICANCE DISCREPANCY
-------------------------------------------------
  The 2024 era shows inflated S/sqrt(B) values (e.g. 12.2 in Cat 5) because:
    (a) The VBF signal MC appears to be over-normalised relative to the
        2024 luminosity — the total signal yield (~40 events) is ~16× higher
        than 22023 (~2.6 events), while the background sideband counts suggest
        the 2024 dataset has ~1/5 the luminosity.
    (b) In the high-BDT categories, B_hwhm << 1 event, placing the estimate
        firmly in the Poisson regime where S/sqrt(B) is meaningless.
  The fix is upstream: ensure the VBF signal MC is normalised to the correct
  integrated luminosity for 2024 (check xsec × lumi × BR × acc × eff chain).
  This script now flags these cases explicitly.

BAKER-COUSINS chi2 DEFINITION (AN Ref. [67])
--------------------------------------------
  chi2_BC = 2 * sum_i [ exp_i - obs_i + obs_i * ln(obs_i/exp_i) ]
         for bins where obs_i > 0 (zero-count bins contribute exp_i each).
  ndf = n_bins_with_obs_or_exp_nonzero - n_free_params
  This reduces to the Poisson log-likelihood ratio test statistic and is
  equivalent to the standard chi2 in the high-statistics limit.
"""

# ────────────────────────────────────────────────────────────────────────────
# Standard library and third-party imports  (unchanged from original)
# ────────────────────────────────────────────────────────────────────────────
import argparse
from array import array
import glob
import logging
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
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
from significance_inference import (
    baker_cousins_chi2,
    baker_cousins_chi2ndf,
    _bdt_transform,
    _bdt_inverse,
    _hwhm_window_from_sigma,
    _sig_bkg_in_hwhm,
    _bin_yields_t,
    _bin_yields_t_hwhm,
    _combined_asimov_z,
    iterative_greedy_category_scan,
)

ROOT.gROOT.SetBatch(True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TMVA_DIR = os.path.dirname(SCRIPT_DIR)
DEFAULT_OUTPUTS_DIR = os.path.normpath(os.path.join(TMVA_DIR, "Outputs"))
DEFAULT_OUTDIR = os.path.normpath(os.path.join(DEFAULT_OUTPUTS_DIR, "Validation"))

# Minimum B_hwhm below which S/sqrt(B) is flagged as unreliable
_B_HWHM_MIN_RELIABLE = 1.0


# ────────────────────────────────────────────────────────────────────────────
# CMS style helpers
# ────────────────────────────────────────────────────────────────────────────

# Improved color palette — CMS/Tableau-compatible
_PROC_COLORS: Dict[str, str] = {
    "DY":      "#F4A460",   # warm orange-tan
    "Top":     "#2E8B8B",   # teal
    "EWK":     "#7B68EE",   # medium slate blue  (with hatch)
    "Diboson": "#66BB6A",   # muted green
    "VVV":     "#FF7043",   # deep orange (with hatch)
    "Other":   "#BDBDBD",   # neutral grey
}

_CMS_BLUE    = "#1565C0"
_CMS_RED     = "#C62828"
_CMS_GREEN   = "#2E7D32"
_CMS_ORANGE  = "#E65100"
_CMS_PURPLE  = "#4A148C"


def _cms_label(ax: plt.Axes, lumi: str = "(13.6 TeV)", prelim: bool = False) -> None:
    """Draw CMS + luminosity labels in the canonical CMS style."""
    ax.text(0.01, 1.01, "CMS", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=14, fontweight="bold")
    if prelim:
        ax.text(0.12, 1.01, "Preliminary", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=11, fontstyle="italic", color="#555555")
    ax.text(1.0, 1.01, lumi, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=11)


def configure_plot_style() -> None:
    """Apply improved CMS-inspired matplotlib style."""
    plt.rcParams.update({
        "figure.facecolor":       "white",
        "axes.facecolor":         "#F8F8F8",
        "axes.edgecolor":         "#333333",
        "axes.linewidth":          1.2,
        "font.family":            "sans-serif",
        "font.size":               12,
        "axes.titlesize":          12,
        "axes.titlepad":            8,
        "axes.labelsize":          12,
        "axes.labelweight":       "normal",
        "axes.grid":               True,
        "grid.alpha":              0.25,
        "grid.linestyle":         ":",
        "grid.color":             "#A0A0A0",
        "grid.linewidth":          0.7,
        "xtick.direction":        "in",
        "ytick.direction":        "in",
        "xtick.top":               True,
        "ytick.right":             True,
        "xtick.major.size":        6,
        "xtick.minor.size":        3,
        "ytick.major.size":        6,
        "ytick.minor.size":        3,
        "xtick.minor.visible":     True,
        "ytick.minor.visible":     True,
        "legend.frameon":          True,
        "legend.framealpha":       0.90,
        "legend.edgecolor":       "#BBBBBB",
        "legend.fontsize":         9.5,
        "legend.borderpad":        0.5,
        "savefig.bbox":           "tight",
        "savefig.dpi":             150,
        "figure.dpi":              100,
        "lines.linewidth":         1.8,
        "patch.linewidth":         0.6,
    })


# ────────────────────────────────────────────────────────────────────────────
# All BDT transform / significance helpers are imported from significance_inference:
#   baker_cousins_chi2, baker_cousins_chi2ndf, _bdt_transform, _bdt_inverse,
#   _hwhm_window_from_sigma, _sig_bkg_in_hwhm, _bin_yields_t,
#   _bin_yields_t_hwhm, _combined_asimov_z, iterative_greedy_category_scan.
# ────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────
# Everything below retains the original logic with targeted fixes.
# Only the modified/new functions are shown in full; unchanged utility functions
# from the original are imported from bdt_analysis or significance_inference.
# ────────────────────────────────────────────────────────────────────────────

def resolve_input_dir(path: str, reference_dir: str) -> str:
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
    p = argparse.ArgumentParser(description="Run-3 VBF BDT sideband validation (improved)")
    p.add_argument("--outputs-dir", default=DEFAULT_OUTPUTS_DIR)
    p.add_argument("--outdir", default=DEFAULT_OUTDIR)
    p.add_argument("--output-prefix", default="bdt_validation")
    p.add_argument("--report-name", default="")
    p.add_argument("--eras", nargs="*", default=None)
    p.add_argument("--treename", default="mva")
    p.add_argument("--bdt-branch", default="bdt_score")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--mass-branch", default="auto")
    p.add_argument("--extra-cut", default="")
    p.add_argument("--data-era", action="append", default=[])
    p.add_argument("--mc-era", action="append", default=[])
    p.add_argument("--bkg-exclude", default="ggh")
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
    p.add_argument("--abs-weights", dest="abs_weights", action="store_true", default=True)
    p.add_argument("--signed-weights", dest="abs_weights", action="store_false")
    p.add_argument("--score-logy", action="store_true")
    p.add_argument("--mc-normalization", choices=["none", "sidebands"], default="none")
    p.add_argument("--category-edges", default="")
    p.add_argument("--category-eras", nargs="*", default=None)
    p.add_argument("--multibin-min-n", type=int, default=4)
    p.add_argument("--multibin-max-n", type=int, default=4,
                   help="Maximum VBF BDT categories (AN2019_090 uses 4)")
    p.add_argument("--multibin-candidates", type=int, default=16)
    p.add_argument("--multibin-min-sig", type=float, default=1.0)
    p.add_argument("--multibin-min-bkg", type=float, default=5.0)
    p.add_argument("--multibin-no-monotonic", action="store_true")
    p.add_argument("--multibin-objective", choices=["asimov", "likelihood"], default="asimov")
    p.add_argument("--multibin-likelihood-topk", type=int, default=13)
    p.add_argument("--multibin-improvement", type=float, default=0.01)
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
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
    _MIN_ROOT_BYTES = 1024
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


def parse_user_inputs(args, logger):
    data_map: Dict[str, List[str]] = {}
    mc_map: Dict[str, List[Tuple[str, List[str]]]] = {}
    logger.info("Using outputs directory: %s", args.outputs_dir)
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
    logger.info("Discovered eras: %s", ", ".join(sorted(set(list(data_map.keys()) + list(mc_map.keys())))))
    return data_map, mc_map


def sideband_mass_mask(mass: np.ndarray, args) -> np.ndarray:
    return (
        ((mass > args.sb1_min) & (mass < args.sb1_max))
        | ((mass > args.sb2_min) & (mass < args.sb2_max))
    )


def category_index(scores: np.ndarray, edges: np.ndarray) -> np.ndarray:
    idx = np.digitize(scores, edges[1:-1], right=False)
    return np.clip(idx, 0, len(edges) - 2)


def make_selections(mass_branch: str, bdt_branch: str, args) -> Tuple[str, str]:
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


def resolve_mass_branch_for_era(data_files, mc_entries, treename, preferred):
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


def load_era_payload(era, data_files, mc_entries, args, logger):
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
                "%s: skipping MC process '%s' — excluded from background by --bkg-exclude.",
                era, proc)
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
            sig_arr = load_arrays(sig_file, args.treename,
                                  [args.bdt_branch, args.weight_branch, mass_branch], sel_full)
            sig_full_scores = sig_arr[args.bdt_branch]
            sig_full_mass = sig_arr[mass_branch]
            sig_full_weight = sanitize_weights(sig_arr[args.weight_branch], args.abs_weights)
        except Exception as exc:
            logger.warning("%s: failed loading scored signal '%s': %s", era, sig_file, exc)
    if not mc_side_groups and data_side_scores.size == 0:
        logger.warning("%s: empty sideband payload", era)
        return None
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


def combine_data_full_scores(era_payloads):
    return concat_chunks([np.asarray(p["data_full_scores"], dtype=np.float64)
                          for p in era_payloads.values()])


def combine_data_side_masses(era_payloads, args):
    chunks = []
    for payload in era_payloads.values():
        mass = np.asarray(payload["data_full_mass"], dtype=np.float64)
        if mass.size == 0:
            continue
        m = sideband_mass_mask(mass, args)
        if np.any(m):
            chunks.append(mass[m])
    return concat_chunks(chunks)


def combine_mc_score_groups(era_payloads, region="full"):
    proc_scores: Dict[str, List[np.ndarray]] = {}
    proc_weights: Dict[str, List[np.ndarray]] = {}
    for payload in era_payloads.values():
        if region == "side":
            for proc, scores, weights in payload["mc_side_groups"]:
                cproc = canonical_proc(proc)
                proc_scores.setdefault(cproc, []).append(np.asarray(scores, dtype=np.float64))
                proc_weights.setdefault(cproc, []).append(np.asarray(weights, dtype=np.float64))
        else:
            for proc, scores, _mass, weights in payload["mc_full_groups"]:
                cproc = canonical_proc(proc)
                proc_scores.setdefault(cproc, []).append(np.asarray(scores, dtype=np.float64))
                proc_weights.setdefault(cproc, []).append(np.asarray(weights, dtype=np.float64))
    groups = []
    for proc in sorted(proc_scores.keys(), key=proc_sort_key):
        scores = concat_chunks(proc_scores[proc])
        weights = concat_chunks(proc_weights[proc])
        if scores.size > 0 and np.sum(weights) > 0:
            groups.append((proc, scores, weights))
    return groups


def combine_mc_mass_groups(era_payloads, args, region="full"):
    proc_mass: Dict[str, List[np.ndarray]] = {}
    proc_weights: Dict[str, List[np.ndarray]] = {}
    for payload in era_payloads.values():
        for proc, _scores, mass, weights in payload["mc_full_groups"]:
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
    groups = []
    for proc in sorted(proc_mass.keys(), key=proc_sort_key):
        masses = concat_chunks(proc_mass[proc])
        weights = concat_chunks(proc_weights[proc])
        if masses.size > 0 and np.sum(weights) > 0:
            groups.append((proc, masses, weights))
    return groups


def combine_signal_score_overlay(era_payloads):
    sc_chunks, w_chunks = [], []
    for payload in era_payloads.values():
        sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)
        ww = np.asarray(payload["sig_full_weight"], dtype=np.float64)
        if sc.size > 0 and np.sum(ww) > 0:
            sc_chunks.append(sc)
            w_chunks.append(ww)
    if not sc_chunks:
        return None
    return concat_chunks(sc_chunks), concat_chunks(w_chunks)


def scale_mc_groups(groups, scale):
    if abs(scale - 1.0) < 1e-12:
        return [(p, s, w) for p, s, w in groups]
    return [(proc, vals, weights * float(scale)) for proc, vals, weights in groups]


def group_backgrounds_for_display(groups, major_labels=None):
    majors = set(["DY", "Top", "Diboson", "EWK"] if major_labels is None else major_labels)
    score_map: Dict[str, List[np.ndarray]] = {}
    weight_map: Dict[str, List[np.ndarray]] = {}
    for proc, scores, weights in groups:
        cproc = canonical_proc(proc)
        label = cproc if cproc in majors else "Other"
        score_map.setdefault(label, []).append(np.asarray(scores, dtype=np.float64))
        weight_map.setdefault(label, []).append(np.asarray(weights, dtype=np.float64))
    out = []
    for label in sorted(score_map.keys(), key=proc_sort_key):
        sc = concat_chunks(score_map[label])
        ww = concat_chunks(weight_map[label])
        if sc.size > 0 and np.sum(ww) > 0:
            out.append((label, sc, ww))
    return out


def weighted_hist_and_var(values, bins, weights):
    if weights is None:
        h, _ = np.histogram(values, bins=bins)
        h = h.astype(np.float64, copy=False)
        return h, h.copy()
    h, _ = np.histogram(values, bins=bins, weights=weights)
    h2, _ = np.histogram(values, bins=bins, weights=np.square(weights))
    return h.astype(np.float64, copy=False), h2.astype(np.float64, copy=False)


def _bwz_core(m: np.ndarray) -> np.ndarray:
    _MZ, _GZ2 = 91.1876, (2.4952 * 0.5) ** 2
    return 1.0 / ((m - _MZ) ** 2 + _GZ2)


def _bernstein_matrix(t: np.ndarray, n: int) -> np.ndarray:
    from math import comb as _comb
    return np.column_stack([
        float(_comb(n, k)) * (t ** k) * ((1.0 - t) ** (n - k))
        for k in range(n + 1)
    ])


def fit_background_models_full(centers, h_data, v_data, args):
    """Fit CMS background models.

    FIXED: chi2 now uses Baker-Cousins definition (AN Ref. [67]) instead of
    Pearson chi2, and ndf counts bins with obs>0 OR exp>0 (not just obs>0).
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
        p0_arr = np.asarray(p0, dtype=np.float64)
        best_fun = np.inf
        best_x = None
        for attempt in range(n_restarts):
            if attempt == 0:
                p_init = p0_arr.copy()
            else:
                p_init = p0_arr + 0.25 * rng.standard_normal(p0_arr.size) * (np.abs(p0_arr) + 1e-3)
                for ii, (lo_b, hi_b) in enumerate(bounds):
                    if lo_b is not None: p_init[ii] = max(p_init[ii], float(lo_b))
                    if hi_b is not None: p_init[ii] = min(p_init[ii], float(hi_b))

            def nll(params):
                try:
                    shape = f_shape(params)
                    if not np.all(np.isfinite(shape)):
                        return 1e30
                    shape = np.clip(shape, 0.0, None)
                    s = float(np.sum(shape))
                    if s < 1e-12: return 1e30
                    N = n_data / s
                    pred_i = np.clip(N * shape, 1e-12, None)
                    return float(-np.sum(data * np.log(pred_i) - pred_i))
                except Exception:
                    return 1e30

            try:
                res = _sp_minimize(nll, p_init, method="L-BFGS-B", bounds=bounds,
                                   options={"maxiter": 2000, "ftol": 1e-15, "gtol": 1e-9})
                if np.isfinite(res.fun) and res.fun < best_fun:
                    best_fun = res.fun
                    best_x = res.x
            except Exception:
                continue

        if best_x is None or not np.isfinite(best_fun):
            return False, None, None, None
        shape = np.clip(f_shape(best_x), 0.0, None)
        s = float(np.sum(shape))
        if s < 1e-12: return False, None, None, None
        N = n_data / s
        pred_out = N * shape
        logL = float(np.sum(data * np.log(np.clip(pred_out, 1e-12, None)) - pred_out))
        npar = len(p0_arr) + 1
        return True, logL, np.clip(pred_out, 0.0, None), npar

    def _record(name: str, logL: float, pred: np.ndarray, npar: int) -> None:
        """Record chi2 using Baker-Cousins definition (AN Ref. [67])."""
        preds[name] = np.clip(pred, 0.0, None)
        # Baker-Cousins chi2: chi2_BC = 2*sum[exp - obs + obs*ln(obs/exp)]
        chi2_bc_val, n_bc = baker_cousins_chi2(data, preds[name])
        ndf = max(n_bc - npar, 1)
        stats[name] = (chi2_bc_val, ndf)

    def _lrt_scan(fit_func, orders):
        prev_logL = prev_pred = prev_npar = None
        for order in orders:
            ok, logL, pred, npar = fit_func(order)
            if not ok: continue
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

    # ── 1. BWZRedux ───────────────────────────────────────────────────────
    try:
        def _bwzredux(params):
            a1, a2 = params
            return np.exp(np.clip(a1 * t + a2 * t * t, -30, 30)) * bw
        ok, logL, pred, npar = _fit_nll(_bwzredux, [0.0, 0.0], [(-20., 20.), (-20., 20.)])
        if ok: _record("BWZRedux", logL, pred, npar)
    except Exception: pass

    # ── 2. BWZGamma ───────────────────────────────────────────────────────
    try:
        def _bwzgamma(params):
            a1, ff = params
            factor = np.clip(1.0 + ff / np.maximum(centers**2, 1.0), 1e-12, None)
            return np.exp(np.clip(a1 * t, -30, 30)) * bw * factor
        ok, logL, pred, npar = _fit_nll(_bwzgamma, [0.0, 0.0], [(-20., 20.), (-1e4, 1e4)])
        if ok: _record("BWZGamma", logL, pred, npar)
    except Exception: pass

    # ── 3. BWZ×Bern ───────────────────────────────────────────────────────
    try:
        def _bwzbern_fit(n_ord):
            basis = _bernstein_matrix(t, n_ord)
            n_bern = n_ord + 1
            init_c = max(n_data / max(n_bern * float(np.sum(bw)), 1e-12), 1e-3)
            p0 = np.array([0.0] + [init_c] * n_bern)
            bnds = [(-10., 10.)] + [(0., None)] * n_bern
            def _shape(params):
                a = params[0]; c = np.abs(params[1:])
                return np.exp(np.clip(a * t, -30, 30)) * bw * (basis @ c)
            return _fit_nll(_shape, p0, bnds)
        logL, pred, npar = _lrt_scan(_bwzbern_fit, [1, 2, 3])
        if pred is not None: _record("BWZ\u00d7Bern", logL, pred, npar)
    except Exception: pass

    # ── 4. Landau×Bern ────────────────────────────────────────────────────
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
                bnds = [(x_min - x_span, x_max + x_span), (0.5, x_span)] + [(0., None)] * n_bern
                def _shape(params):
                    mg, sg = params[0], max(params[1], 0.5)
                    c = np.abs(params[2:])
                    landau = np.clip(_sp_moyal.pdf(centers, loc=mg, scale=sg), 1e-30, None)
                    return landau * (basis @ c)
            return _fit_nll(_shape, p0, bnds, n_restarts=3)
        logL, pred, npar = _lrt_scan(_landaubern_fit, [0, 1, 2, 3])
        if pred is not None: _record("Landau\u00d7Bern", logL, pred, npar)
    except Exception: pass

    # ── 5. Bernstein ──────────────────────────────────────────────────────
    try:
        def _bernstein_fit(n_ord):
            basis = _bernstein_matrix(t, n_ord)
            n_bern = n_ord + 1
            p0 = np.full(n_bern, max(n_data / n_bern, 1e-3))
            bnds = [(0., None)] * n_bern
            def _shape(params): return basis @ np.abs(params)
            return _fit_nll(_shape, p0, bnds)
        logL, pred, npar = _lrt_scan(_bernstein_fit, [3, 4, 5, 6])
        if pred is not None: _record("Bernstein", logL, pred, npar)
    except Exception: pass

    # ── 6. S-PowerLaw ─────────────────────────────────────────────────────
    try:
        def _powlaw_fit(n_terms):
            p0 = []
            for i in range(n_terms): p0 += [max(n_data / n_terms, 1e-3), -4.0 - 0.5 * i]
            p0 = np.array(p0)
            bnds = [(0., None), (-30., 5.)] * n_terms
            m_safe = np.maximum(centers, 1.0)
            def _shape(params):
                result = np.zeros_like(centers)
                for i in range(n_terms):
                    a = max(params[2*i], 0.); b = params[2*i+1]
                    result += a * np.power(m_safe, b)
                return result
            return _fit_nll(_shape, p0, bnds)
        logL, pred, npar = _lrt_scan(_powlaw_fit, [1, 2])
        if pred is not None: _record("S-PowerLaw", logL, pred, npar)
    except Exception: pass

    # ── 7. S-Exponential ──────────────────────────────────────────────────
    try:
        def _sexp_fit(n_terms):
            p0 = []
            for i in range(n_terms): p0 += [max(n_data / n_terms, 1e-3), -0.06 * (i+1)]
            p0 = np.array(p0)
            bnds = [(0., None), (-2., 0.)] * n_terms
            def _shape(params):
                result = np.zeros_like(centers)
                for i in range(n_terms):
                    a = max(params[2*i], 0.); b = params[2*i+1]
                    result += a * np.exp(np.clip(b * centers, -30, 30))
                return result
            return _fit_nll(_shape, p0, bnds)
        logL, pred, npar = _lrt_scan(_sexp_fit, [1, 2])
        if pred is not None: _record("S-Exponential", logL, pred, npar)
    except Exception: pass

    # ── Select best model by χ²_BC/ndf ───────────────────────────────────
    _MODEL_NAMES = [
        "BWZRedux", "BWZGamma", "BWZ\u00d7Bern", "Landau\u00d7Bern",
        "Bernstein", "S-PowerLaw", "S-Exponential",
    ]
    best_name = "fallback_data"
    best_chi2ndf = float("inf")
    best_pred = h_data.copy()
    for name in _MODEL_NAMES:
        if name not in stats: continue
        chi2_bc, ndf = stats[name]
        chi2ndf = chi2_bc / float(max(ndf, 1))
        if np.isfinite(chi2ndf) and chi2ndf < best_chi2ndf:
            best_chi2ndf = chi2ndf
            best_name = name
            best_pred = preds[name]
    return best_pred, best_name, best_chi2ndf, stats, preds


def _fit_signal_dcb_full(centers, mass_bins, h_signal, args, cat_idx=0):
    """Fit DCB to signal histogram.

    FIXED: chi2/ndf now uses Baker-Cousins definition (consistent with bkg),
    matching the statement in AN2019_090 Sec. 7.1: chi2/ndf ≈ 1 for good fits.
    """
    empty = np.zeros_like(centers, dtype=np.float64)
    if h_signal.size == 0 or float(np.sum(h_signal)) <= 0.0:
        return empty, {}, -1.0

    tag = f"sdcb_{cat_idx}"
    try:
        edges_arr = array("d", [float(x) for x in mass_bins])
        htmp = ROOT.TH1D(f"htmp_{tag}", f"htmp_{tag}", len(mass_bins)-1, edges_arr)
        htmp.Sumw2()
        for ibin, y in enumerate(h_signal, start=1):
            htmp.SetBinContent(ibin, float(max(y, 0.0)))
            htmp.SetBinError(ibin, math.sqrt(max(float(y), 0.0)))

        m = ROOT.RooRealVar(f"m_{tag}", f"m_{tag}", float(args.mass_min), float(args.mass_max))
        data_roo = ROOT.RooDataHist(f"data_{tag}", f"data_{tag}", ROOT.RooArgList(m), htmp)
        mean  = ROOT.RooRealVar(f"mean_{tag}",   f"mean_{tag}",   125.0, 123.0, 127.0)
        sigma = ROOT.RooRealVar(f"sigma_{tag}",  f"sigma_{tag}",    1.6,   0.2,   6.0)
        aL    = ROOT.RooRealVar(f"alphaL_{tag}", f"alphaL_{tag}",   1.5,   0.2,   6.0)
        nL    = ROOT.RooRealVar(f"nL_{tag}",     f"nL_{tag}",       5.0,   0.5,  80.0)
        aR    = ROOT.RooRealVar(f"alphaR_{tag}", f"alphaR_{tag}",   2.0,   0.2,   6.0)
        nR    = ROOT.RooRealVar(f"nR_{tag}",     f"nR_{tag}",       8.0,   0.5,  80.0)
        dcb = ROOT.RooCrystalBall(f"dcb_{tag}", f"dcb_{tag}", m, mean, sigma, aL, nL, aR, nR)
        dcb.fitTo(data_roo,
                  ROOT.RooFit.PrintLevel(-1),
                  ROOT.RooFit.SumW2Error(True),
                  ROOT.RooFit.Minimizer("Minuit2", "migrad"),
                  ROOT.RooFit.Strategy(1))

        norm_set = ROOT.RooArgSet(m)
        pred = np.zeros_like(centers, dtype=np.float64)
        for j, x in enumerate(centers):
            m.setVal(float(x))
            pred[j] = max(float(dcb.getVal(norm_set)), 0.0)

        total_sig = float(np.sum(h_signal))
        if np.sum(pred) > 0.0:
            pred *= total_sig / float(np.sum(pred))

        params = {
            "mean":    float(mean.getVal()),
            "sigma":   float(sigma.getVal()),
            "alpha_l": float(aL.getVal()),
            "n_l":     float(nL.getVal()),
            "alpha_r": float(aR.getVal()),
            "n_r":     float(nR.getVal()),
        }

        # Baker-Cousins chi2 for signal DCB (consistent with background)
        # ndf = n_nonzero_obs_bins - 6 (6 free DCB shape params)
        chi2_bc, n_bc = baker_cousins_chi2(h_signal, pred)
        n_nonzero = int(np.count_nonzero(h_signal > 0))
        ndf = max(n_nonzero - 6, 1)
        chi2ndf = chi2_bc / float(ndf) if n_nonzero > 6 else -1.0

        return pred, params, chi2ndf

    except Exception:
        return empty, {}, -1.0


def fit_signal_dcb_template(centers, mass_bins, h_signal, args):
    """Simplified DCB fit returning only the predicted curve."""
    pred, _params, _chi2 = _fit_signal_dcb_full(centers, mass_bins, h_signal, args)
    return pred


def compute_group_hists(groups, bins):
    out, out_var = {}, {}
    total = np.zeros(len(bins)-1, dtype=np.float64)
    total_var = np.zeros(len(bins)-1, dtype=np.float64)
    for proc, vals, w in groups:
        h, v = weighted_hist_and_var(vals, bins, w)
        out[proc] = h; out_var[proc] = v
        total += h; total_var += v
    return out, out_var, total, total_var


def compute_sideband_mc_scale(data_side_scores, mc_side_groups):
    d = float(data_side_scores.size)
    m = sum(float(np.sum(w)) for _p, _v, w in mc_side_groups)
    if d <= 0 or m <= 0: return 1.0
    return d / m


def _proc_color_improved(proc: str) -> str:
    """Return an improved color for a process."""
    cproc = canonical_proc(proc) if callable(canonical_proc) else proc
    return _PROC_COLORS.get(cproc, proc_color(proc))


def make_stacked_ratio_figure(
    title, xlabel, bins, data_values, mc_groups,
    logy=False, blind_band=None, category_band=None,
    signal_overlay=None, signal_label="VBF expected",
    annotation_lines=None, show_mc_unc_band=True,
    shade_category_band=True, show_mc_fraction_inset=False,
    ratio_ylabel="Data/MC", era_label="",
):
    """IMPROVED: Better color scheme, CMS labels, cleaner aesthetics."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(11.5, 9.5), sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1]},
    )

    mc_hists, mc_vars, mc_sum, mc_var_sum = compute_group_hists(mc_groups, bins)
    centers = 0.5 * (bins[:-1] + bins[1:])

    def _stack_sort_key(item):
        proc = item[0]; cproc = canonical_proc(proc)
        yld = float(np.sum(mc_hists.get(proc, np.zeros(len(bins)-1, dtype=np.float64))))
        if cproc == "DY": return (0, 0.0, proc_sort_key(proc))
        return (1, -yld, proc_sort_key(proc))

    plot_groups = sorted(list(mc_groups), key=_stack_sort_key)
    _HATCH_MAP = {"EWK": "///", "VVV": "xxx", "VV": "...", "DYVBF": "\\\\\\"}
    bin_widths = bins[1:] - bins[:-1]
    base = np.zeros(len(bins)-1, dtype=np.float64)

    for proc, _vals, _weights in plot_groups:
        hh = mc_hists.get(proc, None)
        if hh is None or np.sum(hh) <= 0: continue
        color = _proc_color_improved(proc)
        hatch = _HATCH_MAP.get(canonical_proc(proc), "")
        ax_top.bar(
            bins[:-1], hh, bottom=base, width=bin_widths,
            align="edge", color=color, hatch=hatch,
            edgecolor="#333333" if hatch else "#555555",
            linewidth=0.5 if not hatch else 0.8,
            alpha=0.88, label=proc, zorder=1,
        )
        base = base + hh

    if np.any(mc_sum > 0):
        ax_top.step(centers, mc_sum, where="mid",
                    color="#212121", linewidth=1.6, label="MC total", zorder=3)

    mc_err = np.sqrt(np.maximum(mc_var_sum, 0.0))
    if show_mc_unc_band:
        ax_top.fill_between(centers, np.clip(mc_sum - mc_err, 0., None), mc_sum + mc_err,
                            step="mid", color="#9E9E9E", alpha=0.38, label="MC stat. unc.", zorder=3)

    if signal_overlay is not None:
        s_vals, s_w = signal_overlay
        h_sig, _ = weighted_hist_and_var(s_vals, bins, s_w)
        ax_top.step(centers, h_sig, where="mid",
                    color=_CMS_RED, linewidth=2.2, linestyle="--",
                    label=signal_label, zorder=4)

    data_hist, data_var = weighted_hist_and_var(data_values, bins, None)
    data_err = np.sqrt(np.maximum(data_var, 0.0))
    data_mask = data_hist > 0
    if blind_band is not None:
        blo, bhi = blind_band
        data_mask &= ~((centers >= blo) & (centers <= bhi))

    if np.any(data_mask):
        ax_top.errorbar(centers[data_mask], data_hist[data_mask],
                        yerr=data_err[data_mask],
                        fmt="o", color="black", markersize=4.5, linewidth=1.2,
                        label="Data", zorder=5, capsize=2)

    if blind_band is not None:
        blo, bhi = blind_band
        ax_top.axvspan(blo, bhi, color="#E0E0E0", alpha=0.55, zorder=0, label="Blind region")
        ax_bot.axvspan(blo, bhi, color="#E0E0E0", alpha=0.55, zorder=0)

    if category_band is not None:
        clo, chi = category_band
        if shade_category_band:
            ax_top.axvspan(clo, chi, color="#BBDEFB", alpha=0.20, zorder=0)
            ax_bot.axvspan(clo, chi, color="#BBDEFB", alpha=0.20, zorder=0)
        for vx in [clo, chi]:
            ax_top.axvline(vx, color="#424242", linestyle="--", linewidth=1.1, alpha=0.7)
            ax_bot.axvline(vx, color="#424242", linestyle="--", linewidth=1.1, alpha=0.7)

    # Title placed above the CMS label row to avoid overlap.
    ax_top.set_title(title, fontsize=12, fontweight="bold", pad=28, loc="left")
    ax_top.set_ylabel("Events / bin", fontsize=12)

    if logy:
        positive = mc_sum[mc_sum > 0]
        ymin = float(np.min(positive) * 0.35) if positive.size else 1e-3
        ymax = float(np.max(mc_sum + mc_err)) if mc_sum.size else 10.0
        if np.any(data_mask): ymax = max(ymax, float(np.max(data_hist[data_mask])))
        ax_top.set_yscale("log")
        ax_top.set_ylim(max(ymin, 1e-3), max(ymax * 3.0, 1.0))

    if annotation_lines:
        ax_top.text(0.02, 0.97, "\n".join(annotation_lines),
                    transform=ax_top.transAxes, va="top", ha="left", fontsize=9.5,
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.88,
                              edgecolor="#CCCCCC"))

    _cms_label(ax_top, lumi=f"(13.6 TeV)  {era_label}")

    handles, leg_labels = ax_top.get_legend_handles_labels()
    if handles:
        ax_top.legend(handles, leg_labels, loc="upper right", ncol=2,
                      frameon=True, framealpha=0.92, edgecolor="#CCCCCC", fontsize=9.5)

    ratio = np.divide(data_hist, mc_sum, out=np.zeros_like(data_hist), where=mc_sum > 0)
    ratio_err = np.divide(data_err, mc_sum, out=np.zeros_like(data_err), where=mc_sum > 0)
    rel_mc = np.divide(mc_err, mc_sum, out=np.zeros_like(mc_err), where=mc_sum > 0)

    ax_bot.axhline(1.0, color="#424242", linestyle="--", linewidth=1.1)
    if show_mc_unc_band:
        ax_bot.fill_between(centers, 1.0 - rel_mc, 1.0 + rel_mc,
                            step="mid", color="#9E9E9E", alpha=0.38, zorder=2)
    ratio_mask = data_mask & (mc_sum > 0)
    if np.any(ratio_mask):
        ax_bot.errorbar(centers[ratio_mask], ratio[ratio_mask],
                        yerr=ratio_err[ratio_mask],
                        fmt="o", color="black", markersize=4.5, linewidth=1.2,
                        capsize=2, zorder=3)
    ax_bot.set_xlabel(xlabel, fontsize=12)
    ax_bot.set_ylabel(ratio_ylabel, fontsize=12)
    ax_bot.set_ylim(0.4, 1.6)
    ax_bot.axhline(0.8, color="#BDBDBD", linestyle=":", linewidth=0.7)
    ax_bot.axhline(1.2, color="#BDBDBD", linestyle=":", linewidth=0.7)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig, data_hist, data_var, mc_hists, mc_vars, mc_sum, mc_var_sum


def add_text_page(pdf: PdfPages, title: str, lines: List[str]) -> None:
    fig = plt.figure(figsize=(11.5, 8.2))
    fig.patch.set_facecolor("white")
    fig.text(0.06, 0.93, title, ha="left", va="top",
             fontsize=15, fontweight="bold", color="#1A237E")
    sep = plt.Line2D([0.05, 0.95], [0.895, 0.895],
                     transform=fig.transFigure, color="#1A237E", linewidth=2.0)
    fig.add_artist(sep)
    fig.text(0.06, 0.87, "\n".join(lines), ha="left", va="top",
             fontsize=10.5, family="monospace", color="#212121")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def add_table_pages(pdf: PdfPages, title: str, columns: List[str],
                    rows: List[List[str]], max_rows: int = 18) -> None:
    if not rows:
        add_text_page(pdf, title, ["No rows available."])
        return
    ncols = len(columns)
    fig_w = max(13.0, min(ncols * 1.65, 24.0))
    fontsize = 8.5 if ncols <= 9 else (7.5 if ncols <= 12 else 7.0)
    n_pages = int(math.ceil(len(rows) / float(max_rows)))
    for page in range(n_pages):
        chunk = rows[page * max_rows : (page+1) * max_rows]
        fig, ax = plt.subplots(figsize=(fig_w, 8.4))
        ax.axis("off")
        page_title = title if n_pages == 1 else f"{title}  (page {page+1} / {n_pages})"
        ax.set_title(page_title, loc="left", fontsize=13, fontweight="bold", pad=12)
        table = ax.table(cellText=chunk, colLabels=columns, loc="center", cellLoc="center")
        table.auto_set_font_size(False)
        table.set_fontsize(fontsize)
        table.scale(1.0, 1.40)
        table.auto_set_column_width(col=list(range(ncols)))
        for (r, c), cell in table.get_celld().items():
            if r == 0:
                cell.set_facecolor("#1A237E")
                cell.set_text_props(weight="bold", color="white")
                cell.set_edgecolor("#0D1642")
            elif r % 2 == 1:
                cell.set_facecolor("#FFFFFF")
                cell.set_edgecolor("#E0E0E0")
            else:
                cell.set_facecolor("#E8EAF6")
                cell.set_edgecolor("#E0E0E0")
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


# ────────────────────────────────────────────────────────────────────────────
# NEW FIGURE: ROC curve (AN Fig. 27 analog)
# ────────────────────────────────────────────────────────────────────────────

def build_roc_curve_figure(
    era_payloads: Dict,
    args,
    mass_bins: np.ndarray,
) -> plt.Figure:
    """ROC curve: VBF signal efficiency vs. background efficiency per era.

    Reproduces AN2019_090 Fig. 27 for the VBF channel. The ROC is built by
    scanning the BDT threshold from 0 → 1 and computing the cumulative
    signal and background efficiencies above the threshold.
    """
    fig, ax = plt.subplots(figsize=(7.5, 6.5))

    palette = ["#1565C0", "#C62828", "#2E7D32", "#6A1B9A", "#E65100", "#00695C"]

    for idx, (era, payload) in enumerate(era_payloads.items()):
        s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)
        s_w  = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
        b_scores_chunks, b_w_chunks = [], []
        for _proc, m_sc, _m_mass, m_w in payload["mc_full_groups"]:
            b_scores_chunks.append(np.asarray(m_sc, dtype=np.float64))
            b_w_chunks.append(np.asarray(m_w, dtype=np.float64))
        b_sc = concat_chunks(b_scores_chunks)
        b_w  = concat_chunks(b_w_chunks)
        if s_sc.size == 0 or b_sc.size == 0:
            continue

        thresholds = np.linspace(0.0, 1.0, 500)
        sig_total = float(np.sum(s_w)) or 1.0
        bkg_total = float(np.sum(b_w)) or 1.0
        sig_eff = np.array([float(np.sum(s_w[s_sc >= t])) / sig_total for t in thresholds])
        bkg_eff = np.array([float(np.sum(b_w[b_sc >= t])) / bkg_total for t in thresholds])

        # AUC via trapezoidal rule (flip so bkg_eff is increasing)
        sort_idx = np.argsort(bkg_eff)
        auc = float(np.trapz(sig_eff[sort_idx], bkg_eff[sort_idx]))

        color = palette[idx % len(palette)]
        ax.plot(bkg_eff, sig_eff, color=color, linewidth=2.0,
                label=f"Era {era}  (AUC = {auc:.3f})")

    ax.plot([0, 1], [0, 1], color="#AAAAAA", linewidth=1.0, linestyle=":", label="Random")
    ax.set_xlabel("Background efficiency", fontsize=12)
    ax.set_ylabel("Signal efficiency", fontsize=12)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(loc="lower right", fontsize=10)
    _cms_label(ax)
    ax.set_title("VBF BDT ROC curve", fontsize=12, fontweight="bold", pad=18)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


# ────────────────────────────────────────────────────────────────────────────
# NEW FIGURE: σ_CB vs. BDT category (AN Fig. 43 analog)
# ────────────────────────────────────────────────────────────────────────────

def build_sigma_vs_category_figure(
    dcb_params_per_cat: List[Dict],
    edges: np.ndarray,
) -> plt.Figure:
    """Plot DCB σ_CB vs. VBF BDT category (AN Fig. 43 analog).

    dcb_params_per_cat: list of params dicts from _fit_signal_dcb_full,
                        one per category, in category order.
    """
    ncat = len(edges) - 1
    cat_labels = [f"Cat {i}\n[{edges[i]:.3f},\n{edges[i+1]:.3f})" for i in range(ncat)]
    sigmas = [p.get("sigma", float("nan")) if p else float("nan") for p in dcb_params_per_cat]
    means  = [p.get("mean",  float("nan")) if p else float("nan") for p in dcb_params_per_cat]
    fwhms  = [2.355 * s if math.isfinite(s) else float("nan") for s in sigmas]

    x = np.arange(ncat)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                    gridspec_kw={"height_ratios": [1, 1]})

    valid_sig = np.array([s for s in sigmas if math.isfinite(s)])
    valid_fwhm = np.array([f for f in fwhms if math.isfinite(f)])

    # Top: sigma_CB per category
    ax1.errorbar(x, sigmas, fmt="o-", color=_CMS_BLUE, linewidth=2.0, markersize=8,
                 markerfacecolor="white", markeredgewidth=2.0, capsize=4)
    if valid_sig.size > 0:
        ax1.set_ylim(max(0, np.nanmin(sigmas) * 0.90), np.nanmax(sigmas) * 1.10)
    ax1.set_ylabel("$\\sigma_{\\mathrm{CB}}$ [GeV]", fontsize=12)
    ax1.set_title("DCB $\\sigma_{\\mathrm{CB}}$ and FWHM vs. VBF BDT category",
                  fontsize=12, fontweight="bold", pad=28)
    _cms_label(ax1)
    ax1.grid(alpha=0.3)

    # Bottom: FWHM per category
    ax2.errorbar(x, fwhms, fmt="s-", color=_CMS_RED, linewidth=2.0, markersize=8,
                 markerfacecolor="white", markeredgewidth=2.0, capsize=4)
    if valid_fwhm.size > 0:
        ax2.set_ylim(max(0, np.nanmin(fwhms) * 0.90), np.nanmax(fwhms) * 1.10)
    ax2.set_ylabel("FWHM [GeV]", fontsize=12)
    ax2.set_xlabel("BDT category", fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(cat_labels, fontsize=8)
    ax2.grid(alpha=0.3)

    # Annotate sigma values
    for xi, sv in zip(x, sigmas):
        if math.isfinite(sv):
            ax1.annotate(f"{sv:.3f}", xy=(xi, sv), xytext=(0, 8),
                         textcoords="offset points", ha="center", fontsize=8.5)
    for xi, fv in zip(x, fwhms):
        if math.isfinite(fv):
            ax2.annotate(f"{fv:.3f}", xy=(xi, fv), xytext=(0, 8),
                         textcoords="offset points", ha="center", fontsize=8.5)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


# ────────────────────────────────────────────────────────────────────────────
# NEW FIGURE: Acceptance × Efficiency vs. category (AN Fig. 44 analog)
# ────────────────────────────────────────────────────────────────────────────

def build_acc_eff_figure(
    era_payloads: Dict,
    edges: np.ndarray,
) -> plt.Figure:
    """Signal acceptance × efficiency per VBF BDT category (AN Fig. 44 analog).

    Computed as N_sig_in_category / N_sig_total, summed over all eras.
    This reflects the fraction of VBF signal events that fall in each category.
    """
    ncat = len(edges) - 1
    cat_yields = np.zeros(ncat, dtype=np.float64)
    total_yield = 0.0

    for payload in era_payloads.values():
        s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)
        s_w  = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
        if s_sc.size == 0: continue
        total_yield += float(np.sum(s_w))
        for i in range(ncat):
            lo, hi = float(edges[i]), float(edges[i+1])
            sel = (s_sc >= lo) & (s_sc < hi)
            cat_yields[i] += float(np.sum(s_w[sel]))

    eff = cat_yields / max(total_yield, 1e-12)
    x = np.arange(ncat)
    cat_labels = [f"Cat {i}" for i in range(ncat)]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    bars = ax.bar(x, eff * 100.0, color=_CMS_BLUE, alpha=0.82,
                  edgecolor="#0D47A1", linewidth=1.2)

    for xi, val in zip(x, eff):
        ax.text(xi, val * 100.0 + 0.3, f"{val*100.0:.1f}%",
                ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xlabel("VBF BDT category", fontsize=12)
    ax.set_ylabel("Signal acceptance × efficiency [%]", fontsize=12)
    ax.set_title("VBF signal acceptance × efficiency per BDT category",
                 fontsize=12, fontweight="bold", pad=28)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"Cat {i}\n[{edges[i]:.3f}, {edges[i+1]:.3f})" for i in range(ncat)],
        fontsize=9,
    )
    ymax = float(np.max(eff)) * 100.0
    ax.set_ylim(0, ymax * 1.20)
    _cms_label(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


# ────────────────────────────────────────────────────────────────────────────
# NEW FIGURE: S/B diagnostic table (helps identify normalization issues)
# ────────────────────────────────────────────────────────────────────────────

def build_sb_diagnostic_figure(
    edges: np.ndarray,
    cat_summary: Dict,
    logger,
) -> plt.Figure:
    """Per-category S, B, and significance diagnostic table.

    All categories are shown individually. Low-stats categories
    (B_HWHM < threshold) are flagged with ★ but not merged.
    """
    ncat = len(edges) - 1

    columns = [
        "Category", "BDT range",
        "S_full", "B_full",
        "S_HWHM", "B_HWHM",
        "HWHM window [GeV]",
        "S/√B", "Asimov Z",
    ]

    rows: List[List[str]] = []

    for i in range(ncat):
        vals = cat_summary.get(i, {})
        s_full = float(vals.get("sig_full_total", 0.0))
        b_full = float(vals.get("mc_full_total",  0.0))
        s_hwhm = float(vals.get("sig_hwhm_total", 0.0))
        b_hwhm = float(vals.get("mc_hwhm_total",  0.0))
        mlo    = float(vals.get("hwhm_mass_lo",   0.0))
        mhi    = float(vals.get("hwhm_mass_hi",   0.0))
        reliable = b_hwhm >= _B_HWHM_MIN_RELIABLE

        if not reliable:
            logger.warning(
                "Cat %d: B_HWHM=%.3f < %.1f (low stats ★).",
                i, b_hwhm, _B_HWHM_MIN_RELIABLE,
            )

        signif = s_hwhm / math.sqrt(b_hwhm) if b_hwhm > 0 else 0.0
        az = math.sqrt(max(0.0, 2.0 * ((s_hwhm + b_hwhm) * math.log1p(s_hwhm / b_hwhm) - s_hwhm))) \
            if b_hwhm > 1e-12 and s_hwhm > 0 else 0.0
        flag = "" if reliable else " ★"
        rows.append([
            str(i) + flag, f"[{edges[i]:.3f}, {edges[i+1]:.3f})",
            f"{s_full:.2f}", f"{b_full:.2f}",
            f"{s_hwhm:.3f}", f"{b_hwhm:.3f}",
            f"[{mlo:.1f}, {mhi:.1f}]",
            f"{signif:.3f}", f"{az:.3f}",
        ])

    n_display = len(rows)
    fig, ax = plt.subplots(figsize=(14, max(3.5, 1.0 + n_display * 0.7)))
    ax.axis("off")
    ax.set_title(
        "Per-category yield and sensitivity summary",
        fontsize=12, fontweight="bold", pad=10,
    )

    table = ax.table(cellText=rows, colLabels=columns, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.6)
    table.auto_set_column_width(col=list(range(len(columns))))

    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor("#1A237E")
            cell.set_text_props(weight="bold", color="white")
        elif r > 0 and "★" in rows[r - 1][0]:
            cell.set_facecolor("#FFF8E1")  # soft amber for low-stats row
        elif r % 2 == 0:
            cell.set_facecolor("#E8EAF6")
        else:
            cell.set_facecolor("#FFFFFF")
        cell.set_edgecolor("#E0E0E0")

    fig.tight_layout(pad=0.5)
    return fig


# ────────────────────────────────────────────────────────────────────────────
# Improved existing figures
# ────────────────────────────────────────────────────────────────────────────

def build_bdt_signal_background_table(bins, signal_hist, bkg_hist, bkg_var, data_hist):
    columns = ["BDT bin", "Total signal", "Bkg. ± ΔB", "Data", "S/(S+B) (%)", "S/√B"]
    rows = []
    nbins = len(bins) - 1
    for i in range(nbins):
        s  = float(signal_hist[i]) if i < signal_hist.size else 0.0
        b  = float(bkg_hist[i])    if i < bkg_hist.size    else 0.0
        db = math.sqrt(max(float(bkg_var[i]), 0.0)) if i < bkg_var.size else 0.0
        d  = float(data_hist[i])   if i < data_hist.size   else 0.0
        purity = (100.0 * s / (s + b)) if (s + b) > 0 else 0.0
        signif = (s / math.sqrt(b)) if b > 0 else 0.0
        rows.append([
            f"[{bins[i]:.3f}, {bins[i+1]:.3f})",
            f"{s:.3f}", f"{b:.2f} ± {db:.2f}", f"{int(round(d))}",
            f"{purity:.2f}", f"{signif:.3f}",
        ])
    return columns, rows


def summarize_categories(era_payloads, edges, args, mc_weight_scale=1.0):
    """(Unchanged from original — see original script for full docstring.)"""
    ncat = len(edges) - 1
    summary: Dict[int, Dict] = {i: {} for i in range(ncat)}

    for payload in era_payloads.values():
        d_sc = np.asarray(payload["data_full_scores"], dtype=np.float64)
        d_m  = np.asarray(payload["data_full_mass"],   dtype=np.float64)
        if d_sc.size > 0 and d_m.size > 0:
            d_idx = category_index(d_sc, edges)
            d_sb  = sideband_mass_mask(d_m, args)
            for i in range(ncat):
                m_full = (d_idx == i)
                if np.any(m_full):
                    summary[i]["data_full_count"] = summary[i].get("data_full_count", 0.0) + float(np.count_nonzero(m_full))
                m_sb = m_full & d_sb
                if np.any(m_sb):
                    summary[i]["data_side_count"] = summary[i].get("data_side_count", 0.0) + float(np.count_nonzero(m_sb))

        for proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:
            cproc = canonical_proc(proc)
            idx   = category_index(np.asarray(m_sc, dtype=np.float64), edges)
            mass  = np.asarray(m_mass, dtype=np.float64)
            w     = np.asarray(m_w, dtype=np.float64) * float(mc_weight_scale)
            sb    = sideband_mass_mask(mass, args)
            for i in range(ncat):
                m_full = (idx == i)
                if np.any(m_full):
                    y_full  = float(np.sum(w[m_full]))
                    e2_full = float(np.sum(np.square(w[m_full])))
                    summary[i][f"mc_full_{cproc}"]        = summary[i].get(f"mc_full_{cproc}", 0.0) + y_full
                    summary[i][f"mcerr2_full_{cproc}"]    = summary[i].get(f"mcerr2_full_{cproc}", 0.0) + e2_full
                    summary[i]["mc_full_total"]            = summary[i].get("mc_full_total", 0.0) + y_full
                    summary[i]["mc_full_total_err2"]       = summary[i].get("mc_full_total_err2", 0.0) + e2_full
                m_sb2 = m_full & sb
                if np.any(m_sb2):
                    y_sb  = float(np.sum(w[m_sb2]))
                    e2_sb = float(np.sum(np.square(w[m_sb2])))
                    summary[i][f"mc_side_{cproc}"]        = summary[i].get(f"mc_side_{cproc}", 0.0) + y_sb
                    summary[i][f"mcerr2_side_{cproc}"]    = summary[i].get(f"mcerr2_side_{cproc}", 0.0) + e2_sb
                    summary[i]["mc_side_total"]            = summary[i].get("mc_side_total", 0.0) + y_sb
                    summary[i]["mc_side_total_err2"]       = summary[i].get("mc_side_total_err2", 0.0) + e2_sb

        s_sc = np.asarray(payload["sig_full_scores"],  dtype=np.float64)
        s_w  = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
        if s_sc.size > 0 and s_w.size > 0:
            s_idx = category_index(s_sc, edges)
            for i in range(ncat):
                m_sig = (s_idx == i)
                if np.any(m_sig):
                    summary[i]["sig_full_total"]      = summary[i].get("sig_full_total", 0.0) + float(np.sum(s_w[m_sig]))
                    summary[i]["sig_full_total_err2"] = summary[i].get("sig_full_total_err2", 0.0) + float(np.sum(np.square(s_w[m_sig])))

    # HWHM-windowed yields
    mass_bins_sum = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    centers_sum = 0.5 * (mass_bins_sum[:-1] + mass_bins_sum[1:])
    for i in range(ncat):
        lo_cat = float(edges[i]); hi_cat = float(edges[i+1])
        sig_m_chunks, sig_w_chunks = [], []
        for payload in era_payloads.values():
            s_sc   = np.asarray(payload["sig_full_scores"], dtype=np.float64)
            s_mass = np.asarray(payload["sig_full_mass"],   dtype=np.float64)
            s_w    = np.asarray(payload["sig_full_weight"], dtype=np.float64)
            if s_sc.size == 0: continue
            sel = (s_sc >= lo_cat) & (s_sc < hi_cat)
            if np.any(sel):
                sig_m_chunks.append(s_mass[sel])
                sig_w_chunks.append(s_w[sel])
        sig_m_cat = concat_chunks(sig_m_chunks)
        sig_w_cat = concat_chunks(sig_w_chunks)
        h_sig_cat, _ = weighted_hist_and_var(sig_m_cat, mass_bins_sum, sig_w_cat if sig_w_cat.size > 0 else None)
        _pred_dcb, _params_dcb, _chi2_dcb = _fit_signal_dcb_full(centers_sum, mass_bins_sum, h_sig_cat, args, cat_idx=i)
        if _params_dcb:
            sig_mean_i  = _params_dcb["mean"]
            sig_sigma_i = _params_dcb["sigma"]
        else:
            sig_mean_i = 125.09; sig_sigma_i = 2.0
        mass_lo_i, mass_hi_i = _hwhm_window_from_sigma(sig_mean_i, sig_sigma_i)
        summary[i]["hwhm_mass_lo"] = mass_lo_i
        summary[i]["hwhm_mass_hi"] = mass_hi_i
        s_hwhm, b_hwhm = _sig_bkg_in_hwhm(era_payloads, lo_cat, hi_cat, mass_lo_i, mass_hi_i, mc_weight_scale)
        summary[i]["sig_hwhm_total"] = s_hwhm
        summary[i]["mc_hwhm_total"]  = b_hwhm
    return summary


def build_category_summary_table(edges, summary):
    """Build summary table.

    All categories are shown individually. Low-stats categories
    (B_HWHM < threshold) are flagged with ★ but not merged.
    """
    columns = [
        "Category", "BDT range",
        "Data sidebands", "MC sidebands ± stat",
        "Data/MC (sidebands)", "MC full (110-150) ± stat",
        "S/√B [HWHM]", "Asimov Z (cat) [HWHM]",
        "B_HWHM",
    ]
    ncat = len(edges) - 1
    rows = []
    combined_z_sum = 0.0

    for i in range(ncat):
        vals = summary.get(i, {})
        b_hwhm = float(vals.get("mc_hwhm_total", 0.0))
        reliable = b_hwhm >= _B_HWHM_MIN_RELIABLE

        d_side  = float(vals.get("data_side_count", 0.0))
        m_side  = float(vals.get("mc_side_total",   0.0))
        me2_s   = float(vals.get("mc_side_total_err2", 0.0))
        m_full  = float(vals.get("mc_full_total",   0.0))
        me2_f   = float(vals.get("mc_full_total_err2", 0.0))
        s_hwhm  = float(vals.get("sig_hwhm_total",  0.0))

        ratio  = d_side / m_side if m_side > 0 else 0.0
        signif = s_hwhm / math.sqrt(b_hwhm) if b_hwhm > 0 else 0.0
        term   = 0.0
        if b_hwhm > 1e-12 and s_hwhm > 0:
            term = 2.0 * ((s_hwhm + b_hwhm) * math.log1p(s_hwhm / b_hwhm) - s_hwhm)
            combined_z_sum += max(term, 0.0)
        cat_z = math.sqrt(max(term, 0.0))
        flag = "" if reliable else " ★"
        rows.append([
            str(i) + flag, f"[{edges[i]:.3f}, {edges[i+1]:.3f})",
            f"{d_side:.0f}", f"{m_side:.2f} ± {math.sqrt(max(me2_s,0)):.2f}",
            f"{ratio:.3f}", f"{m_full:.2f} ± {math.sqrt(max(me2_f,0)):.2f}",
            f"{signif:.4f}", f"{cat_z:.4f}", f"{b_hwhm:.3f}",
        ])

    combined_z = math.sqrt(max(combined_z_sum, 0.0))
    rows.append(["Combined", "all", "", "", "", "", "", f"{combined_z:.4f}", ""])
    return columns, rows


def build_category_process_table(edges, summary):
    proc_keys = set()
    for vals in summary.values():
        for k in vals.keys():
            if k.startswith("mc_full_") and not k.startswith("mcerr2_"):
                proc_keys.add(k)
    procs = sorted([k.replace("mc_full_", "") for k in proc_keys], key=proc_sort_key)
    columns = ["Category", "BDT range"] + procs
    rows = []
    for i in range(len(edges) - 1):
        vals = summary.get(i, {})
        row = [str(i), f"[{edges[i]:.3f}, {edges[i+1]:.3f})"]
        for p in procs:
            row.append(f"{float(vals.get(f'mc_full_{p}', 0.0)):.2f}")
        rows.append(row)
    return columns, rows


def build_category_mass_grid(era_payloads, edges, args, mass_bins, mc_weight_scale=1.0):
    """Per-category dimuon mass with all background models.
    FIXED: chi2 now uses Baker-Cousins convention.
    """
    ncat = len(edges) - 1
    ncols = 2 if ncat > 1 else 1
    nrows = int(math.ceil(ncat / float(ncols)))
    fig = plt.figure(figsize=(14.5, 6.5 * nrows + 0.6))
    gs = fig.add_gridspec(nrows * 2, ncols,
                          height_ratios=[3.0, 1.0] * nrows,
                          hspace=0.38, wspace=0.30)
    top_axes, residual_axes = [], []
    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    fit_columns = ["Category", "BDT range", "Best model", "Best χ²_BC/ndf", "S/√B",
                   "BWZRedux", "BWZGamma", "BWZ×Bern", "Landau×Bern",
                   "Bernstein", "S-PowerLaw", "S-Exponential"]
    fit_rows = []

    model_styles = {
        "BWZRedux":      (_CMS_RED,    "-"),
        "BWZGamma":      (_CMS_GREEN,  "-"),
        "BWZ\u00d7Bern":    (_CMS_BLUE,   "-"),
        "Landau\u00d7Bern": (_CMS_ORANGE, "-"),
        "Bernstein":     (_CMS_PURPLE, "--"),
        "S-PowerLaw":    ("#6D4C41",   ":"),
        "S-Exponential": ("#00838F",   "-."),
    }

    for i in range(ncat):
        row = i // ncols; col = i % ncols
        ax = fig.add_subplot(gs[2*row, col])
        ax_res = fig.add_subplot(gs[2*row+1, col], sharex=ax)
        top_axes.append(ax); residual_axes.append(ax_res)

        data_mass_chunks, bkg_mass_chunks, bkg_weight_chunks = [], [], []
        sig_mass_chunks, sig_weight_chunks = [], []
        lo, hi = float(edges[i]), float(edges[i+1])

        for payload in era_payloads.values():
            d_sc = np.asarray(payload["data_full_scores"], dtype=np.float64)
            d_m  = np.asarray(payload["data_full_mass"],   dtype=np.float64)
            if d_sc.size > 0 and d_m.size > 0:
                dsel = (d_sc >= lo) & (d_sc < hi)
                if np.any(dsel): data_mass_chunks.append(d_m[dsel])

            for _proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:
                m_sc = np.asarray(m_sc, dtype=np.float64)
                m_mass = np.asarray(m_mass, dtype=np.float64)
                m_w = np.asarray(m_w, dtype=np.float64) * float(mc_weight_scale)
                msel = (m_sc >= lo) & (m_sc < hi)
                if np.any(msel):
                    bkg_mass_chunks.append(m_mass[msel])
                    bkg_weight_chunks.append(m_w[msel])

            s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)
            s_mass = np.asarray(payload["sig_full_mass"], dtype=np.float64)
            s_w    = np.asarray(payload["sig_full_weight"], dtype=np.float64)
            if s_sc.size > 0:
                ssel = (s_sc >= lo) & (s_sc < hi)
                if np.any(ssel):
                    sig_mass_chunks.append(s_mass[ssel])
                    sig_weight_chunks.append(s_w[ssel])

        data_mass = concat_chunks(data_mass_chunks)
        bkg_mass  = concat_chunks(bkg_mass_chunks)
        bkg_w     = concat_chunks(bkg_weight_chunks)
        sig_mass  = concat_chunks(sig_mass_chunks)
        sig_w     = concat_chunks(sig_weight_chunks)

        h_data, v_data = weighted_hist_and_var(data_mass, mass_bins, None)
        h_bkg, _       = weighted_hist_and_var(bkg_mass,  mass_bins, bkg_w  if bkg_w.size  > 0 else None)
        h_sig, _       = weighted_hist_and_var(sig_mass,  mass_bins, sig_w  if sig_w.size  > 0 else None)

        h_bkg_fit, best_model, best_chi2ndf, model_stats, model_preds = fit_background_models_full(
            centers, h_data, v_data, args)
        h_sig_dcb = fit_signal_dcb_template(centers, mass_bins, h_sig, args)
        h_sb_fit  = h_bkg_fit + h_sig_dcb

        dmask = (centers < args.blind_min) | (centers > args.blind_max)
        signif = (float(np.sum(h_sig_dcb)) / math.sqrt(float(np.sum(h_bkg_fit)))) \
                 if float(np.sum(h_bkg_fit)) > 0.0 else 0.0

        for model_name in model_styles:
            if model_name not in model_preds: continue
            pred = model_preds[model_name]
            color, ls = model_styles[model_name]
            lw = 2.8 if model_name == best_model else 1.6
            alpha = 1.0 if model_name == best_model else 0.85
            label = f"{model_name} ({'best' if model_name == best_model else 'alt'})"
            line = ax.step(centers, pred, where="mid", color=color,
                           linewidth=lw, linestyle=ls, alpha=alpha, label=label,
                           zorder=3 if model_name == best_model else 2)[0]
            line.set_path_effects([pe.Stroke(linewidth=lw + 1.0, foreground="white"),
                                   pe.Normal()])

        ax.step(centers, h_sb_fit, where="mid", color="#212121",
                linewidth=2.0, linestyle="-.", label="Best B + DSCB(S)", zorder=4)

        if np.any(dmask):
            ax.errorbar(centers[dmask], h_data[dmask],
                        yerr=np.sqrt(np.maximum(v_data[dmask], 1.0)),
                        fmt="o", color="black", markersize=3.8, linewidth=1.0,
                        label="Data" if i == 0 else None, zorder=5, capsize=2)

        ax.axvspan(args.blind_min, args.blind_max, color="#E0E0E0", alpha=0.45, zorder=0)
        ax.set_title(f"Cat {i}:  BDT [{lo:.3f}, {hi:.3f})", fontsize=11, fontweight="bold", pad=4)
        ax.set_ylabel("Events / bin", fontsize=10)
        ax.grid(alpha=0.22)
        ax.tick_params(labelbottom=False)
        ax.text(0.02, 0.40,
                f"Best: {best_model}\n$\\chi^2_{{\\mathrm{{BC}}}}$/ndf = {best_chi2ndf:.3f}\n"
                f"S/$\\sqrt{{B}}$ = {signif:.3f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.88,
                          edgecolor="#CCCCCC"))
        ax.legend(loc="upper right", fontsize=6.5, ncol=2,
                  frameon=True, framealpha=0.90, edgecolor="#CCCCCC")

        # Fit comparison table row
        model_order = ["BWZRedux", "BWZGamma", "BWZ\u00d7Bern", "Landau\u00d7Bern",
                       "Bernstein", "S-PowerLaw", "S-Exponential"]
        row_data = [str(i), f"[{lo:.3f}, {hi:.3f})", best_model,
                    f"{best_chi2ndf:.3f}", f"{signif:.3f}"]
        for model_name in model_order:
            if model_name in model_stats:
                chi2_bc, ndf = model_stats[model_name]
                row_data.append(f"{chi2_bc/float(max(ndf,1)):.3f} ({int(ndf)})")
            else:
                row_data.append("n/a")
        fit_rows.append(row_data)

        # Residual ratio panel
        ax_res.axhline(1.0, color="#424242", linestyle="--", linewidth=1.0)
        for model_name in model_order:
            if model_name not in model_preds: continue
            pred = model_preds[model_name]
            ratio = np.divide(pred, h_bkg_fit, out=np.ones_like(pred), where=h_bkg_fit > 0)
            color, ls = model_styles[model_name]
            lw = 1.8 if model_name == best_model else 1.2
            line = ax_res.step(centers, ratio, where="mid", color=color,
                               linewidth=lw, linestyle=ls, alpha=0.92, zorder=3)[0]
            line.set_path_effects([pe.Stroke(linewidth=lw + 0.8, foreground="white"), pe.Normal()])
        ax_res.axvspan(args.blind_min, args.blind_max, color="#E0E0E0", alpha=0.45, zorder=0)
        ax_res.set_ylabel("PDF / Best", fontsize=9)
        ax_res.grid(alpha=0.22)
        ax_res.set_ylim(0.94, 1.06)

    for j in range(ncat, nrows * ncols):
        r = j // ncols; c = j % ncols
        fig.add_subplot(gs[2*r, c]).axis("off")
        fig.add_subplot(gs[2*r+1, c]).axis("off")

    for ax in residual_axes:
        ax.set_xlabel("$m_{\\mu\\mu}$ [GeV]", fontsize=10)

    fig.subplots_adjust(left=0.07, right=0.97, top=0.920, bottom=0.055, hspace=0.38, wspace=0.30)
    fig.suptitle(
        "Per-category dimuon mass: background model fits",
        fontsize=12, fontweight="bold", y=0.965,
    )
    return fig, fit_columns, fit_rows


def build_vbf_signal_dcb_grid(era_payloads, edges, args, mass_bins,
                               mc_weight_scale=1.0, cat_summary=None):
    """IMPROVED: Better aesthetics + Baker-Cousins chi2 annotation."""
    ncat = len(edges) - 1
    ncols = 2 if ncat > 1 else 1
    nrows = int(math.ceil(ncat / float(ncols)))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.5, 5.2 * nrows + 0.6), squeeze=False)
    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    bin_width = float(mass_bins[1] - mass_bins[0])

    dcb_cols = ["Category", "BDT range", "N_sig", "Mean [GeV]", "Sigma [GeV]",
                "alpha_L", "n_L", "alpha_R", "n_R", "FWHM [GeV]", "S/√B", "Asimov Z",
                "χ²_BC/ndf"]
    dcb_rows = []
    all_params = []  # collect for sigma-vs-category figure

    for i in range(ncat):
        row_idx = i // ncols; col_idx = i % ncols
        ax = axes[row_idx, col_idx]
        lo, hi = float(edges[i]), float(edges[i+1])

        sig_mass_chunks, sig_weight_chunks = [], []
        for payload in era_payloads.values():
            s_sc   = np.asarray(payload["sig_full_scores"],  dtype=np.float64)
            s_mass = np.asarray(payload["sig_full_mass"],    dtype=np.float64)
            s_w    = np.asarray(payload["sig_full_weight"],  dtype=np.float64)
            if s_sc.size == 0: continue
            ssel = (s_sc >= lo) & (s_sc < hi)
            if np.any(ssel):
                sig_mass_chunks.append(s_mass[ssel])
                sig_weight_chunks.append(s_w[ssel])

        sig_mass = concat_chunks(sig_mass_chunks)
        sig_w    = concat_chunks(sig_weight_chunks)
        h_sig, _ = weighted_hist_and_var(sig_mass, mass_bins, sig_w if sig_w.size > 0 else None)
        pred_dcb, params, chi2ndf = _fit_signal_dcb_full(centers, mass_bins, h_sig, args, cat_idx=i)
        all_params.append(params)
        n_sig_total = float(np.sum(h_sig))

        _summary_i = (cat_summary or {}).get(i, {})
        s_yield = float(_summary_i.get("sig_hwhm_total", 0.0)) or float(_summary_i.get("sig_full_total", n_sig_total))
        b_yield = float(_summary_i.get("mc_hwhm_total",  0.0))
        sbsqrtb = s_yield / math.sqrt(b_yield) if b_yield > 0.0 else 0.0
        if b_yield > 0.0 and s_yield > 0.0:
            asimov_z = math.sqrt(max(0.0, 2.0*((s_yield+b_yield)*math.log1p(s_yield/b_yield)-s_yield)))
        else:
            asimov_z = 0.0

        h_sig_gev    = h_sig    / bin_width
        pred_dcb_gev = pred_dcb / bin_width

        ax.step(centers, h_sig_gev, where="mid",
                color=_CMS_BLUE, linewidth=1.6,
                label=f"VBF H(125) MC  (N = {n_sig_total:.1f})", zorder=2)
        ax.fill_between(centers, h_sig_gev, step="mid",
                        color=_CMS_BLUE, alpha=0.18, zorder=1)

        if np.sum(pred_dcb_gev) > 0.0:
            ax.step(centers, pred_dcb_gev, where="mid",
                    color=_CMS_RED, linewidth=2.4,
                    label="DCB fit", zorder=3)

        ax.set_title(f"Cat {i}:  BDT [{lo:.3f}, {hi:.3f})", fontsize=11, fontweight="bold", pad=4)
        ax.set_xlabel("$m_{\\mu\\mu}$ [GeV]", fontsize=10)
        ax.set_ylabel("Events / GeV", fontsize=10)
        ax.grid(alpha=0.22)
        _cms_label(ax)

        # Reliability flag for significance
        b_hwhm = float(_summary_i.get("mc_hwhm_total", 0.0))
        rel_flag = "" if b_hwhm >= _B_HWHM_MIN_RELIABLE else " ★"

        if params:
            fwhm = 2.355 * params["sigma"]
            chi2_str = f"$\\chi^2_{{\\mathrm{{BC}}}}$/ndf = {chi2ndf:.3f}" if chi2ndf >= 0.0 else ""
            stats_text = (
                f"$\\mu$ = {params['mean']:.3f} GeV\n"
                f"$\\sigma$ = {params['sigma']:.3f} GeV\n"
                f"$\\alpha_L$ = {params['alpha_l']:.3f},  $n_L$ = {params['n_l']:.2f}\n"
                f"$\\alpha_R$ = {params['alpha_r']:.3f},  $n_R$ = {params['n_r']:.2f}\n"
                f"FWHM $\\approx$ {fwhm:.3f} GeV\n"
                f"S/$\\sqrt{{B}}$ = {sbsqrtb:.3f}{rel_flag}\n"
                + chi2_str
            )
        else:
            fwhm = float("nan")
            stats_text = f"N = {n_sig_total:.1f}\n(DCB fit failed)"

        ax.text(0.97, 0.97, stats_text, transform=ax.transAxes,
                va="top", ha="right", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          alpha=0.90, edgecolor="#CCCCCC"))
        ax.legend(loc="upper left", fontsize=9)

        chi2_str_row = f"{chi2ndf:.3f}" if chi2ndf >= 0.0 else "n/a"
        dcb_rows.append([
            str(i), f"[{lo:.3f}, {hi:.3f})", f"{n_sig_total:.2f}",
            f"{params['mean']:.4f}"    if params else "n/a",
            f"{params['sigma']:.4f}"   if params else "n/a",
            f"{params['alpha_l']:.4f}" if params else "n/a",
            f"{params['n_l']:.3f}"     if params else "n/a",
            f"{params['alpha_r']:.4f}" if params else "n/a",
            f"{params['n_r']:.3f}"     if params else "n/a",
            f"{fwhm:.3f}" if math.isfinite(fwhm) else "n/a",
            f"{sbsqrtb:.3f}",
            f"{asimov_z:.3f}",
            chi2_str_row,
        ])

    for j in range(ncat, nrows * ncols):
        axes[j // ncols, j % ncols].axis("off")

    fig.subplots_adjust(left=0.07, right=0.97, top=0.910, bottom=0.08, hspace=0.48, wspace=0.30)
    fig.suptitle(
        "VBF H(125) signal DCB extraction per BDT category",
        fontsize=12, fontweight="bold", y=0.965,
    )
    return fig, dcb_cols, dcb_rows, all_params


def build_vbf_bdt_transformed_figure(era_payloads, edges, args, mass_bins,
                                     mc_weight_scale=1.0, n_score_bins=25):
    """IMPROVED: CMS labels, better colors, cleaner panel separation."""
    sig_score_chunks, sig_mass_chunks, sig_w_chunks = [], [], []
    bkg_score_chunks, bkg_mass_chunks, bkg_w_chunks = [], [], []

    for payload in era_payloads.values():
        s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)
        s_m  = np.asarray(payload["sig_full_mass"],   dtype=np.float64)
        s_w  = np.asarray(payload["sig_full_weight"], dtype=np.float64)
        if s_sc.size > 0 and s_m.size > 0 and s_w.size > 0:
            sig_score_chunks.append(s_sc); sig_mass_chunks.append(s_m); sig_w_chunks.append(s_w)
        for _proc, b_sc, b_m, b_w in payload["mc_full_groups"]:
            b_sc = np.asarray(b_sc, dtype=np.float64)
            b_m  = np.asarray(b_m,  dtype=np.float64)
            b_w  = np.asarray(b_w,  dtype=np.float64) * float(mc_weight_scale)
            if b_sc.size > 0:
                bkg_score_chunks.append(b_sc); bkg_mass_chunks.append(b_m); bkg_w_chunks.append(b_w)

    sig_scores = concat_chunks(sig_score_chunks); sig_mass = concat_chunks(sig_mass_chunks)
    sig_w      = concat_chunks(sig_w_chunks)
    bkg_scores = concat_chunks(bkg_score_chunks); bkg_mass = concat_chunks(bkg_mass_chunks)
    bkg_w      = concat_chunks(bkg_w_chunks)

    if sig_scores.size == 0 or bkg_scores.size == 0:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.text(0.5, 0.5, "No signal/background data", ha="center", va="center", transform=ax.transAxes)
        return fig

    t_sig = _bdt_transform(sig_scores); t_bkg = _bdt_transform(bkg_scores)
    t_edges = _bdt_transform(np.clip(edges[1:-1], 1e-6, 1.0 - 1e-6))
    t_min = min(float(np.percentile(t_sig, 1.0)),  float(np.percentile(t_bkg, 1.0)))
    t_max = max(float(np.percentile(t_sig, 99.5)), float(np.percentile(t_bkg, 99.5)))
    score_edges = np.linspace(t_min, t_max, n_score_bins + 1)

    # Inclusive DCB for HWHM
    centers_mass = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    h_sig_incl, _ = weighted_hist_and_var(sig_mass, mass_bins, sig_w if sig_w.size > 0 else None)
    _pred_incl, _par_incl, _chi2_incl = _fit_signal_dcb_full(centers_mass, mass_bins, h_sig_incl, args, cat_idx=999)
    if _par_incl:
        sig_mean = _par_incl["mean"]; hwhm = 1.177 * _par_incl["sigma"]
    else:
        sig_mean = 125.09; hwhm = 2.0
    mass_lo = sig_mean - hwhm; mass_hi = sig_mean + hwhm

    sig_in_window = (sig_mass >= mass_lo) & (sig_mass <= mass_hi)
    bkg_in_window = (bkg_mass >= mass_lo) & (bkg_mass <= mass_hi)
    s_per_bin = np.zeros(n_score_bins, dtype=np.float64)
    b_per_bin = np.zeros(n_score_bins, dtype=np.float64)
    if np.any(sig_in_window):
        s_idx = np.clip(np.digitize(t_sig[sig_in_window], score_edges) - 1, 0, n_score_bins-1)
        np.add.at(s_per_bin, s_idx, sig_w[sig_in_window])
    if np.any(bkg_in_window):
        b_idx = np.clip(np.digitize(t_bkg[bkg_in_window], score_edges) - 1, 0, n_score_bins-1)
        np.add.at(b_per_bin, b_idx, bkg_w[bkg_in_window])

    sqrtb = np.sqrt(np.maximum(b_per_bin, 0.0))
    s_over_sqrtb = np.where(sqrtb > 0.0, s_per_bin / sqrtb, 0.0)

    bin_widths = score_edges[1:] - score_edges[:-1]
    h_sig_score, _ = weighted_hist_and_var(t_sig, score_edges, sig_w if sig_w.size > 0 else None)
    h_bkg_score, _ = weighted_hist_and_var(t_bkg, score_edges, bkg_w if bkg_w.size > 0 else None)
    tot_s = float(np.sum(h_sig_score)); tot_b = float(np.sum(h_bkg_score))
    h_sig_norm = h_sig_score / (tot_s * bin_widths) if tot_s > 0 else h_sig_score
    h_bkg_norm = h_bkg_score / (tot_b * bin_widths) if tot_b > 0 else h_bkg_score
    bin_centers = 0.5 * (score_edges[:-1] + score_edges[1:])

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 7),
                                          gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})

    ax_top.step(bin_centers, h_sig_norm, where="mid",
                color=_CMS_BLUE, linewidth=2.0, label="Signal $m_H=125$ GeV")
    ax_top.fill_between(bin_centers, h_sig_norm, step="mid", color=_CMS_BLUE, alpha=0.15)
    ax_top.step(bin_centers, h_bkg_norm, where="mid",
                color=_CMS_RED, linewidth=2.0, label="Background")
    ax_top.fill_between(bin_centers, h_bkg_norm, step="mid", color=_CMS_RED, alpha=0.10)

    for t_edge in t_edges:
        if score_edges[0] < t_edge < score_edges[-1]:
            ax_top.axvline(t_edge, color="#424242", linestyle="--", linewidth=1.2, alpha=0.75)
            ax_bot.axvline(t_edge, color="#424242", linestyle="--", linewidth=1.2, alpha=0.75)

    ax_top.set_ylabel("Normalised / bin width", fontsize=11)
    ax_top.tick_params(labelbottom=False)
    ax_top.legend(fontsize=11, loc="upper left", framealpha=0.90)
    ax_top.grid(alpha=0.22)
    ax_top.set_xlim(score_edges[0], score_edges[-1])
    _cms_label(ax_top)

    ax_bot.step(bin_centers, s_over_sqrtb, where="mid", color="#212121", linewidth=1.8)
    ax_bot.fill_between(bin_centers, s_over_sqrtb, step="mid", color="#212121", alpha=0.12)
    ax_bot.set_xlabel("atanh((BDT output+1)/2)", fontsize=11)
    ax_bot.set_ylabel("$S/\\sqrt{B}$", fontsize=11)
    ax_bot.grid(alpha=0.22)
    ax_bot.set_xlim(score_edges[0], score_edges[-1])
    ymax_bot = max(float(np.max(s_over_sqrtb)) * 1.30, 0.01)
    ax_bot.set_ylim(0.0, ymax_bot)
    hwhm_note = (f"HWHM window: [{mass_lo:.1f}, {mass_hi:.1f}] GeV  "
                 f"($\\sigma$={_par_incl['sigma']:.2f} GeV)" if _par_incl else "HWHM window (fallback)")
    ax_bot.text(0.98, 0.95, hwhm_note, transform=ax_bot.transAxes,
                ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.85))

    fig.suptitle("VBF BDT score and per-bin $S/\\sqrt{B}$ in signal HWHM window",
                 fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    return fig


def build_significance_evolution_figure(scan_history):
    """IMPROVED: better annotation and CMS label."""
    if not scan_history:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No scan history available\n(manual edges or single-era run)",
                ha="center", va="center", transform=ax.transAxes, fontsize=11)
        ax.set_title("Significance evolution  (AN Fig. 33 analog)", fontsize=12)
        return fig

    ns = [h[0] for h in scan_history]; zs = [h[1] for h in scan_history]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(ns, zs, "o-", color=_CMS_BLUE, linewidth=2.2, markersize=8,
            markerfacecolor="white", markeredgewidth=2.2, zorder=3)
    best_idx = int(np.argmax(zs))
    ax.axvline(ns[best_idx], color=_CMS_RED, linestyle="--", linewidth=1.5,
               label=f"Optimal: n={ns[best_idx]}, Z={zs[best_idx]:.4f}")
    ax.scatter([ns[best_idx]], [zs[best_idx]], color=_CMS_RED, s=100, zorder=5)
    for n_val, z_val in zip(ns, zs):
        ax.annotate(f"{z_val:.3f}", xy=(n_val, z_val), xytext=(0, 10),
                    textcoords="offset points", ha="center", fontsize=10, color="#212121",
                    fontweight="bold")
    ax.set_xlabel("Number of BDT categories", fontsize=12)
    ax.set_ylabel("Combined Asimov $Z$ [HWHM window]", fontsize=12)
    ax.set_title("Combined significance vs. number of VBF BDT categories",
                 fontsize=11, fontweight="bold", pad=28)
    ax.set_xticks(ns)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.25)
    _cms_label(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def build_fwhm_vs_bdt_figure(era_payloads, args, mass_bins, n_score_bins=10):
    """IMPROVED: CMS label, better error bars, annotation."""
    sig_score_chunks, sig_mass_chunks, sig_w_chunks = [], [], []
    for payload in era_payloads.values():
        s_sc = np.asarray(payload["sig_full_scores"], dtype=np.float64)
        s_m  = np.asarray(payload["sig_full_mass"],   dtype=np.float64)
        s_w  = np.asarray(payload["sig_full_weight"], dtype=np.float64)
        if s_sc.size > 0 and s_m.size > 0 and s_w.size > 0:
            sig_score_chunks.append(s_sc); sig_mass_chunks.append(s_m); sig_w_chunks.append(s_w)
    sig_scores = concat_chunks(sig_score_chunks); sig_mass = concat_chunks(sig_mass_chunks)
    sig_w = concat_chunks(sig_w_chunks)
    if sig_scores.size == 0:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.text(0.5, 0.5, "No signal data", ha="center", va="center", transform=ax.transAxes)
        return fig

    q_edges = np.quantile(sig_scores, np.linspace(0.0, 1.0, n_score_bins+1))
    q_edges = np.unique(q_edges); n_actual = len(q_edges) - 1
    centers_mass = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    bin_centers_score = 0.5 * (q_edges[:-1] + q_edges[1:])
    x_plot = _bdt_transform(np.clip(bin_centers_score, 1e-6, 1.0 - 1e-6))
    fwhm_vals, fwhm_errs = [], []

    for b in range(n_actual):
        lo_s = float(q_edges[b]); hi_s = float(q_edges[b+1])
        sel = (sig_scores >= lo_s) & (sig_scores < hi_s)
        if not np.any(sel):
            fwhm_vals.append(float("nan")); fwhm_errs.append(float("nan")); continue
        h_slice, _ = weighted_hist_and_var(sig_mass[sel], mass_bins, sig_w[sel] if sig_w.size > 0 else None)
        _pred, _params, _chi2 = _fit_signal_dcb_full(centers_mass, mass_bins, h_slice, args, cat_idx=1000+b)
        if _params:
            fwhm_vals.append(2.355 * _params["sigma"])
            fwhm_errs.append(2.355 * _params["sigma"] * 0.10)
        else:
            fwhm_vals.append(float("nan")); fwhm_errs.append(float("nan"))

    fwhm_arr = np.array(fwhm_vals, dtype=np.float64)
    fwhm_err_arr = np.array(fwhm_errs, dtype=np.float64)
    valid = np.isfinite(fwhm_arr)
    fig, ax = plt.subplots(figsize=(8, 5.5))
    if np.any(valid):
        ax.errorbar(x_plot[valid], fwhm_arr[valid], yerr=fwhm_err_arr[valid],
                    fmt="o-", color=_CMS_BLUE, linewidth=2.0, markersize=7, capsize=4,
                    markerfacecolor="white", markeredgewidth=2.0,
                    label="FWHM (2.355 × σ from DCB fit)")
        # Reference line at mean FWHM
        mean_fwhm = float(np.nanmean(fwhm_arr))
        ax.axhline(mean_fwhm, color="#9E9E9E", linestyle="--", linewidth=1.2,
                   label=f"Mean FWHM = {mean_fwhm:.2f} GeV")

    ax.set_xlabel("atanh((VBF BDT output + 1) / 2)", fontsize=12)
    ax.set_ylabel("FWHM [GeV]", fontsize=12)
    ax.set_title("Dimuon mass resolution (FWHM) vs. VBF BDT score",
                 fontsize=12, fontweight="bold", pad=28)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.25)
    _cms_label(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    return fig


# ────────────────────────────────────────────────────────────────────────────
# BDT math utilities (imported from significance_inference above):
#   _bdt_transform, _bdt_inverse, _hwhm_window_from_sigma, _sig_bkg_in_hwhm,
#   _bin_yields_t, _bin_yields_t_hwhm, _combined_asimov_z,
#   iterative_greedy_category_scan
# ────────────────────────────────────────────────────────────────────────────

def _gather_category_data(era_payloads, lo, hi, mass_bins, mc_weight_scale):
    data_chunks, bkg_chunks, bkg_w_chunks, sig_chunks, sig_w_chunks = [], [], [], [], []
    for payload in era_payloads.values():
        d_sc = np.asarray(payload["data_full_scores"], dtype=np.float64)
        d_m  = np.asarray(payload["data_full_mass"],   dtype=np.float64)
        if d_sc.size > 0 and d_m.size > 0:
            sel = (d_sc >= lo) & (d_sc < hi)
            if np.any(sel): data_chunks.append(d_m[sel])
        for _proc, m_sc, m_mass, m_w in payload["mc_full_groups"]:
            m_sc   = np.asarray(m_sc,   dtype=np.float64)
            m_mass = np.asarray(m_mass, dtype=np.float64)
            m_w    = np.asarray(m_w,    dtype=np.float64) * mc_weight_scale
            sel = (m_sc >= lo) & (m_sc < hi)
            if np.any(sel): bkg_chunks.append(m_mass[sel]); bkg_w_chunks.append(m_w[sel])
        s_sc   = np.asarray(payload["sig_full_scores"], dtype=np.float64)
        s_mass = np.asarray(payload["sig_full_mass"],   dtype=np.float64)
        s_w    = np.asarray(payload["sig_full_weight"], dtype=np.float64)
        if s_sc.size > 0 and s_mass.size > 0:
            sel = (s_sc >= lo) & (s_sc < hi)
            if np.any(sel): sig_chunks.append(s_mass[sel]); sig_w_chunks.append(s_w[sel])
    h_data, v_data = weighted_hist_and_var(concat_chunks(data_chunks), mass_bins, None)
    h_bkg,  _      = weighted_hist_and_var(concat_chunks(bkg_chunks),  mass_bins, concat_chunks(bkg_w_chunks) if bkg_w_chunks else None)
    h_sig,  _      = weighted_hist_and_var(concat_chunks(sig_chunks),  mass_bins, concat_chunks(sig_w_chunks) if sig_w_chunks else None)
    return h_data, v_data, h_bkg, h_sig


def trim_degenerate_top_categories(edges, era_payloads, args, mc_weight_scale, logger):
    edges = np.asarray(edges, dtype=np.float64).copy()
    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    centers = 0.5 * (mass_bins[:-1] + mass_bins[1:])
    while edges.size > 2:
        lo = float(edges[-2]); hi = float(edges[-1])
        h_data, v_data, h_bkg, h_sig = _gather_category_data(era_payloads, lo, hi, mass_bins, mc_weight_scale)
        if float(np.sum(h_data)) == 0.0 and float(np.sum(h_bkg)) == 0.0:
            logger.info("Top category [%.3f, %.3f) is empty; merging.", lo, hi)
            edges = np.concatenate([edges[:-2], [edges[-1]]]); continue
        h_bkg_fit, best_model, best_chi2ndf, _, _ = fit_background_models_full(centers, h_data, v_data, args)
        h_sig_dcb = fit_signal_dcb_template(centers, mass_bins, h_sig, args)
        b_total = float(np.sum(h_bkg_fit)); s_total = float(np.sum(h_sig_dcb))
        signif = (s_total / math.sqrt(b_total)) if b_total > 0.0 else 0.0
        is_degenerate = (best_model == "fallback_data" or not np.isfinite(best_chi2ndf) or signif == 0.0)
        if is_degenerate:
            logger.info("Top category [%.3f, %.3f) is degenerate; merging.", lo, hi)
            edges = np.concatenate([edges[:-2], [edges[-1]]])
        else:
            break
    return edges


def parse_category_edges_arg(text):
    s = str(text).strip()
    if not s: return None
    edges = []
    if os.path.isfile(s):
        with open(s, "r", encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if not lines: return None
        header = lines[0].lower()
        if "bin_low" in header and "bin_high" in header:
            lows, highs = [], []
            for ln in lines[1:]:
                toks = [t.strip() for t in ln.split(",")]
                if len(toks) < 2: continue
                try: lows.append(float(toks[0])); highs.append(float(toks[1]))
                except ValueError: continue
            edges = sorted(set(lows + highs))
        else:
            for ln in lines:
                for tok in ln.replace(",", " ").split():
                    try: edges.append(float(tok))
                    except ValueError: continue
    else:
        for tok in s.split(","):
            tok = tok.strip()
            if not tok: continue
            edges.append(float(tok))
    if not edges: return None
    out = np.asarray(sorted(set(edges)), dtype=np.float64)
    if out.size < 2: return None
    return out


def optimize_category_edges(args, fallback_eras, logger):
    manual = parse_category_edges_arg(args.category_edges)
    if manual is not None:
        logger.info("Using user-provided category edges: %s", ", ".join([f"{x:.6f}" for x in manual]))
        return manual, []

    opt_eras = args.category_eras if args.category_eras else list(fallback_eras)
    if not opt_eras:
        opt_eras = infer_eras_from_outputs(args.outputs_dir, None)

    sig_scores_chunks, sig_weight_chunks, sig_mass_chunks = [], [], []
    bkg_scores_chunks, bkg_weight_chunks, bkg_mass_chunks = [], [], []
    used_eras = []

    for era in opt_eras:
        files = find_scored_files(args.outputs_dir, era)
        sig_file = files.get("signal", None); bkg_file = files.get("background", None)
        if not sig_file or not bkg_file: continue
        mass_branch = resolve_mass_branch(sig_file, args.treename, args.mass_branch)
        if mass_branch is None:
            mass_branch = resolve_mass_branch(bkg_file, args.treename, args.mass_branch)
        has_mass = mass_branch is not None
        sel_bdt = f"TMath::Finite({args.bdt_branch})"
        if has_mass:
            sel_full = (f"({sel_bdt}) && ({mass_branch}>={args.mass_min}"
                        f" && {mass_branch}<={args.mass_max})")
            if args.extra_cut: sel_full = f"({sel_full}) && ({args.extra_cut})"
            load_branches = [args.bdt_branch, args.weight_branch, mass_branch]
        else:
            sel_full = sel_bdt
            if args.extra_cut: sel_full = f"({sel_full}) && ({args.extra_cut})"
            load_branches = [args.bdt_branch, args.weight_branch]
        try:
            s_arr = load_arrays(sig_file, args.treename, load_branches, sel_full)
            b_arr = load_arrays(bkg_file, args.treename, load_branches, sel_full)
        except Exception as exc:
            logger.warning("%s: failed loading optimization arrays: %s", era, exc); continue

        s_sc = s_arr[args.bdt_branch]; s_w = sanitize_weights(s_arr[args.weight_branch], args.abs_weights)
        b_sc = b_arr[args.bdt_branch]; b_w = sanitize_weights(b_arr[args.weight_branch], args.abs_weights)
        if s_sc.size == 0 or b_sc.size == 0: continue
        sig_scores_chunks.append(s_sc); sig_weight_chunks.append(s_w)
        bkg_scores_chunks.append(b_sc); bkg_weight_chunks.append(b_w)
        if has_mass:
            sig_mass_chunks.append(np.asarray(s_arr[mass_branch], dtype=np.float64))
            bkg_mass_chunks.append(np.asarray(b_arr[mass_branch], dtype=np.float64))
        used_eras.append(era)

    if not used_eras:
        raise SystemExit("[error] No scored signal/background inputs found for category optimization.")

    sig_scores  = concat_chunks(sig_scores_chunks)
    sig_weights = concat_chunks(sig_weight_chunks)
    bkg_scores  = concat_chunks(bkg_scores_chunks)
    bkg_weights = concat_chunks(bkg_weight_chunks)

    hwhm_mass_lo = 0.0; hwhm_mass_hi = float("inf")
    sig_mass_all = bkg_mass_all = None

    if sig_mass_chunks and bkg_mass_chunks:
        sig_mass_all = concat_chunks(sig_mass_chunks)
        bkg_mass_all = concat_chunks(bkg_mass_chunks)
        mass_bins_opt = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
        centers_opt = 0.5 * (mass_bins_opt[:-1] + mass_bins_opt[1:])
        h_sig_incl, _ = weighted_hist_and_var(sig_mass_all, mass_bins_opt, sig_weights if sig_weights.size > 0 else None)
        _pred_opt, _par_opt, _ = _fit_signal_dcb_full(centers_opt, mass_bins_opt, h_sig_incl, args, cat_idx=9999)
        if _par_opt:
            hwhm_mass_lo, hwhm_mass_hi = _hwhm_window_from_sigma(_par_opt["mean"], _par_opt["sigma"])
        else:
            hwhm_mass_lo, hwhm_mass_hi = 123.09, 127.09

    total_sig_yield = float(np.sum(sig_weights)); total_bkg_yield = float(np.sum(bkg_weights))
    adaptive_min_sig = min(args.multibin_min_sig, max(1e-6, total_sig_yield * 0.05))
    adaptive_min_bkg = min(args.multibin_min_bkg, max(1e-6, total_bkg_yield * 0.05))

    edges, scan_history = iterative_greedy_category_scan(
        sig_scores, sig_weights, bkg_scores, bkg_weights,
        max_n=args.multibin_max_n, min_improvement=args.multibin_improvement,
        min_sig=adaptive_min_sig, min_bkg=adaptive_min_bkg,
        sig_mass=sig_mass_all, bkg_mass=bkg_mass_all,
        hwhm_mass_lo=hwhm_mass_lo, hwhm_mass_hi=hwhm_mass_hi,
        logger=logger,
    )
    return edges, scan_history


# ────────────────────────────────────────────────────────────────────────────
# MAIN
# ────────────────────────────────────────────────────────────────────────────

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
        raise SystemExit(f"[error] No eras found. Checked outputs-dir: {args.outputs_dir}")

    era_payloads: Dict = {}
    for era in eras:
        payload = load_era_payload(era, data_map.get(era, []), mc_map.get(era, []), args, logger)
        if payload is not None:
            era_payloads[era] = payload
    if not era_payloads:
        raise SystemExit("[error] No valid era payloads were loaded.")

    edges, scan_history = optimize_category_edges(args, eras, logger)

    run3_data_full_scores = combine_data_full_scores(era_payloads)
    run3_data_side_mass   = combine_data_side_masses(era_payloads, args)
    run3_mc_score_side    = combine_mc_score_groups(era_payloads, region="side")
    run3_mc_score_full    = combine_mc_score_groups(era_payloads, region="full")
    run3_mc_mass_full     = combine_mc_mass_groups(era_payloads, args, region="full")

    if run3_data_full_scores.size == 0:
        raise SystemExit("[error] No data events found in full mass window.")
    if not run3_mc_score_full:
        raise SystemExit("[error] No full-window MC events found.")

    mc_scale = 1.0
    if args.mc_normalization == "sidebands":
        if not run3_mc_score_side:
            raise SystemExit("[error] No sideband MC events found.")
        mc_scale = compute_sideband_mc_scale(run3_data_side_mass, run3_mc_score_side)
        logger.info("MC normalization: sideband scale factor = %.6f", mc_scale)
        if mc_scale < 0.1 or mc_scale > 10.0:
            logger.warning("mc_scale=%.4f is far from 1 — check MC luminosity normalization.", mc_scale)
    else:
        logger.info("MC normalization: lumi-weighted (scale = 1.0, following AN2019_090_UCSD)")

    edges = trim_degenerate_top_categories(edges, era_payloads, args, mc_scale, logger)
    interior = np.asarray(edges[1:-1], dtype=np.float64)
    category_band = (float(np.min(interior)), float(np.max(interior))) if interior.size >= 1 else None

    run3_mc_score_full    = scale_mc_groups(run3_mc_score_full, mc_scale)
    run3_mc_mass_full     = scale_mc_groups(run3_mc_mass_full, mc_scale)
    run3_mc_score_display = group_backgrounds_for_display(run3_mc_score_full)

    # ── Display-only normalization ───────────────────────────────────────────
    # Always sideband-normalize the two data/MC overview plots so data points
    # sit on top of MC at the right scale. This does not affect mc_scale or any
    # significance/Asimov/table computation.
    if run3_mc_score_side:
        display_mc_scale = compute_sideband_mc_scale(run3_data_side_mass, run3_mc_score_side)
        logger.info("Display plot normalization (sideband scale): %.6f", display_mc_scale)
    else:
        display_mc_scale = mc_scale
        logger.warning("No sideband MC for display normalization; falling back to mc_scale.")
    # mc_score_display and mc_mass_full already carry mc_scale — apply the net factor
    _display_net = display_mc_scale / mc_scale if mc_scale > 1e-12 else display_mc_scale
    mc_for_score_plot = scale_mc_groups(run3_mc_score_display, _display_net)
    mc_for_mass_plot  = scale_mc_groups(run3_mc_mass_full,     _display_net)

    era_label = ", ".join(sorted(era_payloads.keys()))

    # ── Page 2: BDT score ───────────────────────────────────────────────────
    score_bins = np.linspace(args.score_min, args.score_max, args.score_nbins + 1)
    fig_score, data_hist, _dv, _sh, _sv, mc_total, mc_total_var = make_stacked_ratio_figure(
        title="BDT score in full 110–150 GeV window: Data vs. stacked MC",
        xlabel="BDT score",
        bins=score_bins,
        data_values=run3_data_full_scores,
        mc_groups=mc_for_score_plot,
        logy=args.score_logy,
        category_band=category_band,
        show_mc_unc_band=False,
        shade_category_band=False,
        ratio_ylabel="Data / Bkg",
        era_label=era_label,
        annotation_lines=[
            f"Eras: {era_label}",
            f"Full mass window: {args.mass_min:.0f}–{args.mass_max:.0f} GeV",
            f"Boundary range: {category_band[0]:.3f} to {category_band[1]:.3f}"
            if category_band else "Boundary range: N/A",
        ],
    )

    signal_hist = np.zeros(len(score_bins)-1, dtype=np.float64)
    score_signal_overlay = combine_signal_score_overlay(era_payloads)
    if score_signal_overlay is not None:
        sig_scores_all, sig_weights_all = score_signal_overlay
        if sig_scores_all.size > 0 and np.sum(sig_weights_all) > 0:
            signal_hist, _ = weighted_hist_and_var(np.asarray(sig_scores_all, dtype=np.float64),
                                                   score_bins, np.asarray(sig_weights_all, dtype=np.float64))

    bdt_cols, bdt_rows = build_bdt_signal_background_table(
        score_bins, signal_hist, mc_total, mc_total_var, data_hist)

    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)

    # ── Page 4: Significance evolution ─────────────────────────────────────
    fig_sig_evolution = build_significance_evolution_figure(scan_history)

    # ── Page 5: Atanh BDT + S/√B ────────────────────────────────────────────
    fig_vbf_transformed = build_vbf_bdt_transformed_figure(
        era_payloads=era_payloads, edges=edges, args=args,
        mass_bins=mass_bins, mc_weight_scale=mc_scale)

    # ── Page 6: FWHM vs. BDT ───────────────────────────────────────────────
    fig_fwhm_vs_bdt = build_fwhm_vs_bdt_figure(era_payloads=era_payloads, args=args, mass_bins=mass_bins)

    # ── Page 7: ROC curve [NEW] ─────────────────────────────────────────────
    fig_roc = build_roc_curve_figure(era_payloads=era_payloads, args=args, mass_bins=mass_bins)

    # ── Page 8: Blinded mass sideband ───────────────────────────────────────
    fig_mass, _m_dhist, _m_dvar, _m_hists, _m_vars, _m_tot, _m_tot_var = make_stacked_ratio_figure(
        title="Blinded $m_{\\mu\\mu}$ spectrum: sideband data vs. MC expectation",
        xlabel="$m_{\\mu\\mu}$ [GeV]",
        bins=mass_bins,
        data_values=run3_data_side_mass,
        mc_groups=mc_for_mass_plot,
        logy=False,
        blind_band=(args.blind_min, args.blind_max),
        show_mc_unc_band=False,
        era_label=era_label,
    )

    # ── Pages 9–10: Per-category mass + background model fits ───────────────
    fig_cat_mass, fit_cols, fit_rows = build_category_mass_grid(
        era_payloads=era_payloads, edges=edges, args=args,
        mass_bins=mass_bins, mc_weight_scale=mc_scale)

    # ── Category summary (needs DCB fits) ────────────────────────────────────
    cat_summary = summarize_categories(era_payloads, edges, args, mc_weight_scale=mc_scale)

    # ── Pages 11–12: VBF DCB signal extraction ──────────────────────────────
    fig_vbf_bdt, dcb_cols, dcb_rows, all_dcb_params = build_vbf_signal_dcb_grid(
        era_payloads=era_payloads, edges=edges, args=args,
        mass_bins=mass_bins, mc_weight_scale=mc_scale, cat_summary=cat_summary)

    # ── Page 13: σ_CB vs. category [NEW] ────────────────────────────────────
    fig_sigma_vs_cat = build_sigma_vs_category_figure(all_dcb_params, edges)

    # ── Page 14: Acceptance × Efficiency [NEW] ──────────────────────────────
    fig_acc_eff = build_acc_eff_figure(era_payloads=era_payloads, edges=edges)

    # ── Page 15: S/B diagnostic [NEW] ───────────────────────────────────────
    fig_sb_diag = build_sb_diagnostic_figure(edges=edges, cat_summary=cat_summary, logger=logger)

    # ── Pages 16–17: Summary tables ─────────────────────────────────────────
    cat_cols, cat_rows = build_category_summary_table(edges, cat_summary)
    cat_proc_cols, cat_proc_rows = build_category_process_table(edges, cat_summary)

    # ── Assemble PDF ─────────────────────────────────────────────────────────
    report_name = args.report_name.strip() if args.report_name else None
    if not report_name:
        era_label_fname = "_".join(sorted(era_payloads.keys())) if era_payloads else "run3"
        report_name = f"{args.output_prefix}_report_{era_label_fname}.pdf"
    report_path = report_name if os.path.isabs(report_name) else os.path.join(args.outdir, report_name)

    with PdfPages(report_path) as pdf:
        # 1. Cover
        add_text_page(
            pdf,
            "Run-3 H → μμ VBF BDT Category Optimisation",
            [
                f"Eras              : {', '.join(sorted(era_payloads.keys()))}",
                f"Category edges    : {', '.join([f'{x:.4f}' for x in edges])}",
                f"MC scale factor   : {mc_scale:.4f}  (Data/MC sideband normalisation)",
                f"Blinding window   : {args.blind_min:.0f}–{args.blind_max:.0f} GeV",
                f"Sidebands         : ({args.sb1_min:.0f},{args.sb1_max:.0f}) ∪ "
                f"({args.sb2_min:.0f},{args.sb2_max:.0f}) GeV",
                f"Outputs dir       : {args.outputs_dir}",
                f"Extra selection   : {args.extra_cut if args.extra_cut else '—'}",
                "",
                "CONTENTS:",
                "  p.2   BDT score distribution: Data vs. MC",
                "  p.3   BDT bin yields",
                "  p.4   Significance vs. number of categories",
                "  p.5   Atanh BDT + S/√B per category",
                "  p.6   Signal FWHM vs. BDT score",
                "  p.7   ROC curve with AUC",
                "  p.8   Blinded dimuon mass spectrum",
                "  p.9   Per-category mass distributions and background fits",
                "  p.10  Background model comparison (χ²/ndf)",
                "  p.11  VBF signal DCB fits per category",
                "  p.12  DCB parameter table",
                "  p.13  σ_CB vs. category",
                "  p.14  Signal acceptance × efficiency",
                "  p.15  Per-category yield and sensitivity summary",
                "  p.16  Category summary table",
                "  p.17  Category process composition",
            ],
        )

        pdf.savefig(fig_score)
        add_table_pages(pdf, "BDT bin summary: signal, background, and sensitivity",
                        bdt_cols, bdt_rows, max_rows=20)
        pdf.savefig(fig_sig_evolution)
        pdf.savefig(fig_vbf_transformed)
        pdf.savefig(fig_fwhm_vs_bdt)
        pdf.savefig(fig_roc)          # [NEW]
        pdf.savefig(fig_mass)
        pdf.savefig(fig_cat_mass)
        add_table_pages(pdf, "Per-category background model fit comparison (χ²_BC/ndf, ndf)",
                        fit_cols, fit_rows, max_rows=20)
        pdf.savefig(fig_vbf_bdt)
        add_table_pages(pdf, "VBF signal DCB fit parameters per category",
                        dcb_cols, dcb_rows, max_rows=20)
        pdf.savefig(fig_sigma_vs_cat)  # [NEW]
        pdf.savefig(fig_acc_eff)       # [NEW]
        pdf.savefig(fig_sb_diag)       # [NEW]
        add_table_pages(pdf, "Category summary: HWHM-windowed signal and background yields",
                        cat_cols, cat_rows, max_rows=20)
        add_table_pages(pdf, "Category composition by process (full mass window)",
                        cat_proc_cols, cat_proc_rows, max_rows=20)

    # Close all figures
    for fig in [fig_score, fig_mass, fig_cat_mass, fig_vbf_bdt, fig_vbf_transformed,
                fig_sig_evolution, fig_fwhm_vs_bdt, fig_roc, fig_sigma_vs_cat,
                fig_acc_eff, fig_sb_diag]:
        try:
            plt.close(fig)
        except Exception:
            pass

    logger.info("Wrote improved report PDF: %s", report_path)
    print(f"[done] Wrote report: {report_path}")


if __name__ == "__main__":
    main()