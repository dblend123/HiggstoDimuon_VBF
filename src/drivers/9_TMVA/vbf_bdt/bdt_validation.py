#!/usr/bin/env python3
"""
VBF BDT per-era validation diagnostic script.

PURPOSE
-------
Produces a detailed per-era diagnostic PDF for a *single* era label supplied
via ``--era``.  It focuses on feature-level checks, fold-KS diagnostics,
ROC curves, and simple-cut acceptance scans.

ROLE IN THE WORKFLOW
--------------------
This script is complementary to ``bdt_event_categorization.py``, which is
the primary multi-era categorization and sideband-validation report.

Directory: vbf_bdt/

  bdt_analysis.py           — shared utility library (import only)
  bdt_event_categorization.py — MAIN script: multi-era sideband validation,
                                 category optimisation, blinded report
  bdt_validation.py         — THIS script: single-era diagnostic report
                               (KS folds, ROC, feature cross-check, …)
  bdt_figures.py            — Thin wrapper that calls bdt_event_categorization
                               functions and emits a figure-only PDF
  significance_inference.py — BDT math, chi2, significance scan utilities

SHARED HELPERS
--------------
All weight-handling, file discovery, and histogram utilities should be imported
from ``bdt_analysis``.  No duplicate implementations should exist here.
"""

import argparse
import glob
import logging
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ROOT

ROOT.gROOT.SetBatch(True)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TMVA_DIR = os.path.dirname(SCRIPT_DIR)
if TMVA_DIR not in sys.path:
    sys.path.insert(0, TMVA_DIR)

try:
    from tmva100_DataPreparation import VARIABLES as TMVA_VARIABLES
except Exception:
    TMVA_VARIABLES = [
        "z", "pt_balance", "mu1_pt_over_m", "mu2_pt_over_m", "mu1_eta", "mu2_eta",
        "y_mumu", "pt_mumu", "jet1_eta", "deta_jj", "dphi_jj", "njets",
        "pt_j1", "pt_j2", "dphi_mumu_jj", "deta_mumu_jj", "phi_CS", "costheta_CS",
    ]

# Shared pure-math utilities from the utility library
from bdt_analysis import (
    weighted_mean_std,
    weighted_corr,
    asimov_q,
    multi_bin_asimov_z,
    ks_2samp,
    compute_roc,
)


PROC_ORDER = ["DY", "DYVBF", "EWK", "TOP", "VV", "VVV", "GGH", "VBF", "BACKGROUND", "SIGNAL", "OTHER"]
PROC_COLOR = {
    "DY": "#f4a261",
    "DYVBF": "#c77dff",
    "EWK": "#6a4c93",
    "TOP": "#2a9d8f",
    "VV": "#457b9d",
    "VVV": "#5f0f40",
    "GGH": "#7f7f7f",
    "VBF": "#2a9d4b",
    "BACKGROUND": "#8093a8",
    "SIGNAL": "#2a9d4b",
    "OTHER": "#8d99ae",
}


def canonical_proc(name: str) -> str:
    s = str(name).strip().lower()
    if "dyvbf" in s:
        return "DYVBF"
    if "dy" in s or s == "z":
        return "DY"
    if "ewk" in s:
        return "EWK"
    if "tt" in s or "top" in s:
        return "TOP"
    if s in ("vvv", "www", "wwz", "wzz", "zzz"):
        return "VVV"
    if "vv" in s or "dibos" in s or s in ("ww", "wz", "zz"):
        return "VV"
    if "ggh" in s:
        return "GGH"
    if "vbf" in s and "dy" not in s:
        return "VBF"
    if "background" in s:
        return "BACKGROUND"
    if "signal" in s:
        return "SIGNAL"
    return "OTHER"


def proc_sort_key(name: str) -> int:
    cn = canonical_proc(name)
    return PROC_ORDER.index(cn) if cn in PROC_ORDER else len(PROC_ORDER)


@dataclass
class Sample:
    proc: str
    files: List[str]
    arrays: Dict[str, np.ndarray]
    is_data: bool = False


REQUESTED_VBF_FEATURES = [
    "jet1_qgl",
    "jet1_eta",
    "jet2_eta",
    "jet1_pt",
    "jet2_pt",
    "mjj",
    "dEta_jj",
    "dPhi_jj",
    "mmj_min_dEta",
    "mmj_min_dPhi",
    "pt_centrality",
    "dimuon_pt",
    "dimuon_rapidity",
    "dimuon_cos_theta_cs",
    "dimuon_phi_cs",
    "dimuon_mass_resolution",
    "mu1_pt_over_mass",
    "mu2_pt_over_mass",
    "mu1_eta",
    "mu2_eta",
    "N_jets",
    "year",
]

REQUESTED_ALIAS_MAP = {
    "jet1_qgl": [("jet1_qgl", "alias", "")],
    "jet1_eta": [("jet1_eta", "alias", ""), ("eta_j1", "alias", "")],
    "jet2_eta": [("jet2_eta", "alias", ""), ("eta_j2", "alias", ""), ("y_j2", "proxy", "rapidity, not eta")],
    "jet1_pt": [("jet1_pt", "alias", ""), ("pt_j1", "alias", "")],
    "jet2_pt": [("jet2_pt", "alias", ""), ("pt_j2", "alias", "")],
    "mjj": [("mjj", "alias", ""), ("m_jj", "alias", ""), ("dijet_mass", "alias", "")],
    "dEta_jj": [("dEta_jj", "alias", ""), ("deta_jj", "alias", "")],
    "dPhi_jj": [("dPhi_jj", "alias", ""), ("dphi_jj", "alias", "")],
    "mmj_min_dEta": [("mmj_min_dEta", "alias", ""), ("deta_mumu_jj", "proxy", "dimuon-dijet deta, not min(mu,jet)")],
    "mmj_min_dPhi": [("mmj_min_dPhi", "alias", ""), ("dphi_mumu_jj", "proxy", "dimuon-dijet dphi, not min(mu,jet)")],
    "pt_centrality": [("pt_centrality", "alias", ""), ("z", "proxy", "zeppenfeld variable, not pt centrality")],
    "dimuon_pt": [("dimuon_pt", "alias", ""), ("pt_mumu", "alias", "")],
    "dimuon_rapidity": [("dimuon_rapidity", "alias", ""), ("y_mumu", "alias", "")],
    "dimuon_cos_theta_cs": [("dimuon_cos_theta_cs", "alias", ""), ("costheta_CS", "alias", "")],
    "dimuon_phi_cs": [("dimuon_phi_cs", "alias", ""), ("phi_CS", "alias", "")],
    "dimuon_mass_resolution": [("dimuon_mass_resolution", "alias", ""), ("mass_resolution", "alias", "")],
    "mu1_pt_over_mass": [("mu1_pt_over_mass", "alias", ""), ("mu1_pt_over_m", "alias", "")],
    "mu2_pt_over_mass": [("mu2_pt_over_mass", "alias", ""), ("mu2_pt_over_m", "alias", "")],
    "mu1_eta": [("mu1_eta", "alias", "")],
    "mu2_eta": [("mu2_eta", "alias", "")],
    "N_jets": [("N_jets", "alias", ""), ("njets", "alias", "")],
    "year": [("year", "derived", "derived from --era")],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="VBF + BDT validation workflow")
    p.add_argument("--era", required=True, help="Era label, e.g. 2024 or 22023")
    p.add_argument("--outputs-dir", default=os.path.join("src", "drivers", "9_TMVA", "Outputs"),
                   help="Directory with scored ROOT files")
    p.add_argument("--output-dir", default=os.path.join("src", "drivers", "9_TMVA", "Outputs", "VBFValidation"),
                   help="Output directory for plots/tables")
    p.add_argument("--treename", default="mva")

    p.add_argument("--data", action="append", default=[],
                   help="Comma-separated data ROOT files")
    p.add_argument("--mc", action="append", default=[],
                   help="PROC:file1.root,file2.root (repeat)")

    p.add_argument("--score-branch", default="auto")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--mass-branch", default="auto")
    p.add_argument("--njets-branch", default="auto")
    p.add_argument("--mjj-branch", default="auto")
    p.add_argument("--deta-jj-branch", default="auto")
    p.add_argument("--dphi-jj-branch", default="auto")
    p.add_argument("--ptj1-branch", default="auto")
    p.add_argument("--fold-branch", default="fold")
    p.add_argument("--njet-gap-branch", default="auto",
                   help="Optional central-gap jet count branch for simple-cut veto")

    p.add_argument("--model-file", default="", help="TMVA model file for fallback scoring if score branch is absent")
    p.add_argument("--model-name", default="myBDT")

    p.add_argument("--signal-procs", default="VBF",
                   help="Comma-separated MC process names treated as signal in ROC/purity")
    p.add_argument("--vbf-procs", default="VBF",
                   help="Comma-separated MC process names treated as VBF in VBF-specific summaries")

    p.add_argument("--selection", default="", help="Extra ROOT filter expression")

    p.add_argument("--sr-min", type=float, default=110.0)
    p.add_argument("--sr-max", type=float, default=150.0)
    p.add_argument("--sb1-min", type=float, default=100.0)
    p.add_argument("--sb1-max", type=float, default=115.0)
    p.add_argument("--sb2-min", type=float, default=135.0)
    p.add_argument("--sb2-max", type=float, default=150.0)

    p.add_argument("--score-nbins", type=int, default=40)
    p.add_argument("--score-min", type=float, default=0.0)
    p.add_argument("--score-max", type=float, default=1.0)

    p.add_argument("--simple-bdt-thresholds", default="0.5,0.7,0.9")
    p.add_argument("--acceptance-bins", default="0,0.2,0.4,0.6,0.8,0.9,0.95,0.99,1.0")

    p.add_argument("--category-fine-bins", type=int, default=100)
    p.add_argument("--category-delta-z", type=float, default=0.01,
                   help="Stop merging when relative dZ from adding next lower fine-bin is below this")
    p.add_argument("--category-min-raw", type=int, default=15)
    p.add_argument("--category-min-sig", type=float, default=0.05)

    p.add_argument("--feature-importance-file", default="",
                   help="Optional feature-importance txt file")
    p.add_argument("--fold-importance-files", default="",
                   help="Comma-separated fold-wise feature-importance txt files")

    # Fold KS diagnostic thresholds (CMS-style loose/nominal/tight defaults)
    p.add_argument(
        "--ks-threshold-labels",
        default="Run2-loose,Run2-nominal,Run2-tight",
        help="Comma-separated labels for fold-KS threshold profiles",
    )
    p.add_argument(
        "--ks-p-thresholds",
        default="0.01,0.05,0.10",
        help="Comma-separated minimum p-values for each KS profile",
    )
    p.add_argument(
        "--ks-d-thresholds",
        default="0.15,0.10,0.07",
        help="Comma-separated maximum KS-D values for each KS profile",
    )

    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args()


def setup_logger(level: str) -> logging.Logger:
    logger = logging.getLogger("vbf_bdt_validation")
    logger.setLevel(getattr(logging, level))
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.handlers.clear()
    logger.addHandler(handler)
    return logger


def parse_csv(text: str) -> List[str]:
    return [x.strip() for x in str(text).split(",") if x.strip()]


def parse_csv_floats(text: str) -> List[float]:
    out = []
    for x in parse_csv(text):
        out.append(float(x))
    return out


def parse_mc_entries(entries: Sequence[str]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for item in entries:
        if ":" not in item:
            raise SystemExit(f"[error] --mc expects PROC:file1,file2, got '{item}'")
        proc, files = item.split(":", 1)
        flist = [f for f in parse_csv(files) if os.path.exists(f)]
        if flist:
            out.setdefault(proc, []).extend(flist)
    return out


def discover_auto_inputs(outputs_dir: str, era: str, logger: logging.Logger) -> Tuple[List[str], Dict[str, List[str]]]:
    data_files: List[str] = []
    mc_map: Dict[str, List[str]] = {}
    tiny_files: List[str] = []

    patterns = [
        os.path.join(outputs_dir, f"{era}_scored_*.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_*.root"),
        os.path.join(outputs_dir, f"{era}_data_scored.root"),
        os.path.join(outputs_dir, era, f"{era}_data_scored.root"),
    ]

    seen = set()
    files = []
    for pat in patterns:
        for path in glob.glob(pat):
            if path not in seen:
                seen.add(path)
                files.append(path)

    for path in files:
        # Skip obviously broken placeholder files before ROOT touches them.
        try:
            fsize = os.path.getsize(path)
        except OSError:
            logger.warning("Skipping file with failed stat(): %s", path)
            continue
        if fsize < 1024:
            tiny_files.append(path)
            continue

        base = os.path.basename(path)
        low = base.lower()
        if "data_scored" in low or low.startswith(f"{era.lower()}_data"):
            data_files.append(path)
            continue

        if "_scored_" in base:
            proc = base.split("_scored_", 1)[1].replace(".root", "")
            proc_u = proc.upper()
            if proc_u == "ALL":
                mc_map.setdefault("BACKGROUND", []).append(path)
            elif proc_u == "BACKGROUND":
                mc_map.setdefault("BACKGROUND", []).append(path)
            elif proc_u == "SIGNAL":
                mc_map.setdefault("SIGNAL", []).append(path)
            elif proc_u == "TEST":
                logger.info("Skipping non-physics scored test file: %s", path)
            else:
                mc_map.setdefault(proc, []).append(path)

    if tiny_files:
        if data_files or mc_map:
            logger.info(
                "Ignored %d tiny file(s) (<1KB); using valid discovered inputs instead. Example: %s",
                len(tiny_files),
                tiny_files[0],
            )
        else:
            for path in tiny_files:
                logger.warning("Skipping tiny file (<1KB), likely corrupt/incomplete: %s", path)

    logger.info("Auto-discovery: data files=%d, MC groups=%d", len(data_files), len(mc_map))
    return data_files, mc_map


def list_branches(path: str, treename: str) -> List[str]:
    if not os.path.exists(path):
        return []
    try:
        if os.path.getsize(path) < 1024:
            return []
    except OSError:
        return []

    f = None
    try:
        f = ROOT.TFile.Open(path)
    except Exception:
        return []

    if not f or f.IsZombie():
        try:
            if f:
                f.Close()
        except Exception:
            pass
        return []

    try:
        t = f.Get(treename)
        if not t:
            f.Close()
            return []
        out = [br.GetName() for br in t.GetListOfBranches()]
        f.Close()
        return out
    except Exception:
        try:
            f.Close()
        except Exception:
            pass
        return []


def union_branches(file_list: Sequence[str], treename: str, logger: Optional[logging.Logger] = None) -> List[str]:
    names = set()
    skipped = 0
    for path in file_list:
        cols = list_branches(path, treename)
        if not cols:
            skipped += 1
            if logger is not None:
                logger.warning("Skipping unreadable/invalid file during branch scan: %s", path)
            continue
        names.update(cols)
    if skipped and logger is not None:
        logger.warning("Branch scan skipped %d unreadable/invalid file(s)", skipped)
    return sorted(names)


def resolve_branch(preferred: str, available: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    avail = set(available)
    if preferred != "auto":
        return preferred if preferred in avail else None
    for c in candidates:
        if c in avail:
            return c
    return None


def resolve_branch_map(args: argparse.Namespace, available: Sequence[str]) -> Dict[str, Optional[str]]:
    branch_map = {
        "score": resolve_branch(args.score_branch, available, ["bdt_score", "score", "bdt"]),
        "weight": resolve_branch(args.weight_branch, available, ["weight", "eventWeight", "evt_w"]),
        "mass": resolve_branch(args.mass_branch, available, ["m_mumu", "m", "mass", "mll"]),
        "njets": resolve_branch(args.njets_branch, available, ["njets", "N_jets", "nJet", "n_jets"]),
        "mjj": resolve_branch(args.mjj_branch, available, ["mjj", "m_jj", "dijet_mass", "mjj_jj"]),
        "deta_jj": resolve_branch(args.deta_jj_branch, available, ["deta_jj", "dEta_jj", "deltaEta_jj"]),
        "dphi_jj": resolve_branch(args.dphi_jj_branch, available, ["dphi_jj", "dPhi_jj", "deltaPhi_jj"]),
        "pt_j1": resolve_branch(args.ptj1_branch, available, ["pt_j1", "jet1_pt", "j1_pt"]),
        "pt_j2": resolve_branch("auto", available, ["pt_j2", "jet2_pt", "j2_pt"]),
        "fold": resolve_branch(args.fold_branch, available, ["fold"]),
        "njet_gap": resolve_branch(args.njet_gap_branch, available, ["n_jets_between", "nJets_between", "n_jets_in_gap", "njet_gap"]),
    }
    return branch_map


def add_note_page(pdf: PdfPages, title: str, lines: Sequence[str]) -> None:
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.axis("off")
    ax.set_title(title)
    ax.text(0.01, 0.98, "\n".join(lines), transform=ax.transAxes, va="top", family="monospace", fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


def safe_legend(ax, **kwargs) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if handles and any(str(lbl).strip() and not str(lbl).startswith("_") for lbl in labels):
        ax.legend(**kwargs)


def build_ks_threshold_profiles(args: argparse.Namespace, logger: logging.Logger) -> List[Tuple[str, float, float]]:
    labels = parse_csv(args.ks_threshold_labels)
    pvals = parse_csv_floats(args.ks_p_thresholds)
    dvals = parse_csv_floats(args.ks_d_thresholds)

    n = min(len(labels), len(pvals), len(dvals))
    if n == 0:
        logger.warning("Invalid KS profile inputs; using fallback profile p>=0.05, D<=0.10")
        return [("Nominal", 0.05, 0.10)]

    if not (len(labels) == len(pvals) == len(dvals)):
        logger.warning(
            "KS profile lengths mismatch labels/p/d = %d/%d/%d; truncating to %d",
            len(labels),
            len(pvals),
            len(dvals),
            n,
        )

    out: List[Tuple[str, float, float]] = []
    for i in range(n):
        out.append((str(labels[i]).strip(), max(float(pvals[i]), 0.0), max(float(dvals[i]), 0.0)))
    return out


def ks_profile_status(d_val: float, p_val: float, p_min: float, d_max: float) -> str:
    if p_val >= p_min and d_val <= d_max:
        return "PASS"
    if p_val >= 0.5 * p_min and d_val <= 1.3 * d_max:
        return "WARN"
    return "FAIL"


# The following pure-math utilities are imported from bdt_analysis above:
# weighted_mean_std, weighted_corr, asimov_q, multi_bin_asimov_z,
# ks_2samp, compute_roc.


def interpolate_tpr_at_fpr(fpr: np.ndarray, tpr: np.ndarray, target: float) -> float:
    if fpr.size == 0:
        return float("nan")
    return float(np.interp(target, fpr, tpr, left=tpr[0], right=tpr[-1]))


def weighted_hist(values: np.ndarray, bins: np.ndarray, weights: Optional[np.ndarray] = None) -> np.ndarray:
    hist, _ = np.histogram(values, bins=bins, weights=weights)
    return hist.astype(np.float64)


def parse_feature_importance_txt(path: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if ":" not in line:
                continue
            name, val = line.split(":", 1)
            name = name.strip()
            sval = val.strip().replace("%", "")
            try:
                out[name] = float(sval)
            except ValueError:
                continue
    return out


def load_tmva_model(path: str, name: str):
    if not path:
        return None
    ROOT.gSystem.Load("libTMVA")
    if ROOT.gSystem.AccessPathName(path, ROOT.kReadPermission):
        raise RuntimeError(f"Model file not readable: {path}")
    return ROOT.TMVA.Experimental.RBDT[""](name, path)


def load_sample(
    proc: str,
    files: Sequence[str],
    treename: str,
    branch_map: Dict[str, Optional[str]],
    selection: str,
    model,
    logger: logging.Logger,
    is_data: bool = False,
) -> Sample:
    fields = ["score", "mass", "weight", "njets", "mjj", "deta_jj", "dphi_jj", "pt_j1", "pt_j2", "fold", "njet_gap"]
    collected = {k: [] for k in fields}

    for path in files:
        cols = set(list_branches(path, treename))
        if not cols:
            logger.warning("Skipping unreadable tree: %s", path)
            continue

        use_map = {}
        for key in fields:
            br = branch_map.get(key)
            if br and br in cols:
                use_map[key] = br

        need_model_score = "score" not in use_map and model is not None

        read_branches = set(use_map.values())
        if need_model_score:
            for v in TMVA_VARIABLES:
                if v in cols:
                    read_branches.add(v)
                else:
                    raise RuntimeError(
                        f"Cannot compute model score in {path}: missing feature branch '{v}'"
                    )

        if not read_branches:
            logger.warning("No readable branches in %s", path)
            continue

        try:
            df = ROOT.RDataFrame(treename, path)
            if selection:
                df = df.Filter(selection)
            arrs = df.AsNumpy(sorted(read_branches))
        except Exception as exc:
            logger.warning("Skipping file due read failure: %s (%s)", path, exc)
            continue

        if len(arrs) == 0:
            continue

        n = len(next(iter(arrs.values())))

        score_vals = None
        if "score" in use_map:
            score_vals = np.asarray(arrs[use_map["score"]], dtype=np.float64)
        elif need_model_score:
            x = np.vstack([np.asarray(arrs[v], dtype=np.float32) for v in TMVA_VARIABLES]).T
            score_vals = np.ravel(model.Compute(x).astype(np.float64, copy=False))

        for key in fields:
            if key == "score":
                vals = score_vals if score_vals is not None else np.full(n, np.nan, dtype=np.float64)
            elif key in use_map:
                vals = np.asarray(arrs[use_map[key]], dtype=np.float64)
            elif key == "weight":
                vals = np.ones(n, dtype=np.float64)
            else:
                vals = np.full(n, np.nan, dtype=np.float64)
            collected[key].append(vals)

    merged = {}
    for key in fields:
        merged[key] = np.concatenate(collected[key]) if collected[key] else np.array([], dtype=np.float64)

    merged["weight"] = sanitize_weights(merged["weight"])

    return Sample(proc=proc, files=list(files), arrays=merged, is_data=is_data)


def merge_samples_by_proc(samples: Sequence[Sample]) -> Dict[str, Dict[str, np.ndarray]]:
    groups: Dict[str, Dict[str, List[np.ndarray]]] = {}
    for s in samples:
        cp = canonical_proc(s.proc)
        if cp not in groups:
            groups[cp] = {k: [] for k in s.arrays.keys()}
        for key, vals in s.arrays.items():
            groups[cp][key].append(vals)
    merged: Dict[str, Dict[str, np.ndarray]] = {}
    for cp, mp in groups.items():
        merged[cp] = {
            k: (np.concatenate(vlist) if vlist else np.array([], dtype=np.float64))
            for k, vlist in mp.items()
        }
    return merged


def combine_group_arrays(group_map: Dict[str, Dict[str, np.ndarray]], proc_set: Sequence[str]) -> Dict[str, np.ndarray]:
    out_keys = ["score", "mass", "weight", "njets", "mjj", "deta_jj", "dphi_jj", "pt_j1", "pt_j2", "fold", "njet_gap"]
    out = {k: [] for k in out_keys}
    canon_set = {canonical_proc(p) for p in proc_set}
    for proc, arrs in group_map.items():
        if proc not in canon_set:
            continue
        for k in out_keys:
            out[k].append(arrs.get(k, np.array([], dtype=np.float64)))
    merged = {}
    for k in out_keys:
        merged[k] = np.concatenate(out[k]) if out[k] else np.array([], dtype=np.float64)
    return merged


def all_mc_arrays(group_map: Dict[str, Dict[str, np.ndarray]]) -> Dict[str, np.ndarray]:
    return combine_group_arrays(group_map, list(group_map.keys()))


def make_mass_masks(arr: Dict[str, np.ndarray], args: argparse.Namespace) -> Dict[str, np.ndarray]:
    m = arr["mass"]
    valid = np.isfinite(m)
    sr = valid & (m >= args.sr_min) & (m < args.sr_max)
    sb = (valid & (m >= args.sb1_min) & (m < args.sb1_max)) | (valid & (m >= args.sb2_min) & (m < args.sb2_max))
    zpeak = valid & (m >= 70.0) & (m < 110.0)
    highsb = valid & (m > 150.0)
    return {"sr": sr, "sb": sb, "zpeak": zpeak, "highsb": highsb}


def crosscheck_requested_features(
    available_branches: Sequence[str],
    era: str,
    output_dir: str,
    logger: logging.Logger,
) -> List[Tuple[str, str, str, str]]:
    avail = set(available_branches)
    rows: List[Tuple[str, str, str, str]] = []

    for req in REQUESTED_VBF_FEATURES:
        if req == "year":
            rows.append((req, "derived", f"year={era}", "from --era"))
            continue

        if req in avail:
            rows.append((req, "exact", req, ""))
            continue

        found = False
        for alias, kind, note in REQUESTED_ALIAS_MAP.get(req, []):
            if alias in avail and kind in ("alias", "proxy"):
                rows.append((req, kind, alias, note))
                found = True
                break
        if not found:
            rows.append((req, "missing", "-", "not found in tree branches"))

    miss = sum(1 for r in rows if r[1] == "missing")
    proxy = sum(1 for r in rows if r[1] == "proxy")
    logger.info("Feature cross-check: missing=%d, proxy=%d", miss, proxy)
    return rows


def build_stack_inputs(
    group_map: Dict[str, Dict[str, np.ndarray]],
    mask_name: str,
    args: argparse.Namespace,
) -> List[Tuple[str, np.ndarray, np.ndarray]]:
    out = []
    for proc in sorted(group_map.keys(), key=proc_sort_key):
        arr = group_map[proc]
        masks = make_mass_masks(arr, args)
        m = masks[mask_name]
        score = arr["score"][m]
        w = arr["weight"][m]
        if score.size == 0:
            continue
        out.append((proc, score, w))
    return out


def stack_and_ratio(
    ax_top,
    ax_bot,
    bins: np.ndarray,
    title: str,
    mc_groups: Sequence[Tuple[str, np.ndarray, np.ndarray]],
    data_scores: Optional[np.ndarray],
    logy: bool = True,
) -> None:
    hist_vals = []
    colors = []
    labels = []
    mc_sum = np.zeros(len(bins) - 1, dtype=np.float64)

    for proc, score, w in mc_groups:
        h = weighted_hist(score, bins=bins, weights=w)
        if np.sum(h) <= 0:
            continue
        hist_vals.append(score)
        colors.append(PROC_COLOR.get(proc, PROC_COLOR["OTHER"]))
        labels.append(proc)
        mc_sum += h

    if hist_vals:
        ax_top.hist(
            hist_vals,
            bins=bins,
            weights=[g[2] for g in mc_groups if weighted_hist(g[1], bins=bins, weights=g[2]).sum() > 0],
            stacked=True,
            color=colors,
            label=labels,
            alpha=0.85,
            histtype="stepfilled",
        )

    if data_scores is not None and data_scores.size > 0:
        h_data = weighted_hist(data_scores, bins=bins, weights=None)
        centers = 0.5 * (bins[:-1] + bins[1:])
        err = np.sqrt(np.maximum(h_data, 1.0))
        ax_top.errorbar(centers, h_data, yerr=err, fmt="o", color="black", label="Data")

        if ax_bot is not None:
            ratio = np.divide(h_data, mc_sum, out=np.zeros_like(h_data, dtype=np.float64), where=mc_sum > 0)
            ratio_err = np.divide(err, mc_sum, out=np.zeros_like(err), where=mc_sum > 0)
            ax_bot.axhline(1.0, color="gray", linestyle="--", linewidth=1)
            ax_bot.errorbar(centers, ratio, yerr=ratio_err, fmt="o", color="black")
            ax_bot.set_ylabel("Data/MC")
            ax_bot.set_ylim(0.4, 1.6)
            ax_bot.grid(alpha=0.3)
            ax_bot.set_xlabel("BDT score")
    else:
        if ax_bot is not None:
            ax_bot.axis("off")

    ax_top.set_title(title)
    ax_top.set_ylabel("Events / bin")
    if logy:
        ax_top.set_yscale("log")
    ax_top.grid(alpha=0.3)
    safe_legend(ax_top, fontsize=8, ncol=2)


def profile_mean(x: np.ndarray, y: np.ndarray, w: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    xc = []
    ym = []
    ye = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = np.isfinite(x) & np.isfinite(y) & (x >= lo) & (x < hi)
        if not np.any(m):
            continue
        xx = x[m]
        yy = y[m]
        ww = w[m]
        sw = ww.sum()
        if sw <= 0:
            continue
        mu = np.sum(ww * yy) / sw
        var = np.sum(ww * (yy - mu) ** 2) / sw
        err = math.sqrt(max(var, 0.0) / max(np.count_nonzero(m), 1))
        xc.append(0.5 * (lo + hi))
        ym.append(mu)
        ye.append(err)
    return np.asarray(xc), np.asarray(ym), np.asarray(ye)


def safe_auc(fpr: np.ndarray, tpr: np.ndarray) -> float:
    if fpr.size == 0:
        return float("nan")
    return float(np.trapz(tpr, fpr))


def optimize_categories_high_to_low(
    sig_scores: np.ndarray,
    sig_w: np.ndarray,
    bkg_scores: np.ndarray,
    bkg_w: np.ndarray,
    fine_bins: int,
    delta_z: float,
    min_raw: int,
    min_sig: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = np.linspace(0.0, 1.0, fine_bins + 1)
    s_hist = weighted_hist(sig_scores, bins=edges, weights=sig_w)
    b_hist = weighted_hist(bkg_scores, bins=edges, weights=bkg_w)
    raw_hist = weighted_hist(np.concatenate([sig_scores, bkg_scores]), bins=edges, weights=None)

    out_edges = [1.0]
    i_hi = fine_bins - 1

    while i_hi >= 0:
        s_acc = 0.0
        b_acc = 0.0
        raw_acc = 0.0
        i = i_hi

        while i >= 0:
            s_new = s_acc + s_hist[i]
            b_new = b_acc + b_hist[i]
            raw_new = raw_acc + raw_hist[i]

            z_old = math.sqrt(max(asimov_q(s_acc, b_acc), 0.0))
            z_new = math.sqrt(max(asimov_q(s_new, b_new), 0.0))

            if z_old > 1e-12:
                rel_gain = (z_new - z_old) / z_old
            else:
                rel_gain = 1.0

            s_acc, b_acc, raw_acc = s_new, b_new, raw_new
            i -= 1

            if raw_acc >= min_raw and s_acc >= min_sig and rel_gain < delta_z:
                break

        low_edge = edges[max(i + 1, 0)]
        out_edges.append(float(low_edge))
        i_hi = i

    out_edges = sorted(set(out_edges))
    if out_edges[0] > 0.0:
        out_edges = [0.0] + out_edges
    if out_edges[-1] < 1.0:
        out_edges = out_edges + [1.0]

    out_edges_arr = np.asarray(out_edges, dtype=np.float64)
    s_bins = weighted_hist(sig_scores, bins=out_edges_arr, weights=sig_w)
    b_bins = weighted_hist(bkg_scores, bins=out_edges_arr, weights=bkg_w)
    return out_edges_arr, s_bins, b_bins


def main() -> None:
    args = parse_args()
    logger = setup_logger(args.log_level)

    os.makedirs(args.output_dir, exist_ok=True)

    data_files = []
    mc_map: Dict[str, List[str]] = {}

    if args.data or args.mc:
        for item in args.data:
            data_files.extend([p for p in parse_csv(item) if os.path.exists(p)])
        mc_map = parse_mc_entries(args.mc)
    else:
        data_files, mc_map = discover_auto_inputs(args.outputs_dir, args.era, logger)

    if not mc_map:
        raise SystemExit("[error] No MC files found. Provide --mc or use a valid --outputs-dir/--era.")

    all_files = data_files + [p for files in mc_map.values() for p in files]
    available = union_branches(all_files, args.treename, logger=logger)
    if not available:
        raise SystemExit("[error] Could not read branches from input trees.")

    branch_map = resolve_branch_map(args, available)

    model = None
    if branch_map["score"] is None and args.model_file:
        model = load_tmva_model(args.model_file, args.model_name)
        logger.info("Loaded model for fallback scoring: %s", args.model_file)

    if branch_map["score"] is None and model is None:
        raise SystemExit("[error] BDT score branch not found and no --model-file was provided.")

    if branch_map["mass"] is None:
        raise SystemExit("[error] Could not resolve dimuon mass branch.")

    # 0) Cross-check requested variables
    feature_rows = crosscheck_requested_features(available, args.era, args.output_dir, logger)

    # Load samples
    data_sample = None
    if data_files:
        data_sample = load_sample("DATA", data_files, args.treename, branch_map, args.selection, model, logger, is_data=True)

    mc_samples: List[Sample] = []
    for proc, files in mc_map.items():
        s = load_sample(proc, files, args.treename, branch_map, args.selection, model, logger, is_data=False)
        if s.arrays["score"].size == 0:
            logger.warning("Skipping empty sample: %s", proc)
            continue
        mc_samples.append(s)

    if not mc_samples:
        raise SystemExit("[error] No usable MC samples after loading.")

    mc_groups = merge_samples_by_proc(mc_samples)

    signal_set = [canonical_proc(x) for x in parse_csv(args.signal_procs)]
    vbf_set = [canonical_proc(x) for x in parse_csv(args.vbf_procs)]

    available_groups = set(mc_groups.keys())
    if not any(sig in available_groups for sig in signal_set):
        if "SIGNAL" in available_groups:
            logger.warning(
                "No groups matched --signal-procs='%s'. Falling back to discovered SIGNAL group.",
                args.signal_procs,
            )
            signal_set = ["SIGNAL"]
        else:
            logger.warning(
                "No groups matched --signal-procs='%s'. Available groups: %s",
                args.signal_procs,
                ", ".join(sorted(available_groups)),
            )

    if not any(vbf in available_groups for vbf in vbf_set):
        if any(sig in available_groups for sig in signal_set):
            logger.warning(
                "No groups matched --vbf-procs='%s'. Falling back to signal groups: %s",
                args.vbf_procs,
                ", ".join(signal_set),
            )
            vbf_set = list(dict.fromkeys(signal_set))
        elif "SIGNAL" in available_groups:
            logger.warning(
                "No groups matched --vbf-procs='%s'. Falling back to discovered SIGNAL group.",
                args.vbf_procs,
            )
            vbf_set = ["SIGNAL"]
        else:
            logger.warning(
                "No groups matched --vbf-procs='%s'. Available groups: %s",
                args.vbf_procs,
                ", ".join(sorted(available_groups)),
            )

    sig_arr = combine_group_arrays(mc_groups, signal_set)
    bkg_arr = combine_group_arrays(mc_groups, [p for p in mc_groups.keys() if p not in set(signal_set)])
    all_arr = all_mc_arrays(mc_groups)

    if sig_arr["score"].size == 0 or bkg_arr["score"].size == 0:
        logger.warning("Signal or background group is empty with current --signal-procs setting.")

    sm = make_mass_masks(sig_arr, args)["sr"]
    bm = make_mass_masks(bkg_arr, args)["sr"]
    s_sc = sig_arr["score"][sm]
    s_w = sig_arr["weight"][sm]
    b_sc = bkg_arr["score"][bm]
    b_w = bkg_arr["weight"][bm]

    score_bins = np.linspace(args.score_min, args.score_max, args.score_nbins + 1)

    pdf_path = os.path.join(args.output_dir, f"{args.era}_vbf_bdt_validation.pdf")

    summary_lines = []
    summary_lines.append(f"Era: {args.era}")
    summary_lines.append(f"Inputs: data_files={len(data_files)} mc_groups={len(mc_groups)}")
    summary_lines.append(f"Signal groups: {', '.join(signal_set) if signal_set else '(none)'}")
    summary_lines.append(f"VBF groups: {', '.join(vbf_set) if vbf_set else '(none)'}")
    summary_lines.append(f"SR weighted yields: S={float(np.sum(s_w)):.3f}, B={float(np.sum(b_w)):.3f}")

    ks_profiles = build_ks_threshold_profiles(args, logger)
    profile_desc = "; ".join([f"{name}: p>={pmin:.3g}, D<={dmax:.3g}" for name, pmin, dmax in ks_profiles])
    summary_lines.append(f"KS profiles: {profile_desc}")

    with PdfPages(pdf_path) as pdf:
        # Header page
        add_note_page(
            pdf,
            f"VBF BDT Validation - {args.era}",
            [
                "Standalone per-era diagnostic script: bdt_validation.py",
                "Logic references: ScanYields / PlotScored / MassSculpt / KS_perfold / DataMC",
                f"Resolved branches: {branch_map}",
            ],
        )

        # Feature cross-check table page
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis("off")
        ax.set_title("Requested VBF feature cross-check")
        col_labels = ["requested", "status", "resolved_as", "notes"]
        table_data = [list(r) for r in feature_rows]
        tbl = ax.table(cellText=table_data, colLabels=col_labels, loc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8)
        tbl.scale(1, 1.3)
        pdf.savefig(fig)
        plt.close(fig)

        # 1) Feature importance ranking
        fi_items = None
        fi_path = ""
        if args.feature_importance_file and os.path.exists(args.feature_importance_file):
            fi_path = args.feature_importance_file
        else:
            pat = os.path.join(args.outputs_dir, "Training", args.era, f"{args.era}_featureimportance_*.txt")
            found = sorted(glob.glob(pat))
            if found:
                fi_path = found[0]

        if fi_path:
            fi_map = parse_feature_importance_txt(fi_path)
            fi_items = sorted(fi_map.items(), key=lambda x: x[1], reverse=True)
            topn = fi_items[:20]
            fig, ax = plt.subplots(figsize=(11, 8))
            names = [x[0] for x in topn][::-1]
            vals = [x[1] for x in topn][::-1]
            ax.barh(names, vals, color="#4c78a8")
            ax.set_xlabel("Feature importance (gain fraction or %)")
            ax.set_title("A) Feature importance ranking")
            ax.grid(alpha=0.3)
            ax.text(0.01, 0.01, f"source: {fi_path}", transform=ax.transAxes, fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)
            summary_lines.append(f"Feature importance source: {fi_path}")
        else:
            logger.info("Skipping feature-importance ranking page: no feature-importance file found.")

        # Optional fold-wise stability
        fold_files = [p for p in parse_csv(args.fold_importance_files) if os.path.exists(p)]
        if len(fold_files) >= 2:
            fold_maps = [parse_feature_importance_txt(p) for p in fold_files]
            feats = sorted(set().union(*[set(m.keys()) for m in fold_maps]))
            means = []
            stds = []
            for feat in feats:
                vals = [m.get(feat, 0.0) for m in fold_maps]
                means.append(float(np.mean(vals)))
                stds.append(float(np.std(vals)))
            order = np.argsort(means)[::-1][:15]
            fig, ax = plt.subplots(figsize=(11, 8))
            y = np.arange(len(order))
            ax.barh(y, [means[i] for i in order], xerr=[stds[i] for i in order], color="#72b7b2")
            ax.set_yticks(y)
            ax.set_yticklabels([feats[i] for i in order])
            ax.invert_yaxis()
            ax.set_xlabel("mean importance")
            ax.set_title("6) Feature importance stability across folds")
            ax.grid(alpha=0.3)
            pdf.savefig(fig)
            plt.close(fig)

        # BDT score vs jet multiplicity
        if np.isfinite(all_arr["njets"]).any():
            dy_like = combine_group_arrays(mc_groups, ["DY", "DYVBF"])
            vbf_like = combine_group_arrays(mc_groups, vbf_set)
            dy_masks = make_mass_masks(dy_like, args)
            vbf_masks = make_mass_masks(vbf_like, args)
            fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
            any_curve = False
            cats = [
                ("Njets = 2", lambda a: np.isfinite(a["njets"]) & (a["njets"] >= 2) & (a["njets"] < 3)),
                ("Njets = 3", lambda a: np.isfinite(a["njets"]) & (a["njets"] >= 3) & (a["njets"] < 4)),
                ("Njets >= 4", lambda a: np.isfinite(a["njets"]) & (a["njets"] >= 4)),
            ]
            for ax, (title, fn) in zip(axes, cats):
                m_dy = dy_masks["sr"] & fn(dy_like)
                m_vbf = vbf_masks["sr"] & fn(vbf_like)

                s_dy = dy_like["score"][m_dy]
                w_dy = dy_like["weight"][m_dy]
                s_vbf = vbf_like["score"][m_vbf]
                w_vbf = vbf_like["weight"][m_vbf]

                hd = weighted_hist(s_dy, bins=score_bins, weights=w_dy)
                hv = weighted_hist(s_vbf, bins=score_bins, weights=w_vbf)
                centers = 0.5 * (score_bins[:-1] + score_bins[1:])

                plotted = False
                if hd.sum() > 0:
                    hd = hd / hd.sum()
                    ax.step(centers, hd, where="mid", color=PROC_COLOR["DY"], linewidth=2, label="DY-like")
                    plotted = True
                    any_curve = True
                if hv.sum() > 0:
                    hv = hv / hv.sum()
                    ax.step(
                        centers,
                        hv,
                        where="mid",
                        color=PROC_COLOR["VBF"],
                        linewidth=2.4,
                        linestyle="--",
                        label="VBF-like",
                    )
                    plotted = True
                    any_curve = True

                if not plotted:
                    ax.text(0.5, 0.5, "No events", transform=ax.transAxes, ha="center", va="center", fontsize=10)
                ax.set_title(title)
                ax.set_xlabel("BDT score")
                ax.grid(alpha=0.3)
                safe_legend(ax, fontsize=8)

            if any_curve:
                axes[0].set_ylabel("Normalized entries")
                safe_legend(axes[0])
                fig.suptitle("BDT score vs jet multiplicity (SR)")
                pdf.savefig(fig)
            else:
                logger.info("Skipping jet-multiplicity page: no non-empty DY/VBF histograms.")
            plt.close(fig)
        else:
            logger.info("Skipping jet-multiplicity page: njets branch not available.")

        # mjj vs score 2D S/(S+B)
        if np.isfinite(sig_arr["mjj"]).any() and np.isfinite(bkg_arr["mjj"]).any():
            s2d, xedges, yedges = np.histogram2d(
                sig_arr["score"][sm],
                sig_arr["mjj"][sm],
                bins=[np.linspace(0.0, 1.0, 41), np.linspace(0.0, 2000.0, 41)],
                weights=sig_arr["weight"][sm],
            )
            b2d, _, _ = np.histogram2d(
                bkg_arr["score"][bm],
                bkg_arr["mjj"][bm],
                bins=[xedges, yedges],
                weights=bkg_arr["weight"][bm],
            )
            frac = np.divide(s2d, s2d + b2d, out=np.zeros_like(s2d), where=(s2d + b2d) > 0)

            fig, ax = plt.subplots(figsize=(9, 7))
            mesh = ax.pcolormesh(xedges, yedges, frac.T, cmap="viridis", vmin=0.0, vmax=1.0)
            cb = fig.colorbar(mesh, ax=ax)
            cb.set_label("Signal fraction S/(S+B)")
            ax.axvline(0.7, color="white", linestyle="--", linewidth=1.5)
            ax.axhline(500.0, color="white", linestyle="--", linewidth=1.5)
            ax.set_xlabel("BDT score")
            ax.set_ylabel("mjj [GeV]")
            ax.set_title("D) mjj vs BDT score (SR)")
            pdf.savefig(fig)
            plt.close(fig)
        else:
            logger.info("Skipping mjj-vs-score 2D page: mjj branch not available.")

        # BDT output by VBF topology (tight/medium/loose)
        if np.isfinite(all_arr["mjj"]).any() and np.isfinite(all_arr["deta_jj"]).any():
            topo = {
                "tight_VBF": (700.0, 4.0),
                "medium_VBF": (400.0, 3.0),
                "loose_VBF": (250.0, 2.5),
            }
            for name, (mjj_cut, deta_cut) in topo.items():
                fig = plt.figure(figsize=(11, 8))
                gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
                ax_t = fig.add_subplot(gs[0, 0])
                ax_b = fig.add_subplot(gs[1, 0], sharex=ax_t)

                mc_topo = []
                for proc in sorted(mc_groups.keys(), key=proc_sort_key):
                    arr = mc_groups[proc]
                    m = make_mass_masks(arr, args)["sr"]
                    m = m & (arr["mjj"] > mjj_cut) & (np.abs(arr["deta_jj"]) > deta_cut)
                    if np.count_nonzero(m) == 0:
                        continue
                    mc_topo.append((proc, arr["score"][m], arr["weight"][m]))

                d_topo = None
                if data_sample is not None:
                    dm = make_mass_masks(data_sample.arrays, args)["sr"]
                    dm = dm & (data_sample.arrays["mjj"] > mjj_cut) & (np.abs(data_sample.arrays["deta_jj"]) > deta_cut)
                    d_topo = data_sample.arrays["score"][dm]

                stack_and_ratio(
                    ax_t,
                    ax_b,
                    score_bins,
                    f"4) {name}: mjj>{mjj_cut:.0f}, |deta_jj|>{deta_cut:.1f}",
                    mc_topo,
                    d_topo,
                    logy=True,
                )

                sig_topo = combine_group_arrays({k: v for k, v in mc_groups.items() if k in set(signal_set)}, signal_set)
                bkg_topo = combine_group_arrays({k: v for k, v in mc_groups.items() if k not in set(signal_set)}, [k for k in mc_groups.keys() if k not in set(signal_set)])
                ms = make_mass_masks(sig_topo, args)["sr"] if sig_topo["score"].size else np.array([], dtype=bool)
                mb = make_mass_masks(bkg_topo, args)["sr"] if bkg_topo["score"].size else np.array([], dtype=bool)
                if ms.size and mb.size:
                    ms = ms & (sig_topo["mjj"] > mjj_cut) & (np.abs(sig_topo["deta_jj"]) > deta_cut)
                    mb = mb & (bkg_topo["mjj"] > mjj_cut) & (np.abs(bkg_topo["deta_jj"]) > deta_cut)
                    ssum = float(sig_topo["weight"][ms].sum())
                    bsum = float(bkg_topo["weight"][mb].sum())
                    purity = ssum / max(ssum + bsum, 1e-12)
                    ax_t.text(
                        0.02,
                        0.94,
                        f"Signal purity = {purity:.3f}",
                        transform=ax_t.transAxes,
                        fontsize=9,
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                    )
                pdf.savefig(fig)
                plt.close(fig)
        else:
            logger.info("Skipping topology page(s): mjj and/or deta_jj unavailable.")

        # BDT category boundaries optimization
        if s_sc.size > 0 and b_sc.size > 0:
            cat_edges, cat_s, cat_b = optimize_categories_high_to_low(
                s_sc,
                s_w,
                b_sc,
                b_w,
                fine_bins=args.category_fine_bins,
                delta_z=args.category_delta_z,
                min_raw=args.category_min_raw,
                min_sig=args.category_min_sig,
            )
            z_tot = multi_bin_asimov_z(cat_s, cat_b)
            summary_lines.append(f"Optimized categories: n={len(cat_edges)-1}, Z_Asimov={z_tot:.4f}")

            fig, ax = plt.subplots(figsize=(11, 7))
            idx = np.arange(len(cat_s))
            ax.bar(idx, cat_b, width=0.62, label="Background", color=PROC_COLOR["DY"], alpha=0.75)

            ax_sig = ax.twinx()
            ax_sig.plot(
                idx,
                cat_s,
                color=PROC_COLOR["VBF"],
                marker="o",
                linewidth=2.0,
                label="Signal",
            )
            labels = [f"[{cat_edges[i]:.2f},{cat_edges[i+1]:.2f})" for i in range(len(cat_s))]
            ax.set_xticks(idx)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            ax.set_ylabel("Background weighted yield")
            ax_sig.set_ylabel("Signal weighted yield")

            b_max = float(np.max(cat_b)) if cat_b.size else 0.0
            s_max = float(np.max(cat_s)) if cat_s.size else 0.0
            ax.set_ylim(0.0, max(1.0, b_max * 1.25))
            ax_sig.set_ylim(0.0, max(1.0e-4, s_max * 1.35))

            ax.set_title("5) Optimized BDT categories (high-to-low merge)")
            ax.grid(alpha=0.3)
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax_sig.get_legend_handles_labels()
            if h1 or h2:
                ax.legend(h1 + h2, l1 + l2, loc="upper right")
            ax.text(0.01, 0.96, f"Total multibin Z_Asimov={z_tot:.3f}", transform=ax.transAxes)
            pdf.savefig(fig)
            plt.close(fig)
        else:
            logger.info("Skipping category-optimization page: signal/background score arrays are empty.")

        # Fold-based KS checks (detailed diagnostic page)
        fold_ok = np.isfinite(sig_arr["fold"]).any() and np.isfinite(bkg_arr["fold"]).any()
        if fold_ok:
            s_fold = sig_arr["fold"].astype(np.int32)
            b_fold = bkg_arr["fold"].astype(np.int32)
            s_mask_sr = make_mass_masks(sig_arr, args)["sr"]
            b_mask_sr = make_mass_masks(bkg_arr, args)["sr"]

            tr_s = s_mask_sr & np.isin(s_fold, [0, 1])
            va_s = s_mask_sr & np.isin(s_fold, [2])
            te_s = s_mask_sr & np.isin(s_fold, [3])
            tr_b = b_mask_sr & np.isin(b_fold, [0, 1])
            va_b = b_mask_sr & np.isin(b_fold, [2])
            te_b = b_mask_sr & np.isin(b_fold, [3])

            d_sig, p_sig = ks_2samp(sig_arr["score"][tr_s], sig_arr["score"][va_s])
            d_bkg, p_bkg = ks_2samp(bkg_arr["score"][tr_b], bkg_arr["score"][va_b])
            summary_lines.append(f"KS(signal train-vs-valid): D={d_sig:.4f}, p={p_sig:.4g}")
            summary_lines.append(f"KS(bkg train-vs-valid): D={d_bkg:.4f}, p={p_bkg:.4g}")

            m_sig_tr, _ = weighted_mean_std(sig_arr["score"][tr_s], sig_arr["weight"][tr_s])
            m_sig_va, _ = weighted_mean_std(sig_arr["score"][va_s], sig_arr["weight"][va_s])
            m_sig_te, _ = weighted_mean_std(sig_arr["score"][te_s], sig_arr["weight"][te_s])
            m_bkg_tr, _ = weighted_mean_std(bkg_arr["score"][tr_b], bkg_arr["weight"][tr_b])
            m_bkg_va, _ = weighted_mean_std(bkg_arr["score"][va_b], bkg_arr["weight"][va_b])
            m_bkg_te, _ = weighted_mean_std(bkg_arr["score"][te_b], bkg_arr["weight"][te_b])

            fig = plt.figure(figsize=(14, 9))
            gs = fig.add_gridspec(2, 2, height_ratios=[2.1, 1.4], hspace=0.25)
            ax_sig = fig.add_subplot(gs[0, 0])
            ax_bkg = fig.add_subplot(gs[0, 1], sharex=ax_sig)
            ax_tbl = fig.add_subplot(gs[1, :])

            centers = 0.5 * (score_bins[:-1] + score_bins[1:])

            def _norm_hist(vals: np.ndarray, ww: np.ndarray) -> np.ndarray:
                h = weighted_hist(vals, score_bins, ww)
                return h / h.sum() if h.sum() > 0 else h

            h_sig_tr = _norm_hist(sig_arr["score"][tr_s], sig_arr["weight"][tr_s])
            h_sig_va = _norm_hist(sig_arr["score"][va_s], sig_arr["weight"][va_s])
            h_sig_te = _norm_hist(sig_arr["score"][te_s], sig_arr["weight"][te_s])

            h_bkg_tr = _norm_hist(bkg_arr["score"][tr_b], bkg_arr["weight"][tr_b])
            h_bkg_va = _norm_hist(bkg_arr["score"][va_b], bkg_arr["weight"][va_b])
            h_bkg_te = _norm_hist(bkg_arr["score"][te_b], bkg_arr["weight"][te_b])

            ax_sig.step(centers, h_sig_tr, where="mid", linewidth=1.9, color="#1f77b4", label="Signal train")
            ax_sig.step(centers, h_sig_va, where="mid", linewidth=1.9, color="#ff7f0e", label="Signal valid")
            ax_sig.step(centers, h_sig_te, where="mid", linewidth=1.6, linestyle="--", color="#2ca02c", label="Signal test")
            ax_sig.set_title("Signal score by fold split")
            ax_sig.set_ylabel("Normalized entries")
            ax_sig.grid(alpha=0.3)
            safe_legend(ax_sig, fontsize=8)

            ax_bkg.step(centers, h_bkg_tr, where="mid", linewidth=1.9, color="#1f77b4", label="Bkg train")
            ax_bkg.step(centers, h_bkg_va, where="mid", linewidth=1.9, color="#ff7f0e", label="Bkg valid")
            ax_bkg.step(centers, h_bkg_te, where="mid", linewidth=1.6, linestyle="--", color="#2ca02c", label="Bkg test")
            ax_bkg.set_title("Background score by fold split")
            ax_bkg.grid(alpha=0.3)
            safe_legend(ax_bkg, fontsize=8)

            for aa in (ax_sig, ax_bkg):
                aa.set_xlabel("BDT score")

            ax_tbl.axis("off")
            status_col_labels = [name for name, _pmin, _dmax in ks_profiles]
            table_cols = [
                "Sample",
                "N train",
                "N valid",
                "N test",
                "Mean(train)",
                "Mean(valid)",
                "Mean(test)",
                "KS D",
                "KS p",
            ] + status_col_labels

            sig_statuses = [
                ks_profile_status(d_sig, p_sig, pmin, dmax)
                for _name, pmin, dmax in ks_profiles
            ]
            bkg_statuses = [
                ks_profile_status(d_bkg, p_bkg, pmin, dmax)
                for _name, pmin, dmax in ks_profiles
            ]

            table_rows = [
                [
                    "Signal",
                    f"{np.count_nonzero(tr_s)}",
                    f"{np.count_nonzero(va_s)}",
                    f"{np.count_nonzero(te_s)}",
                    f"{m_sig_tr:.4f}",
                    f"{m_sig_va:.4f}",
                    f"{m_sig_te:.4f}",
                    f"{d_sig:.4f}",
                    f"{p_sig:.4g}",
                ] + sig_statuses,
                [
                    "Background",
                    f"{np.count_nonzero(tr_b)}",
                    f"{np.count_nonzero(va_b)}",
                    f"{np.count_nonzero(te_b)}",
                    f"{m_bkg_tr:.4f}",
                    f"{m_bkg_va:.4f}",
                    f"{m_bkg_te:.4f}",
                    f"{d_bkg:.4f}",
                    f"{p_bkg:.4g}",
                ] + bkg_statuses,
            ]
            tbl = ax_tbl.table(cellText=table_rows, colLabels=table_cols, loc="center", cellLoc="center")
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8.2)
            tbl.scale(1.0, 1.35)
            ax_tbl.set_title("Fold-based KS checks (train vs validation)", fontsize=12, pad=10)

            status_start_col = len(table_cols) - len(ks_profiles)
            for r in range(1, 1 + len(table_rows)):
                for c in range(status_start_col, len(table_cols)):
                    cell = tbl[(r, c)]
                    val = str(table_rows[r - 1][c])
                    if val == "PASS":
                        cell.set_facecolor("#d9f2d9")
                    elif val == "WARN":
                        cell.set_facecolor("#fff4cc")
                    else:
                        cell.set_facecolor("#f8d7da")
                    cell.set_text_props(weight="bold")

            for c in range(len(table_cols)):
                hcell = tbl[(0, c)]
                hcell.set_facecolor("#e9ecef")
                hcell.set_text_props(weight="bold")

            ax_tbl.text(
                0.01,
                0.96,
                "KS thresholds: " + profile_desc,
                transform=ax_tbl.transAxes,
                fontsize=8.5,
                ha="left",
                va="top",
            )

            summary_lines.append("KS profile checks (Signal): " + ", ".join([f"{nm}={st}" for (nm, _, _), st in zip(ks_profiles, sig_statuses)]))
            summary_lines.append("KS profile checks (Background): " + ", ".join([f"{nm}={st}" for (nm, _, _), st in zip(ks_profiles, bkg_statuses)]))

            if np.isfinite(sig_arr["njets"]).any() and np.isfinite(bkg_arr["njets"]).any():
                corr_sig = weighted_corr(sig_arr["njets"][s_mask_sr], sig_arr["score"][s_mask_sr], sig_arr["weight"][s_mask_sr])
                corr_bkg = weighted_corr(bkg_arr["njets"][b_mask_sr], bkg_arr["score"][b_mask_sr], bkg_arr["weight"][b_mask_sr])
                summary_lines.append(f"Corr(score, njets) signal={corr_sig:.3f}, background={corr_bkg:.3f}")
                ax_tbl.text(
                    0.01,
                    0.02,
                    f"Corr(score, njets): signal={corr_sig:.3f}, background={corr_bkg:.3f}",
                    transform=ax_tbl.transAxes,
                    fontsize=9,
                    ha="left",
                    va="bottom",
                )

            pdf.savefig(fig)
            plt.close(fig)
        else:
            logger.info("Skipping fold-based KS checks: fold branch not available.")

        # Control-region response checks (always show zpeak/highSB/SR, annotate empty)
        region_defs = [
            ("zpeak", "Z-peak (70-110)"),
            ("highsb", "High sideband (>150)"),
            ("sr", "SR (110-150)"),
        ]
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
        centers = 0.5 * (score_bins[:-1] + score_bins[1:])

        for ax, (rname, rtitle) in zip(axes, region_defs):
            m_mc = make_mass_masks(all_arr, args)[rname]
            h_mc = weighted_hist(all_arr["score"][m_mc], score_bins, all_arr["weight"][m_mc])
            mc_yield = float(np.sum(all_arr["weight"][m_mc]))

            h_data: Optional[np.ndarray] = None
            data_yield = 0.0
            if data_sample is not None:
                m_d = make_mass_masks(data_sample.arrays, args)[rname]
                h_data = weighted_hist(data_sample.arrays["score"][m_d], score_bins, None)
                data_yield = float(np.sum(h_data))

            plotted = False
            if np.sum(h_mc) > 0:
                h_mc = h_mc / h_mc.sum()
                ax.step(centers, h_mc, where="mid", color="tab:blue", linewidth=2, label="MC total")
                plotted = True

            if h_data is not None and np.sum(h_data) > 0:
                h_data = h_data / h_data.sum()
                ax.step(centers, h_data, where="mid", color="black", linewidth=1.5, label="Data")
                plotted = True

            if not plotted:
                ax.text(0.5, 0.5, "No events", transform=ax.transAxes, ha="center", va="center", fontsize=10)

            ax.text(
                0.03,
                0.96,
                f"MC={mc_yield:.1f}\nData={data_yield:.0f}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
            )
            ax.set_title(rtitle)
            ax.set_xlabel("BDT score")
            ax.grid(alpha=0.3)
            safe_legend(ax, fontsize=8)

        axes[0].set_ylabel("Normalized entries")
        fig.suptitle("Control-region BDT response")
        pdf.savefig(fig)
        plt.close(fig)

        # BDT score vs VBF kinematics profiles
        prof_done = False
        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        if np.isfinite(sig_arr["mjj"]).any() and np.isfinite(bkg_arr["mjj"]).any():
            sb = make_mass_masks(sig_arr, args)["sr"]
            bb = make_mass_masks(bkg_arr, args)["sr"]
            bins_mjj = np.array([0, 200, 400, 600, 800, 1200, 1800, 2500], dtype=np.float64)
            xc, ym, ye = profile_mean(sig_arr["mjj"][sb], sig_arr["score"][sb], sig_arr["weight"][sb], bins_mjj)
            axes[0].errorbar(xc, ym, yerr=ye, fmt="o-", color=PROC_COLOR["VBF"], label="Signal")
            xc, ym, ye = profile_mean(bkg_arr["mjj"][bb], bkg_arr["score"][bb], bkg_arr["weight"][bb], bins_mjj)
            axes[0].errorbar(xc, ym, yerr=ye, fmt="s-", color=PROC_COLOR["DY"], label="Background")
            axes[0].set_xlabel("mjj [GeV]")
            axes[0].set_ylabel("<BDT score>")
            axes[0].grid(alpha=0.3)
            safe_legend(axes[0], fontsize=8)
            prof_done = True
        else:
            axes[0].axis("off")

        if np.isfinite(sig_arr["deta_jj"]).any() and np.isfinite(bkg_arr["deta_jj"]).any():
            sb = make_mass_masks(sig_arr, args)["sr"]
            bb = make_mass_masks(bkg_arr, args)["sr"]
            bins_de = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=np.float64)
            xc, ym, ye = profile_mean(np.abs(sig_arr["deta_jj"][sb]), sig_arr["score"][sb], sig_arr["weight"][sb], bins_de)
            axes[1].errorbar(xc, ym, yerr=ye, fmt="o-", color=PROC_COLOR["VBF"], label="Signal")
            xc, ym, ye = profile_mean(np.abs(bkg_arr["deta_jj"][bb]), bkg_arr["score"][bb], bkg_arr["weight"][bb], bins_de)
            axes[1].errorbar(xc, ym, yerr=ye, fmt="s-", color=PROC_COLOR["DY"], label="Background")
            axes[1].set_xlabel("|deta_jj|")
            axes[1].set_ylabel("<BDT score>")
            axes[1].grid(alpha=0.3)
            safe_legend(axes[1], fontsize=8)
            prof_done = True
        else:
            axes[1].axis("off")

        if np.isfinite(sig_arr["pt_j1"]).any() and np.isfinite(bkg_arr["pt_j1"]).any():
            sb = make_mass_masks(sig_arr, args)["sr"]
            bb = make_mass_masks(bkg_arr, args)["sr"]
            bins_pt = np.array([0, 40, 60, 80, 120, 160, 220, 320, 500], dtype=np.float64)
            xc, ym, ye = profile_mean(sig_arr["pt_j1"][sb], sig_arr["score"][sb], sig_arr["weight"][sb], bins_pt)
            axes[2].errorbar(xc, ym, yerr=ye, fmt="o-", color=PROC_COLOR["VBF"], label="Signal")
            xc, ym, ye = profile_mean(bkg_arr["pt_j1"][bb], bkg_arr["score"][bb], bkg_arr["weight"][bb], bins_pt)
            axes[2].errorbar(xc, ym, yerr=ye, fmt="s-", color=PROC_COLOR["DY"], label="Background")
            axes[2].set_xlabel("jet1 pt [GeV]")
            axes[2].set_ylabel("<BDT score>")
            axes[2].grid(alpha=0.3)
            safe_legend(axes[2], fontsize=8)
            prof_done = True
        else:
            axes[2].axis("off")

        if prof_done:
            fig.suptitle("BDT score vs VBF kinematics")
            pdf.savefig(fig)
            plt.close(fig)
        else:
            plt.close(fig)
            logger.info("Skipping BDT-vs-kinematics page: missing mjj/deta_jj/pt_j1 branches.")

        # VBF acceptance and purity by BDT score
        vbf_arr = combine_group_arrays(mc_groups, vbf_set)
        non_vbf_arr = combine_group_arrays(mc_groups, [k for k in mc_groups.keys() if k not in set(vbf_set)])
        if vbf_arr["score"].size > 0:
            bins_acc = np.array(parse_csv_floats(args.acceptance_bins), dtype=np.float64)
            bins_acc = np.unique(np.clip(bins_acc, 0.0, 1.0))
            if bins_acc[0] > 0.0:
                bins_acc = np.insert(bins_acc, 0, 0.0)
            if bins_acc[-1] < 1.0:
                bins_acc = np.append(bins_acc, 1.0)

            m_vbf = make_mass_masks(vbf_arr, args)["sr"]
            m_non = make_mass_masks(non_vbf_arr, args)["sr"] if non_vbf_arr["score"].size else np.array([], dtype=bool)

            s_hist = weighted_hist(vbf_arr["score"][m_vbf], bins_acc, vbf_arr["weight"][m_vbf])
            b_hist = weighted_hist(non_vbf_arr["score"][m_non], bins_acc, non_vbf_arr["weight"][m_non]) if m_non.size else np.zeros(len(bins_acc)-1)

            s_tot = float(s_hist.sum())

            simple_mask_vbf = np.zeros_like(m_vbf, dtype=bool)
            if np.isfinite(vbf_arr["mjj"]).any() and np.isfinite(vbf_arr["deta_jj"]).any() and np.isfinite(vbf_arr["njets"]).any():
                simple_mask_vbf = (
                    m_vbf
                    & (vbf_arr["mjj"] > 500.0)
                    & (np.abs(vbf_arr["deta_jj"]) > 3.5)
                    & (vbf_arr["njets"] >= 2)
                )

            rows = []
            for i, (lo, hi) in enumerate(zip(bins_acc[:-1], bins_acc[1:])):
                m_bin = m_vbf & (vbf_arr["score"] >= lo) & (vbf_arr["score"] < hi)
                s_bin = float(vbf_arr["weight"][m_bin].sum())
                s_bin_simple = float(vbf_arr["weight"][m_bin & simple_mask_vbf].sum()) if simple_mask_vbf.size else float("nan")
                b_bin = float(b_hist[i])
                frac_vbf = s_bin / max(s_tot, 1e-12)
                frac_simple = s_bin_simple / max(s_bin, 1e-12) if np.isfinite(s_bin_simple) else float("nan")
                purity = s_bin / max(s_bin + b_bin, 1e-12)
                rows.append([lo, hi, frac_vbf, frac_simple, purity])

            fig, ax = plt.subplots(figsize=(11, 6))
            labels = [f"[{r[0]:.2f},{r[1]:.2f})" for r in rows]
            x = np.arange(len(rows))
            frac_vbf = np.asarray([r[2] for r in rows], dtype=np.float64)
            frac_simple = np.asarray([r[3] for r in rows], dtype=np.float64)
            purity = np.asarray([r[4] for r in rows], dtype=np.float64)

            ax.plot(x, frac_vbf, marker="o", linewidth=2.0, label="Fraction of VBF signal", color="#1f77b4")

            ax_zoom = ax.twinx()
            ax_zoom.plot(x, frac_simple, marker="s", linewidth=2.0, label="VBF passing simple cuts", color="#ff7f0e")
            ax_zoom.plot(x, purity, marker="^", linewidth=2.0, label="Purity S/(S+B)", color="#2ca02c")
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
            ax.set_ylim(0.0, max(1.0, float(np.nanmax(frac_vbf)) * 1.15))
            small_max = float(np.nanmax(np.nan_to_num(np.concatenate([frac_simple, purity]), nan=0.0)))
            ax_zoom.set_ylim(0.0, max(0.02, small_max * 1.3))

            ax.set_ylabel("Fraction of VBF signal")
            ax_zoom.set_ylabel("Small-fraction metrics")
            ax.set_title("10) VBF signal acceptance and purity by BDT score")
            ax.grid(alpha=0.3)

            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax_zoom.get_legend_handles_labels()
            if h1 or h2:
                ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)

            high_bin_purity = rows[-1][4] if rows else float("nan")
            summary_lines.append(f"Highest BDT-bin purity={high_bin_purity:.3f}")
        else:
            logger.info("Skipping VBF-acceptance page: no VBF-labeled sample after process mapping.")

        add_note_page(pdf, "Validation Summary", summary_lines)

    logger.info("Wrote: %s", pdf_path)


if __name__ == "__main__":
    main()
