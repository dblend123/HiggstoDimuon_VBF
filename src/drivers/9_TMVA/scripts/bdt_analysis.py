#!/usr/bin/env python3

import argparse
import glob
import os
import sys
import logging
import itertools
import math

import numpy as np
import ROOT
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Add parent directory to path so we can import tmva100_DataPreparation
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tmva100_DataPreparation import VARIABLES
from significance_inference import (
    build_shape_inputs,
    compute_shape_expected_significance,
    run_shape_toy_ensembles,
    run_bias_test_matrix,
)

ROOT.gROOT.SetBatch(True)

# Reference luminosity (fb^-1) for significance scaling across eras
# Use this to normalize multi-bin yields to a common scale
REF_LUMI_FB = 100.0


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
    p.add_argument("--alpha", type=float, default=1e-3, help="KS p-value threshold")
    p.add_argument("--dthr-pct", type=float, default=0.20, help="KS D(%) threshold")
    p.add_argument("--abs-weights", dest="abs_weights", action="store_true", default=True,
                   help="Use |weight| for yields/significance (default: enabled)")
    p.add_argument("--signed-weights", dest="abs_weights", action="store_false",
                   help="Disable absolute-value transform (negative/invalid weights are still removed)")
    p.add_argument("--output", default="bdt_analysis.pdf")
    p.add_argument("--multibin-min-n", type=int, default=4,
                   help="Minimum number of BDT categories for multi-bin optimization")
    p.add_argument("--multibin-max-n", type=int, default=6,
                   help="Maximum number of BDT categories for multi-bin optimization")
    p.add_argument("--multibin-candidates", type=int, default=16,
                   help="Number of interior candidate boundary points for scanning")
    p.add_argument("--multibin-min-sig", type=float, default=1.0,
                   help="Minimum weighted signal yield required in each BDT bin")
    p.add_argument("--multibin-min-bkg", type=float, default=5.0,
                   help="Minimum weighted background yield required in each BDT bin")
    p.add_argument("--multibin-no-monotonic", action="store_true",
                   help="Disable strict monotonic S/B constraint across ordered BDT bins")
    p.add_argument("--multibin-objective", choices=["asimov", "likelihood"], default="likelihood",
                   help="Objective for selecting BDT bin edges: ZA(total) or profile-likelihood Zexp")
    p.add_argument("--multibin-likelihood-topk", type=int, default=20,
                   help="For likelihood objective, evaluate full profile-likelihood only for top-K ZA candidates")
    p.add_argument("--shape-likelihood", action="store_true",
                   help="Run simultaneous shape-based likelihood across optimized BDT bins")
    p.add_argument("--shape-mass-nbins", type=int, default=40,
                   help="Number of mass bins for simultaneous likelihood")
    p.add_argument("--shape-poly-order", type=int, default=1,
                   help="Per-bin polynomial modifier order for background (0-2)")
    p.add_argument("--shape-bkg-norm-unc", type=float, default=0.20,
                   help="Gaussian prior width (fractional) on per-bin background normalization")
    p.add_argument("--shape-poly-unc", type=float, default=0.30,
                   help="Gaussian prior width for per-bin background shape coefficients")
    p.add_argument("--shape-sig-norm-unc", type=float, default=0.05,
                   help="Gaussian prior width (fractional) on correlated signal normalization")
    p.add_argument("--shape-sig-mass-shift-unc", type=float, default=0.20,
                   help="1σ correlated signal mass shift in GeV for template morphing")
    p.add_argument("--shape-sig-resolution-unc", type=float, default=0.10,
                   help="1σ correlated relative signal resolution nuisance for template morphing")
    p.add_argument("--shape-mu-max", type=float, default=3.0,
                   help="Upper bound of mu grid in profiled likelihood scan")
    p.add_argument("--shape-mu-points", type=int, default=121,
                   help="Number of mu scan points for profile likelihood")
    p.add_argument("--shape-data-source", choices=["asimov", "mc-splusb", "mc-bonly"], default="asimov",
                   help="Data used in shape-likelihood fit: Asimov S+B, MC S+B histogram, or MC B-only histogram")
    p.add_argument("--shape-core-poly-order", type=int, default=2,
                   help="Polynomial order used for core-PDF per-BDT deformation (1-2)")
    p.add_argument("--shape-core-use-lowmed", action="store_true",
                   help="Build core background template from low/medium BDT bins only (exclude highest bin)")
    p.add_argument("--shape-core-sb-threshold", type=float, default=0.02,
                   help="Maximum S/B allowed for a BDT bin to contribute to core template when low/med core is enabled")
    p.add_argument("--shape-eps-prior-scale", type=float, default=1.0,
                   help="Scale factor for epsilon deformation Gaussian priors around data-driven anchors")
    p.add_argument("--dy-ge2-unc", type=float, default=0.05,
                   help="Log-normal 1σ fractional uncertainty for DY component with nGenJets>=2")
    p.add_argument("--genjet-branch", default="auto",
                   help="Generator-jet multiplicity branch for DY split (auto tries common names)")
    p.add_argument("--shape-toys", type=int, default=0,
                   help="Number of toy ensembles to run (0 disables)")
    p.add_argument("--shape-toy-seed", type=int, default=12345,
                   help="Random seed for toy ensembles")
    p.add_argument("--shape-toy-mu", default="0,1",
                   help="Comma-separated injected signal strengths for toy tests (e.g. '0,1')")
    p.add_argument("--shape-bias-tests", action="store_true",
                   help="Run generation-vs-fit background-function bias matrix and report pass/fail")
    p.add_argument("--shape-bias-mu-true", type=float, default=1.0,
                   help="Injected signal strength for bias tests (mu_true)")
    p.add_argument("--shape-bias-threshold", type=float, default=0.20,
                   help="Bias pass threshold on |(mu_fit-mu_true)/mu_true|")
    p.add_argument("--shape-bias-toys", type=int, default=30,
                   help="Number of pseudo-datasets per generation/fit model pair for bias tests")
    p.add_argument("--shape-bias-param-grid", type=int, default=21,
                   help="Grid points used for internal parameter scan in background-model fits")

    # Optional Data/MC summary (same format as tmva_DataMC_BDT.py)
    p.add_argument("--data-era", action="append", default=[],
                   help="ERA:/path/data1.root[,/path/data2.root]")
    p.add_argument("--mc-era", action="append", default=[],
                   help="ERA:PROC:/path/mc1.root[,/path/mc2.root]  (repeat)")
    p.add_argument("--datamc-nbins", type=int, default=40)
    p.add_argument("--datamc-logy", action="store_true")
    p.add_argument("--lumi", type=float, default=None,
                   help="Luminosity in fb^-1; if set, MC is rescaled to match data counts. For 2024, the value is 104.80585 fb^-1.")

    # Optional model-based scan (tmva_ScanYields.py logic)
    p.add_argument("--model-file", default="tmva101.root")
    p.add_argument("--model-name", default="myBDT")
    p.add_argument("--scan-sig", nargs="+", default=[], help="Signal inputs (files/dirs/globs)")
    p.add_argument("--scan-bkg", nargs="+", default=[], help="Background inputs (files/dirs/globs)")
    p.add_argument("--scan-cuts", default="0.00:1.00:401",
                   help="scan range as start:end:npoints, inclusive")

    # Optional retrain-per-fold KS (tmva_KS_perfold.py logic)
    p.add_argument("--retrain-ks", action="store_true")
    p.add_argument("--ks-dir", default="DataPreparation_Files",
                   help="Dir containing *_all_signal.root and *_all_background.root")

    # Luminosity normalization for significance calculations
    p.add_argument("--sig-lumi", type=float, default=None,
                   help="Signal luminosity in fb^-1 for normalizing yields in multi-bin optimization. If not set, raw weights are used.")
    p.add_argument("--sig-xsec", type=float, default=None,
                   help="Signal cross-section in pb (used with --sig-lumi to normalize). If set with --sig-lumi, signal weights are rescaled.")

    # Combined era fit
    p.add_argument("--combined-eras", nargs="+", default=None,
                   help="Combine multiple eras into single fit (e.g., 12022 22022). Each era is lumi-scaled then concatenated.")
    p.add_argument("--no-lumi-scaling", action="store_true",
                   help="Disable era luminosity rescaling to reference luminosity")

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


def resolve_mass_branch(path, treename, preferred):
    cols = set(list_branches(path, treename))
    if preferred != "auto" and preferred in cols:
        return preferred
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


def find_scored_process_file(outputs_dir, era, proc):
    proc = str(proc)
    candidates = [
        os.path.join(outputs_dir, f"{era}_scored_{proc}.root"),
        os.path.join(outputs_dir, era, f"{era}_scored_{proc}.root"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    return None


def resolve_genjet_branch(path, treename, preferred="auto"):
    cols = set(list_branches(path, treename))
    if preferred and preferred != "auto" and preferred in cols:
        return preferred
    candidates = [
        "nGenJet", "nGenJets", "n_genjet", "nGenJ", "nGenJetAK4",
        "nJetGen", "n_jets_gen",
    ]
    for c in candidates:
        if c in cols:
            return c
    if "njets" in cols:
        return "njets"
    return None


def sanitize_weights(w, use_abs=False):
    w = np.asarray(w, dtype=np.float64)
    w[~np.isfinite(w)] = 0.0
    if use_abs:
        w = np.abs(w)
    w[w < 0.0] = 0.0
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


def asimov_q(S, B):
    if S <= 0 or B <= 0:
        return 0.0
    return float(2.0 * ((S + B) * np.log(1.0 + S / max(B, 1e-12)) - S))


def multi_bin_asimov_z(s_bins, b_bins):
    qsum = 0.0
    for s, b in zip(s_bins, b_bins):
        qsum += asimov_q(float(s), float(b))
    return float(np.sqrt(max(qsum, 0.0)))


def format_ratio(val):
    if not np.isfinite(val):
        return "inf"
    return f"{val:.3g}"


def parse_group_map(items, flagname):
    out = {}
    for s in items:
        if ":" not in s:
            raise SystemExit(f"[error] {flagname} expects key:val, got '{s}'")
        k, v = s.split(":", 1)
        out.setdefault(k, []).append(v)
    return out


def canonical_proc(name):
    s = name.lower()
    if "dy" in s or "z" in s:
        return "DY"
    if "ewk" in s or "z+ewk" in s:
        return "EWK"
    if "top" in s or "tt" in s:
        return "Top"
    if "vv" in s or "dibos" in s:
        return "Diboson"
    if "ggh" in s:
        return "ggH"
    if "vbf" in s:
        return "VBF"
    return "Other"


def proc_sort_key(name):
    order = ["DY", "EWK", "Top", "Diboson", "Other", "ggH", "VBF"]
    cn = canonical_proc(name)
    return order.index(cn) if cn in order else len(order)


def proc_color(name):
    colors = {
        "DY": "#f4a261",
        "EWK": "#6a4c93",
        "Top": "#2a9d8f",
        "Diboson": "#457b9d",
        "Other": "#8d99ae",
        "ggH": "#7f7f7f",
        "VBF": "#4caf50",
    }
    return colors.get(canonical_proc(name), "#7f7f7f")


def expand_inputs(items):
    files = []
    for it in items:
        if os.path.isdir(it):
            files += sorted(glob.glob(os.path.join(it, "*.root")))
        elif it.endswith(".root") and os.path.isfile(it):
            files.append(it)
        else:
            files += sorted(glob.glob(it))
    return files


def parse_cuts(spec):
    a, b, n = spec.split(":")
    return np.linspace(float(a), float(b), int(n))


def parse_float_csv(spec, default_vals=None):
    if spec is None:
        return list(default_vals or [])
    txt = str(spec).strip()
    if not txt:
        return list(default_vals or [])
    out = []
    for token in txt.split(","):
        token = token.strip()
        if not token:
            continue
        try:
            out.append(float(token))
        except ValueError:
            continue
    return out if out else list(default_vals or [])


def resolve_shape_data_override(shape_inputs, shape_data_source):
    src = str(shape_data_source).strip().lower()
    if src == "mc-splusb":
        return np.asarray(shape_inputs["s_hist"] + shape_inputs["b_hist"], dtype=np.float64)
    if src == "mc-bonly":
        return np.asarray(shape_inputs["b_hist"], dtype=np.float64)
    return None


def compute_scores_model(model, xs):
    X = np.vstack(xs).T
    return model.Compute(X.astype(np.float32, copy=False)).astype(np.float64, copy=False)


def compute_bin_yields(scores, weights, edges):
    vals = np.asarray(scores, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    out = []
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        if i == len(edges) - 2:
            mask = (vals >= lo) & (vals <= hi)
        else:
            mask = (vals >= lo) & (vals < hi)
        out.append(float(w[mask].sum()))
    return np.asarray(out, dtype=np.float64)


def candidate_boundaries(sig_scores, sig_w, bkg_scores, bkg_w, n_candidates):
    n_candidates = max(int(n_candidates), 4)
    all_scores = np.concatenate([sig_scores, bkg_scores])
    all_w = np.concatenate([
        sanitize_weights(sig_w, use_abs=True),
        sanitize_weights(bkg_w, use_abs=True),
    ])
    qs = np.linspace(0.0, 1.0, n_candidates + 2)[1:-1]
    cands = []
    for q in qs:
        v = weighted_quantile(all_scores, all_w, float(q))
        if np.isfinite(v) and (0.0 < v < 1.0):
            cands.append(float(v))
    cands = sorted(set([round(x, 6) for x in cands]))
    return cands


def evaluate_multibin_config(sig_scores, sig_w, bkg_scores, bkg_w, edges,
                             min_sig, min_bkg, enforce_monotonic=True):
    s_bins = compute_bin_yields(sig_scores, sig_w, edges)
    b_bins = compute_bin_yields(bkg_scores, bkg_w, edges)

    if np.any(s_bins < float(min_sig)) or np.any(b_bins < float(min_bkg)):
        return None

    ratios = np.divide(s_bins, b_bins, out=np.full_like(s_bins, np.inf, dtype=np.float64), where=b_bins > 0)
    if enforce_monotonic and np.any(np.diff(ratios) <= 0.0):
        return None

    zapprox_bins = np.divide(s_bins, np.sqrt(np.maximum(b_bins, 1e-12)))
    ztot = multi_bin_asimov_z(s_bins, b_bins)

    return {
        "edges": [float(x) for x in edges],
        "s_bins": s_bins,
        "b_bins": b_bins,
        "sb_bins": ratios,
        "zapprox_bins": zapprox_bins,
        "za_total": float(ztot),
        "n_bins": int(len(edges) - 1),
    }


def optimize_multibin_boundaries(sig_scores, sig_w, bkg_scores, bkg_w,
                                 min_n=4, max_n=6, n_candidates=16,
                                 min_sig=1.0, min_bkg=5.0, enforce_monotonic=True,
                                 objective="asimov",
                                 sig_mass=None, bkg_mass=None,
                                 mass_min=110.0, mass_max=150.0, mass_nbins=40,
                                 poly_order=1,
                                 bkg_norm_unc=0.20,
                                 poly_unc=0.30,
                                 sig_norm_unc=0.05,
                                 sig_mass_shift_unc=0.20,
                                 sig_resolution_unc=0.10,
                                 mu_max=3.0,
                                 mu_points=121,
                                 likelihood_topk=20,
                                 shape_data_source="asimov"):
    min_n = max(2, int(min_n))
    max_n = max(min_n, int(max_n))
    cands = candidate_boundaries(sig_scores, sig_w, bkg_scores, bkg_w, n_candidates)

    objective = str(objective).strip().lower()
    if objective not in ("asimov", "likelihood"):
        objective = "asimov"

    use_likelihood = (
        objective == "likelihood"
        and sig_mass is not None
        and bkg_mass is not None
    )

    metric_key = "za_total"
    metric_label = "ZA(total)"
    if use_likelihood:
        metric_key = "zexp_profile"
        metric_label = "Zexp(profile)"

    best_global = None
    best_by_n = {}
    scan_stats = {}

    for n_bins in range(min_n, max_n + 1):
        n_inner = n_bins - 1
        if len(cands) < n_inner:
            scan_stats[n_bins] = {"scanned": 0, "valid": 0}
            continue

        scanned = 0
        valid = 0
        best_local = None
        valid_results = []
        for interior in itertools.combinations(cands, n_inner):
            scanned += 1
            edges = [0.0] + list(interior) + [1.0]
            res = evaluate_multibin_config(
                sig_scores, sig_w, bkg_scores, bkg_w, edges,
                min_sig=min_sig, min_bkg=min_bkg, enforce_monotonic=enforce_monotonic,
            )
            if res is None:
                continue
            valid += 1
            valid_results.append(res)

        if use_likelihood and valid_results:
            topk = max(1, int(likelihood_topk))
            valid_results = sorted(valid_results, key=lambda r: r["za_total"], reverse=True)
            shortlisted = valid_results[:min(topk, len(valid_results))]
            evaluated = []
            for res in shortlisted:
                shape_inputs = build_shape_inputs(
                    sig_mass,
                    sig_scores,
                    sig_w,
                    bkg_mass,
                    bkg_scores,
                    bkg_w,
                    res["edges"],
                    mass_min,
                    mass_max,
                    mass_nbins,
                )
                shape_fit = compute_shape_expected_significance(
                    shape_inputs,
                    poly_order=poly_order,
                    bkg_norm_unc=bkg_norm_unc,
                    poly_unc=poly_unc,
                    sig_norm_unc=sig_norm_unc,
                    sig_mass_shift_unc=sig_mass_shift_unc,
                    sig_resolution_unc=sig_resolution_unc,
                    mu_max=mu_max,
                    mu_points=mu_points,
                    data_override=resolve_shape_data_override(shape_inputs, shape_data_source),
                )
                out = dict(res)
                out["zexp_profile"] = float(shape_fit["zexp"])
                out["q0_profile"] = float(shape_fit["q0"])
                out["mu_hat_profile"] = float(shape_fit["mu_hat"])
                evaluated.append(out)
            if evaluated:
                best_local = max(evaluated, key=lambda r: r[metric_key])
        elif valid_results:
            best_local = max(valid_results, key=lambda r: r["za_total"])

        scan_stats[n_bins] = {"scanned": int(scanned), "valid": int(valid)}
        if best_local is not None:
            best_local["metric_key"] = metric_key
            best_local["metric_label"] = metric_label
            best_local["opt_metric"] = float(best_local[metric_key])
            best_by_n[n_bins] = best_local
            if best_global is None or best_local["opt_metric"] > best_global["opt_metric"]:
                best_global = best_local

    return {
        "best": best_global,
        "best_by_n": best_by_n,
        "scan_stats": scan_stats,
        "candidates": cands,
        "objective_requested": objective,
        "objective_used": "likelihood" if use_likelihood else "asimov",
        "metric_key": metric_key,
        "metric_label": metric_label,
    }


def plot_shape_fit_overlays(pdf, title, mass_edges, bdt_edges, asimov, model_b, model_sb):
    n_bdt = asimov.shape[0]
    ncols = 2 if n_bdt > 1 else 1
    nrows = int(np.ceil(n_bdt / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.7 * nrows), sharex=True)
    axes = np.atleast_1d(axes).flatten()
    centers = 0.5 * (mass_edges[:-1] + mass_edges[1:])

    for i in range(len(axes)):
        ax = axes[i]
        if i >= n_bdt:
            ax.axis("off")
            continue
        lo, hi = bdt_edges[i], bdt_edges[i + 1]
        ax.errorbar(centers, asimov[i], yerr=np.sqrt(np.maximum(asimov[i], 1.0)),
                    fmt="o", color="black", label="Asimov")
        ax.plot(centers, model_b[i], color="tab:orange", linewidth=2, label="B-only fit")
        ax.plot(centers, model_sb[i], color="tab:blue", linewidth=2, linestyle="--", label="S+B best fit")
        ax.set_title(f"Bin {i + 1}: {lo:.3f}-{hi:.3f}")
        ax.set_xlabel("mμμ [GeV]")
        ax.set_ylabel("Expected events")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8)

    fig.suptitle(title, fontsize=14)
    pdf.savefig(fig)
    plt.close(fig)


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


def plot_mass_sculpt_norm(pdf, title, mass, score, weight, bdt_bins, corr_val):
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(bdt_bins) - 1))
    for i, (lo, hi) in enumerate(zip(bdt_bins[:-1], bdt_bins[1:])):
        mask = (score >= lo) & (score < hi)
        if not np.any(mask):
            continue
        m = mass[mask]
        w = weight[mask]
        ax.hist(m, bins=40, weights=w, histtype="step",
                linewidth=2, color=colors[i], density=True,
                label=f"{lo:.2f} <= BDT < {hi:.2f}")
    ax.set_xlabel("Mass [GeV]")
    ax.set_ylabel("A.U.")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    if np.isfinite(corr_val):
        ax.text(
            0.02, 0.95,
            f"Weighted corr(mass, BDT) = {corr_val:.3f}",
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


def plot_datamc(pdf, title, data_scores, mc_groups, bins, logy):
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(11, 9), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )

    mc_hists = []
    mc_labels = []
    mc_colors = []
    mc_sum = np.zeros(len(bins) - 1, dtype=np.float64)

    for name, scores, weights in mc_groups:
        hist, _ = np.histogram(scores, bins=bins, weights=weights)
        mc_hists.append(hist)
        mc_labels.append(canonical_proc(name))
        mc_colors.append(proc_color(name))
        mc_sum += hist

    ax_top.hist(
        [bins[:-1]] * len(mc_hists),
        bins=bins,
        weights=mc_hists,
        stacked=True,
        color=mc_colors,
        label=mc_labels,
        histtype="stepfilled",
        alpha=0.85,
    )
    data_hist, _ = np.histogram(data_scores, bins=bins)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    ax_top.errorbar(bin_centers, data_hist, yerr=np.sqrt(np.maximum(data_hist, 1.0)),
                    fmt="o", color="black", label="Data")
    ax_top.set_ylabel("Events / bin")
    ax_top.set_title(title)
    ax_top.grid(alpha=0.3)
    if logy:
        ax_top.set_yscale("log")
    ax_top.legend(fontsize=8, ncol=2)

    ratio = np.divide(data_hist, mc_sum, out=np.zeros_like(data_hist, dtype=np.float64), where=mc_sum > 0)
    ratio_err = np.divide(np.sqrt(np.maximum(data_hist, 1.0)), mc_sum, out=np.zeros_like(data_hist, dtype=np.float64), where=mc_sum > 0)
    ax_bot.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    ax_bot.errorbar(bin_centers, ratio, yerr=ratio_err, fmt="o", color="black")
    ax_bot.set_xlabel("BDT score")
    ax_bot.set_ylabel("Data/MC")
    ax_bot.set_ylim(0.5, 1.5)
    ax_bot.grid(alpha=0.3)

    pdf.savefig(fig)
    plt.close(fig)


def plot_scan_curves(pdf, title, cuts, z_vals, best_z):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(cuts, z_vals, linewidth=2, label="Z = S/sqrt(S+B)")
    ax.axvline(best_z["cut"], color="tab:red", linestyle="--", linewidth=1.5)
    ax.set_xlabel("BDT score cut")
    ax.set_ylabel("Significance")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    ax.text(
        0.02, 0.95,
        "Best cut\n"
        f"Z: cut={best_z['cut']:.3f}, Z={best_z['val']:.3f}",
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    pdf.savefig(fig)
    plt.close(fig)


def rescale_mc_to_data(data_count, mc_groups, logger):
    """
    Auto-rescale MC weights to match data counts.
    Assumes MC has NOT been normalized by cross-section × luminosity.
    """
    if data_count <= 0 or not mc_groups:
        return mc_groups, 1.0
    
    mc_total_weight = sum(float(w.sum()) for _, _, w in mc_groups)
    if mc_total_weight <= 0:
        logger.warning("MC total weight is 0; cannot rescale")
        return mc_groups, 1.0
    
    scale = data_count / mc_total_weight
    logger.info("MC rescale: data_count=%.0f, mc_total_weight=%.1f → scale=%.6f",
                data_count, mc_total_weight, scale)
    
    scaled = []
    for proc, scores, weights in mc_groups:
        scaled.append((proc, scores, weights * scale))
    
    return scaled, scale


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


def get_era_luminosity_fb(era):
    """
    Return approximate integrated luminosity (fb^-1) for given era.
    Maps common CMS era codes to known luminosities.
    """
    lumi_map = {
        "12022": 7.9805,   
        "12023": 17.794,    
        "2024": 104.80585,    
        "22022": 26.672,   
        "22023": 9.4510,    
    }
    return lumi_map.get(era, None)

def apply_lumi_scaling(sig_w, bkg_w, era_lumi_fb, ref_lumi_fb=None, logger=None):
    """
    Rescale signal and background weights to a reference luminosity.
    If era_lumi_fb is None, no rescaling is applied.
    Returns rescaled weights and the scale factor used.
    """
    if era_lumi_fb is None or era_lumi_fb <= 0:
        return sig_w, bkg_w, 1.0
    
    if ref_lumi_fb is None:
        ref_lumi_fb = REF_LUMI_FB
    
    scale = float(ref_lumi_fb / era_lumi_fb)
    if logger:
        logger.info("Luminosity scaling: era_lumi=%.1f fb^-1 → ref_lumi=%.1f fb^-1 → scale=%.4f",
                   era_lumi_fb, ref_lumi_fb, scale)
    
    return sig_w * scale, bkg_w * scale, scale


def combined_era_fit(outputs_dir, combined_eras, treename, bdt_branch, weight_branch,
                     mass_branch, mass_min, mass_max, abs_weights, logger,
                     multibin_min_n, multibin_max_n, multibin_candidates,
                     multibin_min_sig, multibin_min_bkg, multibin_no_monotonic,
                     shape_likelihood, shape_mass_nbins, shape_poly_order,
                     shape_bkg_norm_unc, shape_poly_unc, shape_sig_norm_unc,
                     shape_sig_mass_shift_unc, shape_sig_resolution_unc,
                     shape_mu_max, shape_mu_points, shape_data_source,
                     shape_core_poly_order, shape_core_use_lowmed,
                     shape_core_sb_threshold, shape_eps_prior_scale, dy_ge2_unc, genjet_branch,
                     no_lumi_scaling,
                     shape_toys, shape_toy_seed, shape_toy_mu_values,
                     shape_bias_tests, shape_bias_mu_true, shape_bias_threshold,
                     shape_bias_toys, shape_bias_param_grid,
                     pdf):
    """
    Combined multi-era fit: load, lumi-scale per-era, concatenate, and fit.
    """
    logger.info("=" * 70)
    logger.info("COMBINED ERA FIT: %s", " + ".join(combined_eras))
    logger.info("=" * 70)
    
    all_sig_scores = []
    all_sig_mass = []
    all_sig_w = []
    all_bkg_scores = []
    all_bkg_mass = []
    all_bkg_w = []
    all_bkg_genjets = []
    all_dy_mass = []
    all_dy_score = []
    all_dy_w = []
    all_dy_genjets = []
    
    for era in combined_eras:
        logger.info("Loading era: %s", era)
        files = find_scored_files(outputs_dir, era)
        sig_file = files["signal"]
        bkg_file = files["background"]
        
        if not sig_file or not bkg_file:
            logger.warning("%s: missing signal/background files; skipping", era)
            continue
        
        mass_branch_era = resolve_mass_branch(sig_file, treename, mass_branch)
        if mass_branch_era is None:
            logger.warning("%s: mass branch not found; skipping", era)
            continue
        
        sel_mass = (
            f"({mass_branch_era}>={mass_min} && {mass_branch_era}<={mass_max} "
            f"&& TMath::Finite({mass_branch_era}))"
        )
        sel_finite = f"TMath::Finite({bdt_branch})"
        selection = f"({sel_finite}) && ({sel_mass})"
        
        branches = [bdt_branch, weight_branch, mass_branch_era]
        sig = load_arrays(sig_file, treename, branches, selection)
        bkg = load_arrays(bkg_file, treename, branches, selection)
        
        sig_scores = sig[bdt_branch]
        sig_mass_vals = sig[mass_branch_era]
        sig_w_raw = sanitize_weights(sig[weight_branch], abs_weights)
        
        bkg_scores = bkg[bdt_branch]
        bkg_mass_vals = bkg[mass_branch_era]
        bkg_w_raw = sanitize_weights(bkg[weight_branch], abs_weights)

        bkg_genjet_branch = resolve_genjet_branch(bkg_file, treename, genjet_branch)
        if bkg_genjet_branch is not None:
            bkg_gen = load_arrays(bkg_file, treename, [bkg_genjet_branch], selection)
            bkg_genjets = np.asarray(bkg_gen[bkg_genjet_branch], dtype=np.float64)
        else:
            bkg_genjets = None
        
        logger.info("%s: S events=%d (sum=%.3f), B events=%d (sum=%.3f) [pre-lumi]",
                   era, sig_scores.size, float(sig_w_raw.sum()),
                   bkg_scores.size, float(bkg_w_raw.sum()))
        
        era_lumi = get_era_luminosity_fb(era)
        if (not bool(no_lumi_scaling)) and era_lumi is not None:
            sig_w, bkg_w, _ = apply_lumi_scaling(sig_w_raw, bkg_w_raw, era_lumi, REF_LUMI_FB, logger)
        else:
            sig_w = sig_w_raw
            bkg_w = bkg_w_raw
            if bool(no_lumi_scaling):
                logger.info("%s: luminosity scaling disabled; no scaling applied", era)
            else:
                logger.warning("%s: luminosity unknown; no scaling applied", era)
        
        logger.info("%s: S sum=%.3f, B sum=%.3f [post-lumi]",
                   era, float(sig_w.sum()), float(bkg_w.sum()))
        
        all_sig_scores.append(sig_scores)
        all_sig_mass.append(sig_mass_vals)
        all_sig_w.append(sig_w)
        all_bkg_scores.append(bkg_scores)
        all_bkg_mass.append(bkg_mass_vals)
        all_bkg_w.append(bkg_w)
        if bkg_genjets is not None:
            all_bkg_genjets.append(bkg_genjets)

        dy_file = find_scored_process_file(outputs_dir, era, "DY")
        if dy_file and mass_branch_era is not None:
            dy_genjet_branch = resolve_genjet_branch(dy_file, treename, genjet_branch)
            dy_cols = [bdt_branch, weight_branch, mass_branch_era]
            if dy_genjet_branch is not None:
                dy_cols.append(dy_genjet_branch)
            dy_arr = load_arrays(dy_file, treename, dy_cols, selection)
            dy_scores = np.asarray(dy_arr[bdt_branch], dtype=np.float64)
            dy_mass = np.asarray(dy_arr[mass_branch_era], dtype=np.float64)
            dy_w_raw = sanitize_weights(dy_arr[weight_branch], abs_weights)
            if (not bool(no_lumi_scaling)) and era_lumi is not None:
                dy_w, _, _ = apply_lumi_scaling(dy_w_raw, dy_w_raw, era_lumi, REF_LUMI_FB, None)
            else:
                dy_w = dy_w_raw
            if dy_genjet_branch is not None:
                dy_genjets = np.asarray(dy_arr[dy_genjet_branch], dtype=np.float64)
                all_dy_mass.append(dy_mass)
                all_dy_score.append(dy_scores)
                all_dy_w.append(dy_w)
                all_dy_genjets.append(dy_genjets)
    
    if not all_sig_scores:
        logger.warning("No valid eras for combined fit")
        return
    
    sig_scores_comb = np.concatenate(all_sig_scores)
    sig_mass_comb = np.concatenate(all_sig_mass)
    sig_w_comb = np.concatenate(all_sig_w)
    bkg_scores_comb = np.concatenate(all_bkg_scores)
    bkg_mass_comb = np.concatenate(all_bkg_mass)
    bkg_w_comb = np.concatenate(all_bkg_w)
    bkg_genjets_comb = np.concatenate(all_bkg_genjets) if all_bkg_genjets else None
    dy_mass_comb = np.concatenate(all_dy_mass) if all_dy_mass else None
    dy_score_comb = np.concatenate(all_dy_score) if all_dy_score else None
    dy_w_comb = np.concatenate(all_dy_w) if all_dy_w else None
    dy_genjets_comb = np.concatenate(all_dy_genjets) if all_dy_genjets else None
    
    logger.info("COMBINED: S events=%d (sum=%.3f), B events=%d (sum=%.3f)",
               sig_scores_comb.size, float(sig_w_comb.sum()),
               bkg_scores_comb.size, float(bkg_w_comb.sum()))
    
    fpr, tpr = compute_roc(sig_scores_comb, sig_w_comb, bkg_scores_comb, bkg_w_comb)
    auc = compute_auc(fpr, tpr)
    plot_roc(pdf, f"Combined ({'+'.join(combined_eras)}): ROC (weighted)", fpr, tpr, auc)
    logger.info("COMBINED: ROC AUC = %.4f", auc)
    
    scan = optimize_multibin_boundaries(
        sig_scores_comb, sig_w_comb, bkg_scores_comb, bkg_w_comb,
        min_n=multibin_min_n, max_n=multibin_max_n,
        n_candidates=multibin_candidates,
        min_sig=multibin_min_sig, min_bkg=multibin_min_bkg,
        enforce_monotonic=(not multibin_no_monotonic),
        objective="asimov",
        sig_mass=sig_mass_comb, bkg_mass=bkg_mass_comb,
        mass_min=mass_min, mass_max=mass_max,
        mass_nbins=shape_mass_nbins,
        poly_order=int(np.clip(shape_poly_order, 0, 2)),
        bkg_norm_unc=float(shape_bkg_norm_unc),
        poly_unc=float(shape_poly_unc),
        sig_norm_unc=float(shape_sig_norm_unc),
        sig_mass_shift_unc=float(shape_sig_mass_shift_unc),
        sig_resolution_unc=float(shape_sig_resolution_unc),
        mu_max=float(shape_mu_max),
        mu_points=int(shape_mu_points),
        likelihood_topk=20,
        shape_data_source=shape_data_source,
    )
    
    best_multi = scan["best"]
    if best_multi is None:
        logger.warning("COMBINED: no valid BDT binning found")
        return
    
    logger.info("COMBINED: Multi-bin optimum: N=%d, ZA=%.3f",
               best_multi["n_bins"], best_multi["za_total"])
    
    if shape_likelihood:
        poly_order = int(np.clip(shape_poly_order, 0, 2))
        shape_inputs = build_shape_inputs(
            sig_mass_comb, sig_scores_comb, sig_w_comb,
            bkg_mass_comb, bkg_scores_comb, bkg_w_comb,
            best_multi["edges"], mass_min, mass_max, shape_mass_nbins,
            bkg_genjets=bkg_genjets_comb,
            dy_mass=dy_mass_comb,
            dy_score=dy_score_comb,
            dy_w=dy_w_comb,
            dy_genjets=dy_genjets_comb,
            core_use_lowmed=shape_core_use_lowmed,
            core_poly_order=shape_core_poly_order,
                core_sb_threshold=shape_core_sb_threshold,
            eps_prior_scale=shape_eps_prior_scale,
        )
        
        shape_fit = compute_shape_expected_significance(
            shape_inputs,
            poly_order=poly_order,
            bkg_norm_unc=float(shape_bkg_norm_unc),
            poly_unc=float(shape_poly_unc),
            sig_norm_unc=float(shape_sig_norm_unc),
            sig_mass_shift_unc=float(shape_sig_mass_shift_unc),
            sig_resolution_unc=float(shape_sig_resolution_unc),
            dy_ge2_unc=float(dy_ge2_unc),
            mu_max=float(shape_mu_max),
            mu_points=int(shape_mu_points),
            data_override=resolve_shape_data_override(shape_inputs, shape_data_source),
        )
        
        logger.info("COMBINED: Shape-likelihood: Z=%.3f, q0=%.3f, mu_hat=%.3f",
                   shape_fit["zexp"], shape_fit["q0"], shape_fit["mu_hat"])
        
        rows_shape = [
            ["Combined Eras", "+".join(combined_eras)],
            ["Expected Z", f"{shape_fit['zexp']:.4f}"],
            ["q0", f"{shape_fit['q0']:.4f}"],
            ["mu_hat", f"{shape_fit['mu_hat']:.4f}"],
            ["N_BDT bins", f"{shape_inputs['n_bdt']}"],
            ["N_mass bins", f"{shape_inputs['n_mass']}"],
            ["poly order", f"{poly_order}"],
            ["DY >=2j unc", f"{dy_ge2_unc:.3f}"],
        ]
        plot_table(pdf,
                  f"Combined ({'+'.join(combined_eras)}): Shape-likelihood summary",
                  ["Quantity", "Value"],
                  rows_shape,
                  col_align=["left", "right"])
        
        rows = []
        for i, (lo, hi) in enumerate(zip(best_multi["edges"][:-1], best_multi["edges"][1:])):
            rows.append([
                str(i + 1),
                f"{lo:.3f}-{hi:.3f}",
                f"{best_multi['s_bins'][i]:.3f}",
                f"{best_multi['b_bins'][i]:.3f}",
                format_ratio(best_multi["sb_bins"][i]),
            ])
        plot_table(pdf,
                  f"Combined ({'+'.join(combined_eras)}): BDT bin yields",
                  ["Bin", "BDT range", "S", "B", "S/B"],
                  rows,
                  col_align=["right", "center", "right", "right", "right"])

        if int(shape_toys) > 0:
            toy_summary = run_shape_toy_ensembles(
                shape_inputs,
                n_toys=int(shape_toys),
                injected_mus=shape_toy_mu_values,
                poly_order=poly_order,
                bkg_norm_unc=float(shape_bkg_norm_unc),
                poly_unc=float(shape_poly_unc),
                sig_norm_unc=float(shape_sig_norm_unc),
                sig_mass_shift_unc=float(shape_sig_mass_shift_unc),
                sig_resolution_unc=float(shape_sig_resolution_unc),
                dy_ge2_unc=float(dy_ge2_unc),
                mu_max=float(shape_mu_max),
                mu_points=int(shape_mu_points),
                seed=int(shape_toy_seed),
            )
            rows_toy = []
            for mu_inj in sorted(toy_summary.keys()):
                item = toy_summary[mu_inj]
                rows_toy.append([
                    f"{mu_inj:.2f}",
                    f"{item['n_toys']}",
                    f"{item['mu_hat_mean']:.4f}",
                    f"{item['mu_hat_std']:.4f}",
                    f"{item['mu_hat_bias']:.4f}",
                    f"{item['z_mean']:.4f}",
                    f"{item['z_std']:.4f}",
                ])
            if rows_toy:
                plot_table(
                    pdf,
                    f"Combined ({'+'.join(combined_eras)}): Toy ensemble validation",
                    ["Injected μ", "N toys", "<μhat>", "σ(μhat)", "Bias", "<Z>", "σ(Z)"],
                    rows_toy,
                    col_align=["right", "right", "right", "right", "right", "right", "right"],
                )

        if bool(shape_bias_tests):
            bias = run_bias_test_matrix(
                shape_inputs,
                mu_true=float(shape_bias_mu_true),
                n_toys=int(shape_bias_toys),
                threshold=float(shape_bias_threshold),
                seed=int(shape_toy_seed),
                param_grid=int(shape_bias_param_grid),
            )
            models = bias["models"]
            rows_bias = []
            for g in models:
                for f in models:
                    item = bias["matrix"][g][f]
                    rows_bias.append([
                        g,
                        f,
                        f"{item['mu_hat_mean']:.4f}",
                        f"{item['mu_hat_std']:.4f}",
                        f"{item['bias']:.4f}",
                        "PASS" if item["pass"] else "FAIL",
                    ])
            plot_table(
                pdf,
                f"Combined ({'+'.join(combined_eras)}): Bias matrix (mu_true={shape_bias_mu_true:.2f}, threshold={shape_bias_threshold:.2f})",
                ["Gen model", "Fit model", "<mu_hat>", "sigma(mu_hat)", "Bias", "Status"],
                rows_bias,
                col_align=["left", "left", "right", "right", "right", "center"],
            )


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(levelname)s] %(message)s",
    )
    logger = logging.getLogger("bdt_analysis")

    args = parse_args()
    logger.info("Starting BDT analysis")
    logger.info("Outputs dir: %s", os.path.abspath(args.outputs_dir))
    logger.info("Mass window: %.1f-%.1f GeV", args.mass_min, args.mass_max)
    logger.info("BDT branch: %s | weight: %s | mass: %s", args.bdt_branch, args.weight_branch, args.mass_branch)
    logger.info("Using abs weights: %s", str(args.abs_weights))
    logger.info("Output PDF: %s", args.output)
    logger.info("KS thresholds: alpha=%g | D<=%.2f%%", args.alpha, args.dthr_pct)
    logger.info(
        "Multi-bin BDT scan: bins=%d-%d | candidates=%d | min(S,B)=(%.3f, %.3f) | monotonic=%s | objective=%s",
        args.multibin_min_n,
        args.multibin_max_n,
        args.multibin_candidates,
        args.multibin_min_sig,
        args.multibin_min_bkg,
        str(not args.multibin_no_monotonic),
        args.multibin_objective,
    )
    logger.info(
        "Shape likelihood: enabled=%s | mass bins=%d | poly order=%d | bkg norm unc=%.3f | poly unc=%.3f | sig norm unc=%.3f | sig mass shift unc=%.3f GeV | sig resolution unc=%.3f",
        str(args.shape_likelihood),
        args.shape_mass_nbins,
        args.shape_poly_order,
        args.shape_bkg_norm_unc,
        args.shape_poly_unc,
        args.shape_sig_norm_unc,
        args.shape_sig_mass_shift_unc,
        args.shape_sig_resolution_unc,
    )
    logger.info("Shape data source: %s", args.shape_data_source)
    logger.info(
        "Core-PDF options: core poly order=%d | low/med core=%s | core S/B threshold=%.4f | eps prior scale=%.3f | DY>=2j unc=%.3f | genjet branch=%s",
        args.shape_core_poly_order,
        str(args.shape_core_use_lowmed),
        args.shape_core_sb_threshold,
        args.shape_eps_prior_scale,
        args.dy_ge2_unc,
        args.genjet_branch,
    )
    toy_mu_values = parse_float_csv(args.shape_toy_mu, default_vals=[0.0, 1.0])
    logger.info("Toy ensembles: n=%d | seed=%d | injected mu=%s",
                args.shape_toys, args.shape_toy_seed, ",".join([f"{x:.2f}" for x in toy_mu_values]))
    logger.info(
        "Bias tests: enabled=%s | mu_true=%.3f | threshold=%.3f | toys=%d | param_grid=%d",
        str(args.shape_bias_tests),
        args.shape_bias_mu_true,
        args.shape_bias_threshold,
        args.shape_bias_toys,
        args.shape_bias_param_grid,
    )
    if args.no_lumi_scaling:
        logger.info("Luminosity scaling: disabled")
    else:
        logger.info("Luminosity scaling (multi-bin yields): reference=%.1f fb^-1", REF_LUMI_FB)
    outputs_dir = os.path.abspath(args.outputs_dir)
    eras = discover_eras(outputs_dir, args.eras)
    if not eras:
        raise SystemExit("No eras discovered. Check --outputs-dir or --eras.")
    logger.info("Discovered eras: %s", ", ".join(eras))

    data_map_raw = parse_group_map(args.data_era, "--data-era") if args.data_era else {}
    mc_map_raw = parse_group_map(args.mc_era, "--mc-era") if args.mc_era else {}
    if data_map_raw or mc_map_raw:
        logger.info("Data/MC inputs provided for eras: %s",
                    ", ".join(sorted(set(list(data_map_raw.keys()) + list(mc_map_raw.keys())))))
        logger.info("Data/MC bins: %d | logy=%s", args.datamc_nbins, str(args.datamc_logy))

    with PdfPages(args.output) as pdf:
        multibin_results_by_era = {}
        for era in eras:
            logger.info("--- Era %s ---", era)
            files = find_scored_files(outputs_dir, era)
            sig_file = files["signal"]
            bkg_file = files["background"]
            all_file = files["all"]

            if not sig_file or not bkg_file:
                print(f"[warn] {era}: missing signal/background scored files; skipping signal-vs-background plots.")
                logger.warning("%s: missing signal/background files", era)

            main_ref = sig_file or bkg_file or all_file
            if not main_ref:
                print(f"[warn] {era}: no scored files found.")
                continue

            mass_branch = resolve_mass_branch(main_ref, args.treename, args.mass_branch)
            if mass_branch is None:
                print(f"[warn] {era}: could not infer mass branch; mass-based plots will be skipped.")
                logger.warning("%s: mass branch not found; skipping mass plots", era)
            else:
                logger.info("%s: mass branch = %s", era, mass_branch)

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
                logger.info(
                    "%s: S events=%d (sum=%.3f), B events=%d (sum=%.3f) [pre-lumi scaling]",
                    era, sig_scores.size, float(sig_w.sum()), bkg_scores.size, float(bkg_w.sum())
                )

                # Apply luminosity scaling to normalize yields across eras for significance calculations
                era_lumi = get_era_luminosity_fb(era)
                if (not bool(args.no_lumi_scaling)) and era_lumi is not None:
                    sig_w, bkg_w, lumi_scale = apply_lumi_scaling(sig_w, bkg_w, era_lumi, REF_LUMI_FB, logger)
                    logger.info(
                        "%s: S events=%d (sum=%.3f), B events=%d (sum=%.3f) [post-lumi scaling]",
                        era, sig_scores.size, float(sig_w.sum()), bkg_scores.size, float(bkg_w.sum())
                    )
                else:
                    if args.no_lumi_scaling:
                        logger.info("%s: luminosity scaling disabled; yields will not be rescaled", era)
                    else:
                        logger.warning("%s: luminosity not known; yields will not be rescaled", era)

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
                logger.info("%s: ROC AUC = %.4f", era, auc)

                # Multi-bin BDT categorization (no hard cut)
                sig_mass_for_opt = sig[mass_branch] if mass_branch else None
                bkg_mass_for_opt = bkg[mass_branch] if mass_branch else None
                scan = optimize_multibin_boundaries(
                    sig_scores,
                    sig_w,
                    bkg_scores,
                    bkg_w,
                    min_n=args.multibin_min_n,
                    max_n=args.multibin_max_n,
                    n_candidates=args.multibin_candidates,
                    min_sig=args.multibin_min_sig,
                    min_bkg=args.multibin_min_bkg,
                    enforce_monotonic=(not args.multibin_no_monotonic),
                    objective=args.multibin_objective,
                    sig_mass=sig_mass_for_opt,
                    bkg_mass=bkg_mass_for_opt,
                    mass_min=args.mass_min,
                    mass_max=args.mass_max,
                    mass_nbins=args.shape_mass_nbins,
                    poly_order=int(np.clip(args.shape_poly_order, 0, 2)),
                    bkg_norm_unc=float(args.shape_bkg_norm_unc),
                    poly_unc=float(args.shape_poly_unc),
                    sig_norm_unc=float(args.shape_sig_norm_unc),
                    sig_mass_shift_unc=float(args.shape_sig_mass_shift_unc),
                    sig_resolution_unc=float(args.shape_sig_resolution_unc),
                    mu_max=float(args.shape_mu_max),
                    mu_points=int(args.shape_mu_points),
                    likelihood_topk=int(args.multibin_likelihood_topk),
                    shape_data_source=args.shape_data_source,
                )
                best_multi = scan["best"]
                multibin_results_by_era[era] = {
                    "best": best_multi,
                    "best_by_n": scan["best_by_n"],
                    "scan_stats": scan["scan_stats"],
                    "metric_label": scan.get("metric_label", "ZA(total)"),
                }

                if scan.get("objective_requested") == "likelihood" and scan.get("objective_used") != "likelihood":
                    logger.warning(
                        "%s: Requested likelihood objective but mass branch is unavailable; falling back to ZA(total)",
                        era,
                    )

                if best_multi is not None:
                    metric_label = best_multi.get("metric_label", scan.get("metric_label", "ZA(total)"))
                    metric_value = float(best_multi.get("opt_metric", best_multi.get("za_total", 0.0)))
                    edges_txt = ", ".join([f"{x:.3f}" for x in best_multi["edges"]])
                    logger.info(
                        "%s: Multi-bin optimum: N=%d, %s=%.3f, edges=[%s]",
                        era,
                        best_multi["n_bins"],
                        metric_label,
                        metric_value,
                        edges_txt,
                    )
                    for n_bins in range(args.multibin_min_n, args.multibin_max_n + 1):
                        st = scan["scan_stats"].get(n_bins, {"scanned": 0, "valid": 0})
                        if n_bins in scan["best_by_n"]:
                            n_best = scan["best_by_n"][n_bins]
                            n_metric = float(n_best.get("opt_metric", n_best.get("za_total", 0.0)))
                            logger.info(
                                "%s: N=%d scanned=%d valid=%d best_%s=%.3f",
                                era,
                                n_bins,
                                st["scanned"],
                                st["valid"],
                                metric_label,
                                n_metric,
                            )
                        else:
                            logger.info(
                                "%s: N=%d scanned=%d valid=%d (no valid monotonic/stat-stable config)",
                                era,
                                n_bins,
                                st["scanned"],
                                st["valid"],
                            )

                    rows = []
                    for i, (lo, hi) in enumerate(zip(best_multi["edges"][:-1], best_multi["edges"][1:])):
                        rows.append([
                            str(i + 1),
                            f"{lo:.3f}-{hi:.3f}",
                            f"{best_multi['s_bins'][i]:.3f}",
                            f"{best_multi['b_bins'][i]:.3f}",
                            format_ratio(best_multi["sb_bins"][i]),
                            f"{best_multi['zapprox_bins'][i]:.3f}",
                        ])
                    plot_table(
                        pdf,
                        f"{era}: Multi-bin BDT categorization (no hard cut)\n"
                        f"Optimized N={best_multi['n_bins']}, {metric_label}={metric_value:.3f}, "
                        f"monotonic S/B={str(not args.multibin_no_monotonic)}",
                        ["Bin", "BDT range", "S", "B", "S/B", "S/sqrt(B)"],
                        rows,
                        col_align=["right", "center", "right", "right", "right", "right"],
                    )

                    if mass_branch:
                        plot_mass_sculpt(
                            pdf,
                            f"{era}: Signal mass by optimized BDT categories",
                            sig[mass_branch],
                            sig_scores,
                            sig_w,
                            best_multi["edges"],
                        )
                        plot_mass_sculpt(
                            pdf,
                            f"{era}: Background mass by optimized BDT categories",
                            bkg[mass_branch],
                            bkg_scores,
                            bkg_w,
                            best_multi["edges"],
                        )
                        plot_mass_sculpt_norm(
                            pdf,
                            f"{era}: Signal normalized mass by optimized BDT categories",
                            sig[mass_branch],
                            sig_scores,
                            sig_w,
                            best_multi["edges"],
                            weighted_corr(sig[mass_branch], sig_scores, sig_w),
                        )
                        plot_mass_sculpt_norm(
                            pdf,
                            f"{era}: Background normalized mass by optimized BDT categories",
                            bkg[mass_branch],
                            bkg_scores,
                            bkg_w,
                            best_multi["edges"],
                            weighted_corr(bkg[mass_branch], bkg_scores, bkg_w),
                        )

                        if args.shape_likelihood:
                            poly_order = int(np.clip(args.shape_poly_order, 0, 2))
                            bkg_genjet_branch = resolve_genjet_branch(bkg_file, args.treename, args.genjet_branch)
                            bkg_genjets = None
                            if bkg_genjet_branch is not None:
                                bkg_gen = load_arrays(bkg_file, args.treename, [bkg_genjet_branch], selection)
                                bkg_genjets = np.asarray(bkg_gen[bkg_genjet_branch], dtype=np.float64)

                            dy_mass_vals = None
                            dy_scores_vals = None
                            dy_w_vals = None
                            dy_genjets_vals = None
                            dy_file = find_scored_process_file(outputs_dir, era, "DY")
                            if dy_file and mass_branch:
                                dy_genjet_branch = resolve_genjet_branch(dy_file, args.treename, args.genjet_branch)
                                dy_cols = [args.bdt_branch, args.weight_branch, mass_branch]
                                if dy_genjet_branch is not None:
                                    dy_cols.append(dy_genjet_branch)
                                dy_arr = load_arrays(dy_file, args.treename, dy_cols, selection)
                                dy_scores_vals = np.asarray(dy_arr[args.bdt_branch], dtype=np.float64)
                                dy_mass_vals = np.asarray(dy_arr[mass_branch], dtype=np.float64)
                                dy_w_vals = sanitize_weights(dy_arr[args.weight_branch], args.abs_weights)
                                if dy_genjet_branch is not None:
                                    dy_genjets_vals = np.asarray(dy_arr[dy_genjet_branch], dtype=np.float64)

                            shape_inputs = build_shape_inputs(
                                sig[mass_branch],
                                sig_scores,
                                sig_w,
                                bkg[mass_branch],
                                bkg_scores,
                                bkg_w,
                                best_multi["edges"],
                                args.mass_min,
                                args.mass_max,
                                args.shape_mass_nbins,
                                bkg_genjets=bkg_genjets,
                                dy_mass=dy_mass_vals,
                                dy_score=dy_scores_vals,
                                dy_w=dy_w_vals,
                                dy_genjets=dy_genjets_vals,
                                core_use_lowmed=args.shape_core_use_lowmed,
                                core_poly_order=args.shape_core_poly_order,
                                core_sb_threshold=args.shape_core_sb_threshold,
                                eps_prior_scale=args.shape_eps_prior_scale,
                            )
                            shape_fit = compute_shape_expected_significance(
                                shape_inputs,
                                poly_order=poly_order,
                                bkg_norm_unc=float(args.shape_bkg_norm_unc),
                                poly_unc=float(args.shape_poly_unc),
                                sig_norm_unc=float(args.shape_sig_norm_unc),
                                sig_mass_shift_unc=float(args.shape_sig_mass_shift_unc),
                                sig_resolution_unc=float(args.shape_sig_resolution_unc),
                                dy_ge2_unc=float(args.dy_ge2_unc),
                                mu_max=float(args.shape_mu_max),
                                mu_points=int(args.shape_mu_points),
                                data_override=resolve_shape_data_override(shape_inputs, args.shape_data_source),
                            )
                            logger.info(
                                "%s: Shape-likelihood Asimov: Z=%.3f, q0=%.3f, mu_hat=%.3f",
                                era,
                                shape_fit["zexp"],
                                shape_fit["q0"],
                                shape_fit["mu_hat"],
                            )

                            rows_shape = [
                                ["Expected Z", f"{shape_fit['zexp']:.4f}"],
                                ["q0", f"{shape_fit['q0']:.4f}"],
                                ["mu_hat", f"{shape_fit['mu_hat']:.4f}"],
                                ["N_BDT bins", f"{shape_inputs['n_bdt']}"],
                                ["N_mass bins", f"{shape_inputs['n_mass']}"],
                                ["poly order", f"{poly_order}"],
                                ["bkg norm unc", f"{args.shape_bkg_norm_unc:.3f}"],
                                ["poly unc", f"{args.shape_poly_unc:.3f}"],
                                ["core poly order", f"{args.shape_core_poly_order}"],
                                ["core from low/med bins", str(bool(args.shape_core_use_lowmed))],
                                ["DY >=2j unc", f"{args.dy_ge2_unc:.3f}"],
                                ["sig norm unc", f"{args.shape_sig_norm_unc:.3f}"],
                                ["sig mass shift unc [GeV]", f"{args.shape_sig_mass_shift_unc:.3f}"],
                                ["sig resolution unc", f"{args.shape_sig_resolution_unc:.3f}"],
                            ]
                            plot_table(
                                pdf,
                                f"{era}: Simultaneous shape-based profile likelihood summary",
                                ["Quantity", "Value"],
                                rows_shape,
                                col_align=["left", "right"],
                            )

                            rows_yields = []
                            sb_model = shape_fit["model_sb"]
                            b_model = shape_fit["model_b"]
                            for i, (lo, hi) in enumerate(zip(best_multi["edges"][:-1], best_multi["edges"][1:])):
                                s_bin = float(shape_inputs["s_hist"][i].sum())
                                b_bin = float(shape_inputs["b_hist"][i].sum())
                                sb_fit = float(sb_model[i].sum())
                                b_fit = float(b_model[i].sum())
                                rows_yields.append([
                                    str(i + 1),
                                    f"{lo:.3f}-{hi:.3f}",
                                    f"{s_bin:.3f}",
                                    f"{b_bin:.3f}",
                                    f"{b_fit:.3f}",
                                    f"{sb_fit:.3f}",
                                ])
                            plot_table(
                                pdf,
                                f"{era}: Shape-likelihood bin yields (Asimov)",
                                ["Bin", "BDT", "S truth", "B truth", "B-only fit", "S+B fit"],
                                rows_yields,
                                col_align=["right", "center", "right", "right", "right", "right"],
                            )

                            profile = np.asarray(shape_fit["profile"], dtype=np.float64)
                            fig, ax = plt.subplots(figsize=(9, 6))
                            ax.plot(profile[:, 0], 2.0 * (profile[:, 1] - np.min(profile[:, 1])), linewidth=2)
                            ax.axvline(shape_fit["mu_hat"], color="tab:red", linestyle="--", linewidth=1.5)
                            ax.set_xlabel("Signal strength μ")
                            ax.set_ylabel("Δ(-2 ln L)")
                            ax.set_title(f"{era}: Profile likelihood scan (simultaneous shape fit)")
                            ax.grid(alpha=0.3)
                            ax.text(
                                0.02,
                                0.95,
                                f"Zexp = {shape_fit['zexp']:.3f}\nq0 = {shape_fit['q0']:.3f}\nμhat = {shape_fit['mu_hat']:.3f}",
                                transform=ax.transAxes,
                                va="top",
                                fontsize=10,
                                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                            )
                            pdf.savefig(fig)
                            plt.close(fig)

                            plot_shape_fit_overlays(
                                pdf,
                                f"{era}: Asimov mass templates by BDT bin (simultaneous fit)",
                                shape_inputs["mass_edges"],
                                best_multi["edges"],
                                shape_inputs["asimov"],
                                shape_fit["model_b"],
                                shape_fit["model_sb"],
                            )

                            if int(args.shape_toys) > 0:
                                toy_summary = run_shape_toy_ensembles(
                                    shape_inputs,
                                    n_toys=int(args.shape_toys),
                                    injected_mus=toy_mu_values,
                                    poly_order=poly_order,
                                    bkg_norm_unc=float(args.shape_bkg_norm_unc),
                                    poly_unc=float(args.shape_poly_unc),
                                    sig_norm_unc=float(args.shape_sig_norm_unc),
                                    sig_mass_shift_unc=float(args.shape_sig_mass_shift_unc),
                                    sig_resolution_unc=float(args.shape_sig_resolution_unc),
                                    dy_ge2_unc=float(args.dy_ge2_unc),
                                    mu_max=float(args.shape_mu_max),
                                    mu_points=int(args.shape_mu_points),
                                    seed=int(args.shape_toy_seed),
                                )
                                rows_toy = []
                                for mu_inj in sorted(toy_summary.keys()):
                                    item = toy_summary[mu_inj]
                                    rows_toy.append([
                                        f"{mu_inj:.2f}",
                                        f"{item['n_toys']}",
                                        f"{item['mu_hat_mean']:.4f}",
                                        f"{item['mu_hat_std']:.4f}",
                                        f"{item['mu_hat_bias']:.4f}",
                                        f"{item['z_mean']:.4f}",
                                        f"{item['z_std']:.4f}",
                                    ])
                                if rows_toy:
                                    plot_table(
                                        pdf,
                                        f"{era}: Toy ensemble validation (core-PDF likelihood)",
                                        ["Injected μ", "N toys", "<μhat>", "σ(μhat)", "Bias", "<Z>", "σ(Z)"],
                                        rows_toy,
                                        col_align=["right", "right", "right", "right", "right", "right", "right"],
                                    )

                            if bool(args.shape_bias_tests):
                                bias = run_bias_test_matrix(
                                    shape_inputs,
                                    mu_true=float(args.shape_bias_mu_true),
                                    n_toys=int(args.shape_bias_toys),
                                    threshold=float(args.shape_bias_threshold),
                                    seed=int(args.shape_toy_seed),
                                    param_grid=int(args.shape_bias_param_grid),
                                )
                                models = bias["models"]
                                rows_bias = []
                                for g in models:
                                    for f in models:
                                        item = bias["matrix"][g][f]
                                        rows_bias.append([
                                            g,
                                            f,
                                            f"{item['mu_hat_mean']:.4f}",
                                            f"{item['mu_hat_std']:.4f}",
                                            f"{item['bias']:.4f}",
                                            "PASS" if item["pass"] else "FAIL",
                                        ])
                                plot_table(
                                    pdf,
                                    f"{era}: Bias matrix (mu_true={args.shape_bias_mu_true:.2f}, threshold={args.shape_bias_threshold:.2f})",
                                    ["Gen model", "Fit model", "<mu_hat>", "sigma(mu_hat)", "Bias", "Status"],
                                    rows_bias,
                                    col_align=["left", "left", "right", "right", "right", "center"],
                                )
                else:
                    logger.warning(
                        "%s: Multi-bin optimization found no valid configuration (check min yields / candidates / monotonicity)",
                        era,
                    )

                # Expected signal yield per optimized BDT bin (mass window)
                if best_multi is None:
                    continue
                bins_yield = best_multi["edges"]
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
                    f"{era}: Expected yield per optimized BDT bin ({args.mass_min:.0f}-{args.mass_max:.0f} GeV)\n"
                    f"S total={s_total:.2f}, B total={b_total:.2f}",
                    ["BDT bin", "S (bin)", "S frac", "S>=lo", "B (bin)", "S/B (>=lo)"],
                    rows,
                    col_align=["center", "right", "right", "right", "right", "right"],
                )
                logger.info("%s: Yield table bins=%d (optimized bins)", era, len(rows))
                for row in rows:
                    logger.info(
                        "%s: bin %s | S=%s | Sfrac=%s | S>=lo=%s | B=%s | S/B=%s",
                        era, row[0], row[1], row[2], row[3], row[4], row[5]
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
                        best_multi["edges"],
                    )
                    plot_mass_sculpt(
                        pdf,
                        f"{era}: Background mass sculpting (GeV)",
                        bkg_mass,
                        bkg_scores,
                        bkg_w,
                        best_multi["edges"],
                    )
                    corr_sig = weighted_corr(sig_mass, sig_scores, sig_w)
                    corr_bkg = weighted_corr(bkg_mass, bkg_scores, bkg_w)
                    plot_mass_sculpt_norm(
                        pdf,
                        f"{era}: Signal mass sculpting (normalized shapes)",
                        sig_mass,
                        sig_scores,
                        sig_w,
                        best_multi["edges"],
                        corr_sig,
                    )
                    plot_mass_sculpt_norm(
                        pdf,
                        f"{era}: Background mass sculpting (normalized shapes)",
                        bkg_mass,
                        bkg_scores,
                        bkg_w,
                        best_multi["edges"],
                        corr_bkg,
                    )
                    if np.isfinite(corr_sig) or np.isfinite(corr_bkg):
                        logger.info(
                            "%s: corr(mass, BDT) signal=%.4f background=%.4f",
                            era, corr_sig, corr_bkg
                        )
                    sig_mean, sig_rms = weighted_mean_std(sig_mass, sig_w)
                    bkg_mean, bkg_rms = weighted_mean_std(bkg_mass, bkg_w)
                    logger.info(
                        "%s: Mass summary S mean=%.2f RMS=%.2f | B mean=%.2f RMS=%.2f",
                        era, sig_mean, sig_rms, bkg_mean, bkg_rms
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
                if best_multi is not None:
                    s_total_bins = float(sig_feat_map["weight"].sum())
                    for lo, hi in zip(best_multi["edges"][:-1], best_multi["edges"][1:]):
                        mask = (sig_feat_map["score"] >= lo) & (sig_feat_map["score"] < hi)
                        s_bin = float(sig_feat_map["weight"][mask].sum())
                        frac = 100.0 * s_bin / s_total_bins if s_total_bins > 0 else 0.0
                        bin_stats.append((f"{lo:.2f}-{hi:.2f}", frac))
                plot_feature_grid(
                    pdf,
                    f"{era}: Signal distributions by BDT bin",
                    feature_grid[:4],
                    sig_feat_map,
                    best_multi["edges"] if best_multi is not None else [0.0, 1.0],
                    mass_branch=mass_branch,
                    bin_stats=bin_stats,
                )
                logger.info("%s: Signal feature grid: %s", era, ", ".join(feature_grid[:4]))
                for feat in feature_grid[:4]:
                    if feat is None:
                        continue
                    corr = weighted_corr(sig_feat_map[feat], sig_feat_map["score"], sig_feat_map["weight"])
                    if np.isfinite(corr):
                        logger.info("%s: %s corr(score)=%.3f", era, feat, corr)

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
                    for row in rows:
                        logger.info(
                            "%s: Fold %s KS=%s N=%s alpha=%s dthr=%s pass=%s",
                            era, row[0], row[1], row[2], row[3], row[4], row[5]
                        )
                else:
                    print(f"[warn] {era}: fold branch missing; KS/CDF skipped.")
                    logger.warning("%s: fold branch missing; KS/CDF skipped", era)

            # Optional Data/MC summary for this era
            if era in data_map_raw and era in mc_map_raw:
                data_files = []
                for v in data_map_raw.get(era, []):
                    data_files.extend([p for p in v.split(",") if p])
                mc_entries = mc_map_raw.get(era, [])
                
                logger.info("%s: Data/MC inputs: %d data file(s), %d MC process(es)",
                           era, len(data_files), len(mc_entries))
                for i, f in enumerate(data_files[:3]):
                    logger.info("%s:   Data file %d: %s", era, i+1, f)
                if len(data_files) > 3:
                    logger.info("%s:   ... and %d more data file(s)", era, len(data_files) - 3)

                if data_files and mc_entries:
                    dm_mass_branch = resolve_mass_branch(data_files[0], args.treename, args.mass_branch)
                    if dm_mass_branch is None:
                        logger.warning("%s: Data/MC mass branch not found; skipping Data/MC plots", era)
                    else:
                        data_cols = set(list_branches(data_files[0], args.treename))
                        logger.info("%s: Data file branches: %s", era, ", ".join(sorted(data_cols)[:5]))
                        if args.bdt_branch not in data_cols:
                            logger.warning(
                                "%s: Data/MC missing bdt branch '%s' in data; skipping Data/MC plots",
                                era, args.bdt_branch
                            )
                            data_cols = None

                        dm_mass_sel = (
                            f"({dm_mass_branch}>={args.mass_min} && {dm_mass_branch}<={args.mass_max} "
                            f"&& TMath::Finite({dm_mass_branch}))"
                        )
                        dm_sel = f"{dm_mass_sel}"
                        if data_cols is not None:
                            dm_sel = f"TMath::Finite({args.bdt_branch}) && {dm_mass_sel}"
                        
                        logger.info("%s: Data/MC selection: %s", era, dm_sel)
                        logger.info("%s: Mass branch: %s [%.0f, %.0f] GeV", era, dm_mass_branch, args.mass_min, args.mass_max)

                        if data_cols is None:
                            continue

                        data_arr = load_arrays(
                            data_files,
                            args.treename,
                            [args.bdt_branch],
                            dm_sel,
                        )
                        data_scores = data_arr[args.bdt_branch]
                        
                        data_mean, data_std = weighted_mean_std(data_scores, np.ones_like(data_scores))
                        logger.info("%s: Data after selection: n_events=%d, mean(BDT)=%.4f, std(BDT)=%.4f",
                                   era, data_scores.size, data_mean, data_std)

                        mc_groups = []
                        mc_full_groups = []
                        for entry in mc_entries:
                            if ":" not in entry:
                                raise SystemExit(f"[error] --mc-era needs ERA:PROC:files, got {entry}")
                            proc, files = entry.split(":", 1)
                            flist = [p for p in files.split(",") if p]
                            logger.info("%s: MC process '%s': %d file(s)", era, proc, len(flist))
                            for i, f in enumerate(flist[:2]):
                                logger.info("%s:   MC file %d: %s", era, i+1, f)
                            if len(flist) > 2:
                                logger.info("%s:   ... and %d more file(s)", era, len(flist) - 2)
                            
                            if flist:
                                mc_cols = set(list_branches(flist[0], args.treename))
                                if args.bdt_branch not in mc_cols:
                                    logger.warning(
                                        "%s: Data/MC missing bdt branch '%s' in MC file %s; skipping Data/MC plots",
                                        era, args.bdt_branch, flist[0]
                                    )
                                    mc_groups = []
                                    break
                            mc_arr = load_arrays(
                                flist,
                                args.treename,
                                [args.bdt_branch, args.weight_branch, dm_mass_branch],
                                dm_sel,
                            )
                            mc_scores = mc_arr[args.bdt_branch]
                            mc_weights = sanitize_weights(mc_arr[args.weight_branch], args.abs_weights)
                            mc_mass = mc_arr[dm_mass_branch]
                            mc_sum = float(mc_weights.sum())
                            mc_mean, mc_std = weighted_mean_std(mc_scores, mc_weights)
                            logger.info("%s: MC '%s' after selection: n_events=%d, sum(weight)=%.3f, mean(BDT)=%.4f, std(BDT)=%.4f",
                                       era, proc, mc_scores.size, mc_sum, mc_mean, mc_std)
                            mc_groups.append((proc, mc_scores, mc_weights))
                            mc_full_groups.append((proc, mc_scores, mc_mass, mc_weights))

                        if not mc_groups:
                            continue

                        # Log pre-rescale statistics
                        total_mc_pre = sum(float(w.sum()) for _, _, w in mc_groups)
                        logger.info("%s: Total MC (pre-rescale): sum(weight)=%.3f, data/MC ratio=%.6f",
                                   era, total_mc_pre, data_scores.size / max(total_mc_pre, 1e-12))

                        # Rescale MC to data counts if requested
                        scale_used = 1.0
                        if args.lumi is not None:
                            mc_groups, scale_used = rescale_mc_to_data(data_scores.size, mc_groups, logger)
                            logger.info("%s: Rescaled MC groups with lumi=%.1f fb^-1; scale=%.6f", 
                                        era, args.lumi, scale_used)
                            total_mc_post = sum(float(w.sum()) for _, _, w in mc_groups)
                            logger.info("%s: Total MC (post-rescale): sum(weight)=%.3f, data/MC ratio=%.6f",
                                       era, total_mc_post, data_scores.size / max(total_mc_post, 1e-12))
                            if scale_used != 1.0:
                                mc_full_groups = [
                                    (proc, scores, masses, weights * scale_used)
                                    for proc, scores, masses, weights in mc_full_groups
                                ]
                        
                        bins = np.linspace(0, 1, args.datamc_nbins + 1)
                        
                        # Compute and log bin-by-bin statistics
                        data_hist, _ = np.histogram(data_scores, bins=bins)
                        mc_hist_total = np.zeros(len(bins) - 1, dtype=np.float64)
                        for _, mc_scores, mc_weights in mc_groups:
                            mc_hist, _ = np.histogram(mc_scores, bins=bins, weights=mc_weights)
                            mc_hist_total += mc_hist
                        
                        logger.info("%s: Bin-by-bin data/MC comparison (first 5 + last 2 bins):", era)
                        for i in range(min(5, len(bins)-1)):
                            ratio = data_hist[i] / max(mc_hist_total[i], 1e-12)
                            logger.info("%s:   bin[%.2f-%.2f]: data=%d, mc=%.1f, ratio=%.3f",
                                       era, bins[i], bins[i+1], int(data_hist[i]), mc_hist_total[i], ratio)
                        if len(bins) > 7:
                            for i in range(len(bins)-3, len(bins)-1):
                                ratio = data_hist[i] / max(mc_hist_total[i], 1e-12)
                                logger.info("%s:   bin[%.2f-%.2f]: data=%d, mc=%.1f, ratio=%.3f",
                                           era, bins[i], bins[i+1], int(data_hist[i]), mc_hist_total[i], ratio)
                        
                        plot_datamc(
                            pdf,
                            f"{era}: Data/MC BDT score ({args.mass_min:.0f}-{args.mass_max:.0f} GeV)",
                            data_scores,
                            mc_groups,
                            bins,
                            args.datamc_logy,
                        )
                        logger.info("%s: Data/MC summary in %.0f-%.0f GeV", era, args.mass_min, args.mass_max)
                        logger.info("%s: Data entries=%d", era, data_scores.size)
                        total_mc = 0.0
                        for proc, mc_scores, mc_weights in mc_groups:
                            proc_sum = float(mc_weights.sum())
                            total_mc += proc_sum
                            proc_ratio = data_scores.size / max(proc_sum, 1e-12)
                            logger.info(
                                "%s: MC %s entries=%d sum=%.3f ratio=%.6f",
                                era, proc, mc_scores.size, proc_sum, proc_ratio
                            )
                        if total_mc > 0:
                            overall_ratio = data_scores.size / total_mc
                            logger.info("%s: Data/MC total ratio = %.6f (data=%.0f, mc_sum=%.3f)",
                                       era, overall_ratio, data_scores.size, total_mc)
            elif data_map_raw or mc_map_raw:
                logger.info("%s: Data/MC not provided for this era; skipping", era)

            # Feature importance table (per-era)
            path, items = read_feature_importance(outputs_dir, era)
            if items:
                rows = [[name, f"{pct:.2f}%"] for name, pct in items[:10]]
                plot_table(
                    pdf,
                    f"{era}: Top-10 feature importance (gain)\n{os.path.basename(path)}",
                    ["Feature", "Gain %"],
                    rows,
                )
                logger.info("%s: Feature importance file: %s", era, os.path.basename(path))
                for name, pct in items[:10]:
                    logger.info("%s: Feature %s gain=%.2f%%", era, name, pct)

            # Optional retrain-per-fold KS (uses DataPreparation_Files)
            if args.retrain_ks:
                sig_all = os.path.join(args.ks_dir, f"{era}_all_signal.root")
                bkg_all = os.path.join(args.ks_dir, f"{era}_all_background.root")
                if os.path.exists(sig_all) and os.path.exists(bkg_all):
                    try:
                        from xgboost import XGBClassifier
                        from tmva101_Training import load_folds
                    except Exception as exc:
                        logger.warning("%s: retrain KS skipped (missing deps): %s", era, exc)
                    else:
                        logger.info("%s: retrain-per-fold KS using %s", era, args.ks_dir)
                        ks_rows = []
                        for k in range(4):
                            folds_tr = ((k + 1) % 4, (k + 2) % 4)
                            folds_va = (k,)
                            folds_te = ((k + 3) % 4,)
                            (x_tr, y_tr, w_tr), (x_va, y_va, w_va), _ = load_folds(
                                [sig_all], [bkg_all],
                                folds_train=folds_tr, folds_valid=folds_va, folds_test=folds_te,
                                sig_extra_weight_branch=None
                            )

                            bdt = XGBClassifier(
                                objective="binary:logistic",
                                eval_metric="logloss",
                                tree_method="hist",
                                missing=np.nan,
                                max_bin=30,
                                max_depth=4,
                                learning_rate=0.1,
                                n_estimators=1000,
                                booster="gbtree",
                                n_jobs=8,
                                gamma=0,
                                min_child_weight=0.03,
                                max_delta_step=0,
                                subsample=0.5,
                                colsample_bytree=1,
                                colsample_bylevel=1,
                                reg_alpha=0,
                                reg_lambda=1,
                                scale_pos_weight=0.7,
                                base_score=0.5,
                                random_state=0,
                                verbosity=0,
                            )
                            bdt.fit(x_tr, y_tr, sample_weight=w_tr)
                            tr_s = bdt.predict_proba(np.asarray(x_tr[y_tr == 1], dtype=np.float32))[:, 1]
                            va_s = bdt.predict_proba(np.asarray(x_va[y_va == 1], dtype=np.float32))[:, 1]

                            D, p = ks_2samp(va_s, tr_s)
                            D_pct = 100.0 * D
                            passed = (p >= args.alpha) and (D_pct <= args.dthr_pct)
                            ks_rows.append([
                                str(k), f"{D_pct:.2f}%", f"{len(va_s):d}",
                                f"{args.alpha:.3g}", f"{args.dthr_pct:.2f}%",
                                "TRUE" if passed else "FALSE",
                            ])
                            logger.info(
                                "%s: retrain KS fold %d D=%.3f%% p=%.3g pass=%s",
                                era, k, D_pct, p, "TRUE" if passed else "FALSE"
                            )
                        plot_table(
                            pdf,
                            f"{era}: Retrain-per-fold KS summary",
                            ["Fold", "KS statistic", "N events", "alpha", "Pass threshold", "KS pass"],
                            ks_rows,
                        )
                else:
                    logger.warning("%s: retrain KS inputs missing in %s", era, args.ks_dir)

        # Optional combined-era fit
        if args.combined_eras and len(args.combined_eras) > 0:
            combined_era_fit(
                outputs_dir=outputs_dir,
                combined_eras=args.combined_eras,
                treename=args.treename,
                bdt_branch=args.bdt_branch,
                weight_branch=args.weight_branch,
                mass_branch=args.mass_branch,
                mass_min=args.mass_min,
                mass_max=args.mass_max,
                abs_weights=args.abs_weights,
                logger=logger,
                multibin_min_n=args.multibin_min_n,
                multibin_max_n=args.multibin_max_n,
                multibin_candidates=args.multibin_candidates,
                multibin_min_sig=args.multibin_min_sig,
                multibin_min_bkg=args.multibin_min_bkg,
                multibin_no_monotonic=args.multibin_no_monotonic,
                shape_likelihood=args.shape_likelihood,
                shape_mass_nbins=args.shape_mass_nbins,
                shape_poly_order=args.shape_poly_order,
                shape_bkg_norm_unc=args.shape_bkg_norm_unc,
                shape_poly_unc=args.shape_poly_unc,
                shape_sig_norm_unc=args.shape_sig_norm_unc,
                shape_sig_mass_shift_unc=args.shape_sig_mass_shift_unc,
                shape_sig_resolution_unc=args.shape_sig_resolution_unc,
                shape_mu_max=args.shape_mu_max,
                shape_mu_points=args.shape_mu_points,
                shape_data_source=args.shape_data_source,
                shape_core_poly_order=args.shape_core_poly_order,
                shape_core_use_lowmed=args.shape_core_use_lowmed,
                shape_core_sb_threshold=args.shape_core_sb_threshold,
                shape_eps_prior_scale=args.shape_eps_prior_scale,
                dy_ge2_unc=args.dy_ge2_unc,
                genjet_branch=args.genjet_branch,
                no_lumi_scaling=args.no_lumi_scaling,
                shape_toys=args.shape_toys,
                shape_toy_seed=args.shape_toy_seed,
                shape_toy_mu_values=toy_mu_values,
                shape_bias_tests=args.shape_bias_tests,
                shape_bias_mu_true=args.shape_bias_mu_true,
                shape_bias_threshold=args.shape_bias_threshold,
                shape_bias_toys=args.shape_bias_toys,
                shape_bias_param_grid=args.shape_bias_param_grid,
                pdf=pdf,
            )

        # Cross-era boundary stability summary (same N bins)
        if len(multibin_results_by_era) >= 2:
            rows = []
            summary_metric_label = None
            for n_bins in range(args.multibin_min_n, args.multibin_max_n + 1):
                edges_list = []
                metric_vals = []
                used_eras = []
                for era in eras:
                    era_res = multibin_results_by_era.get(era, {})
                    by_n = era_res.get("best_by_n", {}) if era_res else {}
                    if n_bins in by_n:
                        used_eras.append(era)
                        metric_vals.append(float(by_n[n_bins].get("opt_metric", by_n[n_bins].get("za_total", 0.0))))
                        edges_list.append(by_n[n_bins]["edges"])
                        if summary_metric_label is None:
                            summary_metric_label = by_n[n_bins].get("metric_label", "ZA(total)")
                if len(edges_list) < 2:
                    continue

                internal = np.asarray([e[1:-1] for e in edges_list], dtype=np.float64)
                stds = np.std(internal, axis=0)
                max_std = float(np.max(stds)) if stds.size else np.nan
                mean_std = float(np.mean(stds)) if stds.size else np.nan
                rows.append([
                    str(n_bins),
                    str(len(used_eras)),
                    f"{np.mean(metric_vals):.3f}",
                    f"{np.min(metric_vals):.3f}",
                    f"{np.max(metric_vals):.3f}",
                    f"{mean_std:.4f}",
                    f"{max_std:.4f}",
                ])

            if rows:
                metric_tag = summary_metric_label if summary_metric_label else "ZA"
                plot_table(
                    pdf,
                    "Cross-era stability of optimized BDT category boundaries\n"
                    "(smaller boundary std = more stable across eras)",
                    ["N bins", "#eras", f"<{metric_tag}>", f"{metric_tag} min", f"{metric_tag} max", "mean std(edge)", "max std(edge)"],
                    rows,
                    col_align=["right", "right", "right", "right", "right", "right", "right"],
                )


    # Optional model-based scan (global, not per-era)
    if args.scan_sig and args.scan_bkg:
        logger.info("Running model-based scan")
        ROOT.gSystem.Load("libTMVA")
        if ROOT.gSystem.AccessPathName(args.model_file, ROOT.kReadPermission):
            raise SystemExit(f"Model file not readable: {args.model_file}")
        model = ROOT.TMVA.Experimental.RBDT[""](args.model_name, args.model_file)

        logger.info("Scan model: %s (method %s)", args.model_file, args.model_name)
        logger.info("Scan cuts: %s", args.scan_cuts)

        sig_files = expand_inputs(args.scan_sig)
        bkg_files = expand_inputs(args.scan_bkg)
        if not sig_files or not bkg_files:
            raise SystemExit("No scan inputs found for signal/background.")

        scan_mass_branch = resolve_mass_branch(sig_files[0], args.treename, args.mass_branch)
        if scan_mass_branch is None and bkg_files:
            scan_mass_branch = resolve_mass_branch(bkg_files[0], args.treename, args.mass_branch)
        scan_sel = ""
        if scan_mass_branch:
            scan_sel = (
                f"TMath::Finite({args.bdt_branch}) && "
                f"({scan_mass_branch}>={args.mass_min} && {scan_mass_branch}<={args.mass_max} "
                f"&& TMath::Finite({scan_mass_branch}))"
            )

        sig_arrs = load_arrays(sig_files, args.treename, VARIABLES + [args.weight_branch], scan_sel)
        bkg_arrs = load_arrays(bkg_files, args.treename, VARIABLES + [args.weight_branch], scan_sel)
        xs_s = [sig_arrs[v] for v in VARIABLES]
        xs_b = [bkg_arrs[v] for v in VARIABLES]
        w_s = sanitize_weights(sig_arrs[args.weight_branch], args.abs_weights)
        w_b = sanitize_weights(bkg_arrs[args.weight_branch], args.abs_weights)

        s_s = np.ravel(compute_scores_model(model, xs_s))
        s_b = np.ravel(compute_scores_model(model, xs_b))
        logger.info("Scan inputs: S entries=%d (sum=%.3f), B entries=%d (sum=%.3f)",
                s_s.size, float(w_s.sum()), s_b.size, float(w_b.sum()))

        cuts = parse_cuts(args.scan_cuts)
        rows = []
        best_z = {"val": -1.0, "cut": 0.0}
        for c in cuts:
            s = float(w_s[s_s >= c].sum())
            b = float(w_b[s_b >= c].sum())
            z = s / np.sqrt(s + b + 1e-12) if (s + b) > 0 else 0.0
            rows.append([f"{c:.4f}", f"{s:.4f}", f"{b:.4f}", f"{z:.4f}"])
            if z > best_z["val"]:
                best_z = {"val": z, "cut": float(c)}

        logger.info("Scan best Z=%.3f at cut=%.4f",
                    best_z["val"], best_z["cut"])
        logger.info("Scan table (coarse sampling):")
        step = max(1, int(len(rows) / 10))
        for r in rows[::step]:
            logger.info("  cut=%s | S=%s | B=%s | Z=%s", r[0], r[1], r[2], r[3])

        plot_scan_curves(
            pdf,
            f"Model-based scan (mass {args.mass_min:.0f}-{args.mass_max:.0f} GeV)",
            cuts,
            [float(r[3]) for r in rows],
            [float(r[3]) for r in rows],
            best_z,
            best_z,
        )
        plot_table(
            pdf,
            "Model-based scan summary",
            ["Cut", "S", "B", "Z"],
            rows[::max(1, int(len(rows) / 20))],
            col_align=["right", "right", "right", "right"],
        )
    elif args.scan_sig or args.scan_bkg:
        logger.warning("Scan skipped: provide both --scan-sig and --scan-bkg")

    logger.info("Wrote %s", args.output)
    print(f"[done] Wrote {args.output}")


if __name__ == "__main__":
    main()
