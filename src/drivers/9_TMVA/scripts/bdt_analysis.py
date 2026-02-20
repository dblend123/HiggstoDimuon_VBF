#!/usr/bin/env python3

import argparse
import glob
import os
import sys
import logging
import itertools

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
    p.add_argument("--shape-mu-max", type=float, default=3.0,
                   help="Upper bound of mu grid in profiled likelihood scan")
    p.add_argument("--shape-mu-points", type=int, default=121,
                   help="Number of mu scan points for profile likelihood")

    # Optional Data/MC summary (same format as tmva_DataMC_BDT.py)
    p.add_argument("--data-era", action="append", default=[],
                   help="ERA:/path/data1.root[,/path/data2.root]")
    p.add_argument("--mc-era", action="append", default=[],
                   help="ERA:PROC:/path/mc1.root[,/path/mc2.root]  (repeat)")
    p.add_argument("--datamc-nbins", type=int, default=40)
    p.add_argument("--datamc-logy", action="store_true")
    p.add_argument("--lumi", type=float, default=None,
                   help="Luminosity in fb^-1; if set, MC is rescaled to match data counts. For 2024, the value is 104.80585 fb^-1.")
    p.add_argument("--ams-cut-offsets", default="0.02,0.05",
                   help="Comma-separated offsets below AMS cut for extra data mass dumps (e.g. 0.01,0.02).")
    p.add_argument("--mass-log-max", type=int, default=0,
                   help="Max number of data masses to log; 0 means log all.")

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


def parse_bins(spec):
    return [float(x) for x in spec.split(",")]


def parse_group_map(items, flagname):
    out = {}
    for s in items:
        if ":" not in s:
            raise SystemExit(f"[error] {flagname} expects key:val, got '{s}'")
        k, v = s.split(":", 1)
        out.setdefault(k, []).append(v)
    return out


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


def compute_scores_model(model, xs):
    X = np.vstack(xs).T
    return model.Compute(X.astype(np.float32, copy=False)).astype(np.float64, copy=False)


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
    za_bins = np.asarray([ams(float(s), float(b)) for s, b in zip(s_bins, b_bins)], dtype=np.float64)
    ztot = multi_bin_asimov_z(s_bins, b_bins)

    return {
        "edges": [float(x) for x in edges],
        "s_bins": s_bins,
        "b_bins": b_bins,
        "sb_bins": ratios,
        "zapprox_bins": zapprox_bins,
        "za_bins": za_bins,
        "za_total": float(ztot),
        "n_bins": int(len(edges) - 1),
    }


def optimize_multibin_boundaries(sig_scores, sig_w, bkg_scores, bkg_w,
                                 min_n=4, max_n=6, n_candidates=16,
                                 min_sig=1.0, min_bkg=5.0, enforce_monotonic=True):
    min_n = max(2, int(min_n))
    max_n = max(min_n, int(max_n))
    cands = candidate_boundaries(sig_scores, sig_w, bkg_scores, bkg_w, n_candidates)

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
            if best_local is None or res["za_total"] > best_local["za_total"]:
                best_local = res

        scan_stats[n_bins] = {"scanned": int(scanned), "valid": int(valid)}
        if best_local is not None:
            best_by_n[n_bins] = best_local
            if best_global is None or best_local["za_total"] > best_global["za_total"]:
                best_global = best_local

    return {
        "best": best_global,
        "best_by_n": best_by_n,
        "scan_stats": scan_stats,
        "candidates": cands,
    }


def orthogonal_poly_terms(x, order):
    terms = []
    if order >= 1:
        terms.append(np.asarray(x, dtype=np.float64))
    if order >= 2:
        terms.append(0.5 * (3.0 * np.asarray(x, dtype=np.float64) ** 2 - 1.0))
    return terms


def build_shape_inputs(sig_mass, sig_score, sig_w,
                       bkg_mass, bkg_score, bkg_w,
                       bdt_edges, mass_min, mass_max, mass_nbins):
    mass_edges = np.linspace(float(mass_min), float(mass_max), int(mass_nbins) + 1)
    mass_centers = 0.5 * (mass_edges[:-1] + mass_edges[1:])
    x = 2.0 * (mass_centers - float(mass_min)) / max(float(mass_max - mass_min), 1e-12) - 1.0

    n_bdt = len(bdt_edges) - 1
    n_mass = len(mass_edges) - 1
    s_hist = np.zeros((n_bdt, n_mass), dtype=np.float64)
    b_hist = np.zeros((n_bdt, n_mass), dtype=np.float64)

    for i, (lo, hi) in enumerate(zip(bdt_edges[:-1], bdt_edges[1:])):
        if i == n_bdt - 1:
            sm = (sig_score >= lo) & (sig_score <= hi)
            bm = (bkg_score >= lo) & (bkg_score <= hi)
        else:
            sm = (sig_score >= lo) & (sig_score < hi)
            bm = (bkg_score >= lo) & (bkg_score < hi)

        s_hist[i], _ = np.histogram(sig_mass[sm], bins=mass_edges, weights=sig_w[sm])
        b_hist[i], _ = np.histogram(bkg_mass[bm], bins=mass_edges, weights=bkg_w[bm])

    s_y = s_hist.sum(axis=1)
    b_y = b_hist.sum(axis=1)
    s_shape = np.divide(s_hist, np.maximum(s_y[:, None], 1e-12))

    core = np.sum(b_hist, axis=0)
    core = np.divide(core, max(float(np.sum(core)), 1e-12))

    asimov = s_hist + b_hist
    return {
        "mass_edges": mass_edges,
        "mass_centers": mass_centers,
        "x": x,
        "s_hist": s_hist,
        "b_hist": b_hist,
        "s_y": s_y,
        "b_y": b_y,
        "s_shape": s_shape,
        "core": core,
        "asimov": asimov,
        "n_bdt": n_bdt,
        "n_mass": n_mass,
    }


def unpack_shape_nuisance(vec, n_bdt, poly_order):
    vec = np.asarray(vec, dtype=np.float64)
    n_norm = n_bdt
    n_poly = n_bdt * max(int(poly_order), 0)
    need = n_norm + n_poly
    if vec.size < need:
        out = np.zeros(need, dtype=np.float64)
        out[:vec.size] = vec
        vec = out
    eta = vec[:n_norm]
    zeta = vec[n_norm:n_norm + n_poly].reshape((n_bdt, max(int(poly_order), 0))) if poly_order > 0 else np.zeros((n_bdt, 0), dtype=np.float64)
    return eta, zeta


def build_bkg_shapes_from_nuisance(core_shape, x, zeta_row, poly_order, poly_unc):
    mod = np.ones_like(core_shape, dtype=np.float64)
    terms = orthogonal_poly_terms(x, poly_order)
    for j in range(min(len(terms), zeta_row.size)):
        mod += float(poly_unc) * float(zeta_row[j]) * terms[j]
    mod = np.clip(mod, 1e-6, None)
    shaped = core_shape * mod
    shaped /= max(float(np.sum(shaped)), 1e-12)
    return shaped


def build_expected_counts(mu, nuisance_vec, shape_inputs,
                          poly_order=1, bkg_norm_unc=0.20, poly_unc=0.30):
    n_bdt = shape_inputs["n_bdt"]
    eta, zeta = unpack_shape_nuisance(nuisance_vec, n_bdt, poly_order)
    lam = np.zeros((shape_inputs["n_bdt"], shape_inputs["n_mass"]), dtype=np.float64)

    for i in range(shape_inputs["n_bdt"]):
        s = max(float(mu), 0.0) * float(shape_inputs["s_y"][i]) * shape_inputs["s_shape"][i]
        b_norm = float(shape_inputs["b_y"][i]) * np.exp(float(bkg_norm_unc) * float(eta[i]))
        b_shape = build_bkg_shapes_from_nuisance(
            shape_inputs["core"],
            shape_inputs["x"],
            zeta[i] if i < zeta.shape[0] else np.zeros(0, dtype=np.float64),
            poly_order,
            poly_unc,
        )
        b = b_norm * b_shape
        lam[i] = s + b
    return np.clip(lam, 1e-12, None)


def poisson_nll(data, model):
    d = np.asarray(data, dtype=np.float64)
    m = np.asarray(model, dtype=np.float64)
    return float(np.sum(m - d * np.log(np.clip(m, 1e-12, None))))


def gaussian_penalty(nuisance_vec):
    v = np.asarray(nuisance_vec, dtype=np.float64)
    return float(0.5 * np.sum(v ** 2))


def total_shape_nll(mu, nuisance_vec, shape_inputs,
                    poly_order=1, bkg_norm_unc=0.20, poly_unc=0.30,
                    data_key="asimov"):
    lam = build_expected_counts(
        mu,
        nuisance_vec,
        shape_inputs,
        poly_order=poly_order,
        bkg_norm_unc=bkg_norm_unc,
        poly_unc=poly_unc,
    )
    data = shape_inputs[data_key]
    return poisson_nll(data, lam) + gaussian_penalty(nuisance_vec)


def profile_nll_fixed_mu(mu, shape_inputs, poly_order, bkg_norm_unc, poly_unc,
                         max_rounds=7, init=None):
    n_params = shape_inputs["n_bdt"] * (1 + max(int(poly_order), 0))
    if init is None:
        vec = np.zeros(n_params, dtype=np.float64)
    else:
        vec = np.asarray(init, dtype=np.float64).copy()
        if vec.size != n_params:
            tmp = np.zeros(n_params, dtype=np.float64)
            tmp[:min(n_params, vec.size)] = vec[:min(n_params, vec.size)]
            vec = tmp

    step_schedule = [1.0, 0.6, 0.35, 0.20, 0.12, 0.07, 0.04]
    step_schedule = step_schedule[:max(1, min(max_rounds, len(step_schedule)))]

    best = total_shape_nll(mu, vec, shape_inputs, poly_order, bkg_norm_unc, poly_unc)
    for step in step_schedule:
        improved = False
        for ip in range(n_params):
            cur = vec[ip]
            cands = np.array([cur - step, cur, cur + step], dtype=np.float64)
            cands = np.clip(cands, -5.0, 5.0)
            vals = []
            for c in cands:
                test = vec.copy()
                test[ip] = c
                vals.append(total_shape_nll(mu, test, shape_inputs, poly_order, bkg_norm_unc, poly_unc))
            best_idx = int(np.argmin(vals))
            if vals[best_idx] + 1e-10 < best:
                vec[ip] = float(cands[best_idx])
                best = float(vals[best_idx])
                improved = True
        if not improved:
            continue

    return best, vec


def compute_shape_expected_significance(shape_inputs,
                                        poly_order=1,
                                        bkg_norm_unc=0.20,
                                        poly_unc=0.30,
                                        mu_max=3.0,
                                        mu_points=121):
    mu_points = max(int(mu_points), 11)
    mu_grid = np.linspace(0.0, float(mu_max), mu_points)

    nll_mu0, vec_mu0 = profile_nll_fixed_mu(
        0.0, shape_inputs, poly_order, bkg_norm_unc, poly_unc, init=None
    )

    prof = []
    vec_init = vec_mu0
    best = None
    for mu in mu_grid:
        nll_mu, vec_init = profile_nll_fixed_mu(
            float(mu),
            shape_inputs,
            poly_order,
            bkg_norm_unc,
            poly_unc,
            init=vec_init,
        )
        prof.append((float(mu), float(nll_mu)))
        if best is None or nll_mu < best[1]:
            best = (float(mu), float(nll_mu), vec_init.copy())

    mu_hat, nll_best, vec_best = best
    q0 = max(0.0, 2.0 * (nll_mu0 - nll_best))
    zexp = float(np.sqrt(q0))

    model_b = build_expected_counts(0.0, vec_mu0, shape_inputs, poly_order, bkg_norm_unc, poly_unc)
    model_sb = build_expected_counts(mu_hat, vec_best, shape_inputs, poly_order, bkg_norm_unc, poly_unc)

    return {
        "zexp": zexp,
        "q0": q0,
        "mu_hat": float(mu_hat),
        "nll_mu0": float(nll_mu0),
        "nll_best": float(nll_best),
        "vec_mu0": vec_mu0,
        "vec_best": vec_best,
        "profile": prof,
        "model_b": model_b,
        "model_sb": model_sb,
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


def plot_text_page(pdf, title, lines):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle(title, fontsize=14, y=0.97)
    text = "\n".join(lines)
    fig.text(0.07, 0.9, text, fontsize=11, va="top")
    pdf.savefig(fig)
    plt.close(fig)


def plot_highcut_mass(pdf, title, data_mass, mc_groups, bins, logy):
    fig, ax = plt.subplots(figsize=(11, 7))

    mc_hists = []
    mc_labels = []
    mc_colors = []

    for name, masses, weights in mc_groups:
        hist, _ = np.histogram(masses, bins=bins, weights=weights)
        mc_hists.append(hist)
        mc_labels.append(canonical_proc(name))
        mc_colors.append(proc_color(name))

    if mc_hists:
        ax.hist(
            [bins[:-1]] * len(mc_hists),
            bins=bins,
            weights=mc_hists,
            stacked=True,
            color=mc_colors,
            label=mc_labels,
            histtype="stepfilled",
            alpha=0.85,
        )

    if data_mass.size:
        data_hist, _ = np.histogram(data_mass, bins=bins)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        ax.errorbar(
            bin_centers,
            data_hist,
            yerr=np.sqrt(np.maximum(data_hist, 1.0)),
            fmt="o",
            color="black",
            label="Data",
        )

    ax.set_xlabel("Mass [GeV]")
    ax.set_ylabel("Events / bin")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=8, ncol=2)

    pdf.savefig(fig)
    plt.close(fig)


def plot_scan_curves(pdf, title, cuts, z_vals, za_vals, best_z, best_za):
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.plot(cuts, z_vals, linewidth=2, label="Z = S/sqrt(S+B)")
    ax.plot(cuts, za_vals, linewidth=2, label="ZA (asymptotic)")
    ax.axvline(best_z["cut"], color="tab:red", linestyle="--", linewidth=1.5)
    ax.axvline(best_za["cut"], color="tab:green", linestyle=":", linewidth=1.5)
    ax.set_xlabel("BDT score cut")
    ax.set_ylabel("Significance")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    ax.text(
        0.02, 0.95,
        "Best cuts\n"
        f"Z:  cut={best_z['cut']:.3f}, Z={best_z['val']:.3f}\n"
        f"ZA: cut={best_za['cut']:.3f}, ZA={best_za['val']:.3f}",
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
        "Multi-bin BDT scan: bins=%d-%d | candidates=%d | min(S,B)=(%.3f, %.3f) | monotonic=%s",
        args.multibin_min_n,
        args.multibin_max_n,
        args.multibin_candidates,
        args.multibin_min_sig,
        args.multibin_min_bkg,
        str(not args.multibin_no_monotonic),
    )
    logger.info(
        "Shape likelihood: enabled=%s | mass bins=%d | poly order=%d | bkg norm unc=%.3f | poly unc=%.3f",
        str(args.shape_likelihood),
        args.shape_mass_nbins,
        args.shape_poly_order,
        args.shape_bkg_norm_unc,
        args.shape_poly_unc,
    )
    outputs_dir = os.path.abspath(args.outputs_dir)
    eras = discover_eras(outputs_dir, args.eras)
    if not eras:
        raise SystemExit("No eras discovered. Check --outputs-dir or --eras.")
    logger.info("Discovered eras: %s", ", ".join(eras))

    bdt_bins = parse_bins(args.bdt_bins)
    logger.info("Base BDT bins: %s", ", ".join([f"{b:.2f}" for b in bdt_bins]))

    data_map_raw = parse_group_map(args.data_era, "--data-era") if args.data_era else {}
    mc_map_raw = parse_group_map(args.mc_era, "--mc-era") if args.mc_era else {}
    if data_map_raw or mc_map_raw:
        logger.info("Data/MC inputs provided for eras: %s",
                    ", ".join(sorted(set(list(data_map_raw.keys()) + list(mc_map_raw.keys())))))
        logger.info("Data/MC bins: %d | logy=%s", args.datamc_nbins, str(args.datamc_logy))

    with PdfPages(args.output) as pdf:
        best_cut_by_era = {}
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
                    "%s: S events=%d (sum=%.3f), B events=%d (sum=%.3f)",
                    era, sig_scores.size, float(sig_w.sum()), bkg_scores.size, float(bkg_w.sum())
                )

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
                best_cut_by_era[era] = best_summary["thr"]
                plot_ams(pdf, f"{era}: AMS vs signal efficiency", effs, ams_vals, best_summary)
                logger.info(
                    "%s: Best AMS=%.3f at eff=%.2f, cut=%.3f, S=%.3f, B=%.3f, S/B=%s",
                    era,
                    best_summary["ams"],
                    best_summary["eff"],
                    best_summary["thr"],
                    best_summary["s"],
                    best_summary["b"],
                    format_ratio(best_summary["s_over_b"]),
                )

                # Multi-bin BDT categorization (no hard cut)
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
                )
                best_multi = scan["best"]
                multibin_results_by_era[era] = {
                    "best": best_multi,
                    "best_by_n": scan["best_by_n"],
                    "scan_stats": scan["scan_stats"],
                }

                if best_multi is not None:
                    edges_txt = ", ".join([f"{x:.3f}" for x in best_multi["edges"]])
                    logger.info(
                        "%s: Multi-bin optimum: N=%d, ZA_total=%.3f, edges=[%s]",
                        era,
                        best_multi["n_bins"],
                        best_multi["za_total"],
                        edges_txt,
                    )
                    for n_bins in range(args.multibin_min_n, args.multibin_max_n + 1):
                        st = scan["scan_stats"].get(n_bins, {"scanned": 0, "valid": 0})
                        if n_bins in scan["best_by_n"]:
                            logger.info(
                                "%s: N=%d scanned=%d valid=%d best_ZA=%.3f",
                                era,
                                n_bins,
                                st["scanned"],
                                st["valid"],
                                scan["best_by_n"][n_bins]["za_total"],
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
                            f"{best_multi['za_bins'][i]:.3f}",
                        ])
                    plot_table(
                        pdf,
                        f"{era}: Multi-bin BDT categorization (no hard cut)\n"
                        f"Optimized N={best_multi['n_bins']}, ZA(total)={best_multi['za_total']:.3f}, "
                        f"monotonic S/B={str(not args.multibin_no_monotonic)}",
                        ["Bin", "BDT range", "S", "B", "S/B", "S/sqrt(B)", "ZA(bin)"],
                        rows,
                        col_align=["right", "center", "right", "right", "right", "right", "right"],
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
                            )
                            shape_fit = compute_shape_expected_significance(
                                shape_inputs,
                                poly_order=poly_order,
                                bkg_norm_unc=float(args.shape_bkg_norm_unc),
                                poly_unc=float(args.shape_poly_unc),
                                mu_max=float(args.shape_mu_max),
                                mu_points=int(args.shape_mu_points),
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
                else:
                    logger.warning(
                        "%s: Multi-bin optimization found no valid configuration (check min yields / candidates / monotonicity)",
                        era,
                    )

                suggested_cuts = pick_top_cuts(thr_vals, ams_vals, max_cuts=3, min_sep=0.02)
                if suggested_cuts:
                    logger.info(
                        "%s: Suggested high-score cuts: %s",
                        era,
                        ", ".join([f"{c:.3f}" for c in suggested_cuts]),
                    )

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
                logger.info("%s: Yield table bins=%d (refined near %.3f)", era, len(rows), best_summary["thr"])
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
                    corr_sig = weighted_corr(sig_mass, sig_scores, sig_w)
                    corr_bkg = weighted_corr(bkg_mass, bkg_scores, bkg_w)
                    plot_mass_sculpt_norm(
                        pdf,
                        f"{era}: Signal mass sculpting (normalized shapes)",
                        sig_mass,
                        sig_scores,
                        sig_w,
                        bdt_bins,
                        corr_sig,
                    )
                    plot_mass_sculpt_norm(
                        pdf,
                        f"{era}: Background mass sculpting (normalized shapes)",
                        bkg_mass,
                        bkg_scores,
                        bkg_w,
                        bdt_bins,
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
                            [args.bdt_branch, dm_mass_branch],
                            dm_sel,
                        )
                        data_scores = data_arr[args.bdt_branch]
                        data_mass = data_arr[dm_mass_branch]
                        
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
                        best_cut = best_cut_by_era.get(era)
                        if best_cut is not None:
                            offsets = []
                            for item in args.ams_cut_offsets.split(","):
                                item = item.strip()
                                if not item:
                                    continue
                                try:
                                    offsets.append(float(item))
                                except ValueError:
                                    logger.warning("%s: Invalid --ams-cut-offsets entry '%s'", era, item)
                            cuts_to_show = []
                            seen = set()
                            for off in [0.0] + [o for o in offsets if o > 0.0]:
                                cut = min(1.0, max(0.0, best_cut - off))
                                key = round(cut, 6)
                                if key in seen:
                                    continue
                                seen.add(key)
                                cuts_to_show.append(cut)

                            mass_bins = np.linspace(args.mass_min, args.mass_max, 21)
                            for cut in cuts_to_show:
                                data_high = data_scores >= cut
                                data_mass_high = data_mass[data_high]
                                data_mass_sorted = np.sort(data_mass_high)
                                mc_high_groups = []
                                mc_high_total = 0.0
                                for proc, mc_scores, mc_mass, mc_weights in mc_full_groups:
                                    mc_mask = mc_scores >= cut
                                    mc_mass_sel = mc_mass[mc_mask]
                                    mc_weight_sel = mc_weights[mc_mask]
                                    mc_high_groups.append((proc, mc_mass_sel, mc_weight_sel))
                                    mc_high_total += float(mc_weight_sel.sum())

                                delta = best_cut - cut
                                if abs(delta) < 1e-6:
                                    cut_label = "AMS suggested cut"
                                else:
                                    cut_label = f"AMS - {delta:.3f}"

                                logger.info(
                                    "%s: High-score data/MC (BDT>=%.3f): data_n=%d, mc_sum=%.3f, ratio=%.3f",
                                    era,
                                    cut,
                                    data_mass_high.size,
                                    mc_high_total,
                                    data_mass_high.size / max(mc_high_total, 1e-12),
                                )
                                if data_mass_sorted.size:
                                    if args.mass_log_max > 0:
                                        max_log = min(args.mass_log_max, data_mass_sorted.size)
                                    else:
                                        max_log = data_mass_sorted.size
                                    mass_log = ", ".join([f"{m:.2f}" for m in data_mass_sorted[:max_log]])
                                    if data_mass_sorted.size > max_log:
                                        mass_log += ", ..."
                                    logger.info(
                                        "%s: Data masses (BDT>=%.3f) [GeV]: %s",
                                        era,
                                        cut,
                                        mass_log,
                                    )

                                mass_lines = [
                                    f"BDT >= {cut:.3f} ({cut_label})",
                                    f"Data events: {data_mass_high.size}",
                                    f"MC sum(weights): {mc_high_total:.3f}",
                                    "Data masses [GeV]:",
                                ]
                                if data_mass_sorted.size:
                                    max_show = 20
                                    shown = ", ".join([f"{m:.2f}" for m in data_mass_sorted[:max_show]])
                                    mass_lines.append(shown)
                                    if data_mass_sorted.size > max_show:
                                        mass_lines.append("...")
                                else:
                                    mass_lines.append("(none)")

                                plot_text_page(
                                    pdf,
                                    f"{era}: Data mass for BDT >= {cut:.3f}",
                                    [line for line in mass_lines if line],
                                )

                                plot_highcut_mass(
                                    pdf,
                                    f"{era}: Mass for BDT >= {cut:.3f}",
                                    data_mass_high,
                                    mc_high_groups,
                                    mass_bins,
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

        # Cross-era boundary stability summary (same N bins)
        if len(multibin_results_by_era) >= 2:
            rows = []
            for n_bins in range(args.multibin_min_n, args.multibin_max_n + 1):
                edges_list = []
                z_vals = []
                used_eras = []
                for era in eras:
                    era_res = multibin_results_by_era.get(era, {})
                    by_n = era_res.get("best_by_n", {}) if era_res else {}
                    if n_bins in by_n:
                        used_eras.append(era)
                        z_vals.append(float(by_n[n_bins]["za_total"]))
                        edges_list.append(by_n[n_bins]["edges"])
                if len(edges_list) < 2:
                    continue

                internal = np.asarray([e[1:-1] for e in edges_list], dtype=np.float64)
                stds = np.std(internal, axis=0)
                max_std = float(np.max(stds)) if stds.size else np.nan
                mean_std = float(np.mean(stds)) if stds.size else np.nan
                rows.append([
                    str(n_bins),
                    str(len(used_eras)),
                    f"{np.mean(z_vals):.3f}",
                    f"{np.min(z_vals):.3f}",
                    f"{np.max(z_vals):.3f}",
                    f"{mean_std:.4f}",
                    f"{max_std:.4f}",
                ])

            if rows:
                plot_table(
                    pdf,
                    "Cross-era stability of optimized BDT category boundaries\n"
                    "(smaller boundary std = more stable across eras)",
                    ["N bins", "#eras", "<ZA>", "ZA min", "ZA max", "mean std(edge)", "max std(edge)"],
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
        best_za = {"val": -1.0, "cut": 0.0}
        for c in cuts:
            s = float(w_s[s_s >= c].sum())
            b = float(w_b[s_b >= c].sum())
            z = s / np.sqrt(s + b + 1e-12) if (s + b) > 0 else 0.0
            za = ams(s, b)
            rows.append([f"{c:.4f}", f"{s:.4f}", f"{b:.4f}", f"{z:.4f}", f"{za:.4f}"])
            if z > best_z["val"]:
                best_z = {"val": z, "cut": float(c)}
            if za > best_za["val"]:
                best_za = {"val": za, "cut": float(c)}

        logger.info("Scan best Z=%.3f at cut=%.4f; best ZA=%.3f at cut=%.4f",
                    best_z["val"], best_z["cut"], best_za["val"], best_za["cut"])
        logger.info("Scan table (coarse sampling):")
        step = max(1, int(len(rows) / 10))
        for r in rows[::step]:
            logger.info("  cut=%s | S=%s | B=%s | Z=%s | ZA=%s", r[0], r[1], r[2], r[3], r[4])

        plot_scan_curves(
            pdf,
            f"Model-based scan (mass {args.mass_min:.0f}-{args.mass_max:.0f} GeV)",
            cuts,
            [float(r[3]) for r in rows],
            [float(r[4]) for r in rows],
            best_z,
            best_za,
        )
        plot_table(
            pdf,
            "Model-based scan summary",
            ["Cut", "S", "B", "Z", "ZA"],
            rows[::max(1, int(len(rows) / 20))],
            col_align=["right", "right", "right", "right", "right"],
        )
    elif args.scan_sig or args.scan_bkg:
        logger.warning("Scan skipped: provide both --scan-sig and --scan-bkg")

    logger.info("Wrote %s", args.output)
    print(f"[done] Wrote {args.output}")


if __name__ == "__main__":
    main()
