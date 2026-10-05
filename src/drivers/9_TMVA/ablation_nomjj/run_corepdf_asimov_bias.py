#!/usr/bin/env python3
"""Screen DY-inspired core PDFs with background-only Asimov spurious-signal tests."""

import argparse
import csv
import math
import os

import numpy as np
import ROOT
from scipy.optimize import brentq, minimize


ROOT.gROOT.SetBatch(True)

MZ = 91.1876
GAMMA_Z = 2.4952


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", default="2024")
    parser.add_argument(
        "--data",
        default="../Outputs_nomjj_newsample/2024/2024_data_scored.root",
    )
    parser.add_argument(
        "--signal",
        default="../Outputs_nomjj_newsample/2024/2024_scored_signal.root",
    )
    parser.add_argument("--treename", default="mva")
    parser.add_argument(
        "--output-dir",
        default="../Outputs_nomjj_newsample/CorePdfBias/2024",
    )
    parser.add_argument("--category-edges", default="0.60,0.85,0.95,1.00")
    parser.add_argument("--mass-min", type=float, default=110.0)
    parser.add_argument("--mass-max", type=float, default=150.0)
    parser.add_argument("--mass-bin-width", type=float, default=0.5)
    parser.add_argument("--sb1-min", type=float, default=110.0)
    parser.add_argument("--sb1-max", type=float, default=115.0)
    parser.add_argument("--sb2-min", type=float, default=135.0)
    parser.add_argument("--sb2-max", type=float, default=150.0)
    parser.add_argument(
        "--core-orders",
        default="0,1,2,3",
        help="Comma-separated Bernstein correction orders to screen",
    )
    parser.add_argument(
        "--max-bias-over-sigma",
        type=float,
        default=0.20,
        help="Provisional screening threshold; replace with group-approved value",
    )
    return parser.parse_args()


def parse_edges(text):
    edges = np.asarray([float(item) for item in text.split(",")], dtype=float)
    if edges.size < 2 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("Category edges must be strictly increasing")
    return edges


def load_arrays(path, tree_name, branches):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    tree = root_file.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Missing tree '{tree_name}' in {path}")
    missing = [branch for branch in branches if not tree.GetBranch(branch)]
    root_file.Close()
    if missing:
        raise RuntimeError(f"Missing branches in {path}: {missing}")
    arrays = ROOT.RDataFrame(tree_name, path).AsNumpy(branches)
    return {key: np.asarray(value, dtype=np.float64) for key, value in arrays.items()}


def normalized_shape(density, widths):
    density = np.asarray(density, dtype=float)
    weighted = np.maximum(density, 1.0e-300) * widths
    total = float(np.sum(weighted))
    if not math.isfinite(total) or total <= 0.0:
        return None
    return weighted / total


def bernstein_value(x, degree, coefficients):
    result = np.zeros_like(x, dtype=float)
    for index, coefficient in enumerate(coefficients):
        result += (
            coefficient
            * math.comb(degree, index)
            * np.power(x, index)
            * np.power(1.0 - x, degree - index)
        )
    return result


def truth_shape(name, parameters, centers, widths, mass_min, mass_max):
    x = (centers - 130.0) / 20.0
    unit_x = (centers - mass_min) / (mass_max - mass_min)
    if name == "exponential":
        density = np.exp(np.clip(parameters[0] * x, -50.0, 50.0))
    elif name == "powerlaw":
        density = np.power(centers / 125.0, -parameters[0])
    elif name == "exp_poly2":
        density = np.exp(
            np.clip(parameters[0] * x + parameters[1] * x * x, -50.0, 50.0)
        )
    elif name.startswith("bernstein"):
        degree = int(name.replace("bernstein", ""))
        coefficients = np.concatenate(([1.0], np.exp(np.asarray(parameters))))
        density = bernstein_value(unit_x, degree, coefficients)
    else:
        raise ValueError(name)
    return normalized_shape(density, widths)


def core_shape(order, parameters, centers, widths, mass_min, mass_max):
    """Modified Z-tail core times a positive Bernstein correction."""
    scaled = (centers - 125.0) / 25.0
    denominator = np.square(centers - MZ) + 0.25 * GAMMA_Z * GAMMA_Z
    core = np.exp(np.clip(parameters[0] * scaled, -50.0, 50.0)) / denominator
    if order > 0:
        unit_x = (centers - mass_min) / (mass_max - mass_min)
        coefficients = np.concatenate(([1.0], np.exp(np.asarray(parameters[1:]))))
        core *= bernstein_value(unit_x, order, coefficients)
    return normalized_shape(core, widths)


def poisson_nll(observed, expected):
    if np.any(~np.isfinite(expected)) or np.any(expected <= 0.0):
        return 1.0e100
    return float(np.sum(expected - observed * np.log(expected)))


def truth_configuration(name):
    if name in ("exponential", "powerlaw"):
        return np.zeros(1), [(-20.0, 20.0)]
    if name == "exp_poly2":
        return np.zeros(2), [(-20.0, 20.0), (-20.0, 20.0)]
    degree = int(name.replace("bernstein", ""))
    return np.zeros(degree), [(-8.0, 8.0)] * degree


def fit_truth_sidebands(name, observed, sb_mask, centers, widths, args):
    shape_start, shape_bounds = truth_configuration(name)
    initial = np.concatenate(([math.log(max(1.0, np.sum(observed) * 2.0))], shape_start))
    bounds = [(math.log(1.0e-3), math.log(1.0e8))] + shape_bounds

    def objective(vector):
        shape = truth_shape(
            name, vector[1:], centers, widths, args.mass_min, args.mass_max
        )
        if shape is None:
            return 1.0e100
        expected = math.exp(vector[0]) * shape[sb_mask]
        return poisson_nll(observed[sb_mask], expected)

    result = minimize(objective, initial, method="L-BFGS-B", bounds=bounds)
    shape = truth_shape(
        name, result.x[1:], centers, widths, args.mass_min, args.mass_max
    )
    return result, math.exp(result.x[0]) * shape


def fit_core_plus_signal(order, observed, signal_counts, centers, widths, args):
    parameter_count = 1 + order
    initial = np.concatenate(
        ([math.log(max(1.0, np.sum(observed)))], np.zeros(parameter_count), [0.0])
    )
    bounds = (
        [(math.log(1.0e-3), math.log(1.0e8))]
        + [(-20.0, 20.0)]
        + [(-8.0, 8.0)] * order
        + [(-20.0, 20.0)]
    )

    def objective(vector):
        shape = core_shape(
            order, vector[1:-1], centers, widths, args.mass_min, args.mass_max
        )
        if shape is None:
            return 1.0e100
        expected = math.exp(vector[0]) * shape + vector[-1] * signal_counts
        return poisson_nll(observed, expected)

    objective.parameter_bounds = bounds

    results = []
    for initial_mu in (-3.0, 0.0, 3.0):
        trial = initial.copy()
        trial[-1] = initial_mu
        results.append(
            minimize(
                objective,
                trial,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 5000, "ftol": 1.0e-12, "gtol": 1.0e-8},
            )
        )
    successful = [result for result in results if result.success]
    result = min(successful or results, key=lambda item: item.fun)
    return result, objective


def numerical_hessian(function, point):
    point = np.asarray(point, dtype=float)
    size = len(point)
    steps = 2.0e-4 * np.maximum(1.0, np.abs(point))
    hessian = np.zeros((size, size), dtype=float)
    f0 = function(point)
    for i in range(size):
        plus = point.copy()
        minus = point.copy()
        plus[i] += steps[i]
        minus[i] -= steps[i]
        hessian[i, i] = (function(plus) - 2.0 * f0 + function(minus)) / steps[i] ** 2
        for j in range(i + 1, size):
            pp = point.copy()
            pm = point.copy()
            mp = point.copy()
            mm = point.copy()
            pp[i] += steps[i]
            pp[j] += steps[j]
            pm[i] += steps[i]
            pm[j] -= steps[j]
            mp[i] -= steps[i]
            mp[j] += steps[j]
            mm[i] -= steps[i]
            mm[j] -= steps[j]
            value = (function(pp) - function(pm) - function(mp) + function(mm)) / (
                4.0 * steps[i] * steps[j]
            )
            hessian[i, j] = value
            hessian[j, i] = value
    return hessian


def estimate_mu_error(objective, result):
    """Return the average 68% profile-likelihood uncertainty on mu.

    A numerical-Hessian estimate is retained only as a fallback.  Profiling is
    more reliable when flexible background coefficients are strongly correlated.
    """
    bounds = getattr(objective, "parameter_bounds", None)
    if bounds is not None and result.success:
        best = np.asarray(result.x, dtype=float)
        minimum = float(result.fun)
        nuisance_bounds = bounds[:-1]
        mu_low, mu_high = bounds[-1]

        def profile_delta(mu_value):
            def nuisance_objective(nuisance):
                return objective(np.concatenate((nuisance, [mu_value])))

            profiled = minimize(
                nuisance_objective,
                best[:-1],
                method="L-BFGS-B",
                bounds=nuisance_bounds,
                options={"maxiter": 5000, "ftol": 1.0e-12, "gtol": 1.0e-8},
            )
            return float(profiled.fun - minimum - 0.5)

        def find_crossing(direction):
            origin = float(best[-1])
            previous = origin
            previous_value = -0.5
            step = 0.25
            for _ in range(20):
                candidate = origin + direction * step
                candidate = min(mu_high, max(mu_low, candidate))
                value = profile_delta(candidate)
                if value >= 0.0:
                    left, right = sorted((previous, candidate))
                    return float(brentq(profile_delta, left, right, xtol=1.0e-5))
                if candidate in (mu_low, mu_high):
                    break
                previous = candidate
                previous_value = value
                step *= 1.7
            return float("nan")

        lower = find_crossing(-1.0)
        upper = find_crossing(+1.0)
        if math.isfinite(lower) and math.isfinite(upper) and upper > lower:
            return 0.5 * (upper - lower)

    # Fallback for an unexpected profile failure.
    try:
        hessian = numerical_hessian(objective, result.x)
        covariance = np.linalg.pinv(hessian, rcond=1.0e-10)
        variance = float(covariance[-1, -1])
        if not math.isfinite(variance) or variance <= 0.0:
            return float("nan")
        return math.sqrt(variance)
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        return float("nan")


def category_mask(scores, low, high, is_last):
    return (scores >= low) & ((scores <= high) if is_last else (scores < high))


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    category_edges = parse_edges(args.category_edges)
    core_orders = tuple(int(item.strip()) for item in args.core_orders.split(","))
    if not core_orders or any(order < 0 for order in core_orders):
        raise ValueError("Core orders must be non-negative integers")

    data = load_arrays(args.data, args.treename, ["m_mumu", "bdt_score"])
    signal = load_arrays(
        args.signal, args.treename, ["m_mumu", "bdt_score", "weight"]
    )

    number_bins = int(round((args.mass_max - args.mass_min) / args.mass_bin_width))
    bin_edges = np.linspace(args.mass_min, args.mass_max, number_bins + 1)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)
    sb_mask = ((centers >= args.sb1_min) & (centers < args.sb1_max)) | (
        (centers >= args.sb2_min) & (centers < args.sb2_max)
    )

    truth_models = (
        "exponential", "powerlaw", "exp_poly2", "bernstein2", "bernstein3", "bernstein4"
    )
    rows = []

    for category_index, (low, high) in enumerate(
        zip(category_edges[:-1], category_edges[1:])
    ):
        category = f"cat{category_index}"
        is_last = category_index == len(category_edges) - 2
        data_selection = category_mask(data["bdt_score"], low, high, is_last)
        data_selection &= (
            ((data["m_mumu"] >= args.sb1_min) & (data["m_mumu"] < args.sb1_max))
            | ((data["m_mumu"] >= args.sb2_min) & (data["m_mumu"] < args.sb2_max))
        )
        observed_sb, _ = np.histogram(data["m_mumu"][data_selection], bins=bin_edges)

        signal_selection = category_mask(signal["bdt_score"], low, high, is_last)
        signal_selection &= (signal["m_mumu"] >= args.mass_min) & (
            signal["m_mumu"] < args.mass_max
        )
        signal_counts, _ = np.histogram(
            signal["m_mumu"][signal_selection],
            bins=bin_edges,
            weights=signal["weight"][signal_selection],
        )
        signal_counts = np.maximum(signal_counts.astype(float), 0.0)
        signal_yield = float(np.sum(signal_counts))
        if signal_yield <= 0.0:
            raise RuntimeError(f"Non-positive signal yield in {category}")

        print(
            f"[{category}] sideband data={int(np.sum(observed_sb))} "
            f"signal(110-150)={signal_yield:.6g}",
            flush=True,
        )

        for truth_name in truth_models:
            truth_result, asimov = fit_truth_sidebands(
                truth_name, observed_sb, sb_mask, centers, widths, args
            )
            for order in core_orders:
                fit_result, objective = fit_core_plus_signal(
                    order, asimov, signal_counts, centers, widths, args
                )
                mu_hat = float(fit_result.x[-1])
                mu_error = estimate_mu_error(objective, fit_result)
                ratio = (
                    abs(mu_hat) / mu_error
                    if math.isfinite(mu_error) and mu_error > 0.0
                    else float("nan")
                )
                passed = bool(
                    truth_result.success
                    and fit_result.success
                    and math.isfinite(ratio)
                    and ratio < args.max_bias_over_sigma
                )
                row = {
                    "category": category,
                    "truth_model": truth_name,
                    "core_order": order,
                    "truth_fit_success": int(truth_result.success),
                    "signal_fit_success": int(fit_result.success),
                    "sideband_events": int(np.sum(observed_sb)),
                    "asimov_full_background": float(np.sum(asimov)),
                    "signal_yield": signal_yield,
                    "mu_hat": mu_hat,
                    "mu_error": mu_error,
                    "spurious_signal_events": mu_hat * signal_yield,
                    "bias_over_sigma": ratio,
                    "threshold": args.max_bias_over_sigma,
                    "passes_screen": int(passed),
                }
                rows.append(row)
                print(
                    f"  truth={truth_name:<12} coreOrder={order} "
                    f"mu={mu_hat:+.4f} +/- {mu_error:.4f} "
                    f"|bias|/sigma={ratio:.3f} pass={int(passed)}",
                    flush=True,
                )

    output_path = os.path.join(
        args.output_dir, f"{args.era}_corepdf_asimov_bias.csv"
    )
    fields = list(rows[0].keys())
    with open(output_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary_path = os.path.join(
        args.output_dir, f"{args.era}_corepdf_asimov_summary.csv"
    )
    summary = []
    for category_index in range(len(category_edges) - 1):
        category = f"cat{category_index}"
        for order in core_orders:
            selected = [
                row for row in rows
                if row["category"] == category and row["core_order"] == order
            ]
            finite = [row["bias_over_sigma"] for row in selected if math.isfinite(row["bias_over_sigma"])]
            summary.append(
                {
                    "category": category,
                    "core_order": order,
                    "maximum_bias_over_sigma": max(finite) if finite else float("nan"),
                    "all_truth_models_pass": int(all(row["passes_screen"] for row in selected)),
                    "failed_truth_models": ";".join(
                        row["truth_model"] for row in selected if not row["passes_screen"]
                    ),
                }
            )
    with open(summary_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    print(f"Wrote detailed results: {output_path}")
    print(f"Wrote summary:          {summary_path}")
    print(
        "This is an Asimov screening study. Passing models still require "
        "pseudoexperiment validation with the analysis-approved criterion."
    )


if __name__ == "__main__":
    main()
