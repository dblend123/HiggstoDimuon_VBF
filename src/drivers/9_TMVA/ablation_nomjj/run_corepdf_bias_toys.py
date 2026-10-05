#!/usr/bin/env python3
"""Run background-only toy bias tests for selected core-PDF orders."""

import argparse
import csv
import math
import os

import numpy as np

import run_corepdf_asimov_bias as core


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
        default="../Outputs_nomjj_newsample/CorePdfBiasToys/2024",
    )
    parser.add_argument("--category-edges", default="0.60,0.85,0.95,1.00")
    parser.add_argument(
        "--selected-orders",
        default="cat0:4,cat1:3,cat2:2",
        help="Comma-separated category:order mapping",
    )
    parser.add_argument("--truth-models", default="exponential,powerlaw,exp_poly2,bernstein2,bernstein3,bernstein4")
    parser.add_argument("--toys", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20241002)
    parser.add_argument("--mass-min", type=float, default=110.0)
    parser.add_argument("--mass-max", type=float, default=150.0)
    parser.add_argument("--mass-bin-width", type=float, default=0.5)
    parser.add_argument("--sb1-min", type=float, default=110.0)
    parser.add_argument("--sb1-max", type=float, default=115.0)
    parser.add_argument("--sb2-min", type=float, default=135.0)
    parser.add_argument("--sb2-max", type=float, default=150.0)
    parser.add_argument("--max-bias-over-sigma", type=float, default=0.20)
    parser.add_argument("--max-fit-failure-rate", type=float, default=0.01)
    return parser.parse_args()


def parse_orders(text):
    result = {}
    for item in text.split(","):
        category, order = item.split(":", 1)
        result[category.strip()] = int(order)
    return result


def quantile(values, probability):
    return float(np.quantile(np.asarray(values, dtype=float), probability))


def main():
    args = parse_args()
    if args.toys <= 0:
        raise ValueError("--toys must be positive")
    os.makedirs(args.output_dir, exist_ok=True)

    category_edges = core.parse_edges(args.category_edges)
    selected_orders = parse_orders(args.selected_orders)
    truth_models = tuple(item.strip() for item in args.truth_models.split(",") if item.strip())
    expected_categories = {f"cat{i}" for i in range(len(category_edges) - 1)}
    if set(selected_orders) != expected_categories:
        raise ValueError(
            f"Selected orders must define exactly {sorted(expected_categories)}; "
            f"received {sorted(selected_orders)}"
        )

    data = core.load_arrays(args.data, args.treename, ["m_mumu", "bdt_score"])
    signal = core.load_arrays(
        args.signal, args.treename, ["m_mumu", "bdt_score", "weight"]
    )

    number_bins = int(round((args.mass_max - args.mass_min) / args.mass_bin_width))
    bin_edges = np.linspace(args.mass_min, args.mass_max, number_bins + 1)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)
    sb_mask = ((centers >= args.sb1_min) & (centers < args.sb1_max)) | (
        (centers >= args.sb2_min) & (centers < args.sb2_max)
    )
    rng = np.random.default_rng(args.seed)

    details = []
    summaries = []

    for category_index, (low, high) in enumerate(
        zip(category_edges[:-1], category_edges[1:])
    ):
        category = f"cat{category_index}"
        order = selected_orders[category]
        is_last = category_index == len(category_edges) - 2

        data_selection = core.category_mask(data["bdt_score"], low, high, is_last)
        data_selection &= (
            ((data["m_mumu"] >= args.sb1_min) & (data["m_mumu"] < args.sb1_max))
            | ((data["m_mumu"] >= args.sb2_min) & (data["m_mumu"] < args.sb2_max))
        )
        observed_sb, _ = np.histogram(data["m_mumu"][data_selection], bins=bin_edges)

        signal_selection = core.category_mask(signal["bdt_score"], low, high, is_last)
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
            f"[{category}] selected core order={order}, toys/truth={args.toys}, "
            f"sideband data={int(np.sum(observed_sb))}",
            flush=True,
        )

        for truth_name in truth_models:
            truth_result, expected_truth = core.fit_truth_sidebands(
                truth_name, observed_sb, sb_mask, centers, widths, args
            )
            asimov_fit, asimov_objective = core.fit_core_plus_signal(
                order, expected_truth, signal_counts, centers, widths, args
            )
            asimov_sigma = core.estimate_mu_error(asimov_objective, asimov_fit)
            asimov_mu = float(asimov_fit.x[-1])

            fitted_mu = []
            failed = 0
            boundary = 0
            for toy_index in range(args.toys):
                toy = rng.poisson(expected_truth).astype(float)
                result, _ = core.fit_core_plus_signal(
                    order, toy, signal_counts, centers, widths, args
                )
                mu_hat = float(result.x[-1])
                at_boundary = abs(mu_hat) > 19.95
                valid = bool(result.success and math.isfinite(mu_hat) and not at_boundary)
                if not result.success or not math.isfinite(mu_hat):
                    failed += 1
                if at_boundary:
                    boundary += 1
                if valid:
                    fitted_mu.append(mu_hat)
                details.append(
                    {
                        "category": category,
                        "truth_model": truth_name,
                        "core_order": order,
                        "toy": toy_index,
                        "fit_success": int(result.success),
                        "at_mu_boundary": int(at_boundary),
                        "mu_hat": mu_hat,
                        "asimov_mu": asimov_mu,
                        "asimov_mu_error": asimov_sigma,
                    }
                )

            failure_rate = failed / args.toys
            boundary_rate = boundary / args.toys
            if fitted_mu:
                mean_mu = float(np.mean(fitted_mu))
                median_mu = float(np.median(fitted_mu))
                rms_mu = float(np.std(fitted_mu, ddof=1)) if len(fitted_mu) > 1 else float("nan")
                q16 = quantile(fitted_mu, 0.16)
                q84 = quantile(fitted_mu, 0.84)
            else:
                mean_mu = median_mu = rms_mu = q16 = q84 = float("nan")

            mean_bias = abs(mean_mu) / asimov_sigma if asimov_sigma > 0.0 else float("nan")
            median_bias = abs(median_mu) / asimov_sigma if asimov_sigma > 0.0 else float("nan")
            worst_bias = max(mean_bias, median_bias)
            pull_width = rms_mu / asimov_sigma if asimov_sigma > 0.0 else float("nan")
            passed = bool(
                truth_result.success
                and asimov_fit.success
                and math.isfinite(worst_bias)
                and worst_bias < args.max_bias_over_sigma
                and failure_rate <= args.max_fit_failure_rate
                and boundary_rate <= args.max_fit_failure_rate
            )
            summaries.append(
                {
                    "category": category,
                    "truth_model": truth_name,
                    "core_order": order,
                    "toys_requested": args.toys,
                    "toys_valid": len(fitted_mu),
                    "fit_failure_rate": failure_rate,
                    "mu_boundary_rate": boundary_rate,
                    "expected_background": float(np.sum(expected_truth)),
                    "signal_yield": signal_yield,
                    "asimov_mu": asimov_mu,
                    "asimov_mu_error": asimov_sigma,
                    "mean_mu": mean_mu,
                    "median_mu": median_mu,
                    "rms_mu": rms_mu,
                    "q16_mu": q16,
                    "q84_mu": q84,
                    "mean_bias_over_sigma": mean_bias,
                    "median_bias_over_sigma": median_bias,
                    "pull_width_using_asimov_sigma": pull_width,
                    "threshold": args.max_bias_over_sigma,
                    "passes_toys": int(passed),
                }
            )
            print(
                f"  truth={truth_name:<12} valid={len(fitted_mu):4d}/{args.toys} "
                f"mean(mu)={mean_mu:+.3f} median={median_mu:+.3f} "
                f"sigmaA={asimov_sigma:.3f} worstBias={worst_bias:.3f} "
                f"pullWidth={pull_width:.3f} pass={int(passed)}",
                flush=True,
            )

    detail_path = os.path.join(args.output_dir, f"{args.era}_corepdf_bias_toys.csv")
    summary_path = os.path.join(
        args.output_dir, f"{args.era}_corepdf_bias_toy_summary.csv"
    )
    with open(detail_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(details[0].keys()))
        writer.writeheader()
        writer.writerows(details)
    with open(summary_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Wrote toy details: {detail_path}")
    print(f"Wrote toy summary: {summary_path}")
    print(
        "The 0.20-sigma criterion is provisional. Final acceptance must use the "
        "analysis-group-approved bias definition and toy count."
    )


if __name__ == "__main__":
    main()
