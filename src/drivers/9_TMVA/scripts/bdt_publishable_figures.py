#!/usr/bin/env python3
"""
Build a publication-style, figure-only PDF from the same inputs used by
bdt_event_categorization.py.

This script reuses the categorization pipeline and emits only polished figures:
- BDT score (sidebands) Data vs MC
- Blinded dimuon mass spectrum
- Per-category dimuon mass with best sideband background fit

Each figure is placed on its own PDF page with legend-aware layout.
"""

import argparse
import os
from typing import Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import bdt_event_categorization as bec


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Publishable figure-only PDF from BDT categorization inputs")

    p.add_argument("--outputs-dir", default=bec.DEFAULT_OUTPUTS_DIR)
    p.add_argument("--outdir", default=bec.DEFAULT_OUTDIR)
    p.add_argument("--output-name", default="bdt_publishable_figures.pdf")

    p.add_argument("--eras", nargs="*", default=None)

    p.add_argument("--treename", default="mva")
    p.add_argument("--bdt-branch", default="bdt_score")
    p.add_argument("--weight-branch", default="weight")
    p.add_argument("--mass-branch", default="auto")
    p.add_argument("--extra-cut", default="")

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

    p.add_argument("--abs-weights", dest="abs_weights", action="store_true", default=True)
    p.add_argument("--signed-weights", dest="abs_weights", action="store_false")

    p.add_argument("--score-logy", action="store_true")
    p.add_argument("--mc-normalization", choices=["none", "sidebands"], default="sidebands")

    p.add_argument("--category-edges", default="")
    p.add_argument("--category-eras", nargs="*", default=None)
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


def apply_publish_style() -> None:
    bec.configure_plot_style()
    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "legend.fontsize": 9,
        "figure.titlesize": 15,
        "savefig.dpi": 200,
    })


def _relayout_legends(fig: plt.Figure, reserve_right: float = 0.80) -> None:
    for ax in fig.axes:
        handles, labels = ax.get_legend_handles_labels()
        if not handles or not labels:
            continue
        ax.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            frameon=True,
            framealpha=0.9,
            ncol=1,
        )

    if fig._suptitle is not None:
        fig._suptitle.set_y(0.992)

    fig.tight_layout(rect=[0.02, 0.02, reserve_right, 0.965])


def _build_payload(args: argparse.Namespace, logger) -> Tuple[dict, float, np.ndarray, tuple]:
    args.outputs_dir = bec.resolve_input_dir(args.outputs_dir, bec.TMVA_DIR)
    if not os.path.isabs(args.outdir):
        args.outdir = os.path.normpath(os.path.abspath(args.outdir))
    os.makedirs(args.outdir, exist_ok=True)

    data_map, mc_map = bec.parse_user_inputs(args, logger)
    eras = sorted(set(list(data_map.keys()) + list(mc_map.keys())))
    if not eras:
        raise SystemExit("[error] No eras found from user input or auto-discovery.")

    era_payloads = {}
    for era in eras:
        payload = bec.load_era_payload(
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

    edges = bec.optimize_category_edges(args, eras, logger)
    interior = np.asarray(edges[1:-1], dtype=np.float64)
    category_band = None
    if interior.size >= 1:
        category_band = (float(np.min(interior)), float(np.max(interior)))

    run3_data_side_scores = bec.combine_data_side_scores(era_payloads)
    run3_data_side_mass = bec.combine_data_side_masses(era_payloads, args)

    run3_mc_score_side = bec.combine_mc_score_groups(era_payloads, region="side")
    run3_mc_mass_full = bec.combine_mc_mass_groups(era_payloads, args, region="full")

    if run3_data_side_scores.size == 0:
        raise SystemExit("[error] No data events found in sidebands for the selected eras.")
    if not run3_mc_score_side:
        raise SystemExit("[error] No sideband MC events found for BDT plotting.")

    mc_scale = 1.0
    if args.mc_normalization == "sidebands":
        mc_scale = bec.compute_sideband_mc_scale(run3_data_side_scores, run3_mc_score_side)
        logger.info("MC normalization (sidebands): scale factor = %.6f", mc_scale)

    run3_mc_score_side = bec.scale_mc_groups(run3_mc_score_side, mc_scale)
    run3_mc_mass_full = bec.scale_mc_groups(run3_mc_mass_full, mc_scale)
    run3_mc_score_display = bec.group_backgrounds_for_display(run3_mc_score_side)

    score_bins = np.linspace(args.score_min, args.score_max, args.score_nbins + 1)
    mass_bins = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)

    return {
        "era_payloads": era_payloads,
        "run3_data_side_scores": run3_data_side_scores,
        "run3_data_side_mass": run3_data_side_mass,
        "run3_mc_score_display": run3_mc_score_display,
        "run3_mc_mass_full": run3_mc_mass_full,
    }, mc_scale, edges, category_band, score_bins, mass_bins


def main() -> None:
    args = parse_args()
    logger = bec.setup_logger(args.log_level)
    apply_publish_style()

    payload, mc_scale, edges, category_band, score_bins, mass_bins = _build_payload(args, logger)

    era_payloads = payload["era_payloads"]

    fig_score, *_ = bec.make_stacked_ratio_figure(
        title="BDT score in sidebands: Data vs stacked MC",
        xlabel="BDT score",
        bins=score_bins,
        data_values=payload["run3_data_side_scores"],
        mc_groups=payload["run3_mc_score_display"],
        logy=args.score_logy,
        category_band=category_band,
        show_mc_unc_band=False,
        shade_category_band=False,
        ratio_ylabel="Data/Bkg",
        annotation_lines=[
            f"Eras: {', '.join(sorted(era_payloads.keys()))}",
            f"Sidebands: ({args.sb1_min:.0f}, {args.sb1_max:.0f}) U ({args.sb2_min:.0f}, {args.sb2_max:.0f}) GeV",
        ],
    )

    fig_mass, *_ = bec.make_stacked_ratio_figure(
        title="Blinded m_mumu spectrum: sideband data points with MC expectation",
        xlabel="m_mumu [GeV]",
        bins=mass_bins,
        data_values=payload["run3_data_side_mass"],
        mc_groups=payload["run3_mc_mass_full"],
        logy=False,
        blind_band=(args.blind_min, args.blind_max),
        show_mc_unc_band=False,
    )

    fig_cat_mass, _fit_cols, _fit_rows = bec.build_category_mass_grid(
        era_payloads=era_payloads,
        edges=edges,
        args=args,
        mass_bins=mass_bins,
        mc_weight_scale=mc_scale,
    )

    _relayout_legends(fig_score, reserve_right=0.80)
    _relayout_legends(fig_mass, reserve_right=0.80)
    _relayout_legends(fig_cat_mass, reserve_right=0.84)

    out_path = args.output_name if os.path.isabs(args.output_name) else os.path.join(args.outdir, args.output_name)
    with PdfPages(out_path) as pdf:
        pdf.savefig(fig_score)
        pdf.savefig(fig_mass)
        pdf.savefig(fig_cat_mass)

    plt.close(fig_score)
    plt.close(fig_mass)
    plt.close(fig_cat_mass)

    logger.info("Wrote publishable figure PDF: %s", out_path)
    print(f"[done] Wrote publishable figures: {out_path}")


if __name__ == "__main__":
    main()
