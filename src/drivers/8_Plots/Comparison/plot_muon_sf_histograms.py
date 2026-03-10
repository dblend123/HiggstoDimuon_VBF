#!/usr/bin/env python3
"""Plot muon SF distributions from a correctionlib JSON file.

This script mirrors the SF usage in:
  src/drivers/process_HiggsAnalysis_wCutsExtended_NoPairing_v0.cc

- ID key:   NUM_MediumID_DEN_TrackerMuons
- ISO key:  NUM_LoosePFIso_DEN_MediumID
- TRG key:  NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium

Notes:
- The C++ code uses abs(eta), clamps eta<=2.39 and pt>=15.01.
- Trigger SF is only applied when requested and pt>=26.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import correctionlib
import matplotlib.pyplot as plt
import numpy as np

DEFAULT_ID_KEY = "NUM_MediumID_DEN_TrackerMuons"
DEFAULT_ISO_KEY = "NUM_LoosePFIso_DEN_MediumID"
DEFAULT_TRG_KEY = "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium"


def evaluate_grid(corr, eta_vals: np.ndarray, pt_vals: np.ndarray) -> np.ndarray:
    out = np.empty(eta_vals.size * pt_vals.size, dtype=np.float64)
    i = 0
    for pt in pt_vals:
        for eta in eta_vals:
            out[i] = float(corr.evaluate(float(eta), float(pt), "nominal"))
            i += 1
    return out


def evaluate_grid_with_pt_cut(
    corr,
    eta_vals: np.ndarray,
    pt_vals: np.ndarray,
    pt_min: float,
) -> np.ndarray:
    kept_pts = pt_vals[pt_vals >= pt_min]
    if kept_pts.size == 0:
        return np.array([], dtype=np.float64)
    return evaluate_grid(corr, eta_vals, kept_pts)


def summarize(values: np.ndarray) -> Dict[str, float]:
    if values.size == 0:
        return {
            "n": 0,
            "min": np.nan,
            "max": np.nan,
            "mean": np.nan,
            "p16": np.nan,
            "p50": np.nan,
            "p84": np.nan,
            "outside_095_105": np.nan,
        }
    outside = np.mean((values < 0.95) | (values > 1.05))
    return {
        "n": int(values.size),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "p16": float(np.percentile(values, 16.0)),
        "p50": float(np.percentile(values, 50.0)),
        "p84": float(np.percentile(values, 84.0)),
        "outside_095_105": float(outside),
    }


def histogram_range(values: np.ndarray, pad_frac: float) -> Tuple[float, float]:
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    if math.isclose(vmin, vmax, rel_tol=1e-12, abs_tol=1e-12):
        delta = max(0.005, abs(vmin) * 0.01)
        return vmin - delta, vmax + delta
    span = vmax - vmin
    pad = max(1e-4, span * pad_frac)
    return vmin - pad, vmax + pad


def dedup_keep_order(values: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def infer_era_label(json_path: Path) -> str:
    stem = json_path.stem  # e.g. muon_Z12023
    prefix = "muon_Z"
    if stem.startswith(prefix) and len(stem) > len(prefix):
        return stem[len(prefix) :]
    return stem


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot muon ID/ISO/TRG SF distributions from correctionlib JSON"
    )
    parser.add_argument(
        "--json",
        default="src/drivers/7_Corrections/MUO/muon_Z2024.json",
        help="Path to muon SF JSON",
    )
    parser.add_argument("--id-key", default=DEFAULT_ID_KEY)
    parser.add_argument("--iso-key", default=DEFAULT_ISO_KEY)
    parser.add_argument("--trg-key", default=DEFAULT_TRG_KEY)
    parser.add_argument(
        "--include-trigger",
        action="store_true",
        help="Attempt to plot trigger SF as in getMuonSF(..., applyTriggerSF=true)",
    )
    parser.add_argument(
        "--all-num-keys",
        action="store_true",
        help="Plot all correction keys that start with NUM_ from the JSON",
    )
    parser.add_argument(
        "--extra-keys",
        nargs="*",
        default=[],
        help="Optional additional correction keys to include",
    )
    parser.add_argument("--eta-max", type=float, default=2.39)
    parser.add_argument("--eta-bins", type=int, default=120)
    parser.add_argument("--pt-min", type=float, default=15.01)
    parser.add_argument("--pt-max", type=float, default=200.0)
    parser.add_argument("--pt-bins", type=int, default=120)
    parser.add_argument(
        "--trg-pt-min",
        type=float,
        default=26.0,
        help="Trigger SF is evaluated only for pt >= this threshold",
    )
    parser.add_argument("--hist-bins", type=int, default=70)
    parser.add_argument(
        "--hist-min",
        type=float,
        default=None,
        help="Optional fixed histogram x-min. If omitted, per-plot min/max is used.",
    )
    parser.add_argument(
        "--hist-max",
        type=float,
        default=None,
        help="Optional fixed histogram x-max. If omitted, per-plot min/max is used.",
    )
    parser.add_argument(
        "--x-pad-frac",
        type=float,
        default=0.08,
        help="Fractional padding around min/max for auto x-range",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=3,
        help="Number of subplot columns when plotting multiple branches",
    )
    parser.add_argument(
        "--max-panels",
        type=int,
        default=0,
        help="Optional cap on number of plotted keys (0 means no cap)",
    )
    parser.add_argument(
        "--out",
        default="src/drivers/8_Plots/Comparison/muon_sf_histograms.png",
        help="Output image path",
    )

    args = parser.parse_args()

    json_path = Path(args.json)
    if not json_path.exists():
        raise FileNotFoundError(f"JSON file not found: {json_path}")

    cset = correctionlib.CorrectionSet.from_file(str(json_path))
    available_keys = set(cset.keys())

    eta_vals = np.linspace(0.0, args.eta_max, args.eta_bins)
    pt_vals = np.linspace(args.pt_min, args.pt_max, args.pt_bins)

    key_specs: List[Tuple[str, str, float | None]] = []

    if args.all_num_keys:
        for k in sorted(k for k in available_keys if k.startswith("NUM_")):
            key_specs.append((k, k, None))
    else:
        key_specs.append(("ID", args.id_key, None))
        key_specs.append(("ISO", args.iso_key, None))
        if args.include_trigger:
            key_specs.append(("TRG", args.trg_key, args.trg_pt_min))

    for k in args.extra_keys:
        key_specs.append((k, k, None))

    key_specs = dedup_keep_order(key_specs)
    if args.max_panels > 0:
        key_specs = key_specs[: args.max_panels]

    sf_values: Dict[str, np.ndarray] = {}
    for label, key, pt_cut in key_specs:
        if key not in available_keys:
            print(f"[WARN] Key not found in JSON: {key}")
            continue
        if pt_cut is None:
            vals = evaluate_grid(cset[key], eta_vals, pt_vals)
        else:
            vals = evaluate_grid_with_pt_cut(cset[key], eta_vals, pt_vals, pt_cut)
        sf_values[label] = vals

    if not sf_values:
        raise RuntimeError("No SF keys were evaluated. Check key names and JSON file.")

    n_panels = len(sf_values)
    ncols = max(1, min(args.ncols, n_panels))
    nrows = int(math.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.8 * ncols, 3.8 * nrows),
        squeeze=False,
    )
    flat_axes = [ax for row in axes for ax in row]

    print("\nMuon SF summary (nominal):")
    print(f"JSON: {json_path}")
    print(
        f"Grid: |eta| in [0, {args.eta_max}] with {args.eta_bins} points, "
        f"pt in [{args.pt_min}, {args.pt_max}] GeV with {args.pt_bins} points"
    )
    print(f"Plotted keys: {n_panels}")

    for ax, (label, values) in zip(flat_axes, sf_values.items()):
        stats = summarize(values)
        if args.hist_min is not None and args.hist_max is not None:
            hist_range = (args.hist_min, args.hist_max)
        else:
            hist_range = histogram_range(values, args.x_pad_frac)

        ax.hist(
            values,
            bins=args.hist_bins,
            range=hist_range,
            color="#2A6EA6",
            alpha=0.8,
            edgecolor="black",
            linewidth=0.4,
        )
        ax.axvline(stats["mean"], color="#111111", linestyle="-", linewidth=1.2)
        ax.axvline(stats["min"], color="#2E86AB", linestyle="--", linewidth=1.1)
        ax.axvline(stats["max"], color="#F18F01", linestyle="--", linewidth=1.1)
        ax.set_title(
            f"{label} SF\n"
            f"mean={stats['mean']:.4f}, min={stats['min']:.4f}, max={stats['max']:.4f}",
            fontsize=10,
        )
        ax.set_xlabel("Scale factor value")
        ax.set_ylabel("Counts on eta-pt grid")
        ax.grid(alpha=0.25)
        ax.set_xlim(hist_range)

        print(
            f"  {label:>3}: n={stats['n']}, mean={stats['mean']:.5f}, "
            f"min={stats['min']:.5f}, max={stats['max']:.5f}, "
            f"p16/p50/p84=({stats['p16']:.5f}/{stats['p50']:.5f}/{stats['p84']:.5f}), "
            f"outside[0.95,1.05]={100.0 * stats['outside_095_105']:.2f}%"
        )

    for ax in flat_axes[n_panels:]:
        ax.set_visible(False)

    era_label = infer_era_label(json_path)
    fig.suptitle(
        f"Muon SF distributions from correctionlib JSON (Era: {era_label})",
        fontsize=12,
    )
    fig.tight_layout()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    print(f"\nSaved plot: {out_path}")


if __name__ == "__main__":
    main()
