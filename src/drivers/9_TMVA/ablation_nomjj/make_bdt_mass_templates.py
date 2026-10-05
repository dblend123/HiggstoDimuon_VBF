#!/usr/bin/env python3
"""Build blinded, category-based dimuon-mass templates for the 2024 VBF BDT analysis."""

import argparse
import csv
import math
import os
from typing import Dict, List, Sequence, Tuple

import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)


PROCESSES = ("VBF", "DY", "EWK", "TT", "VV", "VVV")

WEIGHT_VARIATIONS = {
    "puUp": "weight_puUp",
    "puDown": "weight_puDown",
    "muSFUp": "weight_muSFUp",
    "muSFDown": "weight_muSFDown",
    "triggerUp": "weight_triggerUp",
    "triggerDown": "weight_triggerDown",
    "btagSFUp": "weight_btagSFUp",
    "btagSFDown": "weight_btagSFDown",
    "btagEffStatUp": "weight_btagEffStatUp",
    "btagEffStatDown": "weight_btagEffStatDown",
    "pdfUp": "weight_pdfUp",
    "pdfDown": "weight_pdfDown",
    "isrUp": "weight_isrUp",
    "isrDown": "weight_isrDown",
    "fsrUp": "weight_fsrUp",
    "fsrDown": "weight_fsrDown",
}

SCALE_BRANCHES = (
    "weight_scale0",
    "weight_scale1",
    "weight_scale3",
    "weight_scale5",
    "weight_scale7",
    "weight_scale8",
)

SHIFTED_VARIATIONS = ("jecUp", "jecDown", "jerUp", "jerDown")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", default="2024")
    parser.add_argument(
        "--input-dir",
        default="../Outputs_nomjj_newsample",
        help="Directory containing nominal scored files and Systematics/<era>",
    )
    parser.add_argument(
        "--output-dir",
        default="../Outputs_nomjj_newsample/Templates/2024",
    )
    parser.add_argument(
        "--category-edges",
        default="0.60,0.85,0.95,1.00",
        help="Comma-separated lower edge followed by all upper edges",
    )
    parser.add_argument("--mass-min", type=float, default=110.0)
    parser.add_argument("--mass-max", type=float, default=150.0)
    parser.add_argument("--mass-nbins", type=int, default=180)
    parser.add_argument("--sb1-min", type=float, default=110.0)
    parser.add_argument("--sb1-max", type=float, default=115.0)
    parser.add_argument("--sb2-min", type=float, default=135.0)
    parser.add_argument("--sb2-max", type=float, default=150.0)
    parser.add_argument(
        "--write-scale-components",
        action="store_true",
        help="Also save the six component scale templates used by the envelope",
    )
    return parser.parse_args()


def parse_edges(text: str) -> np.ndarray:
    edges = np.asarray([float(x.strip()) for x in text.split(",")], dtype=float)
    if edges.size < 2 or not np.all(np.diff(edges) > 0.0):
        raise ValueError("Category edges must be strictly increasing")
    if edges[0] < 0.0 or edges[-1] > 1.0:
        raise ValueError("Category edges must stay within [0,1]")
    return edges


def require_file(path: str) -> None:
    if not os.path.isfile(path):
        raise FileNotFoundError(path)


def load_arrays(path: str, tree_name: str, branches: Sequence[str]) -> Dict[str, np.ndarray]:
    require_file(path)
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {path}")
    tree = root_file.Get(tree_name)
    if not tree:
        root_file.Close()
        raise RuntimeError(f"Missing tree '{tree_name}' in {path}")
    present = {branch.GetName() for branch in tree.GetListOfBranches()}
    missing = [branch for branch in branches if branch not in present]
    root_file.Close()
    if missing:
        raise RuntimeError(f"Missing branches in {path}:{tree_name}: {missing}")
    arrays = ROOT.RDataFrame(tree_name, path).AsNumpy(list(branches))
    return {key: np.asarray(value, dtype=np.float64) for key, value in arrays.items()}


def nominal_path(input_dir: str, era: str, process: str) -> str:
    stem = "signal" if process == "VBF" else process
    return os.path.join(input_dir, era, f"{era}_scored_{stem}.root")


def shifted_path(input_dir: str, era: str, process: str, variation: str) -> str:
    return os.path.join(
        input_dir,
        "Systematics",
        era,
        process,
        f"{era}_scored_{process}_{variation}.root",
    )


def category_mask(score: np.ndarray, low: float, high: float, final: bool) -> np.ndarray:
    if final:
        return (score >= low) & (score <= high)
    return (score >= low) & (score < high)


def histogram_arrays(
    mass: np.ndarray,
    weight: np.ndarray,
    mask: np.ndarray,
    mass_edges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    finite = mask & np.isfinite(mass) & np.isfinite(weight)
    contents, _ = np.histogram(mass[finite], bins=mass_edges, weights=weight[finite])
    sumw2, _ = np.histogram(
        mass[finite], bins=mass_edges, weights=weight[finite] * weight[finite]
    )
    return contents.astype(float), sumw2.astype(float)


def make_hist(
    name: str,
    contents: np.ndarray,
    sumw2: np.ndarray,
    mass_min: float,
    mass_max: float,
) -> ROOT.TH1D:
    hist = ROOT.TH1D(name, name, len(contents), mass_min, mass_max)
    hist.SetDirectory(0)
    hist.Sumw2()
    for index, (content, variance) in enumerate(zip(contents, sumw2), start=1):
        hist.SetBinContent(index, float(content))
        hist.SetBinError(index, math.sqrt(max(float(variance), 0.0)))
    return hist


def record_diagnostic(
    rows: List[Dict[str, object]],
    category: str,
    process: str,
    variation: str,
    hist: ROOT.TH1,
) -> None:
    contents = np.asarray(
        [hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)], dtype=float
    )
    rows.append(
        {
            "category": category,
            "process": process,
            "variation": variation,
            "integral": float(np.sum(contents)),
            "negative_bins": int(np.count_nonzero(contents < 0.0)),
            "zero_bins": int(np.count_nonzero(contents == 0.0)),
            "minimum_bin": float(np.min(contents)) if contents.size else float("nan"),
            "maximum_bin": float(np.max(contents)) if contents.size else float("nan"),
        }
    )


def write_hist(
    output_file: ROOT.TFile,
    directory: ROOT.TDirectory,
    hist: ROOT.TH1,
    diagnostics: List[Dict[str, object]],
    category: str,
    process: str,
    variation: str,
) -> None:
    directory.cd()
    hist.Write(hist.GetName(), ROOT.TObject.kOverwrite)
    record_diagnostic(diagnostics, category, process, variation, hist)
    output_file.cd()


def main() -> None:
    args = parse_args()
    category_edges = parse_edges(args.category_edges)
    mass_edges = np.linspace(args.mass_min, args.mass_max, args.mass_nbins + 1)
    os.makedirs(args.output_dir, exist_ok=True)

    output_path = os.path.join(args.output_dir, f"{args.era}_bdt_mass_templates.root")
    diagnostics_path = os.path.join(
        args.output_dir, f"{args.era}_bdt_mass_template_diagnostics.csv"
    )

    nominal_branches = ["m_mumu", "bdt_score", "weight"]
    nominal_branches.extend(WEIGHT_VARIATIONS.values())
    nominal_branches.extend(SCALE_BRANCHES)
    nominal_branches.extend(("weight_zptUp", "weight_zptDown"))

    nominal_arrays = {
        process: load_arrays(
            nominal_path(args.input_dir, args.era, process), "mva", nominal_branches
        )
        for process in PROCESSES
    }

    data_path = os.path.join(args.input_dir, args.era, f"{args.era}_data_scored.root")
    data_arrays = load_arrays(data_path, "mva", ["m_mumu", "bdt_score"])

    shifted_arrays: Dict[Tuple[str, str], Dict[str, np.ndarray]] = {}
    for process in PROCESSES:
        for variation in SHIFTED_VARIATIONS:
            shifted_arrays[(process, variation)] = load_arrays(
                shifted_path(args.input_dir, args.era, process, variation),
                f"mva_{variation}",
                ["m_mumu", "bdt_score", "weight"],
            )

    output_file = ROOT.TFile(output_path, "RECREATE")
    if not output_file or output_file.IsZombie():
        raise RuntimeError(f"Cannot create {output_path}")

    output_file.cd()
    ROOT.TNamed("analysis_era", args.era).Write()
    ROOT.TNamed("category_edges", ",".join(map(str, category_edges))).Write()
    ROOT.TNamed("mass_binning", f"{args.mass_nbins},{args.mass_min},{args.mass_max}").Write()
    ROOT.TNamed(
        "data_blinding",
        f"DATA filled only in [{args.sb1_min},{args.sb1_max}) and "
        f"[{args.sb2_min},{args.sb2_max}); central region omitted",
    ).Write()

    diagnostics: List[Dict[str, object]] = []

    for category_index, (low, high) in enumerate(
        zip(category_edges[:-1], category_edges[1:])
    ):
        category = f"cat{category_index}"
        directory = output_file.mkdir(category)
        final_category = category_index == len(category_edges) - 2

        data_score_mask = category_mask(
            data_arrays["bdt_score"], low, high, final_category
        )
        data_mass = data_arrays["m_mumu"]
        blinded_mass_mask = (
            ((data_mass >= args.sb1_min) & (data_mass < args.sb1_max))
            | ((data_mass >= args.sb2_min) & (data_mass < args.sb2_max))
        )
        data_contents, data_sumw2 = histogram_arrays(
            data_mass,
            np.ones_like(data_mass),
            data_score_mask & blinded_mass_mask,
            mass_edges,
        )
        data_hist = make_hist(
            "data_obs_blinded", data_contents, data_sumw2, args.mass_min, args.mass_max
        )
        write_hist(
            output_file,
            directory,
            data_hist,
            diagnostics,
            category,
            "data_obs_blinded",
            "nominal",
        )

        for process in PROCESSES:
            arrays = nominal_arrays[process]
            score_mask = category_mask(
                arrays["bdt_score"], low, high, final_category
            )
            mass_mask = (
                (arrays["m_mumu"] >= args.mass_min)
                & (arrays["m_mumu"] < args.mass_max)
            )
            selection = score_mask & mass_mask

            nominal_contents, nominal_sumw2 = histogram_arrays(
                arrays["m_mumu"], arrays["weight"], selection, mass_edges
            )
            nominal_hist = make_hist(
                process,
                nominal_contents,
                nominal_sumw2,
                args.mass_min,
                args.mass_max,
            )
            write_hist(
                output_file,
                directory,
                nominal_hist,
                diagnostics,
                category,
                process,
                "nominal",
            )

            process_weight_variations = dict(WEIGHT_VARIATIONS)
            if process == "DY":
                process_weight_variations.update(
                    {"zptUp": "weight_zptUp", "zptDown": "weight_zptDown"}
                )

            for variation, branch in process_weight_variations.items():
                contents, sumw2 = histogram_arrays(
                    arrays["m_mumu"], arrays[branch], selection, mass_edges
                )
                hist = make_hist(
                    f"{process}_{variation}",
                    contents,
                    sumw2,
                    args.mass_min,
                    args.mass_max,
                )
                write_hist(
                    output_file,
                    directory,
                    hist,
                    diagnostics,
                    category,
                    process,
                    variation,
                )

            scale_contents = []
            scale_sumw2 = []
            for branch in SCALE_BRANCHES:
                contents, sumw2 = histogram_arrays(
                    arrays["m_mumu"], arrays[branch], selection, mass_edges
                )
                scale_contents.append(contents)
                scale_sumw2.append(sumw2)
                if args.write_scale_components:
                    component = branch.replace("weight_", "")
                    hist = make_hist(
                        f"{process}_{component}",
                        contents,
                        sumw2,
                        args.mass_min,
                        args.mass_max,
                    )
                    write_hist(
                        output_file,
                        directory,
                        hist,
                        diagnostics,
                        category,
                        process,
                        component,
                    )

            # Include nominal in the binwise envelope so scaleUp is never below
            # nominal and scaleDown is never above it.
            scale_contents_array = np.vstack([nominal_contents] + scale_contents)
            scale_sumw2_array = np.vstack([nominal_sumw2] + scale_sumw2)
            up_indices = np.argmax(scale_contents_array, axis=0)
            down_indices = np.argmin(scale_contents_array, axis=0)
            bin_indices = np.arange(args.mass_nbins)

            for variation, indices in (
                ("scaleUp", up_indices),
                ("scaleDown", down_indices),
            ):
                contents = scale_contents_array[indices, bin_indices]
                sumw2 = scale_sumw2_array[indices, bin_indices]
                hist = make_hist(
                    f"{process}_{variation}",
                    contents,
                    sumw2,
                    args.mass_min,
                    args.mass_max,
                )
                write_hist(
                    output_file,
                    directory,
                    hist,
                    diagnostics,
                    category,
                    process,
                    variation,
                )

            for variation in SHIFTED_VARIATIONS:
                shifted = shifted_arrays[(process, variation)]
                shifted_score_mask = category_mask(
                    shifted["bdt_score"], low, high, final_category
                )
                shifted_mass_mask = (
                    (shifted["m_mumu"] >= args.mass_min)
                    & (shifted["m_mumu"] < args.mass_max)
                )
                contents, sumw2 = histogram_arrays(
                    shifted["m_mumu"],
                    shifted["weight"],
                    shifted_score_mask & shifted_mass_mask,
                    mass_edges,
                )
                hist = make_hist(
                    f"{process}_{variation}",
                    contents,
                    sumw2,
                    args.mass_min,
                    args.mass_max,
                )
                write_hist(
                    output_file,
                    directory,
                    hist,
                    diagnostics,
                    category,
                    process,
                    variation,
                )

    output_file.Write()
    output_file.Close()

    fieldnames = [
        "category",
        "process",
        "variation",
        "integral",
        "negative_bins",
        "zero_bins",
        "minimum_bin",
        "maximum_bin",
    ]
    with open(diagnostics_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(diagnostics)

    print(f"Wrote templates:   {output_path}")
    print(f"Wrote diagnostics: {diagnostics_path}")
    print(f"Categories: {category_edges.tolist()}")
    print("DATA signal region [115,135) remains empty by construction.")


if __name__ == "__main__":
    main()
