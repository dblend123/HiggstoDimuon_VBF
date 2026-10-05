#!/usr/bin/env python3
"""Fit candidate analytic background models to blinded BDT-category sidebands."""

import argparse
import csv
import math
import os
from array import array

import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", default="2024")
    parser.add_argument(
        "--data",
        default="../Outputs_nomjj_newsample/2024/2024_data_scored.root",
    )
    parser.add_argument("--treename", default="mva")
    parser.add_argument(
        "--output-dir",
        default="../Outputs_nomjj_newsample/BackgroundFits/2024",
    )
    parser.add_argument("--category-edges", default="0.60,0.85,0.95,1.00")
    parser.add_argument("--mass-min", type=float, default=110.0)
    parser.add_argument("--mass-max", type=float, default=150.0)
    parser.add_argument("--sb1-min", type=float, default=110.0)
    parser.add_argument("--sb1-max", type=float, default=115.0)
    parser.add_argument("--sb2-min", type=float, default=135.0)
    parser.add_argument("--sb2-max", type=float, default=150.0)
    parser.add_argument("--plot-bin-width", type=float, default=1.0)
    return parser.parse_args()


def parse_edges(text):
    values = np.asarray([float(item) for item in text.split(",")], dtype=float)
    if len(values) < 2 or np.any(np.diff(values) <= 0.0):
        raise ValueError("Category edges must be strictly increasing")
    return values


def load_data(path, treename):
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    tree = root_file.Get(treename)
    if not tree:
        raise RuntimeError(f"Missing tree {treename} in {path}")
    required = ("m_mumu", "bdt_score")
    missing = [name for name in required if not tree.GetBranch(name)]
    root_file.Close()
    if missing:
        raise RuntimeError(f"Missing branches: {missing}")
    arrays = ROOT.RDataFrame(treename, path).AsNumpy(list(required))
    return {name: np.asarray(values, dtype=np.float64) for name, values in arrays.items()}


def make_dataset(name, mass, values):
    dataset = ROOT.RooDataSet(name, name, ROOT.RooArgSet(mass))
    for value in values:
        mass.setVal(float(value))
        dataset.add(ROOT.RooArgSet(mass))
    return dataset


def make_candidates(mass, category):
    """Return (name, pdf, floating_parameters, owned_objects) candidates."""
    candidates = []

    slope = ROOT.RooRealVar(f"exp_slope_{category}", "exponential slope", -0.03, -1.0, 0.0)
    exponential = ROOT.RooExponential(
        f"exp_{category}", "exponential", mass, slope
    )
    candidates.append(("exponential", exponential, [slope], [slope, exponential]))

    exponent = ROOT.RooRealVar(
        f"power_exponent_{category}", "power-law exponent", 4.0, 0.01, 30.0
    )
    power = ROOT.RooGenericPdf(
        f"power_{category}", "pow(@0,-@1)", ROOT.RooArgList(mass, exponent)
    )
    candidates.append(("powerlaw", power, [exponent], [exponent, power]))

    centered_mass = ROOT.RooFormulaVar(
        f"centered_mass_{category}", "(@0-130.0)/20.0", ROOT.RooArgList(mass)
    )
    ep_a = ROOT.RooRealVar(f"exppoly_a_{category}", "linear coefficient", -0.5, -10.0, 10.0)
    ep_b = ROOT.RooRealVar(f"exppoly_b_{category}", "quadratic coefficient", 0.0, -10.0, 10.0)
    exp_poly = ROOT.RooGenericPdf(
        f"exppoly2_{category}",
        "exp(@1*@0+@2*@0*@0)",
        ROOT.RooArgList(centered_mass, ep_a, ep_b),
    )
    candidates.append(
        ("exp_poly2", exp_poly, [ep_a, ep_b], [centered_mass, ep_a, ep_b, exp_poly])
    )

    for degree in (2, 3, 4):
        coefficients = []
        owned = []
        first = ROOT.RooRealVar(
            f"bern{degree}_c0_{category}", "fixed first coefficient", 1.0
        )
        first.setConstant(True)
        coefficients.append(first)
        owned.append(first)
        floating = []
        for index in range(1, degree + 1):
            coefficient = ROOT.RooRealVar(
                f"bern{degree}_c{index}_{category}",
                f"Bernstein coefficient {index}",
                1.0,
                1.0e-5,
                100.0,
            )
            coefficients.append(coefficient)
            floating.append(coefficient)
            owned.append(coefficient)
        pdf = ROOT.RooBernstein(
            f"bernstein{degree}_{category}",
            f"Bernstein degree {degree}",
            mass,
            ROOT.RooArgList(*coefficients),
        )
        owned.append(pdf)
        candidates.append((f"bernstein{degree}", pdf, floating, owned))

    return candidates


def pdf_value(pdf, mass, value):
    mass.setVal(float(value))
    return max(0.0, float(pdf.getVal(ROOT.RooArgSet(mass))))


def expected_bin_counts(pdf, mass, bin_edges, sideband_mask, observed_sideband):
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)
    density = np.asarray([pdf_value(pdf, mass, x) for x in centers])
    norm = float(np.sum(density[sideband_mask] * widths[sideband_mask]))
    if not math.isfinite(norm) or norm <= 0.0:
        return np.zeros_like(centers)
    return observed_sideband * density * widths / norm


def fit_category(category, values, args, output_root, metrics):
    mass = ROOT.RooRealVar(
        f"m_mumu_{category}", "m_{#mu#mu}", args.mass_min, args.mass_max, "GeV"
    )
    mass.setRange("lowSB", args.sb1_min, args.sb1_max)
    mass.setRange("highSB", args.sb2_min, args.sb2_max)
    mass.setRange("full", args.mass_min, args.mass_max)

    dataset = make_dataset(f"data_sidebands_{category}", mass, values)
    candidates = make_candidates(mass, category)

    number_bins = int(round((args.mass_max - args.mass_min) / args.plot_bin_width))
    bin_edges = np.linspace(args.mass_min, args.mass_max, number_bins + 1)
    observed, _ = np.histogram(values, bins=bin_edges)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    sideband_mask = ((centers >= args.sb1_min) & (centers < args.sb1_max)) | (
        (centers >= args.sb2_min) & (centers < args.sb2_max)
    )

    canvas = ROOT.TCanvas(f"canvas_{category}", category, 900, 750)
    frame = mass.frame(ROOT.RooFit.Title(f"{category}: blinded sideband candidates"))
    dataset.plotOn(frame, ROOT.RooFit.Binning(number_bins), ROOT.RooFit.Name("data"))

    colors = (ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kOrange + 7, ROOT.kCyan + 2)
    fitted_objects = []

    for index, (model_name, pdf, parameters, owned) in enumerate(candidates):
        result = pdf.fitTo(
            dataset,
            ROOT.RooFit.Range("lowSB,highSB"),
            ROOT.RooFit.Save(True),
            ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.Strategy(1),
            ROOT.RooFit.SumW2Error(False),
        )
        fitted_objects.extend(owned)
        fitted_objects.append(result)

        expected = expected_bin_counts(
            pdf, mass, bin_edges, sideband_mask, float(len(values))
        )
        valid = sideband_mask & (expected > 0.0)
        chi2 = float(np.sum((observed[valid] - expected[valid]) ** 2 / expected[valid]))
        parameter_count = len(parameters)
        populated_bins = int(np.count_nonzero(valid))
        ndf = max(1, populated_bins - parameter_count)
        p_value = float(ROOT.TMath.Prob(chi2, ndf))
        minimum_nll = float(result.minNll())
        aic = 2.0 * parameter_count + 2.0 * minimum_nll
        bic = parameter_count * math.log(max(1, len(values))) + 2.0 * minimum_nll

        metrics.append(
            {
                "category": category,
                "model": model_name,
                "sideband_events": len(values),
                "parameters": parameter_count,
                "fit_status": int(result.status()),
                "covariance_quality": int(result.covQual()),
                "minimum_nll": minimum_nll,
                "aic": aic,
                "bic": bic,
                "chi2": chi2,
                "ndf": ndf,
                "chi2_ndf": chi2 / ndf,
                "p_value": p_value,
            }
        )

        pdf.plotOn(
            frame,
            ROOT.RooFit.Range("full"),
            ROOT.RooFit.NormRange("lowSB,highSB"),
            ROOT.RooFit.LineColor(colors[index]),
            ROOT.RooFit.LineWidth(2),
            ROOT.RooFit.Name(model_name),
        )

        output_root.cd()
        result.Write(f"fitresult_{category}_{model_name}")

    frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    frame.GetYaxis().SetTitle(f"Events / {args.plot_bin_width:g} GeV")
    frame.Draw()

    blind_box = ROOT.TBox(args.sb1_max, 0.0, args.sb2_min, frame.GetMaximum())
    blind_box.SetFillColorAlpha(ROOT.kGray + 1, 0.22)
    blind_box.SetLineColor(ROOT.kGray + 2)
    blind_box.Draw("same")

    legend = ROOT.TLegend(0.56, 0.57, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(frame.findObject("data"), "Observed sidebands", "lep")
    for model_name, _, _, _ in candidates:
        legend.AddEntry(frame.findObject(model_name), model_name, "l")
    legend.Draw()

    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextSize(0.035)
    label.DrawLatex(0.14, 0.92, f"{args.era} {category}  (115-135 GeV blinded)")

    pdf_path = os.path.join(args.output_dir, f"{args.era}_{category}_sideband_models.pdf")
    canvas.SaveAs(pdf_path)
    output_root.cd()
    canvas.Write()
    dataset.Write()

    print(f"[{category}] sideband events = {len(values)}")
    for row in metrics:
        if row["category"] == category:
            print(
                f"  {row['model']:<12} status={row['fit_status']} "
                f"cov={row['covariance_quality']} chi2/ndf={row['chi2_ndf']:.3f} "
                f"p={row['p_value']:.4g} AIC={row['aic']:.2f}"
            )


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    category_edges = parse_edges(args.category_edges)
    arrays = load_data(args.data, args.treename)

    output_root_path = os.path.join(
        args.output_dir, f"{args.era}_sideband_background_fits.root"
    )
    output_root = ROOT.TFile(output_root_path, "RECREATE")
    if not output_root or output_root.IsZombie():
        raise RuntimeError(f"Cannot create {output_root_path}")

    metrics = []
    mass_values = arrays["m_mumu"]
    scores = arrays["bdt_score"]
    sideband = ((mass_values >= args.sb1_min) & (mass_values < args.sb1_max)) | (
        (mass_values >= args.sb2_min) & (mass_values < args.sb2_max)
    )

    for index, (low, high) in enumerate(zip(category_edges[:-1], category_edges[1:])):
        category = f"cat{index}"
        if index == len(category_edges) - 2:
            score_selection = (scores >= low) & (scores <= high)
        else:
            score_selection = (scores >= low) & (scores < high)
        selected = mass_values[score_selection & sideband]
        if len(selected) == 0:
            raise RuntimeError(f"No sideband data in {category}")
        fit_category(category, selected, args, output_root, metrics)

    output_root.Close()

    metrics_path = os.path.join(
        args.output_dir, f"{args.era}_sideband_background_fit_metrics.csv"
    )
    fields = [
        "category", "model", "sideband_events", "parameters", "fit_status",
        "covariance_quality", "minimum_nll", "aic", "bic", "chi2", "ndf",
        "chi2_ndf", "p_value",
    ]
    with open(metrics_path, "w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(metrics)

    print(f"Wrote fit results: {output_root_path}")
    print(f"Wrote metrics:     {metrics_path}")
    print("Candidate selection is not final until bias/spurious-signal tests pass.")


if __name__ == "__main__":
    main()
