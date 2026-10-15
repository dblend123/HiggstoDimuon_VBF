#!/usr/bin/env python3
"""Build the blinded 2024 VBF BDT RooWorkspace and a Combine datacard.

The pre-unblinding data_obs object is hybrid: real data in the sidebands and
the fitted background Asimov expectation in [115,135) GeV.  This avoids
interpreting deliberately empty blinded bins as observations of zero events.
"""

import argparse
import math
import os
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(100000)
ROOT.Math.MinimizerOptions.SetDefaultMaxIterations(100000)

CATEGORIES = ("cat0", "cat1", "cat2")
BACKGROUND_MODELS = {
    "cat0": "SumExp",
    # FEWZxBernstein has an unidentifiable p0 direction in Combine.
    # LandaUxBernstein passed all seven 500-toy truth-model tests.
    "cat1": "LandaUxBernstein",
    # Both PowerLaw and SumExp passed all seven 500-toy generation models in
    # the highest-score category.  Use SumExp here because the three-parameter
    # PowerLaw is nearly non-identifiable over 110--150 GeV and produced a
    # forced-positive-definite covariance in Combine.
    "cat2": "SumExp",
}
SHAPE_NUISANCES = (
    "pu", "muSF", "trigger", "btagSF", "btagEffStat",
    "pdf", "isr", "fsr", "scale", "jec", "jer",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--era", default="2024")
    parser.add_argument(
        "--templates",
        default="../Outputs_nomjj_newsample/Templates/2024/2024_bdt_mass_templates.root",
    )
    parser.add_argument(
        "--output-dir",
        default="../Outputs_nomjj_newsample/Workspace/2024",
    )
    parser.add_argument("--mass-min", type=float, default=110.0)
    parser.add_argument("--mass-max", type=float, default=150.0)
    parser.add_argument("--blind-min", type=float, default=115.0)
    parser.add_argument("--blind-max", type=float, default=135.0)
    parser.add_argument(
        "--signal-rate-scale",
        type=float,
        default=1.0,
        help="Optional scale applied by the datacard to all VBF rates.",
    )
    return parser.parse_args()


def ws_import(workspace, obj):
    # Plain import is the most portable PyROOT form across the LCG ROOT builds.
    # RooWorkspace reuses already-present observables with matching names.
    getattr(workspace, "import")(obj)


def require_hist(root_file, path):
    hist = root_file.Get(path)
    if not hist or not hist.InheritsFrom("TH1"):
        raise RuntimeError(f"Missing histogram: {path}")
    clone = hist.Clone(path.replace("/", "__"))
    clone.SetDirectory(0)
    return clone


def make_background_pdf(mass, category, model_name):
    """Return the selected positive background PDF and owned dependencies."""
    owned = []
    t = ROOT.RooFormulaVar(
        f"t_{category}",
        f"(@0-{mass.getMin():.17g})/{(mass.getMax()-mass.getMin()):.17g}",
        ROOT.RooArgList(mass),
    )
    x = ROOT.RooFormulaVar(
        f"x_{category}", "2*@0-1", ROOT.RooArgList(t)
    )
    owned.extend((t, x))

    if model_name == "SumExp":
        p1 = ROOT.RooRealVar(f"bkg_p1_{category}", "p1", -1.2, -8.0, 2.0)
        p2 = ROOT.RooRealVar(f"bkg_p2_{category}", "p2", -0.2, -8.0, 2.0)
        mix = ROOT.RooRealVar(f"bkg_mix_{category}", "logit mixture", 0.0, -8.0, 8.0)
        pdf = ROOT.RooGenericPdf(
            f"bkg_{category}",
            "(1/(1+exp(-@3)))*exp(@1*@0)+(1-1/(1+exp(-@3)))*exp(@2*@0)",
            ROOT.RooArgList(x, p1, p2, mix),
        )
        owned.extend((p1, p2, mix, pdf))
    elif model_name == "FEWZxBernstein":
        p0 = ROOT.RooRealVar(f"bkg_p0_{category}", "power", -1.0, -10.0, 4.0)
        b1 = ROOT.RooRealVar(f"bkg_b1_{category}", "log c1", 0.0, -5.0, 5.0)
        b2 = ROOT.RooRealVar(f"bkg_b2_{category}", "log c2", 0.0, -5.0, 5.0)
        pdf = ROOT.RooGenericPdf(
            f"bkg_{category}",
            "pow(@0/91.1876,@2)*((1-@1)*(1-@1)+2*exp(@3)*@1*(1-@1)+exp(@4)*@1*@1)",
            ROOT.RooArgList(mass, t, p0, b1, b2),
        )
        owned.extend((p0, b1, b2, pdf))
    elif model_name == "LandaUxBernstein":
        mpv = ROOT.RooRealVar(
            f"bkg_mpv_{category}",
            "Landau MPV",
            0.5 * (mass.getMin() + mass.getMax()),
            mass.getMin(),
            mass.getMax(),
        )
        # The three-parameter fit drives the Landau scale to its
        # upper boundary and cannot determine an upper likelihood crossing.
        # Fix it at the fitted boundary, leaving MPV and slope free.
        log_scale = ROOT.RooConstVar(
            f"bkg_logscale_{category}",
            "fixed log Landau scale",
            math.log(10.0),
        )
        b1 = ROOT.RooRealVar(
            f"bkg_b1_{category}",
            "exponential slope",
            0.0,
            -6.0,
            6.0,
        )

        # The positive floor prevents log(0) when the Landau tail
        # underflows numerically in the sidebands. This matches the
        # clipping used by the validated bias-toy implementation.
        pdf = ROOT.RooGenericPdf(
            f"bkg_{category}",
            (
                "(TMath::Landau(@0,@1,exp(@2),0)+1e-12)"
                "*exp(@3*@4)"
            ),
            ROOT.RooArgList(mass, mpv, log_scale, b1, x),
        )
        owned.extend((mpv, log_scale, b1, pdf))
    elif model_name == "PowerLaw":
        p0 = ROOT.RooRealVar(f"bkg_p0_{category}", "power", -1.5, -10.0, 4.0)
        p1 = ROOT.RooRealVar(f"bkg_p1_{category}", "linear", 0.0, -5.0, 5.0)
        p2 = ROOT.RooRealVar(f"bkg_p2_{category}", "quadratic", 0.0, -5.0, 5.0)
        pdf = ROOT.RooGenericPdf(
            f"bkg_{category}",
            f"pow(@0/{mass.getMin():.17g},@2)*exp(@3*@1+@4*@1*@1)",
            ROOT.RooArgList(mass, x, p0, p1, p2),
        )
        owned.extend((p0, p1, p2, pdf))
    else:
        raise ValueError(f"Unknown model {model_name}")
    return pdf, owned


def make_roodatahist(name, mass, hist):
    return ROOT.RooDataHist(name, name, ROOT.RooArgList(mass), hist)


def range_fraction(pdf, mass, range_name):
    integral = pdf.createIntegral(
        ROOT.RooArgSet(mass),
        ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)),
        ROOT.RooFit.Range(range_name),
    )
    return float(integral.getVal())


def make_hybrid_data(source_hist, pdf, mass, blind_min, blind_max, name):
    hybrid = source_hist.Clone(name + "_hist")
    hybrid.SetDirectory(0)
    sideband_total = float(source_hist.Integral())
    sb_fraction = range_fraction(pdf, mass, "lowSB") + range_fraction(pdf, mass, "highSB")
    if not math.isfinite(sb_fraction) or sb_fraction <= 0.0:
        raise RuntimeError("Invalid fitted sideband fraction")
    full_yield = sideband_total / sb_fraction

    norm_set = ROOT.RooArgSet(mass)
    for index in range(1, hybrid.GetNbinsX() + 1):
        low = hybrid.GetXaxis().GetBinLowEdge(index)
        high = hybrid.GetXaxis().GetBinUpEdge(index)
        if low >= blind_min and high <= blind_max:
            range_name = f"blindbin_{index}"
            mass.setRange(range_name, low, high)
            fraction = range_fraction(pdf, mass, range_name)
            expected = full_yield * fraction
            hybrid.SetBinContent(index, expected)
            hybrid.SetBinError(index, math.sqrt(max(expected, 0.0)))
    return make_roodatahist(name, mass, hybrid), full_yield


def validate_signal_systematics(root_file):
    missing = []
    for category in CATEGORIES:
        require_hist(root_file, f"{category}/VBF")
        for nuisance in SHAPE_NUISANCES:
            for direction in ("Up", "Down"):
                path = f"{category}/VBF_{nuisance}{direction}"
                if not root_file.Get(path):
                    missing.append(path)
    if missing:
        raise RuntimeError("Missing VBF systematic templates:\n  " + "\n  ".join(missing))


def write_datacard(path, workspace_name, rates, rate_errors, rate_scale):
    bins = list(CATEGORIES)
    lines = [
        "imax 3",
        "jmax 1",
        "kmax *",
        "------------",
    ]
    for category in bins:
        lines.append(
            f"shapes data_obs {category} {workspace_name} w:data_obs_{category}"
        )
        lines.append(
            f"shapes bkg {category} {workspace_name} w:bkg_{category}"
        )
        lines.append(
            f"shapes VBF {category} {workspace_name} "
            f"w:VBF_{category} w:VBF_{category}_$SYSTEMATIC"
        )
    lines.extend([
        "------------",
        "bin " + " ".join(bins),
        "observation " + " ".join("-1" for _ in bins),
        "------------",
        "bin " + " ".join(item for category in bins for item in (category, category)),
        "process " + " ".join(item for _ in bins for item in ("VBF", "bkg")),
        "process " + " ".join(item for _ in bins for item in ("0", "1")),
        "rate " + " ".join(
            item
            for category in bins
            for item in (f"{rates[category] * rate_scale:.10g}", "1")
        ),
        "------------",
    ])
    for nuisance in SHAPE_NUISANCES:
        lines.append(
            f"{nuisance:<14} shape " + " ".join(item for _ in bins for item in ("1", "-"))
        )
    # The analytic-background channels cannot use autoMCStats together with
    # the RooDataHist signal templates in current Combine.  The VBF samples
    # have O(2e5) effective events per category, so retain their small total
    # MC statistical uncertainty as independent category-level lnN terms.
    for target_category in bins:
        relative = (
            rate_errors[target_category] / abs(rates[target_category])
            if rates[target_category] != 0.0
            else 0.0
        )
        factor = 1.0 + relative
        effects = []
        for category in bins:
            effects.extend(
                (f"{factor:.8g}", "-")
                if category == target_category
                else ("-", "-")
            )
        lines.append(
            f"VBF_mcstat_{target_category:<4} lnN " + " ".join(effects)
        )
    lines.extend([
        "bkg_p1_cat0 flatParam",
        "bkg_p2_cat0 flatParam",
        "bkg_mix_cat0 flatParam",
        "bkg_cat0_norm flatParam",
        "bkg_mpv_cat1 flatParam",
        "bkg_b1_cat1 flatParam",
        "bkg_cat1_norm flatParam",
        "bkg_p1_cat2 flatParam",
        "bkg_p2_cat2 flatParam",
        "bkg_mix_cat2 flatParam",
        "bkg_cat2_norm flatParam",
    ])
    lines.extend([
        "# Add the collaboration-approved luminosity lnN and any remaining",
        "# normalization nuisances here before producing final results.",
        "# autoMCStats is intentionally disabled for this mixed model:",
        "# VBF is a RooDataHist template while bkg is an analytic RooAbsPdf.",
        "# VBF MC statistics are represented by the category-level lnN terms.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    templates = ROOT.TFile.Open(args.templates)
    if not templates or templates.IsZombie():
        raise RuntimeError(f"Cannot open {args.templates}")
    validate_signal_systematics(templates)

    workspace = ROOT.RooWorkspace("w", "2024 VBF BDT blinded workspace")
    mass = ROOT.RooRealVar("m_mumu", "m_{#mu#mu}", args.mass_min, args.mass_max, "GeV")
    reference_signal = require_hist(templates, "cat0/VBF")
    mass.setBins(reference_signal.GetNbinsX())
    mass.setRange("lowSB", args.mass_min, args.blind_min)
    mass.setRange("highSB", args.blind_max, args.mass_max)
    mass.setRange("full", args.mass_min, args.mass_max)
    ws_import(workspace, mass)

    keepalive = [mass]
    rates = {}
    rate_errors = {}
    summary = []
    for category in CATEGORIES:
        data_hist = require_hist(templates, f"{category}/data_obs_blinded")
        signal_hist = require_hist(templates, f"{category}/VBF")
        rates[category] = float(signal_hist.Integral())
        rate_errors[category] = math.sqrt(
            sum(
                float(signal_hist.GetBinError(index)) ** 2
                for index in range(1, signal_hist.GetNbinsX() + 1)
            )
        )
        sideband_data = make_roodatahist(f"data_sidebands_{category}", mass, data_hist)

        # Put the signal nominal and every systematic template in the same
        # RooFit observable space as the analytic background.  Keep these as
        # RooDataHist objects (rather than unit-normalized RooHistPdfs) so
        # Combine retains both the bin-shape change and the variation's total
        # normalization.  Importing them into this workspace also avoids the
        # axis mismatch seen when external TH1 templates were mixed with the
        # physical-mass background PDFs.
        signal_datahist = make_roodatahist(
            f"VBF_{category}", mass, signal_hist
        )
        ws_import(workspace, signal_datahist)
        keepalive.extend([signal_hist, signal_datahist])

        for nuisance in SHAPE_NUISANCES:
            for direction in ("Up", "Down"):
                variation = f"{nuisance}{direction}"
                variation_hist = require_hist(
                    templates, f"{category}/VBF_{variation}"
                )
                variation_datahist = make_roodatahist(
                    f"VBF_{category}_{variation}",
                    mass,
                    variation_hist,
                )
                ws_import(workspace, variation_datahist)
                keepalive.extend([variation_hist, variation_datahist])

        model_name = BACKGROUND_MODELS[category]
        pdf, owned = make_background_pdf(mass, category, model_name)
        fit_result = pdf.fitTo(
            sideband_data,
            ROOT.RooFit.Range("lowSB,highSB"),
            ROOT.RooFit.Save(True),
            ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.Strategy(1),
            ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
            ROOT.RooFit.Offset(True),
        )
        if int(fit_result.status()) != 0 or int(fit_result.covQual()) < 2:
            retry = pdf.fitTo(
                sideband_data,
                ROOT.RooFit.Range("lowSB,highSB"),
                ROOT.RooFit.Save(True),
                ROOT.RooFit.PrintLevel(-1),
                ROOT.RooFit.Strategy(2),
                ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                ROOT.RooFit.Offset(True),
                ROOT.RooFit.Hesse(True),
            )
            first_quality = (
                int(fit_result.status()) == 0,
                int(fit_result.covQual()),
                -float(fit_result.minNll()),
            )
            retry_quality = (
                int(retry.status()) == 0,
                int(retry.covQual()),
                -float(retry.minNll()),
            )
            if retry_quality > first_quality:
                fit_result = retry

        # Ensure the live PDF carries the parameters from the selected result
        # (the retry may have been evaluated but not selected).
        live_parameters = pdf.getParameters(sideband_data)
        final_parameters = fit_result.floatParsFinal()
        for parameter in final_parameters:
            live_parameter = live_parameters.find(parameter.GetName())
            if live_parameter:
                live_parameter.setVal(parameter.getVal())
                live_parameter.setError(parameter.getError())

        # Status 0 is ideal.  Status 1 is retained only with a good/full
        # covariance matrix; this is common for strongly correlated analytic
        # background parameters and remains suitable for the blinded pilot.
        if int(fit_result.status()) not in (0, 1) or int(fit_result.covQual()) < 2:
            raise RuntimeError(
                f"{category} {model_name} fit failed: status={fit_result.status()} "
                f"covQual={fit_result.covQual()}"
            )
        if int(fit_result.status()) != 0:
            print(
                f"WARNING: {category} {model_name} retained with "
                f"status={fit_result.status()} covQual={fit_result.covQual()} "
                f"edm={fit_result.edm():.6g}",
                flush=True,
            )

        data_obs, full_bkg = make_hybrid_data(
            data_hist, pdf, mass, args.blind_min, args.blind_max,
            f"data_obs_{category}",
        )
        bkg_norm = ROOT.RooRealVar(
            f"bkg_{category}_norm", "background normalization",
            full_bkg, max(1.0e-6, 0.05 * full_bkg), 20.0 * full_bkg,
        )

        ws_import(workspace, pdf)
        ws_import(workspace, bkg_norm)
        ws_import(workspace, data_obs)
        ws_import(workspace, sideband_data)
        # RooFitResult is a TObject rather than a RooAbsArg/RooAbsData.  It is
        # diagnostic-only and does not belong in the statistical workspace;
        # retain it in Python until the workspace has been written instead.
        keepalive.extend(owned + [pdf, fit_result, bkg_norm, data_obs, sideband_data])
        summary.append(
            f"{category}: model={model_name} sideband={data_hist.Integral():.0f} "
            f"full_bkg={full_bkg:.6g} signal={rates[category]:.6g} "
            f"fit_status={fit_result.status()} covQual={fit_result.covQual()}"
        )

    workspace_path = output_dir / f"{args.era}_vbf_bdt_workspace.root"
    output = ROOT.TFile(str(workspace_path), "RECREATE")
    workspace.Write()
    output.Close()
    templates.Close()

    datacard_path = output_dir / f"{args.era}_vbf_bdt_datacard.txt"
    write_datacard(
        datacard_path,
        workspace_path.name,
        rates,
        rate_errors,
        args.signal_rate_scale,
    )
    summary_path = output_dir / f"{args.era}_vbf_bdt_workspace_summary.txt"
    summary_path.write_text(
        "PRE-UNBLINDING HYBRID ASIMOV WORKSPACE\n"
        "Real data are used only outside [115,135) GeV.\n"
        "The blinded interval contains fitted-background Asimov values.\n\n"
        + "\n".join(summary) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote workspace: {workspace_path}")
    print(f"Wrote datacard:  {datacard_path}")
    print(f"Wrote summary:   {summary_path}")
    print("This is an expected/pre-unblinding workspace, not an observed-results workspace.")


if __name__ == "__main__":
    main()

