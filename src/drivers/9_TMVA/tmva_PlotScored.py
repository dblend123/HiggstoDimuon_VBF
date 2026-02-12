#!/usr/bin/env python3
# Auto-plot BDT score shapes and significance per era + combined.
# Supports either:
#   (A) discovering test files in a directory (default behavior), or
#   (B) explicit --signal and --backgrounds file paths.
#
# Outputs:
#   <out_prefix>_shapes.(pdf|png)
#   <out_prefix>_Z.(pdf|png)
#   <out_prefix>_ZA.(pdf|png)

import os, glob, numpy as np, ROOT
import argparse

# ----------------- CONFIG DEFAULTS (match training) -----------------
from tmva100_DataPreparation import VARIABLES

DEFAULT_TREENAME       = "mva"
DEFAULT_WEIGHT_BRANCH  = "weight"

DEFAULT_PREP_DIR       = "DataPreparation_Files"   # where *_test_{signal,background}.root live
DEFAULT_MODEL_FILE     = "tmva101.root"
DEFAULT_MODEL_NAME     = "myBDT"

NBINS   = 60
XRANGE  = (0.0, 1.0)
# -------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot BDT score shapes and significance; discover files by era or take explicit inputs."
    )
    # Discovery-mode directory
    p.add_argument("--prep-dir", default=DEFAULT_PREP_DIR,
                   help=f"Directory containing *_test_signal.root and *_test_background.root pairs (default: {DEFAULT_PREP_DIR})")
    # Explicit files (override discovery mode if provided)
    p.add_argument("--signal", help="Explicit path to a *_test_signal.root file")
    p.add_argument("--backgrounds", nargs="+",
                   help="One or more explicit background ROOT files (merged logically)")

    # Model / tree / weight
    p.add_argument("--model-file", default=DEFAULT_MODEL_FILE,
                   help=f"TMVA model file (default: {DEFAULT_MODEL_FILE})")
    p.add_argument("--model-name", default=DEFAULT_MODEL_NAME,
                   help=f"TMVA method name inside the model file (default: {DEFAULT_MODEL_NAME})")
    p.add_argument("--treename", default=DEFAULT_TREENAME,
                   help=f"TTree name (default: {DEFAULT_TREENAME})")
    p.add_argument("--weight-branch", default=DEFAULT_WEIGHT_BRANCH,
                   help=f"Weight branch name (default: {DEFAULT_WEIGHT_BRANCH})")

    # Output
    p.add_argument("--out-prefix", default="ALL_tmva",
                   help="Prefix for output files (default: ALL_tmva)")
    return p.parse_args()

def load_test_pairs(prep_dir):
    # Returns pairs = [(era, sig_file, bkg_file), ...]
    pairs = []
    for sig in sorted(glob.glob(os.path.join(prep_dir, "*_test_signal.root"))):
        era = os.path.basename(sig).replace("_test_signal.root", "")
        bkg = os.path.join(prep_dir, f"{era}_test_background.root")
        if os.path.isfile(bkg):
            pairs.append((era, sig, bkg))
    if not pairs:
        raise FileNotFoundError(f"No test pairs found in {prep_dir}")
    return pairs

def arrays_from_file(fname, treename, branches):
    df = ROOT.RDataFrame(treename, fname)
    cols = df.AsNumpy(branches)
    return [np.asarray(cols[b]) for b in branches]

def sanitize_weights(w):
    w = np.asarray(w, dtype=np.float64)
    w[~np.isfinite(w)] = 0.0
    w = np.abs(w)  # for plotting & quick significance
    w[w <= 1e-9] = 1e-9
    return w

def compute_scores(model, xs):
    X = np.vstack(xs).T.astype(np.float32, copy=False)
    return np.ravel(model.Compute(X).astype(np.float64, copy=False))

def effs_vs_cut(scores, weights, cuts):
    # efficiency for score >= cut
    order = np.argsort(scores)
    s = scores[order]; w = weights[order]
    csum = np.cumsum(w)                  # cumulative sum from low to high
    total = csum[-1] if csum.size else 1.0
    idx = np.searchsorted(s, cuts, side="left")
    kept = total - np.concatenate(([0.0], csum))[idx]
    return kept / max(total, 1e-12)

def ZA_asymptotic(S, B):
    if S <= 0.0 or B <= 0.0: return 0.0
    return float(np.sqrt(2.0 * ((S+B) * np.log(1.0 + S/max(B,1e-12)) - S)))

def plot_shapes(sig_scores, sig_w, bkg_scores, bkg_w, outprefix):
    xmin, xmax = XRANGE
    hS = ROOT.TH1F("hS",";BDT score;Normalized entries", NBINS, xmin, xmax)
    hB = ROOT.TH1F("hB",";BDT score;Normalized entries", NBINS, xmin, xmax)
    for x, w in zip(sig_scores, sig_w): hS.Fill(float(x), float(w))
    for x, w in zip(bkg_scores, bkg_w): hB.Fill(float(x), float(w))
    if hS.Integral()>0: hS.Scale(1.0/hS.Integral())
    if hB.Integral()>0: hB.Scale(1.0/hB.Integral())
    hS.SetLineColor(ROOT.kRed);   hS.SetLineWidth(3)
    hB.SetLineColor(ROOT.kBlue);  hB.SetLineWidth(3)
    c = ROOT.TCanvas("c_shapes","",750,650)
    hS.Draw("hist"); hB.Draw("hist same")
    leg = ROOT.TLegend(0.15,0.78,0.45,0.90); leg.SetBorderSize(0)
    leg.AddEntry(hS,"Signal (shape)","l")
    leg.AddEntry(hB,"Background (shape)","l")
    leg.Draw()
    c.SaveAs(outprefix + "_shapes.pdf"); c.SaveAs(outprefix + "_shapes.png")

def plot_significance(sig_scores, sig_w, bkg_scores, bkg_w, outprefix):
    xmin, xmax = XRANGE
    cuts = np.linspace(xmin, xmax, 400)
    epsS = effs_vs_cut(sig_scores, sig_w, cuts)
    epsB = effs_vs_cut(bkg_scores, bkg_w, cuts)
    S0 = float(sig_w.sum()); B0 = float(bkg_w.sum())
    S = S0 * epsS; B = B0 * epsB
    Z  = np.where(S+B>0, S/np.sqrt(S+B), 0.0)
    ZA = np.array([ZA_asymptotic(s, b) for s, b in zip(S, B)], dtype="f8")

    iZ  = int(np.argmax(Z));  cZ  = cuts[iZ];  Zbest  = Z[iZ]
    iZA = int(np.argmax(ZA)); cZA = cuts[iZA]; ZAbest = ZA[iZA]

    # Z plot
    c1 = ROOT.TCanvas("cZ","",750,650)
    gZ = ROOT.TGraph(len(cuts), cuts.astype("f8"), Z.astype("f8"))
    gZ.SetTitle(";BDT score cut; S / #sqrt{S+B}")
    gZ.SetLineWidth(3); gZ.Draw("AL")
    l1 = ROOT.TLine(float(cZ), 0.0, float(cZ), float(Zbest)); l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.Draw()
    lab = ROOT.TLatex(); lab.SetNDC(True); lab.SetTextSize(0.035)
    lab.DrawLatex(0.15,0.92, f"Best Z cut = {cZ:.3f}, Z = {Zbest:.2f}")
    c1.SaveAs(outprefix + "_Z.pdf"); c1.SaveAs(outprefix + "_Z.png")

    # ZA plot
    c2 = ROOT.TCanvas("cZA","",750,650)
    gZA = ROOT.TGraph(len(cuts), cuts.astype("f8"), ZA.astype("f8"))
    gZA.SetTitle(";BDT score cut; Z_{A}")
    gZA.SetLineWidth(3); gZA.Draw("AL")
    l2 = ROOT.TLine(float(cZA), 0.0, float(cZA), float(ZAbest)); l2.SetLineColor(ROOT.kRed); l2.SetLineStyle(2); l2.Draw()
    lab.DrawLatex(0.15,0.92, f"Best ZA cut = {cZA:.3f}, ZA = {ZAbest:.2f}")
    c2.SaveAs(outprefix + "_ZA.pdf"); c2.SaveAs(outprefix + "_ZA.png")

    print(f"[{outprefix}] S0={S0:.2f}, B0={B0:.2f} | best Z: cut={cZ:.3f}, Z={Zbest:.2f} | best ZA: cut={cZA:.3f}, ZA={ZAbest:.2f}")

def discover_or_use_explicit(args):
    """
    Returns:
      - If explicit files provided: ([(era, sig, bkg)], out_prefix)
      - Else discovery pairs from prep-dir: ([(era, sig, bkg), ...], out_prefix)
    """
    if args.signal or args.backgrounds:
        if not (args.signal and args.backgrounds):
            raise SystemExit("If you use explicit file mode, provide BOTH --signal and --backgrounds.")
        # Use a virtual 'era' label for explicit mode
        era = "EXPLICIT"
        # If multiple backgrounds are passed, we will read them all and concatenate arrays later
        return ([(era, args.signal, args.backgrounds)], args.out_prefix)

    # Discovery mode: find *_test_signal.root and *_test_background.root in prep-dir
    pairs = load_test_pairs(args.prep_dir)
    return (pairs, args.out_prefix)

if __name__ == "__main__":
    args = parse_args()

    ROOT.gSystem.Load("libTMVA")
    if ROOT.gSystem.AccessPathName(args.model_file, ROOT.kReadPermission):
        raise SystemExit(f"Model file not readable: {args.model_file}")

    # Load the trained TMVA model
    bdt = ROOT.TMVA.Experimental.RBDT[""](args.model_name, args.model_file)

    pairs, outprefix = discover_or_use_explicit(args)

    # Accumulate ALL eras (or explicit set)
    all_sig_scores, all_sig_w = [], []
    all_bkg_scores, all_bkg_w = [], []

    for era, sig, bkg in pairs:
        # ----- signal -----
        sig_arrays = arrays_from_file(sig, args.treename, VARIABLES + [args.weight_branch])
        xs_sig, w_sig = sig_arrays[:-1], sanitize_weights(sig_arrays[-1])
        s_sig = compute_scores(bdt, xs_sig)
        all_sig_scores.append(s_sig); all_sig_w.append(w_sig)

        # ----- background(s) -----
        if isinstance(bkg, (list, tuple)):  # explicit multi-file background
            era_bkg_scores, era_bkg_w = [], []
            for bfile in bkg:
                bkg_arrays = arrays_from_file(bfile, args.treename, VARIABLES + [args.weight_branch])
                xs_bkg, w_bkg = bkg_arrays[:-1], sanitize_weights(bkg_arrays[-1])
                s_bkg = compute_scores(bdt, xs_bkg)
                era_bkg_scores.append(s_bkg); era_bkg_w.append(w_bkg)
            all_bkg_scores.append(np.concatenate(era_bkg_scores))
            all_bkg_w.append(np.concatenate(era_bkg_w))
        else:
            bkg_arrays = arrays_from_file(bkg, args.treename, VARIABLES + [args.weight_branch])
            xs_bkg, w_bkg = bkg_arrays[:-1], sanitize_weights(bkg_arrays[-1])
            s_bkg = compute_scores(bdt, xs_bkg)
            all_bkg_scores.append(s_bkg); all_bkg_w.append(w_bkg)

    # Single combined arrays
    sS = np.concatenate(all_sig_scores); wS = np.concatenate(all_sig_w)
    sB = np.concatenate(all_bkg_scores); wB = np.concatenate(all_bkg_w)

    plot_shapes(sS, wS, sB, wB, outprefix)
    plot_significance(sS, wS, sB, wB, outprefix)

    print(f"[done] Wrote {outprefix}_shapes.(pdf|png) and {outprefix}_{{Z,ZA}}.(pdf|png)")