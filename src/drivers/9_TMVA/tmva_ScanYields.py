#!/usr/bin/env python3
# Scan BDT score cuts over many ROOT files and compute S, B, Z, ZA.
# Works directly from your TMVA RBDT model and "mva" trees.

import argparse, os, glob, numpy as np, ROOT

# ----------------- config: must match training order -----------------
from tmva100_DataPreparation import VARIABLES
TREENAME = "mva"
WEIGHT_BRANCH = "weight"     # your physics weight (already scaled to lumi)
# --------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Scan BDT cuts and compute yields/significance.")
    p.add_argument("--model-file", default="tmva101.root", help="TMVA model file (from training)")
    p.add_argument("--model-name", default="myBDT", help="Model name used in SaveXGBoost")
    p.add_argument("--sig", nargs="+", required=True,
                   help="Signal inputs (files or directories). You can pass many.")
    p.add_argument("--bkg", nargs="+", required=True,
                   help="Background inputs (files or directories). You can pass many.")
    p.add_argument("--treename", default=TREENAME)
    p.add_argument("--weight",   default=WEIGHT_BRANCH)
    p.add_argument("--cuts",     default="0.00:1.00:401",
                   help="scan range as start:end:npoints, e.g. 0.00:1.00:401 (inclusive grid)")
    p.add_argument("--out-prefix", default="tmva_fullRun",
                   help="prefix for CSV and plots")
    p.add_argument("--abs-weights", action="store_true",
                   help="Use |weight| for S,B (useful for quick sanity if negative weights are present).")
    return p.parse_args()

def expand_inputs(items):
    files = []
    for it in items:
        if os.path.isdir(it):
            files += sorted(glob.glob(os.path.join(it, "*.root")))
        elif it.endswith(".root") and os.path.isfile(it):
            files.append(it)
        else:
            # allow globs directly
            files += sorted(glob.glob(it))
    return files

def load_arrays(file_list, treename, vars_needed, weight_name):
    flist = ROOT.std.vector("string")()
    for f in file_list: flist.push_back(f)
    if len(flist) == 0:
        return (np.empty((0,), dtype=np.float32),) * (len(vars_needed) + 1)

    df = ROOT.RDataFrame(treename, flist)
    cols = set(df.GetColumnNames())
    missing = [v for v in vars_needed if v not in cols]
    if missing:
        raise RuntimeError("Missing variables in input tree: " + ", ".join(missing))
    if weight_name not in cols:
        raise RuntimeError("Weight branch '{}' not found".format(weight_name))

    arrays = df.AsNumpy(vars_needed + [weight_name])
    xs = [np.asarray(arrays[v], dtype=np.float32) for v in vars_needed]
    w  = np.asarray(arrays[weight_name], dtype=np.float64)  # keep precision
    return xs + [w]

def compute_scores(model, xs):
    X = np.vstack(xs).T  # shape (N, nvars)
    return model.Compute(X.astype(np.float32, copy=False)).astype(np.float64, copy=False)

def yields_above_cut(scores, weights, cut, use_abs=False):
    if use_abs:
        weights = np.abs(weights)
    mask = scores >= cut
    return float(weights[mask].sum())

def ZA_asymptotic(S, B):
    if B <= 0 or S <= 0:
        return 0.0
    return float(np.sqrt(2.0*((S+B)*np.log(1.0 + S/max(B,1e-12)) - S)))

def parse_cuts(spec):
    a, b, n = spec.split(":")
    a, b, n = float(a), float(b), int(n)
    return np.linspace(a, b, n)

if __name__ == "__main__":
    args = parse_args()

    # Ensure TMVA is loaded and model is accessible
    ROOT.gSystem.Load("libTMVA")
    if ROOT.gSystem.AccessPathName(args.model_file, ROOT.kReadPermission):
        raise SystemExit("Model file not readable: {}".format(args.model_file))
    model = ROOT.TMVA.Experimental.RBDT[""](args.model_name, args.model_file)

    # Gather files
    sig_files = expand_inputs(args.sig)
    bkg_files = expand_inputs(args.bkg)
    if not sig_files: raise SystemExit("No signal files found.")
    if not bkg_files: raise SystemExit("No background files found.")

    print("[info] signal files:", len(sig_files))
    print("[info] background files:", len(bkg_files))

    # Load arrays and compute scores
    *xs_s, w_s = load_arrays(sig_files, args.treename, VARIABLES, args.weight)
    *xs_b, w_b = load_arrays(bkg_files, args.treename, VARIABLES, args.weight)
    print("[info] entries: S={}, B={}".format(w_s.size, w_b.size))

    s_s = compute_scores(model, xs_s)
    s_b = compute_scores(model, xs_b)

    s_s = np.ravel(s_s)      # make (N,) instead of (N,1)
    s_b = np.ravel(s_b)
    w_s = np.ravel(w_s)
    w_b = np.ravel(w_b)

    # scan
    cuts = parse_cuts(args.cuts)
    rows = []
    best_Z = (-1.0, None)
    best_ZA = (-1.0, None)
    for c in cuts:
        S = yields_above_cut(s_s, w_s, c, use_abs=args.abs_weights)
        B = yields_above_cut(s_b, w_b, c, use_abs=args.abs_weights)
        Z  = S/np.sqrt(S+B+1e-12) if (S+B)>0 else 0.0
        ZA = ZA_asymptotic(S, B)
        rows.append((c, S, B, Z, ZA))
        if Z  > best_Z[0]:   best_Z  = (Z,  c)
        if ZA > best_ZA[0]:  best_ZA = (ZA, c)

    # print summary
    print("\n=== Best cuts ===")
    print(" by  Z  (S/sqrt(S+B)):  cut = {:.4f}  Z = {:.3f}".format(best_Z[1], best_Z[0]))
    print(" by ZA (asymptotic)  :  cut = {:.4f}  ZA= {:.3f}".format(best_ZA[1], best_ZA[0]))

    # write CSV
    csv = args.out_prefix + "_yields.csv"
    with open(csv, "w") as f:
        f.write("cut,S,B,Z,ZA\n")
        for c,S,B,Z,ZA in rows:
            f.write("{:.6f},{:.6f},{:.6f},{:.6f},{:.6f}\n".format(c,S,B,Z,ZA))
    print("[out] wrote", csv)

    # quick plots
    import ROOT
    xs = np.array([r[0] for r in rows], dtype="f8")
    Zs = np.array([r[3] for r in rows], dtype="f8")
    ZAs= np.array([r[4] for r in rows], dtype="f8")

    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1","",750,650)
    gZ = ROOT.TGraph(len(xs), xs, Zs);   gZ.SetLineWidth(3); gZ.SetTitle(";BDT score cut; S / #sqrt{S+B}")
    gZ.Draw("AL")
    l1 = ROOT.TLine(best_Z[1], 0.0, best_Z[1], best_Z[0]); l1.SetLineColor(ROOT.kRed); l1.SetLineStyle(2); l1.Draw()
    lab = ROOT.TLatex(); lab.SetNDC(True); lab.SetTextSize(0.035)
    lab.DrawLatex(0.15,0.92, "Best Z cut = {:.3f}, Z = {:.2f}".format(best_Z[1], best_Z[0]))
    c1.SaveAs(args.out_prefix + "_Z.pdf"); c1.SaveAs(args.out_prefix + "_Z.png")

    c2 = ROOT.TCanvas("c2","",750,650)
    gZA= ROOT.TGraph(len(xs), xs, ZAs);  gZA.SetLineWidth(3); gZA.SetTitle(";BDT score cut; Z_{A}")
    gZA.Draw("AL")
    l2 = ROOT.TLine(best_ZA[1], 0.0, best_ZA[1], best_ZA[0]); l2.SetLineColor(ROOT.kRed); l2.SetLineStyle(2); l2.Draw()
    lab.DrawLatex(0.15,0.92, "Best ZA cut = {:.3f}, ZA = {:.2f}".format(best_ZA[1], best_ZA[0]))
    c2.SaveAs(args.out_prefix + "_ZA.pdf"); c2.SaveAs(args.out_prefix + "_ZA.png")

    print("[out] wrote", args.out_prefix + "_Z.(pdf|png)", "and", args.out_prefix + "_ZA.(pdf|png)")

