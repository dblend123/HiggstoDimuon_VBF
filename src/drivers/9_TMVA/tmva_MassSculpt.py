#!/usr/bin/env python3
# Mass sculpting check: plot m (e.g. dimuon mass) shape in BDT score slices
# Writes plots to Outputs/ and prints simple correlations.
# python3 tmva_MassSculpt.py   --input scored_all.root   --treename mva   --mass-branch m   --mass-min 115 --mass-max 135   --bdt-bins 0.00,0.15,0.30,0.45,0.60,0.75,1.00   --outdir Outputs

import os, argparse, numpy as np, ROOT

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, 
        help="ROOT file with bdt_score branch (e.g. scored_all.root)")
    ap.add_argument("--treename", default="mva")
    ap.add_argument("--mass-branch", default="m", 
        help="Mass variable branch name (e.g. m or m_mumu)")
    ap.add_argument("--is-signal", action="store_true",
        help="If set, plot signal; otherwise defaults to background")
    ap.add_argument("--weight", default="", 
        help="Optional event weight branch name")
    ap.add_argument("--mass-min", type=float, default=110.0)
    ap.add_argument("--mass-max", type=float, default=150.0)
    ap.add_argument("--nbins-mass", type=int, default=40)
    ap.add_argument("--bdt-bins", default="0.00,0.15,0.30,0.45,0.60,0.75,1.00",
        help="Comma-separated BDT edges (inclusive low, exclusive high)")
    ap.add_argument("--outdir", default="Outputs")
    ap.add_argument("--year", default="22023", 
        help="Year tag for output filename")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        raise SystemExit(f"[error] cannot open {args.input}")
    t = f.Get(args.treename)
    if not t:
        raise SystemExit(f"[error] tree '{args.treename}' not found in {args.input}")

    cols = set([br.GetName() for br in t.GetListOfBranches()])
    need = {args.mass_branch, "bdt_score"}
    if not need.issubset(cols):
        raise SystemExit(f"[error] missing required branches: {need - cols}")
    has_y = "y" in cols
    has_w = (args.weight in cols) if args.weight else False

    # Select class (signal vs background) if label exists; otherwise use all entries
    # Expect label 1 for signal, 0 for background (same convention as training).
    sel_class = ""
    if has_y:
        sel_class = "(y==1)" if args.is_signal else "(y==0)"

    # Mass window and finite filters
    mass_cut = f"({args.mass_branch}>={args.mass_min} && {args.mass_branch}<={args.mass_max})"
    finite_cut = f"finite({args.mass_branch}) && finite(bdt_score)"
    base_cut = " && ".join([c for c in [sel_class, mass_cut, finite_cut] if c])

    # Build RDataFrame and read arrays
    df = ROOT.RDataFrame(args.treename, args.input).Filter(base_cut)
    if has_w:
        df = df.Define("evt_w", args.weight)
    else:
        df = df.Define("evt_w", "1.0")

    # Grab numpy arrays
    m_vals = np.array(df.AsNumpy([args.mass_branch])[args.mass_branch], dtype=np.float64)
    s_vals = np.array(df.AsNumpy(["bdt_score"])["bdt_score"], dtype=np.float64)
    w_vals = np.array(df.AsNumpy(["evt_w"])["evt_w"], dtype=np.float64)

    # Quick diagnostics: Pearson corr(m, bdt) in window
    def wmean(a, w): return np.average(a, weights=w) if a.size else np.nan
    def wcov(x, y, w):
        mu_x, mu_y = wmean(x, w), wmean(y, w)
        return np.average((x-mu_x)*(y-mu_y), weights=w)
    def wvar(x, w): return wcov(x, x, w)

    if m_vals.size:
        cov = wcov(m_vals, s_vals, w_vals)
        r = cov / np.sqrt(wvar(m_vals, w_vals) * wvar(s_vals, w_vals))
        print(f"[info] Weighted Pearson corr(mass, bdt) in [{args.mass_min},{args.mass_max}] = {r:.4f}")

    # Prepare BDT slices
    edges = [float(x) for x in args.bdt_bins.split(",")]
    edges = sorted(edges)
    labels = [f"{edges[i]:.2f} ≤ BDT < {edges[i+1]:.2f}" for i in range(len(edges)-1)]

    # Make stacked CDF-like overlay plot of normalized mass shapes per BDT slice
    c = ROOT.TCanvas("c_mass", "", 800, 650)
    leg = ROOT.TLegend(0.62, 0.62, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    # Draw empty frame
    frame = ROOT.TH1F("frame","", args.nbins_mass, args.mass_min, args.mass_max)
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    frame.GetYaxis().SetTitle("A.U.")
    frame.GetYaxis().SetRangeUser(0, 1.05)
    frame.Draw("AXIS")

    # Color palette
    colors = [ROOT.kAzure+2, ROOT.kCyan+2, ROOT.kTeal+4, ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kRed+1]

    hists = []
    for i in range(len(edges)-1):
        lo, hi = edges[i], edges[i+1]
        mask = (s_vals >= lo) & (s_vals < hi)
        if not np.any(mask):
            continue
        m_i = m_vals[mask]; w_i = w_vals[mask]
        h = ROOT.TH1F(f"h{i}", "", args.nbins_mass, args.mass_min, args.mass_max)
        # fill with weights
        for mv, wv in zip(m_i, w_i):
            h.Fill(mv, wv)
        # normalize to unit area (shape only)
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0/integral)
        h.SetLineWidth(3)
        h.SetLineColor(colors[i % len(colors)])
        hists.append((h, labels[i]))

    # auto y-range from max of shapes
    ymax = 0.0
    for h,_ in hists:
        ymax = max(ymax, h.GetMaximum())
    frame.GetYaxis().SetRangeUser(0, ymax*1.25 if ymax>0 else 0.05)
    frame.Draw("AXIS")
    for h, lab in hists:
        h.Draw("HIST SAME")
        leg.AddEntry(h, lab, "l")
    leg.Draw()

    what = "signal" if args.is_signal else "background"
    title = ROOT.TLatex()
    title.SetNDC(True); title.SetTextSize(0.04)
    title.DrawLatex(0.15, 0.93,
        f"Mass sculpting check ({what}), {args.mass_min:.0f} < m_{{#mu#mu}} < {args.mass_max:.0f} GeV")

    out_png = os.path.join(args.outdir, f"{args.year}_tmva_mass_sculpt_{what}.png")
    out_pdf = os.path.join(args.outdir, f"{args.year}_tmva_mass_sculpt_{what}.pdf")
    c.SaveAs(out_pdf); c.SaveAs(out_png)
    print(f"[info] wrote {out_pdf} and {out_png}")

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    main()

