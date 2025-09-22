#!/usr/bin/env python3
import argparse, numpy as np, ROOT

ROOT.gStyle.SetOptStat(0)

def parse():
    p = argparse.ArgumentParser(description="Plot bdt_score for signal/background and find an optimal cut.")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--combined", nargs=3, metavar=("SCORED_ALL.root","TEST_SIG.root","TEST_BKG.root"),
                   help="Use one scored file created from [test_signal, test_background] in that order, "
                        "plus the original test files to determine the split point.")
    g.add_argument("--separate", nargs=2, metavar=("SCORED_SIG.root","SCORED_BKG.root"),
                   help="Use two scored files (signal and background) produced separately.")
    p.add_argument("--treename", default="mva")
    p.add_argument("--score", default="bdt_score")
    p.add_argument("--weight", default="weight")
    p.add_argument("--bins", type=int, default=60)
    p.add_argument("--out", default="tmva104_bdt_score")
    return p.parse_args()

def arrays_from_file(fname, treename, branches):
    df = ROOT.RDataFrame(treename, fname)
    cols = df.AsNumpy(branches)
    return [np.asarray(cols[b]) for b in branches]

def get_sig_bkg_arrays(args):
    b = [args.score, args.weight]
    if args.combined:
        scored_all, test_sig, test_bkg = args.combined
        nsig = int(ROOT.TChain(args.treename).Add(test_sig) or ROOT.TChain(args.treename).GetEntries())
        # A trick to get entries without rebuilding chains twice
        ch_s = ROOT.TChain(args.treename); ch_s.Add(test_sig); nsig = ch_s.GetEntries()
        ch_b = ROOT.TChain(args.treename); ch_b.Add(test_bkg); nbkg = ch_b.GetEntries()

        s_all, w_all = arrays_from_file(scored_all, args.treename, b)
        assert s_all.size == nsig + nbkg, "Entry mismatch with combined file."

        s_sig, w_sig = s_all[:nsig], w_all[:nsig]
        s_bkg, w_bkg = s_all[nsig:], w_all[nsig:]
        return (s_sig, w_sig), (s_bkg, w_bkg)

    else:
        scored_sig, scored_bkg = args.separate
        s_sig, w_sig = arrays_from_file(scored_sig, args.treename, b)
        s_bkg, w_bkg = arrays_from_file(scored_bkg, args.treename, b)
        return (s_sig, w_sig), (s_bkg, w_bkg)

def effs_vs_cut(scores, weights, cuts):
    # Efficiency for scores >= cut
    order = np.argsort(scores)
    s = scores[order]; w = weights[order]
    csum = np.cumsum(w)                        # increasing with score
    total = csum[-1] if csum.size else 1.0
    # map each cut to index where score >= cut
    idx = np.searchsorted(s, cuts, side="left")
    kept = (total - np.concatenate(([0.0], csum)))[idx] if idx.size else np.array([])
    return kept / max(total, 1e-9)

if __name__ == "__main__":
    args = parse()
    (s_sig, w_sig), (s_bkg, w_bkg) = get_sig_bkg_arrays(args)

    # --- histograms (normalized shapes) ---
    xmin = 0.0; xmax = 1.0
    nb   = args.bins
    hS = ROOT.TH1F("hS",";BDT score;Normalized entries", nb, xmin, xmax)
    hB = ROOT.TH1F("hB",";BDT score;Normalized entries", nb, xmin, xmax)
    for x, w in zip(s_sig, w_sig): hS.Fill(float(x), float(abs(w)))
    for x, w in zip(s_bkg, w_bkg): hB.Fill(float(x), float(abs(w)))
    if hS.Integral()>0: hS.Scale(1.0/hS.Integral())
    if hB.Integral()>0: hB.Scale(1.0/hB.Integral())

    hS.SetLineColor(ROOT.kRed);   hS.SetLineWidth(3)
    hB.SetLineColor(ROOT.kBlue);  hB.SetLineWidth(3)

    c1 = ROOT.TCanvas("c1","",750,650)
    hS.Draw("hist"); hB.Draw("hist same")
    leg = ROOT.TLegend(0.15,0.78,0.45,0.90); leg.SetBorderSize(0)
    leg.AddEntry(hS,"Signal (shape)","l")
    leg.AddEntry(hB,"Background (shape)","l")
    leg.Draw()
    c1.Update(); c1.SaveAs(args.out+"_shapes.pdf"); c1.SaveAs(args.out+"_shapes.png")

    # --- find optimal cut on test set (maximize S/sqrt(S+B)) ---
    # Use absolute weights for training-style figure of merit
    cuts = np.linspace(xmin, xmax, 400)
    epsS = effs_vs_cut(s_sig, np.abs(w_sig), cuts)
    epsB = effs_vs_cut(s_bkg, np.abs(w_bkg), cuts)
    S0 = np.sum(np.abs(w_sig)); B0 = np.sum(np.abs(w_bkg))
    S = S0 * epsS; B = B0 * epsB
    Z = S / np.sqrt(S + B + 1e-9)
    iopt = int(np.argmax(Z))
    best_cut = float(cuts[iopt]); best_Z = float(Z[iopt])
    print(f"[tmva104] Best cut (maximize S/sqrt(S+B)) = {best_cut:.4f}  →  Z = {best_Z:.3f}")
    print(f"[tmva104] Efficiencies at best cut:  epsS = {epsS[iopt]:.3f},  epsB = {epsB[iopt]:.3f}")

    # Plot Z(cut)
    c2 = ROOT.TCanvas("c2","",750,650)
    grZ = ROOT.TGraph(len(cuts), np.array(cuts, dtype="f8"), np.array(Z, dtype="f8"))
    grZ.SetTitle(";BDT score cut; S / #sqrt{S+B}")
    grZ.SetLineWidth(3); grZ.SetLineColor(ROOT.kBlack)
    grZ.Draw("AL")
    line = ROOT.TLine(best_cut, 0.0, best_cut, best_Z); line.SetLineColor(ROOT.kRed); line.SetLineStyle(2); line.Draw()
    lab  = ROOT.TLatex(); lab.SetNDC(True); lab.SetTextSize(0.035)
    lab.DrawLatex(0.15,0.92, f"Best cut = {best_cut:.3f},  Z = {best_Z:.2f}")
    c2.Update(); c2.SaveAs(args.out+"_significance.pdf"); c2.SaveAs(args.out+"_significance.png")

    print(f"[tmva104] Wrote {args.out}_shapes.(pdf|png) and {args.out}_significance.(pdf|png)")

