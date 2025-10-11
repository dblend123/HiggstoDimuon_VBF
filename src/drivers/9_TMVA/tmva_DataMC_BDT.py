
#!/usr/bin/env python3
import os, argparse, ROOT
ROOT.gROOT.SetBatch(True)

PROC_ORDER = ["DY", "EWK", "Top", "Diboson", "Other", "ggH", "VBF"]
PROC_COLOR = {
    "DY":       ROOT.kOrange+7,
    "EWK":      ROOT.kViolet-5,
    "Top":      ROOT.kAzure+1,
    "Diboson":  ROOT.kTeal+4,
    "Other":    ROOT.kPink+6,
    "ggH":      ROOT.kGray+1,
    "VBF":      ROOT.kGreen+2,
}
def canonical_proc(name: str) -> str:
    s = name.lower()
    if "dy" in s or "z" in s:           return "DY"
    if "ewk" in s or "z+ewk" in s:       return "EWK"
    if "top" in s or "tt" in s:          return "Top"
    if "vv" in s or "dibos" in s:        return "Diboson"
    if "ggh" in s:                       return "ggH"
    if "vbf" in s:                       return "VBF"
    return "Other"
def proc_sort_key(name: str) -> int:
    cn = canonical_proc(name);  return PROC_ORDER.index(cn) if cn in PROC_ORDER else len(PROC_ORDER)
def proc_color(name: str) -> int:
    cn = canonical_proc(name);  return PROC_COLOR.get(cn, ROOT.kGray+2)
def parse_group_map(items, flagname):
    out = {}
    for s in items:
        if ":" not in s:
            raise SystemExit(f"[error] {flagname} expects key:val, got '{s}'")
        k, v = s.split(":", 1)
        out.setdefault(k, []).append(v)
    return out
def sum_h(hlist):
    if not hlist: return None
    hsum = hlist[0].Clone(hlist[0].GetName()+"_sum")
    for h in hlist[1:]: hsum.Add(h)
    return hsum
def hist_data(files, treename, bdt_branch, sel, nbins, xmin, xmax):
    hlist = []
    for path in files:
        if not os.path.exists(path): raise SystemExit(f"[error] data file not found: {path}")
        df = ROOT.RDataFrame(treename, path)
        if sel: df = df.Filter(sel)
        h = df.Histo1D((f"hData_{os.path.basename(path)}","",nbins,xmin,xmax), bdt_branch).GetPtr().Clone()
        h.Sumw2(); hlist.append(h)
    h = sum_h(hlist)
    if not h: raise SystemExit("[error] empty data histogram")
    h.SetMarkerStyle(20); h.SetLineColor(ROOT.kBlack); h.SetMarkerSize(0.9)
    return h
def hist_mc_group(name, files, treename, bdt_branch, wbranch, sel, nbins, xmin, xmax, color):
    hlist = []
    for path in files:
        if not os.path.exists(path): raise SystemExit(f"[error] mc file not found: {path}")
        df = ROOT.RDataFrame(treename, path)
        if sel: df = df.Filter(sel)
        h = df.Histo1D((f"hMC_{name}_{os.path.basename(path)}","",nbins,xmin,xmax), bdt_branch, wbranch).GetPtr().Clone()
        h.Sumw2(); hlist.append(h)
    h = sum_h(hlist)
    if not h: return None
    h.SetFillColor(color); h.SetLineColor(ROOT.kBlack); h.SetLineWidth(1); h.SetFillStyle(1001)
    return h
def make_panel(era, h_data, group_hists, nbins, xmin, xmax, title_top, logy, outdir):
    h_mc_tot = sum_h([h.Clone(f"{h.GetName()}_c") for h in group_hists]) if group_hists else None
    if not h_mc_tot: raise SystemExit(f"[error] no MC hist for era {era}")
    stack = ROOT.THStack(f"stack_{era}", "")
    for h in group_hists: stack.Add(h)

    can = ROOT.TCanvas(f"can_{era}", f"Data/MC {era}", 900, 900)
    p_top = ROOT.TPad(f"p_top_{era}", "", 0, 0.28, 1, 1.0)
    p_bot = ROOT.TPad(f"p_bot_{era}", "", 0, 0.00, 1, 0.28)

    # margins – make Y label & legend visible
    for p in (p_top, p_bot):
        p.SetLeftMargin(0.14); p.SetRightMargin(0.04)
    p_top.SetBottomMargin(0.03)
    p_bot.SetTopMargin(0.05); p_bot.SetBottomMargin(0.38)
    if logy: p_top.SetLogy()
    p_top.Draw(); p_bot.Draw()

    # --- Top pad ---
    p_top.cd()
    stack.Draw("HIST")
    # Ensure axes titles exist on the THStack frame
    frame = stack.GetHistogram()
    frame.GetXaxis().SetTitle("BDT score")
    frame.GetYaxis().SetTitle("Events / bin")
    frame.GetXaxis().SetLabelSize(0)  # hide labels on top pad
    frame.GetYaxis().SetTitleSize(0.055); frame.GetYaxis().SetTitleOffset(1.0)

    if logy: stack.SetMinimum(0.01)
    for h in group_hists: h.Draw("HIST SAME")   # reinforce visibility
    h_data.Draw("E1 SAME")

    # Legend
    leg = ROOT.TLegend(0.66, 0.62, 0.90, 0.90)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.035)
    leg.AddEntry(h_data, "Data", "ep")
    for h in reversed(group_hists):
        leg.AddEntry(h, h.GetName().replace("MC_",""), "f")
    leg.Draw()

    # Title (ASCII only)
    tl = ROOT.TLatex(); tl.SetNDC(True); tl.SetTextSize(0.045)
    tl.DrawLatex(0.14, 0.93, title_top)
    ROOT.gPad.RedrawAxis()

    # --- Bottom pad: ratio ---
    p_bot.cd()
    ratio = h_data.Clone(f"ratio_{era}")
    ratio.Divide(h_mc_tot)
    ratio.SetMinimum(0.5); ratio.SetMaximum(1.5)
    ratio.GetYaxis().SetTitle("Data / MC")
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().SetTitleSize(0.11); ratio.GetYaxis().SetLabelSize(0.11); ratio.GetYaxis().SetTitleOffset(0.5)
    ratio.GetXaxis().SetTitle("BDT score")
    ratio.GetXaxis().SetTitleSize(0.13); ratio.GetXaxis().SetLabelSize(0.12)

    band = h_mc_tot.Clone(f"band_{era}")
    for b in range(1, band.GetNbinsX()+1):
        mc = h_mc_tot.GetBinContent(b); err = h_mc_tot.GetBinError(b)
        band.SetBinContent(b, 1.0 if mc>0 else 0.0)
        band.SetBinError(b, (err/mc) if mc>0 else 0.0)
    band.SetFillColorAlpha(ROOT.kGray+2, 0.35); band.SetMarkerSize(0)

    ratio.Draw("E1")
    band.Draw("E2 SAME")
    ratio.Draw("E1 SAME")

    # Diagnostic: MC-only snapshot
    can_mc = ROOT.TCanvas(f"can_mc_{era}", "MC only", 800, 600)
    stack.Draw("HIST")
    for h in group_hists: h.Draw("HIST SAME")
    stack.GetXaxis().SetTitle("BDT score"); stack.GetYaxis().SetTitle("Events / bin")
    out_mc = os.path.join(outdir, f"tmva_bdt_datamc_{era}_MCONLY.png")
    can_mc.SaveAs(out_mc)

    return can, h_mc_tot

def main():
    ap = argparse.ArgumentParser("BDT output: Data/MC per era with ratio")
    ap.add_argument("--data-era", action="append", default=[],
                    help="ERA:/path/data1.root[,/path/data2.root]")
    ap.add_argument("--mc-era", action="append", default=[],
                    help="ERA:PROC:/path/mc1.root[,/path/mc2.root]  (repeat)")
    ap.add_argument("--treename", default="mva")
    ap.add_argument("--bdt-branch", default="bdt_score")
    ap.add_argument("--mc-weight", default="weight")
    ap.add_argument("--mass-branch", default="m_mumu")
    ap.add_argument("--mass-min", type=float, default=110.0)
    ap.add_argument("--mass-max", type=float, default=150.0)
    ap.add_argument("--extra-cut", default="", help='extra RDataFrame filter, e.g. "njets>=1"')
    ap.add_argument("--nbins", type=int, default=40)
    ap.add_argument("--xmin", type=float, default=0.0)
    ap.add_argument("--xmax", type=float, default=1.0)
    ap.add_argument("--logy", action="store_true")
    ap.add_argument("--outdir", default="Outputs")
    ap.add_argument("--title", default="BDT output in signal fit mass range")
    args = ap.parse_args()

    print(f"[debug] Required branches: {{'{args.bdt_branch}', '{args.mc_weight}', '{args.mass_branch}'}}")
    os.makedirs(args.outdir, exist_ok=True)
    data_map_raw = parse_group_map(args.data_era, "--data-era")
    mc_map_raw   = parse_group_map(args.mc_era,   "--mc-era")
    eras = sorted(set(list(data_map_raw.keys()) + list(mc_map_raw.keys())))
    if not eras: raise SystemExit("[error] no eras provided")
    mass_sel = f"({args.mass_branch}>={args.mass_min} && {args.mass_branch}<={args.mass_max})"
    sel = mass_sel if not args.extra_cut else f"({mass_sel}) && ({args.extra_cut})"

    canvases = []
    for era in eras:
        data_files = []
        for v in data_map_raw.get(era, []):
            data_files.extend([p for p in v.split(",") if p])
        if not data_files: raise SystemExit(f"[error] no data files for era {era}")

        mc_entries = mc_map_raw.get(era, [])
        groups = []
        for entry in mc_entries:
            if ":" not in entry: raise SystemExit(f"[error] --mc-era needs ERA:PROC:files, got {entry}")
            proc, files = entry.split(":", 1)
            groups.append((proc, [p for p in files.split(",") if p]))

        h_data = hist_data(data_files, args.treename, args.bdt_branch, sel,
                           args.nbins, args.xmin, args.xmax)
        groups_sorted = sorted(groups, key=lambda x: proc_sort_key(x[0]))
        group_hists = []
        print(f"[groups] era {era}:")
        for proc, flist in groups_sorted:
            color = proc_color(proc)
            h = hist_mc_group(proc, flist, args.treename, args.bdt_branch,
                              args.mc_weight, sel, args.nbins, args.xmin, args.xmax,
                              color)
            if h:
                h.SetName(f"MC_{canonical_proc(proc)}")
                print(f"  - {proc}: files={len(flist)}  integral={h.Integral():.3f}  maxbin={h.GetMaximum():.3f}")
                group_hists.append(h)
            else:
                print(f"  - {proc}: (empty)")
        if not group_hists: raise SystemExit("[error] no MC histograms were built; check inputs.")

        h_mc_tot = group_hists[0].Clone("hMC_total"); h_mc_tot.Reset("ICES")
        for h in group_hists: h_mc_tot.Add(h)
        print(f"[yield] era {era}: MCtot = {h_mc_tot.Integral():.3f} ; DATA = {h_data.Integral():.3f}")

        title_top = f"{args.title}  --  {args.mass_min:.0f} < {args.mass_branch} < {args.mass_max:.0f} GeV"
        can, _ = make_panel(era, h_data, group_hists, args.nbins, args.xmin, args.xmax, title_top, args.logy, args.outdir)
        canvases.append((era, can))

    if len(canvases) == 1:
        era, c = canvases[0]
        out = os.path.join(args.outdir, f"tmva_bdt_datamc_{era}")
        c.SaveAs(out+".pdf"); c.SaveAs(out+".png")
        print(f"[done] {out}.(pdf|png)")
    else:
        grid = ROOT.TCanvas("grid","", 2400, 800)
        grid.Divide(len(canvases),1,0,0)
        for i,(_,c) in enumerate(canvases, start=1):
            grid.cd(i); c.DrawClonePad()
        out = os.path.join(args.outdir, "tmva_bdt_datamc_eras")
        grid.SaveAs(out+".pdf"); grid.SaveAs(out+".png")
        print(f"[done] {out}.(pdf|png)")

if __name__ == "__main__":
    main()
