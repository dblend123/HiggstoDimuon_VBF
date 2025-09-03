void Overlay_JetVars_WithRatio() {
    // ---------- Open inputs ----------
    TFile* f_data   = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/Data/DATA/Muon_FullDataRun22CD.root");
    TFile* f_DY     = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/DY/DY_FullMCRun22CD.root");
    TFile* f_TT     = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/TT/TT_FullMCRun22CD.root");
    TFile* f_WZ     = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/WZ/WZto2L2Q_FullMCRunCD.root");
    TFile* f_WZ3L   = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/WZto3LNu/WZto3LNu_FullMCRun22CD.root");
    TFile* f_Z2L2Q  = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/ZZto2L2Q/ZZto2L2Q_FullMCRun22CD.root");
    TFile* f_Z2L2N  = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/ZZto2L2Nu/ZZTo2L2Nu_FullMCRun22CD.root");
    TFile* f_signal = TFile::Open("/eos/user/d/dblend/Analyzed_SkimmedRun22CD/MC/VBF/VBF_MCRun22CD.root");

    if (!f_data || !f_DY || !f_TT || !f_WZ || !f_WZ3L || !f_Z2L2Q || !f_Z2L2N || !f_signal) {
        std::cerr << "[ERROR] One or more files couldn't be opened!" << std::endl;
        return;
    }

    // ---------- Helper to fetch and clone ----------
    auto get = [](TFile* f, const std::string& path) -> TH1* {
        auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
        return h ? (TH1*)h->Clone() : nullptr;
    };

    // ---------- Vars to overlay ----------
    std::vector<std::string> vars = {
        "DiJetPt", "Jet1Pt", "Jet2Pt",
        "DiJetdeta", "Jet1eta", "Jet2eta", "DiJetMass"
    };

    const std::string dir = "Jet_PreDiJet_VBF/";

    // ---------- Loop over variables ----------
    for (const auto& var : vars) {
        // Pull histograms
        TH1* h_data   = get(f_data,   dir + var);
        TH1* h_DY     = get(f_DY,     dir + var);
        TH1* h_TT     = get(f_TT,     dir + var);
        TH1* h_WZ     = get(f_WZ,     dir + var);
        TH1* h_WZ3L   = get(f_WZ3L,   dir + var);
        TH1* h_Z2L2Q  = get(f_Z2L2Q,  dir + var);
        TH1* h_Z2L2N  = get(f_Z2L2N,  dir + var);
        TH1* h_signal = get(f_signal, dir + var);

        if (!h_data || !h_DY || !h_TT || !h_WZ || !h_WZ3L || !h_Z2L2Q || !h_Z2L2N || !h_signal) {
            std::cerr << "[WARN] Missing hist(s) for " << var << " — skipping.\n";
            continue;
        }

        // Build VV = WZ + WZ3L + ZZ(2L2Q) + ZZ(2L2Nu)
        TH1* h_VV = (TH1*)h_WZ->Clone(("h_VV_"+var).c_str());
        h_VV->Add(h_WZ3L);
        h_VV->Add(h_Z2L2Q);
        h_VV->Add(h_Z2L2N);

        // ---------- Style (overlay lines + data points) ----------
        h_data->SetMarkerStyle(20); h_data->SetMarkerSize(1.2);
        h_data->SetMarkerColor(kBlack); h_data->SetLineColor(kBlack);

        h_DY->SetLineColor(kOrange+3);  h_DY->SetLineWidth(3);
        h_TT->SetLineColor(kRed+1);     h_TT->SetLineWidth(3);
        h_VV->SetLineColor(kAzure+2);   h_VV->SetLineWidth(3);
        h_signal->SetLineColor(kMagenta+1); h_signal->SetLineWidth(3);

        // ---------- Canvas with ratio pad ----------
        TCanvas* c = new TCanvas(("c_"+var).c_str(), var.c_str(), 800, 800);
        TPad* pad1 = new TPad(("pad1_"+var).c_str(), "Top pad",    0, 0.32, 1, 1);
        TPad* pad2 = new TPad(("pad2_"+var).c_str(), "Ratio pad",  0, 0.02, 1, 0.32);
        pad1->SetBottomMargin(0.03);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.30);
        pad1->Draw(); pad2->Draw();

        // ---------- Top pad: overlays ----------
        pad1->cd();
        // Optional: comment out next line if you prefer linear scale
         pad1->SetLogy();

        // Pick a sensible y-max to fit data and MC comfortably
        double ymax = std::max({
            h_data->GetMaximum(),
            h_DY->GetMaximum(),
            h_TT->GetMaximum(),
            h_VV->GetMaximum(),
            h_signal->GetMaximum()
        });
        ymax *= 1.6;

        // Use data as the "frame" to inherit axis titles/bins; then draw others
        h_data->SetTitle("");
        h_data->GetYaxis()->SetTitle("Events");
        h_data->GetYaxis()->SetTitleSize(0.05);
        h_data->GetYaxis()->SetLabelSize(0.04);
        h_data->GetYaxis()->SetTitleOffset(1.2);
        h_data->SetMaximum(ymax);
        h_data->SetMinimum(0.1);
        h_data->Draw("E1");              // data points
        h_DY->Draw("HIST SAME");
        h_TT->Draw("HIST SAME");
        h_VV->Draw("HIST SAME");
        h_signal->Draw("HIST SAME");

        // Legend
        TLegend* leg = new TLegend(0.60, 0.60, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_signal, "VBF H#rightarrow#mu#mu", "l");
        leg->AddEntry(h_TT,     "t#bar{t} + Jets",        "l");
        leg->AddEntry(h_VV,     "VV + Jets",              "l");
        leg->AddEntry(h_DY,     "DY + Jets",              "l");
        leg->AddEntry(h_data,   "Data",                   "lep");
        leg->Draw();

        // CMS prelim + lumi
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.13, 0.93, "CMS #scale[0.8]{#it{Preliminary}}");
        latex.DrawLatex(0.55, 0.93, "#sqrt{s} = 13 TeV, 7.98 fb^{-1}");

        // ---------- Bottom pad: Data/MC ratio ----------
        pad2->cd();

        // MC total for ratio = DY + TT + VV (exclude signal)
        TH1* h_mctot = (TH1*)h_DY->Clone(("h_mctot_"+var).c_str());
        h_mctot->Add(h_TT);
        h_mctot->Add(h_VV);

        TH1* h_ratio = (TH1*)h_data->Clone(("h_ratio_"+var).c_str());
        // Avoid divide-by-zero bins
        for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
            double mc = h_mctot->GetBinContent(i);
            if (mc <= 0) {
                h_ratio->SetBinContent(i, 0);
                h_ratio->SetBinError(i, 0);
            }
        }
        h_ratio->Divide(h_mctot);

        // Ratio style
        h_ratio->GetYaxis()->SetTitle("Data/MC");
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetYaxis()->SetNdivisions(505);

        h_ratio->GetXaxis()->SetLabelSize(0.12);
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        // Keep the underlying x-title from the input histogram
        h_ratio->SetMinimum(0.0);
        h_ratio->SetMaximum(2.0);
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetLineColor(kBlack);

        h_ratio->Draw("E1");

        // Optional: draw a line at 1
        TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0,
                                h_ratio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineStyle(2);
        line->Draw("SAME");

        // ---------- Save ----------
        TString outdir = "/eos/user/d/dblend/StackedPlots/Jets/";
        TString base = "OverlayRatio_" + TString(var.c_str());
        c->SaveAs(outdir + base + ".pdf");
        c->SaveAs(outdir + base + ".png");
    }

    // ---------- Cleanup ----------
    f_data->Close(); f_DY->Close(); f_TT->Close(); f_signal->Close();
    f_WZ->Close(); f_WZ3L->Close(); f_Z2L2Q->Close(); f_Z2L2N->Close();
}

