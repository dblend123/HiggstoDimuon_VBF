void StackedPlot_AllMC_WithBlinding() {
    // Input files
    TFile* f_data = TFile::Open("Muon_DataRun22CD.root");
    TFile* f_DY   = TFile::Open("DY_MCRun22CD.root");
    TFile* f_TT   = TFile::Open("TT_MCRun22CD.root");
    TFile* f_WZ   = TFile::Open("WZ_MCRun22CD.root");
    TFile* f_WZ3L = TFile::Open("WZto3LNu_MCRun22CD.root");
    TFile* f_Z2L2Q= TFile::Open("ZiZTo2L2Q_MCRun22CD.root");
    TFile* f_Z2L2N= TFile::Open("ZZTo2L2Nu_MCRun22CD.root");
    TFile* f_signal = TFile::Open("VBF_Run22BCD.root");

    if (!f_data || !f_DY || !f_TT || !f_WZ || !f_WZ3L || !f_Z2L2Q || !f_Z2L2N || !f_signal) {
        std::cerr << "One or more files couldn't be opened!" << std::endl;
        return;
    }

    // Get histogram helper
    auto get = [](TFile* f, const std::string& path) -> TH1* {
        auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
        return h ? (TH1*)h->Clone() : nullptr;
    };

    // Get histograms
    TH1* h_data  = get(f_data,  "DM_Pre-Selection/DiMuonMass");
    TH1* h_DY    = get(f_DY,    "corr_DM_VBF/DiMuonMass");
    TH1* h_TT    = get(f_TT,    "corr_DM_VBF/DiMuonMass");
    TH1* h_WZ    = get(f_WZ,    "corr_DM_VBF/DiMuonMass");
    TH1* h_WZ3L  = get(f_WZ3L,  "corr_DM_VBF/DiMuonMass");
    TH1* h_Z2L2Q = get(f_Z2L2Q, "corr_DM_VBF/DiMuonMass");
    TH1* h_Z2L2N = get(f_Z2L2N, "corr_DM_VBF/DiMuonMass");
    TH1* h_signal = get(f_signal, "corr_DM_VBF/DiMuonMass");

    if (!h_data || !h_DY || !h_TT || !h_WZ || !h_WZ3L || !h_Z2L2Q || !h_Z2L2N || !h_signal) {
        std::cerr << "One or more histograms are missing!" << std::endl;
        return;
    }

    // Combine VV
    TH1* h_VV = (TH1*)h_WZ->Clone("h_VV");
    h_VV->Add(h_WZ3L);
    h_VV->Add(h_WZ);
    h_VV->Add(h_Z2L2Q);
    h_VV->Add(h_Z2L2N);

    // Style
    h_data->SetMarkerStyle(20); h_data->SetMarkerColor(kBlack); h_data->SetLineColor(kBlack); h_data->SetMarkerSize(1.2);
    h_DY->SetFillColor(kOrange-2); h_DY->SetLineColor(kOrange+3);
    h_TT->SetFillColor(kRed-4); h_TT->SetLineColor(kRed+1);
    h_VV->SetFillColor(kAzure-9); h_VV->SetLineColor(kAzure+2);
    h_signal->SetFillColor(kMagenta+1); h_signal->SetLineWidth(3); 

    // Stack
    THStack* hs = new THStack("hs", ";;Events / 1 [GeV]");
    hs->Add(h_signal, "HIST");
    hs->Add(h_TT, "HIST");
    hs->Add(h_DY, "HIST");
    hs->Add(h_VV, "HIST");

    // Canvas and pads
    TCanvas* c = new TCanvas("c", "CMS Data vs MC", 800, 800);
    TPad* pad1 = new TPad("pad1", "Top pad", 0, 0.32, 1, 1);
    TPad* pad2 = new TPad("pad2", "Bottom pad", 0, 0.02, 1, 0.32);
    pad1->SetBottomMargin(0.03);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.3);
    pad1->Draw(); pad2->Draw();

    // --- Top Pad ---
    pad1->cd();
    pad1->SetLogy();
    hs->SetMinimum(0.1);
    hs->Draw("HIST");

    double max_data = h_data->GetMaximum();
    double max_stack = hs->GetMaximum();
    hs->SetMaximum(std::max(max_data, max_stack) * 2.0);
    hs->GetYaxis()->SetTitle("Events / 1 [GeV]");
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetLabelSize(0.04);
    hs->GetYaxis()->SetTitleOffset(1.2);

    // Blinded data
    TH1* h_data_blinded = (TH1*)h_data->Clone("h_data_blinded");
    int bin_min_data = h_data->FindBin(110.0);
    int bin_max_data = h_data->FindBin(150.0);
    for (int i = bin_min_data; i <= bin_max_data; ++i) {
        h_data_blinded->SetBinContent(i, 0);
        h_data_blinded->SetBinError(i, 0);
    }
    h_data_blinded->Draw("E1 SAME");

    // Vertical lines at blinded region
    TLine* l1 = new TLine(110, 0.1, 110, hs->GetMaximum());
    TLine* l2 = new TLine(150, 0.1, 150, hs->GetMaximum());
    l1->SetLineStyle(2); l1->SetLineWidth(2);
    l2->SetLineStyle(2); l2->SetLineWidth(2);
    l1->Draw(); l2->Draw();

    // Legend
    TLegend* leg = new TLegend(0.60, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_signal, "VBF H#rightarrow#mu#mu", "f");
    leg->AddEntry(h_TT,   "t#bar{t} + Jets", "f");
    leg->AddEntry(h_VV,   "VV + Jets", "f");
    leg->AddEntry(h_DY,   "DY + Jets", "f");
    leg->AddEntry(h_data, "Data", "lep");
    leg->Draw();

    // CMS text
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.13, 0.93, "CMS #scale[0.8]{#it{Preliminary}}");
    latex.DrawLatex(0.55, 0.93, "#sqrt{s} = 13 TeV, 7.98 fb^{-1}");

    // --- Ratio Pad ---
    pad2->cd();

    // Clone data and MC histograms
    TH1* h_ratio = (TH1*)h_data->Clone("h_ratio");
    TH1* h_mc_total = (TH1*)h_VV->Clone("h_mc_total");
    h_mc_total->Add(h_TT);
    h_mc_total->Add(h_DY);

    // Calculate Data/MC ratio
    h_ratio->Divide(h_mc_total);

    // --- Blind ratio in 110–150 GeV ---
    int bin_min = h_ratio->FindBin(110.0);
    int bin_max = h_ratio->FindBin(150.0);
    for (int i = bin_min; i <= bin_max; ++i) {
        h_ratio->SetBinContent(i, 0);
        h_ratio->SetBinError(i, 0);
    }

    // Style
    h_ratio->GetYaxis()->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetTitleSize(0.12);
    h_ratio->GetYaxis()->SetTitleOffset(0.4);
    h_ratio->GetYaxis()->SetLabelSize(0.1);
    h_ratio->GetYaxis()->SetNdivisions(505);

    h_ratio->GetXaxis()->SetLabelSize(0.12);
    h_ratio->GetXaxis()->SetTitleSize(0.12);
    h_ratio->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}}^{H} [GeV]");

    h_ratio->SetMinimum(0.0);
    h_ratio->SetMaximum(2.0);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerStyle(20);

    // Draw
    h_ratio->Draw("E1");

    // Save
    c->SaveAs("Stacked_DiMuonMass_blinded_withRatio.pdf");
    c->SaveAs("Stacked_DiMuonMass_blinded_withRatio.png");

    // Cleanup
    f_data->Close(); f_DY->Close(); f_TT->Close(); f_signal->Close();
    f_WZ->Close(); f_WZ3L->Close(); f_Z2L2Q->Close(); f_Z2L2N->Close();
}

