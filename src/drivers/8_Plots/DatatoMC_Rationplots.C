void DatatoMC_Rationplots(const char* filename_mc = "output_Muon_FullDYRun22_v5-v2.root", const char* filename_data = "output_Muon_FullDataRun22CD.root") {
    // Open MC and Data files
    TFile* file_mc   = TFile::Open(filename_mc);
    TFile* file_data = TFile::Open(filename_data);

    if (!file_mc || !file_data || file_mc->IsZombie() || file_data->IsZombie()) {
        std::cerr << "Could not open one of the input files." << std::endl;
        return;
    }

/*    // Get histograms
    auto* h_uncorr_Z     = dynamic_cast<TH1*>(file->Get("DM_VBF/DiMuonMass_Z"));
    auto* h_corr_Z       = dynamic_cast<TH1*>(file->Get("corr_DM_VBF/DiMuonMass_Z"));
    auto* h_uncorr_Higgs = dynamic_cast<TH1*>(file->Get("DM_VBF/DiMuonMass_SB_SR1"));
    auto* h_corr_Higgs   = dynamic_cast<TH1*>(file->Get("corr_DM_VBF/DiMuonMass_SB_SR1"));

    if (!h_uncorr_Z || !h_corr_Z || !h_uncorr_Higgs || !h_corr_Higgs) {
        std::cerr << "Missing histograms in file." << std::endl;
        return;
    }

    // Clone to detach from file
    h_uncorr_Z     = (TH1*)h_uncorr_Z->Clone("h_uncorr_Z");
    h_corr_Z       = (TH1*)h_corr_Z->Clone("h_corr_Z");
    h_uncorr_Higgs = (TH1*)h_uncorr_Higgs->Clone("h_uncorr_Higgs");
    h_corr_Higgs   = (TH1*)h_corr_Higgs->Clone("h_corr_Higgs");

    // Normalize
    h_uncorr_Z->Scale(1.0 / h_uncorr_Z->Integral());
    h_corr_Z->Scale(1.0 / h_corr_Z->Integral());
    h_uncorr_Higgs->Scale(1.0 / h_uncorr_Higgs->Integral());
    h_corr_Higgs->Scale(1.0 / h_corr_Higgs->Integral());

    // Style
    h_uncorr_Z->SetLineColor(kBlue);
    h_corr_Z->SetLineColor(kRed);
    h_uncorr_Z->SetLineWidth(2);
    h_corr_Z->SetLineWidth(2);
    h_uncorr_Z->SetTitle("DiMuon Mass (Z region);M(#mu#mu) [GeV];A.U.");

    h_uncorr_Higgs->SetLineColor(kBlue);
    h_corr_Higgs->SetLineColor(kRed);
    h_uncorr_Higgs->SetLineWidth(2);
    h_corr_Higgs->SetLineWidth(2);
    h_uncorr_Higgs->SetTitle("DiMuon Mass (Higgs region);M(#mu#mu) [GeV];A.U.");

    // ===== Z REGION PLOT =====
    TCanvas* c = new TCanvas("c", "Z Region", 800, 800);

    TPad* pad1 = new TPad("pad1", "Main Z", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    pad1->SetGrid();
    h_uncorr_Z->Draw("HIST");
    h_corr_Z->Draw("HIST SAME");

    TLegend* leg_Z = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg_Z->AddEntry(h_uncorr_Z, "Uncorrected", "l");
    leg_Z->AddEntry(h_corr_Z, "Rochester Corrected", "l");
    leg_Z->Draw();
    c->cd();

    TPad* pad2 = new TPad("pad2", "Ratio Z", 0, 0.0, 1, 0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.3);
    pad2->Draw();
    pad2->cd();
    pad2->SetGridy();

    auto* ratio = (TH1*)h_corr_Z->Clone("ratio_Z");
    ratio->Divide(h_uncorr_Z);
    ratio->SetTitle(";M(#mu#mu) [GeV];Corr / Uncorr");
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetYaxis()->SetLabelSize(0.08);
    ratio->GetXaxis()->SetLabelSize(0.08);
    ratio->GetXaxis()->SetTitleSize(0.1);
    ratio->GetXaxis()->SetTitleOffset(1.0);
    ratio->SetLineColor(kBlack);
    ratio->SetLineWidth(2);
    ratio->SetMinimum(0.8);
    ratio->SetMaximum(1.2);
    ratio->Draw("HIST");

    // Save Z canvas
    //c->SaveAs("Signal_Z_dimuon_mass_comparison.pdf");
    //c->SaveAs("Signal_Z_dimuon_mass_comparison.png");
    c->SaveAs("DY_Z_dimuon_mass_comparison.pdf");
    c->SaveAs("DY_Z_dimuon_mass_comparison.png");

    // ===== HIGGS REGION PLOT =====
    TCanvas* c_Higgs = new TCanvas("c_Higgs", "Higgs Region", 800, 800);

    TPad* padH1 = new TPad("padH1", "Main Higgs", 0, 0.3, 1, 1.0);
    padH1->SetBottomMargin(0.02);
    padH1->Draw();
    padH1->cd();
    padH1->SetGrid();

    // Optional: zoom to Higgs mass window
    h_uncorr_Higgs->GetXaxis()->SetRangeUser(110, 150);
    h_corr_Higgs->GetXaxis()->SetRangeUser(110, 150);

    h_uncorr_Higgs->Draw("HIST");
    h_corr_Higgs->Draw("HIST SAME");

    TLegend* leg_Higgs = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg_Higgs->AddEntry(h_uncorr_Higgs, "Uncorrected Higgs", "l");
    leg_Higgs->AddEntry(h_corr_Higgs, "Rochester Corrected Higgs", "l");
    leg_Higgs->Draw();
    c_Higgs->cd();

    TPad* padH2 = new TPad("padH2", "Ratio Higgs", 0, 0.0, 1, 0.3);
    padH2->SetTopMargin(0.05);
    padH2->SetBottomMargin(0.3);
    padH2->Draw();
    padH2->cd();
    padH2->SetGridy();

    auto* ratio_Higgs = (TH1*)h_corr_Higgs->Clone("ratio_Higgs");
    ratio_Higgs->Divide(h_uncorr_Higgs);
    ratio_Higgs->SetTitle(";M(#mu#mu) [GeV];Corr / Uncorr");
    ratio_Higgs->GetYaxis()->SetTitleSize(0.1);
    ratio_Higgs->GetYaxis()->SetTitleOffset(0.4);
    ratio_Higgs->GetYaxis()->SetLabelSize(0.08);
    ratio_Higgs->GetXaxis()->SetLabelSize(0.08);
    ratio_Higgs->GetXaxis()->SetTitleSize(0.1);
    ratio_Higgs->GetXaxis()->SetTitleOffset(1.0);
    ratio_Higgs->SetLineColor(kBlack);
    ratio_Higgs->SetLineWidth(2);
    ratio_Higgs->SetMinimum(0.8);
    ratio_Higgs->SetMaximum(1.2);
    ratio_Higgs->Draw("HIST");

    // Save Higgs canvas
    //c_Higgs->SaveAs("Signal_Higgs_dimuon_mass_comparison.pdf");
    //c_Higgs->SaveAs("Signal_Higgs_dimuon_mass_comparison.png");
    c_Higgs->SaveAs("DY_Higgs_dimuon_mass_comparison.pdf");
    c_Higgs->SaveAs("DY_Higgs_dimuon_mass_comparison.png");

    file->Close();*/
    
    // ===== DATA / MC COMPARISON IN Z REGION =====
    //TFile* file_data = TFile::Open(filename_data);
    //TFile* file_mc   = TFile::Open(filename_mc);

    if (!file_data || !file_mc) {
        std::cerr << "Could not open one of the files for Data/MC comparison." << std::endl;
        return;
    }

    //auto* h_data_Z = dynamic_cast<TH1*>(file_data->Get("corr_DM_VBF/DiMuonMass_Z"));
    //auto* h_mc_Z   = dynamic_cast<TH1*>(file_mc->Get("corr_DM_VBF/DiMuonMass_Z"));
    auto* h_data_Z = dynamic_cast<TH1*>(file_data->Get("DM_VBF/DiMuonMass_Z"));
    auto* h_mc_Z   = dynamic_cast<TH1*>(file_mc->Get("DM_VBF/DiMuonMass_Z"));

    if (!h_data_Z || !h_mc_Z) {
        std::cerr << "Missing corr_DM_VBF/DiMuonMass_Z in data or MC." << std::endl;
        return;
    }

    // Clone and normalize (optional; or comment these lines if using absolute scale)
    h_data_Z = (TH1*)h_data_Z->Clone("h_data_Z");
    h_mc_Z   = (TH1*)h_mc_Z->Clone("h_mc_Z");
    //h_data_Z->Scale(1.0 / h_data_Z->Integral());
    //h_mc_Z->Scale(1.0 / h_mc_Z->Integral());

    // Styling
    h_data_Z->SetLineColor(kBlack);
    h_data_Z->SetMarkerStyle(20);
    h_data_Z->SetMarkerColor(kBlack);
    h_data_Z->SetTitle("Data vs MC: DiMuon Mass (Z);M(#mu#mu) [GeV];A.U.");
    h_data_Z->SetMaximum( std::max(h_data_Z->GetMaximum(), h_mc_Z->GetMaximum()) * 1.2 );

    h_mc_Z->SetLineColor(kGreen+2);
    h_mc_Z->SetLineWidth(2);

    // Create canvas
    TCanvas* c_ratio = new TCanvas("c_ratio", "Data / MC Z", 800, 800);
    TPad* pad_top = new TPad("pad_top", "Top", 0, 0.3, 1, 1.0);
    pad_top->SetBottomMargin(0.02);
    pad_top->Draw();
    pad_top->cd();
    pad_top->SetGrid();

    h_data_Z->Draw("E1");
    h_mc_Z->Draw("HIST SAME");

    TLegend* leg_dm = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg_dm->AddEntry(h_data_Z, "Data", "lep");
    leg_dm->AddEntry(h_mc_Z, "MC (DY)", "l");
    leg_dm->Draw();

    c_ratio->cd();

    // Ratio pad
    TPad* pad_bot = new TPad("pad_bot", "Bottom", 0, 0.0, 1, 0.3);
    pad_bot->SetTopMargin(0.05);
    pad_bot->SetBottomMargin(0.3);
    pad_bot->Draw();
    pad_bot->cd();
    pad_bot->SetGridy();

    auto* ratio_dm = (TH1*)h_data_Z->Clone("ratio_dm");
    ratio_dm->Divide(h_mc_Z);
    ratio_dm->SetTitle(";M(#mu#mu) [GeV];Data / MC");
    ratio_dm->GetYaxis()->SetTitleSize(0.1);
    ratio_dm->GetYaxis()->SetTitleOffset(0.4);
    ratio_dm->GetYaxis()->SetLabelSize(0.08);
    ratio_dm->GetXaxis()->SetLabelSize(0.08);
    ratio_dm->GetXaxis()->SetTitleSize(0.1);
    ratio_dm->GetXaxis()->SetTitleOffset(1.0);
    ratio_dm->SetLineColor(kBlack);
    ratio_dm->SetLineWidth(2);
    ratio_dm->SetMinimum(0.5); // adjust as needed
    ratio_dm->SetMaximum(1.5); // adjust as needed
    ratio_dm->Draw("E1");

    // Save plot
    c_ratio->SaveAs("DataMC_Z_dimuon_mass_comparison.pdf");
    c_ratio->SaveAs("DataMC_Z_dimuon_mass_comparison.png");
    //c_ratio->SaveAs("NoScale_DataMC_Z_dimuon_mass_comparison.pdf");
    //c_ratio->SaveAs("NoScale_DataMC_Z_dimuon_mass_comparison.png");

    file_data->Close();
    file_mc->Close();
}

