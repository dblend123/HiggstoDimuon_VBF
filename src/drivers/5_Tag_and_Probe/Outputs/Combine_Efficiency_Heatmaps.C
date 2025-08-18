void Combine_Efficiency_Heatmaps(const char* inputFile = "TagAndProbe_Heatmaps_Muon_FullDataRun22.root", const char* outputFile = "output_TagAndProbe_Heatmaps_Muon_FullDataRun22.root") {
    // Load input file
    TFile* infile = TFile::Open(inputFile, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: Cannot open input file.\n";
        return;
    }

    // Get histograms
    auto h2_pt_eta_pass = (TH2D*)infile->Get("h2_pt_eta_pass");
    auto h2_pt_eta_fail = (TH2D*)infile->Get("h2_pt_eta_fail");
    if (!h2_pt_eta_pass || !h2_pt_eta_fail) {
        std::cerr << "Error: Missing histograms in input file.\n";
        return;
    }

    // Clone total and compute efficiency
    auto h2_pt_eta_total = (TH2D*)h2_pt_eta_pass->Clone("h2_pt_eta_total");
    h2_pt_eta_total->Add(h2_pt_eta_fail);

    auto h2_pt_eta_eff = (TH2D*)h2_pt_eta_pass->Clone("h2_pt_eta_eff");
    h2_pt_eta_eff->Divide(h2_pt_eta_total);

    //h2_pt_eta_eff->SetTitle("DY-Run22 p_{T} vs. #eta Efficiency Map;#eta;p_{T} [GeV]");
    h2_pt_eta_eff->SetTitle("Data-Run22 p_{T} vs. #eta Efficiency Map;#eta;p_{T} [GeV]");
    h2_pt_eta_eff->SetMinimum(0.5);
    h2_pt_eta_eff->SetMaximum(1.0);

    gStyle->SetOptStat(0);

    // Canvas: Passing Probes
    TCanvas* c1 = new TCanvas("c1", "Passing Probes", 800, 600);
    h2_pt_eta_pass->SetTitle("Passing Probes;#eta;p_{T} [GeV]");
    h2_pt_eta_pass->Draw("COLZ");

    // Canvas: Failing Probes
    TCanvas* c2 = new TCanvas("c2", "Failing Probes", 800, 600);
    h2_pt_eta_fail->SetTitle("Failing Probes;#eta;p_{T} [GeV]");
    h2_pt_eta_fail->Draw("COLZ");

    // Canvas: Efficiency
    TCanvas* c3 = new TCanvas("c3", "Efficiency", 800, 600);
    h2_pt_eta_eff->Draw("COLZ");

    // Draw bin values + binomial errors
    TLatex latex;
    latex.SetTextSize(0.03);
    latex.SetTextAlign(22); // Centered
    
    for (int ix = 1; ix <= h2_pt_eta_eff->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h2_pt_eta_eff->GetNbinsY(); ++iy) {
            double pass  = h2_pt_eta_pass->GetBinContent(ix, iy);
            double total = h2_pt_eta_total->GetBinContent(ix, iy);
            if (total < 1e-6) continue; // avoid division by zero

            double eff = pass / total;
            double err = std::sqrt(eff * (1.0 - eff) / total);

            double x = h2_pt_eta_eff->GetXaxis()->GetBinCenter(ix);
            double y = h2_pt_eta_eff->GetYaxis()->GetBinCenter(iy);
            latex.DrawLatex(x, y, Form("%.2f#pm%.2f", eff, err));
        }
    }
    
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.13, 0.91, "CMS Preliminary");
    latex.DrawLatex(0.62, 0.91, "Data (Tight ID)");
    //latex.DrawLatex(0.62, 0.91, "Simulation (Tight ID)");

    // Save to output
    TFile* outfile = new TFile(outputFile, "RECREATE");
    h2_pt_eta_pass->Write();
    h2_pt_eta_fail->Write();
    h2_pt_eta_total->Write();
    h2_pt_eta_eff->Write();
    outfile->Close();

    std::string baseName = std::string(outputFile);
    if (baseName.find(".root") != std::string::npos)
        baseName = baseName.substr(0, baseName.find(".root"));

    //c3->SaveAs((baseName + "_eff.pdf").c_str());
    c3->SaveAs((baseName + "_eff.png").c_str());

    std::cout << "Efficiency heatmaps written to " << outputFile << std::endl;
}

