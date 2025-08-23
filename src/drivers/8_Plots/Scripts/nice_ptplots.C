#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>
#include <vector>

void DrawCMSLabel() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.25, 0.82, "#font[62]{CMS}");
    latex.DrawLatex(0.25, 0.78, "#it{Simulation}");
    latex.DrawLatex(0.25, 0.74, "#it{Preliminary}");
    latex.DrawLatex(0.22, 0.65, "#font[62]{Signal MC 2022}");
}

void nice_ptplots() {
    // Disable the stat box
    gStyle->SetOptStat(0);

    TFile *file = new TFile("output_Muon_FullSignalRun22.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open file output_Muon_FullSignalRun22.root" << std::endl;
        return;
    }

    TDirectory *DMdir = (TDirectory*)file->Get("DM_VBF");
    if (!DMdir) {
        std::cerr << "Error: Directory 'DM_VBF' not found in file." << std::endl;
        return;
    }

    // Bin edges for mass histograms
    std::vector<int> ptBins = {20, 25, 30, 35, 40, 50, 60, 80, 120}; // Last bin is upper edge

    for (size_t i = 0; i < ptBins.size() - 1; ++i) {
        int ptLow = ptBins[i];
        int ptHigh = ptBins[i + 1];
        std::string histName = "DiMuonMass_pt" + std::to_string(ptLow);
        TH1F *hist = (TH1F*)DMdir->Get(histName.c_str());

        if (!hist) {
            std::cerr << "Error: Histogram '" << histName << "' not found in directory 'DM_VBF'." << std::endl;
            continue;
        }

        hist->GetXaxis()->SetTitle("Invariant Mass [GeV]");
        hist->GetYaxis()->SetTitle("Events");
        hist->SetTitle(Form("%d < m_{#mu#mu} < %d GeV", ptLow, ptHigh));
        hist->SetFillColor(kGreen - 9);

        TCanvas *canvas = new TCanvas(("canvas_" + histName).c_str(), (histName + " Canvas").c_str(), 800, 600);
        canvas->SetTickx();
        canvas->SetTicky();
        canvas->SetGridx();
        canvas->SetGridy();
        gPad->SetLeftMargin(0.15);

        Double_t totalEvents = 0.0;
        Int_t numBins = hist->GetNbinsX();
        for (Int_t j = 1; j <= numBins; ++j) {
            totalEvents += hist->GetBinContent(j);
        }
        hist->SetEntries(totalEvents);
        hist->Draw("HIST");
        DrawCMSLabel();

        canvas->SaveAs((histName + ".png").c_str());
    }

    // DiMuonPtVsEta 2D histogram
    TH2F *hist2D = (TH2F*)DMdir->Get("DiMuonPtVsEta");
    if (!hist2D) {
        std::cerr << "Error: Histogram 'DiMuonPtVsEta' not found in directory 'DM_VBF'." << std::endl;
        return;
    }

    hist2D->GetXaxis()->SetTitle("Dimuon #eta");
    hist2D->GetYaxis()->SetTitle("Dimuon p_{T} [GeV]");
    hist2D->SetTitle("Dimuon p_{T} vs #eta");

    TCanvas *canvas2D = new TCanvas("canvas_DiMuonPtVsEta", "DiMuonPtVsEta Canvas", 800, 600);
    canvas2D->SetTickx();
    canvas2D->SetTicky();
    canvas2D->SetGridx();
    canvas2D->SetGridy();
    gPad->SetLeftMargin(0.15);

    hist2D->Draw("COLZ");
    DrawCMSLabel();

    canvas2D->SaveAs("DiMuonPtVsEta.png");
}

