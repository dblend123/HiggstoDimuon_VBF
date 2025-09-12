#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include <iostream>
#include <algorithm>
#include <string>

void StackedPlot_AllMC_WithBlinding(TString era = "12022") {
    // Build base paths with era inserted right after Analyzed_SkimmedFiles/
    const TString base = "/eos/user/d/dblend/Analyzed_SkimmedFiles/" + era + "/";

    // -------- Open inputs --------
    TFile* f_data       = TFile::Open(base + "MC/DATA/DATA_all.root");
    TFile* f_DY         = TFile::Open(base + "MC/DY/DY_all.root");
    TFile* f_TT         = TFile::Open(base + "MC/TT/TT_all.root");

    // Diboson (VV)
    TFile* f_WZto2L2Q   = TFile::Open(base + "MC/WZto2L2Q/WZto2L2Q_all.root");
    TFile* f_WZto3LNu   = TFile::Open(base + "MC/WZto3LNu/WZto3LNu_all.root");
    TFile* f_ZZto2L2Nu  = TFile::Open(base + "MC/ZZto2L2Nu/ZZto2L2Nu_all.root");
    TFile* f_ZZto2L2Q   = TFile::Open(base + "MC/ZZto2L2Q/ZZto2L2Q_all.root");
    TFile* f_ZZto4L     = TFile::Open(base + "MC/ZZto4L/ZZto4L_all.root");

    // Triboson (VVV)
    TFile* f_WWW        = TFile::Open(base + "MC/WWW/WWW_all.root");
    TFile* f_WWZ        = TFile::Open(base + "MC/WWZ/WWZ_all.root");
    TFile* f_WZZ        = TFile::Open(base + "MC/WZZ/WZZ_all.root");
    TFile* f_ZZZ        = TFile::Open(base + "MC/ZZZ/ZZZ_all.root");

    // Signal
    TFile* f_signal     = TFile::Open(base + "MC/VBF/VBF_all.root");

    if (!f_data || !f_DY || !f_TT ||
        !f_WZto2L2Q || !f_WZto3LNu || !f_ZZto2L2Nu || !f_ZZto2L2Q || !f_ZZto4L ||
        !f_WWW || !f_WWZ || !f_WZZ || !f_ZZZ ||
        !f_signal) {
        std::cerr << "[ERROR] One or more files couldn't be opened!" << std::endl;
        return;
    }

    auto get = [](TFile* f, const std::string& path) -> TH1* {
        auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
        return h ? (TH1*)h->Clone() : nullptr;
    };

    const std::string hpath = "corr_DM_VBF/DiMuonMass";

    // -------- Get histograms --------
    TH1* h_data     = get(f_data,      hpath);
    TH1* h_DY       = get(f_DY,        hpath);
    TH1* h_TT       = get(f_TT,        hpath);

    // Diboson
    TH1* h_WZ2L2Q   = get(f_WZto2L2Q,  hpath);
    TH1* h_WZ3LNu   = get(f_WZto3LNu,  hpath);
    TH1* h_ZZ2L2Nu  = get(f_ZZto2L2Nu, hpath);
    TH1* h_ZZ2L2Q   = get(f_ZZto2L2Q,  hpath);
    TH1* h_ZZ4L     = get(f_ZZto4L,    hpath);

    // Triboson
    TH1* h_WWW      = get(f_WWW,       hpath);
    TH1* h_WWZ      = get(f_WWZ,       hpath);
    TH1* h_WZZ      = get(f_WZZ,       hpath);
    TH1* h_ZZZ      = get(f_ZZZ,       hpath);

    // Signal
    TH1* h_signal   = get(f_signal,    hpath);

    if (!h_data || !h_DY || !h_TT ||
        !h_WZ2L2Q || !h_WZ3LNu || !h_ZZ2L2Nu || !h_ZZ2L2Q || !h_ZZ4L ||
        !h_WWW || !h_WWZ || !h_WZZ || !h_ZZZ ||
        !h_signal) {
        std::cerr << "[ERROR] One or more histograms are missing!" << std::endl;
        return;
    }

    // -------- Build groups --------
    // VV (diboson)
    TH1* h_VV = (TH1*)h_WZ2L2Q->Clone("h_VV");
    h_VV->Add(h_WZ3LNu);
    h_VV->Add(h_ZZ2L2Q);
    h_VV->Add(h_ZZ2L2Nu);
    h_VV->Add(h_ZZ4L);

    // VVV (triboson)
    TH1* h_VVV = (TH1*)h_WWW->Clone("h_VVV");
    h_VVV->Add(h_WWZ);
    h_VVV->Add(h_WZZ);
    h_VVV->Add(h_ZZZ);

    // -------- Style --------
    h_data->SetMarkerStyle(20); h_data->SetMarkerColor(kBlack); h_data->SetLineColor(kBlack); h_data->SetMarkerSize(1.2);
    h_DY->SetFillColor(kOrange-2); h_DY->SetLineColor(kOrange+3);
    h_TT->SetFillColor(kRed-4);    h_TT->SetLineColor(kRed+1);
    h_VV->SetFillColor(kAzure-9);  h_VV->SetLineColor(kAzure+2);
    h_VVV->SetFillColor(kGreen-7); h_VVV->SetLineColor(kGreen+3);
    h_signal->SetFillColor(kMagenta+1); h_signal->SetLineWidth(3);

    // -------- Stack (signal on top for visibility) --------
    THStack* hs = new THStack("hs", ";;Events / 1 [GeV]");
    hs->Add(h_DY, "HIST");
    hs->Add(h_TT, "HIST");
    hs->Add(h_VV, "HIST");
    hs->Add(h_VVV, "HIST");
    hs->Add(h_signal, "HIST");

    // -------- Canvas & pads --------
    TCanvas* c = new TCanvas("c", "CMS Data vs MC", 800, 800);
    TPad* pad1 = new TPad("pad1", "Top pad",    0, 0.32, 1, 1);
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

    double max_data  = h_data->GetMaximum();
    double max_stack = hs->GetMaximum();
    hs->SetMaximum(std::max(max_data, max_stack) * 2.0);
    hs->GetYaxis()->SetTitle("Events / 1 [GeV]");
    hs->GetYaxis()->SetTitleSize(0.05);
    hs->GetYaxis()->SetLabelSize(0.04);
    hs->GetYaxis()->SetTitleOffset(1.2);

    // Blinded data in 110–150 GeV
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
    TLegend* leg = new TLegend(0.60, 0.50, 0.88, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
    leg->AddEntry(h_signal, "VBF H#rightarrow#mu#mu", "f");
    leg->AddEntry(h_TT,     "t#bar{t} + Jets",       "f");
    leg->AddEntry(h_VVV,    "VVV + Jets",            "f");
    leg->AddEntry(h_VV,     "VV + Jets",             "f");
    leg->AddEntry(h_DY,     "DY + Jets",             "f");
    leg->AddEntry(h_data,   "Data",                  "lep");
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

    // MC total for ratio excludes signal = DY + TT + VV + VVV
    TH1* h_mc_total = (TH1*)h_DY->Clone("h_mc_total");
    h_mc_total->Add(h_TT);
    h_mc_total->Add(h_VV);
    h_mc_total->Add(h_VVV);

    TH1* h_ratio = (TH1*)h_data->Clone("h_ratio");
    h_ratio->Divide(h_mc_total);

    // Blind ratio in 110–150 GeV
    int bin_min = h_ratio->FindBin(110.0);
    int bin_max = h_ratio->FindBin(150.0);
    for (int i = bin_min; i <= bin_max; ++i) {
        h_ratio->SetBinContent(i, 0);
        h_ratio->SetBinError(i, 0);
    }

    // Ratio style
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
    h_ratio->Draw("E1");

    // -------- Save --------
    TString outdir = "/eos/user/d/dblend/StackedPlots/" + era + "/Muons/";
    gSystem->Exec("mkdir -p " + outdir); // ensure outdir exists
    c->SaveAs(outdir + "Stacked_DiMuonMass_blinded_withRatio.pdf");
    c->SaveAs(outdir + "Stacked_DiMuonMass_blinded_withRatio.png");

    // -------- Cleanup --------
    f_data->Close(); f_DY->Close(); f_TT->Close(); f_signal->Close();
    f_WZto2L2Q->Close(); f_WZto3LNu->Close(); f_ZZto2L2Q->Close(); f_ZZto2L2Nu->Close(); f_ZZto4L->Close();
    f_WWW->Close(); f_WWZ->Close(); f_WZZ->Close(); f_ZZZ->Close();
}

