
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

static TH1* getHist(TFile* f, const std::string& path) {
    if (!f) return nullptr;
    auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
    return h ? (TH1*)h->Clone() : nullptr;
}

void Jet_Overlay(TString era = "12022") {
    gStyle->SetOptStat(0);

    const TString base = "/eos/user/d/dblend/Analyzed_SkimmedFiles/" + era + "/";

    // ---------- Open inputs ----------
    TFile* f_data       = TFile::Open(base + "MC/DATA/DATA_all.root");
    TFile* f_DY         = TFile::Open(base + "MC/DY/DY_all.root");
    TFile* f_TT         = TFile::Open(base + "MC/TT/TT_all.root");

    // Diboson
    TFile* f_WZto2L2Q   = TFile::Open(base + "MC/WZto2L2Q/WZto2L2Q_all.root");
    TFile* f_WZto3LNu   = TFile::Open(base + "MC/WZto3LNu/WZto3LNu_all.root");
    TFile* f_ZZto2L2Nu  = TFile::Open(base + "MC/ZZto2L2Nu/ZZto2L2Nu_all.root");
    TFile* f_ZZto2L2Q   = TFile::Open(base + "MC/ZZto2L2Q/ZZto2L2Q_all.root");
    TFile* f_ZZto4L     = TFile::Open(base + "MC/ZZto4L/ZZto4L_all.root");

    // Triboson
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

    // ---------- Vars to overlay ----------
    std::vector<std::string> vars = {
        "DiJetPt", "Jet1Pt", "Jet2Pt",
        "DiJetdeta", "DiJetdphi", "Jet1eta", "Jet2eta", "DiJetMass"
    };
    const std::string dir = "corr_Jet_PreDiJet_VBF/";

    for (const auto& var : vars) {
        // Pull histograms
        TH1* h_data   = getHist(f_data,   dir + var);
        TH1* h_DY     = getHist(f_DY,     dir + var);
        TH1* h_TT     = getHist(f_TT,     dir + var);

        // Diboson
        TH1* h_WZ2L2Q = getHist(f_WZto2L2Q,  dir + var);
        TH1* h_WZ3L   = getHist(f_WZto3LNu,  dir + var);
        TH1* h_ZZ2L2Q = getHist(f_ZZto2L2Q,  dir + var);
        TH1* h_ZZ2L2N = getHist(f_ZZto2L2Nu, dir + var);
        TH1* h_ZZ4L   = getHist(f_ZZto4L,    dir + var);

        // Triboson
        TH1* h_WWW    = getHist(f_WWW,       dir + var);
        TH1* h_WWZ    = getHist(f_WWZ,       dir + var);
        TH1* h_WZZ    = getHist(f_WZZ,       dir + var);
        TH1* h_ZZZ    = getHist(f_ZZZ,       dir + var);

        // Signal
        TH1* h_signal = getHist(f_signal, dir + var);

        if (!h_data || !h_DY || !h_TT ||
            !h_WZ2L2Q || !h_WZ3L || !h_ZZ2L2Q || !h_ZZ2L2N || !h_ZZ4L ||
            !h_WWW || !h_WWZ || !h_WZZ || !h_ZZZ ||
            !h_signal) {
            std::cerr << "[WARN] Missing hist(s) for " << var << " — skipping.\n";
            continue;
        }

        // Group VV
        TH1* h_VV = (TH1*)h_WZ2L2Q->Clone(("h_VV_"+var).c_str());
        h_VV->Add(h_WZ3L);
        h_VV->Add(h_ZZ2L2Q);
        h_VV->Add(h_ZZ2L2N);
        h_VV->Add(h_ZZ4L);

        // Group VVV
        TH1* h_VVV = (TH1*)h_WWW->Clone(("h_VVV_"+var).c_str());
        h_VVV->Add(h_WWZ);
        h_VVV->Add(h_WZZ);
        h_VVV->Add(h_ZZZ);

        // ---------- Style ----------
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(1.2);
        h_data->SetMarkerColor(kBlack);
        h_data->SetLineColor(kBlack);

        h_DY->SetFillColor(kOrange-2); h_DY->SetLineColor(kOrange+3); h_DY->SetLineWidth(2);
        h_TT->SetFillColor(kRed-4);    h_TT->SetLineColor(kRed+1);    h_TT->SetLineWidth(2);
        h_VV->SetFillColor(kAzure-9);  h_VV->SetLineColor(kAzure+2);  h_VV->SetLineWidth(2);
        h_VVV->SetFillColor(kGreen-7); h_VVV->SetLineColor(kGreen+3); h_VVV->SetLineWidth(2);
        h_signal->SetFillColor(kMagenta+1); h_signal->SetLineColor(kMagenta+3); h_signal->SetLineWidth(3);

        // ---------- Canvas with ratio pad ----------
        TCanvas* c = new TCanvas(("c_"+var).c_str(), var.c_str(), 800, 800);
        TPad* pad1 = new TPad(("pad1_"+var).c_str(), "Top pad",    0, 0.32, 1, 1);
        TPad* pad2 = new TPad(("pad2_"+var).c_str(), "Ratio pad",  0, 0.02, 1, 0.32);
        pad1->SetBottomMargin(0.03);
        pad2->SetTopMargin(0.03);
        pad2->SetBottomMargin(0.30);
        pad1->Draw(); pad2->Draw();

        // ---------- Top pad ----------
        pad1->cd();
        pad1->SetLogy();

        double ymax = std::max({
            h_data->GetMaximum(),
            h_DY->GetMaximum(),
            h_TT->GetMaximum(),
            h_VV->GetMaximum(),
            h_VVV->GetMaximum(),
            h_signal->GetMaximum()
        }) * 1.6;

        // Use data as frame
        h_data->SetTitle("");
        h_data->GetYaxis()->SetTitle("Events");
        h_data->GetYaxis()->SetTitleSize(0.05);
        h_data->GetYaxis()->SetLabelSize(0.04);
        h_data->GetYaxis()->SetTitleOffset(1.2);
        h_data->SetMaximum(ymax);
        h_data->SetMinimum(0.1);
        h_data->Draw("E1");

        // Draw MC
        h_DY->Draw("HIST SAME");
        h_TT->Draw("HIST SAME");
        h_VV->Draw("HIST SAME");
        h_VVV->Draw("HIST SAME");
        h_signal->Draw("HIST SAME");
        h_data->Draw("E1 SAME");

        // Legend
        TLegend* leg = new TLegend(0.60, 0.50, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(h_signal, "VBF H#rightarrow#mu#mu", "f");
        leg->AddEntry(h_TT,     "t#bar{t} + Jets",        "f");
        leg->AddEntry(h_DY,     "DY + Jets",              "f");
        leg->AddEntry(h_VV,     "VV + Jets",              "f");
        leg->AddEntry(h_VVV,    "VVV + Jets",             "f");
        leg->AddEntry(h_data,   "Data",                   "lep");
        leg->Draw();

        // CMS prelim + lumi
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.13, 0.93, "CMS #scale[0.8]{#it{Preliminary}}");
        latex.DrawLatex(0.55, 0.93, "#sqrt{s} = 13 TeV, 7.98 fb^{-1}");

        // ---------- Bottom pad: ratio ----------
        pad2->cd();

        // MC total for ratio excludes signal
        TH1* h_mctot = (TH1*)h_DY->Clone(("h_mctot_"+var).c_str());
        h_mctot->Add(h_TT);
        h_mctot->Add(h_VV);
        h_mctot->Add(h_VVV);

        TH1* h_ratio = (TH1*)h_data->Clone(("h_ratio_"+var).c_str());
        for (int i = 1; i <= h_ratio->GetNbinsX(); ++i) {
            if (h_mctot->GetBinContent(i) <= 0) {
                h_ratio->SetBinContent(i, 0);
                h_ratio->SetBinError(i, 0);
            }
        }
        h_ratio->Divide(h_mctot);

        h_ratio->GetYaxis()->SetTitle("Data/MC");
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.4);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetXaxis()->SetLabelSize(0.12);
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->SetMinimum(0.0);
        h_ratio->SetMaximum(2.0);
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetLineColor(kBlack);
        h_ratio->Draw("E1");

        // Horizontal line at 1
        TLine* line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0,
                                h_ratio->GetXaxis()->GetXmax(), 1.0);
        line->SetLineStyle(2);
        line->Draw("SAME");

        // ---------- Save ----------
        TString outdir = "/eos/user/d/dblend/StackedPlots/" + era + "/Jets/";
        gSystem->Exec("mkdir -p " + outdir);
        TString baseName = "OverlayRatio_" + TString(var.c_str());
        c->SaveAs(outdir + baseName + ".pdf");
        c->SaveAs(outdir + baseName + ".png");

        delete c;
    }

    // ---------- Cleanup ----------
    f_data->Close(); f_DY->Close(); f_TT->Close(); f_signal->Close();
    f_WZto2L2Q->Close(); f_WZto3LNu->Close(); f_ZZto2L2Q->Close(); f_ZZto2L2Nu->Close(); f_ZZto4L->Close();
    f_WWW->Close(); f_WWZ->Close(); f_WZZ->Close(); f_ZZZ->Close();
}
