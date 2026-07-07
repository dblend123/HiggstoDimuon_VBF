#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

static TH1* getHist(TFile* f, const std::string& path) {
    if (!f) return nullptr;
    auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
    if (!h) return nullptr;
    auto* c = (TH1*)h->Clone();
    c->SetDirectory(nullptr);
    return c;
}

static TFile* open_or_warn(const TString& p) {
    TFile* f = TFile::Open(p);
    if (!f) std::cerr << "[ERROR] Could not open: " << p << "\n";
    return f;
}


static std::string jetVarToBDTSystName(const std::string& var) {
    if (var == "Jet1Pt")         return "Jet1Pt";
    if (var == "Jet2Pt")         return "Jet2Pt";
    if (var == "Jet1eta")        return "Jet1Eta";
    if (var == "DiJetdeta")      return "DijetDeltaEta";
    if (var == "DiJetdphi")      return "DijetDeltaPhi";
    if (var == "hdphi_mumu_jj")  return "DeltaPhiDimuonDijet";
    if (var == "hdeta_mumu_jj")  return "DeltaEtaDimuonDijet";
    if (var == "nJets")          return "NJets";

    // These are not currently booked in BDTTheorySyst_VBF from your BDT input list.
    // If you want orange bands for them too, add matching SystHist1D objects in the analyzer.
    return "";
}

static TH1* makeVariationTotal(
    const std::string& outName,
    const std::string& systDir,
    const std::string& systVar,
    const std::string& suffix,
    const std::vector<std::pair<TFile*, TH1*>>& samples
) {
    TH1* total = nullptr;

    for (const auto& item : samples) {
        TFile* f = item.first;
        TH1* h_nom = item.second;
        if (!h_nom) continue;

        TH1* h_var = nullptr;
        if (!systVar.empty()) {
            h_var = getHist(f, systDir + systVar + suffix);
        }

        // If a variation histogram is missing, use nominal for this sample.
        // This gives zero theory deviation for that sample instead of crashing.
        TH1* h_to_add = h_var ? h_var : (TH1*)h_nom->Clone();
        h_to_add->SetDirectory(nullptr);

        if (!total) {
            total = (TH1*)h_to_add->Clone(outName.c_str());
            total->SetDirectory(nullptr);
        } else {
            total->Add(h_to_add);
        }

        delete h_to_add;
    }

    return total;
}

static TH1* makeOrangeTheoryBand(
    const std::string& name,
    TH1* h_nom,
    TH1* h_scaleUp, TH1* h_scaleDown,
    TH1* h_pdfUp,   TH1* h_pdfDown,
    TH1* h_isrUp,   TH1* h_isrDown,
    TH1* h_fsrUp,   TH1* h_fsrDown
) {
    if (!h_nom) return nullptr;

    TH1* band = (TH1*)h_nom->Clone(name.c_str());
    band->SetDirectory(nullptr);

    auto dev = [](TH1* h_up, TH1* h_down, TH1* h_nominal, int i) {
        const double nom = h_nominal->GetBinContent(i);
        const double up = h_up ? h_up->GetBinContent(i) : nom;
        const double dn = h_down ? h_down->GetBinContent(i) : nom;
        return std::max(std::abs(up - nom), std::abs(dn - nom));
    };

    for (int i = 1; i <= band->GetNbinsX(); ++i) {
        const double nom = h_nom->GetBinContent(i);

        const double dScale = dev(h_scaleUp, h_scaleDown, h_nom, i);
        const double dPDF   = dev(h_pdfUp,   h_pdfDown,   h_nom, i);
        const double dISR   = dev(h_isrUp,   h_isrDown,   h_nom, i);
        const double dFSR   = dev(h_fsrUp,   h_fsrDown,   h_nom, i);

        // Match the current muon plotting convention: theory band only,
        // with MC statistical uncertainty shown separately by h_mcBand.
        const double totalErr = std::sqrt(dScale*dScale + dPDF*dPDF + dISR*dISR + dFSR*dFSR);

        const double relErrRaw = (nom > 0.0 ? totalErr / nom : 0.0);
        const double relErr    = std::min(relErrRaw, 0.45);  // display cap

        band->SetBinContent(i, 1.0);
        band->SetBinError(i, relErr);
    }

    band->SetFillColor(kOrange-3);
    band->SetFillStyle(1001);
    band->SetMarkerSize(0);
    band->SetLineColor(kOrange+7);
    band->SetLineWidth(1);
    return band;
}

static TH1* makeAzureExperimentalBand(
    const std::string& name,
    TH1* h_nom,
    TH1* h_mcStatBand,
    TH1* h_puUp,     TH1* h_puDown,
    TH1* h_muonUp,   TH1* h_muonDown,
    TH1* h_trigUp,   TH1* h_trigDown,
    TH1* h_btagUp,   TH1* h_btagDown,
    TH1* h_zptUp,    TH1* h_zptDown
) {
    if (!h_nom) return nullptr;

    TH1* band = (TH1*)h_nom->Clone(name.c_str());
    band->SetDirectory(nullptr);

    auto dev = [](TH1* h_up, TH1* h_down, TH1* h_nominal, int i) {
        const double nom = h_nominal->GetBinContent(i);
        const double up = h_up ? h_up->GetBinContent(i) : nom;
        const double dn = h_down ? h_down->GetBinContent(i) : nom;
        return std::max(std::abs(up - nom), std::abs(dn - nom));
    };

    for (int i = 1; i <= band->GetNbinsX(); ++i) {
        const double nom = h_nom->GetBinContent(i);

        const double dPU   = dev(h_puUp,   h_puDown,   h_nom, i);
        const double dMuon = dev(h_muonUp, h_muonDown, h_nom, i);
        const double dTrig = dev(h_trigUp, h_trigDown, h_nom, i);
        const double dBTag = dev(h_btagUp, h_btagDown, h_nom, i);
        const double dZPt  = dev(h_zptUp,  h_zptDown,  h_nom, i);

        // Match the current muon plotting convention: experimental band only,
        // with MC statistical uncertainty shown separately by h_mcBand.
        const double totalErr = std::sqrt(dPU*dPU +
                                          dMuon*dMuon +
                                          dTrig*dTrig 
                                          //dBTag*dBTag +
                                          //dZPt*dZPt
                                          );

        const double relErrRaw = (nom > 0.0 ? totalErr / nom : 0.0);
        const double relErr    = std::min(relErrRaw, 0.45);  // display cap

        band->SetBinContent(i, 1.0);
        band->SetBinError(i, relErr);
    }

    band->SetFillColor(kAzure+1);
    band->SetFillStyle(1001);
    band->SetMarkerSize(0);
    band->SetLineColor(kAzure+1);
    return band;
}

static double getEraLumi_pb(const TString& era) {
    // Values (pb^-1) – kept consistent with Muon_Overlay_AllPlots.C
    if      (era == "12022") return  7980.454145;
    else if (era == "22022") return  26671.60971;
    else if (era == "12023") return  18704.13341;
    else if (era == "22023") return  9693.130053;
    else if (era == "2024")  return  104805.9585;
    else if (era == "2025")  return  114040.38;
    std::cerr << "[WARN] Unknown era \"" << era << "\"; returning 0 pb^-1\n";
    return 0.0;
}

void Jet_Overlay_AllPlots(TString basePath, TString era = "12022") {
    gStyle->SetOptStat(0);

    // Ensure trailing slash
    if (!basePath.EndsWith("/")) basePath += "/";

    // ----- Era Lumi (pb^-1 and fb^-1) -----
    const double eraLumi_pb = getEraLumi_pb(era);
    const double eraLumi_fb = eraLumi_pb / 1000.0; // convert pb^-1 -> fb^-1

    // This is the directory that actually holds DATA/, DY/, etc.
    const TString base = basePath;

    // --------------------------------------------------
    // Parse CorrType and ConfigName from basePath:
    //   .../<era>/<CorrType>/<ConfigName>/
    // e.g. CorrType = "KIT", ConfigName = "noPU_noScaleSmear_..."
    // --------------------------------------------------
    TString tmp = base;
    if (tmp.EndsWith("/")) tmp.Chop();          // strip trailing '/'
    Int_t last1 = tmp.Last('/');                // split off ConfigName
    TString configName = "UnknownConfig";
    TString upOne      = tmp;
    if (last1 >= 0) {
        configName = tmp(last1+1, tmp.Length() - last1 - 1);
        upOne      = tmp(0, last1);
    }
    Int_t last2 = upOne.Last('/');              // split off CorrType
    TString corrType = "UnknownCorr";
    if (last2 >= 0) {
        corrType = upOne(last2+1, upOne.Length() - last2 - 1);
    }
    // ---------- Open inputs ----------
TFile* f_data       = open_or_warn(base + "DATA/DATA_all.root");
TFile* f_DY         = open_or_warn(base + "DY/DY_all.root");
TFile* f_EWK        = open_or_warn(base + "EWK/EWK_all.root");
TFile* f_TT         = open_or_warn(base + "TT/TT_all.root");

// Diboson
TFile* f_WZto2L2Q   = open_or_warn(base + "WZto2L2Q/WZto2L2Q_all.root");
TFile* f_WZto3LNu   = open_or_warn(base + "WZto3LNu/WZto3LNu_all.root");
TFile* f_ZZto2L2Nu  = open_or_warn(base + "ZZto2L2Nu/ZZto2L2Nu_all.root");
TFile* f_ZZto2L2Q   = open_or_warn(base + "ZZto2L2Q/ZZto2L2Q_all.root");
TFile* f_ZZto4L     = open_or_warn(base + "ZZto4L/ZZto4L_all.root");

// Triboson
TFile* f_WWW        = open_or_warn(base + "WWW/WWW_all.root");
TFile* f_WWZ        = open_or_warn(base + "WWZ/WWZ_all.root");
TFile* f_WZZ        = open_or_warn(base + "WZZ/WZZ_all.root");
TFile* f_ZZZ        = open_or_warn(base + "ZZZ/ZZZ_all.root");

// Signal
TFile* f_signal     = open_or_warn(base + "VBF/VBF_all.root");

    if (!f_data || !f_DY || !f_EWK || !f_TT ||
        !f_WZto2L2Q || !f_WZto3LNu || !f_ZZto2L2Nu || !f_ZZto2L2Q || !f_ZZto4L ||
        !f_WWW || !f_WWZ || !f_WZZ || !f_ZZZ ||
        !f_signal) {
        std::cerr << "[ERROR] One or more files couldn't be opened!" << std::endl;
        return;
    }

    // ---------- Vars to overlay ----------
    std::vector<std::string> vars = {
        "DiJetPt", "Jet1Pt", "Jet2Pt",
        "DiJetdeta", "DiJetdphi", "Jet1eta", 
        "Jet2eta", "DiJetMass", "hdphi_mumu_jj", 
        "hdeta_mumu_jj", "nJets"
    };
    const std::string dir = "corr_Jet_VBF_SR/";
    
    // ---------- Output directory + summary file ----------
TString outdir = "/eos/user/d/dblend/StackedPlots/" + era + "/Jets/" + corrType + "/" + configName + "/";
gSystem->Exec("mkdir -p " + outdir);
std::ofstream plotLog((outdir + "OverlayRatios_All.txt").Data(), std::ios::trunc);
if (!plotLog) {
    std::cerr << "[ERROR] Could not open " << (outdir + "OverlayRatios_All.txt") << " for writing\n";
    return;
}

for (const auto& var : vars) {
        // Pull histograms
        TH1* h_data   = getHist(f_data,   dir + var);
        TH1* h_DY     = getHist(f_DY,     dir + var);
        TH1* h_EWK    = getHist(f_EWK,    dir + var);
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

        if (!h_data || !h_DY || !h_EWK || !h_TT ||
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

        // ---------- Totals for logging ----------
        const double y_data = h_data->Integral(1, h_data->GetNbinsX());
        const double y_DY   = h_DY->Integral(1, h_DY->GetNbinsX());
        const double y_EWK  = h_EWK->Integral(1, h_EWK->GetNbinsX());
        const double y_TT   = h_TT->Integral(1, h_TT->GetNbinsX());
        const double y_VV   = h_VV->Integral(1, h_VV->GetNbinsX());
        const double y_VVV  = h_VVV->Integral(1, h_VVV->GetNbinsX());
        const double y_mc   = y_DY + y_EWK + y_TT + y_VV + y_VVV;

        const double r_data_mc = (y_mc > 0.0 ? y_data / y_mc : 0.0);
        const double eff_data_lumi_fb = r_data_mc * eraLumi_fb;

        const char* xTitle =
            (var=="DiJetPt")       ? "p_{T}(jj) [GeV]" :
            (var=="Jet1Pt")        ? "p_{T}(j_{1}) [GeV]" :
            (var=="Jet2Pt")        ? "p_{T}(j_{2}) [GeV]" :
            (var=="DiJetdeta")     ? "#Delta#eta(jj)" :
            (var=="DiJetdphi")     ? "#Delta#phi(jj)" :
            (var=="Jet1eta")       ? "#eta(j_{1})" :
            (var=="Jet2eta")       ? "#eta(j_{2})" :
            (var=="DiJetMass")     ? "m_{jj} [GeV]" :
            (var=="hdphi_mumu_jj") ? "#Delta#phi(#mu#mu,jj)" :
            (var=="hdeta_mumu_jj") ? "#Delta#eta(#mu#mu,jj)" :
            (var=="nJets")         ? "N_{jets}" :
                              var.c_str();

        // Emit summary lines like your example
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "[PLOT] dir=" << dir
                  << " var=" << var
                  << " data=" << y_data
                  << " DY="   << y_DY
                  << " EWK="  << y_EWK
                  << " TT="   << y_TT
                  << " VV="   << y_VV
                  << " VVV="  << y_VVV
                  << " mc="   << y_mc
                  << " data/mc=" << r_data_mc
                  << std::endl;

        std::cout << "[PLOT] eff_data_lumi \xE2\x89\x88 "  // "≈" in UTF-8
                  << eff_data_lumi_fb << " fb^-1 (using eraLumi="
                  << eraLumi_fb << " fb^-1)"
                  << std::endl;
        
        // Also write the same info into the log file
        plotLog << std::fixed << std::setprecision(6);
        plotLog << "Variable: " << var << "\n";
        plotLog << "Data: " << y_data << "\n";
        plotLog << "MC: "   << y_mc
                << "  (DY=" << y_DY
                << ", EWK=" << y_EWK
                << ", TT=" << y_TT
                << ", VV=" << y_VV
                << ", VVV=" << y_VVV << ")\n";
        plotLog << "Data/MC: " << r_data_mc << "\n";
        plotLog << "Era lumi (fb^-1): " << eraLumi_fb << "\n";
        plotLog << "Effective data lumi (fb^-1): " << eff_data_lumi_fb << "\n";
        plotLog << "------------------------------------------\n";

        // ---------- Style ----------
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(1.2);
        h_data->SetMarkerColor(kBlack);
        h_data->SetLineColor(kBlack);

        h_DY->SetFillColor(kOrange-2); h_DY->SetLineColor(kOrange+3);
        h_TT->SetFillColor(kRed-4);    h_TT->SetLineColor(kRed+1);
        h_VV->SetFillColor(kAzure-9);  h_VV->SetLineColor(kAzure+2);
        h_VVV->SetFillColor(kGreen-7); h_VVV->SetLineColor(kGreen+3);
        h_EWK->SetFillColor(kGray+1);  h_EWK->SetLineColor(kGray+3);
        h_signal->SetFillColor(kMagenta+1); h_signal->SetLineColor(kMagenta+3); h_signal->SetLineWidth(3);

        // ---------- Canvas with ratio pad ----------
        TCanvas* c = new TCanvas(("c_"+var).c_str(), var.c_str(), 850, 800);

        TPad* pad1 = new TPad(("pad1_"+var).c_str(), "Top",    0, 0.30, 1, 1);
        TPad* pad2 = new TPad(("pad2_"+var).c_str(), "Bottom", 0, 0.02, 1, 0.30);

        pad1->SetLeftMargin(0.14);
        pad1->SetRightMargin(0.04);
        pad1->SetBottomMargin(0.02);

        pad2->SetLeftMargin(0.14);
        pad2->SetRightMargin(0.04);
        pad2->SetTopMargin(0.02);
        pad2->SetBottomMargin(0.32);
        
        pad1->Draw(); 
        pad2->Draw();

        // ---------- Top pad ----------
        pad1->cd();
        pad1->SetLogy();

        // ---- Build a stacked MC (backgrounds only) ----
        THStack* hs = new THStack(("hs_"+var).c_str(), "");

// Add BIGGEST FIRST so it sits on the bottom of the stack
hs->Add(h_EWK);
hs->Add(h_VVV);
hs->Add(h_VV);
hs->Add(h_TT);
hs->Add(h_DY);

// Empty frame owns the axes
TH1* h_frame = (TH1*)h_data->Clone(("h_frame_"+var).c_str());
h_frame->Reset("ICES");
h_frame->SetTitle("");
h_frame->SetMinimum(0.1);
h_frame->SetMaximum(1.0e6);

h_frame->GetYaxis()->SetTitle("Events");
h_frame->GetYaxis()->SetTitleSize(0.055);
h_frame->GetYaxis()->SetTitleOffset(1.05);
h_frame->GetYaxis()->SetLabelSize(0.040);

h_frame->GetXaxis()->SetLabelSize(0);
h_frame->GetXaxis()->SetTitleSize(0);

h_frame->SetStats(0);
h_frame->SetLineColor(0);
h_frame->SetMarkerColor(0);
h_frame->SetFillStyle(0);

h_frame->Draw("HIST");

hs->Draw("HIST SAME");
h_signal->Draw("HIST SAME");
h_data->Draw("E1 SAME");

gPad->RedrawAxis();

        // Legend
TLegend* leg = new TLegend(0.62, 0.74, 0.97, 0.89);

leg->SetNColumns(2);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->SetTextFont(42);
leg->SetTextSize(0.026);      // bigger text
leg->SetMargin(0.12);
leg->SetEntrySeparation(0.01);

        leg->AddEntry(h_data,   "Data",       "lep");
        leg->AddEntry(h_DY,     "DY",         "f");
        leg->AddEntry(h_EWK,    "EWK",        "f");

        leg->AddEntry(h_TT,     "Top",        "f");
        leg->AddEntry(h_VV,     "Diboson",    "f");
        leg->AddEntry(h_VVV,    "Other Bkg.", "f");

        leg->AddEntry(h_signal, "VBF",        "f");

        leg->Draw();

        // CMS prelim + lumi
        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(62);
        latex.DrawLatex(0.12, 0.93, "CMS");

        latex.SetTextFont(52); // italic
        latex.DrawLatex(0.19, 0.93, "Preliminary");

        latex.SetTextFont(62); // bold
        latex.DrawLatex(0.62, 0.93, Form("%.2f fb^{-1} (13.6 TeV)", eraLumi_pb/1000.0));

        // ---------- Bottom pad: ratio ----------
        pad2->cd();

        // MC total for ratio excludes signal
        TH1* h_mctot = (TH1*)h_DY->Clone(("h_mctot_"+var).c_str());
        h_mctot->Add(h_EWK);
        h_mctot->Add(h_TT);
        h_mctot->Add(h_VV);
        h_mctot->Add(h_VVV);

        // ---------- Orange theory band: MC stat ⊕ scale/PDF/ISR/FSR ----------
        const std::string systDir = "BDTTheorySyst_VBF/";
        const std::string systVar = jetVarToBDTSystName(var);

        std::vector<std::pair<TFile*, TH1*>> bkgSamples = {
            {f_DY, h_DY},
            {f_EWK, h_EWK},
            {f_TT, h_TT},
            {f_WZto2L2Q, h_WZ2L2Q},
            {f_WZto3LNu, h_WZ3L},
            {f_ZZto2L2Nu, h_ZZ2L2N},
            {f_ZZto2L2Q, h_ZZ2L2Q},
            {f_ZZto4L, h_ZZ4L},
            {f_WWW, h_WWW},
            {f_WWZ, h_WWZ},
            {f_WZZ, h_WZZ},
            {f_ZZZ, h_ZZZ}
        };

        TH1* h_orangeBand = nullptr;
        TH1* h_scaleUpTot = nullptr;
        TH1* h_scaleDownTot = nullptr;
        TH1* h_pdfUpTot = nullptr;
        TH1* h_pdfDownTot = nullptr;
        TH1* h_isrUpTot = nullptr;
        TH1* h_isrDownTot = nullptr;
        TH1* h_fsrUpTot = nullptr;
        TH1* h_fsrDownTot = nullptr;

        if (!systVar.empty()) {
            h_scaleUpTot   = makeVariationTotal("h_scaleUpTot_"   + var, systDir, systVar, "_scaleUp",   bkgSamples);
            h_scaleDownTot = makeVariationTotal("h_scaleDownTot_" + var, systDir, systVar, "_scaleDown", bkgSamples);
            h_pdfUpTot     = makeVariationTotal("h_pdfUpTot_"     + var, systDir, systVar, "_pdfUp",     bkgSamples);
            h_pdfDownTot   = makeVariationTotal("h_pdfDownTot_"   + var, systDir, systVar, "_pdfDown",   bkgSamples);
            h_isrUpTot     = makeVariationTotal("h_isrUpTot_"     + var, systDir, systVar, "_isrUp",     bkgSamples);
            h_isrDownTot   = makeVariationTotal("h_isrDownTot_"   + var, systDir, systVar, "_isrDown",   bkgSamples);
            h_fsrUpTot     = makeVariationTotal("h_fsrUpTot_"     + var, systDir, systVar, "_fsrUp",     bkgSamples);
            h_fsrDownTot   = makeVariationTotal("h_fsrDownTot_"   + var, systDir, systVar, "_fsrDown",   bkgSamples);

            h_orangeBand = makeOrangeTheoryBand("h_orangeTheoryBand_" + var, h_mctot,
                                                h_scaleUpTot, h_scaleDownTot,
                                                h_pdfUpTot,   h_pdfDownTot,
                                                h_isrUpTot,   h_isrDownTot,
                                                h_fsrUpTot,   h_fsrDownTot);
        } else {
            std::cerr << "[WARN] No BDTTheorySyst_VBF mapping for jet variable " << var
                      << "; drawing gray MC-stat band only.\n";
        }

TH1* h_mcBand = (TH1*)h_mctot->Clone(("h_mcBand_"+var).c_str());
h_mcBand->SetDirectory(nullptr);

for (int i = 1; i <= h_mcBand->GetNbinsX(); ++i) {
    double mc  = h_mctot->GetBinContent(i);
    double err = h_mctot->GetBinError(i);

    if (mc > 0) {
        h_mcBand->SetBinContent(i, 1.0);
        h_mcBand->SetBinError(i, err / mc);
    } else {
        h_mcBand->SetBinContent(i, 1.0);
        h_mcBand->SetBinError(i, 0.0);
    }
}

h_mcBand->SetFillColor(kGray+2);
h_mcBand->SetFillStyle(1001);
h_mcBand->SetMarkerSize(0);
h_mcBand->SetLineColor(kGray+2);

        // ---------- Azure experimental band: MC stat ⊕ PU/muon/trigger/btag/ZpT ----------
        const std::string expSystDir = "BDTExperimentalSyst_VBF/";

        TH1* h_azureBand = nullptr;
        TH1* h_puUpTot = nullptr;
        TH1* h_puDownTot = nullptr;
        TH1* h_muonSFUpTot = nullptr;
        TH1* h_muonSFDownTot = nullptr;
        TH1* h_trigUpTot = nullptr;
        TH1* h_trigDownTot = nullptr;
        TH1* h_btagUpTot = nullptr;
        TH1* h_btagDownTot = nullptr;
        TH1* h_zptUpTot = nullptr;
        TH1* h_zptDownTot = nullptr;

        if (!systVar.empty()) {
            h_puUpTot       = makeVariationTotal("h_puUpTot_"       + var, expSystDir, systVar, "_puUp",       bkgSamples);
            h_puDownTot     = makeVariationTotal("h_puDownTot_"     + var, expSystDir, systVar, "_puDown",     bkgSamples);
            h_muonSFUpTot   = makeVariationTotal("h_muonSFUpTot_"   + var, expSystDir, systVar, "_muonSFUp",   bkgSamples);
            h_muonSFDownTot = makeVariationTotal("h_muonSFDownTot_" + var, expSystDir, systVar, "_muonSFDown", bkgSamples);
            h_trigUpTot     = makeVariationTotal("h_trigUpTot_"     + var, expSystDir, systVar, "_triggerUp",     bkgSamples);
            h_trigDownTot   = makeVariationTotal("h_trigDownTot_"   + var, expSystDir, systVar, "_triggerDown",   bkgSamples);
            h_btagUpTot     = makeVariationTotal("h_btagUpTot_"     + var, expSystDir, systVar, "_btagUp",     bkgSamples);
            h_btagDownTot   = makeVariationTotal("h_btagDownTot_"   + var, expSystDir, systVar, "_btagDown",   bkgSamples);
            h_zptUpTot      = makeVariationTotal("h_zptUpTot_"      + var, expSystDir, systVar, "_zptUp",      bkgSamples);
            h_zptDownTot    = makeVariationTotal("h_zptDownTot_"    + var, expSystDir, systVar, "_zptDown",    bkgSamples);

            h_azureBand = makeAzureExperimentalBand("h_azureExperimentalBand_" + var, h_mctot, h_mcBand,
                                                    h_puUpTot,       h_puDownTot,
                                                    h_muonSFUpTot,   h_muonSFDownTot,
                                                    h_trigUpTot,     h_trigDownTot,
                                                    h_btagUpTot,     h_btagDownTot,
                                                    h_zptUpTot,      h_zptDownTot);
        }

        TH1* h_ratio = (TH1*)h_data->Clone(("h_ratio_"+var).c_str());
        h_ratio->SetTitle("");
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
        h_ratio->GetXaxis()->SetTitle(xTitle);
        h_ratio->SetMinimum(0.5);
        h_ratio->SetMaximum(1.5);
        h_ratio->SetMarkerStyle(20);
        h_ratio->SetLineColor(kBlack);
        h_ratio->SetStats(0);
        h_ratio->Draw("AXIS");

        // Match Muon_Overlay_AllPlots aesthetics exactly:
        // draw widest/colored bands first, then gray MC band, then data points.
        if (h_azureBand)  h_azureBand->Draw("E2 SAME");
        if (h_orangeBand) h_orangeBand->Draw("E2 SAME");
        h_mcBand->Draw("E2 SAME");
        h_ratio->Draw("E1 SAME");

        // Horizontal dotted guide lines every 0.2 in Data/MC
        pad2->Update();
        {
            double rxmin = gPad->GetUxmin();
            double rxmax = gPad->GetUxmax();

            for (double y = 0.6; y <= 1.4001; y += 0.2) {
                TLine* line = new TLine(rxmin, y, rxmax, y);
                line->SetLineStyle(3);      // dotted
                line->SetLineWidth(1);
                line->SetLineColor(kBlack);
                line->Draw("SAME");
            }
        }

        // ---------- Save ----------
TString baseName = "OverlayRatio_" + TString(var.c_str());
c->SaveAs(outdir + baseName + ".pdf");
c->SaveAs(outdir + baseName + ".png");
        
        delete h_frame;
        delete c;
        delete h_mcBand;
        delete h_orangeBand;
        delete h_scaleUpTot;
        delete h_scaleDownTot;
        delete h_pdfUpTot;
        delete h_pdfDownTot;
        delete h_isrUpTot;
        delete h_isrDownTot;
        delete h_fsrUpTot;
        delete h_fsrDownTot;
        delete h_azureBand;
        delete h_puUpTot;
        delete h_puDownTot;
        delete h_muonSFUpTot;
        delete h_muonSFDownTot;
        delete h_trigUpTot;
        delete h_trigDownTot;
        delete h_btagUpTot;
        delete h_btagDownTot;
        delete h_zptUpTot;
        delete h_zptDownTot;
    }

    // ---------- Cleanup ----------
    f_data->Close(); f_DY->Close(); f_EWK->Close(); f_TT->Close(); f_signal->Close();
    f_WZto2L2Q->Close(); f_WZto3LNu->Close(); f_ZZto2L2Q->Close(); f_ZZto2L2Nu->Close(); f_ZZto4L->Close();
    f_WWW->Close(); f_WWZ->Close(); f_WZZ->Close(); f_ZZZ->Close(); plotLog.close();
}
