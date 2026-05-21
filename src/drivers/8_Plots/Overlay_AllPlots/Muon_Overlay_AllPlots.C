// ================== Muon_Overlay.C (clean two-mode version) ==================
// Muon_Overlay_AllPlots("/eos/user/d/dblend/Analyzed_SkimmedFiles/2024/KIT/noPU_noScaleSmear_noMuonSF_noJERC_noJetVeto_noGolden_noFlags_noBtagSF_noJetCriteria", "2024");
//

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <memory>
#include <algorithm>

// ---------- helpers ----------
static TH1* getHist(TFile* f, const std::string& path) {
  if (!f) return nullptr;
  auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
  if (!h) return nullptr;
  auto* c = (TH1*)h->Clone();
  c->SetDirectory(nullptr);
  return c;
}

static double integralExcluding(const TH1* h, double xmin, double xmax) {
  if (!h) return 0.0;
  const int nb = h->GetNbinsX();
  double sum = 0.0;
  for (int i = 1; i <= nb; ++i) {
    const double lo = h->GetXaxis()->GetBinLowEdge(i);
    const double hi = h->GetXaxis()->GetBinUpEdge(i);
    const bool overlap = (hi > xmin && lo < xmax);
    if (!overlap) sum += h->GetBinContent(i);
  }
  return sum;
}

static TFile* open_or_warn(const TString& p) {
  TFile* f = TFile::Open(p);
  if (!f) std::cerr << "[ERROR] Could not open: " << p << "\n";
  return f;
}

// ---------- draw one variable ----------
static void draw_one(const TString& outdir, const std::string& dir, const std::string& var,
                     TFile* f_data, TFile* f_DY, TFile* f_EWK, TFile* f_TT,
                     TFile* f_WZ2L2Q, TFile* f_WZ3LNu, TFile* f_ZZ2L2Nu, TFile* f_ZZ2L2Q, TFile* f_ZZ4L,
                     TFile* f_WWW, TFile* f_WWZ, TFile* f_WZZ, TFile* f_ZZZ,
                     TFile* f_signal, const char* xTitle, bool doBlindMass, double eraLumiPb,
                     std::ostream* summary = nullptr)
{
  const std::string hpath = dir + var;

  // read
  TH1* h_data   = getHist(f_data,   hpath);
  TH1* h_DY     = getHist(f_DY,     hpath);
  TH1* h_EWK    = getHist(f_EWK,    hpath);
  TH1* h_TT     = getHist(f_TT,     hpath);
  TH1* h_WZ2L2Q = getHist(f_WZ2L2Q, hpath);
  TH1* h_WZ3L   = getHist(f_WZ3LNu, hpath);
  TH1* h_ZZ2L2N = getHist(f_ZZ2L2Nu,hpath);
  TH1* h_ZZ2L2Q = getHist(f_ZZ2L2Q, hpath);
  TH1* h_ZZ4L   = getHist(f_ZZ4L,   hpath);
  TH1* h_WWW    = getHist(f_WWW,    hpath);
  TH1* h_WWZ    = getHist(f_WWZ,    hpath);
  TH1* h_WZZ    = getHist(f_WZZ,    hpath);
  TH1* h_ZZZ    = getHist(f_ZZZ,    hpath);
  TH1* h_signal = getHist(f_signal, hpath);

  if (!h_data || !h_DY || !h_TT || !h_EWK ||
      !h_WZ2L2Q || !h_WZ3L || !h_ZZ2L2N || !h_ZZ2L2Q || !h_ZZ4L ||
      !h_WWW || !h_WWZ || !h_WZZ || !h_ZZZ || !h_signal) {
    std::cerr << "[WARN] Missing hist(s) for " << var << " — skipping.\n";
    // clean
    delete h_data; delete h_DY; delete h_EWK; delete h_TT; delete h_WZ2L2Q; delete h_WZ3L;
    delete h_ZZ2L2N; delete h_ZZ2L2Q; delete h_ZZ4L; delete h_WWW; delete h_WWZ; delete h_WZZ; delete h_ZZZ; delete h_signal;
    return;
  }

  // groups
  TH1* h_VV  = (TH1*)h_WZ2L2Q->Clone(("h_VV_"+var).c_str());
  h_VV->SetDirectory(nullptr);
  h_VV->Add(h_WZ3L); h_VV->Add(h_ZZ2L2Q); h_VV->Add(h_ZZ2L2N); h_VV->Add(h_ZZ4L);

  TH1* h_VVV = (TH1*)h_WWW->Clone(("h_VVV_"+var).c_str());
  h_VVV->SetDirectory(nullptr);
  h_VVV->Add(h_WWZ); h_VVV->Add(h_WZZ); h_VVV->Add(h_ZZZ);

  for (TH1* h : {h_data,h_DY,h_EWK,h_TT,h_VV,h_VVV,h_signal})
    if (h && h->GetSumw2N()==0) h->Sumw2(kTRUE);

  // style
  h_data->SetMarkerStyle(20); h_data->SetMarkerColor(kBlack); h_data->SetLineColor(kBlack); h_data->SetMarkerSize(1.2);
  h_DY->SetFillColor(kOrange-2); h_DY->SetLineColor(kOrange+3);
  h_TT->SetFillColor(kRed-4);    h_TT->SetLineColor(kRed+1);
  h_VV->SetFillColor(kAzure-9);  h_VV->SetLineColor(kAzure+2);
  h_VVV->SetFillColor(kGreen-7); h_VVV->SetLineColor(kGreen+3);
  h_EWK->SetFillColor(kGray+1);  h_EWK->SetLineColor(kGray+3);
  h_signal->SetFillColor(kMagenta+1); h_signal->SetLineColor(kMagenta+3); h_signal->SetLineWidth(3);

  // canvas/pads
  std::unique_ptr<TCanvas> c(new TCanvas(("c_"+var).c_str(), var.c_str(), 800, 800));
  std::unique_ptr<TPad> pad1(new TPad(("pad1_"+var).c_str(), "Top",    0, 0.32, 1, 1));
  std::unique_ptr<TPad> pad2(new TPad(("pad2_"+var).c_str(), "Bottom", 0, 0.02, 1, 0.32));
  pad1->SetBottomMargin(0.03); pad2->SetTopMargin(0.03); pad2->SetBottomMargin(0.3);
  pad1->Draw(); pad2->Draw();

  // top
  pad1->cd(); 
  pad1->SetLogy();

  TH1* h_data_draw = (TH1*)h_data->Clone(("h_data_draw_"+TString(var.c_str())).Data());
  h_data_draw->SetDirectory(nullptr);

  // Blinding windows:
  //   - Inclusive mass plot: hide the full Higgs search region [110,150] GeV.
  //   - Dedicated Higgs plot: hide only the signal window [115,135] GeV,
  //     leaving the sidebands [110,115) and (135,150] visible.
  const bool blind = (doBlindMass && (var=="DiMuonMass" || var=="DiMuonMass_Higgs"));
  const double blindMin = (var=="DiMuonMass_Higgs" ? 115.0 : 110.0);
  const double blindMax = (var=="DiMuonMass_Higgs" ? 135.0 : 150.0);

  if (blind) {
    const int bmin = h_data_draw->FindBin(blindMin);
    const int bmax = h_data_draw->FindBin(blindMax);
    for (int i=bmin;i<=bmax;++i){ h_data_draw->SetBinContent(i,0); h_data_draw->SetBinError(i,0); }
  }

  // integrals (exclude blinded window for both Data/MC when blinding is active)
  auto iNoUF = [](TH1* h){ return h ? h->Integral(1,h->GetNbinsX()) : 0.0; };
  auto iMC   = [&](TH1* h){ return blind ? integralExcluding(h,blindMin,blindMax) : iNoUF(h); };

  const double data_int = blind ? integralExcluding(h_data,blindMin,blindMax) : iNoUF(h_data);
  const double dy_int   = iMC(h_DY);
  const double ewk_int  = iMC(h_EWK);
  const double tt_int   = iMC(h_TT);
  const double vv_int   = iMC(h_VV);
  const double vvv_int  = iMC(h_VVV);
  const double mc_int   = dy_int + ewk_int + tt_int + vv_int + vvv_int; // include EWK!

  const double dataOverMc = (mc_int>0.0 ? data_int/mc_int : 0.0);
  const double effLumiPb  = dataOverMc * eraLumiPb;

  std::cout << std::fixed
            << "[PLOT] dir=" << dir << " var=" << var
            << " Data=" << data_int << " DY=" << dy_int << " EWK=" << ewk_int
            << " TT=" << tt_int << " VV=" << vv_int << " VVV=" << vvv_int
            << " MC=" << mc_int << " Data/MC=" << dataOverMc << "\n";

  if (summary && *summary) {
    (*summary) << std::fixed << std::setprecision(6)
               << "Variable: " << var << "\n"
               << "Blinded: " << (blind ? Form("yes [%.0f,%.0f] GeV", blindMin, blindMax) : "no") << "\n"
               << "Data: " << data_int << "\n"
               << "MC: "   << mc_int
               << "  (DY=" << dy_int << ", EWK=" << ewk_int << ", TT=" << tt_int
               << ", VV=" << vv_int << ", VVV=" << vvv_int << ")\n"
               << "Data/MC: " << dataOverMc << "\n"
               << "Era lumi (fb^-1): " << (eraLumiPb/1000.0) << "\n"
               << "Effective data lumi (fb^-1): " << (effLumiPb/1000.0) << "\n"
               << "------------------------------------------\n";
  }

  // ---------- STACKED MC IN TOP PAD ----------
  // Build a stack for the MC backgrounds (order = bottom -> top in the stack)
  THStack* hs = new THStack(("hs_"+var).c_str(), "");

  // Choose stacking order to taste; here: EWK, VVV, VV, TT, DY
  hs->Add(h_EWK);
  hs->Add(h_VVV);
  hs->Add(h_VV);
  hs->Add(h_TT);
  hs->Add(h_DY);

  // Y max: consider stack, data, and signal
  const double ymax = std::max({
    h_data_draw->GetMaximum(),
    hs->GetMaximum("nostack"),
    h_signal->GetMaximum()
  }) * 1.8;

  // Draw stack first so it owns the axes
  hs->Draw("HIST");
  hs->SetMaximum(ymax);
  hs->SetMinimum(0.1);

  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetTitleOffset(1.2);
  hs->GetXaxis()->SetTitle(xTitle);

  // Then draw signal and data
  h_signal->Draw("HIST SAME");
  h_data_draw->Draw("E1 SAME");

  // Blinded window vertical lines (using pad's current y range so they always show up)
  if (blind) {
    pad1->Update();
    double y0 = gPad->GetUymin();
    double y1 = gPad->GetUymax();

    TLine* v1 = new TLine(blindMin, y0, blindMin, y1);
    TLine* v2 = new TLine(blindMax, y0, blindMax, y1);
    v1->SetLineStyle(2); v2->SetLineStyle(2);
    v1->SetLineWidth(3); v2->SetLineWidth(3);
    v1->SetLineColor(kGray+2); v2->SetLineColor(kGray+2);
    v1->Draw("SAME"); v2->Draw("SAME");
  }

  TLegend leg(0.60, 0.50, 0.88, 0.88);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.035);
  leg.AddEntry(h_signal, "VBF H#rightarrow#mu#mu", "f");
  leg.AddEntry(h_TT,     "t#bar{t} + Jets",        "f");
  leg.AddEntry(h_DY,     "DY + Jets",              "f");
  leg.AddEntry(h_EWK,    "EWK + Jets",             "f");
  leg.AddEntry(h_VV,     "VV + Jets",              "f");
  leg.AddEntry(h_VVV,    "VVV + Jets",             "f");
  leg.AddEntry(h_data,   "Data",                   "lep");
  leg.Draw();

  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.04);
  latex.DrawLatex(0.13, 0.93, "CMS #scale[0.8]{#it{Preliminary}}");
  latex.DrawLatex(0.52, 0.93, Form("#sqrt{s} = 13.6 TeV, %.2f fb^{-1}", eraLumiPb/1000.0));

  // bottom ratio
  pad2->cd();
  TH1* h_mctot = (TH1*)h_DY->Clone(("h_mctot_"+var).c_str());
  h_mctot->SetDirectory(nullptr);
  h_mctot->Add(h_EWK); h_mctot->Add(h_TT); h_mctot->Add(h_VV); h_mctot->Add(h_VVV);

  TH1* h_ratio = (TH1*)h_data->Clone(("h_ratio_"+var).c_str());
  h_ratio->SetDirectory(nullptr);
  for (int i=1;i<=h_ratio->GetNbinsX();++i){
    if (h_mctot->GetBinContent(i)<=0){ h_ratio->SetBinContent(i,0); h_ratio->SetBinError(i,0); }
  }
  h_ratio->Divide(h_mctot);
  if (blind) {
    const int bmin = h_ratio->FindBin(blindMin), bmax = h_ratio->FindBin(blindMax);
    for (int i=bmin;i<=bmax;++i){ h_ratio->SetBinContent(i,0); h_ratio->SetBinError(i,0); }
  }

  h_ratio->GetYaxis()->SetTitle("Data/MC");
  h_ratio->GetYaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetTitleOffset(0.4);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetNdivisions(505);
  h_ratio->GetXaxis()->SetLabelSize(0.12);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetXaxis()->SetTitle(xTitle);
  h_ratio->SetMinimum(0.6);
  h_ratio->SetMaximum(1.4);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(20);
  h_ratio->Draw("E1");

  // Horizontal line at Data/MC = 1 with proper x-range
  pad2->Update();
  {
    double rxmin = gPad->GetUxmin();
    double rxmax = gPad->GetUxmax();
    TLine* line = new TLine(rxmin, 1.0, rxmax, 1.0);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->Draw("SAME");
  }

  // Blinded region vertical lines in ratio pad
  if (blind) {
    pad2->Update();
    double rymin = gPad->GetUymin();
    double rymax = gPad->GetUymax();

    TLine* v1r = new TLine(blindMin, rymin, blindMin, rymax);
    TLine* v2r = new TLine(blindMax, rymin, blindMax, rymax);
    v1r->SetLineStyle(2); v2r->SetLineStyle(2);
    v1r->SetLineWidth(3); v2r->SetLineWidth(3);
    v1r->SetLineColor(kGray+2); v2r->SetLineColor(kGray+2);
    v1r->Draw("SAME"); v2r->Draw("SAME");
  }

  gSystem->Exec("mkdir -p " + outdir);
  const TString baseName = "OverlayRatio_" + TString(var.c_str());
  c->SaveAs(outdir + baseName + ".pdf");
  c->SaveAs(outdir + baseName + ".png");

  // clean local clones (file-owned hists were cloned; delete our clones)
  delete h_data; delete h_DY; delete h_EWK; delete h_TT;
  delete h_WZ2L2Q; delete h_WZ3L; delete h_ZZ2L2N; delete h_ZZ2L2Q; delete h_ZZ4L;
  delete h_WWW; delete h_WWZ; delete h_WZZ; delete h_ZZZ; delete h_signal;
  delete h_VV; delete h_VVV; delete h_data_draw; delete h_mctot; delete h_ratio;
  delete hs;
}

// ---------- driver ----------
// Run on exactly one directory (the one you pass).
// basePath should be the directory that contains DATA/, DY/, EWK/, TT/, etc.
void Muon_Overlay_AllPlots(TString basePath, TString era = "12022") {
  // Ensure trailing slash
  if (!basePath.EndsWith("/")) basePath += "/";

  // Lumi lookup (pb) – same as your original values
  double eraLumiPb = 0.0;
  if      (era == "12022") eraLumiPb =  7980.454145;
  else if (era == "22022") eraLumiPb =  26671.60971;
  else if (era == "12023") eraLumiPb =  18704.13341;
  else if (era == "22023") eraLumiPb =  9693.130053;
  else if (era == "2024")  eraLumiPb =  104805.9585;
  else if (era == "2025")  eraLumiPb =  114040.38;
  else {
    std::cerr << "[WARN] Unknown era \"" << era << "\"; setting eraLumiPb = 0\n";
  }

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

  const std::vector<std::string> vars = {
    "DiMuonMass","DiMuonMass_Z","DiMuonMass_Higgs","DiMuonpt","DiMuoneta","DiMuonrapidity",
    "Muon1pt","Muon2pt","Muon1eta","Muon2eta",
    "Muon1PTdDM","Muon2PTdDM","CollinsSoper","cosCollinsSoper","PtBalance","Zeppenfeld"
  };
  //const std::string dir = "corr_DM_VBF/";
  const std::string dir = "corr_DM_Pre-Selection/";

  // Open data
  std::unique_ptr<TFile> f_data(TFile::Open(base + "DATA/DATA_all.root"));
  if (!f_data) {
    std::cerr << "[ERROR] Missing " << (base + "DATA/DATA_all.root") << "\n";
    return;
  }

  // Open MC samples
  std::unique_ptr<TFile> f_DY       ( open_or_warn(base + "DY/DY_all.root") );
  std::unique_ptr<TFile> f_EWK      ( open_or_warn(base + "EWK/EWK_all.root") );
  std::unique_ptr<TFile> f_TT       ( open_or_warn(base + "TT/TT_all.root") );
  std::unique_ptr<TFile> f_WZ2L2Q   ( open_or_warn(base + "WZto2L2Q/WZto2L2Q_all.root") );
  std::unique_ptr<TFile> f_WZ3LNu   ( open_or_warn(base + "WZto3LNu/WZto3LNu_all.root") );
  std::unique_ptr<TFile> f_ZZ2L2Nu  ( open_or_warn(base + "ZZto2L2Nu/ZZto2L2Nu_all.root") );
  std::unique_ptr<TFile> f_ZZ2L2Q   ( open_or_warn(base + "ZZto2L2Q/ZZto2L2Q_all.root") );
  std::unique_ptr<TFile> f_ZZ4L     ( open_or_warn(base + "ZZto4L/ZZto4L_all.root") );
  std::unique_ptr<TFile> f_WWW      ( open_or_warn(base + "WWW/WWW_all.root") );
  std::unique_ptr<TFile> f_WWZ      ( open_or_warn(base + "WWZ/WWZ_all.root") );
  std::unique_ptr<TFile> f_WZZ      ( open_or_warn(base + "WZZ/WZZ_all.root") );
  std::unique_ptr<TFile> f_ZZZ      ( open_or_warn(base + "ZZZ/ZZZ_all.root") );
  std::unique_ptr<TFile> f_signal   ( open_or_warn(base + "VBF/VBF_all.root") );

  if (!f_DY || !f_EWK || !f_TT || !f_WZ2L2Q || !f_WZ3LNu || !f_ZZ2L2Nu ||
      !f_ZZ2L2Q || !f_ZZ4L || !f_WWW || !f_WWZ || !f_WZZ || !f_ZZZ || !f_signal) {
    std::cerr << "[ERROR] Missing one or more MC files under base path:\n  " << base << "\n";
    return;
  }

  // ---------------------------------------------------------
  // Save in original-like structure under fixed era directory:
  //   /eos/user/d/dblend/StackedPlots/12022/Muons/<CorrType>/<ConfigName>/
  // ---------------------------------------------------------
  TString outdir = "/eos/user/d/dblend/StackedPlots/" + era + "/Muons/" + corrType + "/" + configName + "/";
  gSystem->Exec("mkdir -p " + outdir);
  std::ofstream summary((outdir + "OverlayRatios_All.txt").Data());

  for (const auto& v : vars) {
    const bool doBlind = (v == "DiMuonMass" || v == "DiMuonMass_Higgs");
    const char* xlab = (v=="DiMuonMass") ? "m_{#mu#mu} [GeV]" :
                       (v=="DiMuonMass_Z") ? "m_{#mu#mu} [GeV]" :
                       (v=="DiMuonMass_Higgs") ? "m_{#mu#mu} [GeV]" :
                       (v=="DiMuonpt") ? "p_{T}(#mu#mu) [GeV]" :
                       (v=="DiMuoneta") ? "#eta(#mu#mu)" :
                       (v=="DiMuonrapidity") ? "y(#mu#mu)" :
                       (v=="Muon1pt") ? "p_{T}(#mu_{1}) [GeV]" :
                       (v=="Muon2pt") ? "p_{T}(#mu_{2}) [GeV]" :
                       (v=="Muon1eta") ? "#eta(#mu_{1})" :
                       (v=="Muon2eta") ? "#eta(#mu_{2})" :
                       (v=="Muon1PTdDM") ? "p_{T}(#mu_{1}) / m_{#mu#mu}" :
                       (v=="Muon2PTdDM") ? "p_{T}(#mu_{2}) / m_{#mu#mu}" :
                       (v=="CollinsSoper") ? "#phi_{CS}" :
                       (v=="cosCollinsSoper") ? "cos#theta_{CS}" :
                       (v=="PtBalance") ? "p_{T}^{balance}" :
                       (v=="Zeppenfeld") ? "z^{*}" : v.c_str();

    draw_one(outdir, dir, v,
             f_data.get(), f_DY.get(), f_EWK.get(), f_TT.get(),
             f_WZ2L2Q.get(), f_WZ3LNu.get(), f_ZZ2L2Nu.get(), f_ZZ2L2Q.get(), f_ZZ4L.get(),
             f_WWW.get(), f_WWZ.get(), f_WZZ.get(), f_ZZZ.get(),
             f_signal.get(), xlab, doBlind, eraLumiPb, &summary);
  }
}
// ============================================================================

