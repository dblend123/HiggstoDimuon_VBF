// ================== Muon_Overlay_DiMuonpt_Full.C ==================
// Full-featured version (matches your previous Muon_Overlay_AllPlots style):
// - Uses ALL backgrounds and the same grouping (VV, VVV, etc.)
// - Includes bottom pad with Data/MC ratio
// - Same blinding logic exists (not used here since var=DiMuonpt)
// - Runs ONLY DiMuonpt
// - Loops over directories:
//    corr_DM_Pre-Selection_Z0J
//    corr_DM_Pre-Selection_Z1J
//    corr_DM_Pre-Selection_Z2J
//    corr_DM_Pre-Selection_wZpt_Z0J
//    corr_DM_Pre-Selection_wZpt_Z1J
//    corr_DM_Pre-Selection_wZpt_Z2J
//
// Usage:
//   root -l -q 'Muon_Overlay_DiMuonpt_Full.C("/eos/.../Analyzed_SkimmedFiles/12022/KIT/<config>/","12022")'

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
static void draw_one(const TString& outdir,
                     const std::string& dir,
                     const std::string& var,
                     const std::string& outTag,
                     TFile* f_data, TFile* f_DY, TFile* f_EWK, TFile* f_TT,
                     TFile* f_WZ2L2Q, TFile* f_WZ3LNu, TFile* f_ZZ2L2Nu, TFile* f_ZZ2L2Q, TFile* f_ZZ4L,
                     TFile* f_WWW, TFile* f_WWZ, TFile* f_WZZ, TFile* f_ZZZ,
                     TFile* f_signal,
                     const char* xTitle,
                     bool doBlindMass,
                     double eraLumiPb,
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
    std::cerr << "[WARN] Missing hist(s) for " << hpath << " — skipping.\n";
    delete h_data; delete h_DY; delete h_EWK; delete h_TT; delete h_WZ2L2Q; delete h_WZ3L;
    delete h_ZZ2L2N; delete h_ZZ2L2Q; delete h_ZZ4L; delete h_WWW; delete h_WWZ; delete h_WZZ; delete h_ZZZ; delete h_signal;
    return;
  }

  // groups
  TH1* h_VV  = (TH1*)h_WZ2L2Q->Clone(("h_VV_"+outTag).c_str());
  h_VV->SetDirectory(nullptr);
  h_VV->Add(h_WZ3L); h_VV->Add(h_ZZ2L2Q); h_VV->Add(h_ZZ2L2N); h_VV->Add(h_ZZ4L);

  TH1* h_VVV = (TH1*)h_WWW->Clone(("h_VVV_"+outTag).c_str());
  h_VVV->SetDirectory(nullptr);
  h_VVV->Add(h_WWZ); h_VVV->Add(h_WZZ); h_VVV->Add(h_ZZZ);

  for (TH1* h : {h_data,h_DY,h_EWK,h_TT,h_VV,h_VVV,h_signal})
    if (h && h->GetSumw2N()==0) h->Sumw2(kTRUE);

  // style (same as your original)
  h_data->SetMarkerStyle(20); h_data->SetMarkerColor(kBlack); h_data->SetLineColor(kBlack); h_data->SetMarkerSize(1.2);
  h_DY->SetFillColor(kOrange-2); h_DY->SetLineColor(kOrange+3);
  h_TT->SetFillColor(kRed-4);    h_TT->SetLineColor(kRed+1);
  h_VV->SetFillColor(kAzure-9);  h_VV->SetLineColor(kAzure+2);
  h_VVV->SetFillColor(kGreen-7); h_VVV->SetLineColor(kGreen+3);
  h_EWK->SetFillColor(kGray+1);  h_EWK->SetLineColor(kGray+3);
  h_signal->SetFillColor(kMagenta+1); h_signal->SetLineColor(kMagenta+3); h_signal->SetLineWidth(3);

  // canvas/pads
  std::unique_ptr<TCanvas> c(new TCanvas(("c_"+outTag).c_str(), outTag.c_str(), 800, 800));
  std::unique_ptr<TPad> pad1(new TPad(("pad1_"+outTag).c_str(), "Top",    0, 0.32, 1, 1));
  std::unique_ptr<TPad> pad2(new TPad(("pad2_"+outTag).c_str(), "Bottom", 0, 0.02, 1, 0.32));
  pad1->SetBottomMargin(0.03); pad2->SetTopMargin(0.03); pad2->SetBottomMargin(0.3);
  pad1->Draw(); pad2->Draw();

  // top
  pad1->cd();
  pad1->SetLogy();

  TH1* h_data_draw = (TH1*)h_data->Clone(("h_data_draw_"+outTag).c_str());
  h_data_draw->SetDirectory(nullptr);
  if (doBlindMass && var=="DiMuonMass") {
    const int bmin = h_data_draw->FindBin(110.0);
    const int bmax = h_data_draw->FindBin(150.0);
    for (int i=bmin;i<=bmax;++i){ h_data_draw->SetBinContent(i,0); h_data_draw->SetBinError(i,0); }
  }

  // integrals (exclude blinded window for both Data/MC if DiMuonMass)
  const bool blind = (doBlindMass && var=="DiMuonMass");
  auto iNoUF = [](TH1* h){ return h ? h->Integral(1,h->GetNbinsX()) : 0.0; };
  auto iMC   = [&](TH1* h){ return blind ? integralExcluding(h,110.0,150.0) : iNoUF(h); };

  const double data_int = blind ? integralExcluding(h_data,110.0,150.0) : iNoUF(h_data);
  const double dy_int   = iMC(h_DY);
  const double ewk_int  = iMC(h_EWK);
  const double tt_int   = iMC(h_TT);
  const double vv_int   = iMC(h_VV);
  const double vvv_int  = iMC(h_VVV);
  const double mc_int   = dy_int + ewk_int + tt_int + vv_int + vvv_int;

  const double dataOverMc = (mc_int>0.0 ? data_int/mc_int : 0.0);
  const double effLumiPb  = dataOverMc * eraLumiPb;

  std::cout << std::fixed
            << "[PLOT] dir=" << dir << " var=" << var
            << " Data=" << data_int << " DY=" << dy_int << " EWK=" << ewk_int
            << " TT=" << tt_int << " VV=" << vv_int << " VVV=" << vvv_int
            << " MC=" << mc_int << " Data/MC=" << dataOverMc << "\n";

  if (summary && *summary) {
    (*summary) << std::fixed << std::setprecision(6)
               << "Dir: " << dir << "\n"
               << "Variable: " << var << "\n"
               << "Blinded: " << (blind ? "yes [110,150] GeV" : "no") << "\n"
               << "Data: " << data_int << "\n"
               << "MC: "   << mc_int
               << "  (DY=" << dy_int << ", EWK=" << ewk_int << ", TT=" << tt_int
               << ", VV=" << vv_int << ", VVV=" << vvv_int << ")\n"
               << "Data/MC: " << dataOverMc << "\n"
               << "Era lumi (fb^-1): " << (eraLumiPb/1000.0) << "\n"
               << "Effective data lumi (fb^-1): " << (effLumiPb/1000.0) << "\n"
               << "------------------------------------------\n";
  }

  // Build total MC (used for ymax and ratio)
  TH1* h_mctot = (TH1*)h_DY->Clone(("h_mctot_"+outTag).c_str());
  h_mctot->SetDirectory(nullptr);
  h_mctot->Add(h_EWK);
  h_mctot->Add(h_TT);
  h_mctot->Add(h_VV);
  h_mctot->Add(h_VVV);

  // ---------- STACKED MC IN TOP PAD ----------
  THStack* hs = new THStack(("hs_"+outTag).c_str(), "");
  // Order: EWK, VVV, VV, TT, DY  (same as your original)
  hs->Add(h_EWK);
  hs->Add(h_VVV);
  hs->Add(h_VV);
  hs->Add(h_TT);
  hs->Add(h_DY);

  // Draw stack first so it owns axes
  hs->Draw("HIST");

  // Y max: use total MC (not hs->GetMaximum("nostack")) so nothing gets clipped
  const double ymax = std::max({h_data_draw->GetMaximum(), h_mctot->GetMaximum(), h_signal->GetMaximum()}) * 1.8;
  hs->SetMaximum(ymax);
  hs->SetMinimum(0.1);

  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetLabelSize(0.04);
  hs->GetYaxis()->SetTitleOffset(1.2);
  hs->GetXaxis()->SetTitle(xTitle);

  // draw signal + data
  h_signal->Draw("HIST SAME");
  h_data_draw->Draw("E1 SAME");

  // blinded window lines (kept for compatibility)
  if (blind) {
    pad1->Update();
    double y0 = gPad->GetUymin();
    double y1 = gPad->GetUymax();
    TLine* v1 = new TLine(110.0, y0, 110.0, y1);
    TLine* v2 = new TLine(150.0, y0, 150.0, y1);
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

  TH1* h_ratio = (TH1*)h_data->Clone(("h_ratio_"+outTag).c_str());
  h_ratio->SetDirectory(nullptr);

  for (int i=1;i<=h_ratio->GetNbinsX();++i){
    if (h_mctot->GetBinContent(i)<=0){ h_ratio->SetBinContent(i,0); h_ratio->SetBinError(i,0); }
  }
  h_ratio->Divide(h_mctot);

  if (blind) {
    const int bmin = h_ratio->FindBin(110.0), bmax = h_ratio->FindBin(150.0);
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
  h_ratio->SetMinimum(0.0);
  h_ratio->SetMaximum(2.0);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(20);
  h_ratio->Draw("E1");

  // horizontal line at 1
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

  // blinded region lines in ratio pad
  if (blind) {
    pad2->Update();
    double rymin = gPad->GetUymin();
    double rymax = gPad->GetUymax();
    TLine* v1r = new TLine(110.0, rymin, 110.0, rymax);
    TLine* v2r = new TLine(150.0, rymin, 150.0, rymax);
    v1r->SetLineStyle(2); v2r->SetLineStyle(2);
    v1r->SetLineWidth(3); v2r->SetLineWidth(3);
    v1r->SetLineColor(kGray+2); v2r->SetLineColor(kGray+2);
    v1r->Draw("SAME"); v2r->Draw("SAME");
  }

  gSystem->Exec("mkdir -p " + outdir);

  const TString baseName = "OverlayRatio_" + TString(outTag.c_str());
  c->SaveAs(outdir + baseName + ".pdf");
  c->SaveAs(outdir + baseName + ".png");

  // clean
  delete h_data; delete h_DY; delete h_EWK; delete h_TT;
  delete h_WZ2L2Q; delete h_WZ3L; delete h_ZZ2L2N; delete h_ZZ2L2Q; delete h_ZZ4L;
  delete h_WWW; delete h_WWZ; delete h_WZZ; delete h_ZZZ; delete h_signal;
  delete h_VV; delete h_VVV; delete h_data_draw; delete h_mctot; delete h_ratio;
  delete hs;
}

// ---------- driver ----------
void Muon_Overlay_DiMuonpt(TString basePath, TString era = "12022") {

  if (!basePath.EndsWith("/")) basePath += "/";

  double eraLumiPb = 0.0;
  if      (era == "12022") eraLumiPb =  7980.454145;
  else if (era == "22022") eraLumiPb =  26671.60971;
  else if (era == "12023") eraLumiPb =  18704.13341;
  else if (era == "22023") eraLumiPb =  9693.130053;
  else if (era == "2024")  eraLumiPb =  104805.9585;
  else std::cerr << "[WARN] Unknown era \"" << era << "\"; eraLumiPb = 0\n";

  // Parse CorrType and ConfigName from basePath: .../<era>/<CorrType>/<ConfigName>/
  TString tmp = basePath;
  if (tmp.EndsWith("/")) tmp.Chop();
  Int_t last1 = tmp.Last('/');
  TString configName = "UnknownConfig";
  TString upOne      = tmp;
  if (last1 >= 0) {
    configName = tmp(last1+1, tmp.Length() - last1 - 1);
    upOne      = tmp(0, last1);
  }
  Int_t last2 = upOne.Last('/');
  TString corrType = "UnknownCorr";
  if (last2 >= 0) corrType = upOne(last2+1, upOne.Length() - last2 - 1);

  // only variable requested
  const std::string var = "DiMuonpt";
  const char* xlab = "p_{T}(#mu#mu) [GeV]";

  // directories requested
  const std::vector<std::string> dirs = {
    "corr_DM_Pre-Selection_Z0j/",
    "corr_DM_Pre-Selection_Z1j/",
    "corr_DM_Pre-Selection_Z2j/",
    "corr_DM_Pre-Selection_wZpt_Z0j/",
    "corr_DM_Pre-Selection_wZpt_Z1j/",
    "corr_DM_Pre-Selection_wZpt_Z2j/"
  };

  // open files
  std::unique_ptr<TFile> f_data   ( open_or_warn(basePath + "DATA/DATA_all.root") );
  if (!f_data) { std::cerr << "[ERROR] Missing DATA file\n"; return; }

  std::unique_ptr<TFile> f_DY       ( open_or_warn(basePath + "DY/DY_all.root") );
  std::unique_ptr<TFile> f_EWK      ( open_or_warn(basePath + "EWK/EWK_all.root") );
  std::unique_ptr<TFile> f_TT       ( open_or_warn(basePath + "TT/TT_all.root") );
  std::unique_ptr<TFile> f_WZ2L2Q   ( open_or_warn(basePath + "WZto2L2Q/WZto2L2Q_all.root") );
  std::unique_ptr<TFile> f_WZ3LNu   ( open_or_warn(basePath + "WZto3LNu/WZto3LNu_all.root") );
  std::unique_ptr<TFile> f_ZZ2L2Nu  ( open_or_warn(basePath + "ZZto2L2Nu/ZZto2L2Nu_all.root") );
  std::unique_ptr<TFile> f_ZZ2L2Q   ( open_or_warn(basePath + "ZZto2L2Q/ZZto2L2Q_all.root") );
  std::unique_ptr<TFile> f_ZZ4L     ( open_or_warn(basePath + "ZZto4L/ZZto4L_all.root") );
  std::unique_ptr<TFile> f_WWW      ( open_or_warn(basePath + "WWW/WWW_all.root") );
  std::unique_ptr<TFile> f_WWZ      ( open_or_warn(basePath + "WWZ/WWZ_all.root") );
  std::unique_ptr<TFile> f_WZZ      ( open_or_warn(basePath + "WZZ/WZZ_all.root") );
  std::unique_ptr<TFile> f_ZZZ      ( open_or_warn(basePath + "ZZZ/ZZZ_all.root") );
  std::unique_ptr<TFile> f_signal   ( open_or_warn(basePath + "VBF/VBF_all.root") );

  if (!f_DY || !f_EWK || !f_TT || !f_WZ2L2Q || !f_WZ3LNu || !f_ZZ2L2Nu ||
      !f_ZZ2L2Q || !f_ZZ4L || !f_WWW || !f_WWZ || !f_WZZ || !f_ZZZ || !f_signal) {
    std::cerr << "[ERROR] Missing one or more MC files under base path:\n  " << basePath << "\n";
    return;
  }

  // output directory (same structure as original)
  TString outdir = "/eos/user/d/dblend/StackedPlots/" + era + "/ZPt_Corrections/" + corrType + "/" + configName + "/";
  gSystem->Exec("mkdir -p " + outdir);
  std::ofstream summary((outdir + "OverlayRatios_DiMuonpt_AllDirs.txt").Data());

  // loop directories
  for (const auto& d : dirs) {
    // build a clean tag for filename
    std::string tag = d;
    if (!tag.empty() && tag.back()=='/') tag.pop_back();
    // unique output tag: <dir>__DiMuonpt
    std::string outTag = tag + "__" + var;

    draw_one(outdir, d, var, outTag,
             f_data.get(), f_DY.get(), f_EWK.get(), f_TT.get(),
             f_WZ2L2Q.get(), f_WZ3LNu.get(), f_ZZ2L2Nu.get(), f_ZZ2L2Q.get(), f_ZZ4L.get(),
             f_WWW.get(), f_WWZ.get(), f_WZZ.get(), f_ZZZ.get(),
             f_signal.get(),
             xlab,
             false, // no blinding for DiMuonpt
             eraLumiPb,
             &summary);
  }

  std::cout << "[DONE] Wrote plots to: " << outdir << "\n";
}
