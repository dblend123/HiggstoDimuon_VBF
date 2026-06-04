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
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <cmath>

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

static TH1* getTotalVariationBand(
    const std::string& name,
    TH1* h_nom,
    TH1* h_up,
    TH1* h_down,
    Color_t color,
    int fillStyle = 1001
) {
    TH1* band = (TH1*)h_nom->Clone(name.c_str());
    band->SetDirectory(nullptr);

    for (int i = 1; i <= band->GetNbinsX(); ++i) {
        const double nom  = h_nom->GetBinContent(i);
        const double up   = h_up   ? h_up->GetBinContent(i)   : nom;
        const double down = h_down ? h_down->GetBinContent(i) : nom;

        const double errUp   = std::abs(up - nom);
        const double errDown = std::abs(nom - down);
        const double err = std::max(errUp, errDown);

        band->SetBinContent(i, 1.0);
        band->SetBinError(i, nom > 0 ? err / nom : 0.0);
    }

    band->SetFillColor(color);
    band->SetFillStyle(fillStyle);
    band->SetMarkerSize(0);
    band->SetLineColor(color);
    return band;
}

// Map nominal overlay variable names to the BDTTheorySyst_VBF histogram names.
// Return "" for variables that do not currently have theory-variation clones.
static std::string theoryVarName(const std::string& var) {
  if (var == "DiMuonMass" || var == "DiMuonMass_Higgs") return "DimuonMass";
  if (var == "DiMuonpt")        return "DimuonPt";
  if (var == "DiMuonrapidity")  return "DimuonRapidity";
  if (var == "Muon1PTdDM")      return "Muon1PtOverM";
  if (var == "Muon2PTdDM")      return "Muon2PtOverM";
  if (var == "Muon1eta")        return "Muon1Eta";
  if (var == "Muon2eta")        return "Muon2Eta";
  if (var == "CollinsSoper")    return "CollinsSoperPhi";
  if (var == "cosCollinsSoper") return "CollinsSoperCosTheta";
  if (var == "PtBalance")       return "PtBalance";
  if (var == "Zeppenfeld")      return "Zeppenfeld";
  return "";
}

static void addHistOrNominal(TH1* total, TFile* f, const std::string& path, TH1* h_nominal) {
  if (!total || !h_nominal) return;

  TH1* h = getHist(f, path);
  if (h) {
    total->Add(h);
    delete h;
  } else {
    // If a sample does not have a variation histogram, keep it at nominal.
    // This lets the band still be built from the samples that do have variations.
    total->Add(h_nominal);
  }
}

static TH1* buildTotalVariation(
    const std::string& name,
    const std::string& systDir,
    const std::string& systVar,
    const std::string& suffix,
    TH1* h_template,
    TFile* f_DY,      TH1* h_DY,
    TFile* f_EWK,     TH1* h_EWK,
    TFile* f_TT,      TH1* h_TT,
    TFile* f_WZ2L2Q,  TH1* h_WZ2L2Q,
    TFile* f_WZ3LNu,  TH1* h_WZ3L,
    TFile* f_ZZ2L2Nu, TH1* h_ZZ2L2N,
    TFile* f_ZZ2L2Q,  TH1* h_ZZ2L2Q,
    TFile* f_ZZ4L,    TH1* h_ZZ4L,
    TFile* f_WWW,     TH1* h_WWW,
    TFile* f_WWZ,     TH1* h_WWZ,
    TFile* f_WZZ,     TH1* h_WZZ,
    TFile* f_ZZZ,     TH1* h_ZZZ
) {
  if (systVar.empty() || !h_template) return nullptr;

  TH1* total = (TH1*)h_template->Clone(name.c_str());
  total->SetDirectory(nullptr);
  total->Reset("ICES");

  const std::string path = systDir + systVar + "_" + suffix;

  addHistOrNominal(total, f_DY,      path, h_DY);
  addHistOrNominal(total, f_EWK,     path, h_EWK);
  addHistOrNominal(total, f_TT,      path, h_TT);
  addHistOrNominal(total, f_WZ2L2Q,  path, h_WZ2L2Q);
  addHistOrNominal(total, f_WZ3LNu,  path, h_WZ3L);
  addHistOrNominal(total, f_ZZ2L2Nu, path, h_ZZ2L2N);
  addHistOrNominal(total, f_ZZ2L2Q,  path, h_ZZ2L2Q);
  addHistOrNominal(total, f_ZZ4L,    path, h_ZZ4L);
  addHistOrNominal(total, f_WWW,     path, h_WWW);
  addHistOrNominal(total, f_WWZ,     path, h_WWZ);
  addHistOrNominal(total, f_WZZ,     path, h_WZZ);
  addHistOrNominal(total, f_ZZZ,     path, h_ZZZ);

  return total;
}

static TH1* buildOrangeTheoryBand(
    const std::string& name,
    TH1* h_nom,
    TH1* h_mcStatBand,
    TH1* h_scaleUp, TH1* h_scaleDown,
    TH1* h_pdfUp,   TH1* h_pdfDown,
    TH1* h_isrUp,   TH1* h_isrDown,
    TH1* h_fsrUp,   TH1* h_fsrDown
) {
  if (!h_nom) return nullptr;

  TH1* band = (TH1*)h_nom->Clone(name.c_str());
  band->SetDirectory(nullptr);

  for (int i = 1; i <= band->GetNbinsX(); ++i) {
    const double nom = h_nom->GetBinContent(i);

    auto oneSyst = [&](TH1* up, TH1* down) {
      if (!up && !down) return 0.0;
      const double u = up   ? std::abs(up->GetBinContent(i)   - nom) : 0.0;
      const double d = down ? std::abs(nom - down->GetBinContent(i)) : 0.0;
      return std::max(u, d);
    };

    const double eScale = oneSyst(h_scaleUp, h_scaleDown);
    const double ePdf   = oneSyst(h_pdfUp,   h_pdfDown);
    const double eIsr   = oneSyst(h_isrUp,   h_isrDown);
    const double eFsr   = oneSyst(h_fsrUp,   h_fsrDown);

    const double statRel = h_mcStatBand ? h_mcStatBand->GetBinError(i) : 0.0;
    const double statAbs = nom * statRel;

    const double totalAbs = std::sqrt(statAbs*statAbs +
                                      eScale*eScale +
                                      ePdf*ePdf +
                                      eIsr*eIsr +
                                      eFsr*eFsr);

    band->SetBinContent(i, 1.0);
    band->SetBinError(i, nom > 0 ? totalAbs / nom : 0.0);
  }

  band->SetFillColor(kOrange-3);
  band->SetFillStyle(1001);
  band->SetMarkerSize(0);
  band->SetLineColor(kOrange-3);
  return band;
}



static TH1* buildAzureExperimentalBand(
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

  for (int i = 1; i <= band->GetNbinsX(); ++i) {
    const double nom = h_nom->GetBinContent(i);

    auto oneSyst = [&](TH1* up, TH1* down) {
      if (!up && !down) return 0.0;
      const double u = up   ? std::abs(up->GetBinContent(i)   - nom) : 0.0;
      const double d = down ? std::abs(nom - down->GetBinContent(i)) : 0.0;
      return std::max(u, d);
    };

    const double ePU   = oneSyst(h_puUp,   h_puDown);
    const double eMuon = oneSyst(h_muonUp, h_muonDown);
    const double eTrig = oneSyst(h_trigUp, h_trigDown);
    const double eBTag = oneSyst(h_btagUp, h_btagDown);
    const double eZPt  = oneSyst(h_zptUp,  h_zptDown);

    const double statRel = h_mcStatBand ? h_mcStatBand->GetBinError(i) : 0.0;
    const double statAbs = nom * statRel;

    const double totalAbs = std::sqrt(statAbs*statAbs +
                                      ePU*ePU +
                                      eMuon*eMuon +
                                      eTrig*eTrig +
                                      eBTag*eBTag +
                                      eZPt*eZPt);

    band->SetBinContent(i, 1.0);
    band->SetBinError(i, nom > 0 ? totalAbs / nom : 0.0);
  }

  band->SetFillColor(kAzure+1);
  band->SetFillStyle(1001);
  band->SetMarkerSize(0);
  band->SetLineColor(kAzure+1);
  return band;
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
  const std::string systDir = "BDTTheorySyst_VBF/";

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

// Empty frame owns the axes
TH1* h_frame = (TH1*)h_data_draw->Clone(("h_frame_"+var).c_str());
h_frame->Reset("ICES");
h_frame->SetTitle("");
h_frame->SetMinimum(0.1);
h_frame->SetMaximum(1.0e8);

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
h_data_draw->Draw("E1 SAME");

gPad->RedrawAxis();

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

  TLatex latex; 
  latex.SetNDC(); 
  latex.SetTextFont(62);
  latex.DrawLatex(0.12, 0.93, "CMS");

  latex.SetTextFont(52); // italic
  latex.DrawLatex(0.19, 0.93, "Preliminary");

  latex.SetTextFont(62); // bold
  latex.DrawLatex(0.62, 0.93, Form("%.2f fb^{-1} (13.6 TeV)", eraLumiPb/1000.0));

  // bottom ratio
  pad2->cd();
  TH1* h_mctot = (TH1*)h_DY->Clone(("h_mctot_"+var).c_str());
  h_mctot->SetDirectory(nullptr);
  h_mctot->Add(h_EWK); 
  h_mctot->Add(h_TT); 
  h_mctot->Add(h_VV); 
  h_mctot->Add(h_VVV);

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

// ---------- Orange theory band: MC stat ⊕ scale/PDF/ISR/FSR ----------
const std::string systVar = theoryVarName(var);
if (systVar.empty()) {
  std::cerr << "[THEORY BAND] No BDTTheorySyst_VBF mapping for " << var
            << "; drawing gray stat band only.\n";
}

TH1* h_scaleUpTot = buildTotalVariation("h_scaleUpTot_"+var, systDir, systVar, "scaleUp", h_mctot,
                                        f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                        f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                        f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                        f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_scaleDownTot = buildTotalVariation("h_scaleDownTot_"+var, systDir, systVar, "scaleDown", h_mctot,
                                          f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                          f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                          f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                          f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_pdfUpTot = buildTotalVariation("h_pdfUpTot_"+var, systDir, systVar, "pdfUp", h_mctot,
                                      f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                      f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                      f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                      f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_pdfDownTot = buildTotalVariation("h_pdfDownTot_"+var, systDir, systVar, "pdfDown", h_mctot,
                                        f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                        f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                        f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                        f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_isrUpTot = buildTotalVariation("h_isrUpTot_"+var, systDir, systVar, "isrUp", h_mctot,
                                      f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                      f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                      f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                      f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_isrDownTot = buildTotalVariation("h_isrDownTot_"+var, systDir, systVar, "isrDown", h_mctot,
                                        f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                        f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                        f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                        f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_fsrUpTot = buildTotalVariation("h_fsrUpTot_"+var, systDir, systVar, "fsrUp", h_mctot,
                                      f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                      f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                      f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                      f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_fsrDownTot = buildTotalVariation("h_fsrDownTot_"+var, systDir, systVar, "fsrDown", h_mctot,
                                        f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                        f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                        f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                        f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);

TH1* h_theoryBand = buildOrangeTheoryBand("h_theoryBand_"+var, h_mctot, h_mcBand,
                                          h_scaleUpTot, h_scaleDownTot,
                                          h_pdfUpTot,   h_pdfDownTot,
                                          h_isrUpTot,   h_isrDownTot,
                                          h_fsrUpTot,   h_fsrDownTot);


// ---------- Azure experimental band: MC stat ⊕ PU/muon/trigger/btag/ZpT ----------
const std::string expSystDir = "BDTExperimentalSyst_VBF/";

TH1* h_puUpTot = buildTotalVariation("h_puUpTot_"+var, expSystDir, systVar, "puUp", h_mctot,
                                     f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                     f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                     f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                     f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_puDownTot = buildTotalVariation("h_puDownTot_"+var, expSystDir, systVar, "puDown", h_mctot,
                                       f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                       f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                       f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                       f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_muonSFUpTot = buildTotalVariation("h_muonSFUpTot_"+var, expSystDir, systVar, "muonSFUp", h_mctot,
                                         f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                         f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                         f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                         f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_muonSFDownTot = buildTotalVariation("h_muonSFDownTot_"+var, expSystDir, systVar, "muonSFDown", h_mctot,
                                           f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                           f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                           f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                           f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_trigUpTot = buildTotalVariation("h_trigUpTot_"+var, expSystDir, systVar, "trigUp", h_mctot,
                                       f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                       f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                       f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                       f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_trigDownTot = buildTotalVariation("h_trigDownTot_"+var, expSystDir, systVar, "trigDown", h_mctot,
                                         f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                         f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                         f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                         f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_btagUpTot = buildTotalVariation("h_btagUpTot_"+var, expSystDir, systVar, "btagUp", h_mctot,
                                       f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                       f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                       f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                       f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_btagDownTot = buildTotalVariation("h_btagDownTot_"+var, expSystDir, systVar, "btagDown", h_mctot,
                                         f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                         f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                         f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                         f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_zptUpTot = buildTotalVariation("h_zptUpTot_"+var, expSystDir, systVar, "zptUp", h_mctot,
                                      f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                      f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                      f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                      f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);
TH1* h_zptDownTot = buildTotalVariation("h_zptDownTot_"+var, expSystDir, systVar, "zptDown", h_mctot,
                                        f_DY,h_DY, f_EWK,h_EWK, f_TT,h_TT,
                                        f_WZ2L2Q,h_WZ2L2Q, f_WZ3LNu,h_WZ3L,
                                        f_ZZ2L2Nu,h_ZZ2L2N, f_ZZ2L2Q,h_ZZ2L2Q, f_ZZ4L,h_ZZ4L,
                                        f_WWW,h_WWW, f_WWZ,h_WWZ, f_WZZ,h_WZZ, f_ZZZ,h_ZZZ);

TH1* h_expBand = buildAzureExperimentalBand("h_expBand_"+var, h_mctot, h_mcBand,
                                            h_puUpTot,     h_puDownTot,
                                            h_muonSFUpTot, h_muonSFDownTot,
                                            h_trigUpTot,   h_trigDownTot,
                                            h_btagUpTot,   h_btagDownTot,
                                            h_zptUpTot,    h_zptDownTot);


  TH1* h_ratio = (TH1*)h_data->Clone(("h_ratio_"+var).c_str());
  h_ratio->SetDirectory(nullptr);
  h_ratio->SetTitle("");
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
  h_ratio->SetMinimum(0.5);
  h_ratio->SetMaximum(1.5);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(20);
  h_ratio->SetStats(0);
h_ratio->Draw("AXIS");

// Draw widest bands first, then gray MC-stat on top.
if (h_expBand)    h_expBand->Draw("E2 SAME");
if (h_theoryBand) h_theoryBand->Draw("E2 SAME");
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
  delete h_scaleUpTot; delete h_scaleDownTot; delete h_pdfUpTot; delete h_pdfDownTot;
  delete h_isrUpTot; delete h_isrDownTot; delete h_fsrUpTot; delete h_fsrDownTot;
  delete h_theoryBand;
  delete h_expBand;
  delete h_puUpTot; delete h_puDownTot;
  delete h_muonSFUpTot; delete h_muonSFDownTot;
  delete h_trigUpTot; delete h_trigDownTot;
  delete h_btagUpTot; delete h_btagDownTot;
  delete h_zptUpTot; delete h_zptDownTot;
  delete hs; delete h_frame; delete h_mcBand;
}

// ---------- driver ----------
// Run on exactly one directory (the one you pass).
// basePath should be the directory that contains DATA/, DY/, EWK/, TT/, etc.
void Muon_Overlay_AllPlots(TString basePath, TString era = "12022") {
  gStyle->SetOptStat(0);

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
  const std::string dir = "corr_DM_VBF/";
  //const std::string dir = "corr_DM_Pre-Selection/";

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

