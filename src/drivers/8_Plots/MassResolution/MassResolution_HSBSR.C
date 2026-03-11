#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TLine.h"
#include "TBox.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TString.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooBreitWigner.h"
#include "RooBernstein.h"
#include "RooVoigtian.h"
#include "RooCrystalBall.h"
#include "RooGenericPdf.h"
#include "RooArgList.h"
#include "RooConstVar.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;

namespace {

static bool isDirectory(const char* p) {
  if (!p) return false;
  void* dir = gSystem->OpenDirectory(p);
  if (dir) { gSystem->FreeDirectory(dir); return true; }
  return false;
}


static TString uniqueName(const char* base) {
  static unsigned long counter = 0;
  TString b = base ? base : "obj";
  b.ReplaceAll("/", "_");
  b.ReplaceAll(" ", "_");
  return Form("%s_%lu", b.Data(), ++counter);
}

static std::string sanitize(std::string s) {
  for (char& c : s) {
    if (c == '/' || c == ' ' || c == '(' || c == ')' || c == ':' || c == ';') c = '_';
  }
  return s;
}

static double eraLumiPb(const std::string& era) {
  if      (era == "12022") return  7980.454145;
  else if (era == "22022") return 26671.60971;
  else if (era == "12023") return 18704.13341;
  else if (era == "22023") return  9693.130053;
  else if (era == "2024")  return 104805.9585;
  std::cerr << "[WARN] Unknown era '" << era << "'; setting lumi=0\n";
  return 0.0;
}

static TH1D* makeWindowedHist(const TH1* hIn, double xmin, double xmax, const char* newname) {
  if (!hIn) return nullptr;
  const int b1 = std::max(1, hIn->GetXaxis()->FindBin(xmin));
  const int b2 = std::min(hIn->GetNbinsX(), hIn->GetXaxis()->FindBin(xmax - 1e-6));
  const int nb = b2 - b1 + 1;
  if (nb <= 0) return nullptr;

  std::vector<double> edges; edges.reserve(nb + 1);
  for (int b = b1; b <= b2 + 1; ++b) edges.push_back(hIn->GetXaxis()->GetBinLowEdge(b));

  TString uname = uniqueName(newname);
  TH1D* h = new TH1D(uname, hIn->GetTitle(), nb, edges.data());
  h->SetDirectory(nullptr);
  h->Sumw2();
  for (int i = 0; i < nb; ++i) {
    const int b = b1 + i;
    h->SetBinContent(i+1, hIn->GetBinContent(b));
    h->SetBinError(i+1, hIn->GetBinError(b));
  }
  h->GetXaxis()->SetTitle(hIn->GetXaxis()->GetTitle());
  h->GetYaxis()->SetTitle(hIn->GetYaxis()->GetTitle());
  return h;
}

static TH1D* makeBlindedDataHist(const TH1* hIn, double blindLo, double blindHi, const char* newname) {
  if (!hIn) return nullptr;
  TString uname = uniqueName(newname);
  TH1D* h = (TH1D*)hIn->Clone(uname);
  h->SetDirectory(nullptr);
  for (int b = 1; b <= h->GetNbinsX(); ++b) {
    const double x = h->GetBinCenter(b);
    if (x >= blindLo && x <= blindHi) {
      h->SetBinContent(b, 0.0);
      h->SetBinError(b, 0.0);
    }
  }
  return h;
}

static TH1D* makeSidebandDataHist(const TH1* hIn, double blindLo, double blindHi, const char* newname) {
  return makeBlindedDataHist(hIn, blindLo, blindHi, newname);
}

static TString guessCorrectionsTag(const std::string& inputPath) {
  TString s = inputPath.c_str();
  if (s.EndsWith("/")) s.Chop();
  if (!isDirectory(s.Data())) s = s(0, s.Last('/'));
  return s(s.Last('/') + 1, s.Length());
}

static TString guessCorrTypeTag(const std::string& inputPath) {
  TString s = inputPath.c_str();
  if (s.EndsWith("/")) s.Chop();
  if (!isDirectory(s.Data())) s = s(0, s.Last('/'));
  TString corrections = s(s.Last('/') + 1, s.Length());
  s = s(0, s.Last('/'));
  return s(s.Last('/') + 1, s.Length());
}

static bool isDataFile(const std::string& p) {
  return (p.find("/DATA/") != std::string::npos) || (p.find("DATA_all.root") != std::string::npos);
}

static bool isSignalFile(const std::string& p) {
  static const std::vector<std::string> tags = {
    "/VBF/", "/GGH/", "/GluGluH/", "/VH/", "/ttH/", "/TTH/", "/QQH/",
    "VBF_all.root", "GGH_all.root", "GluGluH_all.root", "QQH_all.root"
  };
  for (const auto& t : tags) if (p.find(t) != std::string::npos) return true;
  return false;
}

static void collectRootFilesFromDirectory(const char* dirPath, std::vector<std::string>& paths) {
  TString base = dirPath;
  if (base.EndsWith("/")) base.Chop();

  void* dir = gSystem->OpenDirectory(base);
  if (!dir) return;

  const char* entry = nullptr;
  while ((entry = gSystem->GetDirEntry(dir))) {
    if (!strcmp(entry, ".") || !strcmp(entry, "..")) continue;
    TString sampleDir = base + "/" + entry;
    if (!isDirectory(sampleDir.Data())) continue;

    TString sampleName = entry;
    TString nominal = sampleDir + "/" + sampleName + "_all.root";
    if (!gSystem->AccessPathName(nominal, kReadPermission)) {
      paths.push_back(nominal.Data());
      continue;
    }

    void* dir2 = gSystem->OpenDirectory(sampleDir);
    const char* e2 = nullptr;
    while (dir2 && (e2 = gSystem->GetDirEntry(dir2))) {
      TString fn = e2;
      if (!fn.EndsWith(".root")) continue;
      TString full = sampleDir + "/" + fn;
      if (!gSystem->AccessPathName(full, kReadPermission)) paths.push_back(full.Data());
    }
    if (dir2) gSystem->FreeDirectory(dir2);
  }
  gSystem->FreeDirectory(dir);
}

static void collectRootFilesFromTxt(const char* txtPath, std::vector<std::string>& paths) {
  std::ifstream fin(txtPath);
  if (!fin.is_open()) return;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    paths.push_back(line);
  }
}

static int chooseColorForSample(const std::string& samplePath, int idx) {
  const std::vector<std::pair<std::string,int>> lut = {
    {"DY", kAzure-9},
    {"TT", kOrange-2},
    {"ST", kOrange+7},
    {"EWK", kGreen+2},
    {"VVV", kMagenta-7},
    {"WW", kMagenta-9},
    {"WZ", kMagenta-6},
    {"ZZ", kMagenta-4},
    {"TW", kOrange+1},
    {"QCD", kGray+1},
    {"WJets", kGreen-7},
    {"ZJets", kAzure-4}
  };
  for (const auto& kv : lut) {
    if (samplePath.find(kv.first) != std::string::npos) return kv.second;
  }
  static const int fallback[] = {kAzure-9, kOrange-2, kGreen+2, kMagenta-7, kCyan-6, kYellow-7, kPink-8};
  return fallback[idx % (sizeof(fallback)/sizeof(int))];
}

static std::string prettySampleName(const std::string& p) {
  std::string s = p;
  auto slash = s.find_last_of('/');
  if (slash != std::string::npos) s = s.substr(0, slash);
  slash = s.find_last_of('/');
  if (slash != std::string::npos) s = s.substr(slash+1);
  return s;
}

} // namespace

void MassResolution_HSBSR(const char* inputPath,
                       const char* era = "12022",
                       const char* histName = "corr_DM_Pre-Selection/DiMuonMass_Higgs") {
  using namespace RooFit;

  if (!inputPath || std::string(inputPath).empty()) {
    std::cerr << "[ERROR] inputPath is empty\n";
    return;
  }

  gStyle->SetOptStat(0);

  const double xmin = 110.0;
  const double xmax = 150.0;
  const double blindLo = 120.0;
  const double blindHi = 130.0;

  std::vector<std::string> paths;
  if (isDirectory(inputPath)) collectRootFilesFromDirectory(inputPath, paths);
  else                        collectRootFilesFromTxt(inputPath, paths);

  if (paths.empty()) {
    std::cerr << "[ERROR] No ROOT files found from input: " << inputPath << "\n";
    return;
  }

  std::unique_ptr<TH1> hDataFull;
  std::unique_ptr<TH1> hBkgFull;
  std::unique_ptr<TH1> hSigFull;
  std::vector<std::unique_ptr<TH1>> bkgComponents;
  std::vector<std::string> bkgNames;

  int compIdx = 0;
  for (const auto& p : paths) {
    std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::cerr << "[WARN] Skipping unreadable file: " << p << "\n";
      continue;
    }

    TH1* h = dynamic_cast<TH1*>(f->Get(histName));
    if (!h) {
      std::cerr << "[WARN] Missing histogram '" << histName << "' in " << p << "\n";
      continue;
    }

    std::unique_ptr<TH1> hc((TH1*)h->Clone((sanitize(prettySampleName(p))+"_clone").c_str()));
    hc->SetDirectory(nullptr);

    if (isDataFile(p)) {
      if (!hDataFull) { hDataFull.reset((TH1*)hc->Clone(uniqueName("hDataFull"))); hDataFull->SetDirectory(nullptr); }
      else hDataFull->Add(hc.get());
    } else if (isSignalFile(p)) {
      if (!hSigFull) { hSigFull.reset((TH1*)hc->Clone(uniqueName("hSigFull"))); hSigFull->SetDirectory(nullptr); }
      else hSigFull->Add(hc.get());
    } else {
      if (!hBkgFull) { hBkgFull.reset((TH1*)hc->Clone(uniqueName("hBkgFull"))); hBkgFull->SetDirectory(nullptr); }
      else hBkgFull->Add(hc.get());

      std::unique_ptr<TH1> hcomp((TH1*)hc->Clone(("hcomp_" + sanitize(prettySampleName(p))).c_str()));
      hcomp->SetDirectory(nullptr);
      hcomp->SetFillColor(chooseColorForSample(p, compIdx));
      hcomp->SetLineColor(kBlack);
      hcomp->SetLineWidth(1);
      bkgNames.push_back(prettySampleName(p));
      bkgComponents.push_back(std::move(hcomp));
      ++compIdx;
    }
  }

  if (!hDataFull) { std::cerr << "[ERROR] No DATA histogram found.\n"; return; }
  if (!hBkgFull)  { std::cerr << "[ERROR] No background histogram found.\n"; return; }
  if (!hSigFull)  { std::cerr << "[WARN] No signal histogram found. Continuing without signal overlay.\n"; }

  std::unique_ptr<TH1D> hData(makeWindowedHist(hDataFull.get(), xmin, xmax, "hData_win"));
  std::unique_ptr<TH1D> hBkg (makeWindowedHist(hBkgFull.get(),  xmin, xmax, "hBkg_win"));
  std::unique_ptr<TH1D> hSig;
  if (hSigFull) hSig.reset(makeWindowedHist(hSigFull.get(), xmin, xmax, "hSig_win"));

  std::vector<std::unique_ptr<TH1D>> bkgCompWins;
  for (size_t i = 0; i < bkgComponents.size(); ++i) {
    auto hw = std::unique_ptr<TH1D>(makeWindowedHist(bkgComponents[i].get(), xmin, xmax, ("bkgCompWin_"+std::to_string(i)).c_str()));
    if (hw) {
      hw->SetFillColor(bkgComponents[i]->GetFillColor());
      hw->SetLineColor(bkgComponents[i]->GetLineColor());
      hw->SetLineWidth(1);
      bkgCompWins.push_back(std::move(hw));
    }
  }

  if (!hData || !hBkg) {
    std::cerr << "[ERROR] Failed to build SR-window histograms.\n";
    return;
  }

  std::unique_ptr<TH1D> hDataBlind(makeBlindedDataHist(hData.get(), blindLo, blindHi, "hDataBlind"));
  std::unique_ptr<TH1D> hDataSB(makeSidebandDataHist(hData.get(), blindLo, blindHi, "hDataSB"));

  RooRealVar m("m", "m_{#mu#mu} [GeV]", xmin, xmax);
  RooDataHist dataSB("dataSB", "dataSB", RooArgList(m), Import(*hDataSB));

  // Background model inspired by the slide's BWZ x Bernstein family.
  RooRealVar mz("mz", "mZ", 91.1876);
  mz.setConstant(true);
  RooRealVar wz("wz", "GammaZ", 2.4952);
  wz.setConstant(true);
  RooBreitWigner bwz("bwz", "Breit-Wigner(Z)", m, mz, wz);

  RooRealVar c1("c1", "c1", 0.30, -5.0, 5.0);
  RooRealVar c2("c2", "c2", 0.20, -5.0, 5.0);
  RooRealVar c3("c3", "c3", 0.10, -5.0, 5.0);
  RooBernstein bern("bern", "bern", m, RooArgList(c1, c2, c3));
  RooProdPdf bkgPdf("bkgPdf", "BWZxBernstein", RooArgSet(bwz, bern));

  std::unique_ptr<RooFitResult> frBkg(
    bkgPdf.fitTo(dataSB, Save(true), PrintLevel(-1), SumW2Error(true),
                 Minimizer("Minuit2", "migrad"), Strategy(1), Range(xmin, xmax))
  );

  // Signal fit from MC with a Voigtian + CB mixture for a realistic Higgs peak.
  std::unique_ptr<RooFitResult> frSig;
  std::unique_ptr<RooPlot> sigFrame;
  double sigChi2 = -1.0;
  double meanVal = 125.0, sigmaVVal = 1.5;
  double cbSigmaVal = 1.7, fracVal = 0.6;

  RooRealVar meanH("meanH", "meanH", 125.0, 123.0, 127.0);
  RooRealVar widthH("widthH", "widthH", 0.00407);
  widthH.setConstant(true);
  RooRealVar sigmaV("sigmaV", "sigmaV", 1.4, 0.2, 5.0);
  RooVoigtian sigVoigt("sigVoigt", "sigVoigt", m, meanH, widthH, sigmaV);

  RooRealVar sigmaCB("sigmaCB", "sigmaCB", 1.8, 0.3, 6.0);
  RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.05, 10.0);
  RooRealVar nL("nL", "nL", 3.0, 0.5, 50.0);
  RooRealVar alphaR("alphaR", "alphaR", 1.8, 0.05, 10.0);
  RooRealVar nR("nR", "nR", 3.0, 0.5, 50.0);
  RooCrystalBall sigCB("sigCB", "sigCB", m, meanH, sigmaCB, alphaL, nL, alphaR, nR);

  RooRealVar frac("frac", "frac", 0.6, 0.0, 1.0);
  RooAddPdf sigPdf("sigPdf", "Voigt + CB", RooArgList(sigVoigt, sigCB), RooArgList(frac));

  if (hSig && hSig->Integral() > 0) {
    RooDataHist sigDH("sigDH", "sigDH", RooArgList(m), Import(*hSig));
    frSig.reset(sigPdf.fitTo(sigDH, Save(true), PrintLevel(-1), SumW2Error(true),
                             Minimizer("Minuit2", "migrad"), Strategy(1), Range(xmin, xmax)));
    sigFrame.reset(m.frame());
    sigDH.plotOn(sigFrame.get(), Name("sigMC"));
    sigPdf.plotOn(sigFrame.get(), Name("sigPdf"), LineWidth(2));
    sigChi2 = sigFrame->chiSquare("sigPdf", "sigMC", 7);
    meanVal = meanH.getVal();
    sigmaVVal = sigmaV.getVal();
    cbSigmaVal = sigmaCB.getVal();
    fracVal = frac.getVal();
  }

  RooPlot* frame = m.frame();
  RooDataHist dataBlindDH("dataBlindDH", "dataBlindDH", RooArgList(m), Import(*hDataBlind));
  dataBlindDH.plotOn(frame, Name("dataBlind"), DataError(RooAbsData::Poisson));
  bkgPdf.plotOn(frame, Name("bkgFit"), LineWidth(3));
  if (hSig && hSig->Integral() > 0) sigPdf.plotOn(frame, Name("sigFit"), LineStyle(kDashed), LineWidth(3));

  const double chi2ndf_bkg = frame->chiSquare("bkgFit", "dataBlind", 3);

  TString corrections = guessCorrectionsTag(inputPath);
  TString corrType    = guessCorrTypeTag(inputPath);
  TString outdir = "/eos/user/d/dblend/StackedPlots/" + TString(era) + "/Mass_Resolutions/" + corrType + "/" + corrections + "/";
  gSystem->Exec("mkdir -p " + outdir);

  TCanvas* c = new TCanvas(uniqueName("c_massres_sr"), "SR mass resolution", 900, 800);
  TPad* pad1 = new TPad(uniqueName("pad1"), "top", 0.0, 0.30, 1.0, 1.0);
  TPad* pad2 = new TPad(uniqueName("pad2"), "bottom", 0.0, 0.00, 1.0, 0.30);
  pad1->SetBottomMargin(0.03);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.30);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);

  THStack* hs = new THStack(uniqueName("hs"), "");
  for (auto& hcomp : bkgCompWins) hs->Add(hcomp.get(), "HIST");
  hs->Draw("HIST");
  hs->SetTitle("");
  hs->GetYaxis()->SetTitle("Events / bin");
  hs->GetYaxis()->SetTitleOffset(1.35);
  hs->GetXaxis()->SetLabelSize(0);

  TH1D* hBkgOutline = (TH1D*)hBkg->Clone(uniqueName("hBkgOutline"));
  hBkgOutline->SetDirectory(nullptr);
  hBkgOutline->SetFillStyle(0);
  hBkgOutline->SetLineColor(kBlack);
  hBkgOutline->SetLineWidth(2);
  hBkgOutline->Draw("HIST SAME");

  TH1D* hSigDraw = nullptr;
  double sigScale = 1.0;
  if (hSig && hSig->Integral() > 0) {
    hSigDraw = (TH1D*)hSig->Clone(uniqueName("hSigDraw"));
    hSigDraw->SetDirectory(nullptr);
    const double maxBkg = std::max(1.0, hBkg->GetMaximum());
    const double maxSig = std::max(1e-9, hSigDraw->GetMaximum());
    sigScale = 0.35 * maxBkg / maxSig;
    hSigDraw->Scale(sigScale);
    hSigDraw->SetLineColor(kRed+1);
    hSigDraw->SetLineWidth(3);
    hSigDraw->SetFillStyle(0);
    hSigDraw->Draw("HIST SAME");
  }

  frame->Draw("SAME");

  const double yMax = std::max({hs->GetMaximum(), hData->GetMaximum(), hBkg->GetMaximum(), 1.0}) * 1.45;
  hs->SetMaximum(yMax);
  hs->SetMinimum(0.0);

  const double blindBoxHeight = yMax * 0.92;
  TBox* blindBox = new TBox(blindLo, 0.0, blindHi, blindBoxHeight);
  blindBox->SetFillColorAlpha(kGray+1, 0.35);
  blindBox->SetLineColor(kGray+2);
  blindBox->SetLineStyle(2);
  blindBox->Draw("SAME");

  TLine* lblind1 = new TLine(blindLo, 0.0, blindLo, blindBoxHeight);
  TLine* lblind2 = new TLine(blindHi, 0.0, blindHi, blindBoxHeight);
  lblind1->SetLineStyle(2); lblind2->SetLineStyle(2);
  lblind1->Draw("SAME");  lblind2->Draw("SAME");

  TLegend* leg = new TLegend(0.58, 0.52, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  TObject* objDataBlind = frame->findObject("dataBlind");
  TObject* objBkgFit    = frame->findObject("bkgFit");
  TObject* objSigFit    = frame->findObject("sigFit");
  if (objDataBlind) leg->AddEntry(objDataBlind, "Data (blinded 120-130 GeV)", "lep");
  leg->AddEntry(hBkgOutline, "Total background MC", "l");
  if (hSigDraw) leg->AddEntry(hSigDraw, Form("Signal MC #times %.1f", sigScale), "l");
  if (objBkgFit) leg->AddEntry(objBkgFit, "BWZ #times Bernstein fit", "l");
  if (hSig && objSigFit) leg->AddEntry(objSigFit, "Signal fit (Voigt + CB)", "l");
  for (size_t i = 0; i < bkgCompWins.size() && i < bkgNames.size(); ++i) {
    leg->AddEntry(bkgCompWins[i].get(), bkgNames[i].c_str(), "f");
  }
  leg->Draw();

  TPaveText* pt = new TPaveText(0.15, 0.60, 0.48, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->AddText(Form("CMS Preliminary   %.2f fb^{-1} (13.6 TeV)", eraLumiPb(era)/1000.0));
  pt->AddText(Form("Background fit #chi^{2}/ndf = %.3f", chi2ndf_bkg));
  pt->AddText("Fit range: 110-150 GeV");
  pt->AddText("Blinded: 120-130 GeV");
  if (hSig && hSig->Integral() > 0) {
    pt->AddText(Form("Signal mean = %.3f GeV", meanVal));
    pt->AddText(Form("Voigt #sigma = %.3f, CB #sigma = %.3f", sigmaVVal, cbSigmaVal));
    pt->AddText(Form("Signal frac(Voigt) = %.3f", fracVal));
    if (sigChi2 > 0) pt->AddText(Form("Signal fit #chi^{2}/ndf = %.3f", sigChi2));
  }
  pt->Draw();

  pad2->cd();
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);

  TH1D* hRatio = (TH1D*)hDataBlind->Clone(uniqueName("hRatio"));
  hRatio->SetDirectory(nullptr);
  for (int b = 1; b <= hRatio->GetNbinsX(); ++b) {
    const double mc = hBkg->GetBinContent(b);
    const double x  = hRatio->GetBinCenter(b);
    if (mc <= 0 || (x >= blindLo && x <= blindHi)) {
      hRatio->SetBinContent(b, 0.0);
      hRatio->SetBinError(b, 0.0);
      continue;
    }
  }
  hRatio->Divide(hBkg.get());
  hRatio->SetTitle("");
  hRatio->GetYaxis()->SetTitle("Data/MC");
  hRatio->GetYaxis()->SetTitleSize(0.11);
  hRatio->GetYaxis()->SetTitleOffset(0.45);
  hRatio->GetYaxis()->SetLabelSize(0.09);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetLabelSize(0.11);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);
  hRatio->Draw("E1");

  TBox* blindRatio = new TBox(blindLo, 0.5, blindHi, 1.5);
  blindRatio->SetFillColorAlpha(kGray+1, 0.35);
  blindRatio->SetLineColor(kGray+2);
  blindRatio->SetLineStyle(2);
  blindRatio->Draw("SAME");

  TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
  l1->SetLineStyle(2);
  l1->SetLineWidth(2);
  l1->Draw("SAME");

  TString baseName = "MassResolution_SR_Blinded";
  c->SaveAs(outdir + baseName + ".pdf");
  c->SaveAs(outdir + baseName + ".png");

  std::cout << "[INFO] Saved plots to: " << outdir << std::endl;
  std::cout << "[INFO] Background model: BWZ x Bernstein" << std::endl;
  if (hSig && hSig->Integral() > 0) {
    std::cout << "[INFO] Signal model: Voigt + CrystalBall" << std::endl;
    std::cout << "[INFO] Signal mean=" << meanVal
              << "  sigmaV=" << sigmaVVal
              << "  sigmaCB=" << cbSigmaVal
              << "  frac=" << fracVal << std::endl;
  }

  // Do not manually delete drawables owned by the canvas/pads.
  // ROOT object ownership can otherwise cause double-deletion / browsing crashes.
  c->Modified();
  c->Update();
}

void MassResolution_SR(const char* inputPath,
                       const char* era = "12022",
                       const char* histName = "corr_DM_Pre-Selection/DiMuonMass_Higgs") {
  MassResolution_HSBSR(inputPath, era, histName);
}
