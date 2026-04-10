// Example:
// root -l
// .L MassResolution_BackgroundEstimation.C
// MassResolution_BackgroundEstimation("","2024")
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TLine.h"
#include "TBox.h"
#include "TStyle.h"
#include "TString.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooBreitWigner.h"
#include "RooExtendPdf.h"
#include "RooArgList.h"
#include "RooAbsReal.h"

#include <algorithm>
#include <fstream>
#include <iostream>
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
    h->SetBinContent(i + 1, hIn->GetBinContent(b));
    h->SetBinError(i + 1, hIn->GetBinError(b));
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
  s = s(0, s.Last('/'));
  return s(s.Last('/') + 1, s.Length());
}

static bool isDataFile(const std::string& p) {
  return (p.find("/DATA/") != std::string::npos) || (p.find("DATA_all.root") != std::string::npos);
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

  const double xmin    = 115.0;
  const double xmax    = 135.0;
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
  for (const auto& p : paths) {
    if (!isDataFile(p)) continue;

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

    std::unique_ptr<TH1> hc((TH1*)h->Clone((sanitize(p)+"_clone").c_str()));
    hc->SetDirectory(nullptr);
    if (!hDataFull) {
      hDataFull.reset((TH1*)hc->Clone(uniqueName("hDataFull")));
      hDataFull->SetDirectory(nullptr);
    } else {
      hDataFull->Add(hc.get());
    }
  }

  if (!hDataFull) {
    std::cerr << "[ERROR] No DATA histogram found.\n";
    return;
  }

  std::unique_ptr<TH1D> hData(makeWindowedHist(hDataFull.get(), xmin, xmax, "hData_win"));
  if (!hData) {
    std::cerr << "[ERROR] Failed to build 110-150 GeV data histogram.\n";
    return;
  }

  std::unique_ptr<TH1D> hDataBlind(makeBlindedDataHist(hData.get(), blindLo, blindHi, "hDataBlind"));

  RooRealVar m("m", "m_{#mu#mu} [GeV]", xmin, xmax);
  m.setRange("leftSB",  xmin,    blindLo);
  m.setRange("rightSB", blindHi, xmax);
  m.setRange("full",    xmin,    xmax);
  m.setRange("blind",   blindLo, blindHi);

  RooDataHist dataSB("dataSB", "data sidebands", RooArgList(m), Import(*hDataBlind));

  RooRealVar meanBW("meanBW", "Breit-Wigner mean", 91.1876, 80.0, 100.0);
  RooRealVar widthBW("widthBW", "Breit-Wigner width", 2.4952, 0.5, 10.0);
  RooBreitWigner bwPdf("bwPdf", "Breit-Wigner background", m, meanBW, widthBW);

  const double nSideband = hDataBlind->Integral();
  RooRealVar nbkg("nbkg", "background yield in sidebands", nSideband, 0.0, 10.0 * std::max(1.0, nSideband));
  RooExtendPdf extBkgPdf("extBkgPdf", "extended BW background", bwPdf, nbkg);

  std::unique_ptr<RooFitResult> frBkg(
    extBkgPdf.fitTo(dataSB,
                    Save(true),
                    Extended(true),
                    PrintLevel(-1),
                    SumW2Error(true),
                    Minimizer("Minuit2", "migrad"),
                    Strategy(1),
                    Range("leftSB,rightSB"))
  );

  std::unique_ptr<RooAbsReal> intLeft (bwPdf.createIntegral(m, NormSet(m), Range("leftSB")));
  std::unique_ptr<RooAbsReal> intRight(bwPdf.createIntegral(m, NormSet(m), Range("rightSB")));
  std::unique_ptr<RooAbsReal> intBlind(bwPdf.createIntegral(m, NormSet(m), Range("blind")));
  std::unique_ptr<RooAbsReal> intFull (bwPdf.createIntegral(m, NormSet(m), Range("full")));

  const double fracLeft  = intLeft  ? intLeft->getVal()  : 0.0;
  const double fracRight = intRight ? intRight->getVal() : 0.0;
  const double fracBlind = intBlind ? intBlind->getVal() : 0.0;
  const double fracSB    = fracLeft + fracRight;

  const double nPredBlind = (fracSB > 0.0) ? nSideband * (fracBlind / fracSB) : 0.0;

  TH1D* hFitBinned = (TH1D*)hData->Clone(uniqueName("hFitBinned"));
  hFitBinned->SetDirectory(nullptr);
  hFitBinned->Reset("ICESM");
  hFitBinned->SetLineColor(kBlue+1);
  hFitBinned->SetLineWidth(3);
  hFitBinned->SetFillStyle(0);

  for (int b = 1; b <= hFitBinned->GetNbinsX(); ++b) {
    const double lo = hFitBinned->GetXaxis()->GetBinLowEdge(b);
    const double hi = hFitBinned->GetXaxis()->GetBinLowEdge(b + 1);
    m.setRange("tmpBin", lo, hi);
    std::unique_ptr<RooAbsReal> iBin(bwPdf.createIntegral(m, NormSet(m), Range("tmpBin")));
    const double fracBin = iBin ? iBin->getVal() : 0.0;
    const double pred = (fracSB > 0.0) ? nSideband * (fracBin / fracSB) : 0.0;
    hFitBinned->SetBinContent(b, pred);
    hFitBinned->SetBinError(b, 0.0);
  }

  RooPlot* frame = m.frame();
  dataSB.plotOn(frame, Name("dataSB"), DataError(RooAbsData::Poisson));
  extBkgPdf.plotOn(frame, Name("bwFit"), Range("full"), NormRange("leftSB,rightSB"), LineWidth(3), LineColor(kBlue+1));
  const double chi2ndf_bkg = frame->chiSquare("bwFit", "dataSB", 3);

  TString corrections = guessCorrectionsTag(inputPath);
  TString corrType    = guessCorrTypeTag(inputPath);
  TString outdir = "/eos/user/d/dblend/StackedPlots/" + TString(era) + "/Mass_Resolution/" + corrType + "/" + corrections + "/";
  gSystem->Exec("mkdir -p " + outdir);

  TCanvas* c = new TCanvas(uniqueName("c_massres_bw_sb"), "Sideband BW background fit", 900, 800);
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

  hDataBlind->SetTitle("");
  hDataBlind->SetMarkerStyle(20);
  hDataBlind->SetMarkerSize(1.0);
  hDataBlind->SetLineColor(kBlack);
  hDataBlind->GetYaxis()->SetTitle("Events / bin");
  hDataBlind->GetYaxis()->SetTitleOffset(1.35);
  hDataBlind->GetXaxis()->SetLabelSize(0);
  hDataBlind->SetMaximum(std::max(1.0, 1.45 * std::max(hData->GetMaximum(), hFitBinned->GetMaximum())));
  hDataBlind->SetMinimum(0.0);
  hDataBlind->Draw("E1");
  hFitBinned->Draw("HIST SAME");

  TBox* blindBox = new TBox(blindLo, 0.0, blindHi, hDataBlind->GetMaximum()*0.92);
  blindBox->SetFillColorAlpha(kGray+1, 0.35);
  blindBox->SetLineColor(kGray+2);
  blindBox->SetLineStyle(2);
  blindBox->Draw("SAME");

  TLine* lblind1 = new TLine(blindLo, 0.0, blindLo, hDataBlind->GetMaximum()*0.92);
  TLine* lblind2 = new TLine(blindHi, 0.0, blindHi, hDataBlind->GetMaximum()*0.92);
  lblind1->SetLineStyle(2);
  lblind2->SetLineStyle(2);
  lblind1->Draw("SAME");
  lblind2->Draw("SAME");

  TLegend* leg = new TLegend(0.58, 0.68, 0.89, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hDataBlind.get(), "Data sidebands", "lep");
  leg->AddEntry(hFitBinned, "Breit-Wigner fit", "l");
  leg->Draw();

  TPaveText* pt = new TPaveText(0.15, 0.58, 0.50, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->AddText(Form("CMS Preliminary   %.2f fb^{-1} (13.6 TeV)", eraLumiPb(era)/1000.0));
  pt->AddText("Fit region: 110-120 #cup 130-150 GeV");
  pt->AddText("Blinded region: 120-130 GeV");
  pt->AddText(Form("BW mean = %.3f GeV", meanBW.getVal()));
  pt->AddText(Form("BW width = %.3f GeV", widthBW.getVal()));
  pt->AddText(Form("Fit #chi^{2}/ndf = %.3f", chi2ndf_bkg));
  pt->AddText(Form("Predicted blinded yield = %.2f", nPredBlind));
  pt->Draw();

  pad2->cd();
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);

  TH1D* hRatio = (TH1D*)hDataBlind->Clone(uniqueName("hRatio"));
  hRatio->SetDirectory(nullptr);
  for (int b = 1; b <= hRatio->GetNbinsX(); ++b) {
    const double fitv = hFitBinned->GetBinContent(b);
    const double x = hRatio->GetBinCenter(b);
    if (fitv <= 0.0 || (x >= blindLo && x <= blindHi)) {
      hRatio->SetBinContent(b, 0.0);
      hRatio->SetBinError(b, 0.0);
      continue;
    }
  }
  hRatio->Divide(hFitBinned);
  hRatio->SetTitle("");
  hRatio->GetYaxis()->SetTitle("Data/Fit");
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

  TString baseName = "MassResolution_HSBSR_BreitWignerOnly";
  c->SaveAs(outdir + baseName + ".pdf");
  c->SaveAs(outdir + baseName + ".png");

  std::cout << "[INFO] Saved plots to: " << outdir << std::endl;
  std::cout << "[INFO] Background model: Breit-Wigner only" << std::endl;
  std::cout << "[INFO] Fit range: 110-120 U 130-150 GeV" << std::endl;
  std::cout << "[INFO] Blinded region: 120-130 GeV" << std::endl;
  std::cout << "[INFO] Fitted BW mean = " << meanBW.getVal() << std::endl;
  std::cout << "[INFO] Fitted BW width = " << widthBW.getVal() << std::endl;
  std::cout << "[INFO] Predicted blinded yield = " << nPredBlind << std::endl;

  c->Modified();
  c->Update();
}

void MassResolution_SR(const char* inputPath,
                       const char* era = "12022",
                       const char* histName = "corr_DM_Pre-Selection/DiMuonMass_Higgs") {
  MassResolution_HSBSR(inputPath, era, histName);
}
