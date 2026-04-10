// Example:
// root -l
// .L MassResolution_ZCR.C
// MassResolution_ZCR("/eos/user/d/dblend/Analyzed_SkimmedFiles/12022/KIT/noJERC_noJetVeto_noBtagSF_noZPtRW/", "12022");
#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TDirectory.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TString.h"
#include "TPad.h"
#include "TLine.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooVoigtian.h"
#include "RooCrystalBall.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

namespace {

static bool isDirectory(const char* p) {
  if (!p) return false;
  void* dir = gSystem->OpenDirectory(p);
  if (dir) { gSystem->FreeDirectory(dir); return true; }
  return false;
}

// --- Lumi lookup (pb) – copied from your macro
static double eraLumiPb(const std::string& era) {
  if      (era == "12022") return  7980.454145;
  else if (era == "22022") return 26671.60971;
  else if (era == "12023") return 18704.13341;
  else if (era == "22023") return  9693.130053;
  else if (era == "2024")  return 104805.9585;
  std::cerr << "[WARN] Unknown era \"" << era << "\"; setting lumi=0\n";
  return 0.0;
}

// Helper: window a TH1 into [xmin, xmax] with identical bin edges inside that window
static TH1D* makeWindowedHist(const TH1* hIn, double xmin, double xmax) {
  if (!hIn) return nullptr;

  const int b1 = std::max(1, hIn->GetXaxis()->FindBin(xmin));
  const int b2 = std::min(hIn->GetNbinsX(), hIn->GetXaxis()->FindBin(xmax));
  const int nb = b2 - b1 + 1;

  std::vector<double> edges; edges.reserve(nb + 1);
  for (int b = b1; b <= b2 + 1; ++b) edges.push_back(hIn->GetXaxis()->GetBinLowEdge(b));

  auto* h = new TH1D((std::string(hIn->GetName()) + "_win").c_str(), hIn->GetTitle(), nb, edges.data());
  h->Sumw2();

  for (int i = 0; i < nb; ++i) {
    const int b = b1 + i;
    h->SetBinContent(i+1, hIn->GetBinContent(b));
    h->SetBinError(i+1,   hIn->GetBinError(b));
  }

  h->GetXaxis()->SetTitle(hIn->GetXaxis()->GetTitle());
  h->GetYaxis()->SetTitle(hIn->GetYaxis()->GetTitle());
  return h;
}

// Numerical quantile from a CDF via bisection
static double quantileFromCdf(RooAbsReal& cdf, RooRealVar& x, double q, double xmin, double xmax) {
  double lo = xmin, hi = xmax;
  for (int it = 0; it < 80; ++it) {
    double mid = 0.5*(lo+hi);
    x.setVal(mid);
    double v = cdf.getVal();
    if (v < q) lo = mid;
    else       hi = mid;
  }
  return 0.5*(lo+hi);
}

// FWHM from a pdf by sampling (same approach as your macro)
static double fwhmFromPdf(RooAbsPdf& pdf, RooRealVar& x, double xmin, double xmax) {
  RooArgSet normSet(x);

  const int N = 5000;
  double maxv = -1;
  double xmaxAt = xmin;
  for (int i = 0; i < N; ++i) {
    double xv = xmin + (xmax-xmin) * (i + 0.5) / N;
    x.setVal(xv);
    double v = pdf.getVal(normSet);
    if (v > maxv) { maxv = v; xmaxAt = xv; }
  }
  if (maxv <= 0) return 0.0;
  const double half = 0.5 * maxv;

  double left = xmin;
  for (int i = 0; i < N; ++i) {
    double xv = xmaxAt - (xmaxAt-xmin) * (i + 0.5) / N;
    x.setVal(xv);
    if (pdf.getVal(normSet) < half) { left = xv; break; }
  }

  double right = xmax;
  for (int i = 0; i < N; ++i) {
    double xv = xmaxAt + (xmax-xmaxAt) * (i + 0.5) / N;
    x.setVal(xv);
    if (pdf.getVal(normSet) < half) { right = xv; break; }
  }

  return std::max(0.0, right - left);
}

// Simple sanitize for filenames
static std::string sanitize(std::string s) {
  std::replace(s.begin(), s.end(), '/', '_');
  std::replace(s.begin(), s.end(), ' ', '_');
  return s;
}

} // namespace

// Usage:
//   root -l -q 'MassResolution_ZCR_Voigt.C("/eos/.../DY/DY_all.root","12022")'
void MassResolution_ZCR(const char* fileListTxt,
                                const char* era = "12022") {  
  using namespace RooFit;

    if (!fileListTxt || std::string(fileListTxt).empty()) {
      std::cerr << "[ERROR] fileListTxt is empty\n";
      return;
    }

    std::vector<std::string> paths;

// Case 1: user passed a directory like .../noZPtRW
if (isDirectory(fileListTxt)) {
  TString base = fileListTxt;
  if (base.EndsWith("/")) base.Chop();

  void* dir = gSystem->OpenDirectory(base);
  const char* entry = nullptr;

  while ((entry = gSystem->GetDirEntry(dir))) {
    if (!strcmp(entry,".") || !strcmp(entry,"..")) continue;

    TString sampleDir = base + "/" + entry;
    if (!isDirectory(sampleDir.Data())) continue;

    // Expect <sampleDir>/<sampleName>_all.root
    TString sampleName = entry;
    TString fpath = sampleDir + "/" + sampleName + "_all.root";

    // Only add if file exists/readable
    if (!gSystem->AccessPathName(fpath, kReadPermission)) {
      paths.push_back(fpath.Data());
    } else {
      // (optional) fallback: accept any *.root inside the sample dir
      void* dir2 = gSystem->OpenDirectory(sampleDir);
      const char* e2 = nullptr;
      while ((e2 = gSystem->GetDirEntry(dir2))) {
        TString fn = e2;
        if (!fn.EndsWith(".root")) continue;
        TString full = sampleDir + "/" + fn;
        if (!gSystem->AccessPathName(full, kReadPermission))
          paths.push_back(full.Data());
      }
      gSystem->FreeDirectory(dir2);
    }
  }
  gSystem->FreeDirectory(dir);

} else {
  // Case 2: treat input as a text file list
  std::ifstream fin(fileListTxt);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Could not open file list: " << fileListTxt << "\n";
    return;
  }
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    paths.push_back(line);
  }
}

if (paths.empty()) {
  std::cerr << "[ERROR] No ROOT files found from input: " << fileListTxt << "\n";
  return;
}

  const std::string eraStr = era ? std::string(era) : std::string("?");
  const double lumi_pb = eraLumiPb(eraStr);
  const double lumi_fb = lumi_pb / 1000.0;

  // ---- reproduce your output-dir logic (same style as signal macro)
  TString tmp = paths[0].c_str();
  if (tmp.EndsWith("/")) tmp.Chop();

  tmp = tmp(0, tmp.Last('/')); // remove filename
  tmp = tmp(0, tmp.Last('/')); // remove sample dir

  TString corrections = tmp(tmp.Last('/') + 1, tmp.Length());
  tmp = tmp(0, tmp.Last('/'));
  TString corrType = tmp(tmp.Last('/') + 1, tmp.Length());

    TString outdir =
    "/eos/user/d/dblend/StackedPlots/" + TString(eraStr.c_str()) +
    "/Mass_Resolutions/ZCR/" + corrType + "/" + corrections + "/";
    gSystem->Exec("mkdir -p " + outdir);

  // ---- Z-window
  const double xmin = 75.0;
  const double xmax = 105.0;

auto fitAndPlot = [&](TH1* hDataIn, TH1* hMcIn, const std::string& histPathTag) {
  if (!hDataIn || !hMcIn) {
    std::cerr << "[WARN] Missing Data/MC hist for tag '" << histPathTag << "'\n";
    return;
  }

  std::unique_ptr<TH1D> hDataWin(makeWindowedHist(hDataIn, xmin, xmax));
  std::unique_ptr<TH1D> hMcWin  (makeWindowedHist(hMcIn,   xmin, xmax));
  if (!hDataWin || !hMcWin) {
    std::cerr << "[ERROR] Failed windowing for tag '" << histPathTag << "'\n";
    return;
  }

  if (hDataWin->GetNbinsX() > 600) { hDataWin->Rebin(2); hMcWin->Rebin(2); }

  // ---- RooFit on DATA (so your “MassResolution” reflects data Z peak)
  // ---- RooFit with DCB on DATA and MC
const double xlo = hDataWin->GetXaxis()->GetXmin();
const double xhi = hDataWin->GetXaxis()->GetXmax();

RooRealVar m("m", "m_{#mu#mu} [GeV]", xlo, xhi);

// Datasets
RooDataHist dataDH("dataDH", "dataDH", RooArgList(m), Import(*hDataWin));
RooDataHist mcDH  ("mcDH",   "mcDH",   RooArgList(m), Import(*hMcWin));

// --------------------
// DATA DCB parameters
// --------------------
RooRealVar meanD ("meanD",  "meanD",  91.1876, 90.5, 91.9);
RooRealVar sigmaD("sigmaD", "sigmaD",  1.3,     0.2,  5.0);

RooRealVar alphaLD("alphaLD","alphaLD", 1.2, 0.05, 10.0);
RooRealVar nLD    ("nLD",    "nLD",     3.0, 0.5,  200.0);
RooRealVar alphaRD("alphaRD","alphaRD", 1.2, 0.05, 10.0);
RooRealVar nRD    ("nRD",    "nRD",     3.0, 0.5,  200.0);

RooCrystalBall dcbData("dcbData", "DCB(Data)", m, meanD, sigmaD, alphaLD, nLD, alphaRD, nRD);

std::unique_ptr<RooFitResult> frD(dcbData.fitTo(
  dataDH, Save(true), PrintLevel(-1), SumW2Error(true),
  Minimizer("Minuit2","migrad"), Strategy(1)
));

// --------------------
// MC DCB parameters
// (separate params so you can print both sets like you want)
// --------------------
RooRealVar meanM ("meanM",  "meanM",  91.1876, 90.5, 91.9);
RooRealVar sigmaM("sigmaM", "sigmaM",  1.1,     0.2,  5.0);

RooRealVar alphaLM("alphaLM","alphaLM", 1.2, 0.05, 10.0);
RooRealVar nLM    ("nLM",    "nLM",     3.0, 0.5,  200.0);
RooRealVar alphaRM("alphaRM","alphaRM", 1.2, 0.05, 10.0);
RooRealVar nRM    ("nRM",    "nRM",     3.0, 0.5,  200.0);

RooCrystalBall dcbMC("dcbMC", "DCB(MC)", m, meanM, sigmaM, alphaLM, nLM, alphaRM, nRM);

std::unique_ptr<RooFitResult> frM(dcbMC.fitTo(
  mcDH, Save(true), PrintLevel(-1), SumW2Error(true),
  Minimizer("Minuit2","migrad"), Strategy(1)
));

// Frames
auto* frame = m.frame(Title(("Z control region fit (DCB): " + histPathTag).c_str()));
dataDH.plotOn(frame, Name("data"));
dcbData.plotOn(frame, LineWidth(2), Name("dcbData"));

// Optional: overlay MC fit curve too
dcbMC.plotOn(frame, LineWidth(2), LineStyle(kDashed), Name("dcbMC"));

// chi2/ndf (count params: mean,sigma,alphaL,nL,alphaR,nR = 6)
const double chi2ndfD = frame->chiSquare("dcbData", "data", 6);

// For DATA resolution metrics
std::unique_ptr<RooAbsReal> cdfD(dcbData.createCdf(m));
const double q16D = quantileFromCdf(*cdfD, m, 0.16, xlo, xhi);
const double q84D = quantileFromCdf(*cdfD, m, 0.84, xlo, xhi);
const double fwhmD = fwhmFromPdf(dcbData, m, xlo, xhi);
const double sigma68D = 0.5 * (q84D - q16D);
const double sigmaFromFwhmD = (fwhmD > 0 ? fwhmD / 2.354820045 : 0.0);

// For MC resolution metrics
std::unique_ptr<RooAbsReal> cdfM(dcbMC.createCdf(m));
const double q16M = quantileFromCdf(*cdfM, m, 0.16, xlo, xhi);
const double q84M = quantileFromCdf(*cdfM, m, 0.84, xlo, xhi);
const double fwhmM = fwhmFromPdf(dcbMC, m, xlo, xhi);
const double sigma68M = 0.5 * (q84M - q16M);
const double sigmaFromFwhmM = (fwhmM > 0 ? fwhmM / 2.354820045 : 0.0);
  
  // ---- 2-pad canvas (pattern copied from overlay) :contentReference[oaicite:4]{index=4}
  auto* c = new TCanvas(("c_"+sanitize(histPathTag)).c_str(), "c", 900, 800);

  TPad* pad1 = new TPad(("pad1_"+sanitize(histPathTag)).c_str(), "Top",    0, 0.32, 1, 1);
  TPad* pad2 = new TPad(("pad2_"+sanitize(histPathTag)).c_str(), "Bottom", 0, 0.02, 1, 0.32);
  pad1->SetBottomMargin(0.03);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.30);
  pad1->Draw(); pad2->Draw();

  // ---- TOP: RooPlot + also draw MC shape behind (optional)
  pad1->cd();
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.05);

  // Draw MC as a histogram in the same pad, scaled to Data integral in the window
  TH1D* hMcDraw = (TH1D*)hMcWin->Clone("hMcDraw");
  hMcDraw->SetDirectory(nullptr);
  const double iData = hDataWin->Integral(1, hDataWin->GetNbinsX());
  const double iMc   = hMcWin->Integral(1, hMcWin->GetNbinsX());
  if (iMc > 0) hMcDraw->Scale(iData / iMc);

  hMcDraw->SetFillStyle(3004);
  hMcDraw->SetLineWidth(2);
  hMcDraw->SetTitle("");
  hMcDraw->GetYaxis()->SetTitle("Events");
  hMcDraw->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  hMcDraw->GetYaxis()->SetTitleOffset(1.35);
  hMcDraw->Draw("HIST");

  frame->Draw("SAME"); // Data points + Voigt curve on top

  auto* pt = new TPaveText(0.58, 0.60, 0.90, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->AddText("=== DATA (DCB) ===");
pt->AddText(Form("mean = %.4f", meanD.getVal()));
pt->AddText(Form("sigma = %.4f", sigmaD.getVal()));
pt->AddText(Form("aL=%.3f nL=%.2f  aR=%.3f nR=%.2f",
                 alphaLD.getVal(), nLD.getVal(), alphaRD.getVal(), nRD.getVal()));
pt->AddText(Form("sigma_{68} = %.4f", sigma68D));
pt->AddText(Form("FWHM = %.4f (#sigma=%.4f)", fwhmD, sigmaFromFwhmD));
pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2ndfD));

pt->AddText("=== MC (DCB) ===");
pt->AddText(Form("mean = %.4f", meanM.getVal()));
pt->AddText(Form("sigma = %.4f", sigmaM.getVal()));
pt->AddText(Form("aL=%.3f nL=%.2f  aR=%.3f nR=%.2f",
                 alphaLM.getVal(), nLM.getVal(), alphaRM.getVal(), nRM.getVal()));
pt->AddText(Form("sigma_{68} = %.4f", sigma68M));
pt->AddText(Form("FWHM = %.4f (#sigma=%.4f)", fwhmM, sigmaFromFwhmM));
  pt->Draw("same");

  TLegend leg(0.15, 0.75, 0.45, 0.88);
  leg.AddEntry(frame->findObject("data"), "Data", "lep");
leg.AddEntry(hMcDraw, "MC (sum, norm. to Data)", "f");
leg.AddEntry(frame->findObject("dcbData"), "DCB fit (Data)", "l");
leg.AddEntry(frame->findObject("dcbMC"), "DCB fit (MC)", "l");
  leg.Draw();

  // ---- BOTTOM: Data/MC ratio using the *same binning window*
  pad2->cd();
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.05);

  TH1D* hRatio = (TH1D*)hDataWin->Clone("hRatio");
  hRatio->SetDirectory(nullptr);

  // Use unscaled MC for ratio (ratio is independent of overall normalization choice)
  for (int i=1; i<=hRatio->GetNbinsX(); ++i) {
    if (hMcWin->GetBinContent(i) <= 0) { hRatio->SetBinContent(i, 0); hRatio->SetBinError(i, 0); }
  }
  hRatio->Divide(hMcWin.get());

  hRatio->SetTitle("");
  hRatio->GetYaxis()->SetTitle("Data/MC");
  hRatio->GetYaxis()->SetTitleSize(0.12);
  hRatio->GetYaxis()->SetTitleOffset(0.40);
  hRatio->GetYaxis()->SetLabelSize(0.10);
  hRatio->GetYaxis()->SetNdivisions(505);

  hRatio->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  hRatio->GetXaxis()->SetTitleSize(0.12);
  hRatio->GetXaxis()->SetLabelSize(0.12);

  hRatio->SetMinimum(0.9);
  hRatio->SetMaximum(1.1);
  hRatio->SetMarkerStyle(20);
  hRatio->Draw("E1");

  pad2->Update();
  double rxmin = gPad->GetUxmin();
  double rxmax = gPad->GetUxmax();
  TLine* one = new TLine(rxmin, 1.0, rxmax, 1.0);
  one->SetLineStyle(2);
  one->SetLineWidth(2);
  one->Draw("SAME");

  const std::string tag = sanitize(histPathTag);
  const TString baseName = "MassResolution_ZCR_" + TString(tag.c_str());
  c->SaveAs(outdir + baseName + ".pdf");
  c->SaveAs(outdir + baseName + ".png");

  delete one;
  delete hRatio;
  delete hMcDraw;
  delete pt;
  delete c;
};

const char* histName = "corr_DM_Pre-Selection/DiMuonMass_Z";

std::unique_ptr<TH1> h_mc;    // total MC background
std::unique_ptr<TH1> h_data;  // data

auto isDataFile = [](const std::string& p){
  // directory mode produces .../DATA/DATA_all.root; this is robust enough
  return (p.find("/DATA/") != std::string::npos) || (p.find("DATA_all.root") != std::string::npos);
};

for (const auto& p : paths) {

  std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "[WARN] Skipping unreadable file: " << p << "\n";
    continue;
  }

  TH1* h = dynamic_cast<TH1*>(f->Get(histName));
  if (!h) {
    std::cerr << "[WARN] Missing " << histName << " in: " << p << "\n";
    continue;
  }

  TH1* hc = dynamic_cast<TH1*>(h->Clone());
  hc->SetDirectory(nullptr);

  if (isDataFile(p)) {
    if (!h_data) {
      h_data.reset(dynamic_cast<TH1*>(hc->Clone("h_data_sum")));
      h_data->SetDirectory(nullptr);
    } else {
      h_data->Add(hc);
    }
  } else {
    if (!h_mc) {
      h_mc.reset(dynamic_cast<TH1*>(hc->Clone("h_mc_sum")));
      h_mc->SetDirectory(nullptr);
    } else {
      h_mc->Add(hc);
    }
  }

  delete hc;
}

if (!h_mc)   { std::cerr << "[ERROR] No MC histograms were loaded.\n"; return; }
if (!h_data) { std::cerr << "[ERROR] No DATA histogram was loaded. Make sure DATA/DATA_all.root is included.\n"; return; }

fitAndPlot(h_data.get(), h_mc.get(), "ZCR_corr_DM_Pre-Selection_DiMuonMass_Z");
}
