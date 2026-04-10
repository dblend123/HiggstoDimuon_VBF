//Example Run:
//root -l
//.L MassResolution_SCR.C
//MassResolution_SCR("/eos/user/d/dblend/Analyzed_SkimmedFiles/2024/KIT/noJERC_noJetVeto_noBtagSF_noZPtRW/VBF/VBF_all.root", "2024")
#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TSystem.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooCrystalBall.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

namespace {

// --- Lumi lookup (pb) – from your original values
static double eraLumiPb(const std::string& era) {
  if      (era == "12022") return  7980.454145;
  else if (era == "22022") return 26671.60971;
  else if (era == "12023") return 18704.13341;
  else if (era == "22023") return  9693.130053;
  else if (era == "2024")  return 104805.9585;
  std::cerr << "[WARN] Unknown era \"" << era << "\"; setting lumi=0\n";
  return 0.0;
}

// Recursively search inside directories for a histogram whose name matches common patterns.
static TH1* findMassHistRecursive(TDirectory* dir, const std::vector<std::string>& preferSubstrings) {
  if (!dir) return nullptr;

  // First pass: collect candidate histograms.
  std::vector<TH1*> cands;
  TIter next(dir->GetListOfKeys());
  while (TKey* key = (TKey*)next()) {
    TObject* obj = key->ReadObj();
    if (!obj) continue;

    if (obj->InheritsFrom(TDirectory::Class())) {
      if (TH1* h = findMassHistRecursive((TDirectory*)obj, preferSubstrings)) return h;
      continue;
    }

    if (!obj->InheritsFrom(TH1::Class())) continue;

    auto* h = (TH1*)obj;
    std::string fullName = std::string(dir->GetPath()) + "/" + h->GetName();

    // Prefer anything that looks like a dimuon mass in VBF.
    // We keep it broad so it works with slightly different naming.
    auto lower = fullName;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    bool looksLikeMass = (lower.find("dimuon") != std::string::npos ||
                          lower.find("mumu")   != std::string::npos ||
                          lower.find("mass")   != std::string::npos);

    if (looksLikeMass) cands.push_back(h);
  }

  if (cands.empty()) return nullptr;

  // Rank candidates by whether they contain preferred substrings.
  auto score = [&](TH1* h) {
    std::string n = std::string(h->GetName());
    std::string ln = n; std::transform(ln.begin(), ln.end(), ln.begin(), ::tolower);
    int s = 0;
    for (size_t i = 0; i < preferSubstrings.size(); ++i) {
      std::string p = preferSubstrings[i];
      std::string lp = p; std::transform(lp.begin(), lp.end(), lp.begin(), ::tolower);
      if (ln.find(lp) != std::string::npos) s += (int)(100 - 10*i);
    }
    // Extra preference: contains "VBF" in path.
    std::string path = std::string(dir->GetPath());
    std::string lpath = path; std::transform(lpath.begin(), lpath.end(), lpath.begin(), ::tolower);
    if (lpath.find("vbf") != std::string::npos) s += 25;
    return s;
  };

  std::sort(cands.begin(), cands.end(), [&](TH1* a, TH1* b){ return score(a) > score(b); });
  return cands.front();
}

static TH1D* makeWindowedHist(const TH1* hIn, double xmin, double xmax) {
  if (!hIn) return nullptr;

  const int b1 = std::max(1, hIn->GetXaxis()->FindBin(xmin));
  const int b2 = std::min(hIn->GetNbinsX(), hIn->GetXaxis()->FindBin(xmax));
  const int nb = b2 - b1 + 1;

  // Create a new histogram with the same binning in [xmin,xmax] (bin edges copied).
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

// Numerical quantile from CDF via bisection.
static double quantileFromCdf(RooAbsReal& cdf, RooRealVar& x, double q, double xmin, double xmax) {
  double lo = xmin;
  double hi = xmax;
  for (int it = 0; it < 80; ++it) {
    double mid = 0.5*(lo+hi);
    x.setVal(mid);
    double v = cdf.getVal();
    if (v < q) lo = mid;
    else       hi = mid;
  }
  return 0.5*(lo+hi);
}

static double fwhmFromPdf(RooAbsPdf& pdf, RooRealVar& x, double xmin, double xmax) {
  // Sample pdf on a fine grid and extract half-max crossings.
  RooArgSet normSet(x); // <-- FIX: persistent normalization set

  const int N = 5000;
  double maxv = -1;
  double xmaxAt = xmin;
  for (int i = 0; i < N; ++i) {
    double xv = xmin + (xmax-xmin) * (i + 0.5) / N;
    x.setVal(xv);
    double v = pdf.getVal(normSet);          // <-- FIXED
    if (v > maxv) { maxv = v; xmaxAt = xv; }
  }
  if (maxv <= 0) return 0.0;
  const double half = 0.5 * maxv;

  // Find left crossing
  double left = xmin;
  for (int i = 0; i < N; ++i) {
    double xv = xmaxAt - (xmaxAt-xmin) * (i + 0.5) / N;
    x.setVal(xv);
    if (pdf.getVal(normSet) < half) { left = xv; break; }   // <-- FIXED
  }

  // Find right crossing
  double right = xmax;
  for (int i = 0; i < N; ++i) {
    double xv = xmaxAt + (xmax-xmaxAt) * (i + 0.5) / N;
    x.setVal(xv);
    if (pdf.getVal(normSet) < half) { right = xv; break; }  // <-- FIXED
  }
  return std::max(0.0, right - left);
}

} // anon

// Usage:
//   root -l -q 'MassResolution_DSCB.C("/eos/.../VBF_all.root","12022")'
void MassResolution_SCR(const char* inputPath,
                         const char* era = "12022") {
  using namespace RooFit;

  if (!inputPath || std::string(inputPath).empty()) {
    std::cerr << "[ERROR] inputPath is empty";
    return;
  }

  std::unique_ptr<TFile> f(TFile::Open(inputPath, "READ"));
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] Could not open file: " << inputPath << " ";
    return;
  }

  const std::string eraStr = era ? std::string(era) : std::string("?");
  const double lumi_pb = eraLumiPb(eraStr);
  const double lumi_fb = lumi_pb / 1000.0;

  TString tmp = inputPath;
  if (tmp.EndsWith("/")) tmp.Chop();

  // remove filename: .../VBF/VBF_all.root -> .../VBF
  tmp = tmp(0, tmp.Last('/'));

  // remove sample dir: .../VBF -> .../noJERC_noJetVeto_noBtagSF_ZPtRW
  tmp = tmp(0, tmp.Last('/'));

  // now tmp ends with .../<era>/<CorrType>/<Corrections>
  TString corrections = tmp(tmp.Last('/') + 1, tmp.Length());
  tmp = tmp(0, tmp.Last('/'));

  TString corrType = tmp(tmp.Last('/') + 1, tmp.Length());

  // output directory
  TString outdir =
    "/eos/user/d/dblend/StackedPlots/" + TString(eraStr.c_str()) +
    "/Mass_Resolutions/SCR/" + corrType + "/" + corrections + "/";

  gSystem->Exec("mkdir -p " + outdir);

  auto sanitize = [](std::string s) {
    std::replace(s.begin(), s.end(), '/', '_');
    std::replace(s.begin(), s.end(), ' ', '_');
    return s;
  };

  auto fitAndPlot = [&](TH1* hIn, const std::string& histPathTag) {
    if (!hIn) {
      std::cerr << "[ERROR] Histogram not found: " << histPathTag << "/DiMuonMass_Higgs";
      return;
    }

    // Window to the Higgs region
    const double xmin = 115.0;
    const double xmax = 135.0;

    std::unique_ptr<TH1D> hWin(makeWindowedHist(hIn, xmin, xmax));
    if (!hWin) {
      std::cerr << "[ERROR] Failed to build windowed histogram for: " << histPathTag << " ";
      return;
    }

    // Optional light rebin if extremely fine
    if (hWin->GetNbinsX() > 400) hWin->Rebin(2);

    // Observable
    RooRealVar m("m", "m_{#mu#mu} [GeV]", xmin, xmax);

    // Data
    RooDataHist data("data", "data", RooArgList(m), Import(*hWin));

    // Double-sided Crystal Ball
    RooRealVar mean("mean", "mean", 125.0, 123.0, 127.0);
    RooRealVar sigma("sigma", "sigma", 1.6, 0.2, 6.0);
    RooRealVar alphaL("alphaL", "alphaL", 1.5, 0.2, 6.0);
    RooRealVar nL("nL", "nL", 5.0, 0.5, 80.0);
    RooRealVar alphaR("alphaR", "alphaR", 2.0, 0.2, 6.0);
    RooRealVar nR("nR", "nR", 8.0, 0.5, 80.0);

    RooCrystalBall dscb("dscb", "Double-sided Crystal Ball", m, mean, sigma, alphaL, nL, alphaR, nR);

    // Fit
    std::unique_ptr<RooFitResult> fr(dscb.fitTo(
      data, Save(true), PrintLevel(-1), SumW2Error(true),
      Minimizer("Minuit2","migrad"), Strategy(1)
    ));

    // Plot
    std::string title = "Higgs-region mass resolution (DSCB): " + histPathTag;
    auto* frame = m.frame(Title(title.c_str()));
    data.plotOn(frame, Name("data"));
    dscb.plotOn(frame, LineWidth(2), Name("dscb"));

    // Goodness of fit
    const int nParams = 6; // mean, sigma, alphaL, nL, alphaR, nR
    double chi2ndf = frame->chiSquare("dscb", "data", nParams);

    // Resolution metrics from fitted pdf
    std::unique_ptr<RooAbsReal> cdf(dscb.createCdf(m));
    double q16 = quantileFromCdf(*cdf, m, 0.16, xmin, xmax);
    double q84 = quantileFromCdf(*cdf, m, 0.84, xmin, xmax);
    double sigma68 = 0.5 * (q84 - q16);

    double fwhm = fwhmFromPdf(dscb, m, xmin, xmax);
    double sigmaFromFwhm = (fwhm > 0 ? fwhm / 2.354820045 : 0.0);

    // Labels
    auto* pt = new TPaveText(0.58, 0.60, 0.90, 0.88, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->AddText(("Era: " + eraStr).c_str());
    pt->AddText(Form("Lumi: %.3f fb^{-1}", lumi_fb));
    pt->AddText(histPathTag.c_str());
    pt->AddText(Form("mean = %.3f GeV", mean.getVal()));
    pt->AddText(Form("sigma (core) = %.3f GeV", sigma.getVal()));
    pt->AddText(Form("sigma_{68} (q84-q16)/2 = %.3f GeV", sigma68));
    pt->AddText(Form("FWHM = %.3f GeV  (#sigma=%.3f)", fwhm, sigmaFromFwhm));
    pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2ndf));

    auto* c = new TCanvas(("c_"+sanitize(histPathTag)).c_str(), "c", 900, 700);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    frame->GetYaxis()->SetTitleOffset(1.35);
    frame->Draw();
    pt->Draw("same");

    // Save
    const std::string tag = sanitize(histPathTag);
    const TString baseName = "MassResolution_DSCB_" + TString(tag.c_str());
    c->SaveAs(outdir + baseName + ".pdf");
    c->SaveAs(outdir + baseName + ".png");

    delete c;
  };

  // --- Your two known histogram locations:
  TH1* h_raw  = dynamic_cast<TH1*>(f->Get("DM_Pre-Selection/DiMuonMass_Higgs"));
  TH1* h_corr = dynamic_cast<TH1*>(f->Get("corr_DM_Pre-Selection/DiMuonMass_Higgs"));

  fitAndPlot(h_raw,  "DM_Pre-Selection");
  fitAndPlot(h_corr, "corr_DM_Pre-Selection");
}
