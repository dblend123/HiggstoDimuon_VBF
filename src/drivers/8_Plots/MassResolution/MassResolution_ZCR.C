// MassResolution_ZCR_Voigt.C
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

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooVoigtian.h"

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

  // ---- helper that fits & saves
  auto fitAndPlot = [&](TH1* hIn, const std::string& histPathTag) {
    if (!hIn) {
      std::cerr << "[WARN] Missing hist for tag '" << histPathTag << "'\n";
      return;
    }

    std::unique_ptr<TH1D> hWin(makeWindowedHist(hIn, xmin, xmax));
    if (!hWin) {
      std::cerr << "[ERROR] Failed windowing for tag '" << histPathTag << "'\n";
      return;
    }

    // Optional light rebin if extremely fine
    if (hWin->GetNbinsX() > 600) hWin->Rebin(2);

    const double xlo = hWin->GetXaxis()->GetXmin();
    const double xhi = hWin->GetXaxis()->GetXmax();

    // (optional) print once for sanity
    std::cout << "[ZCR] Using RooFit range [" << xlo << "," << xhi << "]\n";

    RooRealVar m("m", "m_{#mu#mu} [GeV]", xlo, xhi);
    RooDataHist data("data", "data", RooArgList(m), Import(*hWin));

    // --- Voigtian: BW (Gamma) ⊗ Gaussian (sigma)
    // mean around Z pole
    RooRealVar mean("mean", "mean", 91.1876, 90.5, 91.9);

    // Natural Z width ~2.4952 GeV.
    // If you want to *float* it, widen the range and remove setConstant(true).
    RooRealVar gamma("gamma", "gamma", 2.4952);
    gamma.setConstant(true);

    // Detector resolution proxy
    RooRealVar sigma("sigma", "sigma", 1.3, 0.2, 5.0);

    RooVoigtian voigt("voigt", "Z lineshape (Voigtian)", m, mean, gamma, sigma);

    std::unique_ptr<RooFitResult> fr(voigt.fitTo(
      data, Save(true), PrintLevel(-1), SumW2Error(true),
      Minimizer("Minuit2","migrad"), Strategy(1)
    ));

    auto* frame = m.frame(Title(("Z control region fit (Voigtian): " + histPathTag).c_str()));
    data.plotOn(frame, Name("data"));
    voigt.plotOn(frame, LineWidth(2), Name("voigt"));

    // chi2/ndf (nParams = mean + sigma; gamma fixed)
    const int nParams = 2;
    const double chi2ndf = frame->chiSquare("voigt", "data", nParams);

    // sigma68 and FWHM from fitted pdf
    std::unique_ptr<RooAbsReal> cdf(voigt.createCdf(m));
    const double q16 = quantileFromCdf(*cdf, m, 0.16, xlo, xhi);
    const double q84 = quantileFromCdf(*cdf, m, 0.84, xlo, xhi);
    const double fwhm = fwhmFromPdf(voigt, m, xlo, xhi);
    
    const double sigma68 = 0.5 * (q84 - q16);

    const double sigmaFromFwhm = (fwhm > 0 ? fwhm / 2.354820045 : 0.0);

    auto* pt = new TPaveText(0.58, 0.60, 0.90, 0.88, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->AddText(("Era: " + eraStr).c_str());
    pt->AddText(Form("Lumi: %.3f fb^{-1}", lumi_fb));
    pt->AddText(histPathTag.c_str());
    pt->AddText(Form("mean = %.4f GeV", mean.getVal()));
    pt->AddText(Form("gamma (fixed) = %.4f GeV", gamma.getVal()));
    pt->AddText(Form("sigma (Gauss) = %.4f GeV", sigma.getVal()));
    pt->AddText(Form("sigma_{68} (q84-q16)/2 = %.4f GeV", sigma68));
    pt->AddText(Form("FWHM = %.4f GeV  (#sigma=%.4f)", fwhm, sigmaFromFwhm));
    pt->AddText(Form("#chi^{2}/ndf = %.3f", chi2ndf));

    auto* c = new TCanvas(("c_"+sanitize(histPathTag)).c_str(), "c", 900, 700);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    frame->GetYaxis()->SetTitleOffset(1.35);
    frame->Draw();
    pt->Draw("same");

    const std::string tag = sanitize(histPathTag);
    const TString baseName = "MassResolution_ZCR_" + TString(tag.c_str());
    c->SaveAs(outdir + baseName + ".pdf");
    c->SaveAs(outdir + baseName + ".png");

    delete c;
  };

// ---- Build summed background histogram from all files

const char* histName = "corr_DM_Pre-Selection/DiMuonMass_Z";

std::unique_ptr<TH1> h_sum;  // total background

for (const auto& p : paths) {

    std::unique_ptr<TFile> f(TFile::Open(p.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "[WARN] Skipping unreadable file: " << p << "\n";
        continue;
    }

    TH1* h = dynamic_cast<TH1*>(f->Get(histName));
    if (!h) {
        std::cerr << "[WARN] Missing " << histName
                  << " in: " << p << "\n";
        continue;
    }

    // Clone so histogram survives after file closes
    TH1* hc = dynamic_cast<TH1*>(h->Clone());
    hc->SetDirectory(nullptr);

    if (!h_sum) {
        h_sum.reset(dynamic_cast<TH1*>(hc->Clone("h_bkg_sum")));
        h_sum->SetDirectory(nullptr);
    } else {
        h_sum->Add(hc);
    }
}

// Safety check
if (!h_sum) {
    std::cerr << "[ERROR] No histograms were loaded.\n";
    return;
}

// Fit the TOTAL stacked background
fitAndPlot(h_sum.get(), "BKG_SUM_corr_DM_Pre-Selection_DiMuonMass_Z");
}
