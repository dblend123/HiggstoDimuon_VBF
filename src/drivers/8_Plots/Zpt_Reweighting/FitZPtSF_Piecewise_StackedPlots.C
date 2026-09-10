// ======================= FitZPtSF_Piecewise_StackedPlots.C =======================
// Option A: read SF_DataOverDY from existing ZPtSF output root file and fit piecewise polynomials.
// Also creates a final combined TF1 per jet bin (piecewise wrapper).
//
// Ranges (fixed):
//   low  : 0   - 60   GeV
//   mid  : 60  - 120  GeV
//   high : 120 - 200  GeV
//
// Outputs are auto-routed to:
//   /eos/user/d/dblend/StackedPlots/<YEAR>/ZPt_Corrections/<MuCorr>/<CfgTag>/Fits/
//
// Usage:
//   root -l -q 'FitZPtSF_Piecewise_StackedPlots.C(".../ZPtSF_*.root")'
// Optional (polynomial orders):
//   root -l -q 'FitZPtSF_Piecewise_StackedPlots.C(".../ZPtSF_*.root", 2, 2, 1)'
//
// Adds:
//  1) chi2/ndof for each segment + overall combined piecewise
//  2) Pull distributions (per segment + full range)
//  3) F-test comparing polN vs pol(N+1) per segment
//
// ==============================================================================

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TMath.h"
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "TMatrixDSym.h"
#include <memory>

// ------------------------- Building Polynomial in u = log(x + shift), order N --------------------------
// f(x) = p0 + p1*u + p2*u^2 + ... + pN*u^N
static TF1* makePolLogTF1(const std::string& name, int order,
                          double xmin, double xmax,
                          double shift = 1.0)
{
  std::string u = Form("log(x+%.6f)", shift);

  std::string expr = Form("([0])");
  for (int i = 1; i <= order; ++i) {
    expr += Form(" + ([%d])*pow(%s,%d)", i, u.c_str(), i);
  }

  TF1* f = new TF1(name.c_str(), expr.c_str(), xmin, xmax);
  f->SetNpx(600);
  f->SetParameter(0, 1.0);
  for (int i = 1; i <= order; ++i) f->SetParameter(i, 0.0);
  return f;
}

// ------------------------- path utils --------------------------
static std::vector<std::string> splitPath(const std::string& p) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : p) {
    if (c == '/') {
      if (!cur.empty()) out.push_back(cur);
      cur.clear();
    } else cur.push_back(c);
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

static bool deriveTagsFromAnyPath(const std::string& path,
                                 std::string& year,
                                 std::string& muCorr,
                                 std::string& cfg)
{
  auto parts = splitPath(path);

  for (int i = 0; i < (int)parts.size(); ++i) {
    if (parts[i] == "StackedPlots") {
      if (i + 5 < (int)parts.size() && parts[i+2] == "ZPt_Corrections") {
        year   = parts[i+1];
        muCorr = parts[i+3];
        cfg    = parts[i+4];
        return true;
      }
    }
  }

  for (int i = 0; i < (int)parts.size(); ++i) {
    if (parts[i] == "Analyzed_SkimmedFiles") {
      if (i + 3 < (int)parts.size()) {
        year   = parts[i+1];
        muCorr = parts[i+2];
        cfg    = parts[i+3];
        return true;
      }
    }
  }

  return false;
}

// ------------------------- hist helper --------------------------
static TH1* getHist(TFile* f, const std::string& path) {
  if (!f) return nullptr;
  auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
  if (!h) {
    std::cerr << "[ERROR] Missing hist: " << path << " in " << f->GetName() << "\n";
    return nullptr;
  }
  auto* c = (TH1*)h->Clone();
  c->SetDirectory(nullptr);
  if (!c->GetSumw2N()) c->Sumw2();
  return c;
}

static void sanitizeSF(TH1* h, double maxAbs = 5.0) {
  if (!h) return;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    const double y = h->GetBinContent(i);
    const double e = h->GetBinError(i);
    if (!std::isfinite(y) || !std::isfinite(e) || y <= 0.0 || y > maxAbs) {
      h->SetBinContent(i, 0.0);
      h->SetBinError(i, 0.0);
    }
  }
}

static void applyRelErrorFloor(TH1* h,
                               double xmin, double xmax,
                               double relFloor = 0.05) // 5% default
{
  if (!h) return;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double x = h->GetBinCenter(i);
    if (x < xmin || x > xmax) continue;

    double y = h->GetBinContent(i);
    double e = h->GetBinError(i);
    if (y == 0 || e <= 0) continue;

    double eFloor = relFloor * std::abs(y);
    double eNew = std::sqrt(e*e + eFloor*eFloor);
    h->SetBinError(i, eNew);
  }
}

static int countUsableBins(const TH1* h, double xmin, double xmax) {
  int n=0;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < xmin || x > xmax) continue;
    if (h->GetBinError(i) > 0 && h->GetBinContent(i) != 0) ++n;
  }
  return n;
}

// ------------------------- stats helpers --------------------------
struct FitInfo {
  TF1* f = nullptr;

  int orderUsed = -1;
  int npar = 0;

  double chi2 = 0.0;
  int ndf = 0;
  double pval = 1.0;
  double ftest_pval = -1.0;

  std::unique_ptr<TMatrixDSym> cov;

  bool useLogX = false;
};

static double fitSigma(
    const FitInfo& fit,
    double x
)
{
    if (!fit.f || !fit.cov)
        return 0.0;

    if (fit.npar <= 0)
        return 0.0;

    if (fit.cov->GetNrows() != fit.npar)
        return 0.0;

    const double base =
        fit.useLogX
            ? std::log(x + 1.0)
            : x;

    std::vector<double> grad(
        fit.npar,
        1.0
    );

    for (int i = 1; i < fit.npar; ++i) {
        grad[i] =
            std::pow(
                base,
                i
            );
    }

    double variance = 0.0;

    for (int i = 0; i < fit.npar; ++i) {
        for (int j = 0; j < fit.npar; ++j) {

            variance +=
                grad[i] *
                (*fit.cov)(i, j) *
                grad[j];
        }
    }

    if (!std::isfinite(variance))
        return 0.0;

    // Tiny negative values can occur from floating-point precision.
    variance =
        std::max(
            0.0,
            variance
        );

    return std::sqrt(variance);
}

static double chi2PValue(double chi2, int ndf) {
  if (ndf <= 0) return 1.0;
  // upper-tail prob: Prob(Chi2 > chi2)
  return TMath::Prob(chi2, ndf);
}

// F-test p-value for nested models: model1 (p1 params, chi1, ndf1), model2 (p2 params, chi2, ndf2)
// Assumes model2 has more params and chi2 <= chi1 typically.
static double fTestPValue(double chi1, int ndf1, int p1,
                          double chi2, int ndf2, int p2)
{
  const int d1 = p2 - p1;
  const int d2 = ndf2;
  if (d1 <= 0 || d2 <= 0) return -1.0;

  const double num = (chi1 - chi2) / (double)d1;
  const double den = (chi2) / (double)d2;
  if (den <= 0) return -1.0;

  const double F = num / den;
  if (!std::isfinite(F) || F < 0) return -1.0;

  // upper-tail p-value for F distribution:
  // p = Prob(F_{d1,d2} > F) = 1 - CDF(F)
  const double cdf = TMath::FDistI(F, d1, d2); // integral from 0..F
  return 1.0 - cdf;
}

static TH1D* makePullHist(const TH1* h, TF1* f, double xmin, double xmax,
                          const std::string& name, int nbins=60, double pullMin=-6, double pullMax=6)
{
  if (!h || !f) return nullptr;
  TH1D* hp = new TH1D(name.c_str(), "", nbins, pullMin, pullMax);
  hp->Sumw2();
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < xmin || x > xmax) continue;
    const double y = h->GetBinContent(i);
    const double e = h->GetBinError(i);
    if (e <= 0 || y == 0) continue;
    const double pull = (y - f->Eval(x)) / e;
    if (std::isfinite(pull)) hp->Fill(pull);
  }
  return hp;
}

// Fit a segment with order, capture fit stats.
// Also optionally compute F-test vs (order+1) if enough bins.
static FitInfo fitSegmentWithStats(
    const TH1* h,
    const std::string& baseName,
    const std::string& category,
    double xmin,
    double xmax,
    int orderRequested,
    bool doFtest=true,
    bool useLogX=false) {

  FitInfo out;
  if (!h) return out;

  out.useLogX = useLogX;

  int usable = countUsableBins(h, xmin, xmax);
  int order = orderRequested;

  if (usable < order + 2) {
    std::cerr << "[WARN] Not enough usable bins for " << baseName
              << " in [" << xmin << "," << xmax << "] (usable=" << usable
              << ", order=" << order << "). Falling back to pol0.\n";
    order = 0;
  }

  out.orderUsed = order;
  out.npar = order + 1;

  if (useLogX) {
    out.f = makePolLogTF1(baseName + "_polLog" + std::to_string(order), order, xmin, xmax, 1.0);
  } else {
    const std::string form = "pol" + std::to_string(order);
    out.f = new TF1((baseName + "_" + form).c_str(), form.c_str(), xmin, xmax);
    out.f->SetParameter(0, 1.0);
    for (int i=1;i<=order;i++) out.f->SetParameter(i, 0.0);
  }

  TFitResultPtr r = const_cast<TH1*>(h)->Fit(out.f, "RSQ0");
  if ((int)r == 0 && r.Get()) { // ROOT uses status==0 for OK
    
    out.cov  = std::make_unique<TMatrixDSym>(r->GetCovarianceMatrix());

    out.chi2 = r->Chi2();
    out.ndf  = r->Ndf();
    out.pval = chi2PValue(
        out.chi2,
        out.ndf);
        
  } else {
    std::cerr << "[WARN] Fit failed for " << baseName << " order=" << order << "\n";
  }

  // F-test against order+1 (nested) if possible
  if (doFtest) {
    const int order2 = order + 1;
    if (usable >= order2 + 2) {

      TF1* f2 = nullptr;
      if (useLogX) {
        f2 = makePolLogTF1(baseName + "_polLog" + std::to_string(order2),
                           order2, xmin, xmax, 1.0);
      } else {
        f2 = new TF1((baseName + "_pol" + std::to_string(order2)).c_str(),
                     ("pol" + std::to_string(order2)).c_str(), xmin, xmax);
        f2->SetParameter(0, 1.0);
        for (int i=1;i<=order2;i++) f2->SetParameter(i, 0.0);
      }

      TFitResultPtr r2 = const_cast<TH1*>(h)->Fit(f2, "RSQ0");
      if ((int)r2 == 0 && r2.Get()) {
        const double chi1 = out.chi2;
        const int    ndf1 = out.ndf;
        const int    p1   = out.npar;

        const double chi2 = r2->Chi2();
        const int    ndf2 = r2->Ndf();
        const int    p2   = order2 + 1;

        out.ftest_pval = fTestPValue(chi1, ndf1, p1, chi2, ndf2, p2);

        // Keep the better (lower chi2) model as the "segment fit"
        const bool preferHigherOrder =
            std::isfinite(out.ftest_pval) &&
            out.ftest_pval >= 0.0 &&
            out.ftest_pval < 0.05;

        if (preferHigherOrder) {

          delete out.f;

          out.f = f2;

          out.orderUsed = order2;
          out.npar = p2;

          out.chi2 = chi2;
          out.ndf  = ndf2;

          out.pval =
              chi2PValue(
                  out.chi2,
                  out.ndf
              );

          // IMPORTANT:
          // Replace the old covariance with the covariance
          // belonging to the newly selected polynomial.
        out.cov =
            std::make_unique<TMatrixDSym>(
                r2->GetCovarianceMatrix()
        );
        } else {

          delete f2;
        }

        } else {
            delete f2;
        }
    }
  }

  // ------------------------------------------------
  // Print FINAL selected fit and its covariance
  // ------------------------------------------------
if (out.f && out.cov) {

  std::cout
      << "\n[FINAL FIT COVARIANCE] "
      << category
      << " / "
      << baseName
      << "\n";

  std::cout
      << "Selected order = pol"
      << out.orderUsed
      << "\n";

  std::cout << "Parameters:\n";

  for (int i = 0; i < out.npar; ++i) {

    std::cout
        << "  p" << i
        << " = "
        << out.f->GetParameter(i)
        << " +/- "
        << std::sqrt(
            std::max(
                0.0,
                (*out.cov)(i, i)
            )
        )
        << "\n";
  }

  std::cout
      << "Covariance matrix:\n";

  out.cov->Print();
}  
  return out;
}

// ------------------------- smooth combined TF1 --------------------------
static TF1* g_low  = nullptr;
static TF1* g_mid  = nullptr;
static TF1* g_high = nullptr;

// Smooth step from 0->1 centered at x0 with width w
static double smoothStep(double x, double x0, double w) {
  if (w <= 0) return (x > x0) ? 1.0 : 0.0;
  return 0.5 * (1.0 + TMath::TanH((x - x0)/w));
}

static double combinedSigma(
    double x,
    const FitInfo& iL,
    const FitInfo& iM,
    const FitInfo& iH
)
{
    const double xcl =
        std::min(
            200.0,
            std::max(
                5.0,
                x
            )
        );

    const double w1 = 8.0;
    const double w2 = 10.0;

    const double s1 =
        smoothStep(
            xcl,
            60.0,
            w1
        );

    const double s2 =
        smoothStep(
            xcl,
            120.0,
            w2
        );

    const double sigL =
        fitSigma(
            iL,
            xcl
        );

    const double sigM =
        fitSigma(
            iM,
            xcl
        );

    const double sigH =
        fitSigma(
            iH,
            xcl
        );

    // Same coefficients used by the nominal smooth blend.
    const double aL =
        (1.0 - s2) *
        (1.0 - s1);

    const double aM =
        (1.0 - s2) *
        s1;

    const double aH =
        s2;

    // Fits are independent because they come from separate
    // pT regions, so combine their uncertainties in quadrature.
    const double variance =
        aL * aL * sigL * sigL +
        aM * aM * sigM * sigM +
        aH * aH * sigH * sigH;

    return std::sqrt(
        std::max(
            0.0,
            variance
        )
    );
}

// Smooth blend: low -> mid around x=60, mid -> high around x=120
double smoothCombinedEval(double* x, double* /*p*/) {
  const double xx = x[0];
  if (!g_low || !g_mid || !g_high) return 1.0;

  // Clamp evaluation range (prevents weird extrapolation)
  const double xcl = std::min(200.0, std::max(5.0, xx));

  // --- widths control how “smooth” the transition is
  // With 4 GeV/bin, a width of ~6–10 GeV is a good start.
  const double w1 = 8.0;   // around 60
  const double w2 = 10.0;  // around 120

  // weights
  const double s1 = smoothStep(xcl, 60.0,  w1);   // 0 (below 60) -> 1 (above 60)
  const double s2 = smoothStep(xcl, 120.0, w2);   // 0 (below 120) -> 1 (above 120)

  // blend low+mid first, then blend that with high
  const double lowmid = (1.0 - s1)*g_low->Eval(xcl) + s1*g_mid->Eval(xcl);
  const double all    = (1.0 - s2)*lowmid          + s2*g_high->Eval(xcl);

  return all;
}

static TF1* makeCombined(const std::string& name, TF1* fL, TF1* fM, TF1* fH) {
  g_low = fL; g_mid = fM; g_high = fH;
  TF1* f = new TF1(name.c_str(), smoothCombinedEval, 0.0, 200.0, 0);
  f->SetNpx(800);
  return f;
}

// compute combined chi2/ndof by summing pulls^2 over [xmin,xmax]
static void computeChi2NdfFromHist(const TH1* h, TF1* f, double xmin, double xmax,
                                  int nparEffective,
                                  double& chi2, int& ndf, int& npts)
{
  chi2 = 0.0; npts = 0;
  if (!h || !f) { ndf = 0; return; }

  for (int i=1; i<=h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    if (x < xmin || x > xmax) continue;
    const double y = h->GetBinContent(i);
    const double e = h->GetBinError(i);
    if (e <= 0 || y == 0) continue;
    const double pull = (y - f->Eval(x)) / e;
    if (!std::isfinite(pull)) continue;
    chi2 += pull*pull;
    npts++;
  }
  ndf = npts - nparEffective;
}

// ------------------------- drawing --------------------------
static TCanvas* drawWithPulls(const std::string& title,
                              TH1* h,
                              TF1* fL, TF1* fM, TF1* fH, TF1* fC,
                              const FitInfo& iL, const FitInfo& iM, const FitInfo& iH,
                              TH1D* pullAll)
{
  if (!h) return nullptr;

  TCanvas* c = new TCanvas(("c_"+title).c_str(), "", 900, 900);

  // two pads: top (data+fit), bottom (pulls)
  TPad* p1 = new TPad(("p1_"+title).c_str(), "", 0, 0.30, 1, 1.0);
  TPad* p2 = new TPad(("p2_"+title).c_str(), "", 0, 0.00, 1, 0.30);
  p1->SetBottomMargin(0.02);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.30);
  p1->Draw(); p2->Draw();

  // --- top pad
  p1->cd();
  h->GetXaxis()->SetRangeUser(5, 200);
  h->GetYaxis()->SetTitle("SF = (Data - non-DY)/DY");
  h->SetMarkerStyle(20);
  h->Draw("E1");

  if (fC) {
    fC->SetLineWidth(3);
    fC->SetLineStyle(1);
    fC->Draw("SAME");
  }

  TLine l1(60, 0, 60, h->GetMaximum()*1.05);
  TLine l2(120,0,120,h->GetMaximum()*1.05);
  l1.SetLineStyle(2); l2.SetLineStyle(2);
  l1.Draw("SAME"); l2.Draw("SAME");

  TLegend leg(0.55, 0.58, 0.90, 0.88);
  leg.SetBorderSize(0); leg.SetFillStyle(0);
  leg.AddEntry(h, "SF hist", "lep");
  leg.AddEntry(fC, "Smooth piecewise fit", "l");
  leg.Draw();

  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.035);
  latex.DrawLatex(0.12, 0.93, title.c_str());

  // show F-test p-values (if computed)
  latex.SetTextSize(0.030);
  if (iL.ftest_pval >= 0) latex.DrawLatex(0.12, 0.85, Form("F-test low (N vs N+1): p = %.3g", iL.ftest_pval));
  if (iM.ftest_pval >= 0) latex.DrawLatex(0.12, 0.81, Form("F-test mid (N vs N+1): p = %.3g", iM.ftest_pval));
  if (iH.ftest_pval >= 0) latex.DrawLatex(0.12, 0.77, Form("F-test high (N vs N+1): p = %.3g", iH.ftest_pval));

  // --- bottom pad: pull distribution histogram (full range)
  p2->cd();
  if (pullAll) {
    pullAll->GetXaxis()->SetTitle("Pull = (data - fit)/#sigma");
    pullAll->GetYaxis()->SetTitle("Entries");
    pullAll->GetXaxis()->SetTitleSize(0.10);
    pullAll->GetYaxis()->SetTitleSize(0.10);
    pullAll->GetXaxis()->SetLabelSize(0.09);
    pullAll->GetYaxis()->SetLabelSize(0.09);
    pullAll->GetYaxis()->SetTitleOffset(0.50);
    pullAll->SetLineWidth(2);
    pullAll->Draw("HIST");
    TLine z(0, 0, 0, pullAll->GetMaximum()*1.05);
    z.SetLineStyle(2);
    z.Draw("SAME");
  }

  c->Update();
  return c;
}

void FitZPtSF_Piecewise_StackedPlots(const char* inSFFile,
                                     int orderLow = 2,
                                     int orderMid = 2,
                                     int orderHigh = 1)
{
  const double xLowMin=5,   xLowMax=60;
  const double xMidMin=60,  xMidMax=120;
  const double xHighMin=120,xHighMax=200;

  std::string year, muCorr, cfg;
  std::string inPath(inSFFile);

  if (!deriveTagsFromAnyPath(inPath, year, muCorr, cfg)) {
    std::cerr << "[WARN] Could not parse tags from input path. Output will go to ./Fits_Unknown/\n";
    year="UnknownYear"; muCorr="UnknownMuCorr"; cfg="UnknownCfg";
  }

  std::string outDir = "/eos/user/d/dblend/StackedPlots/" + year +
                       "/ZPt_Corrections/" + muCorr + "/" + cfg + "/Fits";
  gSystem->Exec(("mkdir -p " + outDir).c_str());

  std::string outRoot = outDir + "/ZPtSF_Fits_" + year + ".root";

  TFile* fin = TFile::Open(inSFFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "[FATAL] Could not open input SF file: " << inSFFile << "\n";
    return;
  }

  TFile* fout = TFile::Open(outRoot.c_str(), "RECREATE");

  const std::vector<std::string> cats = {
    "corr_DM_Pre-Selection_Z0j",
    "corr_DM_Pre-Selection_Z1j",
    "corr_DM_Pre-Selection_Z2j"
  };

  for (const auto& cat : cats) {
    const std::string hpath = cat + "/SF_DataMinusNonDY_OverDY";
    TH1* h = getHist(fin, hpath);
    if (!h) continue;

    sanitizeSF(h, 5.0);
    //applyRelErrorFloor(h, xLowMin, xLowMax, 0.05);  // 5% floor in low region
    // Use the propagated statistical uncertainties from
    // (Data - nonDY) / DY for the nominal goodness-of-fit test.
    // Do not add an arbitrary error floor here.
    // applyRelErrorFloor(h, xLowMin, xLowMax, 0.05);


    // --- Fit segments w/ stats (+ F-test vs N+1)
    FitInfo iL = fitSegmentWithStats(
        h, "f_low", cat,
        xLowMin, xLowMax,
        orderLow, true, true
    );

    FitInfo iM = fitSegmentWithStats(
        h, "f_mid", cat,
        xMidMin, xMidMax,
        orderMid, true, false
    );

    FitInfo iH = fitSegmentWithStats(
        h, "f_high", cat,
        xHighMin, xHighMax,
        orderHigh, true, false
    );

    TF1* fL = iL.f;
    TF1* fM = iM.f;
    TF1* fH = iH.f;

    // Final combined TF1 (piecewise function)
    TF1* fC = makeCombined("f_combined", fL, fM, fH);
    // ------------------------------------------------------------
    // Goodness of fit of the actual combined function
    // ------------------------------------------------------------
    const int nparEffective =
        iL.npar +
        iM.npar +
        iH.npar;

    double combinedChi2 = 0.0;
    int combinedNdf = 0;
    int combinedNpts = 0;

    computeChi2NdfFromHist(
        h,
        fC,
        0.0,
        200.0,
        nparEffective,
        combinedChi2,
        combinedNdf,
        combinedNpts
    );

    const double combinedPValue =
        combinedNdf > 0
            ? TMath::Prob(combinedChi2, combinedNdf)
            : 0.0;

// ------------------------------------------------------------
// Use the measured, background-subtracted binned SF directly.
// This preserves the variable-width 10/20 GeV bins.
// ------------------------------------------------------------
TH1* hZPtNom =
    (TH1*)h->Clone("zptWeightNom");

TH1* hZPtUp =
    (TH1*)h->Clone("zptWeightUp");

TH1* hZPtDown =
    (TH1*)h->Clone("zptWeightDown");

hZPtNom->SetDirectory(nullptr);
hZPtUp->SetDirectory(nullptr);
hZPtDown->SetDirectory(nullptr);

for (int ibin = 1;
     ibin <= hZPtNom->GetNbinsX();
     ++ibin)
{
    const double nominal =
        h->GetBinContent(ibin);

    const double error =
        h->GetBinError(ibin);

    const double up =
        nominal + error;

    const double down =
        std::max(
            0.0,
            nominal - error
        );

    hZPtNom->SetBinContent(ibin, nominal);
    hZPtNom->SetBinError(ibin, error);

    hZPtUp->SetBinContent(ibin, up);
    hZPtUp->SetBinError(ibin, 0.0);

    hZPtDown->SetBinContent(ibin, down);
    hZPtDown->SetBinError(ibin, 0.0);
}

// Give the underflow and overflow bins safe frozen values.
const int firstBin = 1;
const int lastBin = hZPtNom->GetNbinsX();

hZPtNom->SetBinContent(
    0,
    hZPtNom->GetBinContent(firstBin)
);
hZPtNom->SetBinContent(
    lastBin + 1,
    hZPtNom->GetBinContent(lastBin)
);

hZPtUp->SetBinContent(
    0,
    hZPtUp->GetBinContent(firstBin)
);
hZPtUp->SetBinContent(
    lastBin + 1,
    hZPtUp->GetBinContent(lastBin)
);

hZPtDown->SetBinContent(
    0,
    hZPtDown->GetBinContent(firstBin)
);
hZPtDown->SetBinContent(
    lastBin + 1,
    hZPtDown->GetBinContent(lastBin)
);

    // --- Pull distributions
    TH1D* pullLow  = makePullHist(h, fL, xLowMin,  xLowMax,  "pull_low");
    TH1D* pullMid  = makePullHist(h, fM, xMidMin,  xMidMax,  "pull_mid");
    TH1D* pullHigh = makePullHist(h, fH, xHighMin, xHighMax, "pull_high");
    TH1D* pullAll  = makePullHist(h, fC, 0.0, 200.0, "pull_all");

    // --- Main diagnostic plot with pull histogram (full range) shown below
    TCanvas* c = drawWithPulls(cat, h, fL, fM, fH, fC, iL, iM, iH, pullAll);
    if (c) {
      c->SaveAs((outDir + "/FitDiag_" + cat + ".png").c_str());
      c->SaveAs((outDir + "/FitDiag_" + cat + ".pdf").c_str());
    }

    // --- Save pull-only canvases too (optional but handy)
    auto savePullCanvas = [&](TH1D* hp, const std::string& tag){
      if (!hp) return (TCanvas*)nullptr;
      TCanvas* cp = new TCanvas(("c_pull_"+tag+"_"+cat).c_str(), "", 700, 600);
      hp->SetLineWidth(2);
      hp->Draw("HIST");
      TLine z(0,0,0,hp->GetMaximum()*1.05);
      z.SetLineStyle(2);
      z.Draw("SAME");
      cp->SaveAs((outDir + "/Pull_" + tag + "_" + cat + ".png").c_str());
      cp->SaveAs((outDir + "/Pull_" + tag + "_" + cat + ".pdf").c_str());
      return cp;
    };

    TCanvas* cpl = savePullCanvas(pullLow,  "low");
    TCanvas* cpm = savePullCanvas(pullMid,  "mid");
    TCanvas* cph = savePullCanvas(pullHigh, "high");
    TCanvas* cpa = savePullCanvas(pullAll,  "all");

    // --- Write outputs into ROOT file
    fout->mkdir(cat.c_str());
    fout->cd(cat.c_str());
    h->Write("SF_DataMinusNonDY_OverDY");

    if (fL) fL->Write("f_low");
    if (fM) fM->Write("f_mid");
    if (fH) fH->Write("f_high");
    if (fC) fC->Write("f_combined");

    hZPtNom->Write();
    hZPtUp->Write();
    hZPtDown->Write();

    // store pull hists
    if (pullLow)  pullLow->Write("pull_low");
    if (pullMid)  pullMid->Write("pull_mid");
    if (pullHigh) pullHigh->Write("pull_high");
    if (pullAll)  pullAll->Write("pull_all");

    // store canvases
    if (c)   c->Write("fit_canvas");
    if (cpl) cpl->Write("pull_canvas_low");
    if (cpm) cpm->Write("pull_canvas_mid");
    if (cph) cph->Write("pull_canvas_high");
    if (cpa) cpa->Write("pull_canvas_all");

    // print summary to stdout (quick scan in batch logs)
    std::cout << "\n=== " << cat << " ===\n";
    std::cout << " low : pol"  << iL.orderUsed << "  chi2/ndf=" << iL.chi2 << "/" << iL.ndf
              << "  p=" << iL.pval << "  Ftest p=" << iL.ftest_pval << "\n";
    //std::cout << "[INFO] Applied low-region rel error floor = 5%\n";
    std::cout << " mid : pol"  << iM.orderUsed << "  chi2/ndf=" << iM.chi2 << "/" << iM.ndf
              << "  p=" << iM.pval << "  Ftest p=" << iM.ftest_pval << "\n";
    std::cout << " high: pol"  << iH.orderUsed << "  chi2/ndf=" << iH.chi2 << "/" << iH.ndf
              << "  p=" << iH.pval << "  Ftest p=" << iH.ftest_pval << "\n";
    std::cout
        << " combined: chi2=" << combinedChi2
        << " npts=" << combinedNpts
        << " npar=" << nparEffective
        << " ndf=" << combinedNdf
        << " chi2/ndf="
        << (combinedNdf > 0
                ? combinedChi2 / combinedNdf
                : 0.0)
        << " p=" << combinedPValue
        << "\n";

std::cout << "\n[COMBINED FIT PULLS] " << cat << "\n";

for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double x = h->GetBinCenter(i);
    const double y = h->GetBinContent(i);
    const double e = h->GetBinError(i);

    if (e <= 0.0 || y == 0.0)
        continue;

    const double prediction = fC->Eval(x);
    const double pull = (y - prediction) / e;

    std::cout
        << " bin=" << i
        << " range=["
        << h->GetXaxis()->GetBinLowEdge(i)
        << ","
        << h->GetXaxis()->GetBinUpEdge(i)
        << "]"
        << " SF=" << y
        << " fit=" << prediction
        << " error=" << e
        << " pull=" << pull
        << "\n";
}

    // cleanup
    delete c;
    delete cpl; delete cpm; delete cph; delete cpa;
    delete pullLow; delete pullMid; delete pullHigh; delete pullAll;
    delete h;
    delete fL; delete fM; delete fH; delete fC;
    delete hZPtNom;
    delete hZPtUp;
    delete hZPtDown;
  }

  fout->Close();
  fin->Close();

  std::cout << "\n[OK] Wrote fit ROOT file:\n  " << outRoot << "\n";
  std::cout << "[OK] Wrote diagnostic plots under:\n  " << outDir << "\n";
}

