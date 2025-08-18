#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooChebychev.h"
#include "RooGenericPdf.h"

void FitWithMultipleModels() {
  using namespace RooFit;

  // Observable
  RooRealVar x("x", "m_{#mu#mu}", 110, 150);
  
  // Load histogram
  TFile Run22File("V3_PartialRun22.root");
  TH1* Run22Hist = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR2");
  RooDataHist data("data", "data", RooArgList(x), Run22Hist);

  // Signal: Gaussian
  RooRealVar mean("mean", "mean", 125.35, 125.2, 125.5);
  RooRealVar sigma("sigma", "sigma", 0.15, 0.05, 0.3);
  RooGaussian gauss("gauss", "signal", x, mean, sigma);
  RooRealVar nsig("nsig", "signal yield", 1000, 0, 3000);

  // ---------- MODEL 1: Polynomial ----------
  RooRealVar p0("p0", "p0", 1, -1e4, 1e4);
  RooRealVar p1("p1", "p1", 0, -500, 500);
  RooRealVar p2("p2", "p2", 0, -100, 100);
  RooPolynomial poly("poly", "poly", x, RooArgList(p0, p1, p2));
  RooRealVar nbkg1("nbkg1", "background yield", 1e5, 0, 1e7);
  RooAddPdf model_poly("model_poly", "gauss+poly", RooArgList(gauss, poly), RooArgList(nsig, nbkg1));

  // ---------- MODEL 2: Exponential ----------
  RooRealVar tau("tau", "tau", -0.1);  // Just set initial value
  tau.setRange(-5, 0);                 // Set limits separately
  RooExponential expo("expo", "expo", x, tau);
  RooRealVar nbkg2("nbkg2", "background yield", 1e5, 0, 1e7);
  RooAddPdf model_expo("model_expo", "gauss+expo", RooArgList(gauss, expo), RooArgList(nsig, nbkg2));

  // ---------- MODEL 3: Chebychev Polynomial ----------
  RooRealVar c1("c1", "c1", 0.1, -1, 1);
  RooRealVar c2("c2", "c2", 0.1, -1, 1);
  RooChebychev cheb("cheb", "chebyshev", x, RooArgList(c1, c2));
  RooRealVar nbkg3("nbkg3", "background yield", 1e5, 0, 1e7);
  RooAddPdf model_cheb("model_cheb", "gauss+cheb", RooArgList(gauss, cheb), RooArgList(nsig, nbkg3));

  // ------------------------
  // Fit all three models
  // ------------------------
  model_poly.fitTo(data, PrintLevel(-1));
  model_expo.fitTo(data, PrintLevel(-1));
  model_cheb.fitTo(data, PrintLevel(-1));

  // ------------------------
  // Plotting
  // ------------------------
  gStyle->SetOptStat(0);
  RooPlot* frame = x.frame(Title("Fit Comparison"));
  data.plotOn(frame);

  // Full models
  model_poly.plotOn(frame, LineColor(kRed), LineStyle(kSolid), Name("poly"));
  model_expo.plotOn(frame, LineColor(kGreen+2), LineStyle(kDashed), Name("expo"));
  model_cheb.plotOn(frame, LineColor(kBlue), LineStyle(kDotted), Name("cheb"));

  // Legend
  auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
  legend->AddEntry(frame->findObject("poly"), "Gauss + Poly", "l");
  legend->AddEntry(frame->findObject("expo"), "Gauss + Expo", "l");
  legend->AddEntry(frame->findObject("cheb"), "Gauss + Cheb", "l");

  // Draw and save
  TCanvas* c = new TCanvas("c", "Fit Comparison", 800, 600);
  frame->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  frame->Draw();
  legend->Draw();

  c->SaveAs("FitWithMultipleModels.pdf");

}
