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

void Gauss_Poly() {
  using namespace RooFit;

  // Observable
  RooRealVar x("x", "m_{#mu#mu}", 110, 150);
  
  // Load histogram
  TFile Run22File("V3_PartialRun22.root");
  TH1* Run22Hist = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR2");
  RooDataHist data("data", "data", RooArgList(x), Run22Hist);
/*
  // Signal: Gaussian
  RooRealVar mean("mean", "mean", 125.35, 125.2, 125.5);
  RooRealVar sigma("sigma", "sigma", 0.15, 0.05, 0.3);
  RooGaussian gauss("gauss", "signal", x, mean, sigma);
  RooRealVar nsig("nsig", "signal yield", 1000, 0, 3000);
*/
  // ---------- MODEL 1: Polynomial ----------
  RooRealVar p0("p0", "p0", 1, -1e4, 1e4);
  RooRealVar p1("p1", "p1", 0, -500, 500);
  RooRealVar p2("p2", "p2", 0, -100, 100);
  RooPolynomial poly("poly", "poly", x, RooArgList(p0, p1, p2));
  RooRealVar nbkg1("nbkg1", "background yield", 1e5, 0, 1e7);
  //RooAddPdf model_poly("model_poly", "gauss+poly", RooArgList(gauss, poly), RooArgList(nsig, nbkg1));

  //model_poly.fitTo(data, PrintLevel(-1));
  poly.fitTo(data, PrintLevel(-1));

  // Plotting
  gStyle->SetOptStat(0);
  RooPlot* frame = x.frame(Title("Gauss + Poly"));
  data.plotOn(frame);

  // Full models
  //model_poly.plotOn(frame, LineColor(kRed), LineStyle(kSolid), Name("poly"));
  poly.plotOn(frame, LineColor(kRed), LineStyle(kSolid), Name("poly"));

  // Legend
  auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
  legend->AddEntry(frame->findObject("poly"), "Gauss + Poly", "l");

  // Draw and save
  TCanvas* c = new TCanvas("c", "Gauss + Poly", 800, 600);
  frame->GetXaxis()->SetTitle("m_{#mu#mu} [GeV]");
  frame->Draw();
  legend->Draw();

  //c->SaveAs("Gauss_Poly.pdf");
  c->SaveAs("Polytest.pdf");

}
