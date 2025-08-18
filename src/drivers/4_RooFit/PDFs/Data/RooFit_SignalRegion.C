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

using namespace RooFit;

void RooFit_SignalRegion(){

	// Observable
	RooRealVar x("MuMu Mass","MuMu Mass",110,150,"GeV/c^{2}");
	x.setBins(25);
	
	// Open ROOT files and load histograms
	TFile Run22File("V3_PartialRun22.root");
	TH1* Run22Hist = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR2");
	RooDataHist data("data", "data", RooArgList(x), Run22Hist);
	
    // Signal model (Gaussian)
    RooRealVar mean("mean", "mean", 125.35, 125.2, 125.5);
    RooRealVar sigma("sigma", "sigma", .15, 0.0, .3);
    RooGaussian gauss("gauss", "gaussian signal", x, mean, sigma);

    // Background model (2nd-order polynomial)
    RooRealVar p0("p0", "p0", 1.0, -1e2, 1e4);
    RooRealVar p1("p1", "p1", 0.0, -500, 500);
    RooRealVar p2("p2", "p2", 0.0, -100, 100);
    RooPolynomial bkg("bkg", "background", x, RooArgList(p0, p1, p2));

	// Yields
    RooRealVar nsig("nsig", "signal yield", 750, 0, 1500);
    RooRealVar nbkg("nbkg", "background yield", 5e5, 1e5, 1e7);

    // Combined model
    RooAddPdf model("model", "signal + background", RooArgList(gauss, bkg), RooArgList(nsig, nbkg));

    // Fit to data
    model.fitTo(data);

    // Overlay fitted models onto another histogram
	RooPlot* frame = x.frame();
	data.plotOn(frame);
    model.plotOn(frame);
    model.plotOn(frame, Components("bkg"), LineStyle(kDashed), LineColor(kRed));
    model.plotOn(frame, Components("gauss"), LineStyle(kDotted), LineColor(kBlue));
    
	// Draw the frame with the overlaid models
	TCanvas canvas("canvas", "canvas", 800, 600);
	frame->Draw();
	canvas.SaveAs("Roofit_PartialRun_SRSB2_Canvas.pdf");

/*
	//Chi Square Calculations	
	RooChi2Var chi2("chi2","chi2", sum, hist, true);
	RooMinuit minuit(chi2);
	minuit.migrad();
	minuit.hesse();
*/
}
