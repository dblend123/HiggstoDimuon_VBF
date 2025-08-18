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
#include "RooCBShape.h"
#include "RooVoigtian.h"
#include "TLegend.h" 
#include "RooHist.h"

void Voigtian_ZMassRegion(){
    using namespace RooFit;

	//Observable
	RooRealVar x("x","m_{#mu#mu}", 85, 97);
	
    // Load histogram 
	TFile *f = new TFile("./test.root","READ");
	
    TH1* Run22Hist = (TH1*)f->Get("NoCats/DiMuonMass_Z");
    
    RooDataHist data("data", "data", RooArgList(x), Run22Hist);
	
    // ---------- MODEL 3: Voigtian DY PDF ----------
	RooRealVar Vmean("Vmean", "Z mass", 91.2, 88.0, 94.0);
	RooRealVar Vwidth("Vwidth", "Z width", 2.495);
    Vwidth.setConstant(true);
 	RooRealVar Vsigma("Vsigma", "Voigtian Sigma", 1.5, 0.1, 5.0);
	RooVoigtian Vmodel("Vmodel", "Voigtian PDF", x, Vmean, Vwidth, Vsigma); 

    // Exponential background: Drell-Yan continuum
    RooRealVar slope("slope", "DY slope", -0.1, -5.0, 0.0);
    RooExponential background("background", "DY background", x, slope);

    // Combine: signal + background
    RooRealVar frac("frac", "Signal fraction", 0.9, 0.0, 1.0);
    RooAddPdf model("model", "Voigtian + background", RooArgList(Vmodel, background), RooArgList(frac));

  // ------------------------
  // Plotting
  // ------------------------
  
    // Full models
    auto fitResult = model.fitTo(data, Save(true), PrintLevel(-1), Minos(true));
    RooPlot* frame = x.frame(Title("DY Fit with Voigtian + Background"));
    data.plotOn(frame, RooFit::Name("datacurve"));
    model.plotOn(frame, VisualizeError(*fitResult, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
    model.plotOn(frame, VisualizeError(*fitResult, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
    model.plotOn(frame, LineColor(kBlue), RooFit::Name("modelcurve"));
    frame->Print("v");

  RooHist* histData = frame->getHist("datacurve");
  RooCurve* curveModel = frame->getCurve("modelcurve");

  double chi2 = 0;
  int ndf = 0;
  for (int i = 0; i < histData->GetN(); ++i) {
      double x, y_data, y_model, err;
      histData->GetPoint(i, x, y_data);
      err = histData->GetErrorY(i);

      y_model = curveModel->Eval(x);
      std::cout << "Bin " << i << ": data = " << y_data << ", model = " << y_model << ", error = " << err << std::endl;

      if (err > 0) {
          chi2 += std::pow((y_data - y_model) / err, 2);
          ++ndf;
      }
  }
  
  int nParams = 4;
  ndf -= nParams;
  std::cout << "Manual Chi²/NDOF = " << chi2 / ndf << std::endl;

  // Observable set
  RooArgSet observables(x);
  //int nParams = 3; // 3 Voigtian parameters

  // Chi²/NDOF from the RooPlot
  double chi2_Voigtian = frame->chiSquare("modelcurve", "datacurve", nParams);

  // NLL from RooDataHist
  auto nll = std::unique_ptr<RooAbsReal>(model.createNLL(data));

  // Print results
  std::cout << "Chi²/NDOF for Voigtian Fit:  " << chi2_Voigtian << std::endl;
  std::cout << "--------------------------------------------" << std::endl;  
  std::cout << "NLL for Voigtian2 Fit:  " << nll->getVal() << std::endl;

  // Legend
  auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
  
  legend->AddEntry(frame->findObject("modelcurve"), "Voigtian + background", "l");

  auto errBand = new TBox();  // dummy object for legend
  errBand->SetFillColor(kGray);
  errBand->SetLineColor(kGray);
  legend->AddEntry(errBand, "1#sigma uncertainty", "f");

  // Canvas and Save
  TCanvas c("c", "Voigtian_Z", 800, 600);
  frame->SetMinimum(0);
  frame->Draw();
  legend->Draw();
  c.SaveAs("Voigtian_ZMassRegion.pdf");
}
