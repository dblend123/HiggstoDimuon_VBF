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
	//RooRealVar x("x","m_{#mu#mu}", 85, 97);
	RooRealVar x("x","m_{#mu#mu}", 70, 115);
	
    // Load histogram 
    TFile *f = new TFile("./TagAndProbe_DimuonInvMass_FullDYRun22.root", "READ");
    TH1* h = (TH1*)f->Get("h_mass_fail_pt_4");
    
    RooDataHist data("data", "data", RooArgList(x), h);
	
    // ---------- MODEL 3: Voigtian DY PDF ----------
	RooRealVar Vmean("Vmean", "Z mass", 91.2, 88.0, 94.0);
	RooRealVar Vwidth("Vwidth", "Z width", 2.495);
    Vwidth.setConstant(true);
 	RooRealVar Vsigma("Vsigma", "Voigtian Sigma", 1.5, 0.1, 5.0);
	RooVoigtian Vmodel("Vmodel", "Voigtian PDF", x, Vmean, Vwidth, Vsigma); 

    // --- Polynomial Modifier (1st order) ---
    //RooRealVar c1("c1", "1st-order shape modifier", 0.0, -0.5, 0.5);
    RooRealVar c1("c1", "1st-order shape modifier", 0.0, -1.0, 1.0);
    RooPolynomial poly("poly", "Shape tweak", x, RooArgList(c1));

    // --- Modified Voigtian: Voigt × Polynomial ---
    RooProdPdf adaptedVoigt("adaptedVoigt", "Adapted signal", RooArgList(Vmodel, poly));

    // Exponential background: Drell-Yan continuum
    //RooRealVar slope("slope", "DY slope", -0.1, -5.0, 0.0);
    RooRealVar slope("slope", "DY slope", 0.0, -.5, 0.5);
    RooExponential background("background", "linear background", x, slope);

    // Combine: signal + background
    //RooRealVar frac("frac", "Signal fraction", 0.9, 0.0, 1.0);
    RooRealVar frac("frac", "Signal fraction", 0.95, 0.0, 1.0);
    RooAddPdf model("model", "Voigtian + background", RooArgList(adaptedVoigt, background), RooArgList(frac));

  // ------------------------
  // Plotting
  // ------------------------
  
    auto fitResult = model.fitTo(data, Save(true), PrintLevel(-1), Minos(true));
    
    // Plot data
    //RooPlot* frame = x.frame(Title("DY Fit: Voigtian x Poly + Background Fit"));
    RooPlot* frame = x.frame(Title("Failing Probes (13.6 TeV 2022)"));
    data.plotOn(frame, RooFit::Name("datacurve"));
    
    // Plot total model with 1σ and 2σ uncertainty bands
    //model.plotOn(frame, VisualizeError(*fitResult, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
    //model.plotOn(frame, VisualizeError(*fitResult, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
    
    // Plot background only
    //model.plotOn(frame, Components("adaptedVoigt"), LineColor(kRed), LineStyle(kDashed), Name("signalcurve"));
    //model.plotOn(frame, Components("background"), LineColor(kGreen + 2), LineStyle(kDotted), Name("backgroundcurve"));

    // Plot signal only (scaled appropriately)
    //model.plotOn(frame, Components("signal"), LineColor(kRed), LineStyle(kDashed), Name("signalcurve"));

    // Plot full model on top
    //model.plotOn(frame, LineColor(kBlue), RooFit::Name("modelcurve"));
    model.plotOn(frame, LineColor(kRed), RooFit::Name("modelcurve"));
    //frame->Print("v");

    // --- Chi² Calculation ---
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
  
  int nParams = 5;
  ndf -= nParams;
  std::cout << "Manual Chi²/NDOF = " << chi2 / ndf << std::endl;

  // Observable set
  //RooArgSet observables(x);
  //int nParams = 3; // 3 Voigtian parameters

  // Chi²/NDOF from the RooPlot
  double chi2_Voigtian = frame->chiSquare("modelcurve", "datacurve", nParams);

  // NLL from RooDataHist
  auto nll = std::unique_ptr<RooAbsReal>(model.createNLL(data));

    // Print results
    std::cout << "Chi²/NDOF for Voigtian + Background Fit:  " << chi2_Voigtian << std::endl;
    std::cout << "--------------------------------------------" << std::endl;  
    std::cout << "NLL for Voigtian + Background Fit:  " << nll->getVal() << std::endl;

    // Legend
    /*auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->AddEntry("modelcurve", "Voigt × poly + bkg", "l");
    legend->AddEntry("signalcurve", "Signal (Voigt × poly)", "l");
    legend->AddEntry("backgroundcurve", "Linear background", "l"); 
    legend->AddEntry(frame->findObject("modelcurve"), "Adapted Voigtian + Linear Background", "l");
    legend->AddEntry(frame->findObject("signalcurve"), "Adapted Voigtian", "l");
    legend->AddEntry(frame->findObject("backgroundcurve"), "Linear Background", "l");*/

    //auto errBand = new TBox();  // dummy object for legend
    TBox* errBand = new TBox();
    errBand->SetFillColor(kGray);
    errBand->SetLineColor(kGray);
    //legend->AddEntry(errBand, "1#sigma uncertainty", "f");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.92, "CMS Simulation Preliminary");

    latex.DrawLatex(0.65, 0.88, "Drell-Yan");
    latex.DrawLatex(0.65, 0.83, "40.0 GeV < p_{T} < 50.0 GeV");

    // Canvas and Save
    TCanvas c("c", "VoigtianAdapted", 800, 600);
    frame->SetMinimum(0);
    frame->Draw();
    //legend->Draw();
    c.SaveAs("lookatthis.pdf");
}
