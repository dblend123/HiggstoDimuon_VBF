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
#include <TLatex.h>

void DrawCMSLabel() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.82, "#font[62]{CMS}");
    latex.DrawLatex(0.20, 0.78, "#it{Simulation}");
    latex.DrawLatex(0.20, 0.74, "#it{Preliminary}");
    latex.DrawLatex(0.62, 0.82, "#font[62]{Signal MC 2022}");
}

void Voigtian_ZMassRegion(TString Tag){
    using namespace RooFit;

	//Observable
	RooRealVar x("x","m_{#mu#mu} (GeV)", 70, 115);
	
    // Load histogram 
    //TFile *f = new TFile("./TagAndProbe_DimuonInvMass_FullDYRun22.root", "READ");
    TFile *f = new TFile("./TagAndProbe_DimuonInvMass_Muon_FullDataRun22.root", "READ");
    
    TString passName = "h_mass_pass_" + Tag;
    TString failName = "h_mass_fail_" + Tag;
    
    TH1* h_pass = (TH1*)f->Get(passName);
    TH1* h_fail = (TH1*)f->Get(failName);
    
    if (!h_pass || !h_fail) {
        std::cerr << "ERROR: Histogram not found: "
              << (!h_pass ? passName : failName) << std::endl;
    return;
    }
    
    TString title_pass = "Passing_" + Tag + "_Probes";
    TString title_fail = "Failing_" + Tag + "_Probes";

    RooDataHist data_pass("data_pass", title_pass.Data(), RooArgList(x), h_pass);
    RooDataHist data_fail("data_fail", title_fail.Data(), RooArgList(x), h_fail);
	
    // ---------- MODEL 3: Voigtian DY PDF ----------
	RooRealVar Vmean("Vmean", "Z mass", 91.2, 88.0, 94.0);
	RooRealVar Vwidth("Vwidth", "Z width", 2.495);
    Vwidth.setConstant(true);
 	RooRealVar Vsigma("Vsigma", "Voigtian Sigma", 1.5, 0.1, 5.0);
	RooVoigtian Vmodel("Vmodel", "Voigtian PDF", x, Vmean, Vwidth, Vsigma); 

    // --- Polynomial Modifier (1st order)  w/ Yields ---
    RooRealVar slopeP("slopeP", "bkg slope (pass)", 0.0, -1.0, 1.0);
    RooPolynomial bkg_pass("bkg_pass", "background pass", x, RooArgList(slopeP));

    RooRealVar slopeF("slopeF", "bkg slope (fail)", 0.0, -1.0, 1.0);
    RooPolynomial bkg_fail("bkg_fail", "background fail", x, RooArgList(slopeF));
   
    RooRealVar nSigP("nSigP", "Signal yield (pass)", 5000, 0, 1e6);
    RooRealVar nBkgP("nBkgP", "Bkg yield (pass)", 1000, 0, 1e6);
    RooAddPdf model_pass("model_pass", "sig + bkg (pass)", RooArgList(Vmodel, bkg_pass), RooArgList(nSigP, nBkgP));

    RooRealVar nSigF("nSigF", "Signal yield (fail)", 3000, 0, 1e6);
    RooRealVar nBkgF("nBkgF", "Bkg yield (fail)", 1000, 0, 1e6);
    RooAddPdf model_fail("model_fail", "sig + bkg (fail)", RooArgList(Vmodel, bkg_fail), RooArgList(nSigF, nBkgF));
/*
    // Exponential background: Drell-Yan continuum
    //RooRealVar slope("slope", "DY slope", -0.1, -5.0, 0.0);
    RooRealVar slope("slope", "DY slope", 0.0, -.5, 0.5);
    RooExponential background("background", "linear background", x, slope);

    // Combine: signal + background
    //RooRealVar frac("frac", "Signal fraction", 0.9, 0.0, 1.0);
    RooRealVar frac("frac", "Signal fraction", 0.95, 0.0, 1.0);
    RooAddPdf model("model", "Voigtian + background", RooArgList(adaptedVoigt, background), RooArgList(frac));
*/
  // ------------------------
  // Fit both categories 
  // ------------------------
  
    auto fr_pass = model_pass.fitTo(data_pass, Save(true), Extended(true), PrintLevel(-1));
    auto fr_fail = model_fail.fitTo(data_fail, Save(true), Extended(true), PrintLevel(-1));

    // ------------------------------
    // Efficiency calculation
    // ------------------------------
    
    double sP = nSigP.getVal();
    double sF = nSigF.getVal();
    double eP = nSigP.getError();
    double eF = nSigF.getError();

    double eff = sP / (sP + sF);
    double eff_err = sqrt(std::pow(eP / (sP + sF), 2) + std::pow(sP * eF / std::pow(sP + sF, 2), 2));

    std::cout << "--------------------------------------------------" << std::endl;
    std::cout << "nSigP = " << sP << " ± " << eP << std::endl;
    std::cout << "nSigF = " << sF << " ± " << eF << std::endl;
    std::cout << "Tag-and-Probe Efficiency = " << eff << " ± " << eff_err << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    // ------------------------------
    // Plotting
    // ------------------------------
    
    RooPlot* frameP = x.frame(Title("Passing Probes (13.6 TeV, 2022)"));
    data_pass.plotOn(frameP, RooFit::Name("pass_datacurve"));
    //model_pass.plotOn(frameP, Components("Vmodel"), LineColor(kRed), LineStyle(kDashed));
    //model_pass.plotOn(frameP, Components("bkg_pass"), LineColor(kGreen + 2), LineStyle(kDotted));
    model_pass.plotOn(frameP, LineColor(kGreen + 2), RooFit::Name("pass_modelcurve"));
    //model_pass.plotOn(frameP, VisualizeError(*fr_pass, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
    //model_pass.plotOn(frameP, VisualizeError(*fr_pass, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));

    RooPlot* frameF = x.frame(Title("Failing Probes (13.6 TeV, 2022)"));
    data_fail.plotOn(frameF, RooFit::Name("fail_datacurve"));
    //model_fail.plotOn(frameF, Components("Vmodel"), LineColor(kRed), LineStyle(kDashed));
    //model_fail.plotOn(frameF, Components("bkg_fail"), LineColor(kGreen + 2), LineStyle(kDotted));
    model_fail.plotOn(frameF, LineColor(kRed), RooFit::Name("fail_modelcurve"));
    //model_fail.plotOn(frameF, VisualizeError(*fr_fail, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
    //model_fail.plotOn(frameF, VisualizeError(*fr_fail, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
    
    // ------------------------------
    // Chi² Calculation
    // ------------------------------
    RooHist* pass_histData = frameP->getHist("pass_datacurve");
    RooCurve* pass_curveModel = frameP->getCurve("pass_modelcurve");

    double pass_chi2 = 0;
    int pass_ndf = 0;
    for (int i = 0; i < pass_histData->GetN(); ++i) {
        double x, y_data, y_model, err;
        pass_histData->GetPoint(i, x, y_data);
        err = pass_histData->GetErrorY(i);
        y_model = pass_curveModel->Eval(x);
        std::cout << "Bin " << i << ": data = " << y_data << ", model = " << y_model << ", error = " << err << std::endl;
        if (err > 0) {
            pass_chi2 += std::pow((y_data - y_model) / err, 2);
            ++pass_ndf;
        }
    }
 
    RooHist* fail_histData = frameF->getHist("fail_datacurve");
    RooCurve* fail_curveModel = frameF->getCurve("fail_modelcurve");

    double fail_chi2 = 0;
    int fail_ndf = 0;
    for (int i = 0; i < fail_histData->GetN(); ++i) {
        double x, y_data, y_model, err;
        fail_histData->GetPoint(i, x, y_data);
        err = fail_histData->GetErrorY(i);
        y_model = fail_curveModel->Eval(x);
        std::cout << "Bin " << i << ": data = " << y_data << ", model = " << y_model << ", error = " << err << std::endl;
        if (err > 0) {
            fail_chi2 += std::pow((y_data - y_model) / err, 2);
            ++fail_ndf;
        }
    }

  int nParams = 5;
  pass_ndf -= nParams;
  fail_ndf -= nParams;
  std::cout << "Manual Chi²/NDOF = " << pass_chi2 / pass_ndf << std::endl;
  std::cout << "Manual Chi²/NDOF = " << fail_chi2 / fail_ndf << std::endl;

  // Chi²/NDOF from the RooPlot
  double chi2_passVoigtian = frameP->chiSquare("pass_modelcurve", "pass_datacurve", nParams);
  double chi2_failVoigtian = frameF->chiSquare("fail_modelcurve", "fail_datacurve", nParams);

  // NLL from RooDataHist
  auto pass_nll = std::unique_ptr<RooAbsReal>(model_pass.createNLL(data_pass));
  auto fail_nll = std::unique_ptr<RooAbsReal>(model_fail.createNLL(data_fail));

    // Print results
    std::cout << "Chi²/NDOF for Passing Voigtian + Background Fit:  " << chi2_passVoigtian << std::endl;
    std::cout << "Chi²/NDOF for Failing Voigtian + Background Fit:  " << chi2_failVoigtian << std::endl;
    std::cout << "--------------------------------------------" << std::endl;  
    std::cout << "NLL for Voigtian + Background Fit:  " << pass_nll->getVal() << std::endl;
    std::cout << "NLL for Voigtian + Background Fit:  " << fail_nll->getVal() << std::endl;

    //auto errBand = new TBox();  // dummy object for legend
    TBox* errBand = new TBox();
    errBand->SetFillColor(kGray);
    errBand->SetLineColor(kGray);
    //legend->AddEntry(errBand, "1#sigma uncertainty", "f");

    // Canvas and Save
    TCanvas cpass("cpass", "pass_VoigtianAdapted", 800, 600);
    frameP->SetMinimum(0);
    frameP->Draw();
    //legend->Draw();
    DrawCMSLabel();
    cpass.SaveAs("Passing_" + Tag + "_Probe_Voigitian_ZMassRegion.png");

    TCanvas cfail("cfail", "fail_VoigtianAdapted", 800, 600);
    frameF->SetMinimum(0);
    frameF->Draw();
    //legend->Draw();
    DrawCMSLabel();
    
    TLatex latex;
    //latex.DrawLatex(0.55, 0.71, "40.0 GeV < muon p_{T} < 50.0 GeV");
    latex.DrawLatex(0.65, 0.71, "-2.4 <  #eta  < -1.6");
    cfail.SaveAs("Failing_" + Tag + "_Probe_Voigitian_ZMassRegion.png");
}
