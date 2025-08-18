#include "RooFormulaVar.h"
#include "TStyle.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
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

void CBG_ZMassRegion(TString Tag) {
    using namespace RooFit;

    // Observable in tight Z window
    RooRealVar x("x", "m_{#mu#mu}", 70, 115);

    // Load histogram
    //TFile *f = new TFile("./TagAndProbe_DimuonInvMass_FullDYRun22.root", "READ");
    TFile *f = new TFile("./TagAndProbe_DimuonInvMass_Muon_FullDataRun22.root", "READ");

    TString passName = "h_mass_pass_" + Tag;
    TString failName = "h_mass_fail_" + Tag;
    
    TH1* h_pass = (TH1*)f->Get(passName);
    TH1* h_fail = (TH1*)f->Get(passName);

    if (!h_pass || !h_fail) {
        std::cerr << "ERROR: Histogram not found: "
              << (!h_pass ? passName : failName) << std::endl;
    return;
    }

    TString title_pass = "Passing_" + Tag + "_Probes";
    TString title_fail = "Failing_" + Tag + "_Probes";
    
    RooDataHist data_pass("data_pass", title_pass.Data(), RooArgList(x), h_pass);
    RooDataHist data_fail("data_fail", title_fail.Data(), RooArgList(x), h_fail);

    // ---------- Signal Shape (e.g. Crystal Ball) ----------
RooRealVar cb_mean("cb_mean", "CB mean", 91.2, 88.0, 94.0);
RooRealVar cb_sigma("cb_sigma", "CB sigma", 1.8, 0.5, 5.0);
RooRealVar cb_alpha("cb_alpha", "CB alpha", 1.2, 0.5, 5.0);
RooRealVar cb_n("cb_n", "CB n", 3.0, 0.1, 10.0);
RooCBShape cb("cb", "Crystal Ball", x, cb_mean, cb_sigma, cb_alpha, cb_n);

RooRealVar cb_sigma2("sigma2", "Width 2", 2.5, 0.5, 5);
RooRealVar cb_alpha2("alpha2", "Alpha 2", -1.2, -5, -0.5);
RooRealVar cb_n2("n2", "n2", 3.0, 0.1, 10);
RooCBShape cb2("cb2", "CB high", x, cb_mean, cb_sigma2, cb_alpha2, cb_n2);

    // ---------- Signal Shape (e.g. Gaussian) ----------
RooRealVar gaus_mean("gaus_mean", "Gauss mean", 91.2, "GeV");  // can tie to cb_mean
RooRealVar gaus_sigma("gaus_sigma", "Gauss sigma", 2.0, 0.5, 5.0);
RooGaussian gaus("gaus", "Gaussian", x, gaus_mean, gaus_sigma);

    // ---------- Fraction of CB in the signal ----------
RooRealVar frac1("frac1", "CB1 fraction", 0.5, 0.0, 1.0);
RooRealVar frac2("frac2", "CB2 fraction", 0.3, 0.0, 1.0);
RooAddPdf sig("sig", "CB1 + CB2 + Gaus", RooArgList(cb, cb2, gaus), RooArgList(frac1, frac2));

  // --------------------- Background Model ---------------------
RooRealVar exp_slope("exp_slope", "Exp Slope", -0.1, -5.0, 0.0);
RooExponential background("background", "exp bkg", x, exp_slope);

// ---------- Background shapes ----------
//RooRealVar slopeP("slopeP", "bkg slope (pass)", 0.0, -1.0, 1.0);
//RooPolynomial bkg_pass("bkg_pass", "background pass", x, RooArgList(slopeP));

//RooRealVar slopeF("slopeF", "bkg slope (fail)", 0.0, -1.0, 1.0);
//RooPolynomial bkg_fail("bkg_fail", "background fail", x, RooArgList(slopeF));

// ---------- Yields ----------
RooRealVar eff("eff", "Efficiency", 0.8, 0.0, 1.0);
RooRealVar nSigTotal("nSigTotal", "Total Signal", 1e5, 0, 1e7);

RooFormulaVar nSigP("nSigP", "@0*@1", RooArgList(eff, nSigTotal));
RooRealVar nBkgP("nBkgP", "Bkg yield (pass)", 1000, 0, 1e6);
RooAddPdf model_pass("model_pass", "sig + bkg (pass)", RooArgList(sig, background), RooArgList(nSigP, nBkgP));

RooFormulaVar nSigF("nSigF", "(1-@0)*@1", RooArgList(eff, nSigTotal));
RooRealVar nBkgF("nBkgF", "Bkg yield (fail)", 1000, 0, 1e6);
RooAddPdf model_fail("model_fail", "sig + bkg (fail)", RooArgList(sig, background), RooArgList(nSigF, nBkgF));

// ---------- Fits ----------
auto fr_pass = model_pass.fitTo(data_pass, Save(true), Extended(true), PrintLevel(-1));
auto fr_fail = model_fail.fitTo(data_fail, Save(true), Extended(true), PrintLevel(-1));

// ---------- Efficiency Calculation ----------
double eff_val = eff.getVal();
double eff_err = eff.getError();

double sP = eff_val * nSigTotal.getVal();
double sF = (1.0 - eff_val) * nSigTotal.getVal();

double eP = eff_err * nSigTotal.getVal();
double eF = eff_err * nSigTotal.getVal();

double eff_calc = sP / (sP + sF);
double eff_calc_err = sqrt(pow(eP / (sP + sF), 2) + pow(sP * eF / pow(sP + sF, 2), 2));

std::cout << "--------------------------------------------------" << std::endl;
std::cout << "nSigP = " << sP << " ± " << eP << std::endl;
std::cout << "nSigF = " << sF << " ± " << eF << std::endl;
std::cout << "Tag-and-Probe Efficiency = " << eff_calc << " ± " << eff_err << std::endl;
std::cout << "--------------------------------------------------" << std::endl;

  // ------------------------
  // Plotting
  // ------------------------

    // Full models
    RooPlot* frameP = x.frame(Title("Passing Probes (13.6 TeV, 2022)"));
    data_pass.plotOn(frameP, Name("pass_datacurve"));
    //model.plotOn(frame, VisualizeError(*fitResult, 2), FillColor(kGray + 2), DrawOption("F"));
    //model.plotOn(frame, VisualizeError(*fitResult, 1), FillColor(kGray), DrawOption("F"));
    model_pass.plotOn(frameP, LineColor(kGreen + 2), Name("pass_modelcurve"));
    //model_pass.plotOn(frameF, VisualizeError(*fr_fail, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
    //model_pass.plotOn(frameF, VisualizeError(*fr_fail, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));

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
    TCanvas cpass("cpass", "pass_CBG", 800, 600);
    frameP->SetMinimum(0);
    frameP->Draw();
    //legend->Draw();
    DrawCMSLabel();
    //cpass.SaveAs("Passing_" + Tag + "_Probe_CBG_ZMassRegion.png");

    TCanvas cfail("cfail", "fail_CBG", 800, 600);
    frameF->SetMinimum(0);
    frameF->Draw();
    //legend->Draw();
    DrawCMSLabel();
    //cfail.SaveAs("Failing_" + Tag + "_Probe_CBG_ZMassRegion.png");
}
