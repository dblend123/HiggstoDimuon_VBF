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

void CrystalBall_ZMassRegion() {
    using namespace RooFit;

    // Observable in tight Z window
    RooRealVar x("x", "m_{#mu#mu}", 70, 115);

    // Load histogram
    //TFile *f = new TFile("./TagAndProbe_DimuonInvMass_FullDYRun22.root", "READ");
    TFile *f = new TFile("./TagAndProbe_DimuonInvMass_Muon_FullDataRun22.root", "READ");
    //TH1* h_pass = (TH1*)f->Get("h_mass_pass_pt_3");
    //TH1* h_fail = (TH1*)f->Get("h_mass_fail_pt_3");
    TH1* h_pass = (TH1*)f->Get("h_mass_pass_eta_3");
    //TH1* h_fail = (TH1*)f->Get("h_mass_fail_eta_3");

    RooDataHist data("data_pass", "Passing Probes", RooArgList(x), h_pass);
    //RooDataHist data_fail("data_fail", "Failing Probes", RooArgList(x), h_fail);

    // ---------- MODEL 3: Crystal-Ball DY PDF ----------
    RooRealVar cb_mean("cb_mean", "CB mean", 91.2, 88.0, 94.0);
    RooRealVar cb_sigma("cb_sigma", "CB sigma", 2.0, 0.5, 5.0);      
    RooRealVar cb_alpha("cb_alpha", "CB alpha", 1.5, 0.5, 10.0);
    RooRealVar cb_n("cb_n", "CB n", 3.0, 0.5, 20.0);
    RooCBShape cb("cb", "Crystal Ball", x, cb_mean, cb_sigma, cb_alpha, cb_n);

    // Linear background: RooPolynomial of degree one
    RooRealVar bg_slope("bg_slope", "Linear slope", 0.0, -0.5, 0.5);  // init at 0 for convergence
    RooPolynomial background("background", "Linear background", x, RooArgList(bg_slope));

    // Combine: singal + background
    RooRealVar frac("frac", "Signal fraction", 0.95, 0.0, 1.0);
    RooAddPdf model("model", "CB + linear bkg", RooArgList(cb, background), RooArgList(frac));

  // ------------------------
  // Plotting
  // ------------------------

    // Full models
    auto fitResult = model.fitTo(data, Save(true), PrintLevel(-1), Minos(true));
    RooPlot* frame = x.frame(Title("Crystal Ball + Linear Background Fit"));
    data.plotOn(frame, Name("datacurve"));
    model.plotOn(frame, VisualizeError(*fitResult, 2), FillColor(kGray + 2), DrawOption("F"));
    model.plotOn(frame, VisualizeError(*fitResult, 1), FillColor(kGray), DrawOption("F"));
    model.plotOn(frame, LineColor(kBlue), Name("modelcurve"));
    frame->Print("v");

    // Chi² calculation
    RooHist* histData = frame->getHist("datacurve");
    RooCurve* curveModel = frame->getCurve("modelcurve");

    double chi2 = 0;
    int ndf = 0;
    for (int i = 0; i < histData->GetN(); ++i) {
        double xVal, y_data, y_model, err;
        histData->GetPoint(i, xVal, y_data);
        err = histData->GetErrorY(i);
        
        y_model = curveModel->Eval(xVal);
        std::cout << "Bin " << i << ": data = " << y_data << ", model = " << y_model << ", error = " << err << std::endl;
        
        if (err > 0) {
            chi2 += std::pow((y_data - y_model) / err, 2);
            ++ndf;
        }
    }

    int nParams = 6; // mean, sigma, alpha, n, slope, frac
    ndf -= nParams;
    std::cout << "Manual Chi²/NDOF = " << chi2 / ndf << std::endl;

    // Observable set
    RooArgSet observables (x);

    // RooFit chi²
    double chi2_roofit = frame->chiSquare("modelcurve", "datacurve", nParams);

    // NLL from RooDataHist
    auto nll = std::unique_ptr<RooAbsReal>(model.createNLL(data));
    
    // Print results

    std::cout << "Chi²/NDOF for Crystal Ball + Signal Fit: " << chi2_roofit << std::endl;
    std::cout << "--------------------------------------------" << std::endl;  
    std::cout << "NLL for Crystal Ball + Signal Fit: " << nll->getVal() << std::endl;

    // Legend
    auto legend = new TLegend(0.6, 0.7, 0.88, 0.88);
    legend->AddEntry(frame->findObject("modelcurve"), "CB + Linear Bkg", "l");

    auto errBand = new TBox();  // dummy for 1σ band
    errBand->SetFillColor(kGray);
    errBand->SetLineColor(kGray);
    legend->AddEntry(errBand, "1#sigma uncertainty", "f");

    // Canvas
    TCanvas c("c", "CB_ZMass", 800, 600);
    frame->SetMinimum(0);
    frame->Draw();
    legend->Draw();
    c.SaveAs("CrystalBall_ZMassRegion.pdf");
}

