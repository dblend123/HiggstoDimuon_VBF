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

void Roofit_ZMassRegion(){
    using namespace RooFit;

	//Observable
	RooRealVar x("x","m_{#mu#mu}", 110, 150);
	
    // Load histogram 
	TFile *f = new TFile("./output_Signal_FullRun22.root","READ");
	
    TH1* Run22Hist2 = (TH1*)f->Get("VBF/DiMuonMass_SB_SR2");
    TH1* Run22Hist25 = (TH1*)f->Get("VBF/DiMuonMass_SB_SR2.5");
    TH1* Run22Hist3 = (TH1*)f->Get("VBF/DiMuonMass_SB_SR3");
    TH1* Run22Hist35 = (TH1*)f->Get("VBF/DiMuonMass_SB_SR3.5");
    
    RooDataHist data2("data2", "data2", RooArgList(x), Run22Hist2);
    RooDataHist data25("data25", "data25", RooArgList(x), Run22Hist25);
    RooDataHist data3("data3", "data3", RooArgList(x), Run22Hist3);
    RooDataHist data35("data35", "data35", RooArgList(x), Run22Hist35);
	
/*	//Build a Gaussian PDF
	RooRealVar sigmean("sigmean", "Z mass", 91.2, 75.0, 105.0);
	RooRealVar sigwidth("sigwidth", "Z width", 2.5, .1, 10.); //Still not sure about these values 
	RooGaussian signalModel("Gaussian","Gaussian PDF", x, sigmean, sigwidth);	

	//Build Crystal-Ball Signal PDF
 	RooRealVar cbMean("cbMean", "CrystalBall Mean", 5.291, 5.0, 5.5);
 	RooRealVar cbSigma("cbSigma", "CrystalBall Sigma", 0.1, 0.01, 1.0);
	RooRealVar cbAlpha("cbAlpha", "CrystalBall Alpha", -1.5, -10.0, 0.0);
 	RooRealVar cbN("cbN", "CrystalBall N", 2.0, 0.1, 10.0);
 	RooCBShape signalModel("crystalBall", "Crystal Ball PDF", x, cbMean, cbSigma, cbAlpha, cbN);
*/
    // ---------- MODEL 3: Voigtian Signal PDF ----------
	RooRealVar Vmean_2("Vmean_2", "Z mass", 125.35, 110.0, 150.0);
	RooRealVar Vwidth_2("Vwidth_2", "Z width", 2.5, .1, 10.); 
 	RooRealVar Vsigma_2("Vsigma_2", "Voigtian Sigma", 2.5, 0.1, 10.0);
	RooVoigtian Vmodel_2("Vmodel_2", "Voigtian PDF", x, Vmean_2, Vwidth_2, Vsigma_2); 

	RooRealVar Vmean_25("Vmean_25", "Z mass", 125.35, 110.0, 150.0);
	RooRealVar Vwidth_25("Vwidth_25", "Z width", 2.5, .1, 10.); 
 	RooRealVar Vsigma_25("Vsigma_25", "Voigtian Sigma", 2.5, 0.1, 10.0);
	RooVoigtian Vmodel_25("Vmodel_25", "Voigtian PDF", x, Vmean_25, Vwidth_25, Vsigma_25); 

	RooRealVar Vmean_3("Vmean_3", "Z mass", 125.35, 110.0, 150.0);
	RooRealVar Vwidth_3("Vwidth_3", "Z width", 2.5, .1, 10.); 
 	RooRealVar Vsigma_3("Vsigma_3", "Voigtian Sigma", 2.5, 0.1, 10.0);
	RooVoigtian Vmodel_3("Vmodel_3", "Voigtian PDF", x, Vmean_3, Vwidth_3, Vsigma_3); 

	RooRealVar Vmean_35("Vmean_35", "Z mass", 125.35, 110.0, 150.0);
	RooRealVar Vwidth_35("Vwidth_35", "Z width", 2.5, .1, 10.); 
 	RooRealVar Vsigma_35("Vsigma_35", "Voigtian Sigma", 2.5, 0.1, 10.0);
	RooVoigtian Vmodel_35("Vmodel_35", "Voigtian PDF", x, Vmean_35, Vwidth_35, Vsigma_35); 

/*	//Build Breit-Wigner background PDF 
	RooRealVar bwpar("bwpar","breit-wigner parameter", -20.0, -100.0, -1.);
	RooBreitWigner background("background", "Breit-Wigner PDF", x, RooConst(91.2), bwpar);	
*/
  
  // ------------------------
  // Fit Models
  // ------------------------
  Vmodel_2.fitTo(data2, PrintLevel(-1));
  Vmodel_25.fitTo(data25, PrintLevel(-1));
  Vmodel_3.fitTo(data3, PrintLevel(-1));
  Vmodel_35.fitTo(data35, PrintLevel(-1));

  // ------------------------
  // Plotting
  // ------------------------
  RooPlot* frame2 = x.frame(Title("Vmodel_2 RooFit"));
  RooPlot* frame25 = x.frame(Title("Vmodel_25 RooFit"));
  RooPlot* frame3 = x.frame(Title("Vmodel_3 RooFit"));
  RooPlot* frame35 = x.frame(Title("Vmodel_35 RooFit"));
  
  data2.plotOn(frame2, RooFit::Name("data2curve"), RooFit::Binning(20));
  data25.plotOn(frame25, RooFit::Name("data25curve"));
  data3.plotOn(frame3, RooFit::Name("data3curve"));
  data35.plotOn(frame35, RooFit::Name("data35curve"));

  // Full models
  auto fitResult2 = Vmodel_2.fitTo(data2, Save(true), PrintLevel(-1), Minos(true));
  Vmodel_2.plotOn(frame2, VisualizeError(*fitResult2, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  Vmodel_2.plotOn(frame2, VisualizeError(*fitResult2, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  Vmodel_2.plotOn(frame2, LineColor(kBlue), RooFit::Name("Vmodel2curve"), RooFit::Binning(20));
  frame2->Print("v");
  auto fitResult25 = Vmodel_25.fitTo(data25, Save(true), PrintLevel(-1), Minos(true));
  Vmodel_25.plotOn(frame25, VisualizeError(*fitResult25, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  Vmodel_25.plotOn(frame25, VisualizeError(*fitResult25, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  Vmodel_25.plotOn(frame25, LineColor(kBlue), RooFit::Name("Vmodel25curve"));
  auto fitResult3 = Vmodel_3.fitTo(data3, Save(true), PrintLevel(-1), Minos(true));

  if (!fitResult3 || fitResult3->covQual() < 2) {
  std::cerr << "fitResult3 has no valid covariance matrix (covQual = "
            << (fitResult3 ? fitResult3->covQual() : -1) << ")." << std::endl;
  } else {
  std::cout << "fitResult3 is OK. covQual = " << fitResult3->covQual() << std::endl;
  fitResult3->Print("v");  // verbose print
  }

  Vmodel_3.plotOn(frame3, VisualizeError(*fitResult3, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  Vmodel_3.plotOn(frame3, VisualizeError(*fitResult3, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  Vmodel_3.plotOn(frame3, LineColor(kBlue), RooFit:: Name("Vmodel3curve"));
  auto fitResult35 =Vmodel_35.fitTo(data35, Save(true), PrintLevel(-1), Minos(true));
  Vmodel_35.plotOn(frame35, VisualizeError(*fitResult35, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  Vmodel_35.plotOn(frame35, VisualizeError(*fitResult35, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  Vmodel_35.plotOn(frame35, LineColor(kBlue), RooFit::Name("Vmodel35curve"));

  RooHist* histData = frame2->getHist("data2curve");
  RooCurve* curveModel = frame2->getCurve("Vmodel2curve");

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
  
  int nParams = 3;
  ndf -= nParams;
  std::cout << "Manual Chi²/NDOF = " << chi2 / ndf << std::endl;

  // Observable set
  RooArgSet observables(x);
  //int nParams = 3; // 3 Voigtian parameters

  // Chi²/NDOF from the RooPlot
  double chi2_Voigtian2 = frame2->chiSquare("Vmodel2curve", "data2curve", nParams);
  double chi2_Voigtian25 = frame25->chiSquare("Vmodel25curve", "data25curve", nParams);
  double chi2_Voigtian3 = frame3->chiSquare("Vmodel3curve", "data3curve", nParams);
  double chi2_Voigtian35 = frame35->chiSquare("Vmodel35curve", "data35curve", nParams);

  // NLL from RooDataHist
  auto nll2 = std::unique_ptr<RooAbsReal>(Vmodel_2.createNLL(data2));
  auto nll25 = std::unique_ptr<RooAbsReal>(Vmodel_25.createNLL(data25));
  auto nll3 = std::unique_ptr<RooAbsReal>(Vmodel_3.createNLL(data3));
  auto nll35 = std::unique_ptr<RooAbsReal>(Vmodel_35.createNLL(data35));

  // Print results
  std::cout << "Chi²/NDOF for Voigtian2 Fit:  " << chi2_Voigtian2 << std::endl;
  std::cout << "Bins in TH1: " << Run22Hist2->GetNbinsX() << std::endl;
  std::cout << "Chi²/NDOF for Voigtian25 Fit:  " << chi2_Voigtian25 << std::endl;
  std::cout << "Chi²/NDOF for Voigtian3 Fit:  " << chi2_Voigtian3 << std::endl;
  std::cout << "Chi²/NDOF for Voigtian35 Fit:  " << chi2_Voigtian35 << std::endl;
  std::cout << "--------------------------------------------" << std::endl;  
  std::cout << "NLL for Voigtian2 Fit:  " << nll2->getVal() << std::endl;
  std::cout << "NLL for Voigtian25 Fit:  " << nll25->getVal() << std::endl;
  std::cout << "NLL for Voigtian3 Fit:  " << nll3->getVal() << std::endl;
  std::cout << "NLL for Voigtian35 Fit:  " << nll35->getVal() << std::endl;

  // Legend
  auto legend2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  auto legend25 = new TLegend(0.6, 0.7, 0.88, 0.88);
  auto legend3 = new TLegend(0.6, 0.7, 0.88, 0.88);
  auto legend35 = new TLegend(0.6, 0.7, 0.88, 0.88);
  
  legend2->AddEntry(frame2->findObject("Vmodel2curve"), "Voigtian", "l");
  legend25->AddEntry(frame25->findObject("Vmodel25curve"), "Voigtian", "l");
  legend3->AddEntry(frame3->findObject("Vmodel3curve"), "Voigtian", "l");
  legend35->AddEntry(frame35->findObject("Vmodel35curve"), "Voigtian", "l");

  auto errBand2 = new TBox();  // dummy object for legend
  errBand2->SetFillColor(kGray);
  errBand2->SetLineColor(kGray);
  legend2->AddEntry(errBand2, "1#sigma uncertainty", "f");

  auto errBand25 = new TBox();  // dummy object for legend
  errBand25->SetFillColor(kGray);
  errBand25->SetLineColor(kGray);
  legend25->AddEntry(errBand25, "1#sigma uncertainty", "f");

  auto errBand3 = new TBox();  // dummy object for legend
  errBand3->SetFillColor(kGray);
  errBand3->SetLineColor(kGray);
  legend3->AddEntry(errBand3, "1#sigma uncertainty", "f");

  auto errBand35 = new TBox();  // dummy object for legend
  errBand35->SetFillColor(kGray);
  errBand35->SetLineColor(kGray);
  legend35->AddEntry(errBand35, "1#sigma uncertainty", "f");

  TString base = "Signal_FullRun22";

  // Canvas and Save
  TCanvas c2("c2", "SR2", 800, 600);
  frame2->SetMinimum(0);
  frame2->Draw();
  legend2->Draw();
  c2.SaveAs(base +  "_SR2.pdf");

  TCanvas c25("c25", "SR2.5", 800, 600);
  frame25->SetMinimum(0);
  frame25->Draw();
  legend25->Draw();
  c25.SaveAs(base + "_SR2_5.pdf");

  TCanvas c3("c3", "SR3", 800, 600);
  frame3->SetMinimum(0);
  frame3->SetMaximum(frame3->GetMaximum() * 1.2);
  frame3->Draw();
  legend3->Draw();
  c3.Update();
  c3.SaveAs(base + "_SR3.pdf");

  TCanvas c35("c35", "SR3.5", 800, 600);
  frame35->SetMinimum(0);
  legend35->Draw();
  frame35->Draw();
  c35.SaveAs(base + "_SR3_5.pdf");
}
