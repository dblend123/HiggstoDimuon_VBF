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
#include "RooHist.h"

void GaussCheb_SignalRegion() {
  using namespace RooFit;

  // Observable
  RooRealVar x("x", "m_{#mu#mu}", 110, 150);
  
  // Load histogram
  TFile Run22File("output_Signal_Run22EE.root");
  
  TH1* Run22Hist2 = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR2");
  TH1* Run22Hist25 = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR2.5");
  TH1* Run22Hist3 = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR3");
  TH1* Run22Hist35 = (TH1*)Run22File.Get("VBF/DiMuonMass_SB_SR3.5");
  
  RooDataHist data2("data2", "data2", RooArgList(x), Run22Hist2);
  RooDataHist data25("data25", "data25", RooArgList(x), Run22Hist25);
  RooDataHist data3("data3", "data3", RooArgList(x), Run22Hist3);
  RooDataHist data35("data35", "data35", RooArgList(x), Run22Hist35);

/*  // Signal: Gaussian
  RooRealVar mean("mean", "mean", 125.35, 125.2, 125.5);
  RooRealVar sigma("sigma", "sigma", 0.15, 0.05, 0.3);
  RooGaussian gauss("gauss", "signal", x, mean, sigma);
  //RooRealVar nsig("nsig", "signal yield", 750, 0, 1500);
*/

  // ---------- MODEL 3: Chebychev Polynomial ----------
  RooRealVar c1_2("c1_2", "c1", 0.1, -1, 1);
  RooRealVar c2_2("c2_2", "c2", 0.1, -1, 1);
  RooRealVar c3_2("c3_2", "c3", 0.1, -1, 1);
  RooChebychev cheb2("cheb2", "chebyshev", x, RooArgList(c1_2, c2_2, c3_2));

  RooRealVar c1_25("c1_25", "c1", 0.1, -1, 1);
  RooRealVar c2_25("c2_25", "c2", 0.1, -1, 1);
  RooRealVar c3_25("c3_25", "c3", 0.1, -1, 1);
  RooChebychev cheb25("cheb25", "chebyshev", x, RooArgList(c1_25, c2_25, c3_25));

  RooRealVar c1_3("c1_3", "c1", 0.1, -1, 1);
  RooRealVar c2_3("c2_3", "c2", 0.1, -1, 1);
  RooRealVar c3_3("c3_3", "c3", 0.1, -1, 1);
  RooChebychev cheb3("cheb3", "chebyshev", x, RooArgList(c1_3, c2_3, c3_3));
  //RooAddPdf cheb3("model_cheb", "gauss+cheb", RooArgList(gauss, cheb33));

  RooRealVar c1_35("c1_35", "c1", 0.1, -1, 1);
  RooRealVar c2_35("c2_35", "c2", 0.1, -1, 1);
  RooRealVar c3_35("c3_35", "c3", 0.1, -1, 1);
  RooChebychev cheb35("cheb35", "chebyshev", x, RooArgList(c1_35, c2_35, c3_35));

  //RooRealVar nbkg3("nbkg3", "background yield", 1e5, 0, 1e7);

  //RooAddPdf model_cheb("model_cheb", "gauss+cheb", RooArgList(gauss, cheb), RooArgList(nsig, nbkg3));
  
  // ------------------------
  // Fit Models
  // ------------------------
  cheb2.fitTo(data2, PrintLevel(-1));
  cheb25.fitTo(data25, PrintLevel(-1));
  cheb3.fitTo(data3, PrintLevel(-1));
  cheb35.fitTo(data35, PrintLevel(-1));

  // ------------------------
  // Plotting
  // ------------------------
  RooPlot* frame2 = x.frame(Title("Cheb2 RooFit"));
  RooPlot* frame25 = x.frame(Title("Cheb25 RooFit"));
  RooPlot* frame3 = x.frame(Title("Cheb3 RooFit"));
  RooPlot* frame35 = x.frame(Title("Cheb35 RooFit"));
  
  data2.plotOn(frame2, Name("data2"));
  data25.plotOn(frame25, Name("data25"));
  data3.plotOn(frame3, Name("data3"));
  data35.plotOn(frame35, Name("data35"));

  // Full models
  auto fitResult2 = cheb2.fitTo(data2, Save(true), PrintLevel(-1), Minos(true));
  cheb2.plotOn(frame2, VisualizeError(*fitResult2, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  cheb2.plotOn(frame2, VisualizeError(*fitResult2, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  cheb2.plotOn(frame2, LineColor(kBlue), Name("cheb2"));
  auto fitResult25 = cheb25.fitTo(data25, Save(true), PrintLevel(-1), Minos(true));
  cheb25.plotOn(frame25, VisualizeError(*fitResult25, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  cheb25.plotOn(frame25, VisualizeError(*fitResult25, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  cheb25.plotOn(frame25, LineColor(kBlue), Name("cheb25"));
  auto fitResult3 = cheb3.fitTo(data3, Save(true), PrintLevel(-1), Minos(true));

  if (!fitResult3 || fitResult3->covQual() < 2) {
  std::cerr << "fitResult3 has no valid covariance matrix (covQual = "
            << (fitResult3 ? fitResult3->covQual() : -1) << ")." << std::endl;
  } else {
  std::cout << "fitResult3 is OK. covQual = " << fitResult3->covQual() << std::endl;
  fitResult3->Print("v");  // verbose print
  }

  cheb3.plotOn(frame3, VisualizeError(*fitResult3, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  cheb3.plotOn(frame3, VisualizeError(*fitResult3, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  cheb3.plotOn(frame3, LineColor(kBlue), Name("cheb3"));
  auto fitResult35 = cheb35.fitTo(data35, Save(true), PrintLevel(-1), Minos(true));
  cheb35.plotOn(frame35, VisualizeError(*fitResult35, 2), FillColor(kGray + 2), FillStyle(1001), DrawOption("F"));
  cheb35.plotOn(frame35, VisualizeError(*fitResult35, 1), FillColor(kGray), FillStyle(1001), DrawOption("F"));
  cheb35.plotOn(frame35, LineColor(kBlue), Name("cheb35"));

  // Observable set
  RooArgSet observables(x);
  int nParams = 3; // 3 cheb parameters

  // Chi²/NDOF from the RooPlot
  double chi2_cheb2 = frame2->chiSquare("cheb2", "data2", nParams);
  double chi2_cheb25 = frame25->chiSquare("cheb25", "data25", nParams);
  double chi2_cheb3 = frame3->chiSquare("cheb3", "data3", nParams);
  double chi2_cheb35 = frame35->chiSquare("cheb35", "data35", nParams);

  // NLL from RooDataHist
  auto nll2 = std::unique_ptr<RooAbsReal>(cheb2.createNLL(data2));
  auto nll25 = std::unique_ptr<RooAbsReal>(cheb25.createNLL(data25));
  auto nll3 = std::unique_ptr<RooAbsReal>(cheb3.createNLL(data3));
  auto nll35 = std::unique_ptr<RooAbsReal>(cheb35.createNLL(data35));

  // Print results
  std::cout << "Chi²/NDOF for Chebyshev2 Fit:  " << chi2_cheb2 << std::endl;
  std::cout << "Chi²/NDOF for Chebyshev25 Fit:  " << chi2_cheb25 << std::endl;
  std::cout << "Chi²/NDOF for Chebyshev3 Fit:  " << chi2_cheb3 << std::endl;
  std::cout << "Chi²/NDOF for Chebyshev35 Fit:  " << chi2_cheb35 << std::endl;
  std::cout << "--------------------------------------------" << std::endl;  
  std::cout << "NLL for Chebyshev2 Fit:  " << nll2->getVal() << std::endl;
  std::cout << "NLL for Chebyshev25 Fit:  " << nll25->getVal() << std::endl;
  std::cout << "NLL for Chebyshev3 Fit:  " << nll3->getVal() << std::endl;
  std::cout << "NLL for Chebyshev35 Fit:  " << nll35->getVal() << std::endl;

  // Legend
  auto legend2 = new TLegend(0.6, 0.7, 0.88, 0.88);
  auto legend25 = new TLegend(0.6, 0.7, 0.88, 0.88);
  auto legend3 = new TLegend(0.6, 0.7, 0.88, 0.88);
  auto legend35 = new TLegend(0.6, 0.7, 0.88, 0.88);
  
  legend2->AddEntry(frame2->findObject("cheb2"), "Cheb", "l");
  legend25->AddEntry(frame25->findObject("cheb25"), "Cheb", "l");
  legend3->AddEntry(frame3->findObject("cheb3"), "Cheb", "l");
  legend35->AddEntry(frame35->findObject("cheb35"), "Cheb", "l");

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

  TString base = "GaussCheb_SignalRegion";

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
  
  // Pulls for SR3
  RooPlot* pullFrame3 = x.frame(Title("Pull SR3"));
  RooHist* pulls3 = frame3->pullHist("data3", "cheb3");
  pullFrame3->addPlotable(pulls3, "P");  // Safe if RooHist is fully declared
  TCanvas cpull3("cpull3", "Pull SR3", 800, 400);
  pullFrame3->Draw();
  cpull3.SaveAs(base + "_SR3.pdf");

}
