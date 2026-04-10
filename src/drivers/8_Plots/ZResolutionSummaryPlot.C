// Example:
// root -l
// .L ZResolutionSummaryPlot.C
// ZResolutionSummaryPlot("/eos/user/d/dblend/StackedPlots/12022/Mass_Resolutions/ZCR/KIT/noJERC_noJetVeto_noBtagSF_noZPtRW/ZPeakSummary_12022.txt", "12022", "sigma68", true);

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct ResPoint {
  std::string region;

  double meanData, errMeanData, meanMC, errMeanMC, ratioPeak, errRatioPeak, diffPeakMeV, errDiffPeakMeV;

  double sigma68Data, errSigma68Data, sigma68MC, errSigma68MC, ratioSigma68, errRatioSigma68, diffSigma68MeV, errDiffSigma68MeV;

  double sigmaCBData, errSigmaCBData, sigmaCBMC, errSigmaCBMC, ratioSigmaCB, errRatioSigmaCB, diffSigmaCBMeV, errDiffSigmaCBMeV;

  double sigmaFwhmData, sigmaFwhmMC, ratioSigmaFwhm, diffSigmaFwhmMeV;
};

void ZResolutionSummaryPlot(const char* txtFile,
                            const char* era = "12022",
                            const char* metric = "sigma68",
                            bool useRatioBottom = true) {
  if (!txtFile || std::string(txtFile).empty()) {
    std::cerr << "[ERROR] txtFile is empty\n";
    return;
  }

  std::ifstream fin(txtFile);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Could not open summary file: " << txtFile << "\n";
    return;
  }

  std::vector<ResPoint> pts;
  std::string first;
  while (fin >> first) {
    if (!first.empty() && first[0] == '#') {
      std::string rest;
      std::getline(fin, rest);
      continue;
    }

    ResPoint p;
    p.region = first;

    fin >> p.meanData >> p.errMeanData >> p.meanMC >> p.errMeanMC >> p.ratioPeak >> p.errRatioPeak >> p.diffPeakMeV >> p.errDiffPeakMeV
        >> p.sigma68Data >> p.errSigma68Data >> p.sigma68MC >> p.errSigma68MC >> p.ratioSigma68 >> p.errRatioSigma68 >> p.diffSigma68MeV >> p.errDiffSigma68MeV
        >> p.sigmaCBData >> p.errSigmaCBData >> p.sigmaCBMC >> p.errSigmaCBMC >> p.ratioSigmaCB >> p.errRatioSigmaCB >> p.diffSigmaCBMeV >> p.errDiffSigmaCBMeV
        >> p.sigmaFwhmData >> p.sigmaFwhmMC >> p.ratioSigmaFwhm >> p.diffSigmaFwhmMeV;

    if (!fin.fail()) pts.push_back(p);
  }
  fin.close();

  if (pts.empty()) {
    std::cerr << "[ERROR] No valid points found\n";
    return;
  }

  const std::string met = metric ? std::string(metric) : "sigma68";
  const int n = (int)pts.size();

  std::vector<double> x(n), ex(n, 0.0);
  std::vector<double> yData(n), eyData(n), yMC(n), eyMC(n), yBot(n), eyBot(n);

  double ymin = 1e9, ymax = -1e9;

  for (int i = 0; i < n; ++i) {
    x[i] = i + 1.0;

    if (met == "sigma68") {
      yData[i] = pts[i].sigma68Data;
      eyData[i] = pts[i].errSigma68Data;
      yMC[i] = pts[i].sigma68MC;
      eyMC[i] = pts[i].errSigma68MC;
      yBot[i] = useRatioBottom ? pts[i].ratioSigma68 : pts[i].diffSigma68MeV;
      eyBot[i] = useRatioBottom ? pts[i].errRatioSigma68 : pts[i].errDiffSigma68MeV;
    } else if (met == "sigmaCB") {
      yData[i] = pts[i].sigmaCBData;
      eyData[i] = pts[i].errSigmaCBData;
      yMC[i] = pts[i].sigmaCBMC;
      eyMC[i] = pts[i].errSigmaCBMC;
      yBot[i] = useRatioBottom ? pts[i].ratioSigmaCB : pts[i].diffSigmaCBMeV;
      eyBot[i] = useRatioBottom ? pts[i].errRatioSigmaCB : pts[i].errDiffSigmaCBMeV;
    } else if (met == "sigmaFwhm") {
      yData[i] = pts[i].sigmaFwhmData;
      eyData[i] = 0.0;
      yMC[i] = pts[i].sigmaFwhmMC;
      eyMC[i] = 0.0;
      yBot[i] = useRatioBottom ? pts[i].ratioSigmaFwhm : pts[i].diffSigmaFwhmMeV;
      eyBot[i] = 0.0;
    } else {
      std::cerr << "[ERROR] metric must be sigma68, sigmaCB, or sigmaFwhm\n";
      return;
    }

    ymin = std::min(ymin, yData[i] - eyData[i]);
    ymin = std::min(ymin, yMC[i] - eyMC[i]);
    ymax = std::max(ymax, yData[i] + eyData[i]);
    ymax = std::max(ymax, yMC[i] + eyMC[i]);
  }

  TCanvas* c = new TCanvas("c_ZResolutionSummary", "ZResolutionSummary", 900, 700);

  TPad* pad1 = new TPad("pad1", "top", 0, 0.25, 1, 1);
  TPad* pad2 = new TPad("pad2", "bottom", 0, 0, 1, 0.25);

  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.04);

  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.32);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.04);

  pad1->Draw();
  pad2->Draw();

  pad1->cd();

  double dy = ymax - ymin;
  if (dy <= 0) dy = 0.1;
  double topMin = ymin - 0.25 * dy;
  double topMax = ymax + 0.35 * dy;

  TH1F* hTop = new TH1F("hTop", "", n, 0.5, n + 0.5);
  hTop->SetMinimum(topMin);
  hTop->SetMaximum(topMax);
  hTop->GetYaxis()->SetTitle("Z-boson resolution [GeV]");
  hTop->GetYaxis()->SetTitleSize(0.05);
  hTop->GetYaxis()->SetLabelSize(0.04);
  hTop->GetYaxis()->SetTitleOffset(1.1);
  hTop->GetXaxis()->SetLabelSize(0.0);

  for (int i = 0; i < n; ++i)
    hTop->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());

  hTop->Draw();

  TGraphErrors* gData = new TGraphErrors(n, x.data(), yData.data(), ex.data(), eyData.data());
  TGraphErrors* gMC   = new TGraphErrors(n, x.data(), yMC.data(), ex.data(), eyMC.data());

  gData->SetMarkerStyle(20);
  gData->SetMarkerSize(1.0);
  gData->SetLineColor(kBlack);
  gData->SetMarkerColor(kBlack);
  gData->SetLineWidth(2);

  gMC->SetMarkerStyle(21);
  gMC->SetMarkerSize(1.0);
  gMC->SetLineColor(kRed);
  gMC->SetMarkerColor(kRed);
  gMC->SetLineWidth(2);

  gData->Draw("P SAME");
  gMC->Draw("P SAME");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.13, 0.92, "CMS");

  TString lumiText;
  if      (std::string(era) == "12022") lumiText = "8.0 fb^{-1} (13.6 TeV)";
  else if (std::string(era) == "22022") lumiText = "26.7 fb^{-1} (13.6 TeV)";
  else if (std::string(era) == "12023") lumiText = "18.7 fb^{-1} (13.6 TeV)";
  else if (std::string(era) == "22023") lumiText = "9.7 fb^{-1} (13.6 TeV)";
  else if (std::string(era) == "2024")  lumiText = "104.8 fb^{-1} (13.6 TeV)";
  else                                   lumiText = TString(era);

  latex.SetTextSize(0.040);
  latex.SetTextAlign(31);
  latex.DrawLatex(0.96, 0.93, lumiText);

  TString metricLabel = (met == "sigma68" ? "#sigma_{68}" : (met == "sigmaCB" ? "DCB #sigma" : "FWHM/2.355"));

  TLegend* leg = new TLegend(0.62, 0.74, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gData, "Data", "lp");
  leg->AddEntry(gMC, "Simulation", "lp");
  leg->AddEntry((TObject*)0, metricLabel.Data(), "");
  leg->Draw();

  pad2->cd();

  TH1F* hBot = new TH1F("hBot", "", n, 0.5, n + 0.5);
  if (useRatioBottom) {
    hBot->SetMinimum(0.8);
    hBot->SetMaximum(1.2);
    hBot->GetYaxis()->SetTitle("Data/MC");
  } else {
    double minB = *std::min_element(yBot.begin(), yBot.end());
    double maxB = *std::max_element(yBot.begin(), yBot.end());
    double dB = maxB - minB;
    if (dB < 20.0) dB = 20.0;
    hBot->SetMinimum(minB - 0.4 * dB);
    hBot->SetMaximum(maxB + 0.4 * dB);
    hBot->GetYaxis()->SetTitle("#Delta res. [MeV]");
  }

  hBot->GetYaxis()->SetTitleSize(0.11);
  hBot->GetYaxis()->SetTitleOffset(0.45);
  hBot->GetYaxis()->SetLabelSize(0.09);
  hBot->GetYaxis()->SetNdivisions(505);

  hBot->GetXaxis()->SetTitle("Region");
  hBot->GetXaxis()->SetTitleSize(0.12);
  hBot->GetXaxis()->SetLabelSize(0.11);

  for (int i = 0; i < n; ++i)
    hBot->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());

  hBot->Draw();

  TGraphErrors* gBot = new TGraphErrors(n, x.data(), yBot.data(), ex.data(), eyBot.data());
  gBot->SetMarkerStyle(20);
  gBot->SetMarkerSize(0.95);
  gBot->SetLineColor(kBlack);
  gBot->SetMarkerColor(kBlack);
  gBot->SetLineWidth(2);
  gBot->Draw("P SAME");

  double ref = useRatioBottom ? 1.0 : 0.0;
  TLine* line = new TLine(0.5, ref, n + 0.5, ref);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("SAME");

  TString outPath = txtFile;
  TString outDir = gSystem->DirName(outPath);
  TString suffix = useRatioBottom ? "_ratio" : "_diffMeV";
  TString outBase = outDir + "/ZResolutionSummaryPlot_" + TString(era) + "_" + metric + suffix;

  c->SaveAs(outBase + ".pdf");
  c->SaveAs(outBase + ".png");

  std::cout << "[INFO] Saved:\n"
            << "  " << outBase << ".pdf\n"
            << "  " << outBase << ".png\n";
}
