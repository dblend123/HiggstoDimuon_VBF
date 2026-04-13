// ================== ZPeakAndResolutionSummaryPlot.C ==================
// Example:
// root -l
// .L ZPeakAndResolutionSummaryPlot.C
// ZPeakAndResolutionSummaryPlot(
//   "/eos/user/d/dblend/StackedPlots/12022/Mass_Resolutions/ZCR/KIT/noJERC_noJetVeto_noBtagSF_noZPtRW/ZResolutionSummary_12022.txt",
//   "12022",
//   "sigma68",
//   true
// );
//
// metric options: "sigma68", "sigmaCB", "sigmaFwhm"
// useRatioBottom = true  -> bottom pads show Data/MC
// useRatioBottom = false -> bottom pads show Data-MC difference in MeV

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

static TString eraLumiText(const char* era) {
  const std::string e = era ? std::string(era) : "";
  if      (e == "12022") return "8.0 fb^{-1} (13.6 TeV)";
  else if (e == "22022") return "26.7 fb^{-1} (13.6 TeV)";
  else if (e == "12023") return "18.7 fb^{-1} (13.6 TeV)";
  else if (e == "22023") return "9.7 fb^{-1} (13.6 TeV)";
  else if (e == "2024")  return "104.8 fb^{-1} (13.6 TeV)";
  return TString(era);
}

void ZPeakAndResolutionSummaryPlot(const char* txtFile,
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
    std::cerr << "[ERROR] No valid points found in " << txtFile << "\n";
    return;
  }

  const std::string met = metric ? std::string(metric) : "sigma68";
  if (met != "sigma68" && met != "sigmaCB" && met != "sigmaFwhm") {
    std::cerr << "[ERROR] metric must be sigma68, sigmaCB, or sigmaFwhm\n";
    return;
  }

  const int n = (int)pts.size();

  std::vector<double> x(n), ex(n, 0.0);

  std::vector<double> yPeakData(n), eyPeakData(n), yPeakMC(n), eyPeakMC(n), yPeakBot(n), eyPeakBot(n);
  std::vector<double> yResData(n),  eyResData(n),  yResMC(n),  eyResMC(n),  yResBot(n),  eyResBot(n);

  double peakMin = 1e9, peakMax = -1e9;
  double resMin  = 1e9, resMax  = -1e9;

  for (int i = 0; i < n; ++i) {
    x[i] = i + 1.0;

    // peak
    yPeakData[i]  = pts[i].meanData;
    eyPeakData[i] = pts[i].errMeanData;
    yPeakMC[i]    = pts[i].meanMC;
    eyPeakMC[i]   = pts[i].errMeanMC;
    yPeakBot[i]   = useRatioBottom ? pts[i].ratioPeak : pts[i].diffPeakMeV;
    eyPeakBot[i]  = useRatioBottom ? pts[i].errRatioPeak : pts[i].errDiffPeakMeV;

    peakMin = std::min(peakMin, yPeakData[i] - eyPeakData[i]);
    peakMin = std::min(peakMin, yPeakMC[i]   - eyPeakMC[i]);
    peakMax = std::max(peakMax, yPeakData[i] + eyPeakData[i]);
    peakMax = std::max(peakMax, yPeakMC[i]   + eyPeakMC[i]);

    // resolution
    if (met == "sigma68") {
      yResData[i]  = pts[i].sigma68Data;
      eyResData[i] = pts[i].errSigma68Data;
      yResMC[i]    = pts[i].sigma68MC;
      eyResMC[i]   = pts[i].errSigma68MC;
      yResBot[i]   = useRatioBottom ? pts[i].ratioSigma68 : pts[i].diffSigma68MeV;
      eyResBot[i]  = useRatioBottom ? pts[i].errRatioSigma68 : pts[i].errDiffSigma68MeV;
    } else if (met == "sigmaCB") {
      yResData[i]  = pts[i].sigmaCBData;
      eyResData[i] = pts[i].errSigmaCBData;
      yResMC[i]    = pts[i].sigmaCBMC;
      eyResMC[i]   = pts[i].errSigmaCBMC;
      yResBot[i]   = useRatioBottom ? pts[i].ratioSigmaCB : pts[i].diffSigmaCBMeV;
      eyResBot[i]  = useRatioBottom ? pts[i].errRatioSigmaCB : pts[i].errDiffSigmaCBMeV;
    } else {
      yResData[i]  = pts[i].sigmaFwhmData;
      eyResData[i] = 0.0;
      yResMC[i]    = pts[i].sigmaFwhmMC;
      eyResMC[i]   = 0.0;
      yResBot[i]   = useRatioBottom ? pts[i].ratioSigmaFwhm : pts[i].diffSigmaFwhmMeV;
      eyResBot[i]  = 0.0;
    }

    resMin = std::min(resMin, yResData[i] - eyResData[i]);
    resMin = std::min(resMin, yResMC[i]   - eyResMC[i]);
    resMax = std::max(resMax, yResData[i] + eyResData[i]);
    resMax = std::max(resMax, yResMC[i]   + eyResMC[i]);
  }

  TCanvas* c = new TCanvas("c_ZPeakAndResolutionSummary", "ZPeakAndResolutionSummary", 1400, 700);

  const double leftMargin   = 0.08;
  const double rightMargin  = 0.03;
  const double midGap       = 0.04;
  const double topFrac      = 0.73;
  const double bottomFrac   = 0.27;

  // left column
  TPad* padPeakTop = new TPad("padPeakTop", "padPeakTop", leftMargin, bottomFrac, 0.5 - midGap/2.0, 1.0);
  TPad* padPeakBot = new TPad("padPeakBot", "padPeakBot", leftMargin, 0.0,        0.5 - midGap/2.0, bottomFrac);

  // right column
  TPad* padResTop  = new TPad("padResTop",  "padResTop",  0.5 + midGap/2.0, bottomFrac, 1.0 - rightMargin, 1.0);
  TPad* padResBot  = new TPad("padResBot",  "padResBot",  0.5 + midGap/2.0, 0.0,        1.0 - rightMargin, bottomFrac);

  padPeakTop->SetBottomMargin(0.02);
  padPeakTop->SetLeftMargin(0.14);
  padPeakTop->SetRightMargin(0.04);

  padPeakBot->SetTopMargin(0.04);
  padPeakBot->SetBottomMargin(0.30);
  padPeakBot->SetLeftMargin(0.14);
  padPeakBot->SetRightMargin(0.04);

  padResTop->SetBottomMargin(0.02);
  padResTop->SetLeftMargin(0.12);
  padResTop->SetRightMargin(0.04);

  padResBot->SetTopMargin(0.04);
  padResBot->SetBottomMargin(0.30);
  padResBot->SetLeftMargin(0.12);
  padResBot->SetRightMargin(0.04);

  padPeakTop->Draw();
  padPeakBot->Draw();
  padResTop->Draw();
  padResBot->Draw();

  // ---------------- Peak top ----------------
  padPeakTop->cd();

  double dPeak = peakMax - peakMin;
  if (dPeak <= 0) dPeak = 0.1;

  TH1F* hPeakTop = new TH1F("hPeakTop", "", n, 0.5, n + 0.5);
  hPeakTop->SetMinimum(peakMin - 0.25 * dPeak);
  hPeakTop->SetMaximum(peakMax + 0.35 * dPeak);
  hPeakTop->GetYaxis()->SetTitle("Z-boson peak [GeV]");
  hPeakTop->GetYaxis()->SetTitleSize(0.055);
  hPeakTop->GetYaxis()->SetLabelSize(0.045);
  hPeakTop->GetYaxis()->SetTitleOffset(1.05);
  hPeakTop->GetXaxis()->SetLabelSize(0.0);
  for (int i = 0; i < n; ++i) hPeakTop->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());
  hPeakTop->Draw();

  TGraphErrors* gPeakData = new TGraphErrors(n, x.data(), yPeakData.data(), ex.data(), eyPeakData.data());
  TGraphErrors* gPeakMC   = new TGraphErrors(n, x.data(), yPeakMC.data(),   ex.data(), eyPeakMC.data());

  gPeakData->SetMarkerStyle(20);
  gPeakData->SetMarkerSize(1.0);
  gPeakData->SetLineColor(kBlack);
  gPeakData->SetMarkerColor(kBlack);
  gPeakData->SetLineWidth(2);

  gPeakMC->SetMarkerStyle(21);
  gPeakMC->SetMarkerSize(1.0);
  gPeakMC->SetLineColor(kRed);
  gPeakMC->SetMarkerColor(kRed);
  gPeakMC->SetLineWidth(2);

  gPeakData->Draw("P SAME");
  gPeakMC->Draw("P SAME");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.05);
  latex.DrawLatex(0.02, 0.93, "CMS");

  latex.SetTextAlign(31);
  latex.SetTextSize(0.043);
  latex.DrawLatex(0.98, 0.93, eraLumiText(era));

  latex.SetTextAlign(13);
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.17, 0.84, "Peak");

  TLegend* legPeak = new TLegend(0.63, 0.74, 0.93, 0.90);
  legPeak->SetBorderSize(0);
  legPeak->SetFillStyle(0);
  legPeak->AddEntry(gPeakData, "Data", "lp");
  legPeak->AddEntry(gPeakMC,   "Simulation", "lp");
  legPeak->Draw();

  // ---------------- Peak bottom ----------------
  padPeakBot->cd();

  TH1F* hPeakBot = new TH1F("hPeakBot", "", n, 0.5, n + 0.5);
  if (useRatioBottom) {
    hPeakBot->SetMinimum(0.995);
    hPeakBot->SetMaximum(1.005);
    hPeakBot->GetYaxis()->SetTitle("Data/MC");
  } else {
    double minB = *std::min_element(yPeakBot.begin(), yPeakBot.end());
    double maxB = *std::max_element(yPeakBot.begin(), yPeakBot.end());
    double dB = maxB - minB;
    if (dB < 5.0) dB = 5.0;
    hPeakBot->SetMinimum(minB - 0.4 * dB);
    hPeakBot->SetMaximum(maxB + 0.4 * dB);
    hPeakBot->GetYaxis()->SetTitle("#Delta peak [MeV]");
  }
  hPeakBot->GetYaxis()->SetTitleSize(0.11);
  hPeakBot->GetYaxis()->SetTitleOffset(0.50);
  hPeakBot->GetYaxis()->SetLabelSize(0.09);
  hPeakBot->GetYaxis()->SetNdivisions(505);
  hPeakBot->GetXaxis()->SetTitle("Region");
  hPeakBot->GetXaxis()->SetTitleSize(0.12);
  hPeakBot->GetXaxis()->SetTitleOffset(1.0);
  hPeakBot->GetXaxis()->SetLabelSize(0.11);
  for (int i = 0; i < n; ++i) hPeakBot->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());
  hPeakBot->Draw();

  TGraphErrors* gPeakBottom = new TGraphErrors(n, x.data(), yPeakBot.data(), ex.data(), eyPeakBot.data());
  gPeakBottom->SetMarkerStyle(20);
  gPeakBottom->SetMarkerSize(0.95);
  gPeakBottom->SetLineColor(kBlack);
  gPeakBottom->SetMarkerColor(kBlack);
  gPeakBottom->SetLineWidth(2);
  gPeakBottom->Draw("P SAME");

  TLine* linePeak = new TLine(0.5, useRatioBottom ? 1.0 : 0.0, n + 0.5, useRatioBottom ? 1.0 : 0.0);
  linePeak->SetLineStyle(2);
  linePeak->SetLineWidth(2);
  linePeak->Draw("SAME");

  // ---------------- Resolution top ----------------
  padResTop->cd();

  double dRes = resMax - resMin;
  if (dRes <= 0) dRes = 0.1;

  TH1F* hResTop = new TH1F("hResTop", "", n, 0.5, n + 0.5);
  hResTop->SetMinimum(resMin - 0.25 * dRes);
  hResTop->SetMaximum(resMax + 0.35 * dRes);
  hResTop->GetYaxis()->SetTitle("Z-boson resolution [GeV]");
  hResTop->GetYaxis()->SetTitleSize(0.055);
  hResTop->GetYaxis()->SetLabelSize(0.045);
  hResTop->GetYaxis()->SetTitleOffset(1.05);
  hResTop->GetXaxis()->SetLabelSize(0.0);
  for (int i = 0; i < n; ++i) hResTop->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());
  hResTop->Draw();

  TGraphErrors* gResData = new TGraphErrors(n, x.data(), yResData.data(), ex.data(), eyResData.data());
  TGraphErrors* gResMC   = new TGraphErrors(n, x.data(), yResMC.data(),   ex.data(), eyResMC.data());

  gResData->SetMarkerStyle(20);
  gResData->SetMarkerSize(1.0);
  gResData->SetLineColor(kBlack);
  gResData->SetMarkerColor(kBlack);
  gResData->SetLineWidth(2);

  gResMC->SetMarkerStyle(21);
  gResMC->SetMarkerSize(1.0);
  gResMC->SetLineColor(kRed);
  gResMC->SetMarkerColor(kRed);
  gResMC->SetLineWidth(2);

  gResData->Draw("P SAME");
  gResMC->Draw("P SAME");

  TString metricLabel = (met == "sigma68" ? "#sigma_{68}" :
                        (met == "sigmaCB" ? "DCB #sigma" : "FWHM/2.355"));

  latex.SetTextAlign(13);
  latex.SetTextSize(0.045);
  latex.DrawLatex(0.17, 0.84, Form("Resolution (%s)", metricLabel.Data()));

  TLegend* legRes = new TLegend(0.56, 0.70, 0.93, 0.90);
  legRes->SetBorderSize(0);
  legRes->SetFillStyle(0);
  legRes->AddEntry(gResData, "Data", "lp");
  legRes->AddEntry(gResMC,   "Simulation", "lp");
  legRes->AddEntry((TObject*)0, metricLabel.Data(), "");
  legRes->Draw();

  // ---------------- Resolution bottom ----------------
  padResBot->cd();

  TH1F* hResBot = new TH1F("hResBot", "", n, 0.5, n + 0.5);
  if (useRatioBottom) {
    hResBot->SetMinimum(0.8);
    hResBot->SetMaximum(1.2);
    hResBot->GetYaxis()->SetTitle("Data/MC");
  } else {
    double minB = *std::min_element(yResBot.begin(), yResBot.end());
    double maxB = *std::max_element(yResBot.begin(), yResBot.end());
    double dB = maxB - minB;
    if (dB < 20.0) dB = 20.0;
    hResBot->SetMinimum(minB - 0.4 * dB);
    hResBot->SetMaximum(maxB + 0.4 * dB);
    hResBot->GetYaxis()->SetTitle("#Delta res. [MeV]");
  }
  hResBot->GetYaxis()->SetTitleSize(0.11);
  hResBot->GetYaxis()->SetTitleOffset(0.50);
  hResBot->GetYaxis()->SetLabelSize(0.09);
  hResBot->GetYaxis()->SetNdivisions(505);
  hResBot->GetXaxis()->SetTitle("Region");
  hResBot->GetXaxis()->SetTitleSize(0.12);
  hResBot->GetXaxis()->SetTitleOffset(1.0);
  hResBot->GetXaxis()->SetLabelSize(0.11);
  for (int i = 0; i < n; ++i) hResBot->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());
  hResBot->Draw();

  TGraphErrors* gResBottom = new TGraphErrors(n, x.data(), yResBot.data(), ex.data(), eyResBot.data());
  gResBottom->SetMarkerStyle(20);
  gResBottom->SetMarkerSize(0.95);
  gResBottom->SetLineColor(kBlack);
  gResBottom->SetMarkerColor(kBlack);
  gResBottom->SetLineWidth(2);
  gResBottom->Draw("P SAME");

  TLine* lineRes = new TLine(0.5, useRatioBottom ? 1.0 : 0.0, n + 0.5, useRatioBottom ? 1.0 : 0.0);
  lineRes->SetLineStyle(2);
  lineRes->SetLineWidth(2);
  lineRes->Draw("SAME");

  TString outPath = txtFile;
  TString outDir = gSystem->DirName(outPath);
  TString suffix = useRatioBottom ? "_ratio" : "_diffMeV";
  TString outBase = outDir + "/ZPeakAndResolutionSummaryPlot_" + TString(era) + "_" + metric + suffix;

  c->SaveAs(outBase + ".pdf");
  c->SaveAs(outBase + ".png");

  std::cout << "[INFO] Saved:\n"
            << "  " << outBase << ".pdf\n"
            << "  " << outBase << ".png\n";
}
