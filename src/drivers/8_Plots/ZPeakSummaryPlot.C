// Example:
// root -l
// .L ZPeakSummaryPlot.C
// ZPeakSummaryPlot("/eos/user/d/dblend/StackedPlots/12022/Mass_Resolutions/ZCR/KIT/noJERC_noJetVeto_noBtagSF_noZPtRW/ZPeakSummary_12022.txt", "12022");

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct PeakPoint {
  std::string region;
  double meanData = 0.0;
  double errData  = 0.0;
  double meanMC   = 0.0;
  double errMC    = 0.0;
  double ratio    = 0.0;
  double errRatio = 0.0;
  double diffMeV  = 0.0;
  double errDiffMeV = 0.0;
};

void ZPeakSummaryPlot(const char* txtFile,
                      const char* era = "12022",
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

  std::vector<PeakPoint> pts;
  std::string region;
  PeakPoint p;

  while (fin >> region) {
    if (region.size() && region[0] == '#') {
      std::string dummy;
      std::getline(fin, dummy);
      continue;
    }

    p.region = region;
    fin >> p.meanData >> p.errData
        >> p.meanMC   >> p.errMC
        >> p.ratio    >> p.errRatio
        >> p.diffMeV  >> p.errDiffMeV;

    if (!fin.fail()) pts.push_back(p);
  }
  fin.close();

  if (pts.empty()) {
    std::cerr << "[ERROR] No valid points found in " << txtFile << "\n";
    return;
  }

  const int n = (int)pts.size();

  std::vector<double> x(n), ex(n, 0.0);
  std::vector<double> yData(n), eyData(n);
  std::vector<double> yMC(n),   eyMC(n);
  std::vector<double> yBot(n),  eyBot(n);

  double ymin = 1e9, ymax = -1e9;

  for (int i = 0; i < n; ++i) {
    x[i] = i + 1.0;

    yData[i] = pts[i].meanData;
    eyData[i] = pts[i].errData;

    yMC[i] = pts[i].meanMC;
    eyMC[i] = pts[i].errMC;

    ymin = std::min(ymin, pts[i].meanData - pts[i].errData);
    ymin = std::min(ymin, pts[i].meanMC   - pts[i].errMC);
    ymax = std::max(ymax, pts[i].meanData + pts[i].errData);
    ymax = std::max(ymax, pts[i].meanMC   + pts[i].errMC);

    if (useRatioBottom) {
      yBot[i]  = pts[i].ratio;
      eyBot[i] = pts[i].errRatio;
    } else {
      yBot[i]  = pts[i].diffMeV;
      eyBot[i] = pts[i].errDiffMeV;
    }
  }

  double padFrac = 0.25;
  TCanvas* c = new TCanvas("c_ZPeakSummary", "ZPeakSummary", 900, 700);

  TPad* pad1 = new TPad("pad1", "top",    0, padFrac, 1, 1);
  TPad* pad2 = new TPad("pad2", "bottom", 0, 0,       1, padFrac);

  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->SetRightMargin(0.04);

  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.32);
  pad2->SetLeftMargin(0.12);
  pad2->SetRightMargin(0.04);

  pad1->Draw();
  pad2->Draw();

  // ---------------- TOP PAD ----------------
  pad1->cd();

  double dy = ymax - ymin;
  if (dy <= 0) dy = 0.1;
  double topMin = ymin - 0.25 * dy;
  double topMax = ymax + 0.35 * dy;

  TH1F* hFrameTop = new TH1F("hFrameTop", "", n, 0.5, n + 0.5);
  hFrameTop->SetMinimum(topMin);
  hFrameTop->SetMaximum(topMax);
  hFrameTop->GetYaxis()->SetTitle("Z-boson peak [GeV]");
  hFrameTop->GetYaxis()->SetTitleSize(0.05);
  hFrameTop->GetYaxis()->SetLabelSize(0.04);
  hFrameTop->GetYaxis()->SetTitleOffset(1.1);
  hFrameTop->GetXaxis()->SetLabelSize(0.0);

  for (int i = 0; i < n; ++i) {
    hFrameTop->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());
  }

  hFrameTop->Draw();

  TGraphErrors* gData = new TGraphErrors(n, x.data(), yData.data(), ex.data(), eyData.data());
  TGraphErrors* gMC   = new TGraphErrors(n, x.data(), yMC.data(),   ex.data(), eyMC.data());

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
  latex.SetTextSize(0.045);
  latex.SetTextFont(42);
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

  TLegend* leg = new TLegend(0.68, 0.78, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gData, "Data", "lp");
  leg->AddEntry(gMC,   "Simulation", "lp");
  leg->Draw();

  // ---------------- BOTTOM PAD ----------------
  pad2->cd();

  TH1F* hFrameBot = new TH1F("hFrameBot", "", n, 0.5, n + 0.5);

  if (useRatioBottom) {
    hFrameBot->SetMinimum(0.995);
    hFrameBot->SetMaximum(1.005);
    hFrameBot->GetYaxis()->SetTitle("Data/MC");
  } else {
    double minDiff = *std::min_element(yBot.begin(), yBot.end());
    double maxDiff = *std::max_element(yBot.begin(), yBot.end());
    double dd = maxDiff - minDiff;
    if (dd < 5.0) dd = 5.0;
    hFrameBot->SetMinimum(minDiff - 0.4 * dd);
    hFrameBot->SetMaximum(maxDiff + 0.4 * dd);
    hFrameBot->GetYaxis()->SetTitle("#Delta peak [MeV]");
  }

  hFrameBot->GetYaxis()->SetTitleSize(0.11);
  hFrameBot->GetYaxis()->SetTitleOffset(0.45);
  hFrameBot->GetYaxis()->SetLabelSize(0.09);
  hFrameBot->GetYaxis()->SetNdivisions(505);

  hFrameBot->GetXaxis()->SetTitle("Region");
  hFrameBot->GetXaxis()->SetTitleSize(0.12);
  hFrameBot->GetXaxis()->SetTitleOffset(1.0);
  hFrameBot->GetXaxis()->SetLabelSize(0.11);

  for (int i = 0; i < n; ++i) {
    hFrameBot->GetXaxis()->SetBinLabel(i + 1, pts[i].region.c_str());
  }

  hFrameBot->Draw();

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

  // output name
  TString outPath = txtFile;
  TString outDir = gSystem->DirName(outPath);
  TString suffix = useRatioBottom ? "_ratio" : "_diffMeV";
  TString outBase = outDir + "/ZPeakSummaryPlot_" + TString(era) + suffix;

  c->SaveAs(outBase + ".pdf");
  c->SaveAs(outBase + ".png");

  std::cout << "[INFO] Saved:\n"
            << "  " << outBase << ".pdf\n"
            << "  " << outBase << ".png\n";
}
