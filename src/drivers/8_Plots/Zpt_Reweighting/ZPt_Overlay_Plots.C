// ======================== ZPt_Overlay_Plots.C =========================
//
// Reads Z pT (DiMuonpt) histograms in 0/1/2+ jet bins:
//   corr_DM_Pre-Selection_Z0j/DiMuonpt
//   corr_DM_Pre-Selection_Z1j/DiMuonpt
//   corr_DM_Pre-Selection_Z2j/DiMuonpt
//
// Your directory layout (example):
//   <baseDir>/DATA/<something>_all.root
//   <baseDir>/DY/<something>_all.root
//   <baseDir>/EWK/<something>_all.root
//   <baseDir>/TT/<something>_all.root
//   <baseDir>/WZto2L2Q/<something>_all.root
//   ... etc.
//
// It groups:
//   VV  = WZto2L2Q, WZto3LNu, ZZto2L2Nu, ZZto2L2Q, ZZto4L
//   VVV = WWW, WWZ, WZZ, ZZZ
//
// It computes per jet bin:
//   SF_DataOverDY = Data/DY (use for Z-pT reweighting)
//   SF_DataOverMC = Data/(DY+EWK+TT+VV+VVV) (diagnostic)
//
// OUTPUT PATH (auto-derived from input path):
//   Input expected to contain ".../Analyzed_SkimmedFiles/<YEAR>/<MuCorr>/<CfgTag>"
//   Plots and SF ROOT file are written to:
//     /eos/user/d/dblend/StackedPlots/<YEAR>/ZPt_Corrections/<MuCorr>/<CfgTag>/
//
// Usage:
//   root -l
//   .L ZPt_Overlay_Plots.C+
//   ZPt_Overlay_Plots("/eos/user/d/dblend/Analyzed_SkimmedFiles/12022/Rochester/noScaleSmear_noBtagSF", "anything", 2);
//
// Note: the second argument is kept for backwards compatibility but is not used for output tags;
//       output tags are taken from the input path.
//
// =====================================================================

#include "TFile.h"
#include "TH1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

// ------------------------- Helper: split path --------------------------
static std::vector<std::string> splitPath(const std::string& p) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : p) {
    if (c == '/') {
      if (!cur.empty()) out.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}

// ------------------ Helper: derive YEAR/MuCorr/Cfg from input path -----
static bool deriveTagsFromInputPath(const std::string& baseDir,
                                   std::string& year,
                                   std::string& muCorrTag,
                                   std::string& cfgTag)
{
  // Expect .../Analyzed_SkimmedFiles/<YEAR>/<MuCorrTag>/<CfgTag>
  auto parts = splitPath(baseDir);

  int idx = -1;
  for (int i = 0; i < (int)parts.size(); ++i) {
    if (parts[i] == "Analyzed_SkimmedFiles") { idx = i; break; }
  }
  if (idx < 0 || idx + 3 >= (int)parts.size()) return false;

  year      = parts[idx + 1];
  muCorrTag = parts[idx + 2];
  cfgTag    = parts[idx + 3];
  return true;
}

// ------------------------- Helper: clone+detach hist -------------------
static TH1* getHist(TFile* f, const std::string& path) {
  if (!f) return nullptr;
  auto* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
  if (!h) return nullptr;
  auto* c = (TH1*)h->Clone();
  c->SetDirectory(nullptr);
  c->Sumw2();
  return c;
}

// ------------------------- Helper: safe divide -------------------------
static void protectDivide(TH1* ratio, const TH1* den) {
  if (!ratio || !den) return;
  for (int i=1;i<=ratio->GetNbinsX();++i) {
    if (den->GetBinContent(i) <= 0) {
      ratio->SetBinContent(i,0);
      ratio->SetBinError(i,0);
    }
  }
}

// ------------------------- Helper: open file ---------------------------
static TFile* openFile(const std::string& path) {
  TFile* f = TFile::Open(path.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] Could not open " << path << "\n";
    return nullptr;
  }
  return f;
}

// -------- Find the first "*_all.root" file inside a folder -------------
static std::string findAllRootInFolder(const std::string& folder) {
  TSystemDirectory dir("d", folder.c_str());
  TList* list = dir.GetListOfFiles();
  if (!list) return "";

  TIter next(list);
  TSystemFile* f = nullptr;
  while ((f = (TSystemFile*)next())) {
    if (f->IsDirectory()) continue;
    std::string name = f->GetName();
    if (name.size() >= 9 && name.substr(name.size()-9) == "_all.root") {
      return folder + "/" + name;
    }
  }
  return "";
}

// ----- Sum histograms over a group of folders (each folder has *_all.root)
static TH1* sumGroupHists(const std::string& baseDir,
                          const std::vector<std::string>& folders,
                          const std::string& histPath,
                          int rebin)
{
  TH1* hsum = nullptr;

  for (const auto& sub : folders) {
    const std::string folder = baseDir + "/" + sub;
    const std::string fpath  = findAllRootInFolder(folder);
    if (fpath.empty()) {
      std::cerr << "[WARN] No *_all.root found in " << folder << "\n";
      continue;
    }

    TFile* f = openFile(fpath);
    if (!f) continue;

    TH1* h = getHist(f, histPath);
    f->Close();
    delete f;

    if (!h) {
      std::cerr << "[WARN] Missing " << histPath << " in " << fpath << "\n";
      continue;
    }

    if (rebin > 1) h->Rebin(rebin);

    if (!hsum) {
      hsum = (TH1*)h->Clone(("SUM_" + sub).c_str());
      hsum->SetDirectory(nullptr);
      hsum->Sumw2();
    } else {
      hsum->Add(h);
    }
    delete h;
  }

  return hsum;
}

// ------------------------- Make one category plot+SF -------------------
static void makeOneCategory(const std::string& outBase,
                            const std::string& catDir,
                            const std::string& histName,
                            TFile* fData, TFile* fDY, TFile* fEWK, TFile* fTT,
                            TH1* hVV, TH1* hVVV,
                            TFile* fout,
                            int rebin)
{
  const std::string histPath = catDir + "/" + histName;

  TH1* hData = getHist(fData, histPath);
  TH1* hDY   = getHist(fDY,   histPath);
  if (!hData || !hDY) { delete hData; delete hDY; return; }

  if (rebin > 1) { hData->Rebin(rebin); hDY->Rebin(rebin); }

  TH1* hEWK = (fEWK ? getHist(fEWK, histPath) : nullptr);
  TH1* hTT  = (fTT  ? getHist(fTT,  histPath) : nullptr);
  if (rebin > 1) {
    if (hEWK) hEWK->Rebin(rebin);
    if (hTT)  hTT->Rebin(rebin);
  }

  // MC total
  TH1* hMC = (TH1*)hDY->Clone(("MCtot_" + catDir).c_str());
  hMC->SetDirectory(nullptr);
  hMC->Sumw2();
  if (hEWK) hMC->Add(hEWK);
  if (hTT)  hMC->Add(hTT);
  if (hVV)  hMC->Add(hVV);
  if (hVVV) hMC->Add(hVVV);

  std::cout << "\n[DEBUG] Category = " << catDir << "\n";
  std::cout << "  hData : bins=" << hData->GetNbinsX()
            << " xmin=" << hData->GetXaxis()->GetXmin()
            << " xmax=" << hData->GetXaxis()->GetXmax() << "\n";
  std::cout << "  hDY   : bins=" << hDY->GetNbinsX()
            << " xmin=" << hDY->GetXaxis()->GetXmin()
            << " xmax=" << hDY->GetXaxis()->GetXmax() << "\n";
  std::cout << "  hMC   : bins=" << hMC->GetNbinsX()
            << " xmin=" << hMC->GetXaxis()->GetXmin()
            << " xmax=" << hMC->GetXaxis()->GetXmax() << "\n";
  if (hEWK) std::cout << "  hEWK  : bins=" << hEWK->GetNbinsX() << "\n";
  if (hTT)  std::cout << "  hTT   : bins=" << hTT->GetNbinsX() << "\n";
  if (hVV)  std::cout << "  hVV   : bins=" << hVV->GetNbinsX() << "\n";
  if (hVVV) std::cout << "  hVVV  : bins=" << hVVV->GetNbinsX() << "\n";

  // SFs
  TH1* sfDY = (TH1*)hData->Clone(("SF_DataOverDY_" + catDir).c_str());
  sfDY->SetDirectory(nullptr);
  protectDivide(sfDY, hDY);
  sfDY->Divide(hDY);

  TH1* sfMC = (TH1*)hData->Clone(("SF_DataOverMC_" + catDir).c_str());
  sfMC->SetDirectory(nullptr);
  protectDivide(sfMC, hMC);
  sfMC->Divide(hMC);

  // ---------------- Plot: stack top, ratio bottom ----------------
  gSystem->Exec(("mkdir -p " + outBase).c_str());

  TCanvas c(("c_" + catDir).c_str(), "", 900, 900);
  TPad p1("p1","",0,0.30,1,1.0);
  TPad p2("p2","",0,0.00,1,0.30);
  p1.SetBottomMargin(0.02);
  p2.SetTopMargin(0.02);
  p2.SetBottomMargin(0.30);
  p1.Draw(); p2.Draw();

  p1.cd();
  THStack hs(("hs_" + catDir).c_str(), "");

  // Stack order: EWK, VVV, VV, TT, DY
  if (hEWK) hs.Add(hEWK);
  if (hVVV) hs.Add(hVVV);
  if (hVV)  hs.Add(hVV);
  if (hTT)  hs.Add(hTT);
  hs.Add(hDY);

  hs.Draw("HIST");
  hs.GetYaxis()->SetTitle("Events");
  hs.GetXaxis()->SetLabelSize(0);
  hs.SetMaximum(std::max(hData->GetMaximum(), hs.GetMaximum()) * 1.8);
  hs.SetMinimum(0.0);

  hData->SetMarkerStyle(20);
  hData->SetLineColor(kBlack);
  hData->Draw("E1 SAME");

  TLegend leg(0.60, 0.55, 0.88, 0.88);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.035);
  if (hTT)  leg.AddEntry(hTT, "t#bar{t} + Jets", "f");
  leg.AddEntry(hDY, "DY + Jets", "f");
  if (hEWK) leg.AddEntry(hEWK,"EWK + Jets","f");
  if (hVV)  leg.AddEntry(hVV, "VV (sum)", "f");
  if (hVVV) leg.AddEntry(hVVV,"VVV (sum)","f");
  leg.AddEntry(hData,"Data","lep");
  leg.Draw();

  TLatex latex; latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.04);
  latex.DrawLatex(0.13, 0.93, "CMS #scale[0.8]{#it{Preliminary}}");
  latex.DrawLatex(0.13, 0.88, (catDir + " (Z peak 70<m_{#mu#mu}<110)").c_str());

  p2.cd();
  TH1* ratio = (TH1*)hData->Clone(("ratio_" + catDir).c_str());
  ratio->SetDirectory(nullptr);
  protectDivide(ratio, hMC);
  ratio->Divide(hMC);

  ratio->SetMinimum(0.0);
  ratio->SetMaximum(2.0);
  ratio->GetYaxis()->SetTitle("Data/MC");
  ratio->GetYaxis()->SetTitleSize(0.12);
  ratio->GetYaxis()->SetTitleOffset(0.40);
  ratio->GetYaxis()->SetLabelSize(0.10);
  ratio->GetYaxis()->SetNdivisions(505);

  ratio->GetXaxis()->SetTitle("p_{T}(Z) [GeV]");
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetLabelSize(0.12);

  ratio->SetMarkerStyle(20);
  ratio->SetLineColor(kBlack);
  ratio->Draw("E1");

  TLine line(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
  line.SetLineStyle(2); line.SetLineWidth(2);
  line.Draw("SAME");

  // Save plots into outBase
  c.SaveAs((outBase + "/ZPtSF_" + catDir + ".png").c_str());
  c.SaveAs((outBase + "/ZPtSF_" + catDir + ".pdf").c_str());

  // Write SF hists into output ROOT file
  fout->mkdir(catDir.c_str());
  fout->cd(catDir.c_str());
  sfDY->Write("SF_DataOverDY");
  sfMC->Write("SF_DataOverMC");

  // cleanup
  delete ratio;
  delete hData; delete hDY; delete hEWK; delete hTT; delete hMC;
  delete sfDY; delete sfMC;
}

// --------------------------- User-facing entry point -------------------
void ZPt_Overlay_Plots(const char* baseDir, const char* /*unusedTag*/, int rebin = 1)
{
  std::string base(baseDir);

  // Derive output tags from baseDir
  std::string year, muCorrTag, cfgTag;
  if (!deriveTagsFromInputPath(base, year, muCorrTag, cfgTag)) {
    std::cerr << "[WARN] Could not parse tags from path (expected .../Analyzed_SkimmedFiles/<YEAR>/<MuCorr>/<CfgTag>)."
              << " Falling back to Unknown tags.\n";
    year = "UnknownYear";
    muCorrTag = "UnknownMuCorr";
    cfgTag = "UnknownCfg";
  }

  // Final output location
  std::string outBase = "/eos/user/d/dblend/StackedPlots/" + year
                      + "/ZPt_Corrections/" + muCorrTag
                      + "/" + cfgTag;

  gSystem->Exec(("mkdir -p " + outBase).c_str());

  // Required: DATA/DY folder each contains one *_all.root
  std::string fDataPath = findAllRootInFolder(base + "/DATA");
  std::string fDYPath   = findAllRootInFolder(base + "/DY");
  std::string fEWKPath  = findAllRootInFolder(base + "/EWK");
  std::string fTTPath   = findAllRootInFolder(base + "/TT");

  if (fDataPath.empty()) std::cerr << "[ERROR] No *_all.root in " << base << "/DATA\n";
  if (fDYPath.empty())   std::cerr << "[ERROR] No *_all.root in " << base << "/DY\n";

  TFile* fData = fDataPath.empty() ? nullptr : openFile(fDataPath);
  TFile* fDY   = fDYPath.empty()   ? nullptr : openFile(fDYPath);
  TFile* fEWK  = fEWKPath.empty()  ? nullptr : openFile(fEWKPath);
  TFile* fTT   = fTTPath.empty()   ? nullptr : openFile(fTTPath);

  if (!fData || !fDY) {
    std::cerr << "[FATAL] Need DATA/*_all.root and DY/*_all.root under " << base << "\n";
    return;
  }

  // VV and VVV groups based on your folder names
  const std::vector<std::string> vvFolders  = {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"};
  const std::vector<std::string> vvvFolders = {"WWW","WWZ","WZZ","ZZZ"};

  // Output SF ROOT file in outBase
  std::string outSF = outBase + "/ZPtSF_" + year + "_" + muCorrTag + "_" + cfgTag + ".root";
  TFile* fout = TFile::Open(outSF.c_str(), "RECREATE");

  const std::vector<std::string> cats = {
    "corr_DM_Pre-Selection_Z0j",
    "corr_DM_Pre-Selection_Z1j",
    "corr_DM_Pre-Selection_Z2j"
  };
  const std::string histName = "DiMuonpt";

  for (const auto& cat : cats) {
    const std::string histPath = cat + "/" + histName;

    // Build grouped VV/VVV hists for THIS category
    TH1* hVV  = sumGroupHists(base, vvFolders,  histPath, rebin);
    TH1* hVVV = sumGroupHists(base, vvvFolders, histPath, rebin);

    makeOneCategory(outBase, cat, histName, fData, fDY, fEWK, fTT, hVV, hVVV, fout, rebin);

    delete hVV;
    delete hVVV;
  }

  fout->Close();
  std::cout << "[OK] Wrote SF ROOT file: " << outSF << "\n";
  std::cout << "[OK] Wrote plots under: " << outBase << "\n";

  // Close inputs
  fData->Close(); fDY->Close();
  if (fEWK) fEWK->Close();
  if (fTT)  fTT->Close();
}

