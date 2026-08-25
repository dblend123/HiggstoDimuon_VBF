// ======================== BTagEff_Maps.C =========================
//
// Builds MC b-tagging efficiency maps from the histograms already
// stored in your analyzer output under:
//
//   BTagEff_Pre-Selection/
//     h_den_b_L, h_num_b_L, h_den_b_M, h_num_b_M
//     h_den_c_L, h_num_c_L, h_den_c_M, h_num_c_M
//     h_den_l_L, h_num_l_L, h_den_l_M, h_num_l_M
//
// The macro sums selected MC samples and computes
//   eff = numerator / denominator
// separately for flavor = b,c,light and WP = L,M.
//
// Output:
//   /eos/user/d/dblend/StackedPlots/<YEAR>/BTagEfficiency/<MuCorr>/<CfgTag>/
//
// Usage:
//   root -l -q 'BTagEff_Maps.C("/eos/user/d/dblend/Analyzed_SkimmedFiles/22022/KIT/noScaleSmear_noBtagSF_noZPtRW")'
//
// =====================================================================

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TString.h"

#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>

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

// -------- Find the first "*_all.root" file inside a folder -------------
static std::string findAllRootInFolder(const std::string& folder) {
  void* dirp = gSystem->OpenDirectory(folder.c_str());
  if (!dirp) return "";

  const char* entry;
  while ((entry = gSystem->GetDirEntry(dirp))) {
    std::string name = entry;
    if (name == "." || name == "..") continue;
    if (name.size() >= 9 && name.substr(name.size()-9) == "_all.root") {
      gSystem->FreeDirectory(dirp);
      return folder + "/" + name;
    }
  }
  gSystem->FreeDirectory(dirp);
  return "";
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

// ------------------------- clone+detach TH2 ----------------------------
static TH2* getHist2(TFile* f, const std::string& path) {
  if (!f) return nullptr;
  auto* h = dynamic_cast<TH2*>(f->Get(path.c_str()));
  if (!h) {
    std::cerr << "[WARN] Missing TH2: " << path << "\n";
    return nullptr;
  }
  auto* c = dynamic_cast<TH2*>(h->Clone());
  c->SetDirectory(nullptr);
  c->Sumw2();
  return c;
}

// ----- Sum histograms over a group of folders (each folder has *_all.root)
static TH2* sumGroupHists2(const std::string& baseDir,
                           const std::vector<std::string>& folders,
                           const std::string& histPath)
{
  TH2* hsum = nullptr;

  for (const auto& sub : folders) {
    const std::string folder = baseDir + "/" + sub;
    const std::string fpath  = findAllRootInFolder(folder);
    if (fpath.empty()) {
      std::cerr << "[WARN] No *_all.root found in " << folder << "\n";
      continue;
    }

    TFile* f = openFile(fpath);
    if (!f) continue;

    TH2* h = getHist2(f, histPath);
    f->Close();
    delete f;

    if (!h) continue;

    if (!hsum) {
      hsum = dynamic_cast<TH2*>(h->Clone(("SUM_" + sub + "_" + histPath).c_str()));
      hsum->SetDirectory(nullptr);
      hsum->Sumw2();
    } else {
      hsum->Add(h);
    }
    delete h;
  }

  return hsum;
}

// ------------------------- safe divide for efficiency -------------------
static TH2* makeEfficiencyMap(TH2* hNum, TH2* hDen, const std::string& name)
{
  if (!hNum || !hDen) return nullptr;

  TH2* hEff = dynamic_cast<TH2*>(hNum->Clone(name.c_str()));
  hEff->SetDirectory(nullptr);
  hEff->Reset();

  for (int ix = 1; ix <= hDen->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hDen->GetNbinsY(); ++iy) {
      const double den = hDen->GetBinContent(ix, iy);
      const double num = hNum->GetBinContent(ix, iy);

        if (den > 0.0) {
            const double eff =
                std::clamp(
                    num / den,
                    0.0,
                    1.0
                );            
            const double sumW   = hDen->GetBinContent(ix, iy);
            const double sumW2  = std::pow(hDen->GetBinError(ix, iy), 2);
            double nEff = 0.0;

            if (sumW2 > 0.0) {
                nEff =
                    (sumW * sumW) /
                    sumW2;
            }

            double err = 0.0;

            if (nEff > 0.0) {
                err =
                    std::sqrt(
                        std::max(
                            0.0,
                            eff * (1.0 - eff) /
                            nEff
                        )
                    );
            }

            hEff->SetBinContent(
                ix,
                iy,
                eff
            );

            hEff->SetBinError(
                ix,
                iy,
                err
            );
        }
        else {

            hEff->SetBinContent(
                ix,
                iy,
                0.0
            );

            hEff->SetBinError(
                ix,
                iy,
                0.0
            );
        }
     }
  }

  return hEff;
}

// ------------------------- draw one map -------------------
static void drawMap(TH2* h, const std::string& outBase, const std::string& tag,
                    const std::string& title, const std::string& era,
                    const std::string& muCorr, const std::string& cfg)
{
  if (!h) return;

  TCanvas c(("c_" + tag).c_str(), "", 900, 750);
  gStyle->SetOptStat(0);
  c.SetRightMargin(0.15);
  c.SetLeftMargin(0.12);
  c.SetBottomMargin(0.12);

  h->SetTitle("");
  h->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h->GetYaxis()->SetTitle("|#eta|");
  h->GetZaxis()->SetTitle("Efficiency");
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetZaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.040);
  h->GetYaxis()->SetLabelSize(0.040);
  h->GetZaxis()->SetLabelSize(0.040);
  h->SetMinimum(0.0);
  h->SetMaximum(1.0);
  h->Draw("COLZ TEXT");

  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(0.040);
  latex.DrawLatex(0.12, 0.94, "CMS #scale[0.8]{#it{Preliminary}}");
  latex.DrawLatex(0.12, 0.89, title.c_str());
  latex.DrawLatex(0.12, 0.84, Form("Era=%s, MuCorr=%s", era.c_str(), muCorr.c_str()));
  latex.DrawLatex(0.12, 0.79, cfg.c_str());

  c.SaveAs((outBase + "/" + tag + ".png").c_str());
  c.SaveAs((outBase + "/" + tag + ".pdf").c_str());
}

// ------------------------- process one flavor/WP pair -------------------
static void processOne(const std::string& baseDir,
                       const std::vector<std::string>& mcFolders,
                       const std::string& dirName,
                       const std::string& denName,
                       const std::string& numName,
                       const std::string& outBase,
                       TFile* fout,
                       const std::string& era,
                       const std::string& muCorr,
                       const std::string& cfg,
                       const std::string& flavorTag,
                       const std::string& wpTag)
{
  const std::string denPath = dirName + "/" + denName;
  const std::string numPath = dirName + "/" + numName;

  TH2* hDen = sumGroupHists2(baseDir, mcFolders, denPath);
  TH2* hNum = sumGroupHists2(baseDir, mcFolders, numPath);

  if (!hDen || !hNum) {
    std::cerr << "[WARN] Missing numerator/denominator for " << flavorTag << " " << wpTag << "\n";
    delete hDen;
    delete hNum;
    return;
  }

  TH2* hEff = makeEfficiencyMap(hNum, hDen, "eff_" + flavorTag + "_" + wpTag);
  if (!hEff) {
    delete hDen;
    delete hNum;
    return;
  }

  const std::string title = "BTag efficiency: flavor=" + flavorTag + ", WP=" + wpTag;
  const std::string stem  = "BTagEff_" + flavorTag + "_" + wpTag;

  drawMap(hEff, outBase, stem, title, era, muCorr, cfg);

  fout->cd();
  hDen->Write(("den_" + flavorTag + "_" + wpTag).c_str());
  hNum->Write(("num_" + flavorTag + "_" + wpTag).c_str());
  hEff->Write(("eff_" + flavorTag + "_" + wpTag).c_str());

  std::cout << "[OK] Built " << stem << "\n";

  delete hDen;
  delete hNum;
  delete hEff;
}

// --------------------------- User-facing entry point -------------------
void BTagEff_Maps(const char* baseDir)
{
  std::string base(baseDir);
  if (!base.empty() && base.back() == '/') base.pop_back();

  std::string year, muCorrTag, cfgTag;
  if (!deriveTagsFromInputPath(base, year, muCorrTag, cfgTag)) {
    std::cerr << "[WARN] Could not parse tags from input path. Falling back to Unknown tags.\n";
    year = "UnknownYear";
    muCorrTag = "UnknownMuCorr";
    cfgTag = "UnknownCfg";
  }

  const std::string outBase = "/eos/user/d/dblend/StackedPlots/" + year
                            + "/BTagEfficiency/" + muCorrTag
                            + "/" + cfgTag;

  gSystem->Exec(("mkdir -p " + outBase).c_str());

  const std::string outRoot = outBase + "/BTagEffMaps_" + year + "_" + muCorrTag + "_" + cfgTag + ".root";
  TFile* fout = TFile::Open(outRoot.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "[FATAL] Could not create output file: " << outRoot << "\n";
    return;
  }

  // MC samples to include in the efficiency maps
  // Adjust this grouping if you want process-specific maps later.
  const std::vector<std::string> mcFolders = {
    "DY", "EWK", "TT",
    "WZto2L2Q", "WZto3LNu",
    "ZZto2L2Nu", "ZZto2L2Q", "ZZto4L",
    "WWW", "WWZ", "WZZ", "ZZZ",
    "VBF"
  };

  const std::string dirName = "BTagEff_Pre-Selection";

  processOne(base, mcFolders, dirName, "h_den_b_L", "h_num_b_L", outBase, fout, year, muCorrTag, cfgTag, "b", "L");
  processOne(base, mcFolders, dirName, "h_den_b_M", "h_num_b_M", outBase, fout, year, muCorrTag, cfgTag, "b", "M");

  processOne(base, mcFolders, dirName, "h_den_c_L", "h_num_c_L", outBase, fout, year, muCorrTag, cfgTag, "c", "L");
  processOne(base, mcFolders, dirName, "h_den_c_M", "h_num_c_M", outBase, fout, year, muCorrTag, cfgTag, "c", "M");

  processOne(base, mcFolders, dirName, "h_den_l_L", "h_num_l_L", outBase, fout, year, muCorrTag, cfgTag, "light", "L");
  processOne(base, mcFolders, dirName, "h_den_l_M", "h_num_l_M", outBase, fout, year, muCorrTag, cfgTag, "light", "M");

  fout->Close();
  delete fout;

  std::cout << "[OK] Wrote efficiency ROOT file: " << outRoot << "\n";
  std::cout << "[OK] Wrote plots under: " << outBase << "\n";
}
