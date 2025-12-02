// Muon_MuCorrComparison.C  (KIT vs NoMuCorr, stamped + Z-fit version, NO THStack)
// root -l -b -q 'Muon_MuCorrComparison.C("22023", "/eos/user/d/dblend/Analyzed_SkimmedFiles/22023/", "noPU_noScaleSmear_noMuonSF_noJERC_noJetVeto_noGolden_noFlags_noFSR_noBtagSF_noJetCriteria")'

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TPad.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TNamed.h>

// ---- global config for subdirectory (e.g. "noPU", "noPU_noScaleSmear_...") ----
static TString gConfig = "";

// ---- helpers ----
static TString modeDir(const TString& mode) {
  // Map mode ("KIT", "NoMuCorr") to directory fragment under base,
  // optionally including the CONFIG subdirectory.
  TString m = mode;
  m.ToLower();
  if (m == "kit") {
    if (gConfig.IsNull() || gConfig == "") return "KIT";
    return "KIT/" + gConfig;
  }
  if (m == "nomucorr") {
    if (gConfig.IsNull() || gConfig == "") return "NoMuCorr";
    return "NoMuCorr/" + gConfig;
  }
  // default (not used anymore)
  return "Rochester";
}

static double eraLumiPb(const TString& era) {
  if (era == "12022") return 7980.454145;
  if (era == "22022") return 26671.60971;
  if (era == "12023") return 18704.13341;
  if (era == "22023") return 9693.130053;
  if (era == "2024")  return 0.0; // put real 2024 lumi here if you want
  std::cerr << "[WARN] Unknown era " << era << ", lumi set to 0\n";
  return 0.0;
}

static TH1* getHist(TFile* f, const std::string& folder, const std::string& name) {
  if (!f) return nullptr;
  TString path = folder + name;
  TH1* h = dynamic_cast<TH1*>(f->Get(path));
  if (!h) {
    TString alt = (folder.size() && folder.back()=='/')
                    ? TString(folder.c_str()) + name.c_str()
                    : TString(folder.c_str()) + "/" + name.c_str();
    h = dynamic_cast<TH1*>(f->Get(alt));
  }
  return h;
}

static TString readStamp(TFile* f, const char* key) {
  if (!f) return "";
  if (auto* n = dynamic_cast<TNamed*>(f->Get(key))) return TString(n->GetTitle());
  return "";
}

static std::map<std::string,TFile*> openModeFiles(const TString& base, const TString& mode,
                                                  const std::vector<TString>& samples) {
  std::map<std::string,TFile*> m;
  const TString sub = modeDir(mode) + "/";
  for (const auto& s : samples) {
    TString path = base + sub + s + "/" + s + "_all.root";
    m[s.Data()] = TFile::Open(path);
    if (!m[s.Data()]) std::cerr << "[WARN] Missing file: " << path << "\n";
  }
  return m;
}

static void closeModeFiles(std::map<std::string,TFile*>& m) {
  for (auto& kv : m) if (kv.second) kv.second->Close();
  m.clear();
}

static TH1* sumGroup(const std::map<std::string,TFile*>& mp,
                     const std::vector<std::string>& names,
                     const std::string& folder,
                     const std::string& hname,
                     const TString& cloneName) {
  TH1* sum = nullptr;
  for (const auto& n : names) {
    auto it = mp.find(n);
    if (it==mp.end() || !it->second) continue;
    TH1* h = getHist(it->second, folder, hname);
    if (!h) continue;
    if (!sum) { sum = (TH1*)h->Clone(cloneName); sum->SetDirectory(nullptr); }
    else      { sum->Add(h); }
  }
  return sum;
}

static TH1* buildTotalMC(const std::map<std::string,TFile*>& fmap,
                         const std::string& folder, const std::string& var,
                         const TString& name) {
  TH1* h_DY  = sumGroup(fmap, {"DY"}, folder, var, name+"_DY");
  TH1* h_EWK = sumGroup(fmap, {"EWK"}, folder, var, name+"_EWK");
  TH1* h_TT  = sumGroup(fmap, {"TT"}, folder, var, name+"_TT");
  TH1* h_VV  = sumGroup(fmap, {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"}, folder, var, name+"_VV");
  TH1* h_VVV = sumGroup(fmap, {"WWW","WWZ","WZZ","ZZZ"}, folder, var, name+"_VVV");
  TH1* tot = nullptr;
  auto add = [&](TH1* h){
    if (!h) return;
    if (!tot) { tot = (TH1*)h->Clone(name); tot->SetDirectory(nullptr); }
    else tot->Add(h);
  };
  add(h_DY); add(h_EWK); add(h_TT); add(h_VV); add(h_VVV);
  if (h_DY)  delete h_DY;
  if (h_EWK) delete h_EWK;
  if (h_TT)  delete h_TT;
  if (h_VV)  delete h_VV;
  if (h_VVV) delete h_VVV;
  return tot;
}

static double integralSafe(TH1* h) { return h ? h->Integral(0, h->GetNbinsX()+1) : 0.0; }

// quick Gaussian fit in Z window
static bool fitZ(TH1* hin, double& mean, double& sigma) {
  if (!hin || hin->GetEntries() == 0) return false;
  auto* h = (TH1*)hin->Clone("fitZ_clone");
  h->SetDirectory(nullptr);
  const double x1 = 80.0, x2 = 100.0;
  int b1 = h->GetXaxis()->FindBin(x1), b2 = h->GetXaxis()->FindBin(x2);
  if (b2 <= b1) { delete h; return false; }
  h->GetXaxis()->SetRange(b1, b2);
  TF1 g("g","gaus", x1, x2);
  g.SetParameters(h->GetMaximum(), 91.0, 2.0);
  int status = h->Fit(&g, "QNR"); // quiet, no draw, use range
  if (status != 0) { delete h; return false; }
  mean  = g.GetParameter(1);
  sigma = g.GetParameter(2);
  delete h;
  return true;
}

static TFile* openDataForMode(const TString& base, const TString& mode) {
  const TString md = modeDir(mode); // "KIT[/CONFIG]" | "NoMuCorr[/CONFIG]"
  const TString path = base + md + "/DATA/DATA_all.root";
  TFile* f = TFile::Open(path);
  if (!f) {
    std::cerr << "[WARN] Missing per-mode Data: " << path
              << "  (falling back to legacy Data/Data_all.root)\n";
    f = TFile::Open(base + "Data/Data_all.root");
  }
  return f;
}

// ---- main ----
void Muon_Overlay_MuCorrComparison(TString era = "2024",
                           TString basePath = "",
                           TString config = "")
{
  gStyle->SetOptStat(0);

  // Set global config for directory structure
  gConfig = config;

  // Start from user-provided base path (non-const so we can tweak it)
  TString base = basePath;
  if (base.IsNull() || base == "") {
    std::cerr << "[INFO] No basePath provided; using default 22023 directory.\n";
    base = "/eos/user/d/dblend/Analyzed_SkimmedFiles/22023/";
  }

  // Ensure base ends with a slash
  if (!base.EndsWith("/")) base += "/";

  const double lumi = eraLumiPb(era);

  std::cerr << "[INFO] Using base path: " << base << "\n";
  if (!gConfig.IsNull() && gConfig != "") {
    std::cerr << "       → Using CONFIG subdir: " << gConfig << "\n";
  }

  // Per-mode Data (fallback to legacy flat Data if missing)
  TFile* f_dataKIT = openDataForMode(base, "KIT");
  TFile* f_dataNOM = openDataForMode(base, "NoMuCorr");

  if (!f_dataKIT && !f_dataNOM) {
    std::cerr << "[ERR] No Data files found for KIT or NoMuCorr (or legacy path).\n";
    return;
  }

  const std::vector<TString> samples = {
    "DY","EWK", "TT","WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L",
    "WWW","WWZ","WZZ","ZZZ",
    "VBF"
  };

  // Only DM preselection regions
  const std::string preDM = "corr_DM_Pre-Selection";
  const std::vector<std::string> regions = {"", "_BB", "_BO", "_BE", "_OO", "_OE", "_EE"};
  std::vector<std::string> folders;
  for (const auto& s : regions) folders.push_back(preDM + s + "/");

  const std::string var = "DiMuonMass_Z";
  const char* xlab = "m_{#mu#mu} [GeV]";

  // Compare only KIT vs NoMuCorr
  const std::vector<TString> modes = {"KIT", "NoMuCorr"};

  // open per-mode file maps
  auto fmapKIT  = openModeFiles(base, "KIT",      samples);
  auto fmapNOM  = openModeFiles(base, "NoMuCorr", samples);

  // show stamps (DY representative)
  std::cerr << "[STAMP] DY KIT="
            << readStamp(fmapKIT.count("DY")? fmapKIT.at("DY"):nullptr,  "muCorrMode")
            << " NoMuCorr="
            << readStamp(fmapNOM.count("DY")? fmapNOM.at("DY"):nullptr,  "muCorrMode")
            << "\n";

  // Loop over regions/folders
  for (const auto& folder : folders) {
    // totals for combined ratio
    TH1* totKIT  = buildTotalMC(fmapKIT,  folder, var, "totKIT");
    TH1* totNOM  = buildTotalMC(fmapNOM,  folder, var, "totNOM");

    // --- check if this folder actually has anything (MC or Data) ---
    bool hasAnyData = false;
    if (f_dataKIT  && getHist(f_dataKIT,  folder, var)) hasAnyData = true;
    if (f_dataNOM  && getHist(f_dataNOM,  folder, var)) hasAnyData = true;

    if (!totKIT && !totNOM && !hasAnyData) {
      // nothing to plot in this folder (e.g. BB, BO, OO, ... that no longer exist)
      if (totKIT) delete totKIT;
      if (totNOM) delete totNOM;
      continue;  // skip to next folder
    }

    // ---------- Per-mode “stacked-style” plots (no THStack) ----------
    for (const auto& MODE : modes) {

      const std::map<std::string,TFile*>& fmap =
        (MODE=="KIT") ? fmapKIT : fmapNOM;

      TFile* f_dataMODE = (MODE=="KIT") ? f_dataKIT : f_dataNOM;

      TH1* hD = getHist(f_dataMODE, folder, var);
      TH1* h_data = nullptr;
      if (hD) {
        h_data = (TH1*)hD->Clone(Form("h_data_%s", MODE.Data()));
        h_data->SetDirectory(nullptr);
      }

      TH1* h_DY   = sumGroup(fmap, {"DY"}, folder, var, "h_DY");
      TH1* h_EWK  = sumGroup(fmap, {"EWK"}, folder, var, "h_EWK");
      TH1* h_TT   = sumGroup(fmap, {"TT"}, folder, var, "h_TT");
      TH1* h_VV   = sumGroup(fmap, {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"}, folder, var, "h_VV");
      TH1* h_VVV  = sumGroup(fmap, {"WWW","WWZ","WZZ","ZZZ"}, folder, var, "h_VVV");
      TH1* h_SIG  = sumGroup(fmap, {"VBF"}, folder, var, "h_SIG");

      // if this MODE has no MC and no Data and no signal, skip drawing
      bool hasMC = (h_DY || h_EWK || h_TT || h_VV || h_VVV);
      if (!hasMC && !h_data && !h_SIG) {
        if (h_DY)  delete h_DY;
        if (h_EWK) delete h_EWK;
        if (h_TT)  delete h_TT;
        if (h_VV)  delete h_VV;
        if (h_VVV) delete h_VVV;
        if (h_SIG) delete h_SIG;
        continue;
      }

      // style
      if (h_SIG) { h_SIG->SetLineColor(kMagenta+2); h_SIG->SetLineWidth(2); h_SIG->SetFillStyle(0); }
      if (h_TT)  { h_TT ->SetFillColor(kRed-4);    h_TT ->SetLineColor(kRed-4);    }
      if (h_DY)  { h_DY ->SetFillColor(kOrange-3); h_DY ->SetLineColor(kOrange-3); }
      if (h_EWK) { h_EWK->SetFillColor(kTeal-6);   h_EWK->SetLineColor(kTeal+2);   }
      if (h_VV)  { h_VV ->SetFillColor(kAzure-9);  h_VV ->SetLineColor(kAzure-9);  }
      if (h_VVV) { h_VVV->SetFillColor(kSpring+4); h_VVV->SetLineColor(kSpring+4); }

      if (h_data) {
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(0.8);
        h_data->SetLineColor(kBlack);
      }

      // choose total MC for this mode (for ratio, ymax decision)
      TH1* totMode = (MODE=="KIT") ? totKIT : totNOM;

      TCanvas c("c","c",900,800);
      TPad pad1("pad1","pad1",0,0.30,1,1);
      TPad pad2("pad2","pad2",0,0.00,1,0.30);
      pad1.SetBottomMargin(0.02);
      pad2.SetTopMargin(0.03);
      pad2.SetBottomMargin(0.35);
      pad1.Draw();
      pad2.Draw();

      // --- Top pad ---
      pad1.cd();
      double ymax = 1.0;
      auto updMax = [&](TH1* h){
        if (!h) return;
        ymax = std::max(ymax, 1.5*h->GetMaximum());
      };
      updMax(totMode);
      updMax(h_data);
      updMax(h_DY);
      updMax(h_EWK);
      updMax(h_TT);
      updMax(h_VV);
      updMax(h_VVV);
      updMax(h_SIG);
      if (ymax <= 0) ymax = 1.0;

      bool useLogY = ((totMode && totMode->Integral() > 0) ||
                      (h_data && h_data->Integral() > 0));
      pad1.SetLogy(useLogY);
      double ymin = useLogY ? 1e-1 : 0.0;

      TH1* frame = nullptr;
      // draw a frame axis first
      if (totMode || h_data || h_SIG || h_DY || h_EWK || h_TT || h_VV || h_VVV) {
        TH1* ref = totMode ? totMode :
                   (h_data ? h_data :
                   (h_DY ? h_DY :
                   (h_EWK ? h_EWK :
                   (h_TT ? h_TT :
                   (h_VV ? h_VV :
                   (h_VVV ? h_VVV : h_SIG))))));
        frame = (TH1*)ref->Clone("frame");
        frame->Reset();
        frame->SetDirectory(nullptr);
        frame->GetYaxis()->SetTitle("Events");
        frame->GetXaxis()->SetTitle("");
        frame->GetYaxis()->SetRangeUser(ymin, ymax);
        frame->Draw("axis");
      }

      // "stacked-style" draw order: VVV, VV, TT, DY, EWK (just drawn one on top of the other)
      if (h_VVV) h_VVV->Draw("hist same");
      if (h_VV)  h_VV ->Draw("hist same");
      if (h_TT)  h_TT ->Draw("hist same");
      if (h_DY)  h_DY ->Draw("hist same");
      if (h_EWK) h_EWK->Draw("hist same");

      if (h_SIG)   h_SIG->Draw("hist same");
      if (h_data)  h_data->Draw("e1 same");

      // Legend
      TLegend leg(0.63,0.60,0.88,0.88);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      if (h_SIG)   leg.AddEntry(h_SIG,  "VBF H#rightarrow#mu#mu", "l");
      if (h_TT)    leg.AddEntry(h_TT,   "t#bar{t} + Jets", "f");
      if (h_DY)    leg.AddEntry(h_DY,   "DY + Jets", "f");
      if (h_EWK)   leg.AddEntry(h_EWK,  "EWK + Jets", "f");
      if (h_VV)    leg.AddEntry(h_VV,   "VV + Jets", "f");
      if (h_VVV)   leg.AddEntry(h_VVV,  "VVV + Jets", "f");
      if (h_data)  leg.AddEntry(h_data, "Data", "lep");
      leg.Draw();

      // CMS-like label
      TPaveText pt(0.12,0.92,0.95,0.99,"NDC");
      pt.SetFillStyle(0);
      pt.SetBorderSize(0);
      pt.AddText(Form("CMS Preliminary   #sqrt{s}=13.6 TeV,  L=%.1f fb^{-1}   Era %s",
                      lumi/1000.0, era.Data()));
      pt.Draw();

      // --- Bottom pad: Data/MC ratio ---
      pad2.cd();
      TH1* h_ratio = nullptr;
      if (h_data && totMode) {
        h_ratio = (TH1*)h_data->Clone("h_ratio");
        h_ratio->SetDirectory(nullptr);
        h_ratio->Divide(totMode);
        h_ratio->GetYaxis()->SetTitle("Data/MC");
        h_ratio->GetXaxis()->SetTitle(xlab);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.11);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.11);
        h_ratio->GetYaxis()->SetRangeUser(0.92, 1.15);
        h_ratio->Draw("e1");
        // dashed unity
        TLine l1(h_ratio->GetXaxis()->GetXmin(),1.0,
                 h_ratio->GetXaxis()->GetXmax(),1.0);
        l1.SetLineStyle(2);
        l1.Draw("same");
        // connecting line through points
        TGraphErrors gr(h_ratio);
        gr.SetLineColor(kGray+2);
        gr.SetLineWidth(2);
        gr.SetMarkerStyle(20);
        gr.SetMarkerSize(0);
        gr.Draw("L same");
      } else {
        TH1* ax = nullptr;
        if (h_data)      ax = (TH1*)h_data->Clone("ax");
        else if (totMode)  ax = (TH1*)totMode->Clone("ax");
        if (ax) {
          ax->SetDirectory(nullptr);
          ax->Reset();
          ax->GetYaxis()->SetTitle("Data/MC");
          ax->GetXaxis()->SetTitle(xlab);
          ax->Draw();
          delete ax;
        }
      }

      // Save
      TString clean = TString(folder).ReplaceAll("/","").ReplaceAll("corr_","");
      TString cleanTag = clean;
      if (!gConfig.IsNull() && gConfig != "") {
        TString cfgSafe = gConfig;
        cfgSafe.ReplaceAll("/","_");
        cleanTag = cfgSafe + "_" + clean;
      }

      TString outdirMode = "/eos/user/d/dblend/StackedPlots/" + era + "/MuCorr_Comparison/" + MODE + "/";
      gSystem->Exec("mkdir -p " + outdirMode);
      TString out = outdirMode + "StackWithRatio_" + cleanTag + ".pdf";
      c.SaveAs(out);
      c.SaveAs(out.ReplaceAll(".pdf",".png"));

      // cleanup for this MODE
      if (frame)  { delete frame;  frame  = nullptr; }
      if (h_ratio){ delete h_ratio; h_ratio = nullptr; }
      if (h_DY)   delete h_DY;
      if (h_EWK)  delete h_EWK;
      if (h_TT)   delete h_TT;
      if (h_VV)   delete h_VV;
      if (h_VVV)  delete h_VVV;
      if (h_SIG)  delete h_SIG;
      if (h_data) { delete h_data; h_data = nullptr; }
    } // per-mode

    // ---------- Combined KIT vs NoMuCorr Data/MC ratio overlay ----------
    TCanvas cR("cR","cR",900,700);
    TPad padR("padR","padR",0,0.00,1,1);
    padR.SetBottomMargin(0.13);
    padR.Draw();
    padR.cd();

    std::vector<std::pair<TString,TH1*>> ratios;

    // Rebuild Data/MC ratio for each mode using that mode's Data
    TH1 *dK=nullptr, *dN=nullptr;

    if (f_dataKIT) {
      if (auto* h = getHist(f_dataKIT, folder, var)) {
        dK = (TH1*)h->Clone("d_kit");
        dK->SetDirectory(nullptr);
      }
    }
    if (f_dataNOM) {
      if (auto* h = getHist(f_dataNOM, folder, var)) {
        dN = (TH1*)h->Clone("d_nom");
        dN->SetDirectory(nullptr);
      }
    }

    if (dK && totKIT)  { dK->Divide(totKIT);  ratios.push_back({"KIT",      dK}); }
    if (dN && totNOM)  { dN->Divide(totNOM);  ratios.push_back({"NoMuCorr", dN}); }

    double yminR = 0.0, ymaxR = 0.2;
    bool firstDrawn = false;
    TLegend legR(0.50,0.75,0.88,0.90);
    legR.SetBorderSize(0);
    legR.SetFillStyle(0);

    for (auto& pr : ratios) {
      TH1* h = pr.second;
      if (!h) continue;
      int col = kBlack;
      if (pr.first=="KIT")      col = kGreen+2;
      if (pr.first=="NoMuCorr") col = kBlue+1;
      h->SetLineColor(col);
      h->SetMarkerColor(col);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(0.7);
      h->GetYaxis()->SetTitle("Data/MC");
      h->GetXaxis()->SetTitle(xlab);
      h->GetYaxis()->SetRangeUser(yminR, ymaxR);
      if (!firstDrawn) {
        h->Draw("e1");
        firstDrawn = true;
      } else {
        h->Draw("e1 same");
      }
      // connecting line
      TGraphErrors gr(h);
      gr.SetLineColor(col);
      gr.SetLineWidth(2);
      gr.SetMarkerSize(0);
      gr.Draw("L same");
      legR.AddEntry(h, pr.first, "lep");
    }

    if (firstDrawn) {
      TH1* href = ratios.front().second;
      TLine l1(href->GetXaxis()->GetXmin(),1.0, href->GetXaxis()->GetXmax(),1.0);
      l1.SetLineStyle(2);
      l1.Draw("same");
    }

    // Build a Data histogram for the Z-fit (use KIT/Data or fallback)
    TH1* hDataR = nullptr;
    {
      TString dataPath = base + modeDir("KIT") + "/DATA/DATA_all.root";
      std::unique_ptr<TFile> fK(TFile::Open(dataPath));
      if (!fK) {
        std::cerr << "[WARN] Missing " << dataPath
                  << " ; falling back to " << (base + "Data/Data_all.root") << "\n";
        fK.reset(TFile::Open(base + "Data/Data_all.root"));
      }
      if (fK) {
        if (auto* h = getHist(fK.get(), folder, var)) {
          hDataR = (TH1*)h->Clone("hDataR");
          hDataR->SetDirectory(nullptr);
        }
      }
    }

    // ---- Z-fit means (vertical lines) ----
    double mData=0, sData=0, mK=0, sK=0, mN=0, sN=0;
    bool okD = fitZ(hDataR,  mData, sData);
    bool okK = fitZ(totKIT,  mK,    sK);
    bool okN = fitZ(totNOM,  mN,    sN);

    if (firstDrawn && (okD||okK||okN)) {
      TH1* href = ratios.front().second;
      double ylow = yminR, yhigh = ymaxR;
      auto drawV = [&](double x, int col, int style){
        if (x<=0 || !href) return;
        TLine ln(x, ylow, x, yhigh);
        ln.SetLineColor(col);
        ln.SetLineStyle(style);
        ln.SetLineWidth(2);
        ln.Draw("same");
      };
      if (okD) drawV(mData, kBlack,     1);
      if (okK) drawV(mK,    kGreen+2,   7);
      if (okN) drawV(mN,    kBlue+1,    9);
    }

    legR.Draw();
    TPaveText ptR(0.12,0.92,0.95,0.99,"NDC");
    ptR.SetFillStyle(0);
    ptR.SetBorderSize(0);
    ptR.AddText(Form("CMS Preliminary   #sqrt{s}=13.6 TeV,  L=%.1f fb^{-1}   Era %s",
                    lumi/1000.0, era.Data()));
    ptR.Draw();

    TString clean2 = TString(folder).ReplaceAll("/","").ReplaceAll("corr_","");
    TString cleanTag2 = clean2;
    if (!gConfig.IsNull() && gConfig != "") {
      TString cfgSafe = gConfig;
      cfgSafe.ReplaceAll("/","_");
      cleanTag2 = cfgSafe + "_" + clean2;
    }

    TString outdirComp = "/eos/user/d/dblend/StackedPlots/" + era + "/MuCorr_Comparison/Combined/";
    gSystem->Exec("mkdir -p " + outdirComp);
    TString outR = outdirComp + "RatioOverlay_" + cleanTag2 + ".pdf";
    cR.SaveAs(outR);
    cR.SaveAs(outR.ReplaceAll(".pdf",".png"));

    // ---------- TXT summary with Z-fit + yields ----------
    TString txtPath = outdirComp + "MuCorr_Summary_" + cleanTag2 + ".txt";
    std::ofstream ofs(txtPath.Data());
    ofs << std::fixed << std::setprecision(3);
    ofs << "# Summary for folder " << folder << " (DiMuonMass_Z)\n";
    double yData = integralSafe(hDataR);
    ofs << "Data: " << yData << "\n";

    auto writeMode = [&](const char* label, const std::map<std::string,TFile*>& fmapX, TH1* totH){
      TH1* c_DY   = sumGroup(fmapX, {"DY"}, folder, var, "c_DY");
      TH1* c_EWK  = sumGroup(fmapX, {"EWK"}, folder, var, "c_EWK");
      TH1* c_TT   = sumGroup(fmapX, {"TT"}, folder, var, "c_TT");
      TH1* c_VV   = sumGroup(fmapX, {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"}, folder, var, "c_VV");
      TH1* c_VVV  = sumGroup(fmapX, {"WWW","WWZ","WZZ","ZZZ"}, folder, var, "c_VVV");
      double yDY   = integralSafe(c_DY);
      double yEWK  = integralSafe(c_EWK);
      double yTT   = integralSafe(c_TT);
      double yVV   = integralSafe(c_VV);
      double yVVV  = integralSafe(c_VVV);
      double yTot  = integralSafe(totH);
      double ratio = (yTot>0 ? yData/yTot : 0.0);
      ofs << "\n[" << label << "]\n";
      ofs << " MC_total: " << yTot << "\n";
      ofs << "   DY:  " << yDY << "\n";
      ofs << "   EWK: " << yEWK << "\n";
      ofs << "   TT:  " << yTT << "\n";
      ofs << "   VV:  " << yVV << "\n";
      ofs << "   VVV: " << yVVV << "\n";
      ofs << " Data/MC: " << ratio << "\n";
      if (c_DY)  delete c_DY;
      if (c_EWK) delete c_EWK;
      if (c_TT)  delete c_TT;
      if (c_VV)  delete c_VV;
      if (c_VVV) delete c_VVV;
    };
    writeMode("KIT",      fmapKIT,  totKIT);
    writeMode("NoMuCorr", fmapNOM,  totNOM);

    ofs << "\n# Z-peak Gaussian fit (80–100 GeV)\n";
    if (okD) ofs << "Data:     mean=" << mData << "  sigma=" << sData << "\n";
    if (okK) ofs << "KIT:      mean=" << mK    << "  sigma=" << sK    << "\n";
    if (okN) ofs << "NoMuCorr: mean=" << mN    << "  sigma=" << sN    << "\n";
    ofs.close();

    // cleanup combined objects
    if (hDataR) { delete hDataR; hDataR = nullptr; }
    if (totKIT) { delete totKIT; totKIT = nullptr; }
    if (totNOM) { delete totNOM; totNOM = nullptr; }
    for (auto& pr : ratios) {
      if (pr.second) delete pr.second;
      pr.second = nullptr;
    }
  } // folders

  closeModeFiles(fmapKIT);
  closeModeFiles(fmapNOM);
  if (f_dataKIT)  f_dataKIT->Close();
  if (f_dataNOM)  f_dataNOM->Close();

  std::cout << "[OK] Wrote KIT vs NoMuCorr plots and summaries for CONFIG = '"
            << gConfig << "'.\n";
}

