// Muon_MuCorrComparison.C  (stamped + Z-fit version)
//
// Usage:
//   root -l -b -q 'Muon_MuCorrComparison.C("12022")'
//
// Inputs it expects:
//   /eos/user/d/dblend/Analyzed_SkimmedFiles/<ERA>/Data/Data_all.root
//   /eos/user/d/dblend/Analyzed_SkimmedFiles/<ERA>/{Rochester,KIT,NoMuCorr}/<Sample>/<Sample>_all.root
//
// Outputs:
//   /eos/user/d/dblend/StackedPlots/<ERA>/MuCorr_Comparison/<MODE>/StackWithRatio_<folder>.{pdf,png}
//   /eos/user/d/dblend/StackedPlots/<ERA>/MuCorr_Comparison/Combined/RatioOverlay_<folder>.{pdf,png}
//   /eos/user/d/dblend/StackedPlots/<ERA>/MuCorr_Comparison/Combined/MuCorr_Summary_<folder>.txt

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
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

// ---- helpers ----
static TString modeDir(const TString& mode) {
  TString m = mode; m.ToLower();
  if (m == "kit")      return "KIT";
  if (m == "nomucorr") return "NoMuCorr";
  return "Rochester";
}
static double eraLumiPb(const TString& era) {
  if (era == "12022") return 7980.454145;
  if (era == "22022") return 26671.60971;
  if (era == "12023") return 18704.13341;
  if (era == "22023") return 9693.130053;
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
  TH1* h_EWK  = sumGroup(fmap, {"EWK"}, folder, var, name+"_EWK");
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
  if (h_DY) delete h_DY;
  if (h_EWK) delete h_EWK;
  if (h_TT) delete h_TT;
  if (h_VV) delete h_VV;
  if (h_VVV) delete h_VVV;
  return tot;
}
static double integralSafe(TH1* h) { return h ? h->Integral(0, h->GetNbinsX()+1) : 0.0; }

// quick Gaussian fit in Z window
static bool fitZ(TH1* hin, double& mean, double& sigma) {
  if (!hin || hin->GetEntries() == 0) return false;
  auto* h = (TH1*)hin->Clone("fitZ_clone"); h->SetDirectory(nullptr);
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
  const TString md = modeDir(mode); // "Rochester" | "KIT" | "NoMuCorr"
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
void Muon_MuCorrComparison(TString era = "12022")
{
  gStyle->SetOptStat(0);
  const TString base = "/eos/user/d/dblend/Analyzed_SkimmedFiles/" + era + "/";
  const double  lumi = eraLumiPb(era);

// Per-mode Data (fallback to legacy flat Data if missing)
TFile* f_dataROCH = openDataForMode(base, "Rochester");
TFile* f_dataKIT  = openDataForMode(base, "KIT");
TFile* f_dataNOM  = openDataForMode(base, "NoMuCorr");
// If all three are null, bail out
if (!f_dataROCH && !f_dataKIT && !f_dataNOM) {
  std::cerr << "[ERR] No Data files found in any mode or legacy path.\n";
  return;
}

  const std::vector<TString> samples = {
    "DY","EWK", "TT","WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L",
    "WWW","WWZ","WZZ","ZZZ",
    "VBF"
  };

  // Only DM preselection & DM_VBF regions
  const std::string preDM = "corr_DM_Pre-Selection";
  //const std::string dmVBF = "corr_DM_VBF";
  const std::vector<std::string> regions = {"", "_BB", "_BO", "_BE", "_OO", "_OE", "_EE"};
  std::vector<std::string> folders;
  for (const auto& s : regions) folders.push_back(preDM + s + "/");
  //for (const auto& s : regions) folders.push_back(dmVBF + s + "/");

  const std::string var = "DiMuonMass_Z";
  const char* xlab = "m_{#mu#mu} [GeV]";
  const std::vector<TString> modes = {"KIT", "Rochester", "NoMuCorr"};

  // open per-mode file maps
  auto fmapKIT  = openModeFiles(base, "KIT",        samples);
  auto fmapROCH = openModeFiles(base, "Rochester",  samples);
  auto fmapNOM  = openModeFiles(base, "NoMuCorr",   samples);

  // show stamps (DY representative)
  std::cerr << "[STAMP] DY Rochester=" << readStamp(fmapROCH.count("DY")? fmapROCH.at("DY"):nullptr, "muCorrMode")
            << " KIT=" << readStamp(fmapKIT.count("DY")? fmapKIT.at("DY"):nullptr, "muCorrMode")
            << " NoMuCorr=" << readStamp(fmapNOM.count("DY")? fmapNOM.at("DY"):nullptr, "muCorrMode") << "\n";

  for (const auto& folder : folders) {
    // totals for combined ratio
    TH1* totKIT  = buildTotalMC(fmapKIT,  folder, var, "totKIT");
    TH1* totROCH = buildTotalMC(fmapROCH, folder, var, "totROCH");
    TH1* totNOM  = buildTotalMC(fmapNOM,  folder, var, "totNOM");

    // ---------- Per-mode stacked plots ----------
    for (const auto& MODE : modes) {
      const std::map<std::string,TFile*>& fmap =
        (MODE=="KIT") ? fmapKIT : (MODE=="Rochester" ? fmapROCH : fmapNOM);

// Mode-matched Data
  TFile* f_dataMODE = (MODE=="KIT") ? f_dataKIT : (MODE=="Rochester" ? f_dataROCH : f_dataNOM);
  TH1* hD = getHist(f_dataMODE, folder, var);
  TH1* h_data = nullptr;
  if (hD) { h_data = (TH1*)hD->Clone(Form("h_data_%s", MODE.Data())); h_data->SetDirectory(nullptr); }

      TH1* h_DY  = sumGroup(fmap, {"DY"}, folder, var, "h_DY");
      TH1* h_EWK  = sumGroup(fmap, {"EWK"}, folder, var, "h_EWK");
      TH1* h_TT  = sumGroup(fmap, {"TT"}, folder, var, "h_TT");
      TH1* h_VV  = sumGroup(fmap, {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"}, folder, var, "h_VV");
      TH1* h_VVV = sumGroup(fmap, {"WWW","WWZ","WZZ","ZZZ"}, folder, var, "h_VVV");
      TH1* h_SIG = sumGroup(fmap, {"VBF"}, folder, var, "h_SIG");

      if (h_SIG) { h_SIG->SetLineColor(kMagenta+2); h_SIG->SetLineWidth(2); h_SIG->SetFillStyle(0); }
      if (h_TT)  { h_TT ->SetFillColor(kRed-4);    h_TT ->SetLineColor(kRed-4);    }
      if (h_DY)  { h_DY ->SetFillColor(kOrange-3); h_DY ->SetLineColor(kOrange-3); }
      if (h_EWK)  { h_EWK ->SetFillColor(kTeal-6); h_EWK ->SetLineColor(kTeal+2); }
      if (h_VV)  { h_VV ->SetFillColor(kAzure-9);  h_VV ->SetLineColor(kAzure-9);  }
      if (h_VVV) { h_VVV->SetFillColor(kSpring+4); h_VVV->SetLineColor(kSpring+4); }

      if (h_data) { h_data->SetMarkerStyle(20); h_data->SetMarkerSize(0.8); h_data->SetLineColor(kBlack); }

      THStack st("st","");
      TH1* h_tot = nullptr;
      auto addToStack = [&](TH1* h){
        if (!h) return;
        st.Add(h);
        if (!h_tot) { h_tot = (TH1*)h->Clone("h_tot"); h_tot->SetDirectory(nullptr); }
        else        { h_tot->Add(h); }
      };
      addToStack(h_VVV); addToStack(h_VV); addToStack(h_TT); addToStack(h_DY); addToStack(h_EWK);
      const bool hasStack = (st.GetHists() != nullptr) && (st.GetNhists() > 0);

      TCanvas c("c","c",900,800);
      TPad pad1("pad1","pad1",0,0.30,1,1);
      TPad pad2("pad2","pad2",0,0.00,1,0.30);
      pad1.SetBottomMargin(0.02);
      pad2.SetTopMargin(0.03);
      pad2.SetBottomMargin(0.35);
      pad1.SetLogy(true);
      pad1.Draw(); pad2.Draw();

      // --- Top pad ---
      pad1.cd();
      double ymax = 1.0;
      if (h_tot)   ymax = std::max(ymax, 1.5*h_tot->GetMaximum());
      if (h_data)  ymax = std::max(ymax, 1.5*h_data->GetMaximum());
      if (ymax <= 0) ymax = 1.0;

      TH1* frame = nullptr;
      if (hasStack) {
        st.Draw("hist");
        st.SetMinimum(1e-1);
        if (st.GetHistogram()) {
          st.GetHistogram()->GetYaxis()->SetTitle("Events");
          st.GetHistogram()->GetXaxis()->SetTitle("");
          st.SetMaximum(ymax);
        }
        if (h_SIG)   h_SIG->Draw("hist same");
        if (h_data)  h_data->Draw("e1 same");
      } else {
        if      (h_data) frame = (TH1*)h_data->Clone("frame");
        else if (h_tot)  frame = (TH1*)h_tot->Clone("frame");
        if (frame) {
          frame->Reset(); frame->SetDirectory(nullptr);
          frame->GetYaxis()->SetTitle("Events");
          frame->GetXaxis()->SetTitle("");
          frame->GetYaxis()->SetRangeUser(1e-1, ymax);
          frame->Draw("axis");
        }
        if (h_SIG)   h_SIG->Draw("hist same");
        if (h_data)  h_data->Draw("e1 same");
      }

      // Legend
      TLegend leg(0.63,0.60,0.88,0.88);
      leg.SetBorderSize(0); leg.SetFillStyle(0);
      if (h_SIG)   leg.AddEntry(h_SIG,  "VBF H#rightarrow#mu#mu", "l");
      if (h_TT)    leg.AddEntry(h_TT,   "t#bar{t} + Jets", "f");
      if (h_DY)    leg.AddEntry(h_DY,   "DY + Jets", "f");
      if (h_EWK)    leg.AddEntry(h_EWK,   "EWK + Jets", "f");
      if (h_VV)    leg.AddEntry(h_VV,   "VV + Jets", "f");
      if (h_VVV)   leg.AddEntry(h_VVV,  "VVV + Jets", "f");
      if (h_data)  leg.AddEntry(h_data, "Data", "lep");
      leg.Draw();

      // CMS-like label
      TPaveText pt(0.12,0.92,0.95,0.99,"NDC");
      pt.SetFillStyle(0); pt.SetBorderSize(0);
      pt.AddText(Form("CMS Preliminary   #sqrt{s}=13.6 TeV,  L=%.1f fb^{-1}   Era %s",
                      lumi/1000.0, era.Data()));
      pt.Draw();

      // --- Bottom pad: Data/MC ratio ---
      pad2.cd();
      TH1* h_ratio = nullptr;
      if (h_data && h_tot) {
        h_ratio = (TH1*)h_data->Clone("h_ratio");
        h_ratio->SetDirectory(nullptr);
        h_ratio->Divide(h_tot);
        h_ratio->GetYaxis()->SetTitle("Data/MC");
        h_ratio->GetXaxis()->SetTitle(xlab);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.11);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.11);
        h_ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
        h_ratio->Draw("e1");
        // dashed unity
        TLine l1(h_ratio->GetXaxis()->GetXmin(),1.0,
                 h_ratio->GetXaxis()->GetXmax(),1.0);
        l1.SetLineStyle(2); l1.Draw("same");
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
        else if (h_tot)  ax = (TH1*)h_tot->Clone("ax");
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
      TString outdirMode = "/eos/user/d/dblend/StackedPlots/" + era + "/MuCorr_Comparison/" + modeDir(MODE) + "/";
      gSystem->Exec("mkdir -p " + outdirMode);
      TString out = outdirMode + "StackWithRatio_" + clean + ".pdf";
      c.SaveAs(out);
      c.SaveAs(out.ReplaceAll(".pdf",".png"));

      // cleanup
      delete frame; delete h_ratio;
      if (h_DY) delete h_DY;
      if (h_EWK) delete h_EWK;
      if (h_TT) delete h_TT;
      if (h_VV) delete h_VV;
      if (h_VVV) delete h_VVV;
      if (h_SIG) delete h_SIG;
    } // per-mode

    // ---------- Combined Data/MC ratio overlay ----------
    TCanvas cR("cR","cR",900,700);
    TPad padR("padR","padR",0,0.00,1,1);
    padR.SetBottomMargin(0.13);
    padR.Draw();
    padR.cd();

std::vector<std::pair<TString,TH1*>> ratios;

// Build a Data/MC ratio for EACH mode using that mode's Data
TH1 *dR=nullptr, *dK=nullptr, *dN=nullptr;
if (f_dataROCH) { if (auto* h = getHist(f_dataROCH, folder, var)) { dR=(TH1*)h->Clone("d_roch"); dR->SetDirectory(nullptr); } }
if (f_dataKIT)  { if (auto* h = getHist(f_dataKIT,  folder, var)) { dK=(TH1*)h->Clone("d_kit");  dK->SetDirectory(nullptr); } }
if (f_dataNOM)  { if (auto* h = getHist(f_dataNOM,  folder, var)) { dN=(TH1*)h->Clone("d_nom");  dN->SetDirectory(nullptr); } }

if (dR && totROCH) { dR->Divide(totROCH); ratios.push_back({"Rochester", dR}); }
if (dK && totKIT)  { dK->Divide(totKIT);  ratios.push_back({"KIT",       dK}); }
if (dN && totNOM)  { dN->Divide(totNOM);  ratios.push_back({"NoMuCorr",  dN}); }

    double ymin = 0.0, ymax = 2.0;
    bool firstDrawn = false;
    TLegend legR(0.50,0.75,0.88,0.90);
    legR.SetBorderSize(0); legR.SetFillStyle(0);

    for (auto& pr : ratios) {
      TH1* h = pr.second; if (!h) continue;
      int col = kBlack;
      if (pr.first=="Rochester") col = kBlue+1;
      if (pr.first=="KIT")       col = kGreen+2;
      if (pr.first=="NoMuCorr")  col = kMagenta+2;
      h->SetLineColor(col);
      h->SetMarkerColor(col);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(0.7);
      h->GetYaxis()->SetTitle("Data/MC");
      h->GetXaxis()->SetTitle(xlab);
      h->GetYaxis()->SetRangeUser(ymin, ymax);
      if (!firstDrawn) { h->Draw("e1"); firstDrawn = true; } else { h->Draw("e1 same"); }
      // connecting line
      TGraphErrors gr(h); gr.SetLineColor(col); gr.SetLineWidth(2); gr.SetMarkerSize(0); gr.Draw("L same");
      legR.AddEntry(h, pr.first, "lep");
    }

    if (firstDrawn) {
      TH1* href = ratios.front().second;
      TLine l1(href->GetXaxis()->GetXmin(),1.0, href->GetXaxis()->GetXmax(),1.0);
      l1.SetLineStyle(2); l1.Draw("same");
    }

// Build a Data histogram for the combined Z-fit (use Rochester/Data)
TH1* hDataR = nullptr;
{
  // try per-mode Data first
  std::unique_ptr<TFile> fR(TFile::Open(base + "Rochester/DATA/DATA_all.root"));
  if (!fR) {
    // fallback to legacy flat Data if per-mode is missing
    std::cerr << "[WARN] Missing " << (base + "Rochester/DATA/DATA_all.root")
              << " ; falling back to " << (base + "Data/Data_all.root") << "\n";
    fR.reset(TFile::Open(base + "Data/Data_all.root"));
  }
  if (fR) {
    if (auto* h = getHist(fR.get(), folder, var)) {
      hDataR = (TH1*)h->Clone("hDataR");
      hDataR->SetDirectory(nullptr);
    }
  }
}

    // ---- Z-fit means (vertical lines) ----
double mData=0, sData=0, mR=0, sR=0, mK=0, sK=0, mN=0, sN=0;
bool okD = fitZ(hDataR,  mData, sData);
bool okR = fitZ(totROCH, mR,    sR);
bool okK = fitZ(totKIT,  mK,    sK);
bool okN = fitZ(totNOM,  mN,    sN);

    if (firstDrawn && (okD||okR||okK||okN)) {
      TH1* href = ratios.front().second;
      double ylow = 0.0, yhigh = 2.0;
      auto drawV = [&](double x, int col, int style){
        if (x<=0 || !href) return;
        TLine ln(x, ylow, x, yhigh);
        ln.SetLineColor(col);
        ln.SetLineStyle(style);
        ln.SetLineWidth(2);
        ln.Draw("same");
      };
      if (okD) drawV(mData, kBlack,     1);
      if (okR) drawV(mR,    kBlue+1,    7);
      if (okK) drawV(mK,    kGreen+2,   9);
      if (okN) drawV(mN,    kMagenta+2, 2);
    }

    legR.Draw();
    TPaveText ptR(0.12,0.92,0.95,0.99,"NDC");
    ptR.SetFillStyle(0); ptR.SetBorderSize(0);
    ptR.AddText(Form("CMS Preliminary   #sqrt{s}=13.6 TeV,  L=%.1f fb^{-1}   Era %s",
                    lumi/1000.0, era.Data()));
    ptR.Draw();

    TString clean = TString(folder).ReplaceAll("/","").ReplaceAll("corr_","");
    TString outdirComp = "/eos/user/d/dblend/StackedPlots/" + era + "/MuCorr_Comparison/Combined/";
    gSystem->Exec("mkdir -p " + outdirComp);
    TString outR = outdirComp + "RatioOverlay_" + clean + ".pdf";
    cR.SaveAs(outR);
    cR.SaveAs(outR.ReplaceAll(".pdf",".png"));

    // ---------- TXT summary with Z-fit + yields ----------
    TString txtPath = outdirComp + "MuCorr_Summary_" + clean + ".txt";
    std::ofstream ofs(txtPath.Data());
    ofs << std::fixed << std::setprecision(3);
    ofs << "# Summary for folder " << folder << " (DiMuonMass_Z)\n";
    double yData = integralSafe(hDataR);
    ofs << "Data: " << yData << "\n";

    auto writeMode = [&](const char* label, const std::map<std::string,TFile*>& fmap, TH1* totH){
      TH1* c_DY  = sumGroup(fmap, {"DY"}, folder, var, "c_DY");
      TH1* c_EWK  = sumGroup(fmap, {"EWK"}, folder, var, "c_EWK");
      TH1* c_TT  = sumGroup(fmap, {"TT"}, folder, var, "c_TT");
      TH1* c_VV  = sumGroup(fmap, {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"}, folder, var, "c_VV");
      TH1* c_VVV = sumGroup(fmap, {"WWW","WWZ","WZZ","ZZZ"}, folder, var, "c_VVV");
      double yDY  = integralSafe(c_DY);
      double yEWK  = integralSafe(c_EWK);
      double yTT  = integralSafe(c_TT);
      double yVV  = integralSafe(c_VV);
      double yVVV = integralSafe(c_VVV);
      double yTot = integralSafe(totH);
      double ratio = (yTot>0 ? yData/yTot : 0.0);
      ofs << "\n[" << label << "]\n";
      ofs << " MC_total: " << yTot << "\n";
      ofs << "   DY:  " << yDY << "\n";
      ofs << "   EWK:  " << yEWK << "\n";
      ofs << "   TT:  " << yTT << "\n";
      ofs << "   VV:  " << yVV << "\n";
      ofs << "   VVV: " << yVVV << "\n";
      ofs << " Data/MC: " << ratio << "\n";
      if (c_DY) delete c_DY;
      if (c_EWK) delete c_EWK;
      if (c_TT) delete c_TT;
      if (c_VV) delete c_VV;
      if (c_VVV) delete c_VVV;
    };
    writeMode("Rochester", fmapROCH, totROCH);
    writeMode("KIT",       fmapKIT,  totKIT);
    writeMode("NoMuCorr",  fmapNOM,  totNOM);

    ofs << "\n# Z-peak Gaussian fit (80–100 GeV)\n";
    if (okD) ofs << "Data:     mean=" << mData << "  sigma=" << sData << "\n";
    if (okR) ofs << "Rochester mean=" << mR    << "  sigma=" << sR    << "\n";
    if (okK) ofs << "KIT:      mean=" << mK    << "  sigma=" << sK    << "\n";
    if (okN) ofs << "NoMuCorr: mean=" << mN    << "  sigma=" << sN    << "\n";
    ofs.close();

    // cleanup
    if (hDataR) { delete hDataR; hDataR = nullptr; }
    if (totKIT) delete totKIT;
    if (totROCH) delete totROCH;
    if (totNOM) delete totNOM;
    for (auto& pr : ratios) { /* clones cleaned above */ }
  }

  closeModeFiles(fmapKIT);
  closeModeFiles(fmapROCH);
  closeModeFiles(fmapNOM);
if (f_dataROCH) f_dataROCH->Close();
if (f_dataKIT)  f_dataKIT->Close();
if (f_dataNOM)  f_dataNOM->Close();

  std::cout << "[OK] Wrote per-mode stacked+ratio plots, combined ratio overlays, and TXT summaries with stamps and Z-fit.\n";
}
