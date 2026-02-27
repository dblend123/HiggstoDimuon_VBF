// root -l -b -q 'Jet_Overlay_Comparison.C("22023","/path/setA/","/path/setB/")'

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

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
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TGraphErrors.h>

// ---------- helpers ----------
static double eraLumiPb(const TString& era) {
  if (era == "12022") return 7980.454145;
  if (era == "22022") return 26671.60971;
  if (era == "12023") return 18704.13341;
  if (era == "22023") return 9693.130053;
  if (era == "2024")  return 0.0; // fill when ready
  std::cerr << "[WARN] Unknown era " << era << ", lumi set to 0\n";
  return 0.0;
}

static TH1* getHist(TFile* f, const std::string& folder, const std::string& name) {
  if (!f) return nullptr;
  TString path = TString(folder.c_str()) + name.c_str();
  TH1* h = dynamic_cast<TH1*>(f->Get(path));
  if (!h) {
    TString alt = TString(folder.c_str());
    if (!alt.EndsWith("/")) alt += "/";
    alt += name.c_str();
    h = dynamic_cast<TH1*>(f->Get(alt));
  }
  return h;
}

static std::map<std::string,TFile*> openFilesFromBase(const TString& base,
                                                      const std::vector<TString>& samples) {
  std::map<std::string,TFile*> m;
  for (const auto& s : samples) {
    TString path = base;
    if (!path.EndsWith("/")) path += "/";
    path += s + "/" + s + "_all.root";
    m[s.Data()] = TFile::Open(path);
    if (!m[s.Data()]) std::cerr << "[WARN] Missing file: " << path << "\n";
  }
  return m;
}

static void closeFiles(std::map<std::string,TFile*>& m) {
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

static TFile* openDataFromBase(const TString& base) {
  TString b = base; if (!b.EndsWith("/")) b += "/";
  // preferred
  TFile* f = TFile::Open(b + "DATA/DATA_all.root");
  if (f) return f;
  // fallback
  f = TFile::Open(b + "Data/Data_all.root");
  return f;
}

static std::vector<std::string> listTH1NamesInFolder(TFile* fAny, const std::string& folderWithSlash) {
  std::vector<std::string> out;
  if (!fAny) return out;

  TDirectory* dir = dynamic_cast<TDirectory*>(fAny->Get(folderWithSlash.c_str()));
  if (!dir) {
    std::string f2 = folderWithSlash;
    while (!f2.empty() && f2.back()=='/') f2.pop_back();
    dir = dynamic_cast<TDirectory*>(fAny->Get(f2.c_str()));
  }
  if (!dir) return out;

  TIter nextkey(dir->GetListOfKeys());
  while (TKey* key = (TKey*)nextkey()) {
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;
    if (cl->InheritsFrom("TH1")) out.push_back(key->GetName());
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// ---------- MAIN ----------
void Jet_Overlay_Comparison(TString era,
                            TString baseSetA,
                            TString baseSetB)
{
  gStyle->SetOptStat(0);

  if (!baseSetA.EndsWith("/")) baseSetA += "/";
  if (!baseSetB.EndsWith("/")) baseSetB += "/";

  const double lumi = eraLumiPb(era);

  const std::vector<TString> samples = {
    "DY","EWK", "TT","WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L",
    "WWW","WWZ","WZZ","ZZZ",
    "VBF"
  };

  auto fmapA = openFilesFromBase(baseSetA, samples);
  auto fmapB = openFilesFromBase(baseSetB, samples);

  TFile* f_dataA = openDataFromBase(baseSetA);
  TFile* f_dataB = openDataFromBase(baseSetB);

  if (!f_dataA && !f_dataB) {
    std::cerr << "[ERR] Could not open Data in either set:\n"
              << "  A: " << baseSetA << "\n"
              << "  B: " << baseSetB << "\n";
    closeFiles(fmapA); closeFiles(fmapB);
    return;
  }

  const std::string jetFolder = "corr_Jet_Pre-Selection/";

  // choose any file to list histograms
  TFile* fList = (fmapA.count("DY") && fmapA["DY"]) ? fmapA["DY"]
               : (fmapB.count("DY") && fmapB["DY"]) ? fmapB["DY"]
               : (f_dataA ? f_dataA : f_dataB);

  auto vars = listTH1NamesInFolder(fList, jetFolder);
  if (vars.empty()) {
    std::cerr << "[ERR] No TH1s found under: " << jetFolder << "\n";
    closeFiles(fmapA); closeFiles(fmapB);
    if (f_dataA) f_dataA->Close();
    if (f_dataB) f_dataB->Close();
    return;
  }

  // output dirs
  TString outA = "/eos/user/d/dblend/StackedPlots/" + era + "/Jet_Comparison/SetA/";
  TString outB = "/eos/user/d/dblend/StackedPlots/" + era + "/Jet_Comparison/SetB/";
  TString outC = "/eos/user/d/dblend/StackedPlots/" + era + "/Jet_Comparison/Combined/";
  gSystem->Exec("mkdir -p " + outA);
  gSystem->Exec("mkdir -p " + outB);
  gSystem->Exec("mkdir -p " + outC);

  // precompute total MC (per histogram, inside loop)
  for (const auto& var : vars) {

    TH1* totA = buildTotalMC(fmapA, jetFolder, var, "totA");
    TH1* totB = buildTotalMC(fmapB, jetFolder, var, "totB");

    auto drawOne = [&](const TString& tag,
                       const std::map<std::string,TFile*>& fmap,
                       TFile* f_data,
                       TH1* totMC,
                       const TString& outdir)
    {
      TH1* hD = getHist(f_data, jetFolder, var);
      TH1* h_data = nullptr;
      if (hD) { h_data = (TH1*)hD->Clone("h_data"); h_data->SetDirectory(nullptr); }

      TH1* h_DY   = sumGroup(fmap, {"DY"}, jetFolder, var, "h_DY");
      TH1* h_EWK  = sumGroup(fmap, {"EWK"}, jetFolder, var, "h_EWK");
      TH1* h_TT   = sumGroup(fmap, {"TT"}, jetFolder, var, "h_TT");
      TH1* h_VV   = sumGroup(fmap, {"WZto2L2Q","WZto3LNu","ZZto2L2Nu","ZZto2L2Q","ZZto4L"}, jetFolder, var, "h_VV");
      TH1* h_VVV  = sumGroup(fmap, {"WWW","WWZ","WZZ","ZZZ"}, jetFolder, var, "h_VVV");
      TH1* h_SIG  = sumGroup(fmap, {"VBF"}, jetFolder, var, "h_SIG");

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

      TCanvas c("c","c",900,800);
      TPad pad1("pad1","pad1",0,0.30,1,1);
      TPad pad2("pad2","pad2",0,0.00,1,0.30);
      pad1.SetBottomMargin(0.02);
      pad2.SetTopMargin(0.03);
      pad2.SetBottomMargin(0.35);
      pad1.Draw(); pad2.Draw();

      pad1.cd();
      double ymax = 1.0;
      auto updMax = [&](TH1* h){ if (h) ymax = std::max(ymax, 1.5*h->GetMaximum()); };
      updMax(totMC); updMax(h_data); updMax(h_DY); updMax(h_EWK); updMax(h_TT); updMax(h_VV); updMax(h_VVV); updMax(h_SIG);
      bool useLogY = ((totMC && totMC->Integral() > 0) || (h_data && h_data->Integral() > 0));
      pad1.SetLogy(useLogY);
      double ymin = useLogY ? 1e-1 : 0.0;

      // frame
      TH1* ref = totMC ? totMC : (h_data ? h_data : (h_DY ? h_DY : (h_EWK ? h_EWK : (h_TT ? h_TT : (h_VV ? h_VV : (h_VVV ? h_VVV : h_SIG))))));
      TH1* frame = nullptr;
      if (ref) {
        frame = (TH1*)ref->Clone("frame");
        frame->Reset();
        frame->SetDirectory(nullptr);
        frame->GetYaxis()->SetTitle("Events");
        frame->GetXaxis()->SetTitle("");
        frame->GetYaxis()->SetRangeUser(ymin, ymax);
        frame->Draw("axis");
      }

      if (h_VVV) h_VVV->Draw("hist same");
      if (h_VV)  h_VV ->Draw("hist same");
      if (h_TT)  h_TT ->Draw("hist same");
      if (h_DY)  h_DY ->Draw("hist same");
      if (h_EWK) h_EWK->Draw("hist same");
      if (h_SIG) h_SIG->Draw("hist same");
      if (h_data) h_data->Draw("e1 same");

      TLegend leg(0.63,0.60,0.88,0.88);
      leg.SetBorderSize(0);
      leg.SetFillStyle(0);
      if (h_SIG)  leg.AddEntry(h_SIG,  "VBF H#rightarrow#mu#mu", "l");
      if (h_TT)   leg.AddEntry(h_TT,   "t#bar{t} + Jets", "f");
      if (h_DY)   leg.AddEntry(h_DY,   "DY + Jets", "f");
      if (h_EWK)  leg.AddEntry(h_EWK,  "EWK + Jets", "f");
      if (h_VV)   leg.AddEntry(h_VV,   "VV + Jets", "f");
      if (h_VVV)  leg.AddEntry(h_VVV,  "VVV + Jets", "f");
      if (h_data) leg.AddEntry(h_data, "Data", "lep");
      leg.Draw();

      TPaveText pt(0.12,0.92,0.95,0.99,"NDC");
      pt.SetFillStyle(0);
      pt.SetBorderSize(0);
      pt.AddText(Form("CMS Preliminary   #sqrt{s}=13.6 TeV,  L=%.1f fb^{-1}   Era %s   (%s)",
                      lumi/1000.0, era.Data(), tag.Data()));
      pt.Draw();

      pad2.cd();
      TH1* h_ratio = nullptr;
      if (h_data && totMC) {
        h_ratio = (TH1*)h_data->Clone("h_ratio");
        h_ratio->SetDirectory(nullptr);
        h_ratio->Divide(totMC);
        h_ratio->GetYaxis()->SetTitle("Data/MC");
        h_ratio->GetXaxis()->SetTitle(var.c_str());
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.11);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetXaxis()->SetTitleSize(0.12);
        h_ratio->GetXaxis()->SetLabelSize(0.11);
        h_ratio->GetYaxis()->SetRangeUser(0.7, 1.3);
        h_ratio->Draw("e1");
        TLine l1(h_ratio->GetXaxis()->GetXmin(),1.0,
                 h_ratio->GetXaxis()->GetXmax(),1.0);
        l1.SetLineStyle(2);
        l1.Draw("same");

        TGraphErrors gr(h_ratio);
        gr.SetLineColor(kGray+2);
        gr.SetLineWidth(2);
        gr.SetMarkerStyle(20);
        gr.SetMarkerSize(0);
        gr.Draw("L same");
      }

      TString out = outdir + "StackWithRatio_" + TString(var) + ".pdf";
      c.SaveAs(out);
      c.SaveAs(out.ReplaceAll(".pdf",".png"));

      if (frame) delete frame;
      if (h_ratio) delete h_ratio;
      if (h_DY) delete h_DY;
      if (h_EWK) delete h_EWK;
      if (h_TT) delete h_TT;
      if (h_VV) delete h_VV;
      if (h_VVV) delete h_VVV;
      if (h_SIG) delete h_SIG;
      if (h_data) delete h_data;
    };

    // per-set plots
    drawOne("SetA", fmapA, f_dataA, totA, outA);
    drawOne("SetB", fmapB, f_dataB, totB, outB);

    // combined overlay: (Data/MC)_A and (Data/MC)_B
    TH1 *dA=nullptr, *dB=nullptr;
    if (f_dataA) if (auto* h = getHist(f_dataA, jetFolder, var)) { dA=(TH1*)h->Clone("dA"); dA->SetDirectory(nullptr); }
    if (f_dataB) if (auto* h = getHist(f_dataB, jetFolder, var)) { dB=(TH1*)h->Clone("dB"); dB->SetDirectory(nullptr); }

    TCanvas cR("cR","cR",900,700);
    bool first = false;
    TLegend legR(0.55,0.78,0.88,0.90);
    legR.SetBorderSize(0);
    legR.SetFillStyle(0);

    if (dA && totA) {
      dA->Divide(totA);
      dA->SetLineColor(kGreen+2);
      dA->SetMarkerColor(kGreen+2);
      dA->SetMarkerStyle(20);
      dA->SetMarkerSize(0.7);
      dA->GetYaxis()->SetTitle("Data/MC");
      dA->GetXaxis()->SetTitle(var.c_str());
      dA->GetYaxis()->SetRangeUser(0.7, 1.3);
      dA->Draw("e1");
      first = true;
      legR.AddEntry(dA, "SetA", "lep");
      TGraphErrors grA(dA); grA.SetLineColor(kGreen+2); grA.SetLineWidth(2); grA.SetMarkerSize(0); grA.Draw("L same");
    }
    if (dB && totB) {
      dB->Divide(totB);
      dB->SetLineColor(kBlue+1);
      dB->SetMarkerColor(kBlue+1);
      dB->SetMarkerStyle(20);
      dB->SetMarkerSize(0.7);
      dB->GetYaxis()->SetTitle("Data/MC");
      dB->GetXaxis()->SetTitle(var.c_str());
      dB->GetYaxis()->SetRangeUser(0.7, 1.3);
      if (!first) { dB->Draw("e1"); first = true; }
      else dB->Draw("e1 same");
      legR.AddEntry(dB, "SetB", "lep");
      TGraphErrors grB(dB); grB.SetLineColor(kBlue+1); grB.SetLineWidth(2); grB.SetMarkerSize(0); grB.Draw("L same");
    }

    if (first) {
      TH1* href = dA ? dA : dB;
      TLine l1(href->GetXaxis()->GetXmin(),1.0, href->GetXaxis()->GetXmax(),1.0);
      l1.SetLineStyle(2);
      l1.Draw("same");
      legR.Draw();

      TPaveText ptR(0.12,0.92,0.95,0.99,"NDC");
      ptR.SetFillStyle(0);
      ptR.SetBorderSize(0);
      ptR.AddText(Form("CMS Preliminary   #sqrt{s}=13.6 TeV,  L=%.1f fb^{-1}   Era %s",
                       lumi/1000.0, era.Data()));
      ptR.Draw();

      TString out = outC + "RatioOverlay_" + TString(var) + ".pdf";
      cR.SaveAs(out);
      cR.SaveAs(out.ReplaceAll(".pdf",".png"));
    }

    if (dA) delete dA;
    if (dB) delete dB;
    if (totA) delete totA;
    if (totB) delete totB;
  }

  closeFiles(fmapA);
  closeFiles(fmapB);
  if (f_dataA) f_dataA->Close();
  if (f_dataB) f_dataB->Close();

  std::cout << "[OK] Done. Wrote to /eos/user/d/dblend/StackedPlots/" << era << "/Jet_Comparison/\n";
}

