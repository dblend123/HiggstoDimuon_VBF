#include <TFile.h>
#include <TH1F.h>
#include <iostream>

void DrawCMSLabel() {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.25, 0.82, "#font[62]{CMS}");
    latex.DrawLatex(0.25, 0.78, "#it{Simulation}");
    latex.DrawLatex(0.25, 0.74, "#it{Preliminary}");
    latex.DrawLatex(0.22, 0.65, "#font[62]{Signal MC 2022}");
}

void nice_plots(){
TFile *file = new TFile("output_Muon_FullSignalRun22.root");
if (!file || file->IsZombie()) {
        std::cerr << "Error: Failed to open file " << "output_Muon_FullSignalRun22.root" << std::endl;
        return;
    }

TDirectory *DMdir = (TDirectory*)file->Get("DM_VBF");
if (!DMdir) {
    std::cerr << "Error: Directory 'DM_VBF' not found in file." << std::endl;
    return;
}

TDirectory *Jetdir = (TDirectory*)file->Get("Jet_VBF");
if (!Jetdir) {
    std::cerr << "Error: Directory 'Jet_VBF' not found in file." << std::endl;
    return;
}

TH1F *DMhist_higgs = (TH1F*)DMdir->Get("DiMuonMass_SB_SR2");
if (!DMhist_higgs) {
    std::cerr << "Error: Histogram 'DiMuonMass_SB_SR2' not found in directory 'DM_VBF'." << std::endl;
    return;
}

TH1F *DMhist_Z = (TH1F*)DMdir->Get("DiMuonMass_Z");
if (!DMhist_Z) {
    std::cerr << "Error: Histogram 'DiMuonMass_Z' not found in directory 'DM_VBF'." << std::endl;
    return;
}

TH1F *Jethistpt = (TH1F*)Jetdir->Get("DiJetPt");
if (!Jethistpt) {
    std::cerr << "Error: Histogram 'DiJetPt' not found in directory 'Jet_VBF'." << std::endl;
    return;
}

TH1F *Jethisteta = (TH1F*)Jetdir->Get("DiJetdeta");
if (!Jethisteta) {
    std::cerr << "Error: Histogram 'DiJetdeta' not found in directory 'Jet_VBF'." << std::endl;
    return;
}

TH1F *Jethistmass = (TH1F*)Jetdir->Get("DiJetMass");
if (!Jethistmass) {
    std::cerr << "Error: Histogram 'DiJetMass' not found in directory 'Jet_VBF'." << std::endl;
    return;
}

DMhist_higgs->GetXaxis()->SetTitle("Invariant Mass [GeV]");
DMhist_higgs->GetYaxis()->SetTitle("Events");
DMhist_higgs->SetTitle("VBF SideBand/Signal Dimuon Mass Region");
DMhist_higgs->SetFillColor(kGreen-9);

DMhist_Z->GetXaxis()->SetTitle("Invariant Mass [GeV]");
DMhist_Z->GetYaxis()->SetTitle("Events");
DMhist_Z->SetTitle("VBF Dimuon Z-Mass Region");
DMhist_Z->SetFillColor(kGreen-9);

Jethistpt->GetXaxis()->SetTitle("DiJet p_{T} [GeV]");
Jethistpt->GetYaxis()->SetTitle("Events");
Jethistpt->SetTitle("VBF DiJet p_{T}");
Jethistpt->SetFillColor(kGreen-9);

Jethisteta->GetXaxis()->SetTitle("#Delta#eta_{jj}");
Jethisteta->GetYaxis()->SetTitle("Events");
Jethisteta->SetTitle("VBF DiJet #Delta#eta_{jj}");
Jethisteta->SetFillColor(kGreen-9);

Jethistmass->GetXaxis()->SetTitle("DiJet Mass [GeV]");
Jethistmass->GetYaxis()->SetTitle("Events");
Jethistmass->SetTitle("VBF DiJet Mass");
Jethistmass->SetFillColor(kGreen-9);

TCanvas *DMcanvas_higgs = new TCanvas("DMcanvas_higgs", "DM_higgs Histogram Canvas", 800, 600);
DMcanvas_higgs->SetTickx();
DMcanvas_higgs->SetTicky();
DMcanvas_higgs->SetGridx();
DMcanvas_higgs->SetGridy();
gPad->SetLeftMargin(0.15);

TCanvas *DMcanvas_Z = new TCanvas("DMcanvas_Z", "DM_Z Histogram Canvas", 800, 600);
DMcanvas_Z->SetTickx();
DMcanvas_Z->SetTicky();
DMcanvas_Z->SetGridx();
DMcanvas_Z->SetGridy();
gPad->SetLeftMargin(0.15);

TCanvas *Jetcanvaspt = new TCanvas("Jetcanvaspt", "Jetpt Histogram Canvas", 800, 600);
Jetcanvaspt->SetTickx();
Jetcanvaspt->SetTicky();
Jetcanvaspt->SetGridx();
Jetcanvaspt->SetGridy();
gPad->SetLeftMargin(0.15);

TCanvas *Jetcanvaseta = new TCanvas("Jetcanvaseta", "Jeteta Histogram Canvas", 800, 600);
Jetcanvaseta->SetTickx();
Jetcanvaseta->SetTicky();
Jetcanvaseta->SetGridx();
Jetcanvaseta->SetGridy();
gPad->SetLeftMargin(0.15);

TCanvas *Jetcanvasmass = new TCanvas("Jetcanvasmass", "Jetmass Histogram Canvas", 800, 600);
Jetcanvasmass->SetTickx();
Jetcanvasmass->SetTicky();
Jetcanvasmass->SetGridx();
Jetcanvasmass->SetGridy();
gPad->SetLeftMargin(0.15);

Int_t DMnumBins_higgs = DMhist_higgs->GetNbinsX();
Double_t DMtotalEvents_higgs = 0.0;
for (Int_t i = 1; i <= DMnumBins_higgs; ++i) {
    Double_t binContent = DMhist_higgs->GetBinContent(i);
    DMtotalEvents_higgs += binContent;
    std::cout << "Bin " << i << ": " << binContent << " events" << std::endl;
}
std::cout << "Total number of events: " << DMtotalEvents_higgs << std::endl;
DMhist_higgs->SetEntries(DMtotalEvents_higgs);
DMcanvas_higgs->cd();
DMhist_higgs->Draw("HIST");
DrawCMSLabel();

Int_t DMnumBins_Z = DMhist_Z->GetNbinsX();
Double_t DMtotalEvents_Z = 0.0;
for (Int_t i = 1; i <= DMnumBins_Z; ++i) {
    Double_t binContent = DMhist_Z->GetBinContent(i);
    DMtotalEvents_Z += binContent;
    std::cout << "Bin " << i << ": " << binContent << " events" << std::endl;
}
std::cout << "Total number of events: " << DMtotalEvents_Z << std::endl;
DMhist_Z->SetEntries(DMtotalEvents_Z);
DMcanvas_Z->cd();
DMhist_Z->Draw("HIST");
DrawCMSLabel();

Int_t JetptnumBins = Jethistpt->GetNbinsX();
Double_t JetpttotalEvents = 0.0;
for (Int_t i = 1; i <= JetptnumBins; ++i) {
    Double_t binContent = Jethistpt->GetBinContent(i);
    JetpttotalEvents += binContent;
    std::cout << "Bin " << i << ": " << binContent << " events" << std::endl;
}
std::cout << "Total number of events: " << JetpttotalEvents << std::endl;
Jethistpt->SetEntries(JetpttotalEvents);
Jetcanvaspt->cd();
Jethistpt->Draw("HIST");
TLatex latex;
latex.SetNDC();
latex.SetTextFont(42);
latex.SetTextSize(0.04);
latex.DrawLatex(0.65, 0.65, "#font[62]{CMS}");
latex.DrawLatex(0.65, 0.61, "#it{Simulation}");
latex.DrawLatex(0.65, 0.56, "#it{Preliminary}");
latex.DrawLatex(0.62, 0.47, "#font[62]{Signal MC 2022}");

Int_t JetetanumBins = Jethisteta->GetNbinsX();
Double_t JetetatotalEvents = 0.0;
for (Int_t i = 1; i <= JetetanumBins; ++i) {
    Double_t binContent = Jethisteta->GetBinContent(i);
    JetetatotalEvents += binContent;
    std::cout << "Bin " << i << ": " << binContent << " events" << std::endl;
}
std::cout << "Total number of events: " << JetetatotalEvents << std::endl;
Jethisteta->SetEntries(JetetatotalEvents);
Jetcanvaseta->cd();
Jethisteta->Draw("HIST");
DrawCMSLabel();

Int_t JetmassnumBins = Jethistmass->GetNbinsX();
Double_t JetmasstotalEvents = 0.0;
for (Int_t i = 1; i <= JetmassnumBins; ++i) {
    Double_t binContent = Jethistmass->GetBinContent(i);
    JetmasstotalEvents += binContent;
    std::cout << "Bin " << i << ": " << binContent << " events" << std::endl;
}
std::cout << "Total number of events: " << JetmasstotalEvents << std::endl;
Jethistmass->SetEntries(JetmasstotalEvents);
Jetcanvasmass->cd();
Jethistmass->Draw("HIST");
DrawCMSLabel();

DMcanvas_higgs->SaveAs("Signal_DiMuonMass_higgs.png"); 
DMcanvas_Z->SaveAs("Signal_DiMuonMass_Z.png"); 
Jetcanvaspt->SaveAs("Signal_Jetpt.png"); 
Jetcanvaseta->SaveAs("Signal_Jeteta.png"); 
Jetcanvasmass->SaveAs("Signal_Jetmass.png"); 
}
