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
    latex.DrawLatex(0.22, 0.65, "#font[62]{Drell-Yan MC 2022}");
}

void nice_plots(){
TFile *file = new TFile("output_Muon_FullDYRun22.root");
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

TH1F *DMhist = (TH1F*)DMdir->Get("DiMuonMass_SB_SR2");
if (!DMhist) {
    std::cerr << "Error: Histogram 'DiMuonMass_SB_SR2' not found in directory 'DM_VBF'." << std::endl;
    return;
}

TH1F *Jethistpt = (TH1F*)Jetdir->Get("DiJetPt");
if (!Jethistpt) {
    std::cerr << "Error: Histogram 'DiJetPt' not found in directory 'DM_VBF'." << std::endl;
    return;
}

TH1F *Jethisteta = (TH1F*)Jetdir->Get("DiJetdeta");
if (!Jethisteta) {
    std::cerr << "Error: Histogram 'DiJetdeta' not found in directory 'DM_VBF'." << std::endl;
    return;
}

TH1F *Jethistmass = (TH1F*)Jetdir->Get("DiJetMass");
if (!Jethistmass) {
    std::cerr << "Error: Histogram 'DiJetMass' not found in directory 'DM_VBF'." << std::endl;
    return;
}

DMhist->GetXaxis()->SetTitle("Invariant Mass [GeV]");
DMhist->GetYaxis()->SetTitle("Events");
DMhist->SetTitle("VBF SideBand/Signal Dimuon Mass Region");
DMhist->SetFillColor(kGreen-9);

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

TCanvas *DMcanvas = new TCanvas("DMcanvas", "DM Histogram Canvas", 800, 600);
DMcanvas->SetTickx();
DMcanvas->SetTicky();
DMcanvas->SetGridx();
DMcanvas->SetGridy();
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

Int_t DMnumBins = DMhist->GetNbinsX();
Double_t DMtotalEvents = 0.0;
for (Int_t i = 1; i <= DMnumBins; ++i) {
    Double_t binContent = DMhist->GetBinContent(i);
    DMtotalEvents += binContent;
    std::cout << "Bin " << i << ": " << binContent << " events" << std::endl;
}
std::cout << "Total number of events: " << DMtotalEvents << std::endl;
DMhist->SetEntries(DMtotalEvents);
DMcanvas->cd();
DMhist->Draw();
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
Jethistpt->Draw();
TLatex latex;
latex.SetNDC();
latex.SetTextFont(42);
latex.SetTextSize(0.04);
latex.DrawLatex(0.65, 0.65, "#font[62]{CMS}");
latex.DrawLatex(0.65, 0.61, "#it{Simulation}");
latex.DrawLatex(0.65, 0.56, "#it{Preliminary}");
latex.DrawLatex(0.62, 0.47, "#font[62]{Drell-Yan MC 2022}");

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
Jethisteta->Draw();
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
Jethistmass->Draw();
DrawCMSLabel();

DMcanvas->SaveAs("DY_DiMuonMass_SR_SB2.png"); 
Jetcanvaspt->SaveAs("DY_Jetpt.png"); 
Jetcanvaseta->SaveAs("DY_Jeteta.png"); 
Jetcanvasmass->SaveAs("DY_Jetmass.png"); 
}
