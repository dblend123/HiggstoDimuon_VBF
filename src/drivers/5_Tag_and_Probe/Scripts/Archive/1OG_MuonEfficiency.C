#include <fstream>  
#include <iostream>
#include <vector>
#include <string>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace ROOT;
using namespace ROOT::VecOps;

void Efficiencies() {
    //ROOT::RDataFrame df("Events", "root://xrootd-cms.infn.it//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/7c887099-39ea-4238-b2f6-32d393c229e9.root");
    std::ifstream infile("test.list");
    std::vector<std::string> fileList;
    std::string line;

    while (std::getline(infile, line)) {
        if (!line.empty()) fileList.push_back(line);
    }

    ROOT::RDataFrame df("Events", fileList);

    //SOMEWHERE HERE NEED TO ASSIGN LEADING/SUBLEADING --> MUONS/JETS  

    //Muons
    auto Muon_total = df.Define("nTotal", "nMuon").Sum("nTotal");
    //auto Muon_leading = df.Define("leadingMuons", "nMuon[nTotal][0]").Sum("leadingMuons");
    //auto Muon_subleading = df.Define("subleadingMuons", "nMuon[nTotal][1]").Sum("subleadingMuons");
    
    auto Muon_global = df.Define("globalMuons", "Sum(Muon_isGlobal)").Sum("globalMuons");
    auto Muon_tracker = df.Define("trackerMuons", "Sum(Muon_isTracker)").Sum("trackerMuons");
    auto Muon_ptcut = df.Define("ptMuons", "Sum(Muon_pt > 29)").Sum("ptMuons");
    auto Muon_etacut = df.Define("etaMuons", "Sum(abs(Muon_eta) < 2.4)").Sum("etaMuons");
    auto Muon_tight = df.Define("tightMuons", "Sum(Muon_tightId)").Sum("tightMuons");
    auto Muon_pass = df.Define("passMuons", "Sum((Muon_isGlobal) && (Muon_isTracker) && (Muon_pt > 29) && (abs(Muon_eta) < 2.4) && (Muon_tightId))").Sum("passMuons");

    //auto Muons_leading = Muon_pass.Define("leadMuon_pt", "Muon_pt[passMuons][0]").Define("subleadMuon_pt", "Muon_pt[passMuons][1]");
    //auto Muon_bothpass = Muon_bothpass.Filter("Sum(passMuons >= 2");
    //auto Muons_goodcounter = Muons_leading.Count();

    //Jets
    auto Jet_total = df.Define("nAllJets", "nJet").Sum("nAllJets");
    auto Jet_ptcut = df.Define("ptJets", "Sum(Jet_pt > 35)").Sum("ptJets");
    auto Jet_etacut = df.Define("etaJets", "Sum(abs(Jet_eta) < 4.7)").Sum("etaJets");
    auto Jet_pass = df.Define("passJets", "Sum((Jet_pt > 35) && (abs(Jet_eta) < 4.7))").Sum("passJets");
    //auto DiJet_eta

    std::cout << "------------------------------ " << "\n";
    std::cout << "Muon Efficiency Information: " << "\n";
    std::cout << "Total muons: " << *Muon_total << "\n";
    std::cout << "Pass isGlobal: " << *Muon_global << "\n";
    std::cout << "Pass isTracker: " << *Muon_tracker << "\n";
    std::cout << "Pass pt > 29: " << *Muon_ptcut << "\n";
    std::cout << "Pass |eta| < 2.4: " << *Muon_etacut << "\n";
    std::cout << "Pass isTight: " << *Muon_tight << "\n";
    std::cout << "Pass all selection: " << *Muon_pass << "\n";
    //std::cout << "Leading Dimuon pair: " << *Muons_goodcounter << "\n";
    std::cout << "------------------------------ " << "\n";
    std::cout << "Jet Efficiency Information: " << "\n";
    std::cout << "Total jets: " << *Jet_total << "\n";
    std::cout << "Pass pt > 35: " << *Jet_ptcut << "\n";
    std::cout << "Pass |eta| < 4.7 " << *Jet_etacut << "\n";
    std::cout << "Pass all selection: " << *Jet_pass << "\n";
    std::cout << "------------------------------ " << "\n";

// Step 1: Only keep events with ≥2 muons
auto df2mu = df.Filter("nMuon >= 2");

// Step 2: Define tag and probe masks && 
// Define tag as a muon that passes tight selection
auto df_with_tags = df2mu.Define("tag_mask", "Muon_pt > 29 && abs(Muon_eta) < 2.4 && Muon_isGlobal && Muon_isTracker && Muon_tightId");
// Probes: all other muons in event
auto df_with_probes = df_with_tags.Define("probe_mask", "!tag_mask");

// Step 3: Build invariant mass of all tag-probe pairs
df_with_probes = df_with_probes.Define("mass_pass", R"cpp(
    std::vector<double> masses;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);

            double m = (tag + probe).M();
            if (m > 60 && m < 120 && Muon_tightId[j]) {
                masses.push_back(m);
            }
        }
    }
    return masses;
)cpp");

// Step 4: Define passing probes (tightId, trigger match, etc.)
df_with_probes = df_with_probes.Define("mass_fail", R"cpp(
        std::vector<double> masses;
        for (size_t i = 0; i < Muon_pt.size(); ++i) {
            if (!tag_mask[i]) continue;
            for (size_t j = 0; j < Muon_pt.size(); ++j) {
                if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
                TLorentzVector tag, probe;
                tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
                probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
                double m = (tag + probe).M();
                if (m > 60 && m < 120 && !Muon_tightId[j]) {
                    masses.push_back(m);
                }
            }
        }
        return masses;
    )cpp");

// Step 5: Plot histograms
    auto h_pass = df_with_probes.Histo1D({"mass_pass", "Passing Probes; m_{#mu#mu} [GeV]; Events", 60, 60, 120}, "mass_pass");
    auto h_fail = df_with_probes.Histo1D({"mass_fail", "Failing Probes; m_{#mu#mu} [GeV]; Events", 60, 60, 120}, "mass_fail");

    // Step 6: Draw
    TCanvas* c1 = new TCanvas("c1", "Tag and Probe Mass", 800, 600);
    gStyle->SetOptStat(0);
    h_pass->SetLineColor(kGreen+2);
    h_pass->SetLineWidth(2);
    h_fail->SetLineColor(kRed+1);
    h_fail->SetLineStyle(2);
    h_fail->SetLineWidth(2);
    h_pass->Draw("HIST");
    h_fail->Draw("HIST SAME");

    auto legend = new TLegend(0.6, 0.7, 0.89, 0.89);
    legend->AddEntry(h_pass.GetPtr(), "Passing probes", "l");
    legend->AddEntry(h_fail.GetPtr(), "Failing probes", "l");
    legend->Draw();

    c1->SaveAs("tag_probe_efficiency.png");

}
