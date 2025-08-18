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

void muonCounter() {
    std::ifstream infile("0DYBackground_Run22.list");
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
    }
