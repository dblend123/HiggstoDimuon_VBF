#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <iostream>
#include <vector>
#include <cmath>

void nEfficiencies() {
    TFile* f = TFile::Open("root://xrootd-cms.infn.it//store/data/Run2022G/Muon/NANOAOD/22Sep2023-v1/30000/7c887099-39ea-4238-b2f6-32d393c229e9.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open file!" << std::endl;
        return;
    }

    TTree* t = (TTree*)f->Get("Events");

    // Define pointers to vectors
    int nMuon;
    std::vector<bool>* Muon_isGlobal = nullptr;
    std::vector<bool>* Muon_isTracker = nullptr;
    std::vector<float>* Muon_pt = nullptr;
    std::vector<float>* Muon_eta = nullptr;
    std::vector<bool>* Muon_tightId = nullptr;

    // Set branch addresses
    t->SetBranchAddress("nMuon", &nMuon);
    t->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal);
    t->SetBranchAddress("Muon_isTracker", &Muon_isTracker);
    t->SetBranchAddress("Muon_pt", &Muon_pt);
    t->SetBranchAddress("Muon_eta", &Muon_eta);
    t->SetBranchAddress("Muon_tightId", &Muon_tightId);

    // Counters
    Long64_t total = 0, nGlobal = 0, nTracker = 0, nPtCut = 0, nEtaCut = 0, nTight = 0;

    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        t->GetEntry(i);
        for (int j = 0; j < nMuon; ++j) {
            ++total;
            if (!Muon_isGlobal->at(j)) continue;
            ++nGlobal;
            if (!Muon_isTracker->at(j)) continue;
            ++nTracker;
            if (Muon_pt->at(j) <= 29) continue;
            ++nPtCut;
            if (fabs(Muon_eta->at(j)) >= 2.4) continue;
            ++nEtaCut;
            if (!Muon_tightId->at(j)) continue;
            ++nTight;
        }
    }

    // Print results
    std::cout << "Total muons: " << total << std::endl;
    std::cout << "Pass isGlobal: " << nGlobal << " (" << 100.0 * nGlobal / total << "%)" << std::endl;
    std::cout << "Pass isTracker: " << nTracker << " (" << 100.0 * nTracker / total << "%)" << std::endl;
    std::cout << "Pass pt > 29: " << nPtCut << " (" << 100.0 * nPtCut / total << "%)" << std::endl;
    std::cout << "Pass |eta| < 2.4: " << nEtaCut << " (" << 100.0 * nEtaCut / total << "%)" << std::endl;
    std::cout << "Pass isTight: " << nTight << " (" << 100.0 * nTight / total << "%)" << std::endl;
}

