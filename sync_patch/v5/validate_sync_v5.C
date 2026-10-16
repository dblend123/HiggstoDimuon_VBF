#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <vector>

void validate_sync_v5(const char* filename) {
    TFile f(filename, "READ");
    if (f.IsZombie()) { std::cerr << "ERROR: cannot open file\n"; return; }
    TTree* t = nullptr;
    f.GetObject("SyncEvents", t);
    if (!t) { std::cerr << "ERROR: SyncEvents not found\n"; return; }
    const char* names[] = {"passJetVetoMap", "nBJetsLoose", "nBJetsMedium",
        "passBJetVeto", "jetPairIndex1", "jetPairIndex2", "jetPassQuality",
        "jetVetoMap", "nJet", "nSelectedJets", "passVBFCategoryApprox",
        "passVBFJetKinematics", "dijetMass", "dijetDEta", "dijetPt",
        "weightGenerator", "weightPileup", "weightMuonID", "weightMuonISO", "weightTrigger"};
    for (const char* name : names) {
        if (!t->GetBranch(name)) { std::cerr << "ERROR: missing branch " << name << "\n"; return; }
    }
    Bool_t passMap=false, passB=false, passApprox=false, passKinematics=false;
    Int_t nLoose=0, nMedium=0, idx1=-1, idx2=-1, nJet=0, nSelected=0;
    Float_t mjj=0, deta=0, ptjj=0;
    Float_t weightGen=0, weightPU=0, weightID=0, weightISO=0, weightTrigger=0;
    std::vector<unsigned char>* quality=nullptr;
    std::vector<unsigned char>* veto=nullptr;
    t->SetBranchAddress("passJetVetoMap", &passMap);
    t->SetBranchAddress("nBJetsLoose", &nLoose);
    t->SetBranchAddress("nBJetsMedium", &nMedium);
    t->SetBranchAddress("passBJetVeto", &passB);
    t->SetBranchAddress("jetPairIndex1", &idx1);
    t->SetBranchAddress("jetPairIndex2", &idx2);
    t->SetBranchAddress("jetPassQuality", &quality);
    t->SetBranchAddress("jetVetoMap", &veto);
    t->SetBranchAddress("nJet", &nJet);
    t->SetBranchAddress("nSelectedJets", &nSelected);
    t->SetBranchAddress("passVBFCategoryApprox", &passApprox);
    t->SetBranchAddress("passVBFJetKinematics", &passKinematics);
    t->SetBranchAddress("dijetMass", &mjj);
    t->SetBranchAddress("dijetDEta", &deta);
    t->SetBranchAddress("dijetPt", &ptjj);
    t->SetBranchAddress("weightGenerator", &weightGen);
    t->SetBranchAddress("weightPileup", &weightPU);
    t->SetBranchAddress("weightMuonID", &weightID);
    t->SetBranchAddress("weightMuonISO", &weightISO);
    t->SetBranchAddress("weightTrigger", &weightTrigger);
    Long64_t vetoMapFailures=0, vbfKinematics=0, approx=0;
    Long64_t badCount=0, badMap=0, badVeto=0, badIndices=0, badWeights=0, badCategory=0;
    for (Long64_t i=0; i<t->GetEntries(); ++i) {
        t->GetEntry(i);
        if (!passMap) ++vetoMapFailures;
        if (passKinematics) ++vbfKinematics;
        if (passApprox) ++approx;
        if (!quality || !veto || static_cast<int>(quality->size())!=nJet ||
            static_cast<int>(veto->size())!=nJet) { ++badCount; continue; }
        int qualified=0; bool anyVeto=false;
        for (int j=0; j<nJet; ++j) {
            if ((*quality)[j]) ++qualified;
            if ((*veto)[j]) anyVeto=true;
        }
        if (qualified!=nSelected) ++badCount;
        // For normal runs without --noJetVeto, any veto-map jet fails the flag.
        if (anyVeto && passMap) ++badMap;
        if (passB != (nMedium==0 && nLoose<=1) || nMedium>nLoose) ++badVeto;
        const bool hasPair=(idx1>=0 && idx2>=0);
        if (hasPair != (nSelected>=2) ||
            (hasPair && (idx1==idx2 || idx1>=nJet || idx2>=nJet ||
                         !(*quality)[idx1] || !(*quality)[idx2]))) ++badIndices;
        if (passApprox && (!passMap || !passB || !passKinematics)) ++badCategory;
        if (!std::isfinite(weightGen) || !std::isfinite(weightPU) ||
            !std::isfinite(weightID) || !std::isfinite(weightISO) ||
            !std::isfinite(weightTrigger)) ++badWeights;
    }
    std::cout << "Entries: " << t->GetEntries() << "\n"
              << "Jet veto map failures retained: " << vetoMapFailures << "\n"
              << "Pass VBF jet kinematics: " << vbfKinematics << "\n"
              << "Pass approximate VBF category: " << approx << "\n"
              << "Bad jet counts: " << badCount << "\n"
              << "Bad jet veto-map bookkeeping: " << badMap << "\n"
              << "Bad b-tag veto bookkeeping: " << badVeto << "\n"
              << "Bad jet-pair indices: " << badIndices << "\n"
              << "Bad approximate category bookkeeping: " << badCategory << "\n"
              << "Nonfinite checked MC weights: " << badWeights << "\n";
}
