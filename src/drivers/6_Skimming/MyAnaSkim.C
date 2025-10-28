//=============================
// MyAnaSkim.C — minimal, standalone NanoAOD skimmer
// Build:  c++ -O2 -std=gnu++17 MyAnaSkim.C -o my_skim $(root-config --cflags --libs)
// Run:    source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh && ./my_skim filelist.txt 20000   # (max events per file; use -1 for all)
//=============================

#define MyAnalysis_cxx

// Make sure 'string' is known before including the header
#include <string>
using std::string;

#include "MyAnalysis.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

using namespace std;

//--------- forward decls ----------
static string FileLocation(const string& initfile);

//--------- output root dir --------
//static string NewPath = "./skims";   // change if you like
//static string NewPath = "/eos/user/d/dblend/SkimmedFiles/12022/DATA/C";   // change if you like
//static string NewPath = "/eos/user/d/dblend/SkimmedFiles/12022/DATA/D";   // change if you like
static string NewPath = "/eos/user/d/dblend/SkimmedFiles/22023/DY/DY_005MCRun22023";   // change if you like

//--------- run stats --------------
static time_t startOfProgram = 0;
static unsigned int fcount = 0;
static unsigned int ievnt  = 0;
static int MaxEvent = -1;        // default cap (use -1 for all)

//--------- cutflow / bookkeeping --
static unsigned int CutEf[50] = {0};
static string eff_fname = "skim";

//--------- kept trees -------------
static TTree* newtree = nullptr;
static TFile* NewTreeFile = nullptr;

// ---- Skim knobs (loosened) ----
static const bool  REQUIRE_OS   = true;   // set false to also keep same-sign
// Vertex
static const float PV_CHI2_MAX  = 10.0f;  // was 4
static const float PV_Z_MAX     = 24.0f;  // was 5
// Muon
static const float MU_PT_MIN    = 5.0f;  // was 23
static const float MU_ETA_MAX   = 4.0f;   // was 2.4
static const float MU_DXY_MAX   = 1.0f;   // was 0.3
static const float MU_DXYBS_MAX = 1.0f;   // was 0.5
static const float MU_DZ_MAX    = 15.0f;  // was 5
static const float MU_ISO_MAX   = 0.50f;  // was 0.30 (tkRelIso)
static const float MU_SIP3D_MAX = 12.0f;   // was 5

//=============================
// MyAnalysis methods (match header signatures)
//=============================
void MyAnalysis::BookHistograms() {
  for (int i=0;i<50;i++) CutEf[i]=0;
}

void MyAnalysis::analyze_event(Long64_t /*entry*/)
{
  // --- PV quality ---
  //if (PV_npvsGood < 1)  return;
  //if (PV_ndof    <= 4)  return;
  //if (PV_chi2    >  PV_CHI2_MAX) return;
  //if (std::fabs(PV_z) > PV_Z_MAX) return;

  // --- muons ---
  if (nMuon < 2) return;

  int ipositive = 0, inegative = 0;
  for (int i=0; i<nMuon; ++i) {
    if (Muon_pt[i]            <= MU_PT_MIN)    continue;
    if (std::fabs(Muon_eta[i])>  MU_ETA_MAX)   continue;
    //if (std::fabs(Muon_dxy[i])>  MU_DXY_MAX)   continue;
    //if (std::fabs(Muon_dxybs[i])> MU_DXYBS_MAX) continue;
    //if (std::fabs(Muon_dz[i]) >  MU_DZ_MAX)    continue;
    //if (Muon_tkRelIso[i]      >  MU_ISO_MAX)   continue;
    //if (Muon_sip3d[i]         >  MU_SIP3D_MAX) continue;

    if (Muon_charge[i] > 0) ++ipositive;
    if (Muon_charge[i] < 0) ++inegative;
  }

  if (REQUIRE_OS) {
    if (ipositive < 1 || inegative < 1) return;
  } else {
    if (ipositive + inegative < 2) return;
  }

  // keep the event
  ++CutEf[0];
  ++ievnt;
  if (newtree) newtree->Fill();
}

void MyAnalysis::Loop()  // NOTE: header says 'virtual void Loop()'
{
  if (fChain == 0) return;

  const Long64_t nentries = fChain->GetEntriesFast();
  for (Long64_t jentry=0; jentry<nentries; ++jentry) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    fChain->GetEntry(jentry);
    analyze_event(ientry);

    if (MaxEvent != -1 && (int)ievnt >= MaxEvent) {
      std::cout << "\nReached requested event cap: " << ievnt << "\n";
      break; // stop this file; main() will decide whether to continue to next file
    }
  }
}

void MyAnalysis::EffPrint(string tag)  // header uses 'string'
{
  std::string filename = "eff_" + tag + ".txt";
  std::ofstream eff_file(filename.c_str());
  if (!eff_file) return;

  eff_file << fcount << " files skimmed\n";
  for (int i=0;i<50;i++){
    if (CutEf[i]!=0) eff_file << "Events after cut " << (i+1) << " : N=" << CutEf[i] << "\n";
  }
  eff_file.close();
}

//=============================
// helpers
//=============================
static string FileLocation(const string& initfile)
{
  const string SK = "Skimmed_";

  // split dir / file
  size_t lastSlash = initfile.find_last_of('/');
  string directory = (lastSlash==string::npos) ? "" : initfile.substr(0,lastSlash);
  string fileName  = (lastSlash==string::npos) ? initfile : initfile.substr(lastSlash+1);

  // extract year if present, else "misc"
  string year = "misc";
  size_t yearPos = directory.find("Run2_DYMC_RootFiles/20");
  if (yearPos != string::npos && yearPos + 24 <= directory.size()) {
    year = directory.substr(yearPos + 20, 4);
  }

  // ensure output dir exists
  //string outDir = NewPath + "/" + year + "/";
  string outDir = NewPath + "/" ;
  //gSystem->mkdir(outDir.c_str(), recursive=true);

  return outDir + SK + fileName;
}

//=============================
// main
//=============================
int main(int argc, char* argv[])
{
  time(&startOfProgram);

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " filelist.txt [maxEvents]\n";
    return 1;
  }

  const string list_file = argv[1];
  if (argc >= 3) MaxEvent = std::stoi(argv[2]);  // -1 for ALL

  // read input list
  std::ifstream list(list_file);
  if (!list) {
    std::cerr << "Cannot open list file: " << list_file << "\n";
    return 2;
  }

  std::vector<string> inputs;
  for (string line; std::getline(list, line); ) {
    if (!line.empty()) inputs.push_back(line);
  }
  if (inputs.empty()) {
    std::cerr << "No inputs found in " << list_file << "\n";
    return 3;
  }

  std::cout << (MaxEvent==-1 ? "Requested ALL events\n" :
                ("Requested " + std::to_string(MaxEvent) + " events\n"));

  // Create one analyzer instance (so counters accumulate across files)
  MyAnalysis Hmm;

  // process files
  for (const auto& infile : inputs) {
    ++fcount;

    // open input file & tree
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "Cannot open input: " << infile << "\n";
      if (f) { f->Close(); delete f; }
      continue;
    }

    TTree* tree = nullptr;
    f->GetObject("Events", tree);
    if (!tree) {
      std::cerr << "No 'Events' tree in: " << infile << "\n";
      f->Close(); delete f;
      continue;
    }

    // init reader on this tree
    Hmm.Init(tree);

    // ---- Slim the Events tree: disable all, then re-enable only what's needed
    tree->SetBranchStatus("*", 0);

    std::vector<const char*> keep = {
      // bookkeeping
      "run", "event", "luminosityBlock",

      // vertices
      "PV_*",

      // *** ADD: event weights & pileup (both namings) ***
      "Generator_*",
      "genWeight",
      "Pileup_*",
      "nTrueInt",

      // muons
      "nMuon",
      "Muon_*", 
      /*"Muon_pt","Muon_eta","Muon_phi","Muon_charge",
      "Muon_isGlobal","Muon_isTracker","Muon_tightId", "Muon_mediumId", 
      "Muon_dxy","Muon_dxybs","Muon_dz","Muon_sip3d","Muon_nTrackerLayers",
      "Muon_pfRelIso04_all","Muon_tkRelIso","Muon_bsConstrainedPt",*/

      // triggers & trig objects
      "HLT_IsoMu24",
      "nTrigObj","TrigObj_pt","TrigObj_eta","TrigObj_phi","TrigObj_id","TrigObj_filterBits",

      // jets
      "nJet",
      "Jet_*",
      /*"Jet_pt","Jet_eta","Jet_phi","Jet_mass","Jet_jetId",
      "Jet_btagDeepFlavB","Jet_area","Jet_chHEF","Jet_neHEF","Jet_chEmEF","Jet_neEmEF","Jet_muEF",
      "Jet_genJetIdx",*/

      // gen jets
      "nGenJet","GenJet_pt","GenJet_eta","GenJet_phi","GenJet_mass",

      // FSR photons
      "nFsrPhoton","FsrPhoton_pt","FsrPhoton_eta","FsrPhoton_phi",
      "FsrPhoton_muonIdx","FsrPhoton_dROverEt2","FsrPhoton_electronIdx",

      // electrons
      "nElectron","Electron_pt","Electron_eta","Electron_phi","Electron_cutBased",

      // Rho
      "Rho_fixedGridRhoFastjetAll",

      // MET
      "MET_pt","MET_phi",
    
      // Flags
      "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter",
      "Flag_BadPFMuonDzFilter", "Flag_eeBadScFilter", "Flag_ecalBadCalibFilter", "Flag_hfNoisyHitsFilter",
    };

    for (auto b : keep) tree->SetBranchStatus(b, 1);

    // open output file
    const string outPath = FileLocation(infile);
    std::cout << "Writing: " << outPath << "\n";
    NewTreeFile = TFile::Open(outPath.c_str(), "RECREATE");
    if (!NewTreeFile || NewTreeFile->IsZombie()) {
      std::cerr << "Cannot create output: " << outPath << "\n";
      if (NewTreeFile) { delete NewTreeFile; NewTreeFile = nullptr; }
      f->Close(); delete f;
      continue;
    }

    // clone ONLY enabled branches
    newtree = tree->CloneTree(0);
    newtree->SetDirectory(NewTreeFile);

    // loop & fill
    Hmm.Loop();

    // write & close
    newtree->Write();
    NewTreeFile->Close();
    delete NewTreeFile; NewTreeFile = nullptr;
    newtree = nullptr;

    // Stop across files if we hit the global cap
    if (MaxEvent != -1 && (int)ievnt >= MaxEvent) break;

    f->Close();
    delete f;
  }

  // print/save cutflow
  //Hmm.EffPrint(eff_fname);

  std::cout << "\nDone. Files skimmed: " << fcount << "  |  Events kept: " << ievnt << "\n";
  return 0;
}

