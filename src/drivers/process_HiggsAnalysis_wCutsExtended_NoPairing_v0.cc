//*******************************************************************
//                   HEADER FILES
//*******************************************************************

// Def Header Files
#include "Muon.h"
#include "Photon.h"
#include "FsrPhoton.h"
#include "Jet.h"
#include "GenJet.h"
#include "Vertex.h"
#include "Event.h"
#include "MET.h"
#include "Constants.h"
#include "Streamer.h"
#include "MetaHiggs.h"
#include "defs/Macros.h"
#include "L1.h"
#include "HLT.h"
#include "TrigObj.h"
#include "Electron.h"
#include "Rho.h"
#include "Flag.h"
#include "correction.h"
#include "KITMuonCorrections.h"


// Other Header Files
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "LumiReweightingStandAlone.h"
#include "HistogramSets.h"
#include <boost/program_options.hpp>
#include <signal.h>
#include <memory>
#include <algorithm>
#include <map>
#include <set>
#include <TVector2.h>
#include <cmath> 
#include <TRandom3.h>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <regex>
#include <utility>
#include <vector>
#include <TLeaf.h>
#include <random>   
#include <cstdint>


//*******************************************************************
//                   HISTOGRAM MACROS
//*******************************************************************

// -------------------------------------
//             VTX 
// -------------------------------------

#define FILL_VTX_HIST(set, x, weight1, weight2) \
    set.h_nPV_nom->Fill(x, weight1); \
    set.h_nPV_pu->Fill(x, weight2);

// -------------------------------------
//             MUONS 
// -------------------------------------

#define FILL_MUON_HIST(set, weight) \
    set.hDiMuonMass->Fill(corr_p4dimuon.M(), weight); \
    set.hDiMuonMass_Z->Fill(corr_p4dimuon.M(), weight); \
    set.hDiMuonMass_Higgs->Fill(corr_p4dimuon.M(), weight); \
    set.hDiMuonpt->Fill(corr_p4dimuon.Pt(), weight); \
    set.hDiMuoneta->Fill(corr_p4dimuon.Eta(), weight); \
    set.hDiMuonrapidity->Fill(corr_p4dimuon.Rapidity(), weight); \
    set.hMuon1pt->Fill(corr_p4m1.Pt(), weight); \
    set.hMuon2pt->Fill(corr_p4m2.Pt(), weight); \
    set.hMuon1eta->Fill(corr_p4m1.Eta(), weight); \
    set.hMuon2eta->Fill(corr_p4m2.Eta(), weight); \
    set.hMuon1PTdDM->Fill(PTdDM1, weight); \
    set.hMuon2PTdDM->Fill(PTdDM2, weight); \
    set.hCollinsSoper->Fill(phiCS, weight); \
    set.hcosCollinsSoper->Fill(cosThetaCS, weight); \

#define FILL_CORR_MUON_HIST(set, weight) FILL_MUON_HIST(set, weight)

// -------------------------------------
//             JETS 
// -------------------------------------

#define FILL_JET_HIST(set, weight) \
    set.hDiJetMass->Fill(dijetmass, weight); \
    set.hDiJetpt->Fill(dijetpt, weight); \
    set.hDiJetdeta->Fill(dijetdeta, weight); \
    set.hDiJetdphi->Fill(dijetdphi, weight); \
    set.hJet1pt->Fill(jet1pt, weight); \
    set.hJet2pt->Fill(jet2pt, weight); \
    set.hJet1eta->Fill(jet1eta, weight); \
    set.hJet2eta->Fill(jet2eta, weight); \
    set.hnJets->Fill(p4jets.size(), weight); \
    set.hdphi_mumu_jj->Fill(g_dphi_mumu_jj, weight); \
    set.hdeta_mumu_jj->Fill(g_deta_mumu_jj, weight); \

#define FILL_CORR_JET_HIST(set, weight) FILL_JET_HIST(set, weight)

#define FILL_BTAG_EFF_HIST(set, jetobj, wpL, wpM, weight)                              \
    do {                                                                                \
        const int flav = ((jetobj)._hadronFlavour == 5 ? 5 :                           \
                         ((jetobj)._hadronFlavour == 4 ? 4 : 0));                      \
        const float pt = (jetobj)._pt;                                                 \
        const float abseta = std::abs((jetobj)._eta);                                  \
        const float disc = (jetobj)._btagPNetB;                                        \
                                                                                        \
        const bool passL = (disc >= (wpL));                                            \
        const bool passM = (disc >= (wpM));                                            \
                                                                                        \
        if (flav == 5) {                                                                \
            (set).h_den_b_L->Fill(pt, abseta, weight);                                 \
            (set).h_den_b_M->Fill(pt, abseta, weight);                                 \
            if (passL) (set).h_num_b_L->Fill(pt, abseta, weight);                      \
            if (passM) (set).h_num_b_M->Fill(pt, abseta, weight);                      \
        } else if (flav == 4) {                                                         \
            (set).h_den_c_L->Fill(pt, abseta, weight);                                 \
            (set).h_den_c_M->Fill(pt, abseta, weight);                                 \
            if (passL) (set).h_num_c_L->Fill(pt, abseta, weight);                      \
            if (passM) (set).h_num_c_M->Fill(pt, abseta, weight);                      \
        } else {                                                                        \
            (set).h_den_l_L->Fill(pt, abseta, weight);                                 \
            (set).h_den_l_M->Fill(pt, abseta, weight);                                 \
            if (passL) (set).h_num_l_L->Fill(pt, abseta, weight);                      \
            if (passM) (set).h_num_l_M->Fill(pt, abseta, weight);                      \
        }                                                                               \
    } while (0)

//*******************************************************************
//                     CONSTANTS
//*******************************************************************
static TRandom3 g_rng(42); // or seed per event if you prefer

struct ZPtReweighter {
  std::unique_ptr<TFile> f;
  std::unique_ptr<TF1> f0;
  std::unique_ptr<TF1> f1;
  std::unique_ptr<TF1> f2;

  double xmin = 0.0;
  double xmax = 200.0;

  bool load(const std::string& rootPath,
            const std::string& dir0 = "corr_DM_Pre-Selection_Z0j",
            const std::string& dir1 = "corr_DM_Pre-Selection_Z1j",
            const std::string& dir2 = "corr_DM_Pre-Selection_Z2j",
            const std::string& fname = "f_combined",
            double xMin = 0.0, double xMax = 200.0)
  {
    xmin = xMin;
    xmax = xMax;

    f.reset(TFile::Open(rootPath.c_str(), "READ"));
    if (!f || f->IsZombie()) {
      std::cerr << "[ZPtRW] ERROR opening " << rootPath << "\n";
      return false;
    }

    auto getClone = [&](const std::string& d) -> std::unique_ptr<TF1> {
      const std::string path = d + "/" + fname;
      TF1* tf = dynamic_cast<TF1*>(f->Get(path.c_str()));
      if (!tf) {
        std::cerr << "[ZPtRW] ERROR missing TF1: " << path << " in " << rootPath << "\n";
        return nullptr;
      }
      // Clone so we can safely keep it after file operations; detach from directory.
      TF1* c = dynamic_cast<TF1*>(tf->Clone((d + "_" + fname + "_clone").c_str()));
      return std::unique_ptr<TF1>(c);
    };

    f0 = getClone(dir0);
    f1 = getClone(dir1);
    f2 = getClone(dir2);

    return (bool)f0 && (bool)f1 && (bool)f2;
  }

  double eval(double zpt, int njets) const {
    if (!f0 || !f1 || !f2) return 1.0;

    // clamp to fit range so you don’t extrapolate wildly
    const double x = std::clamp(zpt, xmin, xmax);

    const TF1* ff = (njets <= 0 ? f0.get() : (njets == 1 ? f1.get() : f2.get()));
    double w = ff->Eval(x);

    // safety: no negative/NaN weights
    if (!std::isfinite(w) || w <= 0.0) w = 1.0;
    return w;
  }
};

// Service Globals
std::string __inputfilename;
std::string __outputfilename;
std::string const NTUPLEMAKER_NAME =  "";

bool __isMC;
bool __useKIT;
bool __useRochester;
bool __noMuCorr;
bool g_disableJetCriteria    = false;  // if true, skip jet kinematic/ID/muon-cleaning cuts
bool g_disableJetJERC        = false;  // if true, ignore JERC JSONs and use raw jet pT
bool g_disableJetVeto        = false;  // if true, ignore jet veto maps
bool g_disableBtagSF         = false;  // if true, ignore b-tag SFs (set them to 1)
bool g_disableMuonScaleSmear = false;  // if true, ignore Scale/Smear json
bool g_disableMuonSF         = false;  // if true, ignore muon SFs
bool g_disableFSRRecovery    = false;  // if true, skip FSR photon recovery
bool g_disablePUWeights      = false;  // if true, PU weight = 1
bool g_disableFlagFilters    = false;  // if true, skip event Flag (MET) filters
bool g_disableZPtRW          = false;

std::string g_zptRWFile      = "None";
static ZPtReweighter g_zptRW;
std::atomic<bool> __continueRunning = true;

static bool g_disableBeamConstraint = false;  // if true, skip Beamspot Constrained Pt 
static inline bool useBeamConstraintPt(const analysis::core::Muon& mu);
static inline float baseMuonPtForCorrections(const analysis::core::Muon& mu);

// MVA tree
TTree* gMvaTree = nullptr;

float g_z             = -99.f;
float g_ptBalance     = -99.f;

float g_mu1_pt_over_m = -99.f;
float g_mu2_pt_over_m = -99.f;
float g_mu1_eta       = -99.f;
float g_mu2_eta       = -99.f;
float g_jet1_eta      = -99.f;
float g_m_mumu        = -99.f;
float g_deta_jj       = -99.f;
float g_dphi_jj       = -99.f;
int   g_njets         = -99.;
float g_dphi_mumu_jj  = -99.f;
float g_deta_mumu_jj  = -99.f;
float g_phi_CS        = -99.f;
float g_costheta_CS   = -99.f;

float g_y_mumu        = -99.f;
float g_y_j1          = -99.f;
float g_y_j2          = -99.f;

float g_pt_mumu       = -99.f;
float g_pt_j1         = -99.f;
float g_pt_j2         = -99.f;   

float g_weight        = 1.f;

// Weights Config
struct SampleConfig {
    double xsec = 0.0;        // pb
    double sumw = 0.0;        // genEventSumw
};

std::string g_weightsFile = "None";  // set via --weights
std::string g_sampleName  = "DATA";  // set via --sample

// --- Jet Debug Counters ---
static unsigned long long jets_in               = 0;
static unsigned long long jets_pass             = 0;
static unsigned long long jets_after_vetoMap    = 0;
static unsigned long long jets_after_id         = 0;
static unsigned long long jets_after_eta        = 0;
static unsigned long long jets_after_dR1        = 0;
static unsigned long long jets_after_dR2        = 0;
static unsigned long long jets_after_pt_default = 0;
static unsigned long long jets_after_pt_hf      = 0;
static unsigned long long jets_after_pt_he      = 0;
static unsigned long long jets_after_pt         = 0;

// Era handling
std::string g_era = "12022";  // default section (can be overridden via --era)
std::unordered_map<std::string, std::unordered_map<std::string, double>> g_sectionLumis;
std::unordered_map<std::string, std::unordered_map<std::string, SampleConfig>> g_sectionSamples;

double g_selected_lumi = 0.0; 
std::string g_selected_lumi_key = "unset";

// Kinematic Cuts  
enum MuonType {kLead=0, kSubLead=1};  // Muon Id --> lead = 0, sublead =1
double _leadmuonPt          = 26.;  // +2 GeV of HLT trigger 
double _subleadmuonPt       = 20.;
double _leadmuonEta         = 2.4;
double _subleadmuonEta      = 2.4;
double _leadmuonIso         = 0.25;
double _subleadmuonIso      = 0.25;
double _leadJetPt           = 35.; //Kind of unsure --> Purdue has 50 GeV w/out Loose PU ID
double _subleadJetPt        = 25.;
double _leadJetEta          = 4.7;
double _subleadJetEta       = 4.7;
double _dijetMass_VBF       = 400;
double _dijetdEta_VBF       = 2.5;

// Namespace 
namespace po = boost::program_options;
using namespace analysis::core;
using namespace analysis::dimuon;
using namespace analysis::processing;
using namespace correction;

//*******************************************************************
//                   CORRECTIONLIB LIBRARY FILES
//*******************************************************************

// -------------------------------------
//             MUONS 
// -------------------------------------

// Muon SF && Momentum Corrections 
std::shared_ptr<correction::CorrectionSet> corrections               = nullptr;     // muon SFs
std::shared_ptr<correction::CorrectionSet> muoncorrections           = nullptr;     // Rochester
std::shared_ptr<correction::CorrectionSet> kitCorrections            = nullptr;     // KIT
std::shared_ptr<correction::CorrectionSet> muonScaleSmearCorrections = nullptr;     // Smearing

// ---------- Era-aware config (edit defaults if your layout changes) ----------
std::string g_muDir             = "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/MUO";
std::string g_muSFJSON          = "auto";   // "auto" => muon_Z<era>.json
std::string g_kitJSON           = "auto";
std::string g_muSFKeyID         = "NUM_MediumID_DEN_TrackerMuons";
std::string g_muSFKeyISO        = "NUM_LoosePFIso_DEN_MediumID";
std::string g_muSFKeyTRG        = "NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium"; 
std::string g_rochJSON          = "auto";   // "auto" => RoccoR<era>.json
static std::string g_rochKeyDT  = "kScaleDT";
static std::string g_rochKeyMC  = "kScaleMC";

static std::string muSFFileForEra(const std::string& era) {
    if (g_muSFJSON != "auto") return g_muSFJSON;
    // Example naming: muon_Z12022.json, muon_Z22022.json, muon_Z12023.json, muon_Z22023.json, muon_Z2024.json
    return g_muDir + "/muon_Z" + era + ".json";
}

static std::string rochFileForEra(const std::string& era) {
    if (g_rochJSON != "auto") return g_rochJSON;
    // Example naming: RoccoR12022.json, RoccoR22022.json, RoccoR12023.json, RoccoR22023.json, RoccoR2024.json
    return g_muDir + "/RoccoR" + era + ".json";
}

static std::string kitFileForEra(const std::string& era) {
    if (g_kitJSON != "auto") return g_kitJSON;
    // Example naming: KIT12022.json, KIT22022.json, KIT12023.json, KIT22023.json, KIT2024.json
    return g_muDir + "/KIT" + era + ".json";
}

static std::string muonScaleSmearFileForEra(const std::string& era) {
    // Example naming: muon_scalesmearing12022.json, muon_scalesmearing22022.json, ...
    std::string baseDir = g_muDir;  // same as for other corrections
    return baseDir + "/muon_scalesmearing" + era + ".json";
}

void loadMuonSF() {
    if (g_disableMuonSF) {
        std::cerr << "[MUON] Skipping muon SF JSON loading (--noMuonSF); muon SFs will be 1.\n";
        corrections = nullptr;
        return;
    }
    try {
        const std::string sfPath = muSFFileForEra(g_era);

        if (!corrections)
            corrections = correction::CorrectionSet::from_file(sfPath);

        std::cerr << "[MUON] Loaded SFs from " << sfPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[MUON SF load ERROR] " << e.what() << std::endl;
        corrections = nullptr;
    }
}

void loadRochesterCorrections() {
    if (__noMuCorr || !__useRochester) {
        std::cerr << "[MUON] Skipping Rochester corrections (mode != roch).\n";
        muoncorrections = nullptr;
        return;
    }
    try {
        const std::string rochPath = rochFileForEra(g_era);

        if (!muoncorrections)
            muoncorrections = correction::CorrectionSet::from_file(rochPath);

        std::cerr << "[MUON] Loaded Rochester from " << rochPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[Rochester load ERROR] " << e.what() << std::endl;
        muoncorrections = nullptr;
    }

    // --- Verify Rochester keys (probe + fallback) ---
    if (!muoncorrections) return;

    try {
        (void) muoncorrections->at(g_rochKeyDT)->evaluate({+1, 30.0, 0.0, 0.0, 0, 0});
        std::cerr << "[MUON] Rochester DT key: " << g_rochKeyDT << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[MUON] WARN: '" << g_rochKeyDT
                  << "' failed probe (" << e.what() << "). Trying fallback 'kScaleDT'.\n";
        try {
            (void) muoncorrections->at("kScaleDT")->evaluate({+1, 30.0, 0.0, 0.0, 0, 0});
            g_rochKeyDT = "kScaleDT";
            std::cerr << "[MUON] Rochester DT key resolved to 'kScaleDT'.\n";
        } catch (const std::exception& e2) {
            std::cerr << "[MUON] ERROR: No valid Rochester DT key found (" << e2.what()
                      << "). DT scale will be 1.0.\n";
            g_rochKeyDT.clear();
        }
    }

    try {
        (void) muoncorrections->at(g_rochKeyMC)->evaluate({+1, 30.0, 0.0, 0.0, 0, 0});
        std::cerr << "[MUON] Rochester MC key: " << g_rochKeyMC << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[MUON] WARN: '" << g_rochKeyMC
                  << "' failed probe (" << e.what() << "). Trying fallback 'kScaleMC'.\n";
        try {
            (void) muoncorrections->at("kScaleMC")->evaluate({+1, 30.0, 0.0, 0.0, 0, 0});
            g_rochKeyMC = "kScaleMC";
            std::cerr << "[MUON] Rochester MC key resolved to 'kScaleMC'.\n";
        } catch (const std::exception& e2) {
            std::cerr << "[MUON] ERROR: No valid Rochester MC key found (" << e2.what()
                      << "). MC scale will be 1.0.\n";
            g_rochKeyMC.clear();
        }
    }
}

void loadKITCorrections() {
    if (__noMuCorr) {
        std::cerr << "[MUON] Skipping KIT corrections (--muCorr none).\n";
        muoncorrections = nullptr;
        return;
    }
    try {
        const std::string kitPath = kitFileForEra(g_era);
        
        // Load (or reload if null). If you want forced reload per era, you can reset to nullptr first.
        if (!kitCorrections)    kitCorrections = correction::CorrectionSet::from_file(kitPath);
        
        std::cerr << "[MUON] Loaded KIT params from " << kitPath << "\n";
        
        (void) kitCorrections->at("cb_params");
        // optional: KIT sets also provide "poly_params" in 2023 and BPix jsons
        // (void) kitCorrections->at("poly_params");
    } catch (const std::exception& e) {
        std::cerr << "[KIT load ERROR] " << e.what() << std::endl;
        kitCorrections = nullptr;
    }
}

void loadMuonScaleSmearCorrections() {
   if (g_disableMuonScaleSmear) { 
        std::cerr << "[SCALESMEARING] Skipping ScaleSmearing JSON loading (noScaleSmear).\n";
        muonScaleSmearCorrections = nullptr;
        return;
    }
    try {
        const std::string path = muonScaleSmearFileForEra(g_era);
        muonScaleSmearCorrections = correction::CorrectionSet::from_file(path);
        std::cerr << "[MUON] Loaded Scale/Smearing corrections from " << path << "\n";
        // Optionally validate
        (void) muonScaleSmearCorrections->at("cb_params");
    } catch (const std::exception& e) {
        std::cerr << "[MUON SCALE/SMEAR load ERROR] " << e.what() << std::endl;
        muonScaleSmearCorrections = nullptr;
    }
}

// -------------------------------------
//           Z-PT REWEIGHTING 
// -------------------------------------

static void initZPtRW() {
    if (g_disableZPtRW) return;

    if (g_zptRWFile == "None") {
        std::cerr << "[ZPtRW] ERROR: enabled but --zptRWFile=None. Disabling.\n";
        g_disableZPtRW = true;
        return;
    }

    if (!g_zptRW.load(g_zptRWFile, "corr_DM_Pre-Selection_Z0j", "corr_DM_Pre-Selection_Z1j", "corr_DM_Pre-Selection_Z2j", "f_combined", 0.0, 200.0)) {
        std::cerr << "[ZPtRW] ERROR: failed to load; disabling\n";
        g_disableZPtRW = true;
        return;
    }

    std::cerr << "[ZPtRW] Loaded ZpT weights from " << g_zptRWFile << "\n";
}

// -------------------------------------
//             JETS 
// -------------------------------------

// Jet Energy Corrections/Resolutions (JEC/JER) 
std::shared_ptr<correction::CorrectionSet> jetCorrections   = nullptr;
std::shared_ptr<correction::CorrectionSet> jetVetoMaps      = nullptr;
std::shared_ptr<correction::CorrectionSet> jetIdCorrections = nullptr;
std::shared_ptr<correction::CorrectionSet> btagCorrections  = nullptr;

// ---------- Era-aware config & helpers ----------
std::string g_jmeDir        = "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/JME";
std::string g_jercJSON      = "auto";       // "auto" => jet_jerc<era>.json
std::string g_jetVetoJSON   = "auto";       // "auto" => jetvetomaps<era>.json
std::string g_jetIdJSON     = "auto";       // "auto" => jetid<era>.json
std::string g_btagJSON      = "auto";       // "auto" => btagging<era>.json
std::string g_btagAlgo      = "deepJet";    // "deepJet", "particleNet", "robustParticleTransformer"

struct JecEraConfig {
  std::string prefixBase;                   // base for JEC keys, we’ll append _MC or _<RunTag>_DATA
  std::string dataRunTag;                   // e.g. RunCD, RunE, RunF, RunG, RunCv123, RunCv4, RunD
  std::string jerScaleKey;                  // exact JER SF key
  std::string jerResKey;                    // exact JER resolution key
  std::string vetoMapKey;                   // veto map key (if available for this era/tag); leave "" to skip
};

static inline const char* pnetCorrName() {
  return (g_era == "12022" || g_era == "22022") ? "particleNet_comb" : "ParticleNet";
}

// Optional: try to infer 2022EE/2023 data run from the input file name
static std::string inferDataRunTag(const std::string& era, const std::string& inputFile) {
  if (era == "22022") {
    // accept both "RunE/F/G" and your list naming "E_...list", "F_...list", "G_...list"
    if (inputFile.find("RunE") != std::string::npos ||
        inputFile.find("/E_")  != std::string::npos ||
         inputFile.find("_E_")  != std::string::npos) return "RunE";
    if (inputFile.find("RunF") != std::string::npos ||
        inputFile.find("/F_")  != std::string::npos ||
        inputFile.find("_F_")  != std::string::npos) return "RunF";
    if (inputFile.find("RunG") != std::string::npos ||
        inputFile.find("/G_")  != std::string::npos ||
        inputFile.find("_G_")  != std::string::npos) return "RunG";
    return "RunE"; // safe default
  }
  if (era == "12023") {
    if (inputFile.find("RunCv4")   != std::string::npos) return "RunCv4";
    if (inputFile.find("RunCv123") != std::string::npos) return "RunCv123";
    return "RunCv123";
  }
  if (era == "22023") return "RunD";
  return "RunCD";
}

static JecEraConfig jecConfigForEra(const std::string& era, const std::string& inputFile) {
    if (era == "12022") {
        return {
            "Summer22_22Sep2023_V3",                               // prefix base
            "RunCD",                                               // DATA tag
            "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",   // JER SF
            "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",  // JER Res
            "Summer22_23Sep2023_RunCD_V1"                          // veto map
        };
    }
    if (era == "22022") {
        const std::string runTag = inferDataRunTag(era, inputFile);
        return {
            "Summer22EE_22Sep2023_V3",
            runTag,
            "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",
            "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
            "Summer22EE_23Sep2023_RunEFG_V1"
        };
    }
    if (era == "12023") {
        const std::string runTag = inferDataRunTag(era, inputFile);
        return {
            "Summer23Prompt23_V2",
            runTag,
            "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi",
            "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi",
            "Summer23Prompt23_RunC_V1"
        };
    }
    if (era == "22023") {
        return {
            "Summer23BPixPrompt23_V3",
            "RunD",
            "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
            "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
            "Summer23BPixPrompt23_RunD_V1" 
        };
    }
    if (era == "2024") {
        return {
            "Summer24Prompt24_V1",
            "RunBCDEFGHI",     
            "Summer24Prompt24_JRV1_MC_ScaleFactor_AK4PFPuppi",
            "Summer24Prompt24_JRV1_MC_PtResolution_AK4PFPuppi",
            "Summer24Prompt24_RunBCDEFGHI_V1"
        };
    }
    return {
        "Summer22_22Sep2023_V2","RunCD",
        "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",
        "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
        "Summer22_23Sep2023_RunCD_V1"
    };
}

static std::string jercFileForEra(const std::string& era) {
  if (g_jercJSON != "auto") return g_jercJSON;
  return g_jmeDir + "/jet_jerc" + era + ".json";
}
static std::string vetoFileForEra(const std::string& era) {
  if (g_jetVetoJSON != "auto") return g_jetVetoJSON;
  return g_jmeDir + "/jetvetomaps" + era + ".json";
}
static std::string jetIdFileForEra(const std::string& era) {
  if (g_jetIdJSON != "auto") return g_jetIdJSON;
  return g_jmeDir + "/jetid" + era + ".json";
}

static std::string btagFileForEra(const std::string& era) {
    if (g_btagJSON != "auto") return g_btagJSON;
    return g_jmeDir + "/btagging" + era + ".json";
}

void loadJetCorrections() {
    if (g_disableJetJERC) {
        std::cerr << "[Jet] JERC disabled via --noJERC; using raw jet pT (no JEC/JER).\n";
        jetCorrections = nullptr;
    } else {
        try {
        if (!jetCorrections)
        jetCorrections = correction::CorrectionSet::from_file(jercFileForEra(g_era));
        } catch (const std::exception& e) {
            std::cerr << "[Jet JERC load ERROR] " << e.what() << std::endl;
            jetCorrections = nullptr;
        }
        }
    if (g_disableJetVeto) {
        std::cerr << "[Jet] Veto maps disabled via --noJetVeto; no jets will be vetoed by maps.\n";
        jetVetoMaps = nullptr;
    } else {
        try {
        if (!jetVetoMaps)
        jetVetoMaps = correction::CorrectionSet::from_file(vetoFileForEra(g_era));
        } catch (const std::exception& e) {
            std::cerr << "[Jet veto map load ERROR] " << e.what() << std::endl;
            jetVetoMaps = nullptr;
        }
        }
  try {
    if (!jetIdCorrections)
      jetIdCorrections = correction::CorrectionSet::from_file(jetIdFileForEra(g_era));
  } catch (const std::exception& e) {
    std::cerr << "[JetID load ERROR] " << e.what() << std::endl;
    jetIdCorrections = nullptr;
  }
  if (g_disableBtagSF) {
    std::cerr << "[BTAG] Skipping b-tag SF JSON loading (--noBtagSF).\n";
    btagCorrections = nullptr;
  } else {
  try {
    if (!btagCorrections)
      btagCorrections = correction::CorrectionSet::from_file(btagFileForEra(g_era));

    if (btagCorrections) {
      auto has = [&](const char* name){
        try { (void)btagCorrections->at(name); return true; } catch (...) { return false; }
      };
      // Probe the names we actually use
      const bool has_wp      = has("particleNet_wp_values") || has("ParticleNet_wp_values");
      const bool has_comb    = has("particleNet_comb")      || has("ParticleNet_comb");
      const bool has_shape   = has("particleNet_shape")     || has("ParticleNet_shape");
      const bool has_light   = has("particleNet_light")     || has("ParticleNet_light");

      std::cerr << "[BTAG] Loaded from " << btagFileForEra(g_era)
                << " | keys: "
                << (has_wp    ? "wp "    : "")
                << (has_comb  ? "comb "  : "")
                << (has_shape ? "shape " : "")
                << (has_light ? "light " : "")
                << "\n";
    }
  } catch (std::exception& e) {
    std::cerr << "[Btag load ERROR] " << e.what() << std::endl;
    btagCorrections = nullptr;
  }
  }
}


// -------------------------------------
//             PU WEIGHTS 
// -------------------------------------

std::shared_ptr<correction::CorrectionSet> pileupCorrections = nullptr;
std::shared_ptr<const correction::Correction> pu_corrections = nullptr;

// Optional CLI-configurable bits (declare these near your other globals if you prefer)
std::string g_puDir   = "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/puWeights";
std::string g_puJSON  = "auto";  // "auto" => use puWeights<era>.json
std::string g_puKey   = "auto";  // "auto" => pick key by era via map below

static std::string puFileForEra(const std::string& era) {
  if (g_puJSON != "auto") return g_puJSON;
  // convention: puWeights<era>.json sitting in g_puDir
  return g_puDir + "/puWeights" + era + ".json";
}

// Fill in the exact names from your JSONs if they differ.
static std::string puKeyForEra(const std::string& era) {
  if (g_puKey != "auto") return g_puKey;
  // Known keys (edit if your JSONs use slightly different strings)
  if (era == "12022") return "Collisions2022_355100_357900_eraBCD_GoldenJson";
  if (era == "22022") return "Collisions2022_359022_362760_eraEFG_GoldenJson";
  if (era == "12023") return "Collisions2023_366403_369802_eraBC_GoldenJson";
  if (era == "22023") return "Collisions2023_369803_370790_eraD_GoldenJson";
  // fallback: common names people use in PU JSONs
  return "weights";  // try a generic name if your files have only one correction
}

void loadPileupCorrections() {
  if (g_disablePUWeights) {
    std::cerr << "[PU] Skipping pileup JSON loading (--noPU); using pileupWeight=1.\n";
    pileupCorrections = nullptr;
    pu_corrections    = nullptr;
    return;
  }

  try {
    const std::string puPath = puFileForEra(g_era);
    const std::string puKey  = puKeyForEra(g_era);

    pileupCorrections = correction::CorrectionSet::from_file(puPath);
    pu_corrections    = pileupCorrections->at(puKey);

    std::cerr << "[PU] Loaded " << puKey << " from " << puPath << "\n";
  } catch (const std::exception& e) {
    std::cerr << "[PU load ERROR] " << e.what() << std::endl;
    std::cerr << "[PU load ERROR] era=" << g_era
              << " file=" << puFileForEra(g_era)
              << " key="  << puKeyForEra(g_era) << "\n";
    pileupCorrections = nullptr;
    pu_corrections    = nullptr;
  }
}

// -------------------------------------
//             GOLDEN JSON
// -------------------------------------

static std::map<unsigned, std::vector<std::pair<unsigned,unsigned>>> g_lumiMask;
std::string g_goldenDir  = "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/Golden_json";
std::string g_goldenJSON = "auto"; 

static void skipWS(const std::string& s, size_t& i) {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
}

static bool consume(const std::string& s, size_t& i, char c) {
    skipWS(s, i);
    if (i < s.size() && s[i] == c) { ++i; return true; }
    return false;
}

static unsigned parseUInt(const std::string& s, size_t& i) {
    skipWS(s, i);
    if (i >= s.size() || !std::isdigit(static_cast<unsigned char>(s[i]))) {
        throw std::runtime_error("Expected unsigned integer while parsing golden JSON");
    }
    unsigned val = 0;
    while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) {
        val = val * 10 + (s[i] - '0');
        ++i;
    }
    return val;
}

static std::string parseQuotedString(const std::string& s, size_t& i) {
    skipWS(s, i);
    if (i >= s.size() || s[i] != '"') {
        throw std::runtime_error("Expected quoted string while parsing golden JSON");
    }
    ++i; // skip opening quote
    std::string out;
    while (i < s.size() && s[i] != '"') {
        out.push_back(s[i]);
        ++i;
    }
    if (i >= s.size()) {
        throw std::runtime_error("Unterminated quoted string in golden JSON");
    }
    ++i; // skip closing quote
    return out;
}

static std::string goldenFileForEra(const std::string& era) {
  if (g_goldenJSON != "auto") return g_goldenJSON;

  if (era == "12022") return g_goldenDir + "/Cert_Collisions2022_355100_362760_Golden.json";
  if (era == "22022") return g_goldenDir + "/Cert_Collisions2022_355100_362760_Golden.json";
  if (era == "12023") return g_goldenDir + "/Cert_Collisions2023_366442_370790_Golden.json";
  if (era == "22023") return g_goldenDir + "/Cert_Collisions2023_366442_370790_Golden.json";
  if (era == "2024")  return g_goldenDir + "/Cert_Collisions2024_378981_386951_Golden.json";

  // Fallback
  return g_goldenDir + "/Cert_Collisions2023_13p6TeV_Golden.json";
}

static void loadGoldenJsonCorrections() {
    g_lumiMask.clear();

    const std::string path = goldenFileForEra(g_era);
    if (path.empty()) {
        std::cerr << "[GOLDEN] No golden JSON configured for era " << g_era << "\n";
        return;
    }

    std::ifstream in(path);
    if (!in) {
        std::cerr << "[GOLDEN] Failed to open file: " << path << "\n";
        return;
    }

    std::string text((std::istreambuf_iterator<char>(in)),
                      std::istreambuf_iterator<char>());

    std::cerr << "[GOLDEN] Raw file size = " << text.size() << " bytes\n";

    try {
        size_t i = 0;
        skipWS(text, i);

        if (!consume(text, i, '{')) {
            throw std::runtime_error("Expected opening '{'");
        }

        while (true) {
            skipWS(text, i);
            if (i < text.size() && text[i] == '}') {
                ++i;
                break;
            }

            std::string runStr = parseQuotedString(text, i);
            unsigned run = static_cast<unsigned>(std::stoul(runStr));

            if (!consume(text, i, ':')) {
                throw std::runtime_error("Expected ':' after run number");
            }
            if (!consume(text, i, '[')) {
                throw std::runtime_error("Expected '[' starting lumi-range list");
            }

            std::vector<std::pair<unsigned,unsigned>> ranges;

            while (true) {
                skipWS(text, i);
                if (i < text.size() && text[i] == ']') {
                    ++i; // end of this run's lumi list
                    break;
                }

                if (!consume(text, i, '[')) {
                    throw std::runtime_error("Expected '[' starting lumi pair");
                }

                unsigned a = parseUInt(text, i);

                if (!consume(text, i, ',')) {
                    throw std::runtime_error("Expected ',' inside lumi pair");
                }

                unsigned b = parseUInt(text, i);

                if (!consume(text, i, ']')) {
                    throw std::runtime_error("Expected ']' ending lumi pair");
                }

                ranges.emplace_back(a, b);

                skipWS(text, i);
                if (i < text.size() && text[i] == ',') {
                    ++i; // next lumi pair
                    continue;
                }
            }

            if (!ranges.empty()) {
                g_lumiMask.emplace(run, std::move(ranges));
            }

            skipWS(text, i);
            if (i < text.size() && text[i] == ',') {
                ++i; // next run block
                continue;
            }
        }

        std::cerr << "[GOLDEN] Loaded " << g_lumiMask.size()
                  << " runs from " << path << "\n";

    } catch (const std::exception& e) {
        std::cerr << "[GOLDEN] Parse error: " << e.what() << "\n";
        std::cerr << "[GOLDEN] File: " << path << "\n";
        g_lumiMask.clear();
    }
}

// Fast membership test (same idea as your PU accessors)
static inline bool passGolden(unsigned run, unsigned lumi) {
  auto it = g_lumiMask.find(run);
  if (it == g_lumiMask.end()) return false;
  for (const auto& rg : it->second) {
    if (lumi >= rg.first && lumi <= rg.second) return true;
  }
  return false;
}

//*******************************************************************
//                      IMPLEMENTING CORRECTIONS
//*******************************************************************

// -------------------------------------
//          MUON PT (ROCHESTER) 
// -------------------------------------

float getMuonRochester(const Muon& mu, bool isMC)
{
    static int errorCount = 0;
    if (!muoncorrections) {
        if (errorCount < 3) {
            std::cerr << "[ERROR] Rochester corrections not loaded!\n";
        }
        ++errorCount;
        return 1.0f;
    }

    // --- Sanitize inputs ---
    int   charge = (mu._charge > 0) ? +1 : (mu._charge < 0 ? -1 : +1); // zero -> +1
    float pt     = baseMuonPtForCorrections(mu);
    float eta    = std::isfinite(mu._eta) ? mu._eta : 0.0f;
    float phi    = std::isfinite(mu._phi) ? mu._phi : 0.0f;

    // keep pt in a safe physical range
    if (!std::isfinite(pt) || pt <= 0.f) return 1.0f;

    // wrap phi into [-pi, pi] if it drifted
    if (phi >  M_PI) phi = std::fmod(phi + M_PI, 2*M_PI) - M_PI;
    if (phi < -M_PI) phi = std::fmod(phi - M_PI, 2*M_PI) + M_PI;

    int s = 0;
    int m = 0;

    try {
        // Primary: numeric leafs (kScaleDT/kScaleMC)
        if (!isMC) {
            static bool printed = false;
            if (!printed) {
                std::cerr << "[ROCH] Applying Rochester correction for first muon event" << std::endl;
                printed = true;
            } 

            if (g_rochKeyDT.empty()) return 1.0f;
            return muoncorrections->at(g_rochKeyDT)->evaluate({charge, pt, eta, phi, s, m});
        } else {
            if (g_rochKeyMC.empty()) return 1.0f;
            return muoncorrections->at(g_rochKeyMC)->evaluate({charge, pt, eta, phi, s, m});
        }
    } catch (const std::exception& e) {
        // Fallback: some JSONs expose a category node requiring a string first (rare).
        try {
            const char* tag = isMC ? "MC" : "DT";
            return muoncorrections->at("RoccoR")->evaluate({std::string(tag), charge, pt, eta, phi, s, m});
        } catch (const std::exception& e2) {
            static int printed = 0;
            if (printed++ < 5) {
                std::cerr << "[ERROR] Rochester failed. primary=" << e.what()
                          << " fallback=" << e2.what()
                          << " (charge=" << charge << ", pt=" << pt
                          << ", eta=" << eta << ", phi=" << phi << ")\n";
            }
            return 1.0f;
        }
    }
}

// -------------------------------------
//          MUON PT (KIT) 
// -------------------------------------

// Contents can be found in KIT.h 

// -------------------------------------
//          MUON SMEARING 
// -------------------------------------

static float getNTrackerLayersOrDefault(const Muon& mu) {
    // use your stored number of tracker layers if available
    if (mu._nTrackerLayers > 0 && std::isfinite(mu._nTrackerLayers))
        return mu._nTrackerLayers;
    return 12.0f; // fallback
}

float getMuonScaleSmear(const Muon& mu, bool isMC) {
    // Smear MC only; keep data untouched (you can change this if your JSON includes data terms)
    if (!isMC) return 1.0f;

    if (!muonScaleSmearCorrections) return 1.0f;

    // Pull CB params (alpha, n, sigma, tail) from the JSONs:
    // inputs: abseta, nTrackerLayers, parameter(0..3)
    // muon_scalesmearing* files follow this schema.
    // e.g. see `cb_params` blocks.  :contentReference[oaicite:3]{index=3}
    const double abseta = std::abs(mu._eta);
    // Prefer your real tracker-layer accessor if available.
    const double nL = getNTrackerLayersOrDefault(mu);

    try {
        auto cb = muonScaleSmearCorrections->at("cb_params");
        const double alpha = cb->evaluate({abseta, nL, 0});
        const double n     = cb->evaluate({abseta, nL, 1});
        const double sigma = cb->evaluate({abseta, nL, 2});
        const double tail  = cb->evaluate({abseta, nL, 3});

        // Minimal, safe smearing using the CB width as a proxy.
        // (If you want true CB-tail sampling, you can replace this with a Crystal Ball draw.)
        static TRandom3 rng(12345);
        const double gaus = rng.Gaus(0.0, std::abs(sigma)*1e-3); // scale down to a relative width
        const double smear = 1.0 + gaus;

        // Keep factors positive
        return smear > 0.0 ? (float)smear : 1.0f;
    } catch (...) {
        // If the set has only polynomial params (some eras), fall back to a light Gaussian
        // based on a quadratic form p0 + p1/pt + p2/pt^2 (if present).
        try {
            auto poly = muonScaleSmearCorrections->at("poly_params"); // present in some files
            const double p0 = poly->evaluate({abseta, nL, 0});
            const double p1 = poly->evaluate({abseta, nL, 1});
            const double p2 = poly->evaluate({abseta, nL, 2});
            const double pt = std::max(1e-3f, mu._pt);
            const double sigma2 = std::max(0.0, p0*p0 + p1*p1/pt + p2*p2/(pt*pt));
            const double sigma = std::sqrt(sigma2);
            static TRandom3 rng2(6789);
            const double smear = rng2.Gaus(1.0, sigma);
            return smear > 0.0 ? (float)smear : 1.0f;
        } catch (...) {
            return 1.0f;
        }
    }
}

// -------------------------------------
//             JEC/JER 
// -------------------------------------
    
// Build the vector<correction::Variable::Type> that correctionlib expects
template <typename... Args>
static inline std::vector<correction::Variable::Type> _corrVals(Args&&... args) {
std::vector<correction::Variable::Type> v;
    v.reserve(sizeof...(Args));
    (v.emplace_back(correction::Variable::Type(std::forward<Args>(args))), ...);
    return v;
}

// Safe evaluate with fallback
template <typename... Args>
static inline double tryEval(const std::string& key, double fallback, Args&&... args) {
    try {
    return jetCorrections->at(key)->evaluate(_corrVals(std::forward<Args>(args)...));
    } catch (...) {
        return fallback;
    }
}

// JET HORN MITIGATION HELPERS

static inline bool isHFHorn(float abseta) {
    return (abseta > 3.0f && abseta < 5.0f);
}

static inline bool isHEHorn(float abseta) {
    return (abseta > 2.5f && abseta < 3.0f);
}

static inline bool isTransitionHorn2024(float abseta) {
    return (abseta > 2.0f && abseta < 2.5f);
}

static inline bool eraIs2022or2023(const std::string& era) {
    return (era == "12022" || era == "22022" ||
            era == "12023" || era == "22023");
}

static inline bool eraNeedsHF50(const std::string& era) {
    // 2022, 2023 => HF requires pT > 50
    return eraIs2022or2023(era);
}

static inline bool eraNeedsHE50(const std::string& era) {
    // 2022, 2023, 2024 => HE requires pT > 50
    return (eraIs2022or2023(era) || era == "2024");
}

static inline float hornMitigatedJetPtCut(float eta, const std::string& era) {
    const float abseta = std::abs(eta);

    if (isHFHorn(abseta) && eraNeedsHF50(era)) return 50.f;
    if (isHEHorn(abseta) && eraNeedsHE50(era)) return 50.f;

    return 25.f; // your current default jet threshold
}

static inline float getDataResidualFactor(float eta, float pt, int run) {
    if (!jetCorrections) return 1.f;

    const JecEraConfig cfg = jecConfigForEra(g_era, __inputfilename);
    const std::string jecBase = cfg.prefixBase;

    float res = (float)tryEval(jecBase + "_DATA_L2L3Residual_AK4PFPuppi",
                             0.0, (double)run, (double)eta, (double)pt);

    if (res == 0.0f)
        res = (float)tryEval(jecBase + "_DATA_L2L3Residual_AK4PFPuppi",
                             0.0, (double)eta, (double)pt);

    if (res == 0.0f)
        res = (float)tryEval(jecBase + "_DATA_L2L3Res_AK4PFPuppi",
                             0.0, (double)run, (double)eta, (double)pt);

    if (res == 0.0f)
        res = (float)tryEval(jecBase + "_DATA_L2L3Res_AK4PFPuppi",
                             1.0, (double)eta, (double)pt);

    return (res > 0.f ? res : 1.f);
}

float getCorrectedJetPt(const Jet& jet, float rho, bool isMC, float genJetPt, bool hasGenJet, int run) {
  if (!jetCorrections) {
    static int n = 0;
    if (n < 3) {
      std::cerr << "[ERROR] Jet corrections not loaded!\n";
      if (++n == 3) std::cerr << "[ERROR] (Further jet correction warnings suppressed...)\n";
    }
    return jet._pt;
  }

  const JecEraConfig cfg = jecConfigForEra(g_era, __inputfilename);

  // Canonical inputs
  const float eta  = std::clamp(jet._eta, -5.191f, 5.191f);
  const float phi  = TVector2::Phi_mpi_pi(jet._phi);
  const float area = std::max(0.f, jet._area);
  const float pt0  = std::max(0.1f, jet._pt);

  const float abseta = std::abs(eta);
  const std::string jecBase = cfg.prefixBase;

  float jec = 1.f;

  if (!isMC) {
    // Prefer compound key
    const std::string compA = jecBase + "_DATA_L1L2L3Res_AK4PFPuppi";
    const std::string compB = jecBase + "_DATA_L1L2L3Residual_AK4PFPuppi";

    double comp = 0.0;
    comp = tryEval(compA, 0.0, (double)area, (double)eta, (double)pt0, (double)rho, (double)phi, (double)run);
    if (comp == 0.0) comp = tryEval(compA, 0.0, (double)area, (double)eta, (double)pt0, (double)rho, (double)run);
    if (comp == 0.0) comp = tryEval(compA, 0.0, (double)area, (double)eta, (double)pt0, (double)rho, (double)phi);
    if (comp == 0.0) comp = tryEval(compA, 0.0, (double)area, (double)eta, (double)pt0, (double)rho);
    if (comp == 0.0) comp = tryEval(compB, 0.0, (double)area, (double)eta, (double)pt0, (double)rho, (double)phi, (double)run);
    if (comp == 0.0) comp = tryEval(compB, 0.0, (double)area, (double)eta, (double)pt0, (double)rho, (double)run);
    if (comp == 0.0) comp = tryEval(compB, 0.0, (double)area, (double)eta, (double)pt0, (double)rho, (double)phi);
    if (comp == 0.0) comp = tryEval(compB, 0.0, (double)area, (double)eta, (double)pt0, (double)rho);

    if (comp > 0.0) {
      jec = (float)comp;
    } else {
      // Piecewise: L1 * L2 * L3 * Residual
      const std::string kL1 = jecBase + "_DATA_L1FastJet_AK4PFPuppi";
      const std::string kL2 = jecBase + "_DATA_L2Relative_AK4PFPuppi";
      const std::string kL3 = jecBase + "_DATA_L3Absolute_AK4PFPuppi";

      const float l1  = (float)tryEval(kL1, 1.0, (double)area, (double)eta, (double)pt0, (double)rho);
      const float l2a = (float)tryEval(kL2, 0.0, (double)eta, (double)phi, (double)pt0);
      const float l2b = (float)tryEval(kL2, 0.0, (double)eta, (double)pt0);
      const float l3a = (float)tryEval(kL3, 0.0, (double)eta, (double)phi, (double)pt0);
      const float l3b = (float)tryEval(kL3, 0.0, (double)eta, (double)pt0);
      const float L2  = (l2a > 0.f ? l2a : (l2b > 0.f ? l2b : 1.f));
      const float L3  = (l3a > 0.f ? l3a : (l3b > 0.f ? l3b : 1.f));

      float res = getDataResidualFactor(eta, pt0, run);
      jec = l1 * L2 * L3 * res;
    }
  } else {
    // MC: L1*L2*L3 (no residual)
    const std::string kL1 = jecBase + "_MC_L1FastJet_AK4PFPuppi";
    const std::string kL2 = jecBase + "_MC_L2Relative_AK4PFPuppi";
    const std::string kL3 = jecBase + "_MC_L3Absolute_AK4PFPuppi";

    const float l1  = (float)tryEval(kL1, 1.0, (double)area, (double)eta, (double)pt0, (double)rho);
    const float l2a = (float)tryEval(kL2, 0.0, (double)eta, (double)phi, (double)pt0);
    const float l2b = (float)tryEval(kL2, 0.0, (double)eta, (double)pt0);
    const float l3a = (float)tryEval(kL3, 0.0, (double)eta, (double)phi, (double)pt0);
    const float l3b = (float)tryEval(kL3, 0.0, (double)eta, (double)pt0);
    const float L2  = (l2a > 0.f ? l2a : (l2b > 0.f ? l2b : 1.f));
    const float L3  = (l3a > 0.f ? l3a : (l3b > 0.f ? l3b : 1.f));

    jec = l1 * L2 * L3;
  }

  float pt = std::max(0.f, pt0 * jec);

  // 2024 mitigation in 2.0 < |eta| < 2.5:
  // if corrected pT falls below 30 GeV, freeze the residual to its value at 30 GeV
  if (!isMC && g_era == "2024" && isTransitionHorn2024(abseta) && pt < 30.f) {
    const float res_at_pt = getDataResidualFactor(eta, std::max(pt, 1.f), run);
    const float res_at_30 = getDataResidualFactor(eta, 30.f, run);

    if (res_at_pt > 0.f && res_at_30 > 0.f) {
      pt *= (res_at_30 / res_at_pt);
    }
  }

  // ---- JER (MC only) ----
  if (isMC) {
    const double sf1   = tryEval(cfg.jerScaleKey, 1.0, (double)eta, (double)pt, std::string("nom"));
    const double sf2   = tryEval(cfg.jerScaleKey, 1.0, (double)eta, std::string("nom"));
    const double jerSF = (sf1 != 1.0 ? sf1 : sf2);

    const bool he_genmatched_only =
        isHEHorn(abseta) &&
        (g_era == "12022" || g_era == "22022" ||
         g_era == "12023" || g_era == "22023" ||
         g_era == "2024");

    if (hasGenJet && genJetPt > 0.f) {
      const double newPt = (double)genJetPt + jerSF * ((double)pt - (double)genJetPt);
      pt = (float)std::max(0.0, newPt);
    } else if (!he_genmatched_only && jerSF > 1.0) {
      const double r3     = tryEval(cfg.jerResKey, 0.0, (double)eta, (double)pt, (double)rho);
      const double r2     = tryEval(cfg.jerResKey, 0.0, (double)eta, (double)pt);
      const double jerRes = (r3 > 0.0 ? r3 : r2);

      if (jerRes > 0.0) {
        const double sigma = std::sqrt(std::max(0.0, jerSF*jerSF - 1.0)) * jerRes;
        pt = std::max(0.f, pt * (float)g_rng.Gaus(1.0, sigma));
      }
    }
  }

  return pt;
}

// -------------------------------------
//              B-TAGGING
// -------------------------------------

// ParticleNet WP threshold from JSON (tries 2023 naming and capitalization)
static inline float pnetWPThreshold(const std::string& wp) {
    if (!btagCorrections) return 1e9f;
    for (const std::string key : {"particleNet_wp_values", "ParticleNet_wp_values"}) {
        try { return (float)btagCorrections->at(key)->evaluate({wp}); }
        catch (...) {}
    }
    std::cerr << "[BTAG] Missing ParticleNet wp_values for wp='" << wp << "'\n";
    return 1e9f;
}

static inline bool isTaggedPNet(const Jet& j, const std::string& wp) {
    const float cut = pnetWPThreshold(wp);
    return j._btagPNetB >= cut;
}

static inline void logPNetWPsOnce() {
    static bool printedWPonce = false;
    if (!printedWPonce && btagCorrections) {
        std::cerr << "[GOLDEN][BTAG] ParticleNet WP: "
                  << "L=" << pnetWPThreshold("L")
                  << " M=" << pnetWPThreshold("M") << "\n";
        printedWPonce = true;
    }
}

// Per-jet WP SF for ParticleNet (MC only; returns 1 on DATA)
static inline double pnetSF_WP_jet(const Jet& j, const std::string& wp, const std::string& syst = "central")
{
    if (g_disableBtagSF) return 1.0;
    if (!__isMC || !btagCorrections) return 1.0;
    const int flav = (j._hadronFlavour == 5 ? 5 : (j._hadronFlavour == 4 ? 4 : 0));
    const double abseta = std::fabs((double)j._eta);
    const double pt     = std::max((double)j._pt, 0.1);

    // choose source based on flavor
    const std::string source = (flav == 5 || flav == 4) ? "comb" : "light";

    // Tables exist in your 2023 JSON: particleNet_comb, particleNet_light
    // signature: (systematic, working_point, flavor, abseta, pt)
    for (const std::string base : {"particleNet_", "ParticleNet_"}) {
        const std::string name = base + source;
        try {
            const auto& c = btagCorrections->at(name);
            return c->evaluate({syst, wp, flav, abseta, pt});
        } catch (...) {
            // try next capitalization
        }
    }
    return 1.0;
}

// -------------------------------------
//              JET VETO 
// -------------------------------------
    
inline bool passesJetVetoMapMinimalSelection(const Jet& j, float corr_pt) {
    if (corr_pt <= 15.f) return false;
    if ((j._jetId & 4) == 0) return false; 
    if ((j._chEmEF + j._neEmEF) >= 0.9f) return false;
    return true;
}

inline bool isVetoedJet(float eta, float phi) {
  if (g_disableJetVeto)
    return false;
  if (!jetVetoMaps) { 
    static int nWarns = 0;
    if (nWarns < 3) {
        std::cerr << "[ERROR] Jet veto corr not loaded!" << std::endl;
        if (++nWarns == 3)
            std::cerr << "[ERROR] (Further jet veto warnings suppressed...)\n";
    }
    return false; 
   }

  const JecEraConfig cfg = jecConfigForEra(g_era, __inputfilename);
  if (cfg.vetoMapKey.empty()) return false; // no maps listed for this era
  try {
    const double val = jetVetoMaps->at(cfg.vetoMapKey)->evaluate({std::string("jetvetomap"), (double)eta, (double)phi});
    return (val != 0.0);
  } catch (...) {
    static bool warned = false;
    if (!warned) {
      std::cerr << "[WARN] Veto map key '" << cfg.vetoMapKey << "' not found for era " << g_era << "\n";
      warned = true;
    }
    return false;
  }
}

// -------------------------------------
//           JET SEPERATION
// -------------------------------------

float jetMuondR(float jeta,float jphi, float meta, float mphi)
{
	TLorentzVector p4j,p4m;
	p4j.SetPtEtaPhiM(10, jeta, jphi, 0);
	p4m.SetPtEtaPhiM(10, meta, mphi, 0);
	return p4j.DeltaR(p4m);
}

//*******************************************************************
//                      CRITERIA CONDITIONS
//*******************************************************************

// -------------------------------------
//              MUON LOCATIONS
// -------------------------------------

bool isBarrel(const Muon& m)  { return std::abs(m._eta) < 0.9; }
bool isOverlap(const Muon& m) { return std::abs(m._eta) >= 0.9 && std::abs(m._eta) < 1.6; }
bool isEndcap(const Muon& m)  { return std::abs(m._eta) >= 1.6 && std::abs(m._eta) < 2.4; }

// -------------------------------------
//              FLAG CRITERIA 
// -------------------------------------

inline bool passFlagFilters(const Flags* fv, bool isData /* = !__isMC */) {
    if (g_disableFlagFilters) return true;
    if (!fv || fv->empty()) return true; // be permissive if not present
    const Flag& f = fv->front();
    return  //f._goodVertices
         f._globalSuperTightHalo2016Filter
         && f._EcalDeadCellTriggerPrimitiveFilter
         && f._BadPFMuonFilter
         && f._BadPFMuonDzFilter
         && f._hfNoisyHitsFilter
         && f._ecalBadCalibFilter
         && ( !isData || f._eeBadScFilter );
}

// -------------------------------------
//          VERTEX CRITERIA 
// -------------------------------------

inline bool passVertex(const Vertices* v) {
    if (!v || v->empty()) return false;
    return v->front()._npvsGood >= 1;
}

// -------------------------------------
//          MUON/LEPTONS CRITERIA 
// -------------------------------------

static inline bool useBeamConstraintPt(const Muon& mu) {

  if (g_disableBeamConstraint) return false;

  // mimic awkward mask: (bsConstrainedChi2 < 30)
  if (!std::isfinite(mu._bsConstrainedChi2)) return false;
  if (mu._bsConstrainedChi2 >= 30.f) return false;

  if (!std::isfinite(mu._bsConstrainedPt)) return false;
  if (mu._bsConstrainedPt <= 0.f) return false;

  return true;
}

static inline float baseMuonPtForCorrections(const Muon& mu) {
  return useBeamConstraintPt(mu) ? mu._bsConstrainedPt : mu._pt;
}

inline float muCorrScale(const Muon& mu, bool isMC, uint32_t evt = 0, uint32_t lumi = 0) {
    if (__noMuCorr) return 1.0f;

    float scale = 1.0f;

    if (__useKIT) {
        if (!kitCorrections) return 1.0f;

        const double pt0 = baseMuonPtForCorrections(mu);
        if (pt0 <= 0.0) return 1.0f;

        const bool is_data = !isMC;
        const int q = (mu._charge > 0 ? +1 : -1);

        const double eta = mu._eta;
        const double phi = mu._phi;

        // Scale (data + MC)
        const double pt_scaled = kit::pt_scale(*kitCorrections, is_data, pt0, eta, phi, q);

        // Smear (MC only), applied as a ratio wrt basePt
        double pt_final = pt_scaled;
        if (isMC) {
            const double pt_smeared_from_base = kit::pt_resol(*kitCorrections, mu, evt, lumi);
            const double smear_ratio = pt_smeared_from_base / pt0;
            pt_final = pt_scaled * smear_ratio;
        }

        return (float)(pt_final / pt0);
    }

    if (__useRochester) {
        scale = getMuonRochester(mu, isMC);
        if (!g_disableMuonScaleSmear) {
            scale *= getMuonScaleSmear(mu, isMC);
        }
        return scale;
    }

    return 1.0f;
}

inline float correctedBsPt(const Muon& mu, bool isMC, uint32_t evt = 0, uint32_t lumi = 0) {
  const float pt0 = baseMuonPtForCorrections(mu);
  return pt0 * muCorrScale(mu, isMC, evt, lumi);
}

bool vetoAdditionalLeptons(Muons* muons, Electrons* electrons, bool isMC)
{
    int nMuons20 = 0;

    // Count muons with pT>20, |eta|<2.4, global||tracker
    for (const auto& mu : *muons) {
        float ptCorr = baseMuonPtForCorrections(mu) * muCorrScale(mu, isMC);

        if (ptCorr <= 20.0f)       continue;
        if (std::abs(mu._eta) >= 2.4f) continue;
        if (!(mu._isGlobal || mu._isTracker)) continue;

        ++nMuons20;
    }

    // Require exactly 2 muons above 20 GeV
    if (nMuons20 != 2) {
        return false;
    }

    // Electron veto: pT>20, |eta|<2.5, ignoring 1.44<|eta|<1.57, passing WP90 MVA-based ID
    for (const auto& el : *electrons) {
        const float aeta = std::abs(el._eta);

        if (el._pt <= 20.0f) continue;
        if (aeta >= 2.5f)    continue;
        if (aeta > 1.44f && aeta < 1.57f) continue;
        if (!el._mvaIso_WP90) continue; 

        return false; // veto event
    }

    return true;
}

bool passMuonHLT(TrigObjs* trig, Muon const& m, bool isMC, uint32_t evt, uint32_t lumi)
{

static int debugCounter = 0;
const int maxDebug = 1;
bool debug = (debugCounter < maxDebug);
++debugCounter;

if (debug) {
std::cout << "\n[passMuonHLT] ===== Call " << debugCounter << " =====\n";
std::cout << "  isMC = " << isMC << "\n";
}

if (!trig) {
if (debug) std::cout << "  [ERROR] trig pointer is null -> return false\n";
return false;
}

if (debug) {
std::cout << "  Number of trigger objects: " << trig->size() << "\n";
}

if (trig->size() == 0) {
if (debug) std::cout << "  [WARN] No trigger objects in event -> return false\n";
return false;
}

    float corrected_pt = correctedBsPt(m, isMC, evt, lumi);
    
if (debug) {
std::cout << "  Muon: pt_corr = " << corrected_pt
  << ", eta = " << m._eta
  << ", phi = " << m._phi << "\n";
}

if (corrected_pt < 26.0f) {
if (debug) std::cout << "  Muon fails pT > 26 GeV requirement -> return false\n";
return false;
}

    int idx = 0;
    for (auto it = trig->begin(); it != trig->end(); ++it, ++idx) {

if (debug) {
std::cout << "    [TrigObj " << idx << "] "
<< "id = " << it->_id
<< ", pt = " << it->_pt
<< ", eta = " << it->_eta
<< ", phi = " << it->_phi
<< ", filterBits(dec) = " << it->_filterBits
<< ", filterBits(hex) = 0x" << std::hex << it->_filterBits << std::dec
<< "\n";
}

        // Only require that this is a muon trigger object.
        if (std::abs(it->_id) != 13) {
            if (debug) std::cout << "      -> skip: abs(id) != 13\n";
            continue;
        }

        // ΔR matching between offline muon and trigger object
        float deta = m._eta - it->_eta;
        float dphi = std::fabs(m._phi - it->_phi);
        if (dphi > M_PI) dphi = 2.0f * M_PI - dphi;
        float dR = std::sqrt(deta * deta + dphi * dphi);

        if (debug) {
            std::cout << "      deta = " << deta
                      << ", dphi = " << dphi
                      << ", dR = " << dR << "\n";
        }

        if (dR < 0.3f) {
            if (debug) {
                std::cout << "      ==> MATCH: dR < 0.3, muon PASSES HLT match\n";
            }
            return true;
        } else if (debug) {
            std::cout << "      -> no match: dR >= 0.3\n";
        }
    }

    if (debug) {
        std::cout << "  No trigger object matched this muon -> return false\n";
    }
    return false;
}

bool passMuon(Muon const& m, float corrected_pt, MuonType id=kLead)
{
    if (id==kLead)
    {
	    if ((m._isGlobal || m._isTracker) 
            && corrected_pt > _leadmuonPt 
            && TMath::Abs(m._eta) < _leadmuonEta 
            && m._mediumId
            && m._pfRelIso04_all <= _leadmuonIso)
		    return true;
    }
    else
    {
	    if ((m._isGlobal || m._isTracker) 
            && corrected_pt > _subleadmuonPt 
            && TMath::Abs(m._eta) < _subleadmuonEta 
            && m._mediumId
            && m._pfRelIso04_all <= _subleadmuonIso)
		    return true;
    }
	return false;
}

bool passMuons(Muon const& m1, Muon const& m2, TrigObjs* trigObjs, uint32_t evt, uint32_t lumi)
{
    float pt1_corr = correctedBsPt(m1, __isMC, evt, lumi);
    float pt2_corr = correctedBsPt(m2, __isMC, evt, lumi);

    const Muon& leadmuon    = (pt1_corr > pt2_corr) ? m1 : m2;
    const Muon& subleadmuon = (pt1_corr > pt2_corr) ? m2 : m1;
    float lead_pt_corr      = std::max(pt1_corr, pt2_corr);
    float sublead_pt_corr   = std::min(pt1_corr, pt2_corr);

    // OS + kinematic/ID/iso
    bool osOK   = (m1._charge != m2._charge);
    bool leadOK = passMuon(leadmuon,    lead_pt_corr,    kLead);
    bool subOK  = passMuon(subleadmuon, sublead_pt_corr, kSubLead);

    if (!osOK || !leadOK || !subOK)
        return false;

    // HLT: require at least one muon matched to trigger with pt>26 (this is inside passMuonHLT)
    bool hlt1OK = passMuonHLT(trigObjs, m1, __isMC, evt, lumi);
    bool hlt2OK = passMuonHLT(trigObjs, m2, __isMC, evt, lumi);
    bool hltOK  = (hlt1OK || hlt2OK);

    if (!hltOK)
        return false;

    return true;
}

// -------------------------------------
//          JET IDENTIFICATION
// -------------------------------------

/*enum class JetAlgo { AK4PUPPI, AK4CHS, AK8PUPPI };
enum class JetIDWP { Tight, TightLeptonVeto };

inline const char* jetIdCorrectionName(JetAlgo algo, JetIDWP wp) {
    switch (algo) {
        case JetAlgo::AK4PUPPI:
            return (wp == JetIDWP::TightLeptonVeto) ? "AK4PUPPI_TightLeptonVeto" : "AK4PUPPI_Tight";
        case JetAlgo::AK4CHS:
            return (wp == JetIDWP::TightLeptonVeto) ? "AK4CHS_TightLeptonVeto" : "AK4CHS_Tight";
        case JetAlgo::AK8PUPPI:
            return (wp == JetIDWP::TightLeptonVeto) ? "AK8PUPPI_TightLeptonVeto" : "AK8PUPPI_Tight";
    }
    return "AK4PUPPI_Tight";
}

bool passJetIDfromJSON(const Jet& jet,
                       JetAlgo algo = JetAlgo::AK4PUPPI,
                       JetIDWP wp   = JetIDWP::Tight)
{
    if (!jetIdCorrections) {
        std::cerr << "[WARN] JetID corrections not loaded, falling back to NanoAOD jetId bits.\n";
        // NanoAOD: 2=loose, 4=tight, 6=tightLepVeto. Keep your previous behavior:
        if (wp == JetIDWP::TightLeptonVeto) return ( (jet._jetId & 4) == 4 );
        return ( (jet._jetId & 2) == 2 ); // Tight
    }

    try {
        const auto& corr = jetIdCorrections->at(jetIdCorrectionName(algo, wp));

        // Inputs expected by the JSON:
        //   eta, chHEF, neHEF, chEmEF, neEmEF, muEF,
        //   chMultiplicity, neMultiplicity, multiplicity (= ch + ne)
        const int chMult = static_cast<int>(jet._chMult);
        const int neMult = static_cast<int>(jet._neMult);
        const int mult   = chMult + neMult;

        // Evaluate → 1.0 if passes, 0.0 otherwise
        const double result = corr->evaluate({
            (double)jet._eta,
            (double)jet._chHEF,
            (double)jet._neHEF,
            (double)jet._chEmEF,
            (double)jet._neEmEF,
            (double)jet._muEF,
            chMult,
            neMult,
            mult
        });

        return (result >= 0.5);
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] JetID evaluate failed: " << e.what()
                  << "  (eta=" << jet._eta << ")\n";
        // Be conservative on error:
        return false;
    }
}*/

//*******************************************************************
//                          SCALE FACTORS
//*******************************************************************

// -------------------------------------
//             MUONS 
// -------------------------------------

float getMuonSF(const Muon& mu, bool isMC, bool applyTriggerSF = false) {
    static int errorCount = 0;
    if (g_disableMuonSF) return 1.0f;
    if (!isMC) {
        if (errorCount < 5) {
            std::cerr << "This is Data!\n";
        }
        ++errorCount;
        return 1.0;
    }

    if (!corrections) {
        if (errorCount < 5) {
            std::cerr << "[ERROR] Correction set not initialized!\n";
        }
        ++errorCount;
        return 1.0;
    }

    float eta = std::abs(mu._eta);
    float pt = mu._pt;

    // Clamp to avoid going outside bin edges
    eta = std::min(eta, 2.39f);          // upper edge is 2.4 (exclusive)
    pt  = std::max(pt, 15.01f);          // lower edge is 15 (exclusive)

    float id_sf = 1.0;
    float iso_sf = 1.0;
    float trig_sf = 1.0;

    try {
        id_sf   = corrections->at(g_muSFKeyID)->evaluate({eta, pt, "nominal"});
        iso_sf  = corrections->at(g_muSFKeyISO)->evaluate({eta, pt, "nominal"});
        if (applyTriggerSF && pt >= 26) {
            trig_sf = corrections->at(g_muSFKeyTRG)->evaluate({eta, pt, "nominal"});
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Failed to evaluate muon SF: " << e.what() << std::endl;
        return 1.0;
    }

    //std::cout << "[DEBUG] Muon SF: eta = " << eta << ", pt = " << pt << ", SF = " << id_sf << std::endl;
    return id_sf * iso_sf * trig_sf;
}

// -------------------------------------
//            B-TAGGED JETS 
// -------------------------------------

static inline double clamp01(double x) {
    if (!std::isfinite(x)) return 0.0;
    return std::max(0.0, std::min(1.0, x));
}

// Replace this with your real efficiency map lookup.
// Must return the MC tagging efficiency for the requested WP
// in your analysis phase space, using the same selected jets.
static inline double getPNetEff_MC(const Jet& j, const std::string& wp)
{
    // TODO:
    // Look up epsilon(flavor, pt, abseta, wp) from your derived MC maps.
    // Example placeholders ONLY:
    const int flav = (j._hadronFlavour == 5 ? 5 : (j._hadronFlavour == 4 ? 4 : 0));
    const double pt = std::max(20.0, (double)j._pt);
    const double abseta = std::fabs((double)j._eta);

    // --- REMOVE THIS PLACEHOLDER ---
    // Just to keep code numerically safe before you plug real maps in.
    if (wp == "M") {
        if (flav == 5) return 0.70;
        if (flav == 4) return 0.20;
        return 0.01;
    }
    if (wp == "L") {
        if (flav == 5) return 0.85;
        if (flav == 4) return 0.35;
        return 0.10;
    }
    return 0.0;
}

static inline double pnetSF_WP_jet_safe(const Jet& j,
                                        const std::string& wp,
                                        const std::string& syst = "central")
{
    const double sf = pnetSF_WP_jet(j, wp, syst);
    if (!std::isfinite(sf) || sf <= 0.0) return 1.0;
    return sf;
}

// -------------------------------------
// Proper fixed-WP event reweighting for a veto using
// Method 1a extended to multiple WPs.
//
// Categories per selected jet:
//   - tagged at M
//   - tagged at L but not M
//   - not tagged at L
//
// Event veto logic is still decided by the reco tags in the event.
// The weight corrects the yield of those selected MC events.
// -------------------------------------
static inline double getBTagSF_PNet_WP_forVeto(const std::vector<const Jet*>& jets,
                                               const std::string& looseWP  = "L",
                                               const std::string& mediumWP = "M",
                                               const std::string& syst     = "central")
{
    if (g_disableBtagSF) return 1.0;
    if (!__isMC || !btagCorrections) return 1.0;

    double w = 1.0;

    for (const Jet* pj : jets) {
        const Jet& j = *pj;

        const bool passM = isTaggedPNet(j, mediumWP);
        const bool passL = isTaggedPNet(j, looseWP);

        const double epsM = clamp01(getPNetEff_MC(j, mediumWP));
        const double epsL = clamp01(getPNetEff_MC(j, looseWP));

        // sanity: by construction loose efficiency must be >= medium efficiency
        const double epsLonly = std::max(0.0, epsL - epsM);

        const double sfM = pnetSF_WP_jet_safe(j, mediumWP, syst);
        const double sfL = pnetSF_WP_jet_safe(j, looseWP,  syst);

        double factor = 1.0;

        if (passM) {
            // tagged at medium
            factor = sfM;
        }
        else if (passL) {
            // tagged at loose but not medium
            // BTV recommended multi-WP factor:
            // (SF_L * epsL - SF_M * epsM) / (epsL - epsM)
            const double denom = epsLonly;
            const double numer = sfL * epsL - sfM * epsM;

            // BTV warns this can become negative in corner cases.
            // Recommended protection: set to unity if numerator < 0
            // to avoid unphysical negative event weights.
            if (!std::isfinite(denom) || denom <= 0.0) {
                factor = 1.0;
            } else if (!std::isfinite(numer) || numer < 0.0) {
                factor = 1.0;
            } else {
                factor = numer / denom;
            }
        }
        else {
            // not tagged at loose
            // (1 - SF_L * epsL) / (1 - epsL)
            const double denom = 1.0 - epsL;
            const double numer = 1.0 - sfL * epsL;

            if (!std::isfinite(denom) || denom <= 0.0) {
                factor = 1.0;
            } else if (!std::isfinite(numer) || numer < 0.0) {
                factor = 1.0;
            } else {
                factor = numer / denom;
            }
        }

        if (!std::isfinite(factor) || factor <= 0.0) {
            factor = 1.0;
        }

        w *= factor;
    }

    if (!std::isfinite(w) || w <= 0.0) return 1.0;
    return w;
}

// wrapper for std::vector<Jet>
static inline double getBTagSF_PNet_WP_forVeto(const std::vector<Jet>& jets,
                                               const std::string& looseWP  = "L",
                                               const std::string& mediumWP = "M",
                                               const std::string& syst     = "central")
{
    std::vector<const Jet*> ptrs;
    ptrs.reserve(jets.size());
    for (const auto& j : jets) ptrs.push_back(&j);
    return getBTagSF_PNet_WP_forVeto(ptrs, looseWP, mediumWP, syst);
}

//*******************************************************************
//                         CATEGORIZE FUNCTION 
//*******************************************************************

void categorize(analysis::core::Jets* jets, analysis::core::GenJets* genjets, const analysis::core::Muon& mu1, 
const analysis::core::Muon& mu2, analysis::core::FsrPhotons* fsrphotons, analysis::core::Photons* photons, 
analysis::core::Vertices* vertices, const analysis::core::MET& met, const analysis::core::Event& event, 
analysis::core::Rhos* rhos, HistSetCollection& hists, float Dataweight, uint32_t evt, uint32_t lumi, 
float w_nom, float w_pu) {

    // -------------------------------------
    //             VERTICES 
    // -------------------------------------
    
 //   int npvsGood = 0;
  //  if (vertices && !vertices->empty()) {
    //    npvsGood = vertices->front()._npvsGood;
   // }

int totalpvsGood = -1;
if (vertices && !vertices->empty()) {
    totalpvsGood = vertices->front()._npvsGood;
}

/*static int dbg3 = 0;
if (dbg3 < 10) {
std::cout << "[STEP3 - CATEGORIZE ENTRY] totalpvsGood="
<< totalpvsGood
<< std::endl;
++dbg3;
}*/

    // -------------------------------------
    //             MUONS 
    // -------------------------------------
    
    // --- KIT, Rochester, noMuCorr Dimuon Set ---
    const float sf1 = muCorrScale(mu1, __isMC, evt, lumi);
    const float sf2 = muCorrScale(mu2, __isMC, evt, lumi); 

    const float pt1_base = (mu1._bsConstrainedPt > 0.f ? mu1._bsConstrainedPt : mu1._pt);
    const float pt2_base = (mu2._bsConstrainedPt > 0.f ? mu2._bsConstrainedPt : mu2._pt);

    // Raw (base) TLVs
    TLorentzVector p4m1, p4m2;
    p4m1.SetPtEtaPhiM(pt1_base, mu1._eta, mu1._phi, PDG_MASS_Mu);
    p4m2.SetPtEtaPhiM(pt2_base, mu2._eta, mu2._phi, PDG_MASS_Mu);

    // Corrected TLVs (consistent with KIT base)
    TLorentzVector corr_p4m1, corr_p4m2;
    corr_p4m1.SetPtEtaPhiM(correctedBsPt(mu1, __isMC, evt, lumi), mu1._eta, mu1._phi, PDG_MASS_Mu);
    corr_p4m2.SetPtEtaPhiM(correctedBsPt(mu2, __isMC, evt, lumi), mu2._eta, mu2._phi, PDG_MASS_Mu);

    TLorentzVector corr_prefsr_p4 = corr_p4m1 + corr_p4m2;
    
    // -------------------------------------
    //             FSR PHOTONS 
    // -------------------------------------
    
    // Add FSR photons to both muons if any pass criteria
    auto applyFSR = []
    (const analysis::core::Muon& mu, TLorentzVector& p4, TLorentzVector& corr_p4, const analysis::core::FsrPhotons* fsrphotons, int muonIdx) {
    const FsrPhoton* bestFsrPhoton = nullptr;
    float min_dROverEt2 = std::numeric_limits<float>::max();

    for (const auto& fsr : *fsrphotons) {
        
        // Basic photon cuts
        if (fsr._pt <= 2.0) continue;
        if (std::abs(fsr._eta) >= 2.5) continue;
        
        // Associated muon/electron index
        if (fsr._muonIdx != muonIdx) continue;
        if (fsr._electronIdx >= 0) continue;
        
        // Isolation and collinearity cuts
        if (fsr._relIso03 / fsr._pt >= 1.8) continue;
        if (fsr._dROverEt2 >= 0.012) continue;
        
        // pT(γ) / pT(μ) cut
        if ((fsr._pt / mu._pt) >= 0.4) continue;

        // Passes all selection — now track best candidate
        if (fsr._dROverEt2 < min_dROverEt2) {
            min_dROverEt2 = fsr._dROverEt2;
            bestFsrPhoton = &fsr;
        }
    }

        if (bestFsrPhoton) {
            TLorentzVector p4fsr;
            p4fsr.SetPtEtaPhiM(bestFsrPhoton->_pt, bestFsrPhoton->_eta, bestFsrPhoton->_phi, 0.0);
            p4 += p4fsr;
            corr_p4 += p4fsr;
            }
    };

    // Apply to mu1 (index 0) and mu2 (index 1)
    if (!g_disableFSRRecovery && fsrphotons && !fsrphotons->empty()) {
    applyFSR(mu1, p4m1, corr_p4m1, fsrphotons, 0);
    applyFSR(mu2, p4m2, corr_p4m2, fsrphotons, 1);
    }

    // With applied FSR corrections 
    TLorentzVector p4dimuon = p4m1 + p4m2;    
    TLorentzVector corr_p4dimuon = corr_p4m1 + corr_p4m2;
    
    // -------------------------------------
    //         COLLINS-SOPER ANGLE 
    // -------------------------------------
    
    double dphi = p4m1.DeltaPhi(p4m2);

    // Collins-Soper Angles
    const bool mu1_is_minus = (mu1._charge < 0);
    const TLorentzVector& muMinus = mu1_is_minus ? corr_p4m1 : corr_p4m2;
    const TLorentzVector& muPlus  = mu1_is_minus ? corr_p4m2 : corr_p4m1;
	
    // ---> Boost everything to the dimuon rest frame
    const TVector3 boostToQ = -corr_p4dimuon.BoostVector();
    TLorentzVector lminus = muMinus; lminus.Boost(boostToQ);
    TLorentzVector lplus  = muPlus;  lplus .Boost(boostToQ);

    // ---> Massless beam four-vectors along ±z (only directions matter)
    TLorentzVector pBeam1(0., 0., +1., 1.);  // +z
    TLorentzVector pBeam2(0., 0., -1., 1.);  // -z
    pBeam1.Boost(boostToQ);
    pBeam2.Boost(boostToQ);

    // Collins–Soper axes in the dilepton rest frame
    const TVector3 b1 = pBeam1.Vect().Unit();
    const TVector3 b2 = pBeam2.Vect().Unit();
    const TVector3 zcs = (b1 - b2).Unit();
    const TVector3 xcs = (b1 + b2).Unit();
    const TVector3 ycs = zcs.Cross(xcs).Unit();

    // Negative muon direction in that frame
    const TVector3 lhat = lminus.Vect().Unit();

    // Angles
    const double cosThetaCS = lhat.Dot(zcs);
    const double phiCS      = std::atan2(lhat.Dot(ycs), lhat.Dot(xcs)); // [-π, π]

    // -------------------------------------
    //              JETS
    // -------------------------------------

    bool pass_bjet_veto = true;
    int n_loose_bjets   = 0;
    int n_medium_bjets  = 0;
    int n_tight_bjets   = 0;
    const float rho     = rhos->at(0)._fixedGridRhoFastjetAll;

    // Cache ParticleNet WPs once per event
    const float pnetL = pnetWPThreshold("L");
    const float pnetM = pnetWPThreshold("M");
    const float pnetT = pnetWPThreshold("T");

    std::vector<TLorentzVector> p4jets;
    std::vector<analysis::core::Jet> btagJets;

    p4jets.reserve(jets->size());
    btagJets.reserve(jets->size());

/*std::cerr
<< "[BTAG WP DBG]"
<< " era=" << g_era
<< " L=" << pnetL
<< " M=" << pnetM
<< " T=" << pnetT
<< "\n";*/

    // -------------------------------------
    //   Event-level jet veto map decision
    // -------------------------------------
    bool eventFailsJetVetoMap = false;

    if (!g_disableJetVeto) {
        for (const auto& j : *jets) {
            float genJetPt_veto = 0.f;
            bool hasGenJet_veto = false;

            if (__isMC && j._genJetIdx >= 0 && j._genJetIdx < (int)genjets->size()) {
                const auto& gj = genjets->at(j._genJetIdx);
                genJetPt_veto = gj._pt;
                hasGenJet_veto = true;
            }

            const int runNum_veto = __isMC ? 0 : static_cast<int>(event._run);
            const float corr_pt_veto = getCorrectedJetPt(
                j, rho, __isMC, genJetPt_veto, hasGenJet_veto, runNum_veto
            );

            if (!passesJetVetoMapMinimalSelection(j, corr_pt_veto))
                continue;

            if (isVetoedJet(j._eta, j._phi)) {
                eventFailsJetVetoMap = true;
                break;
            }
        }
    }

    if (eventFailsJetVetoMap) return;

    // -------------------------------------
    //       Main selected-jet loop
    // -------------------------------------
    for (const auto& j : *jets) {
        ++jets_in;

        // Gen matching for JER
        float genJetPt = 0.f;
        bool hasGenJet = false;

        if (__isMC && j._genJetIdx >= 0 && j._genJetIdx < (int)genjets->size()) {
            const auto& gj = genjets->at(j._genJetIdx);
            genJetPt = gj._pt;
            hasGenJet = true;
        }

        const int runNum = __isMC ? 0 : static_cast<int>(event._run);
        const float corr_pt = getCorrectedJetPt(
            j, rho, __isMC, genJetPt, hasGenJet, runNum
        );

        const float raw_pt_safe = std::max(1e-6f, j._pt);
        const float corr_mass   = j._mass * (corr_pt / raw_pt_safe);

        // This counter means the event survived the earlier event-level veto stage
        ++jets_after_vetoMap;

        if (!g_disableJetCriteria) {
            if (j._jetId < 4)
                continue;
            ++jets_after_id;

            const float abseta = std::abs(j._eta);

            // Default pt cut
            if (corr_pt <= 25.f)
                continue;
            ++jets_after_pt_default;

            // HF horn mitigation
            if (isHFHorn(abseta) && eraNeedsHF50(g_era)) {
                if (corr_pt <= 50.f)
                    continue;
                ++jets_after_pt_hf;
            }

            // HE horn mitigation
            if (isHEHorn(abseta) && eraNeedsHE50(g_era)) {
                if (corr_pt <= 50.f)
                    continue;
                ++jets_after_pt_he;
            }

            ++jets_after_pt;

            if (abseta >= 4.7f)
                continue;
            ++jets_after_eta;

            if (jetMuondR(j._eta, j._phi, mu1._eta, mu1._phi) <= 0.4)
                continue;
            ++jets_after_dR1;

            if (jetMuondR(j._eta, j._phi, mu2._eta, mu2._phi) <= 0.4)
                continue;
            ++jets_after_dR2;
        }

        ++jets_pass;

        // Save corrected four-vector for later dijet kinematics
        TLorentzVector p4;
        p4.SetPtEtaPhiM(corr_pt, j._eta, j._phi, corr_mass);
        p4jets.push_back(p4);

        // Save corrected jet copy for b-tag efficiency / veto logic
        analysis::core::Jet jSF = j;
        jSF._pt            = corr_pt;
        jSF._mass          = corr_mass;
        jSF._eta           = j._eta;
        jSF._phi           = j._phi;
        jSF._btagPNetB     = j._btagPNetB;
        jSF._hadronFlavour = j._hadronFlavour;

        btagJets.push_back(jSF);
    }

    // -------------------------------------
    //   Fill b-tag efficiency maps + count tags
    // -------------------------------------
    for (const auto& jb : btagJets) {
        FILL_BTAG_EFF_HIST(hists.BTagEff_PreSelection, jb, pnetL, pnetM, Dataweight);

        const float disc = jb._btagPNetB;
        if (disc >= pnetL) ++n_loose_bjets;
        if (disc >= pnetM) ++n_medium_bjets;
        if (disc >= pnetT) ++n_tight_bjets;
    }

    pass_bjet_veto = !(n_medium_bjets >= 1 || n_loose_bjets >= 2);
    if (!pass_bjet_veto) return;

    // -------------------------------------
    //      Sort jets / derive jet vars
    // -------------------------------------
    std::sort(
        p4jets.begin(),
        p4jets.end(),
        [](const TLorentzVector& a, const TLorentzVector& b) {
            return a.Pt() > b.Pt();
        }
    );

    const bool is0jet = (p4jets.size() == 0);
    const bool is1jet = (p4jets.size() == 1);
    const bool is2jet = (p4jets.size() >= 2);

    const bool has1jet = (p4jets.size() >= 1);
    const bool has2jet = (p4jets.size() >= 2);

    float jet1eta   = -99.f, jet2eta   = -99.f;
    float jet1pt    = -1.f,  jet2pt    = -1.f;
    float dijetpt   = -1.f,  dijetdeta = -1.f;
    float dijetmass = -1.f,  dijetdphi = -1.f;

    TLorentzVector p4lead, p4sub, dijet;

    if (has1jet) {
        p4lead  = p4jets[0];
        jet1eta = p4lead.Eta();
        jet1pt  = p4lead.Pt();
    }

    if (has2jet) {
        p4sub   = p4jets[1];
        jet2eta = p4sub.Eta();
        jet2pt  = p4sub.Pt();

        dijet     = p4lead + p4sub;
        dijetpt   = dijet.Pt();
        dijetdeta = TMath::Abs(p4lead.Eta() - p4sub.Eta());
        dijetmass = dijet.M();
        dijetdphi = std::fabs(p4lead.DeltaPhi(p4sub));
    }

    float PTdDM1 = corr_p4m1.Pt() / corr_p4dimuon.M();
    float PTdDM2 = corr_p4m2.Pt() / corr_p4dimuon.M();

    // Choosing Largest Dijet Mass Pair
    /*double best_mjj = -1.0;
    int i_best = -1, j_best = -1;

    for (size_t i = 0; i < p4jets.size(); ++i) {
        for (size_t j = i + 1; j < p4jets.size(); ++j) {
            double mjj = (p4jets[i] + p4jets[j]).M();
            if (mjj > best_mjj) {
                best_mjj = mjj;
                i_best = static_cast<int>(i);
                j_best = static_cast<int>(j);
            }
        }
    }

    if (i_best < 0 || j_best < 0) return;
    TLorentzVector jA = p4jets[i_best];
    TLorentzVector jB = p4jets[j_best];
    TLorentzVector p4lead = (jA.Pt() >= jB.Pt()) ? jA : jB;
    TLorentzVector p4sub  = (jA.Pt()  <  jB.Pt()) ? jA : jB;*/
    
    // -------------------------------------
    //          SF && WEIGHTS  
    // -------------------------------------
   
    float mu1_SF            = getMuonSF(mu1, __isMC);                                          // Leading Muon Scale-Factors
    float mu2_SF            = getMuonSF(mu2, __isMC);                                          // Sub-Leading Muon Scale-Factors
    float btagjet_SF        = getBTagSF_PNet_WP_forVeto(btagJets, "L", "M", "central");        // B-Tagged Jets Scale-Factors 
    float totalWeight       = Dataweight * mu1_SF * mu2_SF * btagjet_SF;                       // Weight (MC)
    float noWeight          = 1.0;                                                             // No Weight (DATA)

    float w_noPU            = w_nom;                                                           // NO pileup Weight
    float w_withPU          = w_pu;                                                            // WITH pileup Weight
    float totalWeight_noPU  = w_noPU * mu1_SF * mu2_SF * btagjet_SF;                           // NO pileup Weight * 
    float totalWeight_PU    = w_withPU * mu1_SF * mu2_SF * btagjet_SF;

    
    // -------------------------------------
    //          Z-PT REWEIGHTING  
    // -------------------------------------
   
    if (corr_p4dimuon.M() > 70 && corr_p4dimuon.M() < 110) {
        if (is0jet) {
            FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_Z0j, totalWeight); 
        }
        else if (is1jet) {
            FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_Z1j, totalWeight); 
        }
        else if (is2jet) {
            FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_Z2j, totalWeight); 
        }
    }
    
    // Only apply to DY MC
    if (!g_disableZPtRW && __isMC && g_sampleName == "DY") {
      const double zpt = corr_p4dimuon.Pt();

      int nj = 2;
      if (is0jet) nj = 0;
      else if (is1jet) nj = 1;
      else nj = 2;

      const double w_zpt = g_zptRW.eval(zpt, nj);
      totalWeight *= w_zpt;
    }
 
    // -------------------------------------
    //          B-TAG REWEIGHTING  
    // -------------------------------------
   
static int DbgEvt = 0;
const int maxDbgEvts = 10;

if (DbgEvt++ < maxDbgEvts) {
std::cerr
<< "[WEIGHT DBG] run=" << event._run
<< " lumi=" << event._lumi
<< " evt=" << event._event
<< " | btagJets.size()=" << btagJets.size()
<< " | mu1_SF=" << mu1_SF
<< " | mu2_SF=" << mu2_SF
<< " | btagjet_SF=" << btagjet_SF
<< " | Dataweight=" << Dataweight
<< " | totalWeight=" << totalWeight
<< " | totalWeight_noPU=" << totalWeight_noPU
<< " | totalWeight_PU=" << totalWeight_PU
<< "\n";
}

    // -------------------------------------
    //        PRE-SELECTION HISTOGRAMS  
    // -------------------------------------
   
    FILL_MUON_HIST(hists.setDM_PreSelection, noWeight);   
    FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection, totalWeight); 
    FILL_JET_HIST(hists.setJet_PreSelection, noWeight); 
    FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection, totalWeight); 

    if (corr_p4dimuon.M() > 70 && corr_p4dimuon.M() < 110) {
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_Z, totalWeight); 
        if (is0jet) {
            FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_wZpt_Z0j, totalWeight); 
        }
        else if (is1jet) {
            FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_wZpt_Z1j, totalWeight); 
        }
        else if (is2jet) {
            FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_wZpt_Z2j, totalWeight); 
        }
    }

    if (corr_p4dimuon.M() > 110 && corr_p4dimuon.M() < 150) {
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_SR, totalWeight); 
        if (!__isMC) {
            FILL_VTX_HIST(hists.PileupReweighting_PreSelection, totalpvsGood, 1.0, 1.0);
        } else {
            FILL_VTX_HIST(hists.PileupReweighting_PreSelection, totalpvsGood, totalWeight_noPU, totalWeight_PU);
        }
    }

    if (isBarrel(mu1) && isBarrel(mu2))
    {
        //	BB
        FILL_MUON_HIST(hists.setDM_PreSelection_BB, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_BB, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_BB, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection_BB, totalWeight); 
    }
    else if ((isBarrel(mu1) && isOverlap(mu2)) ||
        (isBarrel(mu2) && isOverlap(mu1)))
    {
        //	BO
        FILL_MUON_HIST(hists.setDM_PreSelection_BO, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_BO, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_BO, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection_BO, totalWeight); 
    }
    else if ((isBarrel(mu1) && isEndcap(mu2)) || 
        (isBarrel(mu2) && isEndcap(mu1)))
    {
        //	BE
        FILL_MUON_HIST(hists.setDM_PreSelection_BE, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_BE, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_BE, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection_BE, totalWeight); 
    }
    else if (isOverlap(mu1) && isOverlap(mu2))
    {
        //	OO
        FILL_MUON_HIST(hists.setDM_PreSelection_OO, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_OO, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_OO, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection_OO, totalWeight); 
    }
    else if ((isOverlap(mu1) && isEndcap(mu2)) ||
        (isOverlap(mu2) && isEndcap(mu1)))
    {
        //	OE
        FILL_MUON_HIST(hists.setDM_PreSelection_OE, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_OE, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_OE, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection_OE, totalWeight); 
    }
    else if (isEndcap(mu1) && isEndcap(mu2))
    {
        //	EE
        FILL_MUON_HIST(hists.setDM_PreSelection_EE, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_EE, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_EE, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection_EE, totalWeight); 
    }
    
    // -------------------------------------
    //          VBF HISTOGRAMS  
    // -------------------------------------
    
    if (p4lead.Pt() < _leadJetPt) return;                       // Leading Jet Cut 
    if (p4sub.Pt() < _subleadJetPt) return;                     // Subleading Jet Cut
    if (n_loose_bjets >= 2 || n_medium_bjets >= 1) return;      // B-Tag Jet Cut 
    
        //  PRE-DIJET CUTS 
        FILL_JET_HIST(hists.setJet_PreDiJet_VBF, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreDiJet_VBF, totalWeight); 

        if (dijetmass > _dijetMass_VBF && TMath::Abs(dijetdeta) > _dijetdEta_VBF) {
            
            // -------------------------------------
            //         MVA INPUT VARIABLES  
            // -------------------------------------
            
            const double m_mumu = corr_p4dimuon.M();

            if(m_mumu > 110 && m_mumu < 150) {

                // Muon Variables
                g_m_mumu             = static_cast<float>(m_mumu);
                g_mu1_eta            = static_cast<float>(corr_p4m1.Eta());
                g_mu2_eta            = static_cast<float>(corr_p4m2.Eta());
                g_mu1_pt_over_m      = (m_mumu > 0.0) ? static_cast<float>(corr_p4m1.Pt() / m_mumu) : 0.f;
                g_mu2_pt_over_m      = (m_mumu > 0.0) ? static_cast<float>(corr_p4m2.Pt() / m_mumu) : 0.f;
                g_y_mumu             = static_cast<float>(corr_p4dimuon.Rapidity());
                g_costheta_CS        = static_cast<float>(cosThetaCS);
                g_phi_CS             = static_cast<float>(TVector2::Phi_mpi_pi(phiCS));

                // Jet Variables 
                g_jet1_eta           = static_cast<float>(p4lead.Eta());
                g_deta_jj            = std::fabs(p4lead.Eta() - p4sub.Eta());
                g_dphi_jj            = std::fabs(TVector2::Phi_mpi_pi(p4lead.Phi() - p4sub.Phi()));
                g_y_j1               = static_cast<float>(p4lead.Rapidity());
                g_y_j2               = static_cast<float>(p4sub.Rapidity());
                g_njets              = static_cast<int>(p4jets.size());
                g_dphi_mumu_jj       = std::fabs(TVector2::Phi_mpi_pi(corr_p4dimuon.Phi() - dijet.Phi()));
                g_deta_mumu_jj       = std::fabs(corr_p4dimuon.Eta() - dijet.Eta());

                // Zeppenfeld variable
                g_z             = std::fabs(g_y_mumu - 0.5f*(g_y_j1 + g_y_j2));

                // Pt-balance ratio
                const double px_sum     =  corr_p4dimuon.Px() + p4lead.Px() + p4sub.Px();
                const double py_sum     =  corr_p4dimuon.Py() + p4lead.Py() + p4sub.Py();
                const double vec_sum    = std::sqrt(px_sum*px_sum + py_sum*py_sum);

                g_pt_mumu       = static_cast<float>(corr_p4dimuon.Pt());
                g_pt_j1         = static_cast<float>(p4lead.Pt());
                g_pt_j2         = static_cast<float>(p4sub.Pt());

                const double scalar_sum = static_cast<double>(g_pt_mumu + g_pt_j1 + g_pt_j2);
                g_ptBalance     = (scalar_sum > 0.0) ? static_cast<float>(vec_sum / scalar_sum) : 0.f;

                // Weight (total per-event weight)
                g_weight = totalWeight;

                // Fill exactly once per selected dimuon+2jet pair
                if (gMvaTree) gMvaTree->Fill();

/*static int vtxDbg = 0;
if (vtxDbg < 20) {
std::cout << "[VTX DBG] totalpvsGood=" << totalpvsGood
<< " totalWeight_noPU=" << totalWeight_noPU
<< " totalWeight_PU=" << totalWeight_PU
<< " genWeight-sign test done"
<< std::endl;
++vtxDbg;
}*/

                //Filling Vertex Histograms 
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_SR, totalWeight)
                if (!__isMC) {
                    FILL_VTX_HIST(hists.PileupReweighting_VBF, totalpvsGood, 1.0, 1.0);
                } else {
                    FILL_VTX_HIST(hists.PileupReweighting_VBF, totalpvsGood, totalWeight_noPU, totalWeight_PU);
                }

            }
            
            if (corr_p4dimuon.M() > 70 && corr_p4dimuon.M() < 110) {
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_Z, totalWeight); 
            }

            //	VBF SELECTION CRITERIA
            FILL_MUON_HIST(hists.setDM_VBF, noWeight)
            FILL_CORR_MUON_HIST(hists.corr_setDM_VBF, totalWeight)
            FILL_JET_HIST(hists.setJet_VBF, noWeight); 
            FILL_CORR_JET_HIST(hists.corr_setJet_VBF, totalWeight); 
            //FILL_BTAG_EFF_HIST(hists.BTagEff_VBF, jSF, pnetL, pnetM, totalWeight);
            
            if (isBarrel(mu1) && isBarrel(mu2))
            {
                //	BB
                FILL_MUON_HIST(hists.setDM_VBF_BB, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_BB, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_BB, noWeight); 
                FILL_CORR_JET_HIST(hists.corr_setJet_VBF_BB, totalWeight); 
            }
            else if ((isBarrel(mu1) && isOverlap(mu2)) ||
                (isBarrel(mu2) && isOverlap(mu1)))
            {
                //	BO
                FILL_MUON_HIST(hists.setDM_VBF_BO, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_BO, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_BO, noWeight); 
                FILL_CORR_JET_HIST(hists.corr_setJet_VBF_BO, totalWeight); 
            }
            else if ((isBarrel(mu1) && isEndcap(mu2)) || 
                (isBarrel(mu2) && isEndcap(mu1)))
            {
                //	BE
                FILL_MUON_HIST(hists.setDM_VBF_BE, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_BE, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_BE, noWeight); 
                FILL_CORR_JET_HIST(hists.corr_setJet_VBF_BE, totalWeight); 
            }
            else if (isOverlap(mu1) && isOverlap(mu2))
            {
                //	OO
                FILL_MUON_HIST(hists.setDM_VBF_OO, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_OO, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_OO, noWeight); 
                FILL_CORR_JET_HIST(hists.corr_setJet_VBF_OO, totalWeight); 
            }
            else if ((isOverlap(mu1) && isEndcap(mu2)) ||
                (isOverlap(mu2) && isEndcap(mu1)))
            {
                //	OE
                FILL_MUON_HIST(hists.setDM_VBF_OE, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_OE, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_OE, noWeight); 
                FILL_CORR_JET_HIST(hists.corr_setJet_VBF_OE, totalWeight); 
            }
            else if (isEndcap(mu1) && isEndcap(mu2))
            {
                //	EE
                FILL_MUON_HIST(hists.setDM_VBF_EE, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_EE, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_EE, noWeight); 
                FILL_CORR_JET_HIST(hists.corr_setJet_VBF_EE, totalWeight); 
            }
		}
}

//*******************************************************************
//                        EVENT COUNTER??? 
//*******************************************************************

float sampleinfo(std::string const& inputname)
{
	Streamer s(inputname, NTUPLEMAKER_NAME+"/Meta");
	s.chainup();

	using namespace analysis::dimuon;
	MetaHiggs *meta=NULL;
	s._chain->SetBranchAddress("Meta", &meta);

	long long int numEvents = 0;
	long long int numEventsWeighted = 0;
	for (int i=0; i<s._chain->GetEntries(); i++)
	{
		s._chain->GetEntry(i);
		numEvents+=meta->_nEventsProcessed;
		numEventsWeighted+=meta->_sumEventWeights;
	}
	std::cout 
		<< "Number of total events processed = " << numEvents << std::endl
		<< "Number of tota levents weighted = " << numEventsWeighted << std::endl;

	return numEventsWeighted;
}


//*******************************************************************
//                       PROGRESS BAR 
//*******************************************************************

void showProgressBar(uint32_t current, uint32_t total) {
    int barWidth = 70; // Width of the progress bar
    float progress = static_cast<float>(current) / total;
    int pos = static_cast<int>(barWidth * progress);

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos)
            std::cout << "=";
        else
            std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(2) << progress * 100 << "%\r";
    std::cout.flush();  // Ensure output is immediately shown
}


//*******************************************************************
//                       READING WEIGHTS 
//*******************************************************************

// Parse an INI with section headers like [12022], [22022], etc.
// Keys inside a section:
//   luminosity_*         = <double>   (any key that starts with "luminosity")
//   <SAMPLE>.xsec        = <double>
//   <SAMPLE>.sumw        = <double>
static void loadWeightsFile(const std::string& path) {
    if (path == "None" || path == "DATA" || path.empty()) return;

    std::ifstream in(path);
    if (!in) {
        std::cerr << "[weights] Could not open file: " << path << "\n";
        return;
    }

    // simple trim lambda
    auto trim = [](std::string& s) {
        size_t b = s.find_first_not_of(" \t\r\n");
        size_t e = s.find_last_not_of(" \t\r\n");
        if (b == std::string::npos) { s.clear(); return; }
        s.assign(s, b, e - b + 1);
    };

    std::string line;
    std::string currentSection = "global";
    std::regex sec_re(R"(^\s*\[([^\]]+)\]\s*$)");
    std::regex kv_re(R"(^\s*([^#][^=]*)\s*=\s*([^\s#]+)\s*$)");

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        // Section header like: [22022]
        {
            std::smatch sm;
            if (std::regex_match(line, sm, sec_re)) {
                currentSection = sm[1].str();
                trim(currentSection);
                continue;
            }
        }

        // key=value
        std::smatch m;
        if (!std::regex_match(line, m, kv_re)) continue;

        std::string key = m[1];
        std::string val = m[2];
        trim(key);
        trim(val);

        try {
            // Any key starting with "luminosity" goes to this section's lumi map
            if (key.rfind("luminosity", 0) == 0) {
                double L = std::stod(val);
                g_sectionLumis[currentSection]["lumi"] = L; // canonical entry
                g_sectionLumis[currentSection][key]    = L; // keep original key
                continue;
            }

            // Expect <SAMPLE>.<field>
            auto dot = key.find('.');
            if (dot == std::string::npos) continue;
            std::string name  = key.substr(0, dot);
            std::string field = key.substr(dot + 1);
            trim(name);
            trim(field);

            double dval = std::stod(val);
            auto& s = g_sectionSamples[currentSection][name]; // per-section sample
            if (field == "xsec")          s.xsec      = dval;
            else if (field == "sumw")     s.sumw      = dval;
            else {
                std::cerr << "[weights] Unknown field '" << field
                          << "' for sample '" << name
                          << "' in section [" << currentSection << "]\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "[weights] Bad line: " << line << " (" << e.what() << ")\n";
        }
    }
}

// Returns the *per-event* base MC weight independent of per-event pu/genWeight.
// i.e. (lumi * xsec) / (sumw)
static double sampleBaseWeight(const std::string& sample) {
    // Look up the sample in the chosen era/section
    auto secIt = g_sectionSamples.find(g_era);
    if (secIt == g_sectionSamples.end()) {
        std::cerr << "[weights] Unknown era '" << g_era << "'. Using 1.0.\n";
        return 1.0;
    }
    const auto& samples = secIt->second;
    auto it = samples.find(sample);
    if (it == samples.end()) {
        std::cerr << "[weights] Unknown sample '" << sample
                  << "' in era '" << g_era << "'. Using 1.0.\n";
        return 1.0;
    }

    const auto& s = it->second;
    const double denom = s.sumw;
    if (denom <= 0.0) {
        std::cerr << "[weights] Non-positive denom for '" << sample
                  << "' in era '" << g_era << "'. Using 1.0.\n";
        return 1.0;
    }

    return (g_selected_lumi * s.xsec) / denom;
}


//*******************************************************************
//                        MAIN PROCESS FUCNTION 
//*******************************************************************

void process()
{
    //	Output Root File
    TFile *outroot = new TFile(__outputfilename.c_str(), "recreate");

    // Also record the exact JSON filenames used (if any)
    TNamed rj("muCorrRochJSON", (!__useRochester || g_rochJSON.empty()) ? "none" : g_rochJSON.c_str());
    rj.Write();
    TNamed kj("muCorrKITJSON", g_kitJSON.empty()  ? "none" : g_kitJSON.c_str());
    kj.Write();

    // Pileup Block
    VtxSet      PileupReweighting_PreSelection(outroot, "PileupReweighting_Pre-Selection");
    VtxSet      PileupReweighting_VBF(outroot, "PileupReweighting_VBF");

    // Z-Reweight Block
    MuonSet     corr_setDM_PreSelection_Z0j(outroot, "corr_DM_Pre-Selection_Z0j");
    MuonSet     corr_setDM_PreSelection_Z1j(outroot, "corr_DM_Pre-Selection_Z1j");
    MuonSet     corr_setDM_PreSelection_Z2j(outroot, "corr_DM_Pre-Selection_Z2j");
    MuonSet     corr_setDM_PreSelection_wZpt_Z0j(outroot, "corr_DM_Pre-Selection_wZpt_Z0j");
    MuonSet     corr_setDM_PreSelection_wZpt_Z1j(outroot, "corr_DM_Pre-Selection_wZpt_Z1j");
    MuonSet     corr_setDM_PreSelection_wZpt_Z2j(outroot, "corr_DM_Pre-Selection_wZpt_Z2j");
   
    // B-Tag Reweight Block 
    BTagEffSet BTagEff_PreSelection(outroot, "BTagEff_Pre-Selection");
    BTagEffSet BTagEff_VBF(outroot, "BTagEff_VBF");

    // Pre-Selection Block
    MuonSet     setDM_PreSelection(outroot, "DM_Pre-Selection");
    MuonSet     corr_setDM_PreSelection(outroot, "corr_DM_Pre-Selection");
    MuonSet     corr_setDM_PreSelection_Z(outroot, "corr_DM_Pre-Selection_Z");
    MuonSet     corr_setDM_PreSelection_SR(outroot, "corr_DM_Pre-Selection_SR");
    JetSet      setJet_PreSelection(outroot, "Jet_Pre-Selection");
    JetSet      corr_setJet_PreSelection(outroot, "corr_Jet_Pre-Selection");
    MuonSet     setDM_PreSelection_BB(outroot, "DM_Pre-Selection_BB");
    MuonSet     corr_setDM_PreSelection_BB(outroot, "corr_DM_Pre-Selection_BB");
    JetSet      setJet_PreSelection_BB(outroot, "Jet_Pre-Selection_BB");
    JetSet      corr_setJet_PreSelection_BB(outroot, "corr_Jet_Pre-Selection_BB");
    MuonSet     setDM_PreSelection_BO(outroot, "DM_Pre-Selection_BO");
    MuonSet     corr_setDM_PreSelection_BO(outroot, "corr_DM_Pre-Selection_BO");
    JetSet      setJet_PreSelection_BO(outroot, "Jet_Pre-Selection_BO");
    JetSet      corr_setJet_PreSelection_BO(outroot, "corr_Jet_Pre-Selection_BO");
    MuonSet     setDM_PreSelection_BE(outroot, "DM_Pre-Selection_BE");
    MuonSet     corr_setDM_PreSelection_BE(outroot, "corr_DM_Pre-Selection_BE");
    JetSet      setJet_PreSelection_BE(outroot, "Jet_Pre-Selection_BE");
    JetSet      corr_setJet_PreSelection_BE(outroot, "corr_Jet_Pre-Selection_BE");
    MuonSet     setDM_PreSelection_OO(outroot, "DM_Pre-Selection_OO");
    MuonSet     corr_setDM_PreSelection_OO(outroot, "corr_DM_Pre-Selection_OO");
    JetSet      setJet_PreSelection_OO(outroot, "Jet_Pre-Selection_OO");
    JetSet      corr_setJet_PreSelection_OO(outroot, "corr_Jet_Pre-Selection_OO");
    MuonSet     setDM_PreSelection_OE(outroot, "DM_Pre-Selection_OE");
    MuonSet     corr_setDM_PreSelection_OE(outroot, "corr_DM_Pre-Selection_OE");
    JetSet      setJet_PreSelection_OE(outroot, "Jet_Pre-Selection_OE");
    JetSet      corr_setJet_PreSelection_OE(outroot, "corr_Jet_Pre-Selection_OE");
    MuonSet     setDM_PreSelection_EE(outroot, "DM_Pre-Selection_EE");
    MuonSet     corr_setDM_PreSelection_EE(outroot, "corr_DM_Pre-Selection_EE");
    JetSet      setJet_PreSelection_EE(outroot, "Jet_Pre-Selection_EE");
    JetSet      corr_setJet_PreSelection_EE(outroot, "corr_Jet_Pre-Selection_EE");
    
    // Pre-DiJet Block
    JetSet      setJet_PreDiJet_VBF(outroot, "Jet_PreDiJet_VBF");
    JetSet      corr_setJet_PreDiJet_VBF(outroot, "corr_Jet_PreDiJet_VBF");
    
    // VBF Block
    MuonSet     setDM_VBF(outroot, "DM_VBF");
    MuonSet     corr_setDM_VBF(outroot, "corr_DM_VBF");
    MuonSet     corr_setDM_VBF_Z(outroot, "corr_DM_VBF_Z");
    MuonSet     corr_setDM_VBF_SR(outroot, "corr_DM_VBF_SR");
    JetSet      setJet_VBF(outroot, "Jet_VBF");
    JetSet      corr_setJet_VBF(outroot, "corr_Jet_VBF");
    MuonSet     setDM_VBF_BB(outroot, "DM_VBF_BB");
    MuonSet     corr_setDM_VBF_BB(outroot, "corr_DM_VBF_BB");
    JetSet      setJet_VBF_BB(outroot, "Jet_VBF_BB");
    JetSet      corr_setJet_VBF_BB(outroot, "corr_Jet_VBF_BB");
    MuonSet     setDM_VBF_BO(outroot, "DM_VBF_BO");
    MuonSet     corr_setDM_VBF_BO(outroot, "corr_DM_VBF_BO");
    JetSet      setJet_VBF_BO(outroot, "Jet_VBF_BO");
    JetSet      corr_setJet_VBF_BO(outroot, "corr_Jet_VBF_BO");
    MuonSet     setDM_VBF_BE(outroot, "DM_VBF_BE");
    MuonSet     corr_setDM_VBF_BE(outroot, "corr_DM_VBF_BE");
    JetSet      setJet_VBF_BE(outroot, "Jet_VBF_BE");
    JetSet      corr_setJet_VBF_BE(outroot, "corr_Jet_VBF_BE");
    MuonSet     setDM_VBF_OO(outroot, "DM_VBF_OO");
    MuonSet     corr_setDM_VBF_OO(outroot, "corr_DM_VBF_OO");
    JetSet      setJet_VBF_OO(outroot, "Jet_VBF_OO");
    JetSet      corr_setJet_VBF_OO(outroot, "corr_Jet_VBF_OO");
    MuonSet     setDM_VBF_OE(outroot, "DM_VBF_OE");
    MuonSet     corr_setDM_VBF_OE(outroot, "corr_DM_VBF_OE");
    JetSet      setJet_VBF_OE(outroot, "Jet_VBF_OE");
    JetSet      corr_setJet_VBF_OE(outroot, "corr_Jet_VBF_OE");
    MuonSet     setDM_VBF_EE(outroot, "DM_VBF_EE");
    MuonSet     corr_setDM_VBF_EE(outroot, "corr_DM_VBF_EE");
    JetSet      setJet_VBF_EE(outroot, "Jet_VBF_EE");
    JetSet      corr_setJet_VBF_EE(outroot, "corr_Jet_VBF_EE");
    
    HistSetCollection hists = {
    
    // Pileup Block
    PileupReweighting_PreSelection, PileupReweighting_VBF,

    // Z-Reweight Block
    corr_setDM_PreSelection_Z0j, corr_setDM_PreSelection_Z1j, corr_setDM_PreSelection_Z2j, 
    corr_setDM_PreSelection_wZpt_Z0j, corr_setDM_PreSelection_wZpt_Z1j, corr_setDM_PreSelection_wZpt_Z2j, 
    
    // B-Tag Reweight Block
    BTagEff_PreSelection, BTagEff_VBF,

    // Pre-Selection Block
    setDM_PreSelection, corr_setDM_PreSelection, corr_setDM_PreSelection_Z, corr_setDM_PreSelection_SR, setJet_PreSelection, corr_setJet_PreSelection,
    setDM_PreSelection_BB, corr_setDM_PreSelection_BB, setJet_PreSelection_BB, corr_setJet_PreSelection_BB,
    setDM_PreSelection_BO, corr_setDM_PreSelection_BO, setJet_PreSelection_BO, corr_setJet_PreSelection_BO,
    setDM_PreSelection_BE, corr_setDM_PreSelection_BE, setJet_PreSelection_BE, corr_setJet_PreSelection_BE,
    setDM_PreSelection_OO, corr_setDM_PreSelection_OO, setJet_PreSelection_OO, corr_setJet_PreSelection_OO,
    setDM_PreSelection_OE, corr_setDM_PreSelection_OE, setJet_PreSelection_OE, corr_setJet_PreSelection_OE,
    setDM_PreSelection_EE, corr_setDM_PreSelection_EE, setJet_PreSelection_EE, corr_setJet_PreSelection_EE,
    
    // Pre-DiJet Block
    setJet_PreDiJet_VBF, corr_setJet_PreDiJet_VBF,

    // VBF Block
    setDM_VBF,    corr_setDM_VBF,    corr_setDM_VBF_Z, corr_setDM_VBF_SR, setJet_VBF, corr_setJet_VBF,
    setDM_VBF_BB, corr_setDM_VBF_BB, setJet_VBF_BB, corr_setJet_VBF_BB, 
    setDM_VBF_BO, corr_setDM_VBF_BO, setJet_VBF_BO, corr_setJet_VBF_BO,
    setDM_VBF_BE, corr_setDM_VBF_BE, setJet_VBF_BE, corr_setJet_VBF_BE,
    setDM_VBF_OO, corr_setDM_VBF_OO, setJet_VBF_OO, corr_setJet_VBF_OO,
    setDM_VBF_OE, corr_setDM_VBF_OE, setJet_VBF_OE, corr_setJet_VBF_OE,
    setDM_VBF_EE, corr_setDM_VBF_EE, setJet_VBF_EE, corr_setJet_VBF_EE
    };

	Streamer streamer(__inputfilename, "Events");
	streamer.chainup();
   
    //MVA Branch outputs 
    gMvaTree = new TTree("mva", "Training features for VBF discriminator");
    gMvaTree->Branch("z",               &g_z);
    gMvaTree->Branch("pt_balance",      &g_ptBalance);
    gMvaTree->Branch("mu1_pt_over_m",   &g_mu1_pt_over_m);
    gMvaTree->Branch("mu2_pt_over_m",   &g_mu2_pt_over_m);
    gMvaTree->Branch("mu1_eta",         &g_mu1_eta);
    gMvaTree->Branch("mu2_eta",         &g_mu2_eta);
    gMvaTree->Branch("m_mumu",          &g_m_mumu);
    gMvaTree->Branch("y_mumu",          &g_y_mumu);
    gMvaTree->Branch("pt_mumu",         &g_pt_mumu);
    
    gMvaTree->Branch("jet1_eta",        &g_jet1_eta);
    gMvaTree->Branch("deta_jj",         &g_deta_jj);
    gMvaTree->Branch("dphi_jj",         &g_dphi_jj);
    gMvaTree->Branch("y_j1",            &g_y_j1);
    gMvaTree->Branch("y_j2",            &g_y_j2);
    gMvaTree->Branch("njets",           &g_njets);
    gMvaTree->Branch("pt_j1",           &g_pt_j1);
    gMvaTree->Branch("pt_j2",           &g_pt_j2);
    gMvaTree->Branch("dphi_mumu_jj",    &g_dphi_mumu_jj, "dphi_mumu_jj/F");
    gMvaTree->Branch("deta_mumu_jj",    &g_deta_mumu_jj, "deta_mumu_jj/F");
    gMvaTree->Branch("phi_CS",          &g_phi_CS,      "phi_CS/F");
    gMvaTree->Branch("costheta_CS",     &g_costheta_CS, "costheta_CS/F");
    
    gMvaTree->Branch("weight",          &g_weight);

    Muons *muons=NULL;
	Muons muons1;
	Muons muons2;
	Jets *jets=NULL;
    GenJets *genjets=NULL;
    Photons *photons=NULL;
    FsrPhotons *fsrphotons=NULL;
	Vertices *vertices=NULL;
	Event *event=NULL;
	EventAuxiliary *aux=NULL;
	MET *met=NULL;
    L1s *l1s=NULL;
    HLTs *hlts=NULL;
    TrigObjs *trigobjs=NULL;
    Electrons *electrons=NULL;
    Rhos *rhos=NULL;
    Flags *flags=NULL;

#define ROOT_CHAIN streamer._chain
    
    //Setting Branch Array for LuminosityBlocks Contents
    SET_BRANCH_UINT(run);
    SET_BRANCH_UINT(luminosityBlock);
    ULong64_t eventNumber = 0;
    ROOT_CHAIN->SetBranchAddress("event", &eventNumber);

    //Setting Branch Array for Muon Properties
    SET_BRANCH_UINT(nMuon);
    SET_BRANCH_FLOAT_ARRAY(Muon_pt);
    SET_BRANCH_FLOAT_ARRAY(Muon_bsConstrainedPt);
    SET_BRANCH_FLOAT_ARRAY(Muon_bsConstrainedChi2);
    SET_BRANCH_FLOAT_ARRAY(Muon_bsConstrainedPtErr);
    SET_BRANCH_FLOAT_ARRAY(Muon_phi);
    SET_BRANCH_FLOAT_ARRAY(Muon_eta);
    SET_BRANCH_INT_ARRAY(Muon_charge);
    SET_BRANCH_BOOL_ARRAY(Muon_isGlobal);
    SET_BRANCH_BOOL_ARRAY(Muon_isTracker);
    SET_BRANCH_UCHAR_ARRAY(Muon_miniIsoId);
    SET_BRANCH_BOOL_ARRAY(Muon_tightId);
    SET_BRANCH_BOOL_ARRAY(Muon_mediumId);
    SET_BRANCH_BOOL_ARRAY(Muon_looseId);
    SET_BRANCH_BOOL_ARRAY(Muon_isPFcand);
    SET_BRANCH_FLOAT_ARRAY(Muon_dxy);
    SET_BRANCH_FLOAT_ARRAY(Muon_dz);
    SET_BRANCH_FLOAT_ARRAY(Muon_pfRelIso04_all);
    SET_BRANCH_UCHAR_ARRAY(Muon_nTrackerLayers);

    //Setting Branch Array for Primary Vertex Properties
    SET_BRANCH_UCHAR(PV_npvs);
    SET_BRANCH_UCHAR(PV_npvsGood);
    SET_BRANCH_FLOAT(PV_ndof);
    SET_BRANCH_FLOAT(PV_chi2);
    SET_BRANCH_FLOAT(PV_score);
    SET_BRANCH_FLOAT(PV_z);
    SET_BRANCH_FLOAT(PV_y);
    SET_BRANCH_FLOAT(PV_x);

    //Setting Branch Array for Missing Transverse Energy Properties
    SET_BRANCH_FLOAT(MET_pt);

    //Setting Branch Array for Jet Properties
    SET_BRANCH_UINT(nJet);
    SET_BRANCH_FLOAT_ARRAY(Jet_eta);
    SET_BRANCH_FLOAT_ARRAY(Jet_mass);
    SET_BRANCH_FLOAT_ARRAY(Jet_pt);
    SET_BRANCH_FLOAT_ARRAY(Jet_phi);
    SET_BRANCH_UCHAR_ARRAY(Jet_jetId);
    SET_BRANCH_FLOAT_ARRAY(Jet_btagDeepFlavB);
    SET_BRANCH_FLOAT_ARRAY(Jet_btagPNetB);
    SET_BRANCH_FLOAT_ARRAY(Jet_btagRobustParTAK4B);
    SET_BRANCH_UCHAR_ARRAY(Jet_hadronFlavour);
    SET_BRANCH_FLOAT_ARRAY(Jet_area);
    SET_BRANCH_FLOAT_ARRAY(Jet_chHEF);
    SET_BRANCH_FLOAT_ARRAY(Jet_neHEF);
    SET_BRANCH_FLOAT_ARRAY(Jet_chEmEF);
    SET_BRANCH_FLOAT_ARRAY(Jet_neEmEF);
    SET_BRANCH_FLOAT_ARRAY(Jet_muEF);
    SET_BRANCH_USHORT_ARRAY(Jet_genJetIdx);

    //Setting Branch Array for GenJet Properties
    SET_BRANCH_UINT(nGenJet);
    SET_BRANCH_FLOAT_ARRAY(GenJet_eta);
    SET_BRANCH_FLOAT_ARRAY(GenJet_mass);
    SET_BRANCH_FLOAT_ARRAY(GenJet_pt);
    SET_BRANCH_FLOAT_ARRAY(GenJet_phi);
    
    //Setting Branch Array for FsrPhoton Properties
    SET_BRANCH_UINT(nFsrPhoton);
    SET_BRANCH_FLOAT_ARRAY(FsrPhoton_pt);
    SET_BRANCH_FLOAT_ARRAY(FsrPhoton_phi);
    SET_BRANCH_FLOAT_ARRAY(FsrPhoton_eta);
    SET_BRANCH_INT_ARRAY(FsrPhoton_muonIdx);
    SET_BRANCH_FLOAT_ARRAY(FsrPhoton_dROverEt2);
    SET_BRANCH_USHORT_ARRAY(FsrPhoton_electronIdx);

    //Setting Branch Array for Photon Properties
    SET_BRANCH_UINT(nPhoton);
    SET_BRANCH_FLOAT_ARRAY(Photon_pt);
    SET_BRANCH_FLOAT_ARRAY(Photon_eta);
    SET_BRANCH_FLOAT_ARRAY(Photon_phi);
    SET_BRANCH_FLOAT_ARRAY(Photon_hoe);
    SET_BRANCH_FLOAT_ARRAY(Photon_mvaID);
    SET_BRANCH_BOOL_ARRAY(Photon_isScEtaEE);
    SET_BRANCH_BOOL_ARRAY(Photon_isScEtaEB);
    SET_BRANCH_BOOL_ARRAY(Photon_pixelSeed);

    //Setting Branch Array for L1 Properties
    SET_BRANCH_BOOL_ARRAY(L1_DoubleMu_12_5);
    SET_BRANCH_BOOL_ARRAY(L1_SingleMu12er1p5);

    //Setting Branch Array for HLT Properties
    SET_BRANCH_BOOL_ARRAY(HLT_IsoMu24);

    //Setting Storage for TrigObj Properties 
    std::vector<Float_t>*   TrigObj_pt          = nullptr;
    std::vector<Float_t>*   TrigObj_eta         = nullptr;
    std::vector<Float_t>*   TrigObj_phi         = nullptr;
    std::vector<UShort_t>*  TrigObj_id          = nullptr;
    std::vector<ULong64_t>* TrigObj_filterBits  = nullptr;

    static const int TO_MAX = 512; // safe cap for leaf-list mode
    static Float_t   TrigObj_pt_arr[TO_MAX];
    static Float_t   TrigObj_eta_arr[TO_MAX];
    static Float_t   TrigObj_phi_arr[TO_MAX];
    static UShort_t  TrigObj_id_arr[TO_MAX];
    static ULong64_t TrigObj_filterBits_arr[TO_MAX];

    UInt_t nTrigObj = 0;
    bool   trigobj_is_vector = false;

    //Setting Branch Array for Electron Properties
    SET_BRANCH_UINT(nElectron);
    SET_BRANCH_FLOAT_ARRAY(Electron_pt);
    SET_BRANCH_FLOAT_ARRAY(Electron_eta);
    SET_BRANCH_FLOAT_ARRAY(Electron_phi);
    SET_BRANCH_UCHAR_ARRAY(Electron_cutBased);
    SET_BRANCH_BOOL_ARRAY(Electron_mvaIso_WP90);

    //Setting Branch Array for Rho Properties
    SET_BRANCH_FLOAT(fixedGridRhoFastjetAll);
   
    //Setting Branch Array for Rho Properties
    SET_BRANCH_BOOL(Flag_goodVertices);
    SET_BRANCH_BOOL(Flag_globalSuperTightHalo2016Filter);
    SET_BRANCH_BOOL(Flag_EcalDeadCellTriggerPrimitiveFilter);
    SET_BRANCH_BOOL(Flag_BadPFMuonFilter);
    SET_BRANCH_BOOL(Flag_BadPFMuonDzFilter);
    SET_BRANCH_BOOL(Flag_eeBadScFilter);
    SET_BRANCH_BOOL(Flag_ecalBadCalibFilter);
    SET_BRANCH_BOOL(Flag_hfNoisyHitsFilter);
   
    // --- Event-level weights from NanoAOD ---
    float Pileup_nTrueInt  = 0.0f;

    float genWeight_branch = 1.0f;
    float Generator_weight_branch = 1.0f;

    if (ROOT_CHAIN->GetBranch("genWeight"))
        ROOT_CHAIN->SetBranchAddress("genWeight", &genWeight_branch);
    if (ROOT_CHAIN->GetBranch("Generator_weight"))
        ROOT_CHAIN->SetBranchAddress("Generator_weight", &Generator_weight_branch);

    // Choose the one that matches Runs/genEventSumw semantics
    auto pickGenW = [&]() -> float {
        if (ROOT_CHAIN->GetBranch("genWeight")) return genWeight_branch;            // preferred
        if (ROOT_CHAIN->GetBranch("Generator_weight")) return Generator_weight_branch; // fallback
        return 1.0f;
    };

    if (ROOT_CHAIN->GetBranch("Pileup_nTrueInt"))
        ROOT_CHAIN->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt);
    else if (ROOT_CHAIN->GetBranch("nTrueInt"))
        ROOT_CHAIN->SetBranchAddress("nTrueInt", &Pileup_nTrueInt);


    muons       = new std::vector<Muon>();
    jets        = new std::vector<Jet>();
    genjets     = new std::vector<GenJet>();
    fsrphotons  = new std::vector<FsrPhoton>();
    photons     = new std::vector<Photon>();
    vertices    = new std::vector<Vertex>();
    met         = new MET();
    l1s         = new std::vector<L1>();
    hlts        = new std::vector<HLT>();
    trigobjs    = new std::vector<TrigObj>();
    electrons   = new std::vector<Electron>();
    rhos        = new std::vector<Rho>();
    flags       = new std::vector<Flag>();
    event       = new Event();


	//	Main Event Loop
	uint32_t numEntries         = streamer._chain->GetEntries();
    uint32_t updateInterval     = numEntries / 100; // Update every 1% of progress
    if (updateInterval == 0) updateInterval = 1; // Handle very small numbers of entries
    
    double sum_w_nominal = 0.0;      // base * genWeight
    double sum_w_pu      = 0.0;           // base * genWeight * PU
    long long n_events   = 0;
    long long n_neg      = 0;

    // Ensure statuses are sane
    ROOT_CHAIN->SetBranchStatus("TrigObj_*", 0);
    ROOT_CHAIN->SetBranchStatus("nTrigObj", 1);
    ROOT_CHAIN->SetBranchStatus("TrigObj_pt", 1);
    ROOT_CHAIN->SetBranchStatus("TrigObj_eta", 1);
    ROOT_CHAIN->SetBranchStatus("TrigObj_phi", 1);
    ROOT_CHAIN->SetBranchStatus("TrigObj_id", 1);
    ROOT_CHAIN->SetBranchStatus("TrigObj_filterBits", 1);

    // Counter is always safe to bind
    ROOT_CHAIN->SetBranchAddress("nTrigObj", &nTrigObj);

    // Detect storage kind
    auto* br_pt = ROOT_CHAIN->GetBranch("TrigObj_pt");
    trigobj_is_vector = br_pt && br_pt->GetClassName()
                      && std::string(br_pt->GetClassName()).find("vector") != std::string::npos;

    if (trigobj_is_vector) {
        ROOT_CHAIN->SetBranchAddress("TrigObj_pt",         &TrigObj_pt);
        ROOT_CHAIN->SetBranchAddress("TrigObj_eta",        &TrigObj_eta);
        ROOT_CHAIN->SetBranchAddress("TrigObj_phi",        &TrigObj_phi);
        ROOT_CHAIN->SetBranchAddress("TrigObj_id",         &TrigObj_id);
        ROOT_CHAIN->SetBranchAddress("TrigObj_filterBits", &TrigObj_filterBits);
        std::cerr << "[INFO] TrigObj binding = VECTOR mode\n";
    } else {
        ROOT_CHAIN->SetBranchAddress("TrigObj_pt",         TrigObj_pt_arr);
        ROOT_CHAIN->SetBranchAddress("TrigObj_eta",        TrigObj_eta_arr);
        ROOT_CHAIN->SetBranchAddress("TrigObj_phi",        TrigObj_phi_arr);
        ROOT_CHAIN->SetBranchAddress("TrigObj_id",         TrigObj_id_arr);
        ROOT_CHAIN->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits_arr);
        std::cerr << "[INFO] TrigObj binding = LEAFLIST mode (counter=nTrigObj, max=" << TO_MAX << ")\n";
    }

std::cout << "[BRANCH CHECK] Flag_goodVertices exists = "
          << (ROOT_CHAIN->GetBranch("Flag_goodVertices") != nullptr) << "\n";
std::cout << "[BRANCH CHECK] PV_npvsGood exists = "
          << (ROOT_CHAIN->GetBranch("PV_npvsGood") != nullptr) << "\n";


// --- DEBUG COUNTERS ---
long long n_afterGolden = 0;
long long n_afterFlags  = 0;
long long n_afterVtx    = 0;
long long n_afterVeto   = 0;
long long n_afterPair   = 0;

    for (uint32_t i=0; i<numEntries && __continueRunning; i++)
	{   
        muons->clear();
        jets->clear();
        genjets->clear();
        photons->clear();
        fsrphotons->clear();
        vertices->clear();
        muons1.clear(); 
        muons2.clear();
		streamer._chain->GetEntry(i);
        l1s->clear();
        hlts->clear();
        trigobjs->clear();
        electrons->clear();
        rhos->clear();

        ROOT_CHAIN->GetEntry(i);
        event->_run  = run;                
        event->_lumi = luminosityBlock;    
        const uint32_t evt  = (uint32_t)(eventNumber & 0xffffffffULL);
        const uint32_t lumi = (uint32_t)luminosityBlock;

        for (size_t i=0; i<nMuon; i++) {
            Muon mu;
            mu._pt                  = Muon_pt[i];
            mu._bsConstrainedPt     = Muon_bsConstrainedPt[i];
            mu._bsConstrainedChi2   = Muon_bsConstrainedChi2[i];
            mu._bsConstrainedPtErr  = Muon_bsConstrainedPtErr[i];
            mu._eta             	= Muon_eta[i];
            mu._phi             	= Muon_phi[i];
            mu._charge          	= Muon_charge[i];
            mu._isGlobal        	= Muon_isGlobal[i];
            mu._isTracker       	= Muon_isTracker[i];
            mu._tightId         	= Muon_tightId[i];
            mu._mediumId        	= Muon_mediumId[i];
            mu._looseId         	= Muon_looseId[i];
            mu._isPFcand        	= Muon_isPFcand[i];
            mu._dxy             	= Muon_dxy[i];
            mu._dz              	= Muon_dz[i];
            mu._pfRelIso04_all  	= Muon_pfRelIso04_all[i];
            mu._nTrackerLayers  	= Muon_nTrackerLayers[i];

            muons->push_back(std::move(mu));
        }

        for (size_t i=0; i<nJet; i++) {
            Jet jet;
            jet._pt 				= Jet_pt[i];
            jet._eta 				= Jet_eta[i];
            jet._phi 				= Jet_phi[i];
            jet._mass 				= Jet_mass[i];
            jet._jetId 				= Jet_jetId[i];
            jet._btagDeepFlavB		= Jet_btagDeepFlavB[i];
            jet._btagPNetB 			= Jet_btagPNetB[i];
            jet._btagRobustParTAK4B	= Jet_btagRobustParTAK4B[i];
            jet._hadronFlavour 		= Jet_hadronFlavour[i];
            jet._area 				= Jet_area[i];
            jet._chHEF 				= Jet_chHEF[i];
            jet._neHEF 				= Jet_neHEF[i];
            jet._chEmEF 			= Jet_chEmEF[i];
            jet._neEmEF 			= Jet_neEmEF[i];
            jet._muEF 				= Jet_muEF[i];
            jet._genJetIdx 			= Jet_genJetIdx[i];

            jets->push_back(std::move(jet));
        }
        
        if (__isMC) {
            for (size_t i=0; i<nGenJet; i++) {
                GenJet genjet;
                genjet._pt = GenJet_pt[i];
                genjet._eta = GenJet_eta[i];
                genjet._phi = GenJet_phi[i];
                genjet._mass = GenJet_mass[i];
                
                genjets->push_back(std::move(genjet));
            }
        }

        for (size_t i=0; i<nFsrPhoton; i++){
            FsrPhoton fsrphoton;
            fsrphoton._pt = FsrPhoton_pt[i];
            fsrphoton._eta = FsrPhoton_eta[i];
            fsrphoton._phi = FsrPhoton_phi[i];
            fsrphoton._muonIdx = FsrPhoton_muonIdx[i];
            fsrphoton._dROverEt2 = FsrPhoton_dROverEt2[i];
            fsrphoton._electronIdx = FsrPhoton_electronIdx[i];

            fsrphotons->push_back(std::move(fsrphoton));
        }

        const size_t nto = trigobj_is_vector
            ? (TrigObj_pt ? TrigObj_pt->size() : 0)
            : std::min<size_t>((size_t)nTrigObj, (size_t)TO_MAX);

        for (size_t i = 0; i < nto; ++i) {
            TrigObj o;
            o._pt  = trigobj_is_vector ? (*TrigObj_pt)[i]         : TrigObj_pt_arr[i];
            o._eta = trigobj_is_vector ? (*TrigObj_eta)[i]        : TrigObj_eta_arr[i];
            o._phi = trigobj_is_vector ? (*TrigObj_phi)[i]        : TrigObj_phi_arr[i];
            o._id  = trigobj_is_vector ? (short)(*TrigObj_id)[i]  : (short)TrigObj_id_arr[i];
            // TrigObj._filterBits is an int; the branch is 64-bit -> implicit narrowing here
            o._filterBits = trigobj_is_vector ? (int)(*TrigObj_filterBits)[i] : (int)TrigObj_filterBits_arr[i];

            trigobjs->push_back(std::move(o));
        }
        
        for (size_t i=0; i<nElectron; i++) {
            Electron electron;
            electron._pt          = Electron_pt[i];
            electron._eta         = Electron_eta[i];
            electron._phi         = Electron_phi[i]; 
            electron._cutBased    = Electron_cutBased[i];
            electron._mvaIso_WP90 = Electron_mvaIso_WP90[i];

            electrons->push_back(std::move(electron));
        }

        Rho rhoObj;
        rhoObj._fixedGridRhoFastjetAll = fixedGridRhoFastjetAll; 
        rhos->push_back(rhoObj);

        float rho = rhos->at(0)._fixedGridRhoFastjetAll;

/*static int dbg1 = 0;
if (dbg1 < 15) {
std::cout << "[STEP1 - BRANCH] PV_npvsGood=" << PV_npvsGood
<< " PV_npvs=" << PV_npvs
<< " evt=" << evt
<< std::endl;
++dbg1;
}*/

        Vertex vertex;
        vertex._chi2        = PV_chi2;
        vertex._ndof        = PV_ndof;
        vertex._score       = PV_score;
        vertex._x           = PV_x;
        vertex._y           = PV_y;
        vertex._z           = PV_z;
        vertex._npvsGood    = static_cast<int>(PV_npvsGood);
        vertex._npvs        = static_cast<int>(PV_npvs);

        vertices->push_back(std::move(vertex));

/*static int dbg2 = 0;
if (dbg2 < 10 && !vertices->empty()) {
std::cout << "[STEP2 - STORED] vertices->_npvsGood="
<< vertices->back()._npvsGood
<< std::endl;
++dbg2;
}*/

        Flag flag;
        flag._goodVertices                          = Flag_goodVertices;
        flag._globalSuperTightHalo2016Filter        = Flag_globalSuperTightHalo2016Filter;
        flag._EcalDeadCellTriggerPrimitiveFilter    = Flag_EcalDeadCellTriggerPrimitiveFilter;
        flag._BadPFMuonFilter                       = Flag_BadPFMuonFilter;
        flag._BadPFMuonDzFilter                     = Flag_BadPFMuonDzFilter;
        flag._eeBadScFilter                         = Flag_eeBadScFilter;       // DATA-only
        flag._ecalBadCalibFilter                    = Flag_ecalBadCalibFilter;
        flag._hfNoisyHitsFilter                     = Flag_hfNoisyHitsFilter;

        flags->push_back(std::move(flag));
       
        // Pileup from correctionlib (JSON expects: NumTrueInteractions, "weights")
        float pileupWeight = 1.0f;
        if (__isMC && pu_corrections && !g_disablePUWeights) {
            try {
                pileupWeight = static_cast<float>(
                    pu_corrections->evaluate({static_cast<double>(Pileup_nTrueInt), "nominal"}));
            } catch (const std::exception& e) {
                std::cerr << "[PU ERROR] evaluate failed: " << e.what() << std::endl;
                pileupWeight = 1.0f;
            }
        }

        float genWeight = __isMC ? pickGenW() : 1.0f;
        double base = (__isMC ? sampleBaseWeight(g_sampleName) : 1.0);
        float Dataweight = static_cast<float>(base) * pileupWeight * genWeight;
        
        const double w_nom = base * genWeight;
        const double w_pu  = w_nom * pileupWeight;
        

sum_w_nominal += w_nom;
sum_w_pu      += w_pu;
n_events++;
if (genWeight < 0) n_neg++;


/*if (i < 50) {
std::cout << std::setprecision(10) 
<< " genWeight=" << genWeight
<< " Generator_weight=" << Generator_weight_branch
<< " Pileup_nTrueInt=" << Pileup_nTrueInt
<< " pileupWeight=" << pileupWeight
<< " base=" << base
<< " total=" << Dataweight << "\n";
}*/

       if (!__isMC && !g_lumiMask.empty()) {
            const unsigned run  = static_cast<unsigned>(event->_run);               // adapt to your names
            const unsigned lumi = static_cast<unsigned>(event->_lumi);              // e.g. event._lumi
            if (!passGolden(run, lumi)) {
                continue; // skip non-certified lumisections
            }
       }
       ++n_afterGolden;

        static int dbg = 0;
if (!__isMC && !g_lumiMask.empty() && dbg < 10) {
unsigned run  = static_cast<unsigned>(event->_run);
unsigned lumi = static_cast<unsigned>(event->_lumi);
std::cout << "[GOLDEN DBG] run=" << run
<< " lumi=" << lumi
<< " pass=" << passGolden(run, lumi) << "\n";
++dbg;
}

/*static int flagDbg = 0;
if (flagDbg < 10 && flags && !flags->empty()) {
const auto& f = flags->front();
std::cout
<< "[FLAG DBG] "
<< "goodVertices=" << f._goodVertices
<< " globalSuperTightHalo2016Filter=" << f._globalSuperTightHalo2016Filter
<< " EcalDeadCellTriggerPrimitiveFilter=" << f._EcalDeadCellTriggerPrimitiveFilter
<< " BadPFMuonFilter=" << f._BadPFMuonFilter
<< " BadPFMuonDzFilter=" << f._BadPFMuonDzFilter
<< " hfNoisyHitsFilter=" << f._hfNoisyHitsFilter
<< " ecalBadCalibFilter=" << f._ecalBadCalibFilter
<< " eeBadScFilter=" << f._eeBadScFilter
<< std::endl;
++flagDbg;
}*/

        //	Selections
        if (!passFlagFilters(flags, !__isMC)) continue;
        ++n_afterFlags;
        if (!passVertex(vertices)) continue;
        ++n_afterVtx;
        if (!vetoAdditionalLeptons(muons, electrons, __isMC)) continue;
        ++n_afterVeto;

        bool anyPair = false;
        // Prepare The Pairs Of Muons
        for (auto it = muons->begin(); it != muons->end(); ++it) {
            for (auto jt = (it + 1); jt != muons->end(); ++jt) {
            
            if (!passMuons(*it, *jt, trigobjs, evt, lumi)) continue;
            anyPair = true;

            categorize(jets, genjets, *it, *jt, fsrphotons, photons, vertices, *met, *event, rhos, hists, Dataweight, evt, lumi, w_nom, w_pu);
            break;
        }
        if (anyPair) break;
    }
    if (!anyPair) continue;
    ++n_afterPair;
    
    // Update the progress bar every `updateInterval` iterations
    if (i % updateInterval == 0 || i == numEntries - 1) {
            showProgressBar(i + 1, numEntries);
        }
    }

std::cout << "\n[CUTFLOW] Total Events = "      << n_events          << " | "
<< " afterGolden = "                     << n_afterGolden     << " | "
<< " afterFlags = "                      << n_afterFlags      << " | "
<< " afterVertex = "                     << n_afterVtx        << " | "
<< " afterVetoAddLeptons = "             << n_afterVeto       << " | "
<< " afterPair = "                       << n_afterPair       
<< std::endl;
   
std::cout << "\n[JET CUT-FLOW] Total Jets = " << jets_in << " | "
<< "pass jetId = " << jets_after_id << " | "
<< "pass pt(>25) = " << jets_after_pt_default << " | "
<< "pass HF horn (pT>50) = " << jets_after_pt_hf << " | "
<< "pass HE horn (pT>50) = " << jets_after_pt_he << " | "
<< "pass pt(final) = " << jets_after_pt << " | "
<< "pass eta = " << jets_after_eta << " | "
<< "pass dR1 = " << jets_after_dR1 << " | "
<< "pass dR2 = " << jets_after_dR2
<< std::endl;

    std::cout << std::endl; // Move to the next line after the progress bar is complete

	outroot->Write();
	outroot->Close();

	return;
}

void sigHandler(int sig)
{
	cout << "### Signal: " << sig << " caughter. Exiting..." << endl;
	__continueRunning = false;
}

int main(int argc, char** argv)
{
	// Register Signals
	signal(SIGABRT, &sigHandler);
	signal(SIGTERM, &sigHandler);
	signal(SIGINT, &sigHandler);

	std::string none;
	std::string muCorrOpt = "roch";

	// Pare Options
	po::options_description desc("Allowed Program Options");
	desc.add_options()
		("help", "produce help messages")
		("input",           po::value<std::string>(),                                           "a file specifying all the ROOT files to process")
		("isMC",            po::value<bool>(),                                                  "type of data: DATA vs MC")
        ("muCorr",          po::value<std::string>(&muCorrOpt)->default_value("roch"),          "'roch' (Rochester), 'kit' (KIT smearing), 'none' (NoMuCorr)")
        ("output",          po::value<std::string>(),                                           "output file name")
        ("noJetCriteria",  po::bool_switch(&g_disableJetCriteria)->default_value(false),        "Disable jet kinematic/ID/muon-cleaning cuts")
        ("noJERC",          po::bool_switch(&g_disableJetJERC)->default_value(false),           "Disable jet JEC/JER from jet_jerc*.json (use raw jet pT)")
        ("noJetVeto",       po::bool_switch(&g_disableJetVeto)->default_value(false),           "Disable jet veto maps (keep all jets; ignore jet_veto_maps JSON)")
        ("noBtagSF",        po::bool_switch(&g_disableBtagSF)->default_value(false),            "Disable b-tagging scale factors (event btag weight = 1)")
        ("noMuonSF",        po::bool_switch(&g_disableMuonSF)->default_value(false),            "Disable muon ID/ISO/trigger scale factors (muon SF = 1)")
        ("noScaleSmear",    po::bool_switch(&g_disableMuonScaleSmear)->default_value(false),    "Disable muon scale/smear JSON factor (only Rochester/KIT left)")
        ("noFlags",         po::bool_switch(&g_disableFlagFilters)->default_value(false),       "Disable event Flag (MET filter) requirements")
        ("noPU",            po::bool_switch(&g_disablePUWeights)->default_value(false),         "Disable pileup reweighting (Pileup_nTrueInt -> weight = 1)")
        ("noFSR",           po::bool_switch(&g_disableFSRRecovery)->default_value(false),       "Disable FSR photon recovery (use pre-FSR muon four-vectors)")
        ("noZPtRW",         po::bool_switch(&g_disableZPtRW)->default_value(false),             "Disable DY Z-pT reweighting")
        ("noBeamConstraint", po::bool_switch(&g_disableBeamConstraint)->default_value(false),   "Disable beam-constrained muon pT usage (bsConstrainedPt)")
        ("leadmuonPt",      po::value<double>(&_leadmuonPt)->default_value(_leadmuonPt),        "Muon Pt Cut")
        ("subleadmuonPt",   po::value<double>(&_subleadmuonPt)->default_value(_subleadmuonPt),  "Muon Pt Cut")
        ("leadmuonEta",     po::value<double>(&_leadmuonEta)->default_value(_leadmuonEta),      "Muon Eta Cut")
        ("subleadmuonEta",  po::value<double>(&_subleadmuonEta)->default_value(_subleadmuonEta),"Muon Eta Cut")
        ("leadmuonIso",     po::value<double>(&_leadmuonIso)->default_value(_leadmuonIso),      "Muon Isolation Cut")
        ("subleadmuonIso",  po::value<double>(&_subleadmuonIso)->default_value(_subleadmuonIso),"Muon Isolation Cut")
        ("leadJetPt",       po::value<double>(&_leadJetPt)->default_value(_leadJetPt),          "Lead Jet Pt Cut")
        ("subleadJetPt",    po::value<double>(&_subleadJetPt)->default_value(_subleadJetPt),    "SubLeading Jet Pt Cut")
        ("leadJetEta",      po::value<double>(&_leadJetEta)->default_value(_leadJetEta),        "Lead Jet Eta Cut")
        ("subleadJetEta",   po::value<double>(&_subleadJetEta)->default_value(_subleadJetEta),  "SubLeading Jet Eta Cut")
        ("dijetMass_VBF",   po::value<double>(&_dijetMass_VBF)->default_value(_dijetMass_VBF),  "DiJet Mass VBF-Category Cut")
        ("dijetdEta_VBF",   po::value<double>(&_dijetdEta_VBF)->default_value(_dijetdEta_VBF),  "DiJet deta VBF-Category Cut")
	    ("zptRWFile",       po::value<std::string>(&g_zptRWFile)->default_value("None"),        "ROOT file containing ZpT TF1s (f_combined)")
        ("weights",         po::value<std::string>(&g_weightsFile)->default_value("None"),      "path to weights config (key=value)")
        ("sample",          po::value<std::string>(&g_sampleName)->default_value("DATA"),       "sample name (e.g. DY, TT, WZ, VBF, DATA)")
        ("era",             po::value<std::string>(&g_era)->default_value("12022"),             "Era section from weights.ini (e.g. 12022, 22022, 12023, 22023)")
        ("puDir",           po::value<std::string>(&g_puDir)->default_value(g_puDir),           "Directory containing PU JSONs")
        ("puJSON",          po::value<std::string>(&g_puJSON)->default_value("auto"),           "PU JSON file path or 'auto' to use puWeights<era>.json")
        ("puKey",           po::value<std::string>(&g_puKey)->default_value("auto"),            "Correction key in PU JSON or 'auto' to use era-based default")
        ("--muDir",         po::value<std::string>(&g_muDir)->default_value(g_muDir),           "Directory containing muon JSONs")
        ("--muSFJSON",      po::value<std::string>(&g_muSFJSON)->default_value("auto"),         "Muon SF JSON (path) or 'auto' => muon_Z<era>.json")
        ("--muSFKeyID",     po::value<std::string>(&g_muSFKeyID)->default_value(g_muSFKeyID),   "Muon ID SF key")
        ("--muSFKeyISO",    po::value<std::string>(&g_muSFKeyISO)->default_value(g_muSFKeyISO), "Muon ISO SF key")
        ("--muSFKeyTRG",    po::value<std::string>(&g_muSFKeyTRG)->default_value(g_muSFKeyTRG), "Muon trigger SF key")
        ("--rochJSON",      po::value<std::string>(&g_rochJSON)->default_value("auto"),         "Rochester JSON (path) or 'auto' => RoccoR<era>.json")
        ("--rochKeyDT",     po::value<std::string>(&g_rochKeyDT)->default_value(g_rochKeyDT),   "Rochester data key")
        ("--rochKeyMC",     po::value<std::string>(&g_rochKeyMC)->default_value(g_rochKeyMC),   "Rochester MC key")
        ("--jmeDir",        po::value<std::string>(&g_jmeDir)->default_value(g_jmeDir),         "Directory containing JME JSONs")
        ("--jercJSON",      po::value<std::string>(&g_jercJSON)->default_value("auto"),         "Path to JERC JSON, or 'auto' => jet_jerc<era>.json")
        ("--jetVetoJSON",   po::value<std::string>(&g_jetVetoJSON)->default_value("auto"),      "Path to veto map JSON, or 'auto' => jetvetomaps<era>.json")
        ("--jetIdJSON",     po::value<std::string>(&g_jetIdJSON)->default_value("auto"),        "Path to jet ID JSON, or 'auto' => jetid<era>.json")
        ("golden",          po::value<std::string>()->default_value("auto"),                    "Golden JSON path or 'auto' for era-based default")
        ("noGolden",        po::bool_switch()->default_value(false),                            "Disable golden JSON masking on DATA");

    ;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
   
    std::cout << "[CONFIG] Beam constraint pT: " 
              << (g_disableBeamConstraint ? "DISABLED" : "ENABLED")
              << std::endl;

    if (vm.count("help") || argc<2)
	{
		std::cout << desc << std::endl;
		return 1;
	}

    if (!vm.count("output")) {
        std::cerr << "Error: --output option is missing!" 
                  << std::endl;
        exit(1);
    }
	
	//	Assign globals
	__inputfilename     = vm["input"].as<std::string>();
	__isMC              = vm["isMC"].as<bool>();
    __outputfilename    = vm["output"].as<std::string>();
    _leadmuonPt         = vm["leadmuonPt"].as<double>();
    _subleadmuonPt      = vm["subleadmuonPt"].as<double>();
    _leadmuonEta        = vm["leadmuonEta"].as<double>();
    _subleadmuonEta     = vm["subleadmuonEta"].as<double>();
    _leadmuonIso        = vm["leadmuonIso"].as<double>();
    _subleadmuonIso     = vm["subleadmuonIso"].as<double>();
    _leadJetPt          = vm["leadJetPt"].as<double>();
    _subleadJetPt       = vm["subleadJetPt"].as<double>();
    _leadJetEta         = vm["leadJetEta"].as<double>();
    _subleadJetEta      = vm["subleadJetEta"].as<double>();
    _dijetMass_VBF      = vm["dijetMass_VBF"].as<double>();
    _dijetdEta_VBF      = vm["dijetdEta_VBF"].as<double>();
    g_goldenJSON        = vm["golden"].as<std::string>();
    const bool has_muCorr = (vm.count("muCorr") > 0);
    bool noGolden       = vm["noGolden"].as<bool>();


    if (!has_muCorr) {
        if (__noMuCorr) {
            muCorrOpt = "none";
        } else if (__useKIT) {
            muCorrOpt = "kit";
        } else {
            muCorrOpt = "roch";
        }
    }

    std::string v = muCorrOpt;
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);

    // default everything off
    __noMuCorr     = false;
    __useKIT       = false;
    __useRochester = false;

    if (v == "none") {
        __noMuCorr = true;
    } else if (v == "kit") {
        __useKIT = true;
    } else if (v == "roch") {
        __useRochester = true;
    } else {
        // unknown value: fall back to Rochester
        std::cerr << "[MUON] WARNING: unknown --muCorr=" << v << " (using Rochester).\n";
        __useRochester = true;
    }

    std::string mode =
        (__noMuCorr     ? "NoMuCorr"   :
         __useKIT       ? "KIT"        :
         __useRochester ? "Rochester"  :
                          "None");

    std::cout << "[MuonCorr] Mode=" << mode
              << "  isMC=" << std::boolalpha << __isMC
              << "  era="  << g_era << std::endl;

    std::cerr << "[MUON] Correction mode = " << mode << std::endl;

    // --- Load weights & pick lumi for the chosen era (only for MC) ---
    if (__isMC && g_weightsFile != "None") {
        loadWeightsFile(g_weightsFile);   // fills g_sectionLumis / g_sectionSamples

    // pick lumi and remember which key it came from
    const auto& lmap = g_sectionLumis[g_era];
    if (!lmap.empty()) {
        // Prefer a specific "luminosity_*" entry if present, else fall back to the canonical "lumi",
        // else just take the first available.
        auto it_specific = std::find_if(lmap.begin(), lmap.end(),
            [](const auto& kv){ return kv.first.rfind("luminosity", 0) == 0; });

        if (it_specific != lmap.end()) {
            g_selected_lumi_key = it_specific->first;
            g_selected_lumi     = it_specific->second;
        } else if (lmap.count("lumi")) {
            g_selected_lumi_key = "lumi";
            g_selected_lumi     = lmap.at("lumi");
        } else {
            g_selected_lumi_key = lmap.begin()->first;
            g_selected_lumi     = lmap.begin()->second;
        }
    } else {
        g_selected_lumi_key = "unset";
        g_selected_lumi     = 0.0;
    }

std::cout << std::fixed
<< "[WEIGHTS] era=" << g_era
<< " sample=" << g_sampleName
<< " picked_lumi_key=" << g_selected_lumi_key
<< " lumi=" << g_selected_lumi
<< " xsec=" << g_sectionSamples[g_era][g_sampleName].xsec
<< " sumw=" << g_sectionSamples[g_era][g_sampleName].sumw
<< " base=" << sampleBaseWeight(g_sampleName)
<< std::endl;
}

    //printCuts();

    // Start Main Loop Process 
    initZPtRW();                                      // load Zpt Root file
    loadMuonSF();                                     // load muon ID/ISO/TRG SF JSON
    if (__useRochester) loadRochesterCorrections();   // Rochester only when requested
    loadMuonScaleSmearCorrections();                  // load MuonScaleSmearing JSONs
    if (__useKIT) loadKITCorrections();               // KIT only when requested   :contentReference[oaicite:6]{index=6}
    if (__isMC) {
        loadJetCorrections();                         // JEC/JER for MC
        loadPileupCorrections();            		  // PU weights for MC
    } else {
        loadJetCorrections();               		  // if you want DATA residuals
    }
	if (!__isMC && !noGolden) {             		  // loads Golden Json files 
        loadGoldenJsonCorrections();
        std::cout << "[GoldenJSON] runs loaded = " << g_lumiMask.size() << std::endl;
    }
    process();

    std::cout << "[MuonCorr] Mode=" << mode
              << "  isMC=" << std::boolalpha << __isMC
              << "  era=" << g_era
              << "\n";

	return 0;
}
