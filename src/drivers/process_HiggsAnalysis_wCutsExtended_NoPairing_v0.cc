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
#include "correction.h"

// ROOT Header Files
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "LumiReweightingStandAlone.h"
#include "HistogramSets.h"
#include "boost/program_options.hpp"
#include "signal.h"
#include "memory"
#include <map>
#include <set>
#include "TVector2.h"
#include <cmath> 

// Reading Weights
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <regex>


//*******************************************************************
//                   HISTOGRAM MACROS
//*******************************************************************

// -------------------------------------
//             MUONS 
// -------------------------------------

#define FILL_MUON_HIST(set, weight) \
    set.hDiMuonMass->Fill(corr_p4dimuon.M(), weight); \
    set.hDiMuonpt->Fill(corr_p4dimuon.Pt(), weight); \
    set.hDiMuoneta->Fill(corr_p4dimuon.Eta(), weight); \
    set.hDiMuonrapidity->Fill(corr_p4dimuon.Rapidity(), weight); \
    set.hMuon1pt->Fill(corr_p4m1.Pt(), weight); \
    set.hMuon2pt->Fill(corr_p4m2.Pt(), weight); \
    set.hMuon1eta->Fill(corr_p4m1.Eta(), weight); \
    set.hMuon2eta->Fill(corr_p4m2.Eta(), weight); \
    set.hMuon1PTdDM->Fill(PTdDM1, weight); \
    set.hMuon2PTdDM->Fill(PTdDM2, weight); \

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

#define FILL_CORR_JET_HIST(set, weight) FILL_JET_HIST(set, weight)


//*******************************************************************
//                     CONSTANTS
//*******************************************************************

// Service Globals
std::string __inputfilename;
std::string __outputfilename;
std::string __puMCfilename;
std::string __puDATAfilename;
std::string const NTUPLEMAKER_NAME =  "";

bool __isMC;
bool __genPUMC;
std::atomic<bool> __continueRunning = true;

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

// Era handling
std::string g_era = "12022";  // default section (can be overridden via --era)
std::unordered_map<std::string, std::unordered_map<std::string, double>> g_sectionLumis;
std::unordered_map<std::string, std::unordered_map<std::string, SampleConfig>> g_sectionSamples;

double g_selected_lumi = 0.0; 


// Kinematic Cuts  
enum MuonType {kLead=0, kSubLead=1};  // Muon Id --> lead = 0, sublead =1
double _leadmuonPt = 29.;
double _subleadmuonPt = 20.;
double _leadmuonEta = 2.4;
double _subleadmuonEta = 2.4;
double _leadmuonIso = 0.15;
double _subleadmuonIso = 0.25;
double _leadJetPt = 35.; //Kind of unsure --> Purdue has 50 GeV w/out Loose PU ID
double _subleadJetPt = 25.;
double _leadJetEta = 4.7;
double _subleadJetEta = 4.7;
double _dijetMass_VBF = 400;
double _dijetdEta_VBF = 2.5;

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
std::shared_ptr<correction::CorrectionSet> corrections = nullptr;     // muon SFs
std::shared_ptr<correction::CorrectionSet> muoncorrections = nullptr; // Rochester

// ---------- Era-aware config (edit defaults if your layout changes) ----------
std::string g_muDir     = "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/MUO";
std::string g_muSFJSON  = "auto";   // "auto" => muon_Z<era>.json
std::string g_muSFKeyID = "NUM_TightID_DEN_TrackerMuons";
std::string g_muSFKeyISO = "NUM_TightPFIso_DEN_TightID";
std::string g_muSFKeyTRG = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";

std::string g_rochJSON  = "auto";   // "auto" => RoccoR<era>.json
static std::string g_rochKeyDT = "kScaleDT";
static std::string g_rochKeyMC = "kScaleMC";

static std::string muSFFileForEra(const std::string& era) {
    if (g_muSFJSON != "auto") return g_muSFJSON;
    // Example naming: muon_Z12022.json, muon_Z22022.json, muon_Z12023.json, muon_Z22023.json
    return g_muDir + "/muon_Z" + era + ".json";
}

static std::string rochFileForEra(const std::string& era) {
    if (g_rochJSON != "auto") return g_rochJSON;
    // Example naming: RoccoR12022.json, RoccoR22022.json, RoccoR12023.json, RoccoR22023.json
    return g_muDir + "/RoccoR" + era + ".json";
}

void loadMuonCorrections() {
    try {
        const std::string sfPath   = muSFFileForEra(g_era);
        const std::string rochPath = rochFileForEra(g_era);

        // Load (or reload if null). If you want forced reload per era, you can reset to nullptr first.
        if (!corrections)     corrections     = correction::CorrectionSet::from_file(sfPath);
        if (!muoncorrections) muoncorrections = correction::CorrectionSet::from_file(rochPath);

        std::cerr << "[MUON] Loaded SFs from " << sfPath << "\n";
        std::cerr << "[MUON] Loaded Rochester from " << rochPath << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[Muon load ERROR] " << e.what() << std::endl;
        corrections = nullptr;
        muoncorrections = nullptr;
    }

    // --- Verify Rochester keys (probe + fallback) ---
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
            std::cerr << "[MUON] ERROR: No valid Rochester DT key found (" << e2.what() << "). DT scale will be 1.0.\n";
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
            std::cerr << "[MUON] ERROR: No valid Rochester MC key found (" << e2.what() << "). MC scale will be 1.0.\n";
            g_rochKeyMC.clear();
        }
    }
}

// -------------------------------------
//             JETS 
// -------------------------------------

// Jet Energy Corrections/Resolutions (JEC/JER) 
std::shared_ptr<correction::CorrectionSet> jetCorrections = nullptr;
std::shared_ptr<correction::CorrectionSet> jetVetoMaps = nullptr;
std::shared_ptr<correction::CorrectionSet> jetIdCorrections = nullptr;

// ---------- Era-aware config & helpers ----------
std::string g_jmeDir        = "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/JME";
std::string g_jercJSON      = "auto"; // "auto" => jet_jerc<era>.json
std::string g_jetVetoJSON   = "auto"; // "auto" => jetvetomaps<era>.json
std::string g_jetIdJSON     = "auto"; // "auto" => jetid<era>.json

struct JecEraConfig {
  std::string prefixBase;   // base for JEC keys, we’ll append _MC or _<RunTag>_DATA
  std::string dataRunTag;   // e.g. RunCD, RunE, RunF, RunG, RunCv123, RunCv4, RunD
  std::string jerScaleKey;  // exact JER SF key
  std::string jerResKey;    // exact JER resolution key
  std::string vetoMapKey;   // veto map key (if available for this era/tag); leave "" to skip
};

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
    // 2022 BCD
    return {
      "Summer22_22Sep2023_V2",        // prefix base
      "RunCD",                        // DATA tag
      "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",   // JER SF
      "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",  // JER Res
      "Summer22_23Sep2023_RunCD_V1"   // veto map
    };
  }
  if (era == "22022") {
    // 2022 E/F/G (DATA is split by run; infer from input)
    const std::string runTag = inferDataRunTag(era, inputFile);
    return {
      "Summer22EE_22Sep2023_V2",
      runTag, // "RunE" or "RunF" or "RunG"
      "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi",
      "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi",
      "Summer22EE_23Sep2023_RunEFG_V1"
    };
  }
  if (era == "12023") {
    // 2023 C (split into RunCv123 vs RunCv4 for DATA; MC uses Prompt set)
    const std::string runTag = inferDataRunTag(era, inputFile); // RunCv123 (default) or RunCv4
    return {
      "Summer23Prompt23_V2",
      runTag,
      // JER MC keys available in your list (both combined and split exist; use combined set):
      "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi",
      "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi",
      "" // veto maps not listed explicitly for 2023 in your dumps; leave blank to skip
    };
  }
  if (era == "22023") {
    // 2023 D (BPixPrompt set)
    return {
      "Summer23BPixPrompt23_V3",
      "RunD",
      "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
      "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
      "" // veto map not listed explicitly; leave blank to skip
    };
  }
  // Fallback to 2022 BCD
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

void loadJetCorrections() {
  try {
    if (!jetCorrections)
      jetCorrections = correction::CorrectionSet::from_file(jercFileForEra(g_era));
  } catch (const std::exception& e) {
    std::cerr << "[Jet JERC load ERROR] " << e.what() << std::endl;
    jetCorrections = nullptr;
  }
  try {
    if (!jetVetoMaps)
      jetVetoMaps = correction::CorrectionSet::from_file(vetoFileForEra(g_era));
  } catch (const std::exception& e) {
    std::cerr << "[Jet veto map load ERROR] " << e.what() << std::endl;
    jetVetoMaps = nullptr;
  }
  try {
    if (!jetIdCorrections)
      jetIdCorrections = correction::CorrectionSet::from_file(jetIdFileForEra(g_era));
  } catch (const std::exception& e) {
    std::cerr << "[JetID load ERROR] " << e.what() << std::endl;
    jetIdCorrections = nullptr;
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
    pu_corrections = nullptr;
  }
}

//*******************************************************************
//                      MOMENTUM CORRECTIONS
//*******************************************************************

// -------------------------------------
//             MUONS 
// -------------------------------------

float getMuonRochester(const Muon& mu, bool isMC)
{
    if (!muoncorrections) {
        std::cerr << "[ERROR] Rochester corrections not loaded!\n";
        return 1.0f;
    }

    // --- Sanitize inputs ---
    int   charge = (mu._charge > 0) ? +1 : (mu._charge < 0 ? -1 : +1); // zero -> +1
    float pt     = std::isfinite(mu._bsConstrainedPt) ? mu._bsConstrainedPt : mu._pt;
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
//             JETS 
// -------------------------------------

float getCorrectedJetPt(const Jet& jet, float rho, bool isMC, float genJetPt, bool hasGenJet, int run) {
    if (!jetCorrections) {
        std::cerr << "[ERROR] Jet corrections not loaded!\n";
        return jet._pt;
    }

    const JecEraConfig cfg = jecConfigForEra(g_era, __inputfilename);

    float eta  = std::clamp(jet._eta, -5.191f, 5.191f);
    float pt   = std::max(jet._pt, 0.1f);
    float phi  = jet._phi;
    float area = jet._area;

    // JEC prefixes (MC vs DATA)
    const std::string prefix = isMC
        ? (cfg.prefixBase + "_MC")
        : (cfg.prefixBase + "_" + cfg.dataRunTag + "_DATA");
  
    static bool jecOnce = false;
    if (!jecOnce) {
        std::cerr << "[JEC] era=" << g_era << " runTag=" << cfg.dataRunTag << " prefix=" << prefix << std::endl;
        jecOnce = true;
    }

// --- L2/L3/Residual adapter ----------------------------------------------
// Some JEC pieces take (eta,phi,pt), others only (eta,pt).
auto evalEtaPhiPtOrEtaPt = [&](const std::string& key, float eta_, float phi_, float pt_) -> float {
    static std::map<std::string,int> mode;
    try {
        if (mode[key] != 2) {
            double v3 = jetCorrections->at(key)->evaluate({eta_, phi_, pt_});
            mode[key] = 3;
            static std::set<std::string> printed;
            if (!printed.count(key)) {
                std::cerr << "[JEC] " << key << " uses 3-arg schema (eta,phi,pt)\n";
                printed.insert(key);
            }
                return v3;
        }
    } catch (...) { /* fall through */ }
    
    double v2 = jetCorrections->at(key)->evaluate({eta_, pt_});
    mode[key] = 2;
    static std::set<std::string> printed;
    if (!printed.count(key)) {
        std::cerr << "[JEC] " << key << " uses 2-arg schema (eta,pt)\n";
        printed.insert(key);
    }
        return v2;
  };

// Insert the run tag (e.g. "RunE") just before the trailing "_V..." in the base,
// e.g. base="Summer22EE_22Sep2023_V2" + "RunE" -> "Summer22EE_22Sep2023_RunE_V2"
auto insertRunBeforeVersion = [](const std::string& base, const std::string& runTag) {
    const auto pos = base.rfind("_V");
    if (pos == std::string::npos) return base + "_" + runTag; // fallback
        return base.substr(0, pos) + "_" + runTag + base.substr(pos);
};

// Build candidate DATA keys for a given step ("L1FastJet_AK4PFPuppi", "L2Relative_AK4PFPuppi", ...)
auto dataKeyCandidates = [&](const std::string& step) -> std::vector<std::string> {
    std::vector<std::string> ks;

    if (g_era == "12022") {
        // Exact from your file: Summer22_22Sep2023_RunCD_V2_DATA_<STEP>
        ks.push_back(insertRunBeforeVersion(cfg.prefixBase, "RunCD") + "_DATA_" + step);
        // Loose fallback: no run tag
        ks.push_back(cfg.prefixBase + "_DATA_" + step);
    }
    else if (g_era == "22022") {
        // Exact per-run keys exist in your file for E/F/G
        ks.push_back(insertRunBeforeVersion(cfg.prefixBase, cfg.dataRunTag) + "_DATA_" + step); // RunE/RunF/RunG
        // Loose fallback: no run tag
        ks.push_back(cfg.prefixBase + "_DATA_" + step);
    }
    else if (g_era == "12023") {
        // Steps have NO run tag in keys: Summer23Prompt23_V2_DATA_<STEP>
        ks.push_back(cfg.prefixBase + "_DATA_" + step);
    }
    else if (g_era == "22023") {
        // Steps have NO run tag in keys: Summer23BPixPrompt23_V3_DATA_<STEP>
        ks.push_back(cfg.prefixBase + "_DATA_" + step);
    } else {
        // Generic fallback
        ks.push_back(cfg.prefixBase + "_DATA_" + step);
    }

    return ks;
};

// --- DATA compound L1L2L3Res adapter ------------------------------------
// Some DATA JSONs expose a compound correction key:
//   <prefixBase>_DATA_L1L2L3Res_AK4PFPuppi
// Inputs can be (area,eta,pt,rho,phi) OR (area,eta,pt,rho,phi,run).
// Try (A,eta,pt,rho,run) first (12023), then other common forms.
// Prefer known signatures per era, then fall back broadly.
auto evalDataCompoundL1L2L3Res = [&](float area_, float eta_, float pt_, float rho_, float phi_, int run_)
    -> std::optional<float>
{
  std::vector<std::string> keys;
  auto addBothSpellings = [&](const std::string& base) {
    keys.push_back(base + "_DATA_L1L2L3Res_AK4PFPuppi");
    keys.push_back(base + "_DATA_L1L2L3Residual_AK4PFPuppi");
  };

  if (g_era == "12022")      addBothSpellings(insertRunBeforeVersion(cfg.prefixBase, "RunCD"));
  else if (g_era == "22022") addBothSpellings(insertRunBeforeVersion(cfg.prefixBase, cfg.dataRunTag));
  else                       addBothSpellings(cfg.prefixBase); // 2023: no run tag in key name

  auto printOnce = [](const std::string& msg) {
    static std::set<std::string> seen;
    if (!seen.count(msg)) { std::cerr << msg << "\n"; seen.insert(msg); }
  };

  for (const auto& key : keys) {
    // --- Era-specific priorities ---
    if (g_era == "22023") {
      // 6-arg (A,eta,pt,rho,phi,run)
      try {
        double v = jetCorrections->at(key)->evaluate({area_, eta_, pt_, rho_, phi_, (double)run_});
        printOnce("[JEC] Using DATA compound key (A,eta,pt,rho,phi,run): " + key);
        return (float)v;
      } catch (...) {}
    }
    if (g_era == "12023") {
      // 5-arg (A,eta,pt,rho,run)
      try {
        double v = jetCorrections->at(key)->evaluate({area_, eta_, pt_, rho_, (double)run_});
        printOnce("[JEC] Using DATA compound key (A,eta,pt,rho,run): " + key);
        return (float)v;
      } catch (...) {}
    }

    // --- Generic fallbacks (cover any payload flavor) ---
    // 6-arg (A,eta,pt,rho,phi,run)
    try {
      double v = jetCorrections->at(key)->evaluate({area_, eta_, pt_, rho_, phi_, (double)run_});
      printOnce("[JEC] Using DATA compound key (A,eta,pt,rho,phi,run): " + key);
      return (float)v;
    } catch (...) {}

    // 5-arg (A,eta,pt,rho,run)
    try {
      double v = jetCorrections->at(key)->evaluate({area_, eta_, pt_, rho_, (double)run_});
      printOnce("[JEC] Using DATA compound key (A,eta,pt,rho,run): " + key);
      return (float)v;
    } catch (...) {}

    // 5-arg (A,eta,pt,rho,phi)
    try {
      double v = jetCorrections->at(key)->evaluate({area_, eta_, pt_, rho_, phi_});
      printOnce("[JEC] Using DATA compound key (A,eta,pt,rho,phi): " + key);
      return (float)v;
    } catch (...) {}

    // 4-arg (A,eta,pt,rho)
    try {
      double v = jetCorrections->at(key)->evaluate({area_, eta_, pt_, rho_});
      printOnce("[JEC] Using DATA compound key (A,eta,pt,rho): " + key);
      return (float)v;
    } catch (...) {}
  }

  return std::nullopt;
};

// 2022 residuals: usually (eta,pt) or (eta,phi,pt)
// 2023 residuals: (run,eta,pt)
auto evalDataResidual = [&](float eta_, float phi_, float pt_, int run_) -> float {
  const std::vector<std::string> steps = {
    "L2L3Residual_AK4PFPuppi",
    "L2L3Res_AK4PFPuppi" // alias, just in case
  };

  auto printOnce = [](const std::string& msg) {
    static std::set<std::string> seen;
    if (!seen.count(msg)) { std::cerr << msg << "\n"; seen.insert(msg); }
  };

  // Helper that tries a list of step names with a given evaluator lambda
  auto tryWith = [&](auto&& evalOne) -> std::optional<float> {
    for (const auto& step : steps) {
      for (const auto& k : dataKeyCandidates(step)) {
        try {
          float v = (float)evalOne(k);
          printOnce("[JEC] Residual matched key: " + k);
          return v;
        } catch (...) {
          // try next candidate
        }
      }
    }
    return std::nullopt;
  };

  if (g_era == "12022" || g_era == "22022") {
    // Prefer (eta,pt)
    if (auto v = tryWith([&](const std::string& k){ return jetCorrections->at(k)->evaluate({eta_, pt_}); }); v) return *v;
    // Fallback (eta,phi,pt)
    if (auto v = tryWith([&](const std::string& k){ return jetCorrections->at(k)->evaluate({eta_, phi_, pt_}); }); v) return *v;
    // Last resort (run,eta,pt) — rarely needed in 2022
    if (auto v = tryWith([&](const std::string& k){ return jetCorrections->at(k)->evaluate({(double)run_, eta_, pt_}); }); v) return *v;

    static bool warned22 = false;
    if (!warned22) { std::cerr << "[JEC WARN] No 2022 DATA residual matched (tried (eta,pt),(eta,phi,pt),(run,eta,pt)). Using 1.0.\n"; warned22 = true; }
    return 1.0f;
  } else {
    // 2023+: prefer (run,eta,pt)
    if (auto v = tryWith([&](const std::string& k){ return jetCorrections->at(k)->evaluate({(double)run_, eta_, pt_}); }); v) return *v;
    // Fallbacks
    if (auto v = tryWith([&](const std::string& k){ return jetCorrections->at(k)->evaluate({eta_, pt_}); }); v) return *v;
    if (auto v = tryWith([&](const std::string& k){ return jetCorrections->at(k)->evaluate({eta_, phi_, pt_}); }); v) return *v;

    static bool warned23 = false;
    if (!warned23) { std::cerr << "[JEC WARN] No 2023 DATA residual matched (tried (run,eta,pt),(eta,pt),(eta,phi,pt)). Using 1.0.\n"; warned23 = true; }
    return 1.0f;
  }
};

// Evaluate a DATA JEC step with fallbacks (handles (eta,phi,pt) vs (eta,pt) internally via your evalEtaPhiPtOrEtaPt)
auto evalDataStep = [&](const std::string& step, float eta_, float phi_, float pt_) -> float {
    for (const auto& k : dataKeyCandidates(step)) {
        try {
            return evalEtaPhiPtOrEtaPt(k, eta_, phi_, pt_);
        } catch (const std::exception&) {
        // try next candidate
        }
    }

static std::set<std::string> warned;
    if (!warned.count(step)) {
        std::cerr << "[JEC WARN] No DATA key matched for step '" << step << "' (era=" << g_era << ", tag=" << cfg.dataRunTag << "). Using 1.0.\n";
        warned.insert(step);
    }
    return 1.0f;
};

// Evaluate DATA L1 with fallbacks (L1 in all your sets is 4-arg: area,eta,pt,rho)
auto evalDataL1 = [&](float area_, float eta_, float pt_, float rho_) -> float {
    const std::string step = "L1FastJet_AK4PFPuppi";
    for (const auto& k : dataKeyCandidates(step)) {
        try {
            return (float)jetCorrections->at(k)->evaluate({area_, eta_, pt_, rho_});
        } catch (const std::exception&) {
        // try next candidate
        }
    }

static bool warned = false;
    if (!warned) {
        std::cerr << "[JEC WARN] No DATA key matched for L1FastJet (era=" << g_era << ", tag=" << cfg.dataRunTag << "). Using 1.0.\n";
        warned = true;
    }
    return 1.0f;
};

float jecFactor = 1.0f;

if (!isMC) {
  if (auto comp = evalDataCompoundL1L2L3Res(area, eta, pt, rho, phi, run); comp.has_value()) {
    jecFactor = comp.value();
  } else {
    float l1  = evalDataL1(area, eta, pt, rho);                                // (A,eta,pt,rho)
    float l2  = evalDataStep("L2Relative_AK4PFPuppi",  eta, phi, pt);          // (eta,phi,pt) or (eta,pt)
    float l3  = evalDataStep("L3Absolute_AK4PFPuppi",  eta, phi, pt);          // (eta,phi,pt) or (eta,pt)
    float res = evalDataResidual(eta, phi, pt, run);  // 2022 → (eta,pt), 2023 → (run,eta,pt)
    jecFactor = l1 * l2 * l3 * res;
  }
} else {
  // MC: keep your current MC path
  float l1  = (float)jetCorrections->at(prefix + "_L1FastJet_AK4PFPuppi")->evaluate({area, eta, pt, rho});
  float l2  = evalEtaPhiPtOrEtaPt(prefix + "_L2Relative_AK4PFPuppi",  eta, phi, pt);
  float l3  = evalEtaPhiPtOrEtaPt(prefix + "_L3Absolute_AK4PFPuppi",  eta, phi, pt);
  float res = 1.0f; // residuals are DATA-only
  jecFactor = l1 * l2 * l3 * res;
}

float corr_pt = pt * jecFactor;

// --- JER adapters -------------------------------------------------------
// Some JER ScaleFactor JSONs take (eta, pt, "nom"), others take (eta, "nom") only.
auto evalJerScale = [&](const std::string& key, double eta_, double pt_) -> double {
    static int mode = 0; // 0=unknown, 3=(eta,pt,"nom"), 2=(eta,"nom")
    try {
        if (mode != 2) {
            double v3 = jetCorrections->at(key)->evaluate({eta_, pt_, std::string("nom")});
            mode = 3;
            static bool printed3=false;
            if (!printed3) { std::cerr << "[JER] " << key << " uses 3-arg schema (eta,pt,\"nom\")\n"; printed3=true; }
                return v3;
        }
    } catch (...) { /* fall through to 2-arg */ }
    
    double v2 = jetCorrections->at(key)->evaluate({eta_, std::string("nom")});
    mode = 2;
    static bool printed2=false;
    if (!printed2) { std::cerr << "[JER] " << key << " uses 2-arg schema (eta,\"nom\")\n"; printed2=true; }
        return v2;
    };

// Some JER Resolution JSONs take (eta, pt, rho), others take (eta, pt).
auto evalJerRes = [&](const std::string& key, double eta_, double pt_, double rho_) -> double {
    static int mode = 0; // 0=unknown, 3=(eta,pt,rho), 2=(eta,pt)
    try {
        if (mode != 2) {
            double v3 = jetCorrections->at(key)->evaluate({eta_, pt_, rho_});
            mode = 3;
            static bool printed3=false;
            if (!printed3) { std::cerr << "[JER] " << key << " uses 3-arg schema (eta,pt,rho)\n"; printed3=true; }
                return v3;
        }
    } catch (...) { /* fall through to 2-arg */ }
    
    double v2 = jetCorrections->at(key)->evaluate({eta_, pt_});
    mode = 2;
    static bool printed2=false;
    if (!printed2) { std::cerr << "[JER] " << key << " uses 2-arg schema (eta,pt)\n"; printed2=true; }
        return v2;
 };

  // ---- JER (MC only) ----
  if (isMC) {
    // Use JRV1 keys from your JSON
    double sf = evalJerScale(cfg.jerScaleKey, (double)eta, (double)corr_pt);

    if (hasGenJet && genJetPt > 0.f) {
      // Truth-anchored smearing
      float smear = 1.0f + (sf - 1.0f) * (corr_pt - genJetPt) / std::max(corr_pt, 1e-3f);
      // Optional safety: avoid negative scaling
      smear = std::max(smear, 0.0f);
      corr_pt *= smear;
    } else {
      // Stochastic smearing
      double reso = evalJerRes(cfg.jerResKey, (double)eta, (double)corr_pt, (double)rho);
      float sigma = static_cast<float>(reso * std::sqrt(std::max(sf*sf - 1.0, 0.0)));
      // Consider a deterministic seed (event, jet index) if you run multithreaded
      static TRandom3 rand(42);
      corr_pt *= (1.0f + rand.Gaus(0, sigma));
      corr_pt  = std::max(corr_pt, 0.0f);
    }
  }

  return corr_pt;
}

inline bool isVetoedJet(float eta, float phi) {
  if (!jetVetoMaps) { std::cerr << "[ERROR] Jet veto corr not loaded!\n"; return false; }
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

float jetMuondR(float jeta,float jphi, float meta, float mphi)
{
	TLorentzVector p4j,p4m;
	p4j.SetPtEtaPhiM(10, jeta, jphi, 0);
	p4m.SetPtEtaPhiM(10, meta, mphi, 0);
	return p4j.DeltaR(p4m);
}

//*******************************************************************
//                      BOOLEAN CONDITIONS
//*******************************************************************

// -------------------------------------
//              MUON LOCATIONS
// -------------------------------------

bool isBarrel(Muon const& m)
{
	return TMath::Abs(m._eta)<0.9;
}

bool isOverlap(Muon const& m)
{
	return TMath::Abs(m._eta)>=0.9 && TMath::Abs(m._eta)<1.6;
}

bool isEndcap(Muon const& m)
{
	return TMath::Abs(m._eta)>=1.6 && TMath::Abs(m._eta)<2.4;
}

// -------------------------------------
//LEPTON CRITERIA 
// -------------------------------------

bool vetoAdditionalLeptons(Muons* muons, Electrons* electrons, bool isMC)
{
    int numPassingMuons = 0;

    // Check muons
    for (const auto& mu : *muons) {
        float rochSF = getMuonRochester(mu, isMC);
        float corrected_pt = mu._pt * rochSF;

        if (corrected_pt > 20.0 && std::abs(mu._eta) < 2.4) {
            if (mu._isGlobal && mu._isTracker && mu._isTight && mu._pfRelIso04_all <= _leadmuonIso) {
                numPassingMuons++;
                if (numPassingMuons > 2) return false; // Veto: more than 2 passing muons
            }
        }
    }

    // Check electrons
    for (const auto& el : *electrons) {
        if (el._pt > 20.0 && std::abs(el._eta) < 2.5 && el._cutBased >= 1) {
            return false; // Veto: found a passing electron
        }
    }

    return true; // Passed: ≤2 good muons AND no good electrons
}

// -------------------------------------
//VERTEX CRITERIA 
// -------------------------------------

bool passVertex(Vertices* v)
{
	if (v->size()==0)
		return false;
	
    int n=0;
	for (Vertices::const_iterator it=v->begin(); it!=v->end() && n<20; ++it)
	{
        double rho = std::sqrt((it->_x * it->_x) + (it->_y * it->_y));

        if (TMath::Abs(it->_z) < 24 
            && it->_ndof > 4 
            && rho <= 0.3)
			return true;
        n++;
	}
	return false;
}

//WILL PROB IMPLEMENT THIS ONE
inline bool passVertex(Vertices* /*v*/, int PV_npvsGood) {
    return PV_npvsGood >= 1;
}

// -------------------------------------
//MUON CRITERIA 
// -------------------------------------

bool passMuonHLT(TrigObjs* trig, Muon const& m, bool isMC)
{
    // If we don't have trigger objects on DATA, fall back to a simple pT check.
    if (trig->size()==0) {
        if (!isMC) {
            float rochSF = getMuonRochester(m, /*isMC=*/false);
            float corrected_pt = m._bsConstrainedPt * rochSF;
            return (corrected_pt >= 26.0f);
        }
         return false;
    }
    
    // Apply Rochester correction
    float rochSF = getMuonRochester(m, isMC);
    float corrected_pt = m._bsConstrainedPt * rochSF;

    if (corrected_pt < 26)
        return false;
	
	for (TrigObjs::const_iterator it=trig->begin(); it!=trig->end(); ++it)
	{
        if (it->_id == 13 && it->_filterBits & (1 << 0)) // Bit 0 = TrkIsoVVL (HLT_IsoMu24) 
        { 
            float deta = m._eta - it->_eta;
            float dphi = std::fabs(m._phi - it->_phi);
            if (dphi > M_PI) dphi = 2 * M_PI - dphi;

            float dR = std::sqrt(deta * deta + dphi * dphi);
            if (dR < 0.3)
                return true;
            }
        }   
	return false;
}

bool passMuon(Muon const& m, float corrected_pt, MuonType id=kLead)
{
    if (id==kLead)
    {
	    if (m._isGlobal 
            && m._isTracker 
            && corrected_pt > _leadmuonPt 
            && TMath::Abs(m._eta) < _leadmuonEta 
            && m._isTight 
            && m._pfRelIso04_all <= _leadmuonIso)
		    return true;
    }
    else
    {
	    if (m._isGlobal 
            && m._isTracker 
            && corrected_pt > _subleadmuonPt 
            && TMath::Abs(m._eta) < _subleadmuonEta 
            && m._isTight 
            && m._pfRelIso04_all <= _subleadmuonIso)
		    return true;
    }
	return false;
}

bool passMuons(Muon const& m1, Muon const& m2, TrigObjs* trigObjs)
{
    // Get Rochester scale factors for both muons
    float rochSF1 = getMuonRochester(m1, __isMC);
    float rochSF2 = getMuonRochester(m2, __isMC);

    float pt1_corr = m1._bsConstrainedPt * rochSF1;
    float pt2_corr = m2._bsConstrainedPt * rochSF2;
    
    const Muon& leadmuon     = (pt1_corr > pt2_corr) ? m1 : m2;
    const Muon& subleadmuon  = (pt1_corr > pt2_corr) ? m2 : m1;
    float lead_pt_corr       = std::max(pt1_corr, pt2_corr);
    float sublead_pt_corr    = std::min(pt1_corr, pt2_corr);
	
    if (m1._charge!=m2._charge && passMuon(leadmuon, lead_pt_corr, kLead) && passMuon(subleadmuon, sublead_pt_corr, kSubLead))
        if (passMuonHLT(trigObjs, m1, __isMC) || passMuonHLT(trigObjs, m2, __isMC))
            return true;

	return false;
}

// -------------------------------------
//JET CRITERIA 
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
        // NanoAOD: 1=loose, 2=tight, 4=tightLepVeto. Keep your previous behavior:
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
    /*if (!isMC) {
        std::cerr << "This is Data!\n";
        return 1.0;
    }*/

    if (!corrections) {
        std::cerr << "[ERROR] Correction set not initialized!\n";
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

//*******************************************************************
//                         CATEGORIZE FUNCTION 
//*******************************************************************

void categorize(analysis::core::Jets* jets, analysis::core::GenJets* genjets, const analysis::core::Muon& mu1, 
const analysis::core::Muon& mu2, analysis::core::FsrPhotons* fsrphotons, analysis::core::Photons* photons, 
analysis::core::Vertices* vertices, const analysis::core::MET& met, const analysis::core::Event& event, 
analysis::core::Rhos* rhos, HistSetCollection& hists, float Dataweight) {

    // -------------------------------------
    //             MUONS 
    // -------------------------------------

    // TLorentz Vectors for Muons 1 & 2
    float roch_sf1 = getMuonRochester(mu1, __isMC);
    float roch_sf2 = getMuonRochester(mu2, __isMC);

    TLorentzVector p4m1, p4m2;
    p4m1.SetPtEtaPhiM(mu1._pt, mu1._eta, mu1._phi, PDG_MASS_Mu);
	p4m2.SetPtEtaPhiM(mu2._pt, mu2._eta, mu2._phi, PDG_MASS_Mu);
    
    TLorentzVector corr_p4m1 = p4m1 * roch_sf1;
    TLorentzVector corr_p4m2 = p4m2 * roch_sf2;
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
    applyFSR(mu1, p4m1, corr_p4m1, fsrphotons, 0);
    applyFSR(mu2, p4m2, corr_p4m2, fsrphotons, 1);

    // With applied FSR corrections 
    TLorentzVector p4dimuon = p4m1 + p4m2;    
    TLorentzVector corr_p4dimuon = corr_p4m1 + corr_p4m2;    
	double dphi = p4m1.DeltaPhi(p4m2);
   
    // -------------------------------------
    //              JETS  
    // -------------------------------------
    
    bool pass_bjet_veto = true;
    int n_loose_bjets = 0;
    int n_medium_bjets = 0;
    float rho = rhos->at(0)._fixedGridRhoFastjetAll;
    
    std::vector<TLorentzVector> p4jets;
    p4jets.reserve(jets->size());

    for (const auto& j : *jets) {
        // veto map
        if (isVetoedJet(j._eta, j._phi)) continue;

        // gen matching for JER
        float genJetPt = 0.f; 
        bool hasGenJet = false;
        
        if (__isMC && j._genJetIdx >= 0 && j._genJetIdx < (int)genjets->size()) {
            const auto& gj = genjets->at(j._genJetIdx);
            genJetPt = gj._pt; 
            hasGenJet = true;
        }

        int runNum = __isMC ? 0 : static_cast<int>(event._run);
        float corr_pt = getCorrectedJetPt(j, rho, __isMC, genJetPt, hasGenJet, runNum);
        //float corr_pt   = getCorrectedJetPt(j, rho, __isMC, genJetPt, hasGenJet);
        float corr_mass = j._mass * (corr_pt / std::max(1e-6f, j._pt));

        if (corr_pt <= 25.f) continue;
        if (std::abs(j._eta) >= 4.7f) continue;
        if ( (j._jetId & 4) != 4 ) continue; // Tight ID
        //const bool passJetID = passJetIDfromJSON(j, JetAlgo::AK4PUPPI, JetIDWP::Tight);
        //if (!passJetID) continue;

        // muon cleaning
        if (jetMuondR(j._eta, j._phi, mu1._eta, mu1._phi) < 0.4) continue;
        if (jetMuondR(j._eta, j._phi, mu2._eta, mu2._phi) < 0.4) continue;

        TLorentzVector p4; p4.SetPtEtaPhiM(corr_pt, j._eta, j._phi, corr_mass);
        p4jets.push_back(p4);

        // b-tag counting should use the same cleaned/accepted jets as above
        if (j._btagDeepFlavB > 0.0583) n_loose_bjets++;
        if (j._btagDeepFlavB > 0.3086) n_medium_bjets++;
    }

    // sort by corrected pT
    std::sort(p4jets.begin(), p4jets.end(), [](const TLorentzVector& a, const TLorentzVector& b){ return a.Pt() > b.Pt(); });

    if (p4jets.size() < 2) return; // no dijet
    
    //  Choosing Largest Dijet Mass Pair
    double best_mjj = -1.0;
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

    // safety (shouldn't happen if size>=2)
    if (i_best < 0 || j_best < 0) return;

    TLorentzVector jA = p4jets[i_best];
    TLorentzVector jB = p4jets[j_best];

    // order the chosen pair by pT (so "lead/sub" still mean pT-ordered)
    TLorentzVector p4lead = (jA.Pt() >= jB.Pt()) ? jA : jB;
    TLorentzVector p4sub  = (jA.Pt()  <  jB.Pt()) ? jA : jB;

    TLorentzVector dijet = p4lead + p4sub;
    
    float dijetdeta      = TMath::Abs(p4lead.Eta() - p4sub.Eta());
    float jet1eta        = p4lead.Eta();
    float jet2eta        = p4sub.Eta();
    
    float dijetmass      = dijet.M();
    
    float dijetpt        = dijet.Pt();
    float jet1pt         = p4lead.Pt();
    float jet2pt         = p4sub.Pt();
   
    float dphi_signed    = p4lead.DeltaPhi(p4sub);   
    float dijetdphi      = std::fabs(dphi_signed); 
    
    float PTdDM1         = (corr_p4m1.Pt())/(dijet.M());
    float PTdDM2         = (corr_p4m2.Pt())/(dijet.M());

    // -------------------------------------
    //          SF && WEIGHTS  
    // -------------------------------------
   
    float mu1_SF         = getMuonSF(mu1, __isMC);
    float mu2_SF         = getMuonSF(mu2, __isMC);
    float totalWeight    = Dataweight * mu1_SF * mu2_SF;
    float noWeight       = 1.0;
    
    // -------------------------------------
    //        PRE-SELECTION HISTOGRAMS  
    // -------------------------------------
   
    //RIGHT NOW THIS INCLUDES THE 1 VBF SELECTION: EXLUDING ADDITINOAL MUONS/ELECTRONS 
    FILL_MUON_HIST(hists.setDM_PreSelection, noWeight);   
    FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection, totalWeight); 
    FILL_JET_HIST(hists.setJet_PreSelection, noWeight); 
    FILL_CORR_JET_HIST(hists.corr_setJet_PreSelection, totalWeight); 

    /*if (isBarrel(mu1) && isBarrel(mu2))
    {
        //	BB
        FILL_MUON_HIST(hists.setDM_PreSelection_BB, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_BB, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_BB, totalWeight); 
    }
    else if ((isBarrel(mu1) && isOverlap(mu2)) ||
        (isBarrel(mu2) && isOverlap(mu1)))
    {
        //	BO
        FILL_MUON_HIST(hists.setDM_PreSelection_BO, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_BO, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_BO, totalWeight); 
    }
    else if ((isBarrel(mu1) & isEndcap(mu2)) || 
        (isBarrel(mu2) && isEndcap(mu1)))
    {
        //	BE
        FILL_MUON_HIST(hists.setDM_PreSelection_BE, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_BE, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_BE, totalWeight); 
    }
    else if (isOverlap(mu1) && isOverlap(mu2))
    {
        //	OO
        FILL_MUON_HIST(hists.setDM_PreSelection_OO, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_OO, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_OO, totalWeight); 
    }
    else if ((isOverlap(mu1) && isEndcap(mu2)) ||
        (isOverlap(mu2) && isEndcap(mu1)))
    {
        //	OE
        FILL_MUON_HIST(hists.setDM_PreSelection_OE, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_OE, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_OE, totalWeight); 
    }
    else if (isEndcap(mu1) && isEndcap(mu2))
    {
        //	EE
        FILL_MUON_HIST(hists.setDM_PreSelection_EE, noWeight);
        FILL_CORR_MUON_HIST(hists.corr_setDM_PreSelection_EE, totalWeight);
        FILL_JET_HIST(hists.setJet_PreSelection_EE, totalWeight); 
    }*/

    if (p4lead.Pt() < _leadJetPt) return;      // leading threshold
    if (p4sub.Pt() < _subleadJetPt) return;   // subleading threshold
    
    // -------------------------------------
    //          VBF HISTOGRAMS  
    // -------------------------------------

    if (n_loose_bjets >= 2 || n_medium_bjets >= 1)
        pass_bjet_veto = false;
    
    if (pass_bjet_veto) {

        //  PRE-JET CUTS 
        FILL_JET_HIST(hists.setJet_PreDiJet_VBF, noWeight); 
        FILL_CORR_JET_HIST(hists.corr_setJet_PreDiJet_VBF, totalWeight); 

        if (dijetmass > _dijetMass_VBF && TMath::Abs(dijetdeta) > _dijetdEta_VBF) {
            
            // -------------------------------------
            //         BDT INPUT VARIABLES  
            // -------------------------------------
            
            const double m_mumu = corr_p4dimuon.M();

            // muon-level vars
            g_m_mumu        = static_cast<float>(m_mumu);
            g_mu1_eta       = static_cast<float>(corr_p4m1.Eta());
            g_mu2_eta       = static_cast<float>(corr_p4m2.Eta());
            g_mu1_pt_over_m = (m_mumu > 0.0) ? static_cast<float>(corr_p4m1.Pt() / m_mumu) : 0.f;
            g_mu2_pt_over_m = (m_mumu > 0.0) ? static_cast<float>(corr_p4m2.Pt() / m_mumu) : 0.f;
            g_y_mumu        = static_cast<float>(corr_p4dimuon.Rapidity());

            // jet-level vars
            g_jet1_eta      = static_cast<float>(p4lead.Eta());
            g_deta_jj       = std::fabs(p4lead.Eta() - p4sub.Eta());
            g_dphi_jj       = std::fabs(TVector2::Phi_mpi_pi(p4lead.Phi() - p4sub.Phi()));
            g_y_j1          = static_cast<float>(p4lead.Rapidity());
            g_y_j2          = static_cast<float>(p4sub.Rapidity());

            // Zeppenfeld variable
            g_z             = std::fabs(g_y_mumu - 0.5f*(g_y_j1 + g_y_j2));

            // pT-balance ratio
            const double px_sum     =  corr_p4dimuon.Px() + p4lead.Px() + p4sub.Px();
            const double py_sum     =  corr_p4dimuon.Py() + p4lead.Py() + p4sub.Py();
            const double vec_sum    = std::sqrt(px_sum*px_sum + py_sum*py_sum);

            g_pt_mumu       = static_cast<float>(corr_p4dimuon.Pt());
            g_pt_j1         = static_cast<float>(p4lead.Pt());
            g_pt_j2         = static_cast<float>(p4sub.Pt());

            const double scalar_sum = static_cast<double>(g_pt_mumu + g_pt_j1 + g_pt_j2);
            g_ptBalance     = (scalar_sum > 0.0) ? static_cast<float>(vec_sum / scalar_sum) : 0.f;

            // weight (total per-event weight)
            g_weight = totalWeight;

            // fill exactly once per selected dimuon+2jet pair
            if (gMvaTree) gMvaTree->Fill();

            //	VBF SELECTION CRITERIA
            FILL_MUON_HIST(hists.setDM_VBF, noWeight)
            FILL_CORR_MUON_HIST(hists.corr_setDM_VBF, totalWeight)
            FILL_JET_HIST(hists.setJet_VBF, noWeight); 
            FILL_CORR_JET_HIST(hists.corr_setJet_VBF, totalWeight); 

            /*if (isBarrel(mu1) && isBarrel(mu2))
            {
                //	BB
                FILL_MUON_HIST(hists.setDM_VBF_BB, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_BB, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_BB, totalWeight); 
            }
            else if ((isBarrel(mu1) && isOverlap(mu2)) ||
                (isBarrel(mu2) && isOverlap(mu1)))
            {
                //	BO
                FILL_MUON_HIST(hists.setDM_VBF_BO, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_BO, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_BO, totalWeight); 
            }
            else if ((isBarrel(mu1) & isEndcap(mu2)) || 
                (isBarrel(mu2) && isEndcap(mu1)))
            {
                //	BE
                FILL_MUON_HIST(hists.setDM_VBF_BE, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_BE, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_BE, totalWeight); 
            }
            else if (isOverlap(mu1) && isOverlap(mu2))
            {
                //	OO
                FILL_MUON_HIST(hists.setDM_VBF_OO, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_OO, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_OO, totalWeight); 
            }
            else if ((isOverlap(mu1) && isEndcap(mu2)) ||
                (isOverlap(mu2) && isEndcap(mu1)))
            {
                //	OE
                FILL_MUON_HIST(hists.setDM_VBF_OE, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_OE, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_OE, totalWeight); 
            }
            else if (isEndcap(mu1) && isEndcap(mu2))
            {
                //	EE
                FILL_MUON_HIST(hists.setDM_VBF_EE, noWeight);
                FILL_CORR_MUON_HIST(hists.corr_setDM_VBF_EE, totalWeight)
                FILL_JET_HIST(hists.setJet_VBF_EE, totalWeight); 
            }*/
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
//                        MAIN PROGRESS FUCNTION 
//*******************************************************************

void process()
{
    //	Output Root File
    TFile *outroot = new TFile(__outputfilename.c_str(), "recreate");

    // Pre-Selection blocks
    MuonSet     setDM_PreSelection(outroot, "DM_Pre-Selection");
    MuonSet     corr_setDM_PreSelection(outroot, "corr_DM_Pre-Selection");
    JetSet      setJet_PreSelection(outroot, "Jet_Pre-Selection");
    JetSet      corr_setJet_PreSelection(outroot, "corr_Jet_Pre-Selection");
    /*MuonSet setDM_PreSelection_BB(outroot, "DM_Pre-Selection_BB");
    MuonSet corr_setDM_PreSelection_BB(outroot, "corr_DM_Pre-Selection_BB");
    JetSet setJet_PreSelection_BB(outroot, "Jet_Pre-Selection_BB");
    MuonSet setDM_PreSelection_BO(outroot, "DM_Pre-Selection_BO");
    MuonSet corr_setDM_PreSelection_BO(outroot, "corr_DM_Pre-Selection_BO");
    JetSet setJet_PreSelection_BO(outroot, "Jet_Pre-Selection_BO");
    MuonSet setDM_PreSelection_BE(outroot, "DM_Pre-Selection_BE");
    MuonSet corr_setDM_PreSelection_BE(outroot, "corr_DM_Pre-Selection_BE");
    JetSet setJet_PreSelection_BE(outroot, "Jet_Pre-Selection_BE");
    MuonSet setDM_PreSelection_OO(outroot, "DM_Pre-Selection_OO");
    MuonSet corr_setDM_PreSelection_OO(outroot, "corr_DM_Pre-Selection_OO");
    JetSet setJet_PreSelection_OO(outroot, "Jet_Pre-Selection_OO");
    MuonSet setDM_PreSelection_OE(outroot, "DM_Pre-Selection_OE");
    MuonSet corr_setDM_PreSelection_OE(outroot, "corr_DM_Pre-Selection_OE");
    JetSet setJet_PreSelection_OE(outroot, "Jet_Pre-Selection_OE");
    MuonSet setDM_PreSelection_EE(outroot, "DM_Pre-Selection_EE");
    MuonSet corr_setDM_PreSelection_EE(outroot, "corr_DM_Pre-Selection_EE");
    JetSet setJet_PreSelection_EE(outroot, "Jet_Pre-Selection_EE");*/
    
    // Pre dijet cuts 
    JetSet setJet_PreDiJet_VBF(outroot, "Jet_PreDiJet_VBF");
    JetSet corr_setJet_PreDiJet_VBF(outroot, "corr_Jet_PreDiJet_VBF");
    
    // VBF blocks
    MuonSet     setDM_VBF(outroot, "DM_VBF");
    MuonSet     corr_setDM_VBF(outroot, "corr_DM_VBF");
    JetSet      setJet_VBF(outroot, "Jet_VBF");
    JetSet      corr_setJet_VBF(outroot, "corr_Jet_VBF");
    /*MuonSet setDM_VBF_BB(outroot, "DM_VBF_BB");
    MuonSet corr_setDM_VBF_BB(outroot, "corr_DM_VBF_BB");
    JetSet setJet_VBF_BB(outroot, "Jet_VBF_BB");
    MuonSet setDM_VBF_BO(outroot, "DM_VBF_BO");
    MuonSet corr_setDM_VBF_BO(outroot, "corr_DM_VBF_BO");
    JetSet setJet_VBF_BO(outroot, "Jet_VBF_BO");
    MuonSet setDM_VBF_BE(outroot, "DM_VBF_BE");
    MuonSet corr_setDM_VBF_BE(outroot, "corr_DM_VBF_BE");
    JetSet setJet_VBF_BE(outroot, "Jet_VBF_BE");
    MuonSet setDM_VBF_OO(outroot, "DM_VBF_OO");
    MuonSet corr_setDM_VBF_OO(outroot, "corr_DM_VBF_OO");
    JetSet setJet_VBF_OO(outroot, "Jet_VBF_OO");
    MuonSet setDM_VBF_OE(outroot, "DM_VBF_OE");
    MuonSet corr_setDM_VBF_OE(outroot, "corr_DM_VBF_OE");
    JetSet setJet_VBF_OE(outroot, "Jet_VBF_OE");
    MuonSet setDM_VBF_EE(outroot, "DM_VBF_EE");
    MuonSet corr_setDM_VBF_EE(outroot, "corr_DM_VBF_EE");
    JetSet setJet_VBF_EE(outroot, "Jet_VBF_EE");*/
    
    HistSetCollection hists = {
    
    // Pre-Selection blocks
    setDM_PreSelection,    corr_setDM_PreSelection,    setJet_PreSelection, corr_setJet_PreSelection,
    /*setDM_PreSelection_BB, corr_setDM_PreSelection_BB, setJet_PreSelection_BB,
    setDM_PreSelection_BO, corr_setDM_PreSelection_BO, setJet_PreSelection_BO,
    setDM_PreSelection_BE, corr_setDM_PreSelection_BE, setJet_PreSelection_BE,
    setDM_PreSelection_OO, corr_setDM_PreSelection_OO, setJet_PreSelection_OO,
    setDM_PreSelection_OE, corr_setDM_PreSelection_OE, setJet_PreSelection_OE,
    setDM_PreSelection_EE, corr_setDM_PreSelection_EE, setJet_PreSelection_EE,*/
    
    // Pre di-jet
    setJet_PreDiJet_VBF, corr_setJet_PreDiJet_VBF,

    // VBF blocks
    setDM_VBF,    corr_setDM_VBF,    setJet_VBF, corr_setJet_VBF
    /*setDM_VBF_BB, corr_setDM_VBF_BB, setJet_VBF_BB,
    setDM_VBF_BO, corr_setDM_VBF_BO, setJet_VBF_BO,
    setDM_VBF_BE, corr_setDM_VBF_BE, setJet_VBF_BE,
    setDM_VBF_OO, corr_setDM_VBF_OO, setJet_VBF_OO,
    setDM_VBF_OE, corr_setDM_VBF_OE, setJet_VBF_OE,
    setDM_VBF_EE, corr_setDM_VBF_EE, setJet_VBF_EE*/
    };

	Streamer streamer(__inputfilename, "Events");
	streamer.chainup();
   
    //MVA Branch outputs 
    gMvaTree = new TTree("mva", "Training features for VBF discriminator");
    gMvaTree->Branch("z",             &g_z);
    gMvaTree->Branch("pt_balance",    &g_ptBalance);
    gMvaTree->Branch("mu1_pt_over_m", &g_mu1_pt_over_m);
    gMvaTree->Branch("mu2_pt_over_m", &g_mu2_pt_over_m);
    gMvaTree->Branch("mu1_eta",       &g_mu1_eta);
    gMvaTree->Branch("mu2_eta",       &g_mu2_eta);
    gMvaTree->Branch("jet1_eta",      &g_jet1_eta);
    gMvaTree->Branch("m_mumu",        &g_m_mumu);
    gMvaTree->Branch("deta_jj",       &g_deta_jj);
    gMvaTree->Branch("dphi_jj",       &g_dphi_jj);
    gMvaTree->Branch("y_mumu",        &g_y_mumu);
    gMvaTree->Branch("y_j1",          &g_y_j1);
    gMvaTree->Branch("y_j2",          &g_y_j2);
    gMvaTree->Branch("pt_mumu",       &g_pt_mumu);
    gMvaTree->Branch("pt_j1",         &g_pt_j1);
    gMvaTree->Branch("pt_j2",         &g_pt_j2);
    gMvaTree->Branch("weight",        &g_weight);
 

    // Loading Muon SF Corrections
    loadMuonCorrections();
    loadJetCorrections();
    loadPileupCorrections();

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

#define ROOT_CHAIN streamer._chain
    
    //Setting Branch Array for LuminosityBlocks Contents
    SET_BRANCH_UINT(run);
    SET_BRANCH_UINT(luminosityBlock);

    //Setting Branch Array for Muon Properties
    SET_BRANCH_UINT(nMuon);
    SET_BRANCH_FLOAT_ARRAY(Muon_pt);
    SET_BRANCH_FLOAT_ARRAY(Muon_bsConstrainedPt);
    SET_BRANCH_FLOAT_ARRAY(Muon_phi);
    SET_BRANCH_FLOAT_ARRAY(Muon_eta);
    SET_BRANCH_INT_ARRAY(Muon_charge);
    SET_BRANCH_BOOL_ARRAY(Muon_isGlobal);
    SET_BRANCH_BOOL_ARRAY(Muon_isTracker);
    SET_BRANCH_UCHAR_ARRAY(Muon_miniIsoId);
    SET_BRANCH_BOOL_ARRAY(Muon_tightId);
    //SET_BRANCH_BOOL_ARRAY(Muon_isPFcand);
    SET_BRANCH_FLOAT_ARRAY(Muon_dxy);
    SET_BRANCH_FLOAT_ARRAY(Muon_dz);
    SET_BRANCH_FLOAT_ARRAY(Muon_pfRelIso04_all);
    SET_BRANCH_UCHAR_ARRAY(Muon_nTrackerLayers);

    //Setting Branch Array for Primary Vertex Properties
    SET_BRANCH_UINT(PV_npvs);
    SET_BRANCH_UINT(PV_npvsGood);
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

    //Setting Branch Array for Trigger Object Properties
    SET_BRANCH_UINT(nTrigObj);
    SET_BRANCH_FLOAT_ARRAY(TrigObj_pt);
    SET_BRANCH_FLOAT_ARRAY(TrigObj_eta);
    SET_BRANCH_FLOAT_ARRAY(TrigObj_phi);
    SET_BRANCH_USHORT_ARRAY(TrigObj_id);
    SET_BRANCH_INT_ARRAY(TrigObj_filterBits);

    //Setting Branch Array for Electron Properties
    SET_BRANCH_UINT(nElectron);
    SET_BRANCH_FLOAT_ARRAY(Electron_pt);
    SET_BRANCH_FLOAT_ARRAY(Electron_eta);
    SET_BRANCH_FLOAT_ARRAY(Electron_phi);
    SET_BRANCH_UCHAR_ARRAY(Electron_cutBased);

    //Setting Branch Array for Rho Properties
    SET_BRANCH_FLOAT(fixedGridRhoFastjetAll);
   
    // --- Event-level weights from NanoAOD ---
    //float Generator_weight = 1.0f;   // falls back to 1 for DATA automatically
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

    // Prefer NanoAODv12 branch names; keep fallbacks for older ntuples
    /*if (ROOT_CHAIN->GetBranch("Generator_weight"))
        ROOT_CHAIN->SetBranchAddress("Generator_weight", &Generator_weight);
    else if (ROOT_CHAIN->GetBranch("genWeight"))
        ROOT_CHAIN->SetBranchAddress("genWeight", &Generator_weight);*/

    if (ROOT_CHAIN->GetBranch("Pileup_nTrueInt"))
        ROOT_CHAIN->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt);
    else if (ROOT_CHAIN->GetBranch("nTrueInt"))
        ROOT_CHAIN->SetBranchAddress("nTrueInt", &Pileup_nTrueInt);


    muons = new std::vector<Muon>();
    jets = new std::vector<Jet>();
    genjets = new std::vector<GenJet>();
    fsrphotons = new std::vector<FsrPhoton>();
    photons = new std::vector<Photon>();
    vertices = new std::vector<Vertex>();
    met = new MET();
    l1s = new std::vector<L1>();
    hlts = new std::vector<HLT>();
    trigobjs = new std::vector<TrigObj>();
    electrons = new std::vector<Electron>();
    rhos = new std::vector<Rho>();
    event = new Event();


	//	Main Event Loop
	uint32_t numEntries = streamer._chain->GetEntries();
    uint32_t updateInterval = numEntries / 100; // Update every 1% of progress
    if (updateInterval == 0) updateInterval = 1; // Handle very small numbers of entries

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

        for (size_t i=0; i<nMuon; i++) {
            Muon mu;
            mu._pt = Muon_pt[i];
            mu._bsConstrainedPt = Muon_bsConstrainedPt[i];
            mu._eta = Muon_eta[i];
            mu._phi = Muon_phi[i];
            mu._charge = Muon_charge[i];
            mu._isGlobal = Muon_isGlobal[i];
            mu._isTracker = Muon_isTracker[i];
            mu._isTight = Muon_tightId[i];
            //mu._isPF = Muon_isPFcand[i];
            mu._dxy = Muon_dxy[i];
            mu._dz = Muon_dz[i];
            mu._pfRelIso04_all = Muon_pfRelIso04_all[i];
            mu._nTrackerLayers = Muon_nTrackerLayers[i];

            muons->push_back(std::move(mu));
        }

        for (size_t i=0; i<nJet; i++) {
            Jet jet;
            jet._pt = Jet_pt[i];
            jet._eta = Jet_eta[i];
            jet._phi = Jet_phi[i];
            jet._mass = Jet_mass[i];
            jet._jetId = Jet_jetId[i];
            jet._btagDeepFlavB = Jet_btagDeepFlavB[i];
            jet._area = Jet_area[i];
            jet._chHEF = Jet_chHEF[i];
            jet._neHEF = Jet_neHEF[i];
            jet._chEmEF = Jet_chEmEF[i];
            jet._neEmEF = Jet_neEmEF[i];
            jet._muEF = Jet_muEF[i];
            jet._genJetIdx = Jet_genJetIdx[i];

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

        for (size_t i=0; i<nTrigObj; i++){
            TrigObj trigobj;
            trigobj._pt = TrigObj_pt[i];
            trigobj._eta = TrigObj_eta[i];
            trigobj._phi = TrigObj_phi[i];
            trigobj._id = TrigObj_id[i];
            trigobj._filterBits = TrigObj_filterBits[i];

            trigobjs->push_back(std::move(trigobj));
        }

        for (size_t i=0; i<nElectron; i++) {
            Electron electron;
            electron._pt = Electron_pt[i];
            electron._eta = Electron_eta[i];
            electron._phi = Electron_phi[i]; 
            electron._cutBased = Electron_cutBased[i];

            electrons->push_back(std::move(electron));
        }

        Rho rhoObj;
        rhoObj._fixedGridRhoFastjetAll = fixedGridRhoFastjetAll; 
        rhos->push_back(rhoObj);

        float rho = rhos->at(0)._fixedGridRhoFastjetAll;

        Vertex vertex;
        vertex._chi2 = PV_chi2;
        vertex._ndof = PV_ndof;
        vertex._score = PV_score;
        vertex._x = PV_x;
        vertex._y = PV_y;
        vertex._z = PV_z;
         
        vertices->push_back(std::move(vertex));
        
        // Pileup from correctionlib (JSON expects: NumTrueInteractions, "weights")
        float pileupWeight = 1.0f;
        if (__isMC && pu_corrections) {
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
        
        // Generator weight straight from the branch
        /*float genWeight = __isMC ? Generator_weight : 1.0f;

        // Final per-event weight
        double base = (__isMC ? sampleBaseWeight(g_sampleName) : 1.0);  // your function
        float Dataweight = static_cast<float>(base) * pileupWeight * genWeight;*/

          /*std::cout << std::setprecision(10) 
          << " genWeight=" << genWeight
          //<< " Generator_weight=" << Generator_weight
          << " Generator_weight=" << Generator_weight_branch
          << " Pileup_nTrueInt=" << Pileup_nTrueInt
          << " pileupWeight=" << pileupWeight
          << " base=" << base
          << " total=" << Dataweight << "\n";*/

        
        //	Selections
		//if (!passVertex(vertices))continue;
        if (!vetoAdditionalLeptons(muons, electrons, __isMC)) continue;

        // Prepare The Pairs Of Muons
        for (auto it = muons->begin(); it != muons->end(); ++it) {
            for (auto jt = (it + 1); jt != muons->end(); ++jt) {
            
            if (!passMuons(*it, *jt, trigobjs)) continue;
            
            categorize(jets, genjets, *it, *jt, fsrphotons, photons, vertices, *met, *event, rhos, hists, Dataweight);
            break;
        }
        break;
    }

    // Update the progress bar every `updateInterval` iterations
    if (i % updateInterval == 0 || i == numEntries - 1) {
            showProgressBar(i + 1, numEntries);
        }
    }

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

void printCuts()
{
        std::cout << endl
        << "Pre-Selection Criteria:" << std::endl
        << "----------------------------------------------------------------------" << std::endl
        << "    Muons:" << std::endl
        << "    -------------" << std::endl
        << "        Muon size >= 2" << std::endl
        << "        Charge of Muon1 != Charge of Muon2" << std::endl
        << "        >=1 muon matched to HLT trigger muon and have pT > 26 GeV" << std::endl
        << std::endl
        << "    Jets:" << std::endl
        << "    -------------" << std::endl
        << "        Jet size >= 2" << std::endl
        << "        pT > 25" << std::endl
        << "        Abs(Eta) < 4.7" << std::endl
        << "        tightId = true" << std::endl
        << "        Loose PU ID or pT > 50 GeV" << std::endl
        << "        Isolated from selected muons by dR(µ, j) > 0.4" << std::endl
        << std::endl
        << "    Electrons:" << std::endl
        << "    -------------" << std::endl
        << "        Electron size = 0" << std::endl
        << std::endl
        << "    Vertices:" << std::endl
        << "    -------------" << std::endl
        << "        Z < 24" << std::endl
        << "        Ndof > 4" << std::endl
        << "        Rho <= 3" << std::endl
        << std::endl
        << "    Leading Muon:" << std::endl
        << "    -------------" << std::endl
        << "    isGlobal && isTracker = true" << std::endl
        << "    tightId = true" << std::endl
        << "    pfRelIso04_all < 0.15" << std::endl
        << "    pT > " << _leadmuonPt << " GeV" << std::endl
        << "    Abs(Eta) < " << _leadmuonEta << std::endl
        << std::endl
        << "    Sub-Leading Muon:" << std::endl
        << "    -----------------" << std::endl
        << "    isGlobal && isTracker = true" << std::endl
        << "    tightId = true" << std::endl
        << "    pfRelIso04_all < 0.25" << std::endl
        << "    pT > " << _subleadmuonPt << " GeV" << std::endl
        << "    Abs(Eta) < " << _subleadmuonEta << std::endl
        << std::endl
        << std::endl
        << "VBF-Selection Criteria:" << std::endl
        << "----------------------------------------------------------------------" << std::endl
        << "    Events:" << std::endl
        << "    -------------" << std::endl
        << "        Excluding events w/ additional muon (electron) with pT > 20 GeV and |η| < 2.4 (2.5)" << std::endl
        << std::endl
        << "    Jets:" << std::endl
        << "    -------------" << std::endl
        << "        At most one b-tagged loose jet and no b-tagged-medium jets" << std::endl
        << std::endl
        << "        Leading Jet:" << std::endl
        << "        ------------" << std::endl
        << "        pT > " << _leadJetPt << " GeV" << std::endl
        << std::endl
        << "        Sub-Leading Jet:" << std::endl
        << "        ----------------" << std::endl
        << "        pT > " << _subleadJetPt << " GeV" << std::endl
        << std::endl
        << "        Di-Jet:" << std::endl
        << "        -------" << std::endl
        << "        Mass > " << _dijetMass_VBF << " GeV" << std::endl
        << "        Abs(Eta) > " << _dijetdEta_VBF << std::endl
        << std::endl;
}

int main(int argc, char** argv)
{
	// Register Signals
	signal(SIGABRT, &sigHandler);
	signal(SIGTERM, &sigHandler);
	signal(SIGINT, &sigHandler);

	std::string none;
	bool genPUMC = false;

	// Pare Options
	po::options_description desc("Allowed Program Options");
	desc.add_options()
		("help", "produce help messages")
		("input", po::value<std::string>(), "a file specifying all the ROOT files to process")
		("isMC", po::value<bool>(), "type of data: DATA vs MC")
		("output", po::value<std::string>(), "output file name")
		("genPUMC", po::value<bool>(&genPUMC)->default_value(false), "true if should generate the MC PU file")
		("puMC", po::value<std::string>(&__puMCfilename)->default_value("None"), "MC PU Reweight file")
		("puDATA", po::value<std::string>(&__puDATAfilename)->default_value("None"), "DATA PU Reweight file")
        ("leadmuonPt", po::value<double>(&_leadmuonPt)->default_value(_leadmuonPt), "Muon Pt Cut")
        ("subleadmuonPt", po::value<double>(&_subleadmuonPt)->default_value(_subleadmuonPt), "Muon Pt Cut")
        ("leadmuonEta", po::value<double>(&_leadmuonEta)->default_value(_leadmuonEta), "Muon Eta Cut")
        ("subleadmuonEta", po::value<double>(&_subleadmuonEta)->default_value(_subleadmuonEta), "Muon Eta Cut")
        ("leadmuonIso", po::value<double>(&_leadmuonIso)->default_value(_leadmuonIso), "Muon Isolation Cut")
        ("subleadmuonIso", po::value<double>(&_subleadmuonIso)->default_value(_subleadmuonIso), "Muon Isolation Cut")
        ("leadJetPt", po::value<double>(&_leadJetPt)->default_value(_leadJetPt), "Lead Jet Pt Cut")
        ("subleadJetPt", po::value<double>(&_subleadJetPt)->default_value(_subleadJetPt), "SubLeading Jet Pt Cut")
        ("leadJetEta", po::value<double>(&_leadJetEta)->default_value(_leadJetEta), "Lead Jet Eta Cut")
        ("subleadJetEta", po::value<double>(&_subleadJetEta)->default_value(_subleadJetEta), "SubLeading Jet Eta Cut")
        ("dijetMass_VBF", po::value<double>(&_dijetMass_VBF)->default_value(_dijetMass_VBF), "DiJet Mass VBF-Category Cut")
        ("dijetdEta_VBF", po::value<double>(&_dijetdEta_VBF)->default_value(_dijetdEta_VBF), "DiJet deta VBF-Category Cut")
	    ("weights", po::value<std::string>(&g_weightsFile)->default_value("None"), "path to weights config (key=value)")
        ("sample",  po::value<std::string>(&g_sampleName)->default_value("DATA"), "sample name (e.g. DY, TT, WZ, VBF, DATA)")
        ("era", po::value<std::string>(&g_era)->default_value("12022"), "Era section from weights.ini (e.g. 12022, 22022, 12023, 22023)")
        ("puDir",  po::value<std::string>(&g_puDir)->default_value(g_puDir), "Directory containing PU JSONs")
        ("puJSON", po::value<std::string>(&g_puJSON)->default_value("auto"), "PU JSON file path or 'auto' to use puWeights<era>.json")
        ("puKey",  po::value<std::string>(&g_puKey)->default_value("auto"),  "Correction key inside PU JSON or 'auto' to use era-based default")
        ("--muDir",     po::value<std::string>(&g_muDir)->default_value(g_muDir), "Directory containing muon JSONs")
        ("--muSFJSON",  po::value<std::string>(&g_muSFJSON)->default_value("auto"), "Muon SF JSON (path) or 'auto' => muon_Z<era>.json")
        ("--muSFKeyID", po::value<std::string>(&g_muSFKeyID)->default_value(g_muSFKeyID),  "Muon ID SF key")
        ("--muSFKeyISO",po::value<std::string>(&g_muSFKeyISO)->default_value(g_muSFKeyISO),"Muon ISO SF key")
        ("--muSFKeyTRG",po::value<std::string>(&g_muSFKeyTRG)->default_value(g_muSFKeyTRG),"Muon trigger SF key")
        ("--rochJSON",  po::value<std::string>(&g_rochJSON)->default_value("auto"),        "Rochester JSON (path) or 'auto' => RoccoR<era>.json")
        ("--rochKeyDT", po::value<std::string>(&g_rochKeyDT)->default_value(g_rochKeyDT),  "Rochester data key")
        ("--rochKeyMC", po::value<std::string>(&g_rochKeyMC)->default_value(g_rochKeyMC),  "Rochester MC key")
        ("--jmeDir",     po::value<std::string>(&g_jmeDir)->default_value(g_jmeDir), "Directory containing JME JSONs")
        ("--jercJSON",   po::value<std::string>(&g_jercJSON)->default_value("auto"), "Path to JERC JSON, or 'auto' => jet_jerc<era>.json")
        ("--jetVetoJSON",po::value<std::string>(&g_jetVetoJSON)->default_value("auto"), "Path to veto map JSON, or 'auto' => jetvetomaps<era>.json")
        ("--jetIdJSON",  po::value<std::string>(&g_jetIdJSON)->default_value("auto"), "Path to jet ID JSON, or 'auto' => jetid<era>.json")
    ;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
    
    if (vm.count("help") || argc<2)
	{
		std::cout << desc << std::endl;
		return 1;
	}

    if (!vm.count("output")) {
    std::cerr << "Error: --output option is missing!" << std::endl;
    exit(1);
    }
	
	//	Assign globals
	__inputfilename = vm["input"].as<std::string>();
	__isMC = vm["isMC"].as<bool>();
	__outputfilename = vm["output"].as<std::string>();
	__genPUMC = vm["genPUMC"].as<bool>();
	__puMCfilename = vm["puMC"].as<std::string>();
	__puDATAfilename = vm["puDATA"].as<std::string>();
    _leadmuonPt = vm["leadmuonPt"].as<double>();
    _subleadmuonPt = vm["subleadmuonPt"].as<double>();
    _leadmuonEta = vm["leadmuonEta"].as<double>();
    _subleadmuonEta = vm["subleadmuonEta"].as<double>();
    _leadmuonIso = vm["leadmuonIso"].as<double>();
    _subleadmuonIso = vm["subleadmuonIso"].as<double>();
    _leadJetPt = vm["leadJetPt"].as<double>();
    _subleadJetPt = vm["subleadJetPt"].as<double>();
    _leadJetEta = vm["leadJetEta"].as<double>();
    _subleadJetEta = vm["subleadJetEta"].as<double>();
    _dijetMass_VBF = vm["dijetMass_VBF"].as<double>();
    _dijetdEta_VBF = vm["dijetdEta_VBF"].as<double>();
    
    // --- Load weights & pick lumi for the chosen era (only for MC) ---
    if (__isMC && g_weightsFile != "None") {
        loadWeightsFile(g_weightsFile);   // fills g_sectionLumis / g_sectionSamples

        // Choose lumi from the section named by --era (e.g. 12022, 22022, 12023, 22023)
        auto sec = g_sectionLumis.find(g_era);
        if (sec != g_sectionLumis.end()) {
            auto &lmap = sec->second;
            // Prefer canonical "lumi" if present; otherwise fall back to any "luminosity_*" in that section
            if (lmap.count("lumi")) {
                g_selected_lumi = lmap.at("lumi");
            } else {
                // fallback: first luminosity_* key if user used specific names like luminosity_BCD/EFG/D
                for (const auto &kv : lmap) {
                    if (kv.first.rfind("luminosity", 0) == 0) { g_selected_lumi = kv.second; break; }
                }
            }
            if (g_selected_lumi <= 0.0) {
                std::cerr << "[weights] No valid lumi found in section [" << g_era
                          << "]; using default " << g_selected_lumi << "\n";
            }
        } else {
            std::cerr << "[weights] Era section [" << g_era
                      << "] not found in weights file; using default lumi = "
                      << g_selected_lumi << "\n";
        }
    }

    printCuts();

    
    // Start Main Loop Process 
	loadMuonCorrections();
    if (__isMC) {
        loadJetCorrections();               // JEC/JER for MC
        loadPileupCorrections();            // PU weights for MC
    } else {
        // Optional: loadJetCorrections() if you want DATA residuals via JSON.
        loadJetCorrections();
    }
    process();
	return 0;
}
