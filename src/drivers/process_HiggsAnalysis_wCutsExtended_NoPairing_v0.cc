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

//*******************************************************************
//                   HISTOGRAM MACROS
//*******************************************************************

// DiMuon set (No Weights) 
#define FILL_DIMUON_HIST(set, weight) \
    set.hDiMuonpt->Fill(corr_p4dimuon.Pt(), weight); \
    set.hDiMuonMass->Fill(corr_p4dimuon.M(), weight); \
    set.hDiMuonMass_Z->Fill(corr_p4dimuon.M(), weight); \
    set.hDiMuoneta->Fill(corr_p4dimuon.Eta(), weight); \
    set.hDiMuondphi->Fill(dphi, weight); \
    set.hDiMuonMass_SB_SR1->Fill(corr_p4dimuon.M(), weight); \
    set.hMuonpt->Fill(mu1._pt, weight); \
    set.hMuonpt->Fill(mu2._pt, weight); \
    set.hMuonDXY->Fill(mu1._dxy, weight); \
    set.hMuonDXY->Fill(mu2._dxy, weight); \
    set.hMuonDZ->Fill(mu1._dz, weight); \
    set.hMuonDZ->Fill(mu2._dz, weight); \
    set.hMuonDXY_vs_MuonDZ->Fill(mu1._dxy, mu1._dz, weight); \
    set.hMuonDXY_vs_MuonDZ->Fill(mu2._dxy, mu2._dz, weight); \
    set.hDimuoneta_vs_DimuonMass->Fill(corr_p4dimuon.M(), corr_p4dimuon.Eta(), weight); \
    set.hDimuonpt_vs_DimuonMass->Fill(corr_p4dimuon.M(), corr_p4dimuon.Pt(), weight); \
    set.hDimuonpt_vs_DimuonMass_JPsi->Fill(corr_p4dimuon.M(), corr_p4dimuon.Pt(), weight); \
    set.hDimuonpt_vs_DimuonMass_Upsilon->Fill(corr_p4dimuon.M(), corr_p4dimuon.Pt(), weight); \
    set.hDimuonpt_vs_DimuonMass_Z->Fill(corr_p4dimuon.M(), corr_p4dimuon.Pt(), weight); \
    set.hDimuonpt_vs_DimuonMass_Higgs->Fill(corr_p4dimuon.M(), corr_p4dimuon.Pt(), weight); \
    if (corr_p4dimuon.M() > 2.9 && corr_p4dimuon.M() < 3.3) { \
        set.hDiMuoneta_JPsi->Fill(corr_p4dimuon.Eta(), weight); \
    } \
    if (corr_p4dimuon.M() > 9.2 && corr_p4dimuon.M() < 10.5) { \
        set.hDiMuoneta_Upsilon->Fill(corr_p4dimuon.Eta(), weight); \
    } \
    if (corr_p4dimuon.M() > 88 && corr_p4dimuon.M() < 95) { \
        set.hDiMuoneta_Z->Fill(corr_p4dimuon.Eta(), weight); \
    } \
    if (corr_p4dimuon.M() > 122 && corr_p4dimuon.M() < 128) { \
        set.hDiMuoneta_Higgs->Fill(corr_p4dimuon.Eta(), weight); \
    }

// DiMuon set (Weights) 
#define FILL_CORR_DIMUON_HIST(set, weight) FILL_DIMUON_HIST(set, weight)

// DiJet set (Weights) 
#define FILL_DIJET_HIST(set, weight) \
    set.hDiJetpt->Fill(dijetpt, weight); \
    set.hDiJetMass->Fill(dijetmass, weight); \
    set.hDiJetdeta->Fill(dijetdeta, weight);

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

// Kinematic Cuts  
enum MuonType {kLead=0, kSubLead=1};  // Muon Id --> lead = 0, sublead =1
double _leadmuonPt = 29.;
double _subleadmuonPt = 20.;
double _leadmuonEta = 2.4;
double _subleadmuonEta = 2.4;
double _leadmuonIso = 0.15;
double _subleadmuonIso = 0.25;
double _leadJetPt = 35.;
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
std::shared_ptr<correction::CorrectionSet> corrections = nullptr;
std::shared_ptr<correction::CorrectionSet> muoncorrections = nullptr;
void loadMuonCorrections() {
    if (!corrections)
        corrections = correction::CorrectionSet::from_file
        ("/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/7_Corrections/muon_Z.json");
    if (!muoncorrections)
        muoncorrections = correction::CorrectionSet::from_file
        ("/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/7_Corrections/RoccoR2022.json");
}

// -------------------------------------
//             JETS 
// -------------------------------------

// Jet Energy Corrections/Resolutions (JEC/JER) 
std::shared_ptr<correction::CorrectionSet> jetCorrections = nullptr;
std::shared_ptr<correction::CorrectionSet> jetidCorrections = nullptr;
std::shared_ptr<correction::CorrectionSet> jetVetoMaps = nullptr;
std::shared_ptr<const correction::Correction> jetVetoCorr = nullptr;
void loadJetCorrections() {
    if (!jetCorrections)
        jetCorrections = correction::CorrectionSet::from_file
        ("/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/7_Corrections/jet_jerc.json");
    if (!jetVetoMaps)
        jetVetoMaps = correction::CorrectionSet::from_file
        ("/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/7_Corrections/jetvetomaps.json");
    if (!jetVetoCorr) {
        jetVetoCorr = jetVetoMaps->at("Summer22_23Sep2023_RunCD_V1");
    }
    //if (!jetidCorrections)
        //jetidCorrections = correction::CorrectionSet::from_file
        //("/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/7_Corrections/jetid.json");
}

// -------------------------------------
//             PU WEIGHTS 
// -------------------------------------

std::shared_ptr<correction::CorrectionSet> pileupCorrections = nullptr;
std::shared_ptr<const correction::Correction> pu_corrections = nullptr;
void loadPileupCorrections() {
    if (!pileupCorrections)
        pileupCorrections = correction::CorrectionSet::from_file
        ("/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/7_Corrections/puWeightsBCD.json");
        pu_corrections = pileupCorrections->at("Collisions2022_355100_357900_eraBCD_GoldenJson");
}

TH1D *hEventWeights = NULL;

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
        return 1.0;
    }

    float pt     = mu._pt;
    float eta    = mu._eta;
    float phi    = mu._phi;
    int   charge = mu._charge;
    int s = 0;
    int m = 0;

    try {
        if (!isMC) {  // DATA
            return muoncorrections->at("kScaleDT")->evaluate({charge, pt, eta, phi, s, m});
        } else {      // MC
            return muoncorrections->at("kScaleMC")->evaluate({charge, pt, eta, phi, s, m});
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Failed to evaluate Rochester correction: " << e.what() << std::endl;
        return 1.0;
    }
}

// -------------------------------------
//             JETS 
// -------------------------------------

float getCorrectedJetPt(const Jet& jet, float rho, bool isMC, float genJetPt, bool hasGenJet) {
    if (!jetCorrections) {
        std::cerr << "[ERROR] Jet corrections not loaded!\n";
        return jet._pt;
    }

    float eta  = std::clamp(jet._eta, -5.191f, 5.191f);
    float pt   = std::max(jet._pt, 0.1f);
    float area = jet._area;

    // JEC prefixes (MC vs DATA)
    std::string era = isMC
    ? "Summer22_22Sep2023_V2_MC"
    : "Summer22_22Sep2023_RunCD_V2_DATA";

    float l1  = jetCorrections->at(era + "_L1FastJet_AK4PFPuppi")->evaluate({area, eta, pt, rho});
    float l2  = jetCorrections->at(era + "_L2Relative_AK4PFPuppi")->evaluate({eta, pt});
    float l3  = jetCorrections->at(era + "_L3Absolute_AK4PFPuppi")->evaluate({eta, pt});
    float res = isMC ? 1.0f : jetCorrections->at(era + "_L2L3Residual_AK4PFPuppi")->evaluate({eta, pt});

    float jecFactor = l1 * l2 * l3 * res;
    float corr_pt   = pt * jecFactor;

    // ---- JER (MC only) ----
    if (isMC) {
        // Use JRV1 keys from your JSON
        float sf = jetCorrections
            ->at("Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi")
            ->evaluate({eta, "nom"}); // "up"/"down" for systematics

        if (hasGenJet && genJetPt > 0.f) {
            // Truth-anchored smearing
            float smearFactor = 1.0f + (sf - 1.0f) * (corr_pt - genJetPt) / corr_pt;
            // Optional safety: avoid negative scaling
            smearFactor = std::max(smearFactor, 0.0f);
            corr_pt *= smearFactor;
        } else {
            // Stochastic smearing
            float resolution = jetCorrections
                ->at("Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi")
                ->evaluate({eta, corr_pt, rho});
            float sigma = resolution * std::sqrt(std::max(sf*sf - 1.0f, 0.0f));
            // Consider a deterministic seed (event, jet index) if you run multithreaded
            static TRandom3 rand(42);
            corr_pt *= (1.0f + rand.Gaus(0, sigma));
            corr_pt  = std::max(corr_pt, 0.0f); // keep physical
        }
    }

    return corr_pt;
}


inline bool isVetoedJet(float eta, float phi) {
    if (!jetVetoCorr) { std::cerr << "[ERROR] Jet veto corr not loaded!\n"; return false; }
    const double val = jetVetoCorr->evaluate({std::string("jetvetomap"), (double)eta, (double)phi});
    return (val != 0.0);
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
//TRIGGER CRITERIA
// -------------------------------------

/*bool passL1(L1s* l1trigger) {
    if (l1trigger->empty() || !(*l1trigger)[0]._SingleMu12er1p5) {
        return false;  // Fail, go to next event
    }
    return true;  // Pass
}

bool passHLT(HLTs* hlttrigger) {
    if (hlttrigger->empty() || !(*hlttrigger)[0]._IsoMu24) {
        return false;  // Fail, go to next event
    }
    return true;  // Pass
}*/

// -------------------------------------
//LEPTON CRITERIA 
// -------------------------------------

bool vetoAdditionalLeptons(Muons* muons, Electrons* electrons)
{
    int numPassingMuons = 0;

    // Check muons
    for (const auto& mu : *muons) {
        float rochSF = getMuonRochester(mu, __isMC);
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

        // Print vertex information
/*        if (rho <= .2 && it->_ndof > 4 && TMath::Abs(it->_z) < 24) {
        std::cout << "Vertex " << n 
                  << ": z = " << it->_z
                  << ", ndof = " << it->_ndof
                  << ", rho = " << rho << " mm" << std::endl;}*/
        
        if (TMath::Abs(it->_z) < 24 
            && it->_ndof > 4 
            && rho <= 0.3)
			return true;
        n++;
	}
	return false;
}

// -------------------------------------
//MUON CRITERIA 
// -------------------------------------

bool passMuonHLT(TrigObjs* trig, Muon const& m, bool isMC)
{
	if (trig->size()==0)
		return false;

    // Apply Rochester correction
    float rochSF = getMuonRochester(m, __isMC);
    float corrected_pt = m._pt * rochSF;

    if (corrected_pt < 26)
        return false;
	
	for (TrigObjs::const_iterator it=trig->begin(); it!=trig->end(); ++it)
	{
        // Print trigger object information
/*        std::cout << "Trigger Object: "
                  << ", id = " << it->_id
                  << ", eta = " << it->_eta
                  << ", pt = " << it->_pt 
                  << ", phi = " << it->_phi << std::endl;*/
        
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

    float pt1_corr = m1._pt * rochSF1;
    float pt2_corr = m2._pt * rochSF2;
    
    const Muon& leadmuon     = (pt1_corr > pt2_corr) ? m1 : m2;
    const Muon& subleadmuon  = (pt1_corr > pt2_corr) ? m2 : m1;
    float lead_pt_corr       = std::max(pt1_corr, pt2_corr);
    float sublead_pt_corr    = std::min(pt1_corr, pt2_corr);
	
    if (m1._charge!=m2._charge && passMuon(leadmuon, lead_pt_corr, kLead) && passMuon(subleadmuon, sublead_pt_corr, kSubLead))
        if (passMuonHLT(trigObjs, m1, !__isMC) || passMuonHLT(trigObjs, m2, !__isMC))
            return true;

	return false;
}

// -------------------------------------
//JET CRITERIA 
// -------------------------------------

/*bool passJetIDfromJSON(const Jet& jet) {
    if (!jetidCorrections) {
        std::cerr << "[ERROR] Jet ID corrections not loaded!\n";
        return false;
    }

    try {
        const auto& corr = jetidCorrections->at("AK4PUPPI_Tight");

        double result = corr->evaluate({
            jet._eta,
            jet._chHEF,
            jet._neHEF,
            jet._chEmEF,
            jet._neEmEF,
            jet._muEF,
            static_cast<int>(jet._chMult),
            static_cast<int>(jet._neMult),
            static_cast<int>(jet._chMult + jet._neMult)
        });

        return (result >= 0.5);  // convention: 1 = pass, 0 = fail
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] Failed to evaluate jet ID: " << e.what() << std::endl;
        return false;
    }
}*/

bool passJet(Jets* jets, GenJets* genjets, float rho)
{
    if (jets->size() < 2) return false;

    std::vector<TLorentzVector> Alljets;
    for (const auto& jet : *jets) {
        
        // VETO MAP (drop the jet if it sits in a bad region)
        if (isVetoedJet(jet._eta, jet._phi)) continue;

        // gen-match locals, JEC+JER
        float genJetPt = 0.0f;
        bool hasGenJet = false;

        if (__isMC && jet._genJetIdx >= 0 && jet._genJetIdx < (int)genjets->size()) {
            const auto& gj = genjets->at(jet._genJetIdx);
            genJetPt = gj._pt;
            hasGenJet = true;
        }

        float corr_pt   = getCorrectedJetPt(jet, rho, __isMC, genJetPt, hasGenJet);
        float corr_mass = jet._mass * (corr_pt / jet._pt);

        if (corr_pt <= 25 || std::abs(jet._eta) >= 4.7) continue;
        if ((jet._jetId & 4) != 4) continue; // Tight ID

        TLorentzVector p4;
        p4.SetPtEtaPhiM(corr_pt, jet._eta, jet._phi, corr_mass);
        Alljets.push_back(p4);
    }

    if (Alljets.size() < 2) return false;
    if (Alljets[0].Pt() < 50) return false;
    return true;
}

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
        id_sf   = corrections->at("NUM_TightID_DEN_TrackerMuons")->evaluate({eta, pt, "nominal"});
        iso_sf  = corrections->at("NUM_TightPFIso_DEN_TightID")->evaluate({eta, pt, "nominal"});
        if (applyTriggerSF && pt >= 26) {
            trig_sf = corrections->at("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight")->evaluate({eta, pt, "nominal"});
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
analysis::core::Rhos* rhos, HistSetCollection& hists, float MCweight) {

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
   
    // -------------------------------------
    //             FSR PHOTONS 
    // -------------------------------------
    
    // Add FSR photons to both muons if any pass criteria
    auto applyFSR = [](const analysis::core::Muon& mu, TLorentzVector& p4, TLorentzVector& corr_p4, const analysis::core::FsrPhotons* fsrphotons, int muonIdx) {
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
   
   //MIGHT BE GOOD TO PLACE A MACRO HERE TO SEE FSR CORRECTIONS
    
    // -------------------------------------
    //              JETS  
    // -------------------------------------
    
    bool pass_bjet_veto = true;
    int n_loose_bjets = 0;
    int n_medium_bjets = 0;
    float rho = rhos->at(0)._fixedGridRhoFastjetAll;
    
    std::vector<TLorentzVector> p4jets;
    for (Jets::const_iterator it=jets->begin(); it!=jets->end(); ++it)
	{
        float genJetPt = -1.0f;
        bool hasGenJet = false;

        // JER matching using NanoAOD index
        if (__isMC && it->_genJetIdx >= 0 && it->_genJetIdx < (int)genjets->size()) {
            const auto& gj = genjets->at(it->_genJetIdx);
            genJetPt = gj._pt;
            hasGenJet = true;
        }
        
        float corr_pt   = getCorrectedJetPt(*it, rho, __isMC, genJetPt, hasGenJet);
        float corr_mass = it->_mass * (corr_pt / it->_pt);

        // Muon cleaning ΔR > 0.4 
        if (jetMuondR(it->_eta, it->_phi, mu1._eta, mu1._phi) < 0.4) continue; 
		if (jetMuondR(it->_eta, it->_phi, mu2._eta, mu2._phi) < 0.4) continue; 
        
        // Count b-tagged jets
        if (it->_btagDeepFlavB > 0.0583) n_loose_bjets++;
        if (it->_btagDeepFlavB > 0.3086) n_medium_bjets++;
        
        // Jet passed all cuts → store for dijet selection
		TLorentzVector p4;
		p4.SetPtEtaPhiM(corr_pt, it->_eta, it->_phi, corr_mass);
		p4jets.push_back(p4);
			
            // Print out the jet's properties
            /*std::cout << "Jet passing selection: " << std::endl;
            std::cout << "  pT:  " << it->_pt << " GeV" << std::endl;
            std::cout << "  eta: " << it->_eta << std::endl;
            std::cout << "  phi: " << it->_phi << std::endl;
            std::cout << "  mass: " << it->_mass << " GeV" << std::endl;*/
    }
   
    // TLorentz Vectors For Leading/Sub-Leading Jets
    if (p4jets.size() < 2) return;
    
    TLorentzVector p4lead = p4jets[0]; 
    TLorentzVector p4sub = p4jets[1];
    TLorentzVector dijet = p4lead + p4sub;
    
    float dijetdeta = TMath::Abs(p4lead.Eta() - p4sub.Eta());
    float dijetmass = dijet.M();
    float dijetpt = dijet.Pt();

    float mu1_SF = getMuonSF(mu1, __isMC);
    float mu2_SF = getMuonSF(mu2, __isMC);
    float totalWeight = MCweight * mu1_SF * mu2_SF;
    float noWeight = 1.0;
    
    // -------------------------------------
    //        PRE-SELECTION HISTOGRAMS  
    // -------------------------------------
    
    FILL_DIMUON_HIST(hists.setDM_PreSelection, noWeight);   
    FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection, totalWeight); 
    FILL_DIJET_HIST(hists.setJet_PreSelection, MCweight); 

    if (isBarrel(mu1) && isBarrel(mu2))
    {
        //	BB
        FILL_DIMUON_HIST(hists.setDM_PreSelection_BB, noWeight);
        FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection_BB, totalWeight);
        FILL_DIJET_HIST(hists.setJet_PreSelection_BB, MCweight); 
    }
    else if ((isBarrel(mu1) && isOverlap(mu2)) ||
        (isBarrel(mu2) && isOverlap(mu1)))
    {
        //	BO
        FILL_DIMUON_HIST(hists.setDM_PreSelection_BO, noWeight);
        FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection_BO, totalWeight);
        FILL_DIJET_HIST(hists.setJet_PreSelection_BO, MCweight); 
    }
    else if ((isBarrel(mu1) & isEndcap(mu2)) || 
        (isBarrel(mu2) && isEndcap(mu1)))
    {
        //	BE
        FILL_DIMUON_HIST(hists.setDM_PreSelection_BE, noWeight);
        FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection_BE, totalWeight);
        FILL_DIJET_HIST(hists.setJet_PreSelection_BE, MCweight); 
    }
    else if (isOverlap(mu1) && isOverlap(mu2))
    {
        //	OO
        FILL_DIMUON_HIST(hists.setDM_PreSelection_OO, noWeight);
        FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection_OO, totalWeight);
        FILL_DIJET_HIST(hists.setJet_PreSelection_OO, MCweight); 
    }
    else if ((isOverlap(mu1) && isEndcap(mu2)) ||
        (isOverlap(mu2) && isEndcap(mu1)))
    {
        //	OE
        FILL_DIMUON_HIST(hists.setDM_PreSelection_OE, noWeight);
        FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection_OE, totalWeight);
        FILL_DIJET_HIST(hists.setJet_PreSelection_OE, MCweight); 
    }
    else if (isEndcap(mu1) && isEndcap(mu2))
    {
        //	EE
        FILL_DIMUON_HIST(hists.setDM_PreSelection_EE, noWeight);
        FILL_CORR_DIMUON_HIST(hists.corr_setDM_PreSelection_EE, totalWeight);
        FILL_DIJET_HIST(hists.setJet_PreSelection_EE, MCweight); 
    }

    // -------------------------------------
    //          VBF HISTOGRAMS  
    // -------------------------------------

    if (n_loose_bjets >= 2 || n_medium_bjets >= 1)
        pass_bjet_veto = false;
	
    if (pass_bjet_veto) {
        if (p4lead.Pt() > _leadJetPt && p4sub.Pt() > _subleadJetPt && dijetmass > _dijetMass_VBF && TMath::Abs(dijetdeta) > _dijetdEta_VBF) {
            
            //	VBF SELECTION CRITERIA
            FILL_DIMUON_HIST(hists.setDM_VBF, noWeight)
            FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF, totalWeight)
            FILL_DIJET_HIST(hists.setJet_VBF, MCweight); 

            if (isBarrel(mu1) && isBarrel(mu2))
            {
                //	BB
                FILL_DIMUON_HIST(hists.setDM_VBF_BB, noWeight);
                FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF_BB, totalWeight)
                FILL_DIJET_HIST(hists.setJet_VBF_BB, MCweight); 
            }
            else if ((isBarrel(mu1) && isOverlap(mu2)) ||
                (isBarrel(mu2) && isOverlap(mu1)))
            {
                //	BO
                FILL_DIMUON_HIST(hists.setDM_VBF_BO, noWeight);
                FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF_BO, totalWeight)
                FILL_DIJET_HIST(hists.setJet_VBF_BO, MCweight); 
            }
            else if ((isBarrel(mu1) & isEndcap(mu2)) || 
                (isBarrel(mu2) && isEndcap(mu1)))
            {
                //	BE
                FILL_DIMUON_HIST(hists.setDM_VBF_BE, noWeight);
                FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF_BE, totalWeight)
                FILL_DIJET_HIST(hists.setJet_VBF_BE, MCweight); 
            }
            else if (isOverlap(mu1) && isOverlap(mu2))
            {
                //	OO
                FILL_DIMUON_HIST(hists.setDM_VBF_OO, noWeight);
                FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF_OO, totalWeight)
                FILL_DIJET_HIST(hists.setJet_VBF_OO, MCweight); 
            }
            else if ((isOverlap(mu1) && isEndcap(mu2)) ||
                (isOverlap(mu2) && isEndcap(mu1)))
            {
                //	OE
                FILL_DIMUON_HIST(hists.setDM_VBF_OE, noWeight);
                FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF_OE, totalWeight)
                FILL_DIJET_HIST(hists.setJet_VBF_OE, MCweight); 
            }
            else if (isEndcap(mu1) && isEndcap(mu2))
            {
                //	EE
                FILL_DIMUON_HIST(hists.setDM_VBF_EE, noWeight);
                FILL_CORR_DIMUON_HIST(hists.corr_setDM_VBF_EE, totalWeight)
                FILL_DIJET_HIST(hists.setJet_VBF_EE, MCweight); 
            }
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
//                        GENERATING PU WEGIHT (NO LONGER NEED) 
//*******************************************************************

/*void generatePUMC()
{
    std::cout << "Generating Pileup Monte-Carlo File..." << std::endl;
    std::cout << "[DEBUG] Writing PU histogram to: " << __puMCfilename << std::endl
    << std::endl;
	TFile *pufile = new TFile(__puMCfilename.c_str(), "recreate");
	TH1D *h = new TH1D("pileup", "pileup", 99, 0, 99);

	Streamer s(__inputfilename, "Events");
	s.chainup();
    
    int pileup_nPU = 0;
    float genWeight = 1.0;
    
    s._chain->SetBranchAddress("Pileup_nPU", &pileup_nPU);
    s._chain->SetBranchAddress("genWeight", &genWeight);

	int numEvents = s._chain->GetEntries();
	std::cout << "Total entries to process for generating PUMC: " << numEvents << std::endl;
    for (uint32_t i=0; i<numEvents; i++)
	{
		s._chain->GetEntry(i);
        if (i < 10) {
        std::cout << "[DEBUG] Event " << i
                  << " Pileup_nPU = " << &pileup_nPU
                  << " genWeight = " << genWeight << std::endl;
    }
		h->Fill(pileup_nPU, genWeight);
	}

	pufile->Write();
	pufile->Close();
}*/

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
//                        MAIN PROGRESS FUCNTION 
//*******************************************************************

void process()
{
    //	Output Root File
    TFile *outroot = new TFile(__outputfilename.c_str(), "recreate");
    hEventWeights = new TH1D("eventWeights", "eventWeights", 1, 0, 1);
	
    DimuonSet setDM_PreSelection(outroot, "DM_Pre-Selection");
    DimuonSet corr_setDM_PreSelection(outroot, "corr_DM_Pre-Selection");
    JetSet setJet_PreSelection(outroot, "Jet_Pre-Selection");
    DimuonSet setDM_PreSelection_BB(outroot, "DM_Pre-Selection_BB");
    DimuonSet corr_setDM_PreSelection_BB(outroot, "corr_DM_Pre-Selection_BB");
    JetSet setJet_PreSelection_BB(outroot, "Jet_Pre-Selection_BB");
    DimuonSet setDM_PreSelection_BO(outroot, "DM_Pre-Selection_BO");
    DimuonSet corr_setDM_PreSelection_BO(outroot, "corr_DM_Pre-Selection_BO");
    JetSet setJet_PreSelection_BO(outroot, "Jet_Pre-Selection_BO");
    DimuonSet setDM_PreSelection_BE(outroot, "DM_Pre-Selection_BE");
    DimuonSet corr_setDM_PreSelection_BE(outroot, "corr_DM_Pre-Selection_BE");
    JetSet setJet_PreSelection_BE(outroot, "Jet_Pre-Selection_BE");
    DimuonSet setDM_PreSelection_OO(outroot, "DM_Pre-Selection_OO");
    DimuonSet corr_setDM_PreSelection_OO(outroot, "corr_DM_Pre-Selection_OO");
    JetSet setJet_PreSelection_OO(outroot, "Jet_Pre-Selection_OO");
    DimuonSet setDM_PreSelection_OE(outroot, "DM_Pre-Selection_OE");
    DimuonSet corr_setDM_PreSelection_OE(outroot, "corr_DM_Pre-Selection_OE");
    JetSet setJet_PreSelection_OE(outroot, "Jet_Pre-Selection_OE");
    DimuonSet setDM_PreSelection_EE(outroot, "DM_Pre-Selection_EE");
    DimuonSet corr_setDM_PreSelection_EE(outroot, "corr_DM_Pre-Selection_EE");
    JetSet setJet_PreSelection_EE(outroot, "Jet_Pre-Selection_EE");
    DimuonSet setDM_VBF(outroot, "DM_VBF");
    DimuonSet corr_setDM_VBF(outroot, "corr_DM_VBF");
    JetSet setJet_VBF(outroot, "Jet_VBF");
    DimuonSet setDM_VBF_BB(outroot, "DM_VBF_BB");
    DimuonSet corr_setDM_VBF_BB(outroot, "corr_DM_VBF_BB");
    JetSet setJet_VBF_BB(outroot, "Jet_VBF_BB");
    DimuonSet setDM_VBF_BO(outroot, "DM_VBF_BO");
    DimuonSet corr_setDM_VBF_BO(outroot, "corr_DM_VBF_BO");
    JetSet setJet_VBF_BO(outroot, "Jet_VBF_BO");
    DimuonSet setDM_VBF_BE(outroot, "DM_VBF_BE");
    DimuonSet corr_setDM_VBF_BE(outroot, "corr_DM_VBF_BE");
    JetSet setJet_VBF_BE(outroot, "Jet_VBF_BE");
    DimuonSet setDM_VBF_OO(outroot, "DM_VBF_OO");
    DimuonSet corr_setDM_VBF_OO(outroot, "corr_DM_VBF_OO");
    JetSet setJet_VBF_OO(outroot, "Jet_VBF_OO");
    DimuonSet setDM_VBF_OE(outroot, "DM_VBF_OE");
    DimuonSet corr_setDM_VBF_OE(outroot, "corr_DM_VBF_OE");
    JetSet setJet_VBF_OE(outroot, "Jet_VBF_OE");
    DimuonSet setDM_VBF_EE(outroot, "DM_VBF_EE");
    DimuonSet corr_setDM_VBF_EE(outroot, "corr_DM_VBF_EE");
    JetSet setJet_VBF_EE(outroot, "Jet_VBF_EE");
    
    HistSetCollection hists = {
    setDM_PreSelection, corr_setDM_PreSelection, setJet_PreSelection,
    setDM_PreSelection_BB, corr_setDM_PreSelection_BB, setJet_PreSelection_BB,
    setDM_PreSelection_BO, corr_setDM_PreSelection_BO, setJet_PreSelection_BO,
    setDM_PreSelection_BE, corr_setDM_PreSelection_BE, setJet_PreSelection_BE,
    setDM_PreSelection_OO, corr_setDM_PreSelection_OO, setJet_PreSelection_OO,
    setDM_PreSelection_OE, corr_setDM_PreSelection_OE, setJet_PreSelection_OE,
    setDM_PreSelection_EE, corr_setDM_PreSelection_EE, setJet_PreSelection_EE,
    setDM_VBF, corr_setDM_VBF, setJet_VBF,
    setDM_VBF_BB, corr_setDM_VBF_BB, setJet_VBF_BB,
    setDM_VBF_BO, corr_setDM_VBF_BO, setJet_VBF_BO,
    setDM_VBF_BE, corr_setDM_VBF_BE, setJet_VBF_BE,
    setDM_VBF_OO, corr_setDM_VBF_OO, setJet_VBF_OO,
    setDM_VBF_OE, corr_setDM_VBF_OE, setJet_VBF_OE,
    setDM_VBF_EE, corr_setDM_VBF_EE, setJet_VBF_EE
    };

    /*reweight::LumiReWeighting* weighter = nullptr;
    if (__isMC && !__genPUMC) {  
    std::cout << "mcPU = " << __puMCfilename << std::endl;
    std::cout << "dataPU = " << __puDATAfilename << std::endl;

    TString mc_pileupfile = __puMCfilename.c_str();
    TString data_pileupfile = __puDATAfilename.c_str();

    weighter = new reweight::LumiReWeighting(
        mc_pileupfile.Data(), data_pileupfile.Data(), "pileup", "pileup");
    } else if (__genPUMC) {
        std::cout << "[INFO] Generating MC PU histogram only. No PU reweighting will be applied.\n";
    }*/

	Streamer streamer(__inputfilename, "Events");
	streamer.chainup();
	
    //	generate the MC Pileup histogram
	//if (__genPUMC && __isMC)
		//generatePUMC();

    auto numEventsWeighted = 100;
    hEventWeights->Fill(0.5, numEventsWeighted);
    
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
    
    //Setting Branch Array for Muon Properties
    SET_BRANCH_UINT(nMuon);
    SET_BRANCH_FLOAT_ARRAY(Muon_pt);
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
    SET_BRANCH_UCHAR_ARRAY(Electron_cutBased);

    //Setting Branch Array for Rho Properties
    SET_BRANCH_FLOAT(fixedGridRhoFastjetAll);
    

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

        for (size_t i=0; i<nMuon; i++) {
            Muon mu;
            mu._pt = Muon_pt[i];
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

        for (size_t i=0; i<nGenJet; i++) {
            GenJet genjet;
            genjet._pt = GenJet_pt[i];
            genjet._eta = GenJet_eta[i];
            genjet._phi = GenJet_phi[i];
            genjet._mass = GenJet_mass[i];
            
            genjets->push_back(std::move(genjet));
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
            electron._phi = Electron_cutBased[i];

            electrons->push_back(std::move(electron));
        }

        Rho rhoObj;
        rhoObj._fixedGridRhoFastjetAll = fixedGridRhoFastjetAll; 
        rhos->push_back(rhoObj);

        float rho = rhos->at(0)._fixedGridRhoFastjetAll;

        //for (size_t i=0; i<PV_npvs; i++){
            Vertex vertex;
            vertex._chi2 = PV_chi2;
            vertex._ndof = PV_ndof;
            vertex._score = PV_score;
            vertex._x = PV_x;
            vertex._y = PV_y;
            vertex._z = PV_z;
         
            vertices->push_back(std::move(vertex));
        //}
        
/*        L1 l1trigger;
        l1trigger._SingleMu12er1p5 = L1_SingleMu12er1p5;
        l1s->push_back(std::move(l1trigger));

        HLT hlttrigger;
        hlttrigger._IsoMu24 = HLT_IsoMu24;
        hlts->push_back(std::move(hlttrigger));
        
        std::cout << "L1 Triggers:\n";
        for (const auto& l1 : *l1s) {
            std::cout << "L1_DoubleMu_12_5: " << l1._DoubleMu_12_5 << std::endl;
        }

        std::cout << "HLT Triggers:\n";
        for (const auto& hlt : *hlts) {
            std::cout << "HLT_IsoMu24: " << hlt._IsoMu24 << std::endl;
        }

        met->_pt = MET_pt;

        if (__isMC && !__genPUMC && weighter && aux)
            pileupWeight = weighter->weight(aux->_nPU);
*/
        float pileupWeight = 1.0;
        if (__isMC && pu_corrections && aux) {
            try {
                pileupWeight = pu_corrections->evaluate({aux->_nPU, "nominal"});
            } catch (const std::exception& e) {
                std::cerr << "[PU ERROR] " << e.what() << std::endl;
                pileupWeight = 1.0;
            }
        }
        float genWeight = (__isMC && aux) ? aux->_genWeight : 1.0;
        
        //----- LUMINOSITY -----

        float luminosity_BCD = 7.980454145e3; 
        
        //----- CROSS-SECTION -----
        
        //Background:
        //-------------------------------------- 
        //float DYcrosssection_BCD = 6.688e3;
        //float EWKcrosssection_BCD = 1.0; //CHANGE LATER 
        //float TTcrosssection_BCD = 762.1;
        //float WZcrosssection_BCD = 7.568;
        //float WZto3LNucrosssection_BCD = 4.924;
        //float ZZto2L2Nucrosssection_BCD = 1.031;
        //float ZZto2L2Qcrosssection_BCD = 6.788;
        
        //Signal:
        //-------------------------------------- 
        //float GGHcrosssection_BCD = 32.69;
        float VBFcrosssection_BCD = 4.18; 
        //float TTHcrosssection_BCD = 0.5742;
        //float ZHcrosssiection_BCD = 0.5478;
        
        //----- TOTAL "genEventSumw" ------
        
        //Background:
        //-------------------------------------- 
        //float DYgenEventSumw_BCD = 1.17668e12;
        //float EWKgenEventSumw_BCD = 8.51310e5;
        //float TTgenEventSumw_BCD = 1.9129e9;
        //float WZgenEventSumw_BCD = 3.1537e7;
        //float WZto3LNugenEventSumw_BCD = 1.3561e7;
        //float ZZto2L2NugenEventSumw_BCD = 1.50043e7;
        //float ZZto2L2QgenEventSumw_BCD = 9.95594e7;
        
        //float Ext_DYgenEventSumw_BCD = 1.73021e12;
        //float Ext_TTgenEventSumw_BCD = 1.93459e9;
        //float Ext_WZgenEventSumw_BCD = 3.07671e7;
        //float Ext_WZto3LNugenEventSumw_BCD = 4.37112e7;
        //float Ext_ZZto2L2NugenEventSumw_BCD = 1.73775e7;
        //float Ext_ZZto2L2QgenEventSumw_BCD = 1.01821e8;
        
        //Signal:
        //-------------------------------------- 
        //float GGHgenEventSumw_BCD = 6.53894e6;
        float VBFgenEventSumw_BCD = 1.66873e6;
        //float TTHgenEventSumw_BCD = 1.13023e5;
        //float ZHgenEventSumw_BCD = 8.09902e5;
        
        //----- FINAL MC-WEIGHT ------
       
        //std::cout << "pileupWeight" << pileupWeight << std::endl;
        //std::cout << "genWeight " << genWeight << std::endl;
        
        //float MCweight  = ((luminosity_BCD * VBFcrosssection_BCD)/(VBFgenEventSumw_BCD + Ext_VBFgenEventSumw_BCD)) * pileupWeight * genWeight;
        float MCweight  = ((luminosity_BCD * VBFcrosssection_BCD)/(VBFgenEventSumw_BCD)) * pileupWeight * genWeight;
        //float MCweight  = 1.0;
 		
        //	Selections
        //if (!passL1(l1s)) continue;  
        //if (!passHLT(hlts)) continue;
        //if (!FsrPhoton(fsrphotons)) continue; HAVEN'T SET THIS UP YET
		if (!passVertex(vertices))continue;
        if (!passJet(jets, genjets, rho)) continue;
        if (!vetoAdditionalLeptons(muons, electrons)) continue;

        // Prepare The Pairs Of Muons
        for (auto it = muons->begin(); it != muons->end(); ++it) {
            for (auto jt = (it + 1); jt != muons->end(); ++jt) {
            
            if (!passMuons(*it, *jt, trigobjs)) continue;
            
            categorize(jets, genjets, *it, *jt, fsrphotons, photons, vertices, *met, *event, rhos, hists, MCweight);
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
        << "        Pt > 25" << std::endl
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
        << "    Pt > " << _leadmuonPt << " GeV" << std::endl
        << "    Abs(Eta) < " << _leadmuonEta << std::endl
        << std::endl
        << "    Sub-Leading Muon:" << std::endl
        << "    -----------------" << std::endl
        << "    isGlobal && isTracker = true" << std::endl
        << "    tightId = true" << std::endl
        << "    pfRelIso04_all < 0.25" << std::endl
        << "    Pt > " << _subleadmuonPt << " GeV" << std::endl
        << "    Abs(Eta) < " << _subleadmuonEta << std::endl
        << std::endl
        << std::endl
        << "VBF-Selection Criteria:" << std::endl
        << "----------------------------------------------------------------------" << std::endl
        << "    Jets:" << std::endl
        << "    -------------" << std::endl
        << "        At most one b-tagged loose jet and no b-tagged-medium jets" << std::endl
        << std::endl
        << "        Leading Jet:" << std::endl
        << "        ------------" << std::endl
        << "        Pt > " << _leadJetPt << " GeV" << std::endl
        << std::endl
        << "        Sub-Leading Jet:" << std::endl
        << "        ----------------" << std::endl
        << "        Pt > " << _subleadJetPt << " GeV" << std::endl
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
    printCuts();

    
    // Start Main Loop Process 
	loadMuonCorrections();
    loadJetCorrections();
    loadPileupCorrections();
    process();
	return 0;
}
