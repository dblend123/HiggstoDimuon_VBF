
//	
#ifdef STANDALONE
#include "Muon.h"
#include "Jet.h"
#include "Electron.h"
#include "Tau.h"
#include "Vertex.h"
#include "Event.h"
#include "MET.h"
#include "Constants.h"
#include "Streamer.h"
#include "MetaHiggs.h"

//	ROOT headers
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "LumiReweightingStandAlone.h"
#include "HistogramSets.h"

#include "boost/program_options.hpp"
#include <signal.h>

#define FILL_JETS(set) \
	set.hDiJetMass->Fill(dijetmass, puweight); \
	set.hDiJetdeta->Fill(TMath::Abs(jetsdeta), puweight); \
	set.hDiMuonpt->Fill(p4dimuon.Pt(), puweight); \
	set.hDiMuonMass->Fill(p4dimuon.M(), puweight); \
	set.hDiMuoneta->Fill(p4dimuon.Eta(), puweight); \
	set.hDiMuondphi->Fill(dimuon_dphi, puweight); \
	set.hMuonpt->Fill(p4m1.Pt(), puweight); \
	set.hMuonpt->Fill(p4m2.Pt(), puweight); \
	set.hMuoneta->Fill(p4m1.Eta(), puweight); \
	set.hMuoneta->Fill(p4m2.Eta(), puweight); \
	set.hMuonphi->Fill(p4m1.Phi(), puweight); \
	set.hMuonphi->Fill(p4m2.Phi(), puweight);

#define FILL_NOJETS(set) \
	set.hDiMuonpt->Fill(p4dimuon.Pt(), puweight); \
	set.hDiMuonMass->Fill(p4dimuon.M(), puweight); \
	set.hDiMuoneta->Fill(p4dimuon.Eta(), puweight); \
	set.hDiMuondphi->Fill(dimuon_dphi, puweight); \
	set.hMuonpt->Fill(p4m1.Pt(), puweight); \
	set.hMuonpt->Fill(p4m2.Pt(), puweight); \
	set.hMuoneta->Fill(p4m1.Eta(), puweight); \
	set.hMuoneta->Fill(p4m2.Eta(), puweight); \
	set.hMuonphi->Fill(p4m1.Phi(), puweight); \
	set.hMuonphi->Fill(p4m2.Phi(), puweight);

/*
 *	Declare/Define all the service globals
 */
std::string __inputfilename;
std::string __outputfilename;
bool __isMC;
bool __genPUMC;
std::string __puMCfilename;
std::string __puDATAfilename;
bool __continueRunning = true;

/*
 *  Define all the Constants
 *  for Muon:
 *  lead - id 0
 *  sublead - is 1
 */
enum MuonType {kLead=0, kSubLead=1};
double __cutmuonHLTMatchedPt = 24.;
double __cutmuonHLTMatchedEta = 2.4;

double __cutleadmuonPt = 10.;
double __cutsubleadmuonPt = 10.;
double __cutleadmuonEta = 2.4;
double __cutsubleadmuonEta = 2.4;
double __cutleadmuonIso = 0.1;
double __cutsubleadmuonIso = 0.1;

double __cutextramuonPt = 10.;
double __cutextramuonEta = 2.4;
double __cutextramuonIso = 0.1;

double __cutleadJetPt = 40.;
double __cutsubleadJetPt = 30.;

double __cutdijetMass_VBFTight = 500;
double __cutjetBTagMediumDiscr = 0.8;

double __cutjetsdeta_VBFTight = 2.5;
double __cutdijetMass_VBFLoose = 250.;
double __cutjetsdeta_VBFLoose = 2.5;
double __cutlowdijetMassVhadH = 60.;
double __cuthighdijetMassVhadH = 110.;
double __cutMETPtZinvH = 40.;
double __cutdimuonPtggFTight = 25.;

double __cutElePt = 10.;
double __cutTauPt = 20.;
double __cutmaxExtraLeptons = 4;

std::string const NTUPLEMAKER_NAME =  "ntuplemaker_H2DiMuonMaker";

namespace po = boost::program_options;
using namespace analysis::core;
using namespace analysis::dimuon;
using namespace analysis::processing;

TH1D *hEventWeights = NULL;
DimuonSet setNoCats("NoCats");

DimuonSet set2Mu1bjets("2Mu1bjets");
DimuonSet set2Mu1bjets0l("2Mu1bjets0l");

DimuonSet set2Mu1bjets2l("2Mu1bjets2l");
DimuonSet set2Mu1bjets2l2e("2Mu1bjets2l2e");
DimuonSet set2Mu1bjets2l1e1mu("2Mu1bjets2l1e1mu");
DimuonSet set2Mu1bjets2l2mu("2Mu1bjets2l2mu");
DimuonSet set2Mu1bjets2lMissed("2Mu1bjets2lMissed");
DimuonSet set2Mu1bjetsMissed("2Mu1bjetsMissed");
DimuonSet set2Mu0bjets("2Mu0bjets");
DimuonSet set2Mu0bjets0l("2Mu0bjets0l");
DimuonSet set2Mu0bjets0l2jets("2Mu0bjets0l2jets");
DimuonSet set2Mu0bjets0l2jetsVBFTight("2Mu0bjets0l2jetsVBFTight");
DimuonSet set2Mu0bjets0l2jetsVBFLoose("2Mu0bjets0l2jetsVBFLoose");
DimuonSet set2Mu0bjets0l2jetsVhadH("2Mu0bjets0l2jetsVhadH");
DimuonSet set2Mu0bjets0l2jetsggF("2Mu0bjets0l2jetsggF");
DimuonSet set2Mu0bjets0l01jets("2Mu0bjets0l01jets");
DimuonSet set2Mu0bjets0l01jetsZinvH("2Mu0bjets0l01jetsZinvH");
DimuonSet set2Mu0bjets0l01jetsggFTight("2Mu0bjets0l01jetsggFTight");
DimuonSet set2Mu0bjets0l01jetsggFTightBB("2Mu0bjets0l01jetsggFTightBB");
DimuonSet set2Mu0bjets0l01jetsggFTightBO("2Mu0bjets0l01jetsggFTightBO");
DimuonSet set2Mu0bjets0l01jetsggFTightBE("2Mu0bjets0l01jetsggFTightBE");
DimuonSet set2Mu0bjets0l01jetsggFTightOO("2Mu0bjets0l01jetsggFTightOO");
DimuonSet set2Mu0bjets0l01jetsggFTightOE("2Mu0bjets0l01jetsggFTightOE");
DimuonSet set2Mu0bjets0l01jetsggFTightEE("2Mu0bjets0l01jetsggFTightEE");
DimuonSet set2Mu0bjets0l01jetsggFLoose("2Mu0bjets0l01jetsggFLoose");
DimuonSet set2Mu0bjets0l01jetsggFLooseBB("2Mu0bjets0l01jetsggFLooseBB");
DimuonSet set2Mu0bjets0l01jetsggFLooseBO("2Mu0bjets0l01jetsggFLooseBO");
DimuonSet set2Mu0bjets0l01jetsggFLooseBE("2Mu0bjets0l01jetsggFLooseBE");
DimuonSet set2Mu0bjets0l01jetsggFLooseOO("2Mu0bjets0l01jetsggFLooseOO");
DimuonSet set2Mu0bjets0l01jetsggFLooseOE("2Mu0bjets0l01jetsggFLooseOE");
DimuonSet set2Mu0bjets0l01jetsggFLooseEE("2Mu0bjets0l01jetsggFLooseEE");
DimuonSet set2Mu0bjets12l("2Mu0bjets12l");
DimuonSet set2Mu0bjets12l1e("2Mu0bjets12l1e");
DimuonSet set2Mu0bjets12l1mu("2Mu0bjets12l1mu");
DimuonSet set2Mu0bjets12l1mu1e("2Mu0bjets12l1mu1e");
DimuonSet set2Mu0bjets12l2e("2Mu0bjets12l2e");
DimuonSet set2Mu0bjets12l2mu("2Mu0bjets12l2mu");

/*
DimuonSet set2Mu1bJets2e("2Mu1bJets2e");
DimuonSet set2Mu1bJets2Mu("2Mu1bJets2Mu");
DimuonSet set2Mu1bJets0e("2Mu1bJets0e");

DimuonSet set2Mu0bJets("2Mu0bJets");
DimuonSet set2Mu0bJets1e("2Mu0bJets1e");
DimuonSet set2Mu0bJets2e("2Mu0bJets2e");
DimuonSet set2Mu0bJets1Mu("2Mu0bJets1Mu");
DimuonSet set2Mu0bJets1e1Mu("2Mu0bJets1e1Mu");
DimuonSet set2Mu0bJets2Mu("2Mu0bJets2Mu");
DimuonSet set2Mu0bJets0e("2Mu0bJets0e");

DimuonSet set2Mu0bJets1e1Tau("2Mu0bJets1e1Tau");
DimuonSet set2Mu0bJets1e0Tau("2Mu0bJets1e0Tau");

DimuonSet set2Mu0bJets1Mu1Tau("2Mu0bJets1Mu1Tau");
DimuonSet set2Mu0bJets1Mu0Tau("2Mu0bJets1Mu0Tau");

DimuonSet set2Mu0bJets0e2Tau("2Mu0bJets0e2Tau");
DimuonSet set2Mu0bJets0e1Tau("2Mu0bJets0e1Tau");
//DimuonSet set2Mu0bJets("2Mu0bJets");

DimuonSet set2Jets("2Jets");
DimuonSet setVBFTight("VBFTight");
DimuonSet setggFTight("ggFTight");
DimuonSet setggFLoose("ggFLoose");
DimuonSet set2Mu0bJets0e2Jets("2Mu0bJets0e2Jets");
DimuonSet set2Mu0bJets0e2JetsVBFTight("2Mu0bJets0e2JetsVBFTight");
DimuonSet set2Mu0bJets0e2JetsggFTight("2Mu0bJets0e2JetsggFTight");
DimuonSet set2Mu0bJets0e2JetsggFLoose("2Mu0bJets0e2JetsggFLoose");

DimuonSet set01Jets("01Jets");
DimuonSet set01JetsTight("01JetsTight");
DimuonSet set01JetsLoose("01JetsLoose");
DimuonSet set01JetsBB("01JetsBB");
DimuonSet set01JetsBO("01JetsBO");
DimuonSet set01JetsBE("01JetsBE");
DimuonSet set01JetsOO("01JetsOO");
DimuonSet set01JetsOE("01JetsOE");
DimuonSet set01JetsEE("01JetsEE");
DimuonSet set01JetsTightBB("01JetsTightBB");
DimuonSet set01JetsTightBO("01JetsTightBO");
DimuonSet set01JetsTightBE("01JetsTightBE");
DimuonSet set01JetsTightOO("01JetsTightOO");
DimuonSet set01JetsTightOE("01JetsTightOE");
DimuonSet set01JetsTightEE("01JetsTightEE");
DimuonSet set01JetsLooseBB("01JetsLooseBB");
DimuonSet set01JetsLooseBO("01JetsLooseBO");
DimuonSet set01JetsLooseBE("01JetsLooseBE");
DimuonSet set01JetsLooseOO("01JetsLooseOO");
DimuonSet set01JetsLooseOE("01JetsLooseOE");
DimuonSet set01JetsLooseEE("01JetsLooseEE");
*/

/*
 *  Declare all the collections that are used in analysis
 */
EventAuxiliary *aux = NULL;
Muons *muons = NULL;
Electrons *electrons = NULL;
Taus *taus = NULL;
Jets *jets = NULL;
Vertices *vertices = NULL;
Event *event = NULL;
MET *met = NULL;
Muons muons_extra;

/*
 *  Define all the auxiliary functions
 */
bool isBarrel(Muon const& m)
{
	return TMath::Abs(m._eta)<0.8;
}

bool isOverlap(Muon const& m)
{
	return TMath::Abs(m._eta)>=0.8 && TMath::Abs(m._eta)<1.6;
}

bool isEndcap(Muon const& m)
{
	return TMath::Abs(m._eta)>=1.6 && TMath::Abs(m._eta)<2.1;
}

bool passVertex(Vertices* v)
{
	if (v->size()==0)
		return false;
	int n=0;
	for (Vertices::const_iterator it=v->begin();
		it!=v->end() && n<20; ++it)
	{
		if (TMath::Abs(it->_z)<24 &&
			it->_ndf>4)
			return true;
		n++;
	}

	return false;
}

bool passMuon(Muon const& m, MuonType id=kLead)
{
    if (id==kLead)
    {
	    if (m._isGlobal && m._isTracker &&
		    m._pt>__cutleadmuonPt && TMath::Abs(m._eta)<__cutleadmuonEta &&
		    m._isTight && (m._trackIsoSumPt/m._pt)<__cutleadmuonIso)
		    return true;
    }
    else
    {
	    if (m._isGlobal && m._isTracker &&
		    m._pt>__cutsubleadmuonPt && TMath::Abs(m._eta)<__cutsubleadmuonEta &&
		    m._isTight && (m._trackIsoSumPt/m._pt)<__cutsubleadmuonIso)
		    return true;
    }
	return false;
}

bool passMuonHLT(Muon const& m)
{
	if ((m._isHLTMatched[1] || m._isHLTMatched[0]) &&
		m._pt>__cutmuonHLTMatchedPt && TMath::Abs(m._eta)<__cutmuonHLTMatchedEta)
		return true;

	return false;
}

bool passMuons(Muon const& m1, Muon const& m2)
{
    Muon const& leadmuon = m1._pt>m2._pt ? m1 : m2;
    Muon const& subleadmuon = m1._pt<m2._pt ? m1 : m2;
	if (m1._charge!=m2._charge &&
		passMuon(leadmuon, kLead) && passMuon(subleadmuon, kSubLead))
		return true;

	return false;
}

float jetMuondR(float jeta,float jphi, float meta, float mphi)
{
	TLorentzVector p4j,p4m;
	p4j.SetPtEtaPhiM(10, jeta, jphi, 0);
	p4m.SetPtEtaPhiM(10, meta, mphi, 0);
	return p4j.DeltaR(p4m);
}

float getdeta(TLorentzVector const& v1, TLorentzVector v2)
{
    return v1.Eta()>v2.Eta() ? v1.Eta()-v2.Eta() : v2.Eta()-v1.Eta();
}

void copy(Muons* muons, Muons& muons_extra, std::pair<int, int> const& p)
{
    int counter = 0;
    for (Muons::const_iterator it=muons->begin(); it!=muons->end(); ++it)
    {
        if (counter!=p.first && counter!=p.second)
            muons_extra.push_back(*it);
        counter++;
    }
}

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
		<< "#events processed total = " << numEvents << std::endl
		<< "#events weighted total = " << numEventsWeighted << std::endl;

	return numEventsWeighted;
}

void generatePUMC()
{
    std::cout << "### Generate PU MC file...." << std::endl;
	TFile *pufile = new TFile(__puMCfilename.c_str(), "recreate");
	TH1D *h = new TH1D("pileup", "pileup", 50, 0, 50);

	Streamer s(__inputfilename, NTUPLEMAKER_NAME+"/Events");
	s.chainup();

	s._chain->SetBranchAddress("EventAuxiliary", &aux);
	int numEvents = s._chain->GetEntries();
	for (uint32_t i=0; i<numEvents; i++)
	{
		s._chain->GetEntry(i);
		h->Fill(aux->_nPU, aux->_genWeight);
	}

	pufile->Write();
	pufile->Close();
}

void process()
{
	//	out ...
	TFile *outroot = new TFile(__outputfilename.c_str(), "recreate");
    hEventWeights = new TH1D("eventWeights", "eventWeights", 1, 0, 1);
    /*
	setNoCats.init();
	set2Jets.init();
	set01Jets.init();
	setVBFTight.init();
	setggFTight.init();
	setggFLoose.init();
	set01JetsTight.init();
	set01JetsLoose.init();
	set01JetsBB.init();
	set01JetsBO.init();
	set01JetsBE.init();
	set01JetsOO.init();
	set01JetsOE.init();
	set01JetsEE.init();
	set01JetsTightBB.init();
	set01JetsTightBO.init();
	set01JetsTightBE.init();
	set01JetsTightOO.init();
	set01JetsTightOE.init();
	set01JetsTightEE.init();
	set01JetsLooseBB.init();
	set01JetsLooseBO.init();
	set01JetsLooseBE.init();
	set01JetsLooseOO.init();
	set01JetsLooseOE.init();
	set01JetsLooseEE.init();
    */

    setNoCats.init();
    set2Mu1bjets.init();
    set2Mu1bjets0l.init();
    set2Mu1bjets2l.init();
    set2Mu1bjets2l2e.init();
    set2Mu1bjets2l1e1mu.init();
    set2Mu1bjets2l2mu.init();
    set2Mu1bjets2lMissed.init();
    set2Mu1bjetsMissed.init();
    set2Mu0bjets.init();
    set2Mu0bjets0l.init();
    set2Mu0bjets0l2jets.init();
    set2Mu0bjets0l2jetsVBFTight.init();
    set2Mu0bjets0l2jetsVBFLoose.init();
    set2Mu0bjets0l2jetsVhadH.init();
    set2Mu0bjets0l2jetsggF.init();
    set2Mu0bjets0l01jets.init();
    set2Mu0bjets0l01jetsZinvH.init();
    set2Mu0bjets0l01jetsggFTight.init();
    set2Mu0bjets0l01jetsggFTightBB.init();
    set2Mu0bjets0l01jetsggFTightBO.init();
    set2Mu0bjets0l01jetsggFTightBE.init();
    set2Mu0bjets0l01jetsggFTightOO.init();
    set2Mu0bjets0l01jetsggFTightOE.init();
    set2Mu0bjets0l01jetsggFTightEE.init();
    set2Mu0bjets0l01jetsggFLoose.init();
    set2Mu0bjets0l01jetsggFLooseBB.init();
    set2Mu0bjets0l01jetsggFLooseBO.init();
    set2Mu0bjets0l01jetsggFLooseBE.init();
    set2Mu0bjets0l01jetsggFLooseOO.init();
    set2Mu0bjets0l01jetsggFLooseOE.init();
    set2Mu0bjets0l01jetsggFLooseEE.init();
    set2Mu0bjets12l.init();
    set2Mu0bjets12l1e.init();
    set2Mu0bjets12l1mu.init();
    set2Mu0bjets12l1mu1e.init();
    set2Mu0bjets12l2e.init();
    set2Mu0bjets12l2mu.init();

	//	get the total events, etc...
	long long int numEventsWeighted = sampleinfo(__inputfilename);
    hEventWeights->Fill(0.5, numEventsWeighted);

	//	generate the MC Pileup histogram
	if (__genPUMC && __isMC)
		generatePUMC();

	Streamer streamer(__inputfilename, NTUPLEMAKER_NAME+"/Events");
	streamer.chainup();

    //  branch
    muons = NULL;
    electrons = NULL;
    taus = NULL;
    jets = NULL;
    vertices = NULL;
    event = NULL;
    aux = NULL;
    met = NULL;
	streamer._chain->SetBranchAddress("Muons", &muons);
    streamer._chain->SetBranchAddress("Electrons", &electrons);
    streamer._chain->SetBranchAddress("Taus", &taus);
	streamer._chain->SetBranchAddress("Jets", &jets);
	streamer._chain->SetBranchAddress("Vertices", &vertices);
	streamer._chain->SetBranchAddress("Event", &event);
	streamer._chain->SetBranchAddress("EventAuxiliary", &aux);
	streamer._chain->SetBranchAddress("MET", &met);

	//	init the PU reweighter
	reweight::LumiReWeighting *weighter = NULL;
	if (__isMC)
	{
		std::cout << "mcPU=" << __puMCfilename << "  dataPU="
			<< __puDATAfilename << std::endl;
		TString mc_pileupfile = __puMCfilename.c_str();
		TString data_pileupfile = __puDATAfilename.c_str();
		weighter = new reweight::LumiReWeighting(
		mc_pileupfile.Data(), data_pileupfile.Data(), "pileup", "pileup");
	}

	//	Event Loop
	uint32_t numEntries = streamer._chain->GetEntries();
	for (uint32_t i=0; i<numEntries && __continueRunning; i++)
	{
        //  
        //  clean the vectors/variables used for event processing
        //
        muons_extra.clear();
        Muon mu1, mu2;

        //
        //  start the event processsing
        //
		streamer._chain->GetEntry(i);
		if (i%1000==0)
			std::cout << "### Event " << i << " / " << numEntries
				<< std::endl;

        //  pu wieght
		float puweight = __isMC ? weighter->weight(aux->_nPU)*aux->_genWeight :
			1.;

		//
		//	Event Selections:
        //	- Primary Vertex Cuts
        //	- HLT Path fired
		//
		if (!passVertex(vertices))
			continue;
		if (!(aux->_hasHLTFired[0] || aux->_hasHLTFired[1]))
			continue;

        //  
        //  1. Find the 2 Muons to be considered as Higgs Candidate
        //  Separate them from the rest of muons.
        //  Skeep the event if there are no such candidates
        //  actual dimuon_higgscandidate is actually set inside of the function
        //
        //  2. Check if there is at least 1 muon that matches HLT and pt>cut
        //  
        std::vector<std::pair<int ,int> > vgoodmus;
        bool matchedHLTMuon = false;
        int icounter = 0;
        int numGoodMus = 0;
        for (Muons::const_iterator it=muons->begin(); it!=muons->end();
            ++it)
        {
            //  HLT matching check
            if (passMuonHLT(*it))
                matchedHLTMuon=true;

            //  
            int jcounter = icounter+1;
            for (Muons::const_iterator jt=(it+1);
                jt!=muons->end(); ++jt)
            {
                if (passMuons(*it, *jt))
                    vgoodmus.push_back(std::make_pair(icounter, jcounter));
                jcounter++;
            }

            icounter++;
            if (passMuon(*it))
                numGoodMus++;
        }

        //  
        //  check that we do have a Higgs Candidate Pair, otherwise next event
        //
        if (vgoodmus.size()==0) continue;
        else if (vgoodmus.size()==1)
        {
            mu1 = muons->at(vgoodmus[0].first);
            mu2 = muons->at(vgoodmus[0].second);
            copy(muons, muons_extra, vgoodmus[0]);
        }
        else // more than 1 pair of Higgs Candidate
        {
            if (numGoodMus==4 && vgoodmus.size()==3)
            {
                continue;
            }
            else if (numGoodMus==4 && vgoodmus.size()==2)
            {
                // 2 different pairs
                std::pair<int, int> p1 = vgoodmus[0];
                std::pair<int, int> p2 = vgoodmus[1];
                TLorentzVector p411, p412, p421, p422, p41, p42;
                p411.SetPtEtaPhiM(muons->at(p1.first)._pt, 
                    muons->at(p1.first)._eta, muons->at(p1.first)._phi,
                    PDG_MASS_Mu);
                p412.SetPtEtaPhiM(muons->at(p1.second)._pt, 
                    muons->at(p1.second)._eta, muons->at(p1.second)._phi,
                    PDG_MASS_Mu);
                p421.SetPtEtaPhiM(muons->at(p2.first)._pt, 
                    muons->at(p2.first)._eta, muons->at(p2.first)._phi,
                    PDG_MASS_Mu);
                p422.SetPtEtaPhiM(muons->at(p2.second)._pt, 
                    muons->at(p2.second)._eta, muons->at(p2.second)._phi,
                    PDG_MASS_Mu);
                p41 = p411 + p412;
                p42 = p421 + p422;
                if (TMath::Abs(p41.M()-PDG_MASS_Z) < TMath::Abs(p42.M() - PDG_MASS_Z))
                {
                    mu1 = muons->at(p2.first);
                    mu2 = muons->at(p2.second);
                    copy(muons, muons_extra, p2);
                }
                else
                {
                    mu1 = muons->at(p1.first);
                    mu2 = muons->at(p1.second);
                    copy(muons, muons_extra, p1);
                }
            }
            else if (numGoodMus==3)
            {
                std::pair<int, int> p1 = vgoodmus[0];
                std::pair<int, int> p2 = vgoodmus[1];
                TLorentzVector p411, p412, p421, p422, p41, p42;
                p411.SetPtEtaPhiM(muons->at(p1.first)._pt, 
                    muons->at(p1.first)._eta, muons->at(p1.first)._phi,
                    PDG_MASS_Mu);
                p412.SetPtEtaPhiM(muons->at(p1.second)._pt, 
                    muons->at(p1.second)._eta, muons->at(p1.second)._phi,
                    PDG_MASS_Mu);
                p421.SetPtEtaPhiM(muons->at(p2.first)._pt, 
                    muons->at(p2.first)._eta, muons->at(p2.first)._phi,
                    PDG_MASS_Mu);
                p422.SetPtEtaPhiM(muons->at(p2.second)._pt, 
                    muons->at(p2.second)._eta, muons->at(p2.second)._phi,
                    PDG_MASS_Mu);
                p41 = p411 + p412;
                p42 = p421 + p422;
                if (p41.Pt()>p42.Pt()) // take with the largest Pt
                {
                    mu1 = muons->at(p1.first);
                    mu2 = muons->at(p1.second);
                    copy(muons, muons_extra, p1);
                }
                else
                {
                    mu1 = muons->at(p2.first);
                    mu2 = muons->at(p2.second);
                    copy(muons, muons_extra, p2);
                }
            }
            else continue;
        }

        //
        //  HLT matching check
        //
        if (!matchedHLTMuon) continue;

        //  
        //  create the dimuons 4vector
        //
        TLorentzVector p4m1, p4m2, p4dimuon;
        p4m1.SetPtEtaPhiM(mu1._pt, mu1._eta, mu1._phi, PDG_MASS_Mu);
        p4m2.SetPtEtaPhiM(mu2._pt, mu2._eta, mu2._phi, PDG_MASS_Mu);
        p4dimuon = p4m1+p4m2;
        double dimuon_dphi = p4m1.DeltaPhi(p4m2);

        //  
        //  Deal with electrons
        //

        //  
        //  Deal with Taus
        //

        //  
        //  Determine the number of extra leptons and make sure the number is less than
        //  max
        //
        int numElectrons = electrons->size();
        int numTaus = taus->size();
        int numExtraMuons = muons->size()-2;
        int numExtraLeptons = numElectrons + numTaus + numExtraMuons;
        if (!(numExtraLeptons<=__cutmaxExtraLeptons)) continue;

        //  
        //  separate jets into 2 vectors:
        //  - bjets
        //  - not bjets
        //
        std::vector<TLorentzVector> p4bjets, p40bjets;
        for (Jets::const_iterator it=jets->begin(); it!=jets->end(); ++it)
        {
            if (it->_pt>30 && TMath::Abs(it->_eta)<2.7 && 
                it->_btag[0]>__cutjetBTagMediumDiscr)
            {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(it->_pt, it->_eta, it->_phi, it->_mass);
                p4bjets.push_back(p4);
            }
            if (it->_pt>30 && TMath::Abs(it->_eta)<4.7)
            {
                TLorentzVector p4;
                p4.SetPtEtaPhiM(it->_pt, it->_eta, it->_phi, it->_mass);
                p40bjets.push_back(p4);
            }
        }

        //
        //  We basically select this event and categorize it now =>
        //  this category is noCats
        //  Start categorization
        //
        FILL_NOJETS(setNoCats);

        if (p4bjets.size()>=1)
        {
            FILL_NOJETS(set2Mu1bjets);
            if (numExtraLeptons==0)
            {
                FILL_NOJETS(set2Mu1bjets0l);
            }
            else if (numExtraLeptons==2)
            {
                FILL_NOJETS(set2Mu1bjets2l);
                if (numElectrons==2 && numExtraMuons==0)
                {
                    if (electrons->at(0)._pt>__cutElePt &&
                        electrons->at(1)._pt>__cutElePt)
                    {
                        FILL_NOJETS(set2Mu1bjets2l2e);
                    }
                }
                else if (numElectrons==1 && numExtraMuons==1)
                {
                    if (electrons->at(0)._pt>__cutElePt &&
                        muons_extra[0]._pt>__cutextramuonPt &&
                        muons_extra[0]._trackIsoSumPt/muons_extra[0]._pt<__cutextramuonIso)
                    {
                        FILL_NOJETS(set2Mu1bjets2l1e1mu);
                    }
                }
                else if (numExtraMuons==2 && numElectrons==0)
                {
                    if (muons_extra[0]._pt>__cutextramuonPt &&
                        muons_extra[1]._pt>__cutextramuonPt &&
                        muons_extra[0]._trackIsoSumPt/muons_extra[0]._pt<__cutextramuonIso
                        &&
                        muons_extra[1]._trackIsoSumPt/muons_extra[1]._pt<__cutextramuonIso)
                    {
                        FILL_NOJETS(set2Mu1bjets2l2mu);
                    }
                }
                else
                {
                    FILL_NOJETS(set2Mu1bjets2lMissed);
                }
            }
            else 
            {
                FILL_NOJETS(set2Mu1bjetsMissed);
            }
        }
        else // no bjets
        {
            FILL_NOJETS(set2Mu0bjets);
            if (numExtraLeptons==0)
            {
                FILL_NOJETS(set2Mu0bjets0l);
                if (p40bjets.size()>=2)
                {
                    //
                    //  2 and more jets
                    //
                    TLorentzVector p4lead = p40bjets[0];
                    TLorentzVector p4sub = p40bjets[1];
                    TLorentzVector p4dijet = p4lead + p4sub;
                    float jetsdeta = getdeta(p4lead, p4sub);
                    float dijetmass = p4dijet.M();
                    float dijetdimuondeta = getdeta(p4dijet, p4dimuon);
                    FILL_JETS(set2Mu0bjets0l2jets);
                    if (dijetmass>__cutdijetMass_VBFTight &&
                        jetsdeta>__cutjetsdeta_VBFTight)
                    {
                        FILL_JETS(set2Mu0bjets0l2jetsVBFTight);
                    }
                    else if (dijetmass>__cutdijetMass_VBFLoose &&
                        jetsdeta>__cutjetsdeta_VBFLoose)
                    {
                        FILL_JETS(set2Mu0bjets0l2jetsVBFLoose);
                    }
                    else if (dijetmass>__cutlowdijetMassVhadH &&
                        dijetmass<__cuthighdijetMassVhadH &&
                        dijetdimuondeta<1.5)
                    {
                        FILL_JETS(set2Mu0bjets0l2jetsVhadH);
                    }
                    else
                    {
                        FILL_JETS(set2Mu0bjets0l2jetsggF);
                    }
                }
                else // 0 leptons 0/1 jets
                {
                    FILL_NOJETS(set2Mu0bjets0l01jets)
                    if (met->_pt>__cutMETPtZinvH)
                    {
                        FILL_NOJETS(set2Mu0bjets0l01jetsZinvH);
                    }
                    else if (met->_pt<=__cutMETPtZinvH && 
                        p4dimuon.Pt()>=__cutdimuonPtggFTight)
                    {
                        FILL_NOJETS(set2Mu0bjets0l01jetsggFTight);
                        if (isBarrel(mu1) && isBarrel(mu2))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFTightBB);
                        }
                        else if ((isBarrel(mu1) && isOverlap(mu2)) ||
                            (isBarrel(mu2) && isOverlap(mu1)))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFTightBO);
                        }
                        else if ((isBarrel(mu1) & isEndcap(mu2)) ||
                            (isBarrel(mu2) && isEndcap(mu1)))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFTightBE);
                        }
                        else if (isOverlap(mu1) && isOverlap(mu2))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFTightOO);
                        }
                        else if ((isOverlap(mu1) && isEndcap(mu2)) ||
                            (isOverlap(mu2) && isEndcap(mu1)))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFTightOE);
                        }
                        else if (isEndcap(mu1) && isEndcap(mu2))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFTightEE);
                        }
                    }
                    else // 01jets ggFLoose
                    {
                        FILL_NOJETS(set2Mu0bjets0l01jetsggFLoose);
                        if (isBarrel(mu1) && isBarrel(mu2))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFLooseBB);
                        }
                        else if ((isBarrel(mu1) && isOverlap(mu2)) ||
                            (isBarrel(mu2) && isOverlap(mu1)))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFLooseBO);
                        }
                        else if ((isBarrel(mu1) & isEndcap(mu2)) ||
                            (isBarrel(mu2) && isEndcap(mu1)))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFLooseBE);
                        }
                        else if (isOverlap(mu1) && isOverlap(mu2))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFLooseOO);
                        }
                        else if ((isOverlap(mu1) && isEndcap(mu2)) ||
                            (isOverlap(mu2) && isEndcap(mu1)))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFLooseOE);
                        }
                        else if (isEndcap(mu1) && isEndcap(mu2))
                        {
                            FILL_NOJETS(set2Mu0bjets0l01jetsggFLooseEE);
                        }
                    }
                }
            }
            else if ((numExtraMuons+numElectrons)==1 ||
                (numExtraMuons+numElectrons)==2)
            {
                FILL_NOJETS(set2Mu0bjets12l);
                if (numElectrons==1 && numExtraMuons==0)
                {
                    if (electrons->at(0)._pt>__cutElePt)
                    {
                        FILL_NOJETS(set2Mu0bjets12l1e);
                    }
                }
                else if (numElectrons==0 && numExtraMuons==1)
                {
                    if (muons_extra[0]._pt>__cutextramuonPt &&
                        muons_extra[0]._trackIsoSumPt/muons_extra[0]._pt<__cutextramuonIso)
                    {
                        FILL_NOJETS(set2Mu0bjets12l1mu);
                    }
                }
                else if (numElectrons==1 && numExtraMuons==1)
                {
                    if (electrons->at(0)._pt>__cutElePt &&
                        muons_extra[0]._pt>__cutextramuonPt &&
                        muons_extra[0]._trackIsoSumPt/muons_extra[0]._pt<__cutextramuonIso)
                    {
                        FILL_NOJETS(set2Mu0bjets12l1mu1e);
                    }
                }
                else if (numElectrons==2 && numExtraMuons==0)
                {
                    if (electrons->at(0)._pt>__cutElePt &&
                        electrons->at(1)._pt>__cutElePt)
                    {
                        FILL_NOJETS(set2Mu0bjets12l2e);
                    }
                }
                else if (numElectrons==0 && numExtraMuons==2)
                {
                    if (muons_extra[0]._pt>__cutextramuonPt &&
                        muons_extra[1]._pt>__cutextramuonPt &&
                        muons_extra[0]._trackIsoSumPt/muons_extra[0]._pt<__cutextramuonIso
                        &&
                        muons_extra[1]._trackIsoSumPt/muons_extra[1]._pt<__cutextramuonIso)
                    {
                        FILL_NOJETS(set2Mu0bjets12l2mu);
                    }
                }
                else continue;
            }
        }
	}

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
    std::cout << "Cuts:" << std::endl
        << "__cutmuonHLTMatchedPt = " << __cutmuonHLTMatchedPt << std::endl
        << "__cutmuonHLTMatchedEta = " << __cutmuonHLTMatchedEta << std::endl
        << "__cutleadmuonPt = " << __cutleadmuonPt << std::endl
        << "__cutsubleadmuonPt = " << __cutsubleadmuonPt << std::endl
        << "__cutleadmuonEta = " << __cutleadmuonEta << std::endl
        << "__cutsubleadmuonEta = " << __cutsubleadmuonEta << std::endl
        << "__cutleadmuonIso = " << __cutleadmuonIso << std::endl
        << "__cutsubleadmuonIso = " << __cutsubleadmuonIso << std::endl
        << "__cutextramuonPt = " <<__cutextramuonPt  << std::endl
        << "__cutextramuonEta = " << __cutextramuonEta << std::endl
        << "__cutextramuonIso = " << __cutextramuonIso << std::endl
        << "__cutleadJetPt = " << __cutleadJetPt << std::endl
        << "__cutsubleadJetPt = " << __cutsubleadJetPt << std::endl
        << "__cutdijetMass_VBFTight = " << __cutdijetMass_VBFTight << std::endl
        << "__cutjetBTagMediumDiscr = " << __cutjetBTagMediumDiscr << std::endl
        << "__cutjetsdeta_VBFTight = " << __cutjetsdeta_VBFTight << std::endl
        << "__cutdijetMass_VBFLoose = " << __cutdijetMass_VBFLoose << std::endl
        << "__cutjetsdeta_VBFLoose = " << __cutjetsdeta_VBFLoose << std::endl
        << "__cutlowdijetMassVhadH = " << __cutlowdijetMassVhadH << std::endl
        << "__cuthighdijetMassVhadH = " << __cuthighdijetMassVhadH << std::endl
        << "__cutMETPtZinvH = " << __cutMETPtZinvH << std::endl
        << "__cutdimuonPtggFTight = " << __cutdimuonPtggFTight << std::endl
        << "__cutElePt = " << __cutElePt << std::endl
        << "__cutTauPt = " << __cutTauPt << std::endl
        << "__cutmaxExtraLeptons = " << __cutmaxExtraLeptons << std::endl;
}
int main(int argc, char** argv)
{
	/*
	 *	Register signals
	 */
	signal(SIGABRT, &sigHandler);
	signal(SIGTERM, &sigHandler);
	signal(SIGINT, &sigHandler);

	std::string none;
	bool genPUMC = false;

	/*
	 *	Pare Options
	 */
	po::options_description desc("Allowed Program Options");
	desc.add_options()
		("help", "produce help messages")
		("input", po::value<std::string>(), "a file specifying all the ROOT files to process")
		("isMC", po::value<bool>(), "type of data: DATA vs MC")
		("output", po::value<std::string>(), "output file name")
		("genPUMC", po::value<bool>(&genPUMC)->default_value(false), "true if should generate the MC PU file")
		("puMC", po::value<std::string>(&none)->default_value("None"), "MC PU Reweight file")
		("puDATA", po::value<std::string>(&none)->default_value("None"), "DATA PU Reweight file")
        ("cutmuonHLTMatchedPt", po::value<double>(&__cutmuonHLTMatchedPt)->default_value(__cutmuonHLTMatchedPt), "Muon Matched Pt Cut")
        ("cutmuonHLTMatchedEta", po::value<double>(&__cutmuonHLTMatchedEta)->default_value(__cutmuonHLTMatchedEta), "Muon Matched Eta Cut")
        ("cutleadmuonPt", po::value<double>(&__cutleadmuonPt)->default_value(__cutleadmuonPt), "Muon Pt Cut")
        ("cutsubleadmuonPt", po::value<double>(&__cutsubleadmuonPt)->default_value(__cutsubleadmuonPt), "Muon Pt Cut")
        ("cutleadmuonEta", po::value<double>(&__cutleadmuonEta)->default_value(__cutleadmuonEta), "Muon Eta Cut")
        ("cutsubleadmuonEta", po::value<double>(&__cutsubleadmuonEta)->default_value(__cutsubleadmuonEta), "Muon Eta Cut")
        ("cutleadmuonIso", po::value<double>(&__cutleadmuonIso)->default_value(__cutleadmuonIso), "Muon Isolation Cut")
        ("cutsubleadmuonIso", po::value<double>(&__cutsubleadmuonIso)->default_value(__cutsubleadmuonIso), "Muon Isolation Cut")

        ("cutextramuonPt", po::value<double>(&__cutextramuonPt)->default_value(__cutextramuonPt), "Pt Cut on Extra Muons")
        ("cutextramuonEta", po::value<double>(&__cutextramuonEta)->default_value(__cutextramuonEta), "Eta cut on Extra Muons")
        ("cutextramuonIso", po::value<double>(&__cutextramuonIso)->default_value(__cutextramuonIso), "Isolation cut on extra muons")

        ("cutleadJetPt", po::value<double>(&__cutleadJetPt)->default_value(__cutleadJetPt), "Lead Jet Pt Cut")
        ("cutsubleadJetPt", po::value<double>(&__cutsubleadJetPt)->default_value(__cutsubleadJetPt), "SubLeading Jet Pt Cut")

        ("cutdijetMass_VBFTight", po::value<double>(&__cutdijetMass_VBFTight)->default_value(__cutdijetMass_VBFTight), "DiJet Mass VBFTight-Category Cut")
        ("cutjetBTagMediumDiscr", po::value<double>(&__cutjetBTagMediumDiscr)->default_value(__cutjetBTagMediumDiscr), "Btag discriminator")

        ("cutjetsdeta_VBFTight", po::value<double>(&__cutjetsdeta_VBFTight)->default_value(__cutjetsdeta_VBFTight), "deta between Jets Cut VBFTight")
        ("cutdijetMass_VBFLoose", po::value<double>(&__cutdijetMass_VBFLoose)->default_value(__cutdijetMass_VBFLoose), "dijet Mass VBFLoose")
        ("cutjetsdeta_VBFLoose", po::value<double>(&__cutjetsdeta_VBFLoose)->default_value(__cutjetsdeta_VBFLoose), "deta between Jets VBFLoose")
        ("cutlowdijetMassVhadH", po::value<double>(&__cutlowdijetMassVhadH)->default_value(__cutlowdijetMassVhadH), "Low dijet Mass VHad H")
        ("cuthighdijetMassVhadH", po::value<double>(&__cuthighdijetMassVhadH)->default_value(__cuthighdijetMassVhadH), "High dijet Mass VHad H")
        ("cutMETPtZinvH", po::value<double>(&__cutMETPtZinvH)->default_value(__cutMETPtZinvH), "MET Pt Zinv H")
        ("cutdimuonPtggFTight", po::value<double>(&__cutdimuonPtggFTight)->default_value(__cutdimuonPtggFTight), "dimuon Pt for ggFTight")
        ("cutElePt", po::value<double>(&__cutElePt)->default_value(__cutElePt), "electron Pt")
        ("cutTauPt", po::value<double>(&__cutTauPt)->default_value(__cutTauPt),
         "Tau Pt")
        ("cutmaxExtraLeptons", po::value<double>(&__cutmaxExtraLeptons)->default_value(__cutmaxExtraLeptons), "Max Number of Extra Leptons")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || argc<2)
	{
		std::cout << desc << std::endl;
		return 1;
	}

	//	Assign globals
	__inputfilename = vm["input"].as<std::string>();
	__isMC = vm["isMC"].as<bool>();
	__outputfilename = vm["output"].as<std::string>();
	__genPUMC = vm["genPUMC"].as<bool>();
	__puMCfilename = vm["puMC"].as<std::string>();
	__puDATAfilename = vm["puDATA"].as<std::string>();

    __cutmuonHLTMatchedPt = vm["cutmuonHLTMatchedPt"].as<double>();
    __cutmuonHLTMatchedEta = vm["cutmuonHLTMatchedEta"].as<double>();

    __cutleadmuonPt = vm["cutleadmuonPt"].as<double>();
    __cutsubleadmuonPt = vm["cutsubleadmuonPt"].as<double>();
    __cutleadmuonEta = vm["cutleadmuonEta"].as<double>();
    __cutsubleadmuonEta = vm["cutsubleadmuonEta"].as<double>();
    __cutleadmuonIso = vm["cutleadmuonIso"].as<double>();
    __cutsubleadmuonIso = vm["cutsubleadmuonIso"].as<double>();

    __cutextramuonPt = vm["cutextramuonPt"].as<double>();
    __cutextramuonEta = vm["cutextramuonEta"].as<double>();
    __cutextramuonIso = vm["cutextramuonIso"].as<double>();

    __cutleadJetPt = vm["cutleadJetPt"].as<double>();
    __cutsubleadJetPt = vm["cutsubleadJetPt"].as<double>();
    __cutdijetMass_VBFTight = vm["cutdijetMass_VBFTight"].as<double>();
    __cutjetBTagMediumDiscr = vm["cutjetBTagMediumDiscr"].as<double>();
    __cutjetsdeta_VBFTight = vm["cutjetsdeta_VBFTight"].as<double>();
    __cutdijetMass_VBFLoose = vm["cutdijetMass_VBFLoose"].as<double>();
    __cutjetsdeta_VBFLoose = vm["cutjetsdeta_VBFLoose"].as<double>();
    __cutlowdijetMassVhadH = vm["cutlowdijetMassVhadH"].as<double>();
    __cuthighdijetMassVhadH = vm["cuthighdijetMassVhadH"].as<double>();
    __cutMETPtZinvH = vm["cutMETPtZinvH"].as<double>();
    __cutdimuonPtggFTight = vm["cutdimuonPtggFTight"].as<double>();

    __cutElePt = vm["cutElePt"].as<double>();
    __cutTauPt = vm["cutTauPt"].as<double>();
    __cutmaxExtraLeptons = vm["cutmaxExtraLeptons"].as<double>();

    printCuts();

	//	start processing
	process();
	return 0;
}

#endif
