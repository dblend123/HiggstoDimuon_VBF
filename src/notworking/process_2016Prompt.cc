
//	
#ifdef STANDALONE
#include "Muon.h"
#include "Jet.h"
#include "Vertex.h"
#include "Event.h"
#include "MET.h"
#include "Constants.h"
#include "Streamer.h"
#include "MetaHiggs.h"
#include "HistogramSets.h"

//	ROOT headers
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1D.h"

std::string const NTUPLEMAKER_NAME =  "ntuplemaker_H2DiMuonMaker";
std::string const CATEGORY_NAMES[5] = {
	"VBFTight", "ggFTight", "VBFLoose", "01JetTight", "01JetLoose" } ;

#define VBFTight 0
#define ggFTight 1
#define VBFLoose 2
#define JET01Tight 3
#define JET01Loose 4
#define NCATEGORIES 5

std::map<std::string, std::vector<std::pair<int, long long int> > > mapcats;
std::vector<std::pair<int, long long int> >  categories[NCATEGORIES];

using namespace analysis::core;
using namespace analysis::dimuon;
using namespace analysis::processing;

DimuonSet setNoCats("_NoCats");
DimuonSet set2Jets("_2Jets");
DimuonSet set01Jets("_01Jets");

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

bool passMuon(Muon const& m)
{
	float muIso = (m._sumChargedHadronPtR04 + max(0., 
		m._sumNeutralHadronEtR04+m._sumPhotonEtR04-0.5*m._sumPUPtR04))/m._pt;
	if (m._isGlobal && m._isTracker &&
		m._pt>10 && TMath::Abs(m._eta)<2.4 &&
		m._isTight && (m._trackIsoSumPt/m._pt)<0.1)
		return true;
	return false;
}

bool passMuonHLT(Muon const& m)
{
	if ((m._isHLTMatched[1] || m._isHLTMatched[0]) &&
		m._pt>20 && TMath::Abs(m._eta)<2.4)
		return true;

	return false;
}

bool passMuons(Muon const& m1, Muon const& m2)
{
	if (m1._charge!=m2._charge &&
		passMuon(m1) && passMuon(m2))
		if (passMuonHLT(m1) || passMuonHLT(m2))
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

void categorize(Jets* jets, Muon const& mu1, Muon const&  mu2, 
	MET const& met, Event const& event)
{
	int run = event._run;
	long long e = event._event;
	std::pair<int, long long int> runevent(run, e);

	TLorentzVector p4m1, p4m2;
	p4m1.SetPtEtaPhiM(mu1._pt, mu1._eta, 
		mu1._phi, PDG_MASS_Mu);
	p4m2.SetPtEtaPhiM(mu2._pt, mu2._eta, 
		mu2._phi, PDG_MASS_Mu);
	TLorentzVector p4dimuon = p4m1 + p4m2;
	if (!(p4dimuon.M()>110 && p4dimuon.M()<160 &&
		mu1._isPF && mu2._isPF))
		return;

	//	Fill the No Categorization Set
	double dphi = p4m1.DeltaPhi(p4m2);
	setNoCats.hDiMuonpt->Fill(p4dimuon.Pt());
	setNoCats.hDiMuonMass->Fill(p4dimuon.M());
	setNoCats.hDiMuoneta->Fill(p4dimuon.Eta());
	setNoCats.hDiMuondphi->Fill(dphi);
	setNoCats.hMuonpt->Fill(p4m1.Pt());
	setNoCats.hMuonpt->Fill(p4m2.Pt());
	setNoCats.hMuoneta->Fill(p4m1.Eta());
	setNoCats.hMuoneta->Fill(p4m2.Eta());
	setNoCats.hMuonphi->Fill(p4m1.Phi());
	setNoCats.hMuonphi->Fill(p4m2.Phi());

	std::vector<TLorentzVector> p4jets;
	for (Jets::const_iterator it=jets->begin(); it!=jets->end(); ++it)
	{
		if (it->_pt>30 && TMath::Abs(it->_eta)<4.7)
		{
			if (!(jetMuondR(it->_eta, it->_phi, mu1._eta, mu1._phi)<0.3) && 
				!(jetMuondR(it->_eta, it->_phi, mu2._eta, mu2._phi)<0.3))
			{
				TLorentzVector p4;
				p4.SetPtEtaPhiM(it->_pt, it->_eta, it->_phi, it->_mass);
				p4jets.push_back(p4);
			}
		}
	}

	if (p4jets.size()>2)
		return;

	bool isPreSelected = false;
	if (p4jets.size()==2)
	{
		TLorentzVector p4lead = p4jets[0]; 
		TLorentzVector p4sub = p4jets[1];
		TLorentzVector dijet = p4lead + p4sub;
		float deta = p4lead.Eta() - p4sub.Eta();
		float dijetmass = dijet.M();
		
		if (p4lead.Pt()>40 && p4sub.Pt()>30 &&
			met._pt<40)
		{
			isPreSelected = true;
	
			set2Jets.hDiJetMass->Fill(dijetmass);
			set2Jets.hDiJetdeta->Fill(TMath::Abs(deta));
			set2Jets.hDiMuonpt->Fill(p4dimuon.Pt());
			set2Jets.hDiMuonMass->Fill(p4dimuon.M());
			set2Jets.hDiMuoneta->Fill(p4dimuon.Eta());
			set2Jets.hDiMuondphi->Fill(dphi);
			set2Jets.hMuonpt->Fill(p4m1.Pt());
			set2Jets.hMuonpt->Fill(p4m2.Pt());
			set2Jets.hMuoneta->Fill(p4m1.Eta());
			set2Jets.hMuoneta->Fill(p4m2.Eta());
			set2Jets.hMuonphi->Fill(p4m1.Phi());
			set2Jets.hMuonphi->Fill(p4m2.Phi());

			if (dijetmass>650 && TMath::Abs(deta)>3.5)
			{
				categories[VBFTight].push_back(runevent);
				mapcats[CATEGORY_NAMES[VBFTight]].push_back(runevent);
				return;}
			if (dijetmass>250 && p4dimuon.Pt()>50)
			{
				categories[ggFTight].push_back(runevent);
				mapcats[CATEGORY_NAMES[ggFTight]].push_back(runevent);
				return;}
			else
			{
				categories[VBFLoose].push_back(runevent);
				mapcats[CATEGORY_NAMES[VBFLoose]].push_back(runevent);
				return;}
		}
	}
	if (!isPreSelected)
	{
		set01Jets.hDiMuonpt->Fill(p4dimuon.Pt());
		set01Jets.hDiMuonMass->Fill(p4dimuon.M());
		set01Jets.hDiMuoneta->Fill(p4dimuon.Eta());
		set01Jets.hDiMuondphi->Fill(dphi);
		set01Jets.hMuonpt->Fill(p4m1.Pt());
		set01Jets.hMuonpt->Fill(p4m2.Pt());
		set01Jets.hMuoneta->Fill(p4m1.Eta());
		set01Jets.hMuoneta->Fill(p4m2.Eta());
		set01Jets.hMuonphi->Fill(p4m1.Phi());
		set01Jets.hMuonphi->Fill(p4m2.Phi());
		if (p4dimuon.Pt()>=10)
		{
			categories[JET01Tight].push_back(runevent);
			mapcats[CATEGORY_NAMES[JET01Tight]].push_back(runevent);
			return;
		}
		else
		{
			categories[JET01Loose].push_back(runevent);
			mapcats[CATEGORY_NAMES[JET01Loose]].push_back(runevent);
			return;
		}
	}
	
	return;
}

void sampleinfo(std::string const& inputname)
{
	Streamer s(inputname, NTUPLEMAKER_NAME+"/Meta");
	s.chainup();

	using namespace analysis::dimuon;
	MetaHiggs *meta=NULL;
	s._chain->SetBranchAddress("Meta", &meta);

	long long int numEvents = 0;
	for (int i=0; i<s._chain->GetEntries(); i++)
	{
		s._chain->GetEntry(i);
		numEvents += meta->_nEventsProcessed;
	}
	std::cout << "#events processed total = " << numEvents << std::endl;
}

void synchronize(std::string const& inputname)
{
	TFile *outroot = new TFile("process_2016Prompt.root", "recreate");
	setNoCats.init();
	set2Jets.init();
	set01Jets.init();

	Streamer streamer(inputname, NTUPLEMAKER_NAME+"/Events");
	streamer.chainup();

	Muons *muons1=NULL;
	Muons *muons2=NULL;
	Jets *jets=NULL;
	Vertices *vertices=NULL;
	Event *event=NULL;
	EventAuxiliary *aux=NULL;
	MET *met=NULL;
	streamer._chain->SetBranchAddress("Muons1", &muons1);
	streamer._chain->SetBranchAddress("Muons2", &muons2);
	streamer._chain->SetBranchAddress("Jets", &jets);
	streamer._chain->SetBranchAddress("Vertices", &vertices);
	streamer._chain->SetBranchAddress("Event", &event);
	streamer._chain->SetBranchAddress("EventAuxiliary", &aux);
	streamer._chain->SetBranchAddress("MET", &met);

	uint32_t numEntries = streamer._chain->GetEntries();
	for (uint32_t i=0; i<numEntries; i++)
	{
		streamer._chain->GetEntry(i);
		if (i%1000==0)
			std::cout << "### Event " << i << " / " << numEntries
				<< std::endl;

		//
		//	Selections
		//
		if (!passVertex(vertices))
			continue;
		if (!(aux->_hasHLTFired[0] || aux->_hasHLTFired[1]))
			continue;
		for (uint32_t im=0; im<muons1->size(); im++)
		{
			if (!passMuons(muons1->at(im), muons2->at(im)))
				continue;
			categorize(jets, muons1->at(im), muons2->at(im), *met, *event);
		}
	}

	outroot->Write();
	outroot->Close();

	return;
}

int main(int argc, char** argv)
{
	if (argc<2 || argc>2)
	{
		std::cout << "Usage:" << std::endl
			<< "./synchronize <input file name>" << std::endl;
		return 1;
	}

	for (int i=0; i<NCATEGORIES; i++)
	{
		std::cout << CATEGORY_NAMES[i] << std::endl;
		mapcats[CATEGORY_NAMES[i]] = std::vector<std::pair<int, long long int>>();
	}

	sampleinfo(argv[1]);
	synchronize(argv[1]);
	return 0;
}

#endif
