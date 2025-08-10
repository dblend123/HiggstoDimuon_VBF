#pragma once

#include "TString.h"
#include "TH1D.h"
#include "TDirectory.h"

namespace analysis
{
	namespace processing
	{
		struct DimuonSet
		{
			DimuonSet(TString const& postfix) 
			{
				_postfix = postfix;
			}

			void init()
			{
				gDirectory->mkdir(_postfix);
				gDirectory->cd(_postfix);
				hDiJetMass = new TH1D("DiJetMass", "DiJetMass",
					20, 0, 1000);
                hDiJetMass->Sumw2(kTRUE);
				hDiJetdeta = new TH1D("DiJetdeta", "DiJetdeta",
					14, 0, 7);
                hDiJetdeta->Sumw2(kTRUE);
				hDiMuonpt = new TH1D("DiMuonpt", "DiMuonpt", 
					100, 0, 200);
                hDiMuonpt->Sumw2(kTRUE);
				hDiMuonMass = new TH1D("DiMuonMass", 
					"DiMuonMass", 7000, 50, 400);
                hDiMuonMass->Sumw2(kTRUE);
				hDiMuoneta = new TH1D("DiMuoneta", "DiMuoneta",
					50, -2.5, 2.5);
                hDiMuoneta->Sumw2(kTRUE);
				hDiMuondphi = new TH1D("DiMuondphi", 
					"DiMoundphi", 18, -3.6, 3.6);
                hDiMuondphi->Sumw2(kTRUE);
				hMuonpt = new TH1D("Muonpt", "Muonpt", 
					50, 0, 100);
                hMuonpt->Sumw2(kTRUE);
				hMuoneta = new TH1D("Muoneta", "Muoneta", 
					50, -2.5, 2.5);
                hMuoneta->Sumw2(kTRUE);
				hMuonphi = new TH1D("Muonphi", "Muonphi", 
					36, -3.6, 3.6);
                hMuonphi->Sumw2(kTRUE);
				gDirectory->cd("../");
			}

            void init2()                                                      //
            {                                                                 //
                gDirectory->mkdir(_postfix);                                  //
                gDirectory->cd(_postfix);                                     //
                hFsrPhotonpt = new TH1F("FsrPhotonpt", "FsrPhotonpt",         //
                    75, 0, 75);                                               //
                //hFsrPhotonpt->Sumw2(kTRUE);                                   //
                hFsrPhotoneta = new TH1F("FsrPhotoneta", "FsrPhotoneta",      //
                    50, -2.5, 2.5);                                           //
                //hFsrPhotoneta->Sumw2(kTRUE);                                  //
                hFsrPhotonphi = new TH1F("FsrPhotonphi", "FsrPhotonphi",      //
                    36, -3.6, 3.6);                                           //
                //hFsrPhotonphi->Sumw2(kTRUE);                                  //
                gDirectory->cd("../");                                        //
            }                                                                 //

			~DimuonSet() {}

            TString _postfix;
			TH1D *hDiJetMass;
			TH1D *hDiJetdeta;
			TH1D *hDiMuonpt;
			TH1D *hDiMuonMass;
			TH1D *hDiMuoneta;
			TH1D *hDiMuondphi;
			TH1D *hMuonpt;
			TH1D *hMuoneta;
			TH1D *hMuonphi;
			TH1D *hNpv;
            TH1F *hFsrPhotonpt;   //
            TH1F *hFsrPhotonphi;  //
            TH1F *hFsrPhotoneta;  //
		};
	}
}
