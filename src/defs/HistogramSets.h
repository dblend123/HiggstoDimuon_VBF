#pragma once

#include "TString.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TH2F.h"

namespace analysis
{
	namespace processing
	{
		struct DimuonSet
		{
			DimuonSet(TDirectory* parentDir, TString const& postfix)
			    : _postfix(postfix)
            {
			    _dir = parentDir->mkdir(postfix);
                _dir->cd();
				
                hDiMuonpt = new TH1F("DiMuonpt", "DiMuonpt", 120, 0, 120);
                //hDiMuonpt->Sumw2(kTRUE);

				hDiMuonMass = new TH1F("DiMuonMass", "DiMuonMass", 150, 50, 200);
                //hDiMuonMass->Sumw2(kTRUE);
				
                /*hDiMuonMass_pt20  = new TH1F("DiMuonMass_pt20",  "DiMuonMass_pt20",  50, 100, 150);
                hDiMuonMass_pt25  = new TH1F("DiMuonMass_pt25",  "DiMuonMass_pt25",  50, 100, 150);
                hDiMuonMass_pt30  = new TH1F("DiMuonMass_pt30",  "DiMuonMass_pt30",  50, 100, 150);
                hDiMuonMass_pt35  = new TH1F("DiMuonMass_pt35",  "DiMuonMass_pt35",  50, 100, 150);
                hDiMuonMass_pt40  = new TH1F("DiMuonMass_pt40",  "DiMuonMass_pt40",  50, 100, 150);
                hDiMuonMass_pt50  = new TH1F("DiMuonMass_pt50",  "DiMuonMass_pt50",  50, 100, 150);
                hDiMuonMass_pt60  = new TH1F("DiMuonMass_pt60",  "DiMuonMass_pt60",  50, 100, 150);
                hDiMuonMass_pt80  = new TH1F("DiMuonMass_pt80",  "DiMuonMass_pt80",  50, 100, 150);
                hDiMuonMass_pt120 = new TH1F("DiMuonMass_pt120", "DiMuonMass_pt120", 50, 100, 150);
                // Define custom bin edges
                double etaBins[] = {-2.4, -1.6, -0.9, 0.9, 1.6, 2.4};
                int nEtaBins = sizeof(etaBins)/sizeof(double) - 1;

                double ptBins[] = {20, 25, 30, 35, 40, 50, 60, 80, 120};
                int nPtBins = sizeof(ptBins)/sizeof(double) - 1;

                // Create TH2F with variable bins
                hDiMuonPt_vs_Eta = new TH2F("DiMuonPtVsEta", "Dimuon p_{T} vs #eta;#eta^{#mu#mu};p_{T}^{#mu#mu} [GeV]",
                                          nEtaBins, etaBins,
                                          nPtBins, ptBins);*/

				hDiMuonMass_Z = new TH1F("DiMuonMass_Z", "DiMuonMass_Z", 90, 70, 115);
                //hDiMuonMass->Sumw2(kTRUE);
				
                hDiMuonMass_SB_SR1 = new TH1F("DiMuonMass_SB_SR1", "DiMuonMass_SB_SR1", 50, 100, 150);
                //hDiMuonMass_SB_SR->Sumw2(kTRUE);
                
                //hDiMuonMass_SB_SR2 = new TH1F("DiMuonMass_SB_SR2", "DiMuonMass_SB_SR2", 25, 100, 150);
                //hDiMuonMass_SB_SR->Sumw2(kTRUE);
                
                //hDiMuonMass_SB_SR25 = new TH1F("DiMuonMass_SB_SR2.5", "DiMuonMass_SB_SR2.5", 20, 100, 150);
                //hDiMuonMass_SB_SR->Sumw25(kTRUE);
                
                //hDiMuonMass_SB_SR3 = new TH1F("DiMuonMass_SB_SR3", "DiMuonMass_SB_SR3", 17, 100, 150);
                //hDiMuonMass_SB_SR->Sumw3(kTRUE);
                
                //hDiMuonMass_SB_SR35 = new TH1F("DiMuonMass_SB_SR3.5", "DiMuonMass_SB_SR3.5", 14, 100, 150);
                //hDiMuonMass_SB_SR->Sumw35(kTRUE);
                
                hDiMuoneta = new TH1F("DiMuoneta", "DiMuoneta", 50, -2.5, 2.5);
                //hDiMuoneta->Sumw2(kTRUE);
				
                hDiMuoneta_JPsi = new TH1F("DiMuoneta_JPsi", "DiMuoneta_JPsi", 50, -2.5, 2.5);
                //hDiMuoneta_JPsi->Sumw2(kTRUE);
				
                hDiMuoneta_Upsilon = new TH1F("DiMuoneta_Upsilon", "DiMuoneta_Upsilon", 50, -2.5, 2.5);
                //hDiMuoneta_Upsilon->Sumw2(kTRUE);
				
                hDiMuoneta_Z = new TH1F("DiMuoneta_Z", "DiMuoneta_Z", 50, -2.5, 2.5);
                //hDiMuoneta_Z->Sumw2(kTRUE);
				
                hDiMuoneta_Higgs = new TH1F("DiMuoneta_Higgs", "DiMuoneta_Higgs", 50, -2.5, 2.5);
                //hDiMuoneta_Higgs->Sumw2(kTRUE);
				
                hDiMuondphi = new TH1F("DiMuondphi","DiMoundphi", 75, -3.2, 3.2);
                //hDiMuondphi->Sumw2(kTRUE);

				hMuonpt = new TH1F("Muonpt", "Muonpt", 100, 0, 120);
                //hMuonpt->Sumw2(kTRUE);

				//hMuoneta = new TH1F("Muoneta", "Muoneta", 50, -2.5, 2.5);
                //hMuoneta->Sumw2(kTRUE);
				
                //hMuonphi = new TH1F("Muonphi", "Muonphi", 36, -3.6, 3.6);
                //hMuonphi->Sumw2(kTRUE);
                
                hMuonDXY = new TH1F("MuonDXY", "MuonDXY", 500, -.025, .025);
                //hMuonDXY->Sumw2(kTRUE);
                
                hMuonDZ = new TH1F("MuonDZ", "MuonDZ", 500, -.035, .035);
                //hMuonDZ->Sumw2(kTRUE);
                
                hMuonDXY_vs_MuonDZ = new TH2F("MuonDXY_vs_MuonDZ", "MuonDXY_vs_MuonDZ", 125, -.05, .05, 125, -.05, .05);
                //hMuonDXY_vs_DZ->Sumw2(kTRUE);
				
                hDimuoneta_vs_DimuonMass = new TH2F("Dimuoneta_vs_DimuonMass", "Dimuoneta_vs_DimuonMass", 150, 0., 150., 50, -8., 8.);  
                //hDimuonpt_vs_DimuonMass->Sumw2(kTRUE);
               
                hDimuonpt_vs_DimuonMass = new TH2F("Dimuonpt_vs_DimuonMass", "Dimuonpt_vs_DimuonMass", 150, 0, 150, 120, 0, 120);  
                //hDimuonpt_vs_DimuonMass->Sumw2(kTRUE);
               
                hDimuonpt_vs_DimuonMass_JPsi = new TH2F("Dimuonpt_vs_DimuonMass_JPsi", "Dimuonpt_vs_DimuonMass_JPsi", 125, 2.6, 3.6, 150, 50, 125);  
                //hDimuonpt_vs_DimuonMass_JPsi->Sumw2(kTRUE);

                hDimuonpt_vs_DimuonMass_Upsilon = new TH2F("Dimuonpt_vs_DimuonMass_Upsilon", "Dimuonpt_vs_DimuonMass_Upsilon", 60, 8.5, 10.5, 75, 50, 125);  
                //hDimuonpt_vs_DimuonMass_Upsilon->Sumw2(kTRUE);

                hDimuonpt_vs_DimuonMass_Z = new TH2F("Dimuonpt_vs_DimuonMass_Z", "Dimuonpt_vs_DimuonMass_Z", 100, 82, 98, 120, 0, 120);  
                //hDimuonpt_vs_DimuonMass_Z->Sumw2(kTRUE);

                hDimuonpt_vs_DimuonMass_Higgs = new TH2F("Dimuonpt_vs_DimuonMass_Higgs", "Dimuonpt_vs_DimuonMass_Higgs", 100, 120, 130, 120, 0, 120);  
                //hDimuonpt_vs_DimuonMass_Higgs->Sumw2(kTRUE);
               
                

                parentDir->cd();
			}
                ~DimuonSet() {}
                
                TDirectory* _dir;
                TString _postfix;
                
                TH1F *hDiJetMass;
                TH1F *hDiJetdeta;
                TH1F *hDiMuonpt;
                TH1F *hDiMuonMass;
                
                /*TH1F *hDiMuonMass_pt20;
                TH1F *hDiMuonMass_pt25;
                TH1F *hDiMuonMass_pt30;
                TH1F *hDiMuonMass_pt35;
                TH1F *hDiMuonMass_pt40;
                TH1F *hDiMuonMass_pt50;
                TH1F *hDiMuonMass_pt60;
                TH1F *hDiMuonMass_pt80;
                TH1F *hDiMuonMass_pt120;
                TH2F *hDiMuonPt_vs_Eta;*/

                TH1F *hDiMuonMass_Z;
                TH1F *hDiMuonMass_SB_SR1;
                /*TH1F *hDiMuonMass_SB_SR2;
                TH1F *hDiMuonMass_SB_SR25;
                TH1F *hDiMuonMass_SB_SR3;
                TH1F *hDiMuonMass_SB_SR35;*/
                TH1F *hDiMuoneta;
                TH1F *hDiMuoneta_JPsi;
                TH1F *hDiMuoneta_Upsilon;
                TH1F *hDiMuoneta_Z;
                TH1F *hDiMuoneta_Higgs;
                TH1F *hDiMuondphi;
                TH1F *hMuonpt;
                //TH1F *hMuoneta;
                //TH1F *hMuonphi;
                TH1F *hNpv;
                TH1F *hMuonDXY;
                TH1F *hMuonDZ;
                TH2F *hMuonDXY_vs_MuonDZ;
                TH2F *hDimuoneta_vs_DimuonMass;
                TH2F *hDimuonpt_vs_DimuonMass;
                TH2F *hDimuonpt_vs_DimuonMass_JPsi;
                TH2F *hDimuonpt_vs_DimuonMass_Upsilon;
                TH2F *hDimuonpt_vs_DimuonMass_Z;
                TH2F *hDimuonpt_vs_DimuonMass_Higgs;
        };
       
        struct JetSet
		{
			JetSet(TDirectory* parentDir, const TString& postfix) 
			    : _postfix(postfix)
            {
			    _dir = parentDir->mkdir(_postfix);
                _dir->cd();
                
                hDiJetpt = new TH1F("DiJetPt", "DiJetPt", 50, 0, 400);
                //hDiJetMass->Sumw2(kTRUE);
				
                hDiJetdeta = new TH1F("DiJetdeta", "DiJetdeta", 40, 0, 8);
                //hDiJetdeta->Sumw2(kTRUE);
                
                hDiJetMass = new TH1F("DiJetMass", "DiJetMass", 50, 0, 1000);
                //hDiJetMass->Sumw2(kTRUE);
                
                parentDir->cd();
            }

			~JetSet() {}
                
            TDirectory* _dir;
            TString _postfix;

            TH1F *hDiJetpt;   
            TH1F *hDiJetdeta;
            TH1F *hDiJetMass;
		};        

        struct FsrSet
		{
			FsrSet(TString const& postfix) 
			{
				_postfix = postfix;
			}

			void init()
            {                                                                 
                gDirectory->mkdir(_postfix);                                  
                gDirectory->cd(_postfix);                                     
                
                hFsrPhotonpt = new TH1F("FsrPhotonpt", "FsrPhotonpt", 75, 0, 75);
                //hFsrPhotonpt->Sumw2(kTRUE);
                
                hFsrPhotoneta = new TH1F("FsrPhotoneta", "FsrPhotoneta",  50, -2.5, 2.5);
                //hFsrPhotoneta->Sumw2(kTRUE);
                
                hFsrPhotonphi = new TH1F("FsrPhotonphi", "FsrPhotonphi", 36, -3.6, 3.6);
                //hFsrPhotonphi->Sumw2(kTRUE);
               
                //X-axis first, then Y-axis second 
                hFsrPhotonpt_vs_FsrPhotoneta = new TH2F("FsrPhotonpt_vs_FsrPhotoneta", "FsrPhotonpt_vs_FsrPhotoneta", 65, -2.5, 2.5, 80, 0, 75);  
                //hFsrPhotonpt_vs_FsrPhotoneta->Sumw2(kTRUE);

                hFsrPhotonpt_vs_FsrPhotonphi = new TH2F("FsrPhotonpt_vs_FsrPhotonphi", "FsrPhotonpt_vs_FsrPhotonphi", 65, -3.6, 3.6, 80, 0, 75);  
                //hFsrPhotonpt_vs_FsrPhotonphi->Sumw2(kTRUE);
                
                hFsrPhoton_Bestpt = new TH1F("FsrPhoton_Bestpt", "FsrPhoton_Bestpt", 75, 0, 75);
                //hFsrPhoton_Bestpt->Sumw2(kTRUE);
                
                hFsrPhoton_Besteta = new TH1F("FsrPhoton_Besteta", "FsrPhoton_Besteta", 60, -2.5, 2.5);
                //hFsrPhoton_Besteta->Sumw2(kTRUE);
                
                hFsrPhoton_Bestphi = new TH1F("FsrPhoton_Bestphi", "FsrPhoton_Bestphi", 60, -3.6, 3.6);
                //hFsrPhoton_Bestphi->Sumw2(kTRUE); 

                //X-axis first, then Y-axis second 
                hFsrPhoton_Bestpt_vs_FsrPhoton_Besteta = new TH2F("FsrPhoton_Bestpt_vs_FsrPhoton_Besteta", "FsrPhoton_Bestpt_vs_FsrPhoton_Besteta", 65, -2.5, 2.5, 80, 0, 75);  
                //hFsrPhoton_Bestpt_vs_FsrPhoton_Besteta->Sumw2(kTRUE);

                hFsrPhoton_Bestpt_vs_FsrPhoton_Bestphi = new TH2F("FsrPhoton_Bestpt_vs_FsrPhoton_Bestphi", "FsrPhoton_Bestpt_vs_FsrPhoton_Bestphi", 65, -3.6, 3.6, 80, 0, 75);  
                //hFsrPhoton_Bestpt_vs_FsrPhoton_Bestphi->Sumw2(kTRUE);

                hFsrPhoton_Dimuonpt = new TH1F("FsrPhoton_Dimuonpt", "FsrPhoton_Dimuonpt", 100, 0, 200);
                //hFsrPhoton_Bestphi->Sumw2(kTRUE);
               
                hFsrPhoton_DimuonMass = new TH1F("FsrPhoton_DimuonMass", "FsrPhoton_DimuonMass", 300, 50, 200);
                //hFsrPhoton_Bestphi->Sumw2(kTRUE);
               
                hFsrPhoton_Dimuoneta = new TH1F("FsrPhoton_Dimuoneta", "FsrPhoton_Dimuoneta", 60, -2.5, 2.5);
                //hFsrPhoton_Dimuoneta->Sumw2(kTRUE);

                hFsrPhoton_Dimuoneta_Z = new TH1F("FsrPhoton_Dimuoneta_Z", "FsrPhoton_Dimuoneta_Z", 60, -2.5, 2.5);
                //hFsrPhoton_Dimuoneta_Z->Sumw2(kTRUE);

                hFsrPhoton_Dimuonphi = new TH1F("FsrPhoton_Dimuonphi", "FsrPhoton_Dimuonphi", 36, -3.6, 3.6);
                //hFsrPhoton_Dimuoniphi->Sumw2(kTRUE);

                //X-axis first, then Y-axis second RENAME TO SOMETHING 
                hFsrPhotonDimuonpt_vs_FsrPhotonDimoneta = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimoneta", "FsrPhotonDimuonpt_vs_FsrPhotonDimoneta", 50, -2.5, 2.5, 75, 0, 75);  
                //hFsrPhotonDimuonpt_vs_FsrPhotonDimoneta->Sumw2(kTRUE);

                hFsrPhotonDimuonpt_vs_FsrPhotonDimuonphi = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimuonphi", "FsrPhotonDimuonpt_vs_FsrPhotonDimuonphi", 50, -3.6, 3.6, 75, 0, 75);  
                //hFsrPhotonDimuonpt_vs_FsrPhotonDimuonphi->Sumw2(kTRUE);

                hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass", "FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass", 200, 0, 200, 200, 0, 200);  
                //hDimuonpt_vs_DimuonMass->Sumw2(kTRUE);
				
                hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_JPsi = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_JPsi", "FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_JPsi", 100, 2.6, 3.6, 150, 50, 125);
                //hDimuonpt_vs_DimuonMass_JPsi->Sumw2(kTRUE);
				
                hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Upsilon = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Upsilon", "FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Upsilon", 40, 8.5, 10.5, 75, 50, 125);
                //hDimuonpt_vs_DimuonMass_Upsilon->Sumw2(kTRUE);
				
                hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Z = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Z", "FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Z", 112, 82, 98, 200, 0, 160);  
                //hDimuonpt_vs_DimuonMass_Z->Sumw2(kTRUE);
				
                hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Higgs = new TH2F("FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Higgs", "FsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Higgs", 50, 120, 130, 150, 0, 150);  
                //hDimuonpt_vs_DimuonMass_Higgs->Sumw2(kTRUE);
				
                hFsrPhotonMass = new TH1F("FsrPhotonMass", "FsrPhotonMass", 300, 50, 400);
                //hFsrPhotonMass->Sumw2(kTRUE);
				
				hFsrPhotonMass_Higgs = new TH1F("FsrPhotonMass_Higgs", "FsrPhotonMass_Higgs", 11, 110, 150);
                //hFsrPhotonMass_Higgs->Sumw2(kTRUE);
                
                gDirectory->cd("../");                                        
            }                                                                 

			~FsrSet() {}

            TString _postfix;
            TH1F *hFsrPhotonpt;   
            TH1F *hFsrPhotonphi;  
            TH1F *hFsrPhotoneta;
            TH2F *hFsrPhotonpt_vs_FsrPhotoneta;
            TH2F *hFsrPhotonpt_vs_FsrPhotonphi;
            TH1F *hFsrPhoton_Bestpt;
            TH1F *hFsrPhoton_Besteta;
            TH1F *hFsrPhoton_Bestphi;
            TH2F *hFsrPhoton_Bestpt_vs_FsrPhoton_Besteta;
            TH2F *hFsrPhoton_Bestpt_vs_FsrPhoton_Bestphi;
            TH1F *hFsrPhoton_Dimuonpt;
            TH1F *hFsrPhoton_DimuonMass;
            TH1F *hFsrPhoton_Dimuoneta;
            TH1F *hFsrPhoton_Dimuoneta_Z;
            TH1F *hFsrPhoton_Dimuonphi;
            TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimoneta;
            TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimuonphi;
		    TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass;
		    TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_JPsi;
		    TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Upsilon;
		    TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Z;
		    TH2F *hFsrPhotonDimuonpt_vs_FsrPhotonDimuonMass_Higgs;
            TH1F *hFsrPhotonMass;
            TH1F *hFsrPhotonMass_Higgs;
        };
        
        struct VertexSet
		{
			VertexSet(TString const& postfix) 
			{
				_postfix = postfix;
			}

			void init()
            {
                gDirectory->mkdir(_postfix);
                gDirectory->cd(_postfix);
                
                hVertexX = new TH1F("VertexX", "VertexX", 68, .16, .185);
                //hVertexX->Sumw2(kTRUE);
                
                hVertexY = new TH1F("VertexY", "VertexY", 35, -.19, -.17);
                //hVertexY->Sumw2(kTRUE);
                
                hVertexZ = new TH1F("VertexZ", "VertexZ", 250, -5, 5);
                //hVertexZ->Sumw2(kTRUE);
                
                hVertexY_vs_VertexX = new TH2F("VertexY_vs_VertexX", "VertexY_vs_VertexX", 125, .145, .205, 125, -.215, -.155);
                //hhVertexY_vs_VertexX(kTRUE);
	        
                hVertexScore = new TH1F("VertexScore", "VertexScore", 125, 0, 35000);
                //hVertexZ->Sumw2(kTRUE);
                
            }

			~VertexSet() {}

            TString _postfix;
            TH1F *hVertexX;  
            TH1F *hVertexY;  
            TH1F *hVertexZ;  
            TH2F *hVertexY_vs_VertexX;
            TH1F *hVertexScore; 
		}; 

        struct PhotonSet
		{
			PhotonSet(TString const& postfix) 
			{
				_postfix = postfix;
			}

			void init()
            {                
                gDirectory->mkdir(_postfix); 
                gDirectory->cd(_postfix);    
                
                hPhotonpt = new TH1F("Photonpt", "Photonpt", 75, 0, 75);
                //hPhotonpt->Sumw2(kTRUE); 
                
                hPhotoneta = new TH1F("Photoneta", "Photoneta", 50, -2.5, 2.5);  
                //hPhotoneta->Sumw2(kTRUE);  
                
                hPhotonphi = new TH1F("Photonphi", "Photonphi", 36, -3.6, 3.6);  
                //hPhotonphi->Sumw2(kTRUE);   
	        }

			~PhotonSet() {}

            TString _postfix;
            TH1F *hPhotonpt;   
            TH1F *hPhotonphi;  
            TH1F *hPhotoneta;  
		};        
        
        struct HistSetCollection {
            analysis::processing::DimuonSet& setDM_PreSelection;
            analysis::processing::DimuonSet& corr_setDM_PreSelection;
            analysis::processing::JetSet& setJet_PreSelection;

            analysis::processing::DimuonSet& setDM_PreSelection_BB;
            analysis::processing::DimuonSet& corr_setDM_PreSelection_BB;
            analysis::processing::JetSet& setJet_PreSelection_BB;

            analysis::processing::DimuonSet& setDM_PreSelection_BO;
            analysis::processing::DimuonSet& corr_setDM_PreSelection_BO;
            analysis::processing::JetSet& setJet_PreSelection_BO;

            analysis::processing::DimuonSet& setDM_PreSelection_BE;
            analysis::processing::DimuonSet& corr_setDM_PreSelection_BE;
            analysis::processing::JetSet& setJet_PreSelection_BE;

            analysis::processing::DimuonSet& setDM_PreSelection_OO;
            analysis::processing::DimuonSet& corr_setDM_PreSelection_OO;
            analysis::processing::JetSet& setJet_PreSelection_OO;

            analysis::processing::DimuonSet& setDM_PreSelection_OE;
            analysis::processing::DimuonSet& corr_setDM_PreSelection_OE;
            analysis::processing::JetSet& setJet_PreSelection_OE;

            analysis::processing::DimuonSet& setDM_PreSelection_EE;
            analysis::processing::DimuonSet& corr_setDM_PreSelection_EE;
            analysis::processing::JetSet& setJet_PreSelection_EE;

            analysis::processing::DimuonSet& setDM_VBF;
            analysis::processing::DimuonSet& corr_setDM_VBF;
            analysis::processing::JetSet& setJet_VBF;

            analysis::processing::DimuonSet& setDM_VBF_BB;
            analysis::processing::DimuonSet& corr_setDM_VBF_BB;
            analysis::processing::JetSet& setJet_VBF_BB;

            analysis::processing::DimuonSet& setDM_VBF_BO;
            analysis::processing::DimuonSet& corr_setDM_VBF_BO;
            analysis::processing::JetSet& setJet_VBF_BO;

            analysis::processing::DimuonSet& setDM_VBF_BE;
            analysis::processing::DimuonSet& corr_setDM_VBF_BE;
            analysis::processing::JetSet& setJet_VBF_BE;

            analysis::processing::DimuonSet& setDM_VBF_OO;
            analysis::processing::DimuonSet& corr_setDM_VBF_OO;
            analysis::processing::JetSet& setJet_VBF_OO;

            analysis::processing::DimuonSet& setDM_VBF_OE;
            analysis::processing::DimuonSet& corr_setDM_VBF_OE;
            analysis::processing::JetSet& setJet_VBF_OE;

            analysis::processing::DimuonSet& setDM_VBF_EE;
            analysis::processing::DimuonSet& corr_setDM_VBF_EE;
            analysis::processing::JetSet& setJet_VBF_EE;
        };
    }
}
