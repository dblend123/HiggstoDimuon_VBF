#pragma once

#include "TString.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TH2F.h"

namespace analysis
{
	namespace processing
	{
		struct MuonSet
		{
			MuonSet(TDirectory* parentDir, TString const& postfix)
			    : _postfix(postfix)
            {
			    _dir = parentDir->mkdir(postfix);
                _dir->cd();
			
                // -----  DiMuon Plots  -----
				hDiMuonMass = new TH1F("DiMuonMass", "DiMuonMass", 150, 50, 200);
                hDiMuonpt = new TH1F("DiMuonpt", "DiMuonpt", 50, 0, 200);
                hDiMuoneta = new TH1F("DiMuoneta", "DiMuoneta", 50, -2.5, 2.5);
                hDiMuonrapidity = new TH1F("DiMuonrapidity", "DiMuonrapidity", 50, -2.5, 2.5);
                //hDiMuonCS = new TH1F("DiMuonCS","DiMounCS", 30, -3.2, 3.2);
                //hDiMuoncosCS = new TH1F("DiMuoncosCS","DiMouncosCS", 30, -3.2, 3.2);

                // -----  Muon Plots  -----
				hMuon1pt = new TH1F("Muon1pt", "Muon1pt", 50, 0, 200);
                hMuon2pt = new TH1F("Muon2pt", "Muon2pt", 50, 0, 200);
				hMuon1eta = new TH1F("Muon1eta", "Muon1eta", 50, -2.5, 2.5);
				hMuon2eta = new TH1F("Muon2eta", "Muon2eta", 50, -2.5, 2.5);
				hMuon1PTdDM = new TH1F("Muon1PTdDM", "Muon1PTdDM", 30, 0, 1.5);
                hMuon2PTdDM = new TH1F("Muon2PTdDM", "Muon2PTdDM", 30, 0, 1.5);
                
                parentDir->cd();
			}
                ~MuonSet() {}
                
                TDirectory* _dir;
                TString _postfix;
                
                // -----  DiMuon  -----
                TH1F *hDiMuonMass;
                TH1F *hDiMuonpt;
                TH1F *hDiMuoneta;
                TH1F *hDiMuonrapidity;
                //TH1F *hDiMuonCS;
                //TH1F *hDiMuoncosCS;
                
                // -----  Muon  -----
                TH1F *hMuon1pt;
                TH1F *hMuon2pt;
                TH1F *hMuon1eta;
                TH1F *hMuon2eta;
                TH1F *hMuon1PTdDM;
                TH1F *hMuon2PTdDM;
        };
       
        struct JetSet
		{
			JetSet(TDirectory* parentDir, const TString& postfix) 
			    : _postfix(postfix)
            {
			    _dir = parentDir->mkdir(_postfix);
                _dir->cd();
                
                // -----  DiJet Plots  -----
                hDiJetMass = new TH1F("DiJetMass", "DiJetMass", 50, 0, 1000);
                hDiJetpt   = new TH1F("DiJetPt", "DiJetPt", 50, 0, 400);
                hDiJetdeta = new TH1F("DiJetdeta", "DiJetdeta", 30, 0, 9);
                hDiJetdphi = new TH1F("DiJetdphi", "DiJetdphi", 30, 0, 3.2);
                
                // -----  Jet Plots  -----
                hJet1pt = new TH1F("Jet1Pt", "Jet1Pt", 50, 0, 200);
                hJet2pt = new TH1F("Jet2Pt", "Jet2Pt", 50, 0, 200);
                hJet1eta = new TH1F("Jet1eta", "Jet1eta", 60, -5, 5);
                hJet2eta = new TH1F("Jet2eta", "Jet2eta", 60, -5, 5);
                //hnJets = new TH1F("nJets", "nJets", 8, 0, 8);
                
                parentDir->cd();
            }

                ~JetSet() {}
                    
                TDirectory* _dir;
                TString _postfix;

                // -----  DiJet  -----
                TH1F *hDiJetMass;
                TH1F *hDiJetpt;   
                TH1F *hDiJetdeta;
                TH1F *hDiJetdphi;
                
                // -----  Jet  -----
                TH1F *hJet1pt;   
                TH1F *hJet2pt;   
                TH1F *hJet1eta;
                TH1F *hJet2eta;
                //TH1F *hnJets;
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
            
            //Fsr blocks
            //analysis::processing::MuonSet& setDM_PreFsr;
            //analysis::processing::MuonSet& setDM_PostFsr;

            //Pre-Selection blocks
            analysis::processing::MuonSet& setDM_PreSelection;
            analysis::processing::MuonSet& corr_setDM_PreSelection;
            analysis::processing::JetSet&  setJet_PreSelection;
            analysis::processing::JetSet&  corr_setJet_PreSelection;

            /*analysis::processing::MuonSet& setDM_PreSelection_BB;
            analysis::processing::MuonSet& corr_setDM_PreSelection_BB;
            analysis::processing::JetSet& setJet_PreSelection_BB;

            analysis::processing::MuonSet& setDM_PreSelection_BO;
            analysis::processing::MuonSet& corr_setDM_PreSelection_BO;
            analysis::processing::JetSet& setJet_PreSelection_BO;

            analysis::processing::MuonSet& setDM_PreSelection_BE;
            analysis::processing::MuonSet& corr_setDM_PreSelection_BE;
            analysis::processing::JetSet& setJet_PreSelection_BE;

            analysis::processing::MuonSet& setDM_PreSelection_OO;
            analysis::processing::MuonSet& corr_setDM_PreSelection_OO;
            analysis::processing::JetSet& setJet_PreSelection_OO;

            analysis::processing::MuonSet& setDM_PreSelection_OE;
            analysis::processing::MuonSet& corr_setDM_PreSelection_OE;
            analysis::processing::JetSet& setJet_PreSelection_OE;

            analysis::processing::MuonSet& setDM_PreSelection_EE;
            analysis::processing::MuonSet& corr_setDM_PreSelection_EE;
            analysis::processing::JetSet& setJet_PreSelection_EE;*/

            // pre-dijet
            analysis::processing::JetSet& setJet_PreDiJet_VBF;
            analysis::processing::JetSet& corr_setJet_PreDiJet_VBF;

            // VBF blocks
            analysis::processing::MuonSet& setDM_VBF;
            analysis::processing::MuonSet& corr_setDM_VBF;
            analysis::processing::JetSet&  setJet_VBF;
            analysis::processing::JetSet&  corr_setJet_VBF;

            /*analysis::processing::MuonSet& setDM_VBF_BB;
            analysis::processing::MuonSet& corr_setDM_VBF_BB;
            analysis::processing::JetSet& setJet_VBF_BB;

            analysis::processing::MuonSet& setDM_VBF_BO;
            analysis::processing::MuonSet& corr_setDM_VBF_BO;
            analysis::processing::JetSet& setJet_VBF_BO;

            analysis::processing::MuonSet& setDM_VBF_BE;
            analysis::processing::MuonSet& corr_setDM_VBF_BE;
            analysis::processing::JetSet& setJet_VBF_BE;

            analysis::processing::MuonSet& setDM_VBF_OO;
            analysis::processing::MuonSet& corr_setDM_VBF_OO;
            analysis::processing::JetSet& setJet_VBF_OO;

            analysis::processing::MuonSet& setDM_VBF_OE;
            analysis::processing::MuonSet& corr_setDM_VBF_OE;
            analysis::processing::JetSet& setJet_VBF_OE;

            analysis::processing::MuonSet& setDM_VBF_EE;
            analysis::processing::MuonSet& corr_setDM_VBF_EE;
            analysis::processing::JetSet& setJet_VBF_EE;*/
        };
    }
}
