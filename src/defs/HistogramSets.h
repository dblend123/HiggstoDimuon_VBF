#pragma once

#include <iostream>
#include <cstdlib>
#include "TString.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TH2F.h"

namespace analysis
{
	namespace processing
	{
        struct SystHist1D
        {
            TH1F* nom = nullptr;

            TH1F* scaleUp = nullptr;
            TH1F* scaleDown = nullptr;

            TH1F* pdfUp = nullptr;
            TH1F* pdfDown = nullptr;

            TH1F* isrUp = nullptr;
            TH1F* isrDown = nullptr;

            TH1F* fsrUp = nullptr;
            TH1F* fsrDown = nullptr;

            void book(const char* name, const char* title, int nbins, double xmin, double xmax)
            {
                nom       = new TH1F(name, title, nbins, xmin, xmax);

                scaleUp   = new TH1F(Form("%s_scaleUp",   name), Form("%s_scaleUp",   title), nbins, xmin, xmax);
                scaleDown = new TH1F(Form("%s_scaleDown", name), Form("%s_scaleDown", title), nbins, xmin, xmax);

                pdfUp     = new TH1F(Form("%s_pdfUp",     name), Form("%s_pdfUp",     title), nbins, xmin, xmax);
                pdfDown   = new TH1F(Form("%s_pdfDown",   name), Form("%s_pdfDown",   title), nbins, xmin, xmax);

                isrUp     = new TH1F(Form("%s_isrUp",     name), Form("%s_isrUp",     title), nbins, xmin, xmax);
                isrDown   = new TH1F(Form("%s_isrDown",   name), Form("%s_isrDown",   title), nbins, xmin, xmax);

                fsrUp     = new TH1F(Form("%s_fsrUp",     name), Form("%s_fsrUp",     title), nbins, xmin, xmax);
                fsrDown   = new TH1F(Form("%s_fsrDown",   name), Form("%s_fsrDown",   title), nbins, xmin, xmax);

                nom->Sumw2();
                scaleUp->Sumw2();   scaleDown->Sumw2();
                pdfUp->Sumw2();     pdfDown->Sumw2();
                isrUp->Sumw2();     isrDown->Sumw2();
                fsrUp->Sumw2();     fsrDown->Sumw2();
            }

            void fill(double x,
                      double w_nom,
                      double w_scaleUp, double w_scaleDown,
                      double w_pdfUp,   double w_pdfDown,
                      double w_isrUp,   double w_isrDown,
                      double w_fsrUp,   double w_fsrDown)
            {

                if (!nom) {
                    std::cerr << "[ERROR] Unbooked histogram in SystHist1D::fill()" << std::endl;
                    abort();
                }

                nom->Fill(x, w_nom);

                scaleUp->Fill(x, w_scaleUp);
                scaleDown->Fill(x, w_scaleDown);

                pdfUp->Fill(x, w_pdfUp);
                pdfDown->Fill(x, w_pdfDown);

                isrUp->Fill(x, w_isrUp);
                isrDown->Fill(x, w_isrDown);

                fsrUp->Fill(x, w_fsrUp);
                fsrDown->Fill(x, w_fsrDown);
            }
        };
		

        struct ExpSystHist1D
        {
            TH1F* nom = nullptr;

            TH1F* puUp = nullptr;
            TH1F* puDown = nullptr;

            TH1F* muonSFUp = nullptr;
            TH1F* muonSFDown = nullptr;

            TH1F* triggerUp = nullptr;
            TH1F* triggerDown = nullptr;

            TH1F* btagSFUp = nullptr;
            TH1F* btagSFDown = nullptr;

            TH1F* btagEffStatUp = nullptr;
            TH1F* btagEffStatDown = nullptr;

            TH1F* zptUp = nullptr;
            TH1F* zptDown = nullptr;

            TH1F* jecUp = nullptr;
            TH1F* jecDown = nullptr;
            
            TH1F* jerUp = nullptr;
            TH1F* jerDown = nullptr;
    
            void book(const char* name, const char* title, int nbins, double xmin, double xmax)
            {
                nom         = new TH1F(name, title, nbins, xmin, xmax);

                puUp       = new TH1F(Form("%s_puUp",       name), Form("%s_puUp",       title), nbins, xmin, xmax);
                puDown     = new TH1F(Form("%s_puDown",     name), Form("%s_puDown",     title), nbins, xmin, xmax);

                muonSFUp   = new TH1F(Form("%s_muonSFUp",   name), Form("%s_muonSFUp",   title), nbins, xmin, xmax);
                muonSFDown = new TH1F(Form("%s_muonSFDown", name), Form("%s_muonSFDown", title), nbins, xmin, xmax);

                triggerUp   = new TH1F(Form("%s_triggerUp",   name), Form("%s_triggerUp",   title), nbins, xmin, xmax);
                triggerDown = new TH1F(Form("%s_triggerDown", name), Form("%s_triggerDown", title), nbins, xmin, xmax);

                btagSFUp =
                    new TH1F(
                        Form("%s_btagSFUp", name),
                        Form("%s_btagSFUp", title),
                        nbins, xmin, xmax
                    );

                btagSFDown =
                    new TH1F(
                        Form("%s_btagSFDown", name),
                        Form("%s_btagSFDown", title),
                        nbins, xmin, xmax
                    );

                btagEffStatUp =
                    new TH1F(
                        Form("%s_btagEffStatUp", name),
                        Form("%s_btagEffStatUp", title),
                        nbins, xmin, xmax
                    );

                btagEffStatDown =
                    new TH1F(
                        Form("%s_btagEffStatDown", name),
                        Form("%s_btagEffStatDown", title),
                        nbins, xmin, xmax
                    );

                zptUp      = new TH1F(Form("%s_zptUp",      name), Form("%s_zptUp",      title), nbins, xmin, xmax);
                zptDown    = new TH1F(Form("%s_zptDown",    name), Form("%s_zptDown",    title), nbins, xmin, xmax);

                jecUp   = new TH1F(Form("%s_jecUp", name),   Form("%s_jecUp", title), nbins, xmin, xmax);
                jecDown = new TH1F(Form("%s_jecDown", name), Form("%s_jecDown", title), nbins, xmin, xmax);
                
                jerUp   = new TH1F(Form("%s_jerUp", name),   Form("%s_jerUp", title), nbins, xmin, xmax);
                jerDown = new TH1F(Form("%s_jerDown", name), Form("%s_jerDown", title), nbins, xmin, xmax);


                nom->Sumw2();
                puUp->Sumw2();          puDown->Sumw2();
                muonSFUp->Sumw2();      muonSFDown->Sumw2();
                triggerUp->Sumw2();     triggerDown->Sumw2();
                btagSFUp->Sumw2();      btagSFDown->Sumw2();
                btagEffStatUp->Sumw2(); btagEffStatDown->Sumw2();
                zptUp->Sumw2();         zptDown->Sumw2();
                jecUp->Sumw2();         jecDown->Sumw2();
                jerUp->Sumw2();         jerDown->Sumw2();
            }
            
            // ---------------------------------------------------------
            // Weight-only overload.
            //
            // This preserves the existing behavior for callers that do
            // not provide shifted JEC/JER values or category weights.
            // ---------------------------------------------------------
            void fill(
                double x,
                double w_nom,
                double w_puUp,          double w_puDown,
                double w_muonSFUp,      double w_muonSFDown,
                double w_triggerUp,     double w_triggerDown,
                double w_btagSFUp,      double w_btagSFDown,
                double w_btagEffStatUp, double w_btagEffStatDown,
                double w_zptUp,         double w_zptDown)
            {
                fill(
                    x,
                    w_nom,

                    w_puUp,
                    w_puDown,

                    w_muonSFUp,
                    w_muonSFDown,

                    w_triggerUp,
                    w_triggerDown,

                    w_btagSFUp,
                    w_btagSFDown,

                    w_btagEffStatUp,
                    w_btagEffStatDown,

                    w_zptUp,
                    w_zptDown,

                    x,
                    w_nom,

                    x,
                    w_nom,

                    x,
                    w_nom,

                    x,
                    w_nom
                );
            }


            // ---------------------------------------------------------
            // Full experimental-systematic overload.
            //
            // JEC and JER variations each receive:
            //   1. their own shifted observable value;
            //   2. their own category-gated event weight.
            // ---------------------------------------------------------
            void fill(
                double x,
                double w_nom,

                double w_puUp,
                double w_puDown,

                double w_muonSFUp,
                double w_muonSFDown,

                double w_triggerUp,
                double w_triggerDown,

                double w_btagSFUp,
                double w_btagSFDown,

                double w_btagEffStatUp,
                double w_btagEffStatDown,

                double w_zptUp,
                double w_zptDown,

                double x_jecUp,
                double w_jecUp,

                double x_jecDown,
                double w_jecDown,

                double x_jerUp,
                double w_jerUp,

                double x_jerDown,
                double w_jerDown)
            {
                if (!nom) {
                    std::cerr
                        << "[ERROR] Unbooked histogram in "
                        << "ExpSystHist1D::fill()\n";

                    std::abort();
                }

                // Do not call TH1::Fill with zero weight.
                // ROOT can still increment the Entries counter for a
                // zero-weight Fill, which obscures category migration.
                auto fillIfWeighted =
                    [](TH1F* hist,
                       double value,
                       double weight)
                    {
                        if (hist && weight != 0.0) {
                            hist->Fill(
                                value,
                                weight
                            );
                        }
                    };

                // Nominal
                fillIfWeighted(
                    nom,
                    x,
                    w_nom
                );

                // Weight-only experimental variations.
                // Their observable value remains nominal.
                fillIfWeighted(
                    puUp,
                    x,
                    w_puUp
                );

                fillIfWeighted(
                    puDown,
                    x,
                    w_puDown
                );

                fillIfWeighted(
                    muonSFUp,
                    x,
                    w_muonSFUp
                );

                fillIfWeighted(
                    muonSFDown,
                    x,
                    w_muonSFDown
                );

                fillIfWeighted(
                    triggerUp,
                    x,
                    w_triggerUp
                );

                fillIfWeighted(
                    triggerDown,
                    x,
                    w_triggerDown
                );

                fillIfWeighted(
                    btagSFUp,
                    x,
                    w_btagSFUp
                );

                fillIfWeighted(
                    btagSFDown,
                    x,
                    w_btagSFDown
                );

                fillIfWeighted(
                    btagEffStatUp,
                    x,
                    w_btagEffStatUp
                );

                fillIfWeighted(
                    btagEffStatDown,
                    x,
                    w_btagEffStatDown
                );

                fillIfWeighted(
                    zptUp,
                    x,
                    w_zptUp
                );

                fillIfWeighted(
                    zptDown,
                    x,
                    w_zptDown
                );

                // Shape + category variations.
                fillIfWeighted(
                    jecUp,
                    x_jecUp,
                    w_jecUp
                );

                fillIfWeighted(
                    jecDown,
                    x_jecDown,
                    w_jecDown
                );

                fillIfWeighted(
                    jerUp,
                    x_jerUp,
                    w_jerUp
                );

                fillIfWeighted(
                    jerDown,
                    x_jerDown,
                    w_jerDown
                );
            }
        };

        struct VtxSet
        {
			VtxSet(TDirectory* parentDir, TString const& postfix)
			    : _postfix(postfix)
            {
			    _dir = parentDir->mkdir(postfix);
                _dir->cd();
			
                h_nPV_nom = new TH1F("h_nPV_nom", ";\nPV;Events", 100, 0, 100);
                h_nPV_pu  = new TH1F("h_nPV_pu",  ";\nPV;Events", 100, 0, 100);
                
                parentDir->cd();
            }
                
                ~VtxSet() {}

                TDirectory* _dir;
                TString _postfix;

                TH1F *h_nPV_nom;
                TH1F *h_nPV_pu;

            };

        struct MuonSet
		{
			MuonSet(TDirectory* parentDir, TString const& postfix)
			    : _postfix(postfix)
            {
			    _dir = parentDir->mkdir(postfix);
                _dir->cd();
			
                // -----  DiMuon Plots  -----
				hDiMuonMass         = new TH1F("DiMuonMass", "DiMuonMass", 130, 70, 200);
				hDiMuonMass_Z       = new TH1F("DiMuonMass_Z", "DiMuonMass_Z", 120, 70, 110);
				hDiMuonMass_Higgs   = new TH1F("DiMuonMass_Higgs", "DiMuonMass_Higgs", 180, 110, 150);
                hDiMuonpt           = new TH1F("DiMuonpt", "DiMuonpt", 60, 0, 200);
                hDiMuoneta          = new TH1F("DiMuoneta", "DiMuoneta", 50, -2.5, 2.5);
                hDiMuonrapidity     = new TH1F("DiMuonrapidity", "DiMuonrapidity", 50, -2.5, 2.5);
                hCollinsSoper       = new TH1F("CollinsSoper","CollinsSoper", 30, -3.2, 3.2);
                hcosCollinsSoper    = new TH1F("cosCollinsSoper","cosCollinsSoper", 10, -1.0, 1.0);
                hPtBalance          = new TH1F("PtBalance", "PtBalance", 30, 0, 1.0);
                hZeppenfeld         = new TH1F("Zeppenfeld", "Zeppenfeld", 20, -5, 5);

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
                TH1F *hDiMuonMass_Z;
                TH1F *hDiMuonMass_Higgs;
                TH1F *hDiMuonpt;
                TH1F *hDiMuoneta;
                TH1F *hDiMuonrapidity;
                TH1F *hCollinsSoper;
                TH1F *hcosCollinsSoper;
                TH1F *hPtBalance;
                TH1F *hZeppenfeld;

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
                hdphi_mumu_jj = new TH1F("hdphi_mumu_jj", "hdphi_mumu_jj", 40, 0, 3.2);
                hdeta_mumu_jj = new TH1F("hdeta_mumu_jj", "hdeta_mumu_jj", 40, 0, 8);

                // -----  Jet Plots  -----
                hJet1pt = new TH1F("Jet1Pt", "Jet1Pt", 50, 0, 200);
                hJet2pt = new TH1F("Jet2Pt", "Jet2Pt", 50, 0, 200);
                hJet1eta = new TH1F("Jet1eta", "Jet1eta", 60, -5, 5);
                hJet2eta = new TH1F("Jet2eta", "Jet2eta", 60, -5, 5);
                hnJets = new TH1F("nJets", "nJets", 8, 0, 8);
                
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
                TH1F *hdphi_mumu_jj;
                TH1F *hdeta_mumu_jj;
                
                // -----  Jet  -----
                TH1F *hJet1pt;   
                TH1F *hJet2pt;   
                TH1F *hJet1eta;
                TH1F *hJet2eta;
                TH1F *hnJets;
		};      
        
        struct BDTTheorySystSet
        {
            BDTTheorySystSet(TDirectory* parentDir, const TString& postfix)
                : _postfix(postfix)
            {
                _dir = parentDir->mkdir(_postfix);
                _dir->cd();

                h_z.book("Zeppenfeld", "Zeppenfeld", 20, -5, 5);

                h_ptBalance.book("PtBalance", "PtBalance", 30, 0, 1.0);

                h_mu1_pt_over_m.book("Muon1PTdDM", "Muon1PTdDM", 30, 0, 1.5);
                h_mu2_pt_over_m.book("Muon2PTdDM", "Muon2PTdDM", 30, 0, 1.5);

                h_mu1_eta.book("Muon1eta", "Muon1eta", 50, -2.5, 2.5);
                h_mu2_eta.book("Muon2eta", "Muon2eta", 50, -2.5, 2.5);

                h_jet1_eta.book("Jet1eta", "Jet1eta", 60, -5, 5);

                h_m_mumu.book("DiMuonMass_Higgs", "DiMuonMass_Higgs", 180, 110, 150);

                h_deta_jj.book("DiJetdeta", "DiJetdeta", 30, 0, 9);
                h_dphi_jj.book("DiJetdphi", "DiJetdphi", 30, 0, 3.2);

                h_njets.book("nJets", "nJets", 8, 0, 8);

                h_dphi_mumu_jj.book("hdphi_mumu_jj", "hdphi_mumu_jj", 40, 0, 3.2);
                h_deta_mumu_jj.book("hdeta_mumu_jj", "hdeta_mumu_jj", 40, 0, 8);

                h_phi_CS.book("CollinsSoper", "CollinsSoper", 30, -3.2, 3.2);
                h_costheta_CS.book("cosCollinsSoper", "cosCollinsSoper", 10, -1, 1);

                h_y_mumu.book("DiMuonrapidity", "DiMuonrapidity", 50, -2.5, 2.5);
                h_y_j1.book("Jet1Rapidity", "Jet1Rapidity", 50, -5, 5);
                h_y_j2.book("Jet2Rapidity", "Jet2Rapidity", 50, -5, 5);

                h_pt_mumu.book("DiMuonpt", "DiMuonpt", 60, 0, 200);

                h_pt_j1.book("Jet1Pt", "Jet1Pt", 50, 0, 200);
                h_pt_j2.book("Jet2Pt", "Jet2Pt", 50, 0, 200);
                
                parentDir->cd();
            }

            ~BDTTheorySystSet() {}

            TDirectory* _dir;
            TString _postfix;

            SystHist1D h_z;
            SystHist1D h_ptBalance;

            SystHist1D h_mu1_pt_over_m;
            SystHist1D h_mu2_pt_over_m;

            SystHist1D h_mu1_eta;
            SystHist1D h_mu2_eta;

            SystHist1D h_jet1_eta;

            SystHist1D h_m_mumu;
            SystHist1D h_deta_jj;
            SystHist1D h_dphi_jj;
            SystHist1D h_njets;

            SystHist1D h_dphi_mumu_jj;
            SystHist1D h_deta_mumu_jj;

            SystHist1D h_phi_CS;
            SystHist1D h_costheta_CS;

            SystHist1D h_y_mumu;
            SystHist1D h_y_j1;
            SystHist1D h_y_j2;

            SystHist1D h_pt_mumu;
            SystHist1D h_pt_j1;
            SystHist1D h_pt_j2;
        };


        struct BDTExperimentalSystSet
        {
            BDTExperimentalSystSet(TDirectory* parentDir, const TString& postfix)
                : _postfix(postfix)
            {
                _dir = parentDir->mkdir(_postfix);
                _dir->cd();

                h_z.book("Zeppenfeld", "Zeppenfeld", 20, -5, 5);

                h_ptBalance.book("PtBalance", "PtBalance", 30, 0, 1.0);

                h_mu1_pt_over_m.book("Muon1PTdDM", "Muon1PTdDM", 30, 0, 1.5);
                h_mu2_pt_over_m.book("Muon2PTdDM", "Muon2PTdDM", 30, 0, 1.5);

                h_mu1_eta.book("Muon1eta", "Muon1eta", 50, -2.5, 2.5);
                h_mu2_eta.book("Muon2eta", "Muon2eta", 50, -2.5, 2.5);

                h_jet1_eta.book("Jet1eta", "Jet1eta", 60, -5, 5);

                h_m_mumu.book("DiMuonMass_Higgs", "DiMuonMass_Higgs", 180, 110, 150);

                h_deta_jj.book("DiJetdeta", "DiJetdeta", 30, 0, 9);
                h_dphi_jj.book("DiJetdphi", "DiJetdphi", 30, 0, 3.2);

                h_njets.book("nJets", "nJets", 8, 0, 8);

                h_dphi_mumu_jj.book("hdphi_mumu_jj", "hdphi_mumu_jj", 40, 0, 3.2);
                h_deta_mumu_jj.book("hdeta_mumu_jj", "hdeta_mumu_jj", 40, 0, 8);

                h_phi_CS.book("CollinsSoper", "CollinsSoper", 30, -3.2, 3.2);
                h_costheta_CS.book("cosCollinsSoper", "cosCollinsSoper", 10, -1, 1);

                h_y_mumu.book("DiMuonrapidity", "DiMuonrapidity", 50, -2.5, 2.5);
                h_y_j1.book("Jet1Rapidity", "Jet1Rapidity", 50, -5, 5);
                h_y_j2.book("Jet2Rapidity", "Jet2Rapidity", 50, -5, 5);

                h_pt_mumu.book("DiMuonpt", "DiMuonpt", 60, 0, 200);

                h_pt_j1.book("Jet1Pt", "Jet1Pt", 50, 0, 200);
                h_pt_j2.book("Jet2Pt", "Jet2Pt", 50, 0, 200);
                
                parentDir->cd();
            }

            ~BDTExperimentalSystSet() {}

            TDirectory* _dir;
            TString _postfix;

            ExpSystHist1D h_z;
            ExpSystHist1D h_ptBalance;

            ExpSystHist1D h_mu1_pt_over_m;
            ExpSystHist1D h_mu2_pt_over_m;

            ExpSystHist1D h_mu1_eta;
            ExpSystHist1D h_mu2_eta;

            ExpSystHist1D h_jet1_eta;

            ExpSystHist1D h_m_mumu;
            ExpSystHist1D h_deta_jj;
            ExpSystHist1D h_dphi_jj;
            ExpSystHist1D h_njets;

            ExpSystHist1D h_dphi_mumu_jj;
            ExpSystHist1D h_deta_mumu_jj;

            ExpSystHist1D h_phi_CS;
            ExpSystHist1D h_costheta_CS;

            ExpSystHist1D h_y_mumu;
            ExpSystHist1D h_y_j1;
            ExpSystHist1D h_y_j2;

            ExpSystHist1D h_pt_mumu;
            ExpSystHist1D h_pt_j1;
            ExpSystHist1D h_pt_j2;
        };

        struct BTagEffSet
        {
            BTagEffSet(TDirectory* parentDir, TString const& postfix)
                : _postfix(postfix)
            {
                _dir = parentDir->mkdir(postfix);
                _dir->cd();

                const int nPtBins = 9;
                double ptBins[nPtBins + 1] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};

                const int nEtaBins = 4;
                double etaBins[nEtaBins + 1] = {0.0, 0.8, 1.6, 2.0, 2.5};

                // --- Loose WP ---
                h_den_b_L = new TH2F("h_den_b_L", "b denom L;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                h_num_b_L = new TH2F("h_num_b_L", "b numer L;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);

                h_den_c_L = new TH2F("h_den_c_L", "c denom L;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                h_num_c_L = new TH2F("h_num_c_L", "c numer L;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);

                h_den_l_L = new TH2F("h_den_l_L", "light denom L;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                h_num_l_L = new TH2F("h_num_l_L", "light numer L;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);

                // --- Medium WP ---
                h_den_b_M = new TH2F("h_den_b_M", "b denom M;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                h_num_b_M = new TH2F("h_num_b_M", "b numer M;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);

                h_den_c_M = new TH2F("h_den_c_M", "c denom M;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                h_num_c_M = new TH2F("h_num_c_M", "c numer M;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);

                h_den_l_M = new TH2F("h_den_l_M", "light denom M;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                h_num_l_M = new TH2F("h_num_l_M", "light numer M;jet p_{T} [GeV];|#eta|",
                                     nPtBins, ptBins, nEtaBins, etaBins);
                
                h_den_b_L->Sumw2();
                h_num_b_L->Sumw2();
                h_den_c_L->Sumw2();
                h_num_c_L->Sumw2();
                h_den_l_L->Sumw2();
                h_num_l_L->Sumw2();

                h_den_b_M->Sumw2();
                h_num_b_M->Sumw2();
                h_den_c_M->Sumw2();
                h_num_c_M->Sumw2();
                h_den_l_M->Sumw2();
                h_num_l_M->Sumw2();
                
                parentDir->cd();
            }

            ~BTagEffSet() {}

            TDirectory* _dir;
            TString _postfix;

            TH2F *h_den_b_L, *h_num_b_L;
            TH2F *h_den_c_L, *h_num_c_L;
            TH2F *h_den_l_L, *h_num_l_L;

            TH2F *h_den_b_M, *h_num_b_M;
            TH2F *h_den_c_M, *h_num_c_M;
            TH2F *h_den_l_M, *h_num_l_M;
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
            
            // Fsr Block
            //analysis::processing::MuonSet& setDM_PreFsr;
            //analysis::processing::MuonSet& setDM_PostFsr;

            // Vertex Block
            analysis::processing::VtxSet& PileupReweighting_PreSelection;
            analysis::processing::VtxSet& PileupReweighting_VBF;

            // Z-Reweighting Block
            analysis::processing::MuonSet& corr_setDM_PreSelection_Z0j;
            analysis::processing::MuonSet& corr_setDM_PreSelection_Z1j;
            analysis::processing::MuonSet& corr_setDM_PreSelection_Z2j;
            analysis::processing::MuonSet& corr_setDM_PreSelection_wZpt_Z0j;
            analysis::processing::MuonSet& corr_setDM_PreSelection_wZpt_Z1j;
            analysis::processing::MuonSet& corr_setDM_PreSelection_wZpt_Z2j;
           
            // B-Tag Reweighting Block
            analysis::processing::BTagEffSet& BTagEff_PreSelection;
            analysis::processing::BTagEffSet& BTagEff_VBF;

            // BDT Theory Syst Block
            analysis::processing::BDTTheorySystSet& BDTTheorySyst_VBF;

            // BDT Experimental Syst Block
            analysis::processing::BDTExperimentalSystSet& BDTExperimentalSyst_VBF;

            // Pre-Selection Block
            analysis::processing::MuonSet& setDM_PreSelection;
            analysis::processing::MuonSet& corr_setDM_PreSelection;
            analysis::processing::MuonSet& corr_setDM_PreSelection_Z;
            analysis::processing::MuonSet& corr_setDM_PreSelection_SR;
            analysis::processing::JetSet&  setJet_PreSelection;
            analysis::processing::JetSet&  corr_setJet_PreSelection;

            analysis::processing::MuonSet& setDM_PreSelection_BB;
            analysis::processing::MuonSet& corr_setDM_PreSelection_BB;
            analysis::processing::JetSet& setJet_PreSelection_BB;
            analysis::processing::JetSet& corr_setJet_PreSelection_BB;

            analysis::processing::MuonSet& setDM_PreSelection_BO;
            analysis::processing::MuonSet& corr_setDM_PreSelection_BO;
            analysis::processing::JetSet& setJet_PreSelection_BO;
            analysis::processing::JetSet& corr_setJet_PreSelection_BO;

            analysis::processing::MuonSet& setDM_PreSelection_BE;
            analysis::processing::MuonSet& corr_setDM_PreSelection_BE;
            analysis::processing::JetSet& setJet_PreSelection_BE;
            analysis::processing::JetSet& corr_setJet_PreSelection_BE;

            analysis::processing::MuonSet& setDM_PreSelection_OO;
            analysis::processing::MuonSet& corr_setDM_PreSelection_OO;
            analysis::processing::JetSet& setJet_PreSelection_OO;
            analysis::processing::JetSet& corr_setJet_PreSelection_OO;

            analysis::processing::MuonSet& setDM_PreSelection_OE;
            analysis::processing::MuonSet& corr_setDM_PreSelection_OE;
            analysis::processing::JetSet& setJet_PreSelection_OE;
            analysis::processing::JetSet& corr_setJet_PreSelection_OE;

            analysis::processing::MuonSet& setDM_PreSelection_EE;
            analysis::processing::MuonSet& corr_setDM_PreSelection_EE;
            analysis::processing::JetSet& setJet_PreSelection_EE;
            analysis::processing::JetSet& corr_setJet_PreSelection_EE;

            // Pre-DiJet Block
            analysis::processing::JetSet& setJet_PreDiJet_VBF;
            analysis::processing::JetSet& corr_setJet_PreDiJet_VBF;

            // VBF Block
            analysis::processing::MuonSet& setDM_VBF;
            analysis::processing::MuonSet& corr_setDM_VBF;
            analysis::processing::MuonSet& corr_setDM_VBF_Z;
            analysis::processing::MuonSet& corr_setDM_VBF_SR;
            analysis::processing::JetSet&  setJet_VBF;
            analysis::processing::JetSet&  corr_setJet_VBF;
            analysis::processing::JetSet&  corr_setJet_VBF_SR;

            analysis::processing::MuonSet& setDM_VBF_BB;
            analysis::processing::MuonSet& corr_setDM_VBF_BB;
            analysis::processing::JetSet& setJet_VBF_BB;
            analysis::processing::JetSet& corr_setJet_VBF_BB;

            analysis::processing::MuonSet& setDM_VBF_BO;
            analysis::processing::MuonSet& corr_setDM_VBF_BO;
            analysis::processing::JetSet& setJet_VBF_BO;
            analysis::processing::JetSet& corr_setJet_VBF_BO;

            analysis::processing::MuonSet& setDM_VBF_BE;
            analysis::processing::MuonSet& corr_setDM_VBF_BE;
            analysis::processing::JetSet& setJet_VBF_BE;
            analysis::processing::JetSet& corr_setJet_VBF_BE;

            analysis::processing::MuonSet& setDM_VBF_OO;
            analysis::processing::MuonSet& corr_setDM_VBF_OO;
            analysis::processing::JetSet& setJet_VBF_OO;
            analysis::processing::JetSet& corr_setJet_VBF_OO;

            analysis::processing::MuonSet& setDM_VBF_OE;
            analysis::processing::MuonSet& corr_setDM_VBF_OE;
            analysis::processing::JetSet& setJet_VBF_OE;
            analysis::processing::JetSet& corr_setJet_VBF_OE;

            analysis::processing::MuonSet& setDM_VBF_EE;
            analysis::processing::MuonSet& corr_setDM_VBF_EE;
            analysis::processing::JetSet& setJet_VBF_EE;
            analysis::processing::JetSet& corr_setJet_VBF_EE;
        };
    }
}
