#include <fstream>  
#include <iostream>
#include <vector>
#include <string>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace ROOT;
using namespace ROOT::VecOps;

void Tag_and_Probe() {
    std::ifstream infile("0DYBackground_Run22EE.list");
    std::vector<std::string> fileList;
    std::string line;

    while (std::getline(infile, line)) {
        if (!line.empty()) fileList.push_back(line);
    }

    ROOT::RDataFrame df("Events", fileList);

// Step 1: Only keep events with ≥2 muons
auto df2mu = df.Filter("nMuon >= 2");

std::cout << "nMuons: " << df2mu.Count().GetValue() << std::endl;

// Step 2: Define tag and probe masks  
// The Tag Muon: muons that pass tight selection
auto df_with_tags = df2mu.Define("tag_mask", R"cpp(
    std::vector<bool> mask(Muon_pt.size(), false);
    if (HLT_IsoMu24) {
        for (size_t i = 0; i < Muon_pt.size(); ++i) {
            if (!Muon_tightId[i]) continue;
            if (Muon_pt[i] <= 27) continue;
            if (Muon_pfRelIso04_all[i] >= 0.15) continue;

            bool matched = false;
            for (size_t j = 0; j < TrigObj_pt.size(); ++j) {
                if (TrigObj_id[j] != 13) continue;
                if (TrigObj_pt[j] < 24) continue;

                float deta = Muon_eta[i] - TrigObj_eta[j];
                float dphi = std::abs(Muon_phi[i] - TrigObj_phi[j]);
                if (dphi > M_PI) dphi = 2 * M_PI - dphi;
                float dr2 = deta*deta + dphi*dphi;
                if (dr2 < 0.01) { // dR < 0.1
                    matched = true;
                    break;
                }
            }
            if (matched) mask[i] = true;
        }
    }
    return mask;
    )cpp");

auto count_tag_quality = df_with_tags.Filter("Sum(tag_mask) > 0").Count();
std::cout << "After tag_mask: " << count_tag_quality.GetValue() << std::endl;

// The Probe Muon: all other muons in event
auto df_with_probes = df_with_tags.Define("probe_mask", R"cpp(
    std::vector<bool> mask;
    for (bool is_tag : tag_mask) {
        mask.push_back(!is_tag);
    }
    return mask;
    )cpp");

auto count_probe_quality = df_with_probes.Filter("Sum(probe_mask) > 0").Count();
std::cout << "After probe_mask: " << count_probe_quality.GetValue() << std::endl;

auto df_with_good_probes = df_with_probes.Define("probe_quality_mask", R"cpp(
    std::vector<bool> mask(Muon_pt.size(), false);
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!probe_mask[i]) continue;
        
        // Require at least two muon stations
        if (Muon_nStations[i] <= 1) continue;

        // Optional: use Muon_highPtId (0=fail, 1=loose, 2=tight)
        if (Muon_highPtId[i] < 1) continue;

        // Optional: require minimum tracker layers
        // if (Muon_nTrackerLayers[i] <= 5) continue;

        mask[i] = true;
    }
    return mask;
)cpp");

auto count_good_probe_quality = df_with_good_probes.Filter("Sum(probe_quality_mask) > 0").Count();
std::cout << "After probe_quality_mask: " << count_good_probe_quality.GetValue() << std::endl;

auto df_with_passing_probes = df_with_good_probes.Define("passing_probe_mask", R"cpp(
    std::vector<bool> mask(Muon_pt.size(), false);
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!probe_quality_mask[i]) continue;

        // Search for a "track-like" matching muon
        float best_dr2 = 999;
        int best_match = -1;

        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (i == j) continue; // skip self

            // Charge match
            //if (Muon_charge[i] == Muon_charge[j]) continue;

            // "Track-like" kinematics
            float p = Muon_pt[j] * std::cosh(Muon_eta[j]);
            //if (!(Muon_pt[j] > 2 && (p > 2 || std::abs(Muon_eta[j]) < 1.0))) continue;
            if (!(Muon_pt[j] > 1 && (p > 1 || std::abs(Muon_eta[j]) < 2.4))) continue;

            // Matching criteria
            float deta = Muon_eta[i] - Muon_eta[j];
            float dphi = std::abs(Muon_phi[i] - Muon_phi[j]);
            if (dphi > M_PI) dphi = 2 * M_PI - dphi;
            float dr2 = deta*deta + dphi*dphi;
            float dpt_rel = std::abs(Muon_pt[i] - Muon_pt[j]) / Muon_pt[i];
            
            //if (dr2 < 0.09 && std::abs(deta) < 0.3 && dpt_rel < 2.0) {
            if (dr2 < 0.25 && std::abs(deta) < 0.5 && dpt_rel < 5.0) {
                if (dr2 < best_dr2) {
                    best_dr2 = dr2;
                    best_match = j;
                }
            }
        }

        //if (best_match >= 0) { comment
            // Require Δvz < 1.0 with *some* tag muon
            //bool matched_to_tag = false; comment
            //for (size_t k = 0; k < Muon_pt.size(); ++k) { comment
                //if (!tag_mask[k]) continue; comment
                //float dvz = std::abs(Muon_vz[best_match] - Muon_vz[k]);
                //if (dvz < 1.0) {
                    //matched_to_tag = true;
                    //break;
                //}
            //} comment

            //if (matched_to_tag) comment
                //mask[i] = true; comment
        //} comment
    }

    return mask;
)cpp");

auto count_passing_probe_quality = df_with_passing_probes.Filter("Sum(passing_probe_mask) > 0").Count();
std::cout << "After passing_probe_mask: " << count_passing_probe_quality.GetValue() << std::endl;

// Step 3: Build invariant mass of PASSING tag-probe pairs
// Passing pt
df_with_passing_probes = df_with_passing_probes.Define("probe_pt_pass", R"cpp(
    std::vector<double> pts;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                pts.push_back(Muon_pt[j]);
            }
        }
    }
    return pts;
)cpp")

// Passing eta
.Define("probe_eta_pass", R"cpp(
    std::vector<double> etas;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                etas.push_back(Muon_eta[j]);
            }
        }
    }
    return etas;
)cpp")

// Passing phi
.Define("probe_phi_pass", R"cpp(
    std::vector<double> phis;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                phis.push_back(Muon_phi[j]);
            }
        }
    }
    return phis;
)cpp")

// Passing dimuon mass
.Define("probe_mass_pass", R"cpp(
    std::vector<double> masses;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                masses.push_back(m);
            }
        }
    }
    return masses;
)cpp");

// Step 4: Build invariant mass of FAILING tag-probe pairs
// Failing pt
df_with_passing_probes = df_with_passing_probes.Define("probe_pt_fail", R"cpp(
    std::vector<double> pts;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                pts.push_back(Muon_pt[j]);
            }
        }
    }
    return pts;
)cpp")

// Failing eta
.Define("probe_eta_fail", R"cpp(
    std::vector<double> etas;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                etas.push_back(Muon_eta[j]);
            }
        }
    }
    return etas;
)cpp")

// Failing phi
.Define("probe_phi_fail", R"cpp(
    std::vector<double> phis;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                phis.push_back(Muon_phi[j]);
            }
        }
    }
    return phis;
)cpp")

// Failing dimuon mass
.Define("probe_mass_fail", R"cpp(
    std::vector<double> masses;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115) {
                masses.push_back(m);
            }
        }
    }
    return masses;
)cpp");

// Step 5: Plot histograms
    // Masses
    auto h_mass_pass = df_with_passing_probes.Histo1D({"probe_mass_pass", "Passing Probes; m_{#mu#mu} [GeV]; Events", 45, 70, 115}, "probe_mass_pass");
    auto h_mass_fail = df_with_passing_probes.Histo1D({"probe_mass_fail", "Failing Probes; m_{#mu#mu} [GeV]; Events", 45, 70, 115}, "probe_mass_fail");
    // Clone passing mass histogram
    auto h_mass_total = (TH1D*) h_mass_pass.GetPtr()->Clone("h_mass_total");
    h_mass_total->Add(h_mass_fail.GetPtr());
    // Create efficiency graph
    auto eff_mass = new TGraphAsymmErrors();
    eff_mass->BayesDivide(h_mass_pass.GetPtr(), h_mass_total);

    // Pts
    auto h_pt_pass  = df_with_passing_probes.Histo1D({"probe_pt_pass",  "Passing Probes; Probe p_{T} [GeV]; Events", 20, 0, 100}, "probe_pt_pass");
    df_with_passing_probes.Display("probe_pt_pass")->Print();
    auto h_pt_fail  = df_with_passing_probes.Histo1D({"probe_pt_fail",  "Failing Probes; Probe p_{T} [GeV]; Events", 20, 0, 100}, "probe_pt_fail");
    // Clone passing pt histogram
    auto h_pt_total = (TH1D*) h_pt_pass.GetPtr()->Clone("h_pt_total");
    h_pt_total->Add(h_pt_fail.GetPtr());
    // Create efficiency graph
    auto eff_pt = new TGraphAsymmErrors();
    eff_pt->BayesDivide(h_pt_pass.GetPtr(), h_pt_total);
    
    // Etas
    auto h_eta_pass = df_with_passing_probes.Histo1D({"probe_eta_pass", "Passing Probes; Probe #eta; Events", 24, -2.4, 2.4}, "probe_eta_pass");
    auto h_eta_fail = df_with_passing_probes.Histo1D({"probe_eta_fail", "Failing Probes; Probe #eta; Events", 24, -2.4, 2.4}, "probe_eta_fail");
    // Clone passing eta histogram
    auto h_eta_total = (TH1D*) h_eta_pass.GetPtr()->Clone("h_eta_total");
    h_eta_total->Add(h_eta_fail.GetPtr());
    // Create efficiency graph
    auto eff_eta = new TGraphAsymmErrors();
    eff_eta->BayesDivide(h_eta_pass.GetPtr(), h_eta_total);
    
    // Phis
    auto h_phi_pass = df_with_passing_probes.Histo1D({"probe_phi_pass", "Passing Probes; Probe #phi; Events", 20, -3.2, 3.2}, "probe_phi_pass");
    auto h_phi_fail = df_with_passing_probes.Histo1D({"probe_phi_fail", "Failing Probes; Probe #phi; Events", 20, -3.2, 3.2}, "probe_phi_fail");
    // Clone passing phi histogram
    auto h_phi_total = (TH1D*) h_phi_pass.GetPtr()->Clone("h_phi_total");
    h_phi_total->Add(h_phi_fail.GetPtr());
    // Create efficiency graph
    auto eff_phi = new TGraphAsymmErrors();
    eff_phi->BayesDivide(h_phi_pass.GetPtr(), h_phi_total);

    // Step 6: Draw
    
    // Mass
    TCanvas* c_mass = new TCanvas("c_mass", "Dimuon Mass", 800, 600);
    h_mass_pass->SetLineColor(kGreen + 2);
    h_mass_pass->SetLineWidth(2);
    h_mass_pass->Draw("HIST");
    h_mass_fail->SetLineColor(kRed + 1);
    h_mass_fail->SetLineStyle(2);
    h_mass_fail->SetLineWidth(2);
    h_mass_fail->Draw("HIST SAME");
    auto leg_mass = new TLegend(0.6, 0.7, 0.88, 0.85);
    leg_mass->AddEntry(h_mass_pass.GetPtr(), "Passing probes", "l");
    leg_mass->AddEntry(h_mass_fail.GetPtr(), "Failing probes", "l");
    leg_mass->Draw();
    c_mass->SaveAs("DY_mass_pass_fail.png");

    auto c_eff_mass = new TCanvas("c_eff_mass", "Efficiency vs Dimuon Mass", 800, 600);
    eff_mass->SetMarkerStyle(20);
    eff_mass->SetMarkerColor(kBlue);
    eff_mass->SetTitle("Efficiency vs m_{#mu#mu};m_{#mu#mu} [GeV];Efficiency");
    eff_mass->Draw("AP");

    c_eff_mass->SaveAs("DY_efficiency_mass.png");

    // pT
    TCanvas* c_pt = new TCanvas("c_pt", "Probe pT", 800, 600);
    h_pt_pass->SetLineColor(kGreen + 2);
    h_pt_pass->SetLineWidth(2);
    h_pt_pass->Draw("HIST");
    h_pt_fail->SetLineColor(kRed + 1);
    h_pt_fail->SetLineStyle(2);
    h_pt_fail->SetLineWidth(2);
    h_pt_fail->Draw("HIST SAME");
    auto leg_pt = new TLegend(0.6, 0.7, 0.88, 0.85);
    leg_pt->AddEntry(h_pt_pass.GetPtr(), "Passing probes", "l");
    leg_pt->AddEntry(h_pt_fail.GetPtr(), "Failing probes", "l");
    leg_pt->Draw();
    c_pt->SaveAs("DY_pt_pass_fail.png");

    auto c_eff_pt = new TCanvas("c_eff_pt", "Efficiency vs Dimuon Mass", 800, 600);
    eff_pt->SetMarkerStyle(20);
    eff_pt->SetMarkerColor(kBlue);
    eff_pt->SetTitle("Efficiency vs Probe p_{T};Probe p_{T} [GeV];Efficiency");
    eff_pt->Draw("AP");

    c_eff_pt->SaveAs("DY_efficiency_pt.png");
    
    // Eta
    TCanvas* c_eta = new TCanvas("c_eta", "Probe Eta", 800, 600);
    h_eta_pass->SetLineColor(kGreen + 2);
    h_eta_pass->SetLineWidth(2);
    h_eta_pass->Draw("HIST");
    h_eta_fail->SetLineColor(kRed + 1);
    h_eta_fail->SetLineStyle(2);
    h_eta_fail->SetLineWidth(2);
    h_eta_fail->Draw("HIST SAME");
    auto leg_eta = new TLegend(0.6, 0.7, 0.88, 0.85);
    leg_eta->AddEntry(h_eta_pass.GetPtr(), "Passing probes", "l");
    leg_eta->AddEntry(h_eta_fail.GetPtr(), "Failing probes", "l");
    leg_eta->Draw();
    c_eta->SaveAs("DY_eta_pass_fail.png");

    auto c_eff_eta = new TCanvas("c_eff_eta", "Efficiency vs Dimuon Mass", 800, 600);
    eff_eta->SetMarkerStyle(20);
    eff_eta->SetMarkerColor(kBlue);
    eff_eta->SetTitle("Efficiency vs Probe #eta;Probe #eta;Efficiency");
    eff_eta->Draw("AP");

    c_eff_eta->SaveAs("DY_efficiency_eta.png");
    
    // Phi
    TCanvas* c_phi = new TCanvas("c_phi", "Probe Phi", 800, 600);
    h_phi_pass->SetLineColor(kGreen + 2);
    h_phi_pass->SetLineWidth(2);
    h_phi_pass->Draw("HIST");
    h_phi_fail->SetLineColor(kRed + 1);
    h_phi_fail->SetLineStyle(2);
    h_phi_fail->SetLineWidth(2);
    h_phi_fail->Draw("HIST SAME");
    auto leg_phi = new TLegend(0.6, 0.7, 0.88, 0.85);
    leg_phi->AddEntry(h_phi_pass.GetPtr(), "Passing probes", "l");
    leg_phi->AddEntry(h_phi_fail.GetPtr(), "Failing probes", "l");
    leg_phi->Draw();
    c_phi->SaveAs("DY_phi_pass_fail.png");

    auto c_eff_phi = new TCanvas("c_eff_phi", "Efficiency vs Dimuon Mass", 800, 600);
    eff_phi->SetMarkerStyle(20);
    eff_phi->SetMarkerColor(kBlue);
    eff_phi->SetTitle("Efficiency vs Probe #phi;Probe #phi;Efficiency");
    eff_phi->Draw("AP");

    c_eff_phi->SaveAs("DY_efficiency_phi.png");
   

}
