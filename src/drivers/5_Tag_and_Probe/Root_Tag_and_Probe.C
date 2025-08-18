#include <fstream>  
#include <iostream>
#include <vector>
#include <string>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <ROOT/RDF/HistoModels.hxx>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TLatex.h>

using namespace ROOT;
using namespace ROOT::VecOps;
using namespace ROOT::RDF;

void Root_Tag_and_Probe(const char* inputList, const char* tag) {
    std::string baseTag(tag);  // Use the passed tag directly

    std::ifstream infile(inputList);
    std::vector<std::string> fileList;
    std::string line;

    while (std::getline(infile, line)) {
        if (!line.empty()) fileList.push_back(line);
    }

    ROOT::RDataFrame df("Events", fileList);

// Step 1: Only keep events with ≥2 muons
auto df2mu = df.Filter("nMuon >= 2");

//std::cout << "nMuons: " << df2mu.Count().GetValue() << std::endl;
auto totalMuons = df2mu.Sum("nMuon");
std::cout << "Total muons in events with >=2 muons: " << totalMuons.GetValue() << std::endl;

// Step 2: Define tag and probe masks  
// The Tag Muon: muons that pass tight selection
auto df_with_tags = df2mu.Define("tag_mask", R"cpp(
     std::vector<bool> mask(Muon_pt.size(), false);
     float max_pt = -1.0;
     int best_idx = -1;
     if (HLT_IsoMu24) {
         for (size_t i = 0; i < Muon_pt.size(); ++i) {
             if (!Muon_tightId[i]) continue;
             if (Muon_pt[i] < 29) continue;
             if (TMath::Abs(Muon_eta[i]) >= 2.4) continue;
             if (Muon_pfRelIso04_all[i] >= 0.15) continue;
             if (Muon_isGlobal[i] == false) continue; 
             if (Muon_isTracker[i] == false) continue; 

             // Trigger Matching
             bool matched = false;
             for (size_t j = 0; j < TrigObj_pt.size(); ++j) {
                 // 13 = Muon
                 if (TrigObj_id[j] != 13) continue;
                 // Only consider triggers with  pT > 24 GeV (Mimics Single-Muon HLT --> HLT_IsoMu24)
                 if (TrigObj_pt[j] < 24) continue;

                 // You calculate ΔR² between the reco muon and the trigger object
                 // ΔR² = Δη² + Δφ² is a common way to match objects spatially
                 float deta = Muon_eta[i] - TrigObj_eta[j];
                 float dphi = std::abs(Muon_phi[i] - TrigObj_phi[j]);
                 if (dphi > M_PI) dphi = 2 * M_PI - dphi;
                 float dr2 = deta*deta + dphi*dphi;
                 //If the angular distance ΔR < 0.1, you consider this muon to be matched to a trigger object (It could have fired the trigger)
                 if (dr2 < 0.01) { // dR < 0.1
                     matched = true;
                     break;
                 }
             }
             // Finally, if any trigger object matched this muon, you mark this muon as a tag
             if (matched && Muon_pt[i] > max_pt) {
             max_pt = Muon_pt[i];
             best_idx = i;
             //if (matched) mask[i] = true;
             }
         }
     }
     if (best_idx >= 0) mask[best_idx] = true;
     return mask;
     )cpp");


//auto count_tag_quality = df_with_tags.Filter("Sum(tag_mask) > 0").Count();
//std::cout << "After tag_mask: " << count_tag_quality.GetValue() << std::endl;
auto nTags = df_with_tags.Define("nTags", "Sum(tag_mask)").Sum("nTags");
std::cout << "Total tag muons: " << nTags.GetValue() << std::endl;

// The Probe Muon: all other muons in event
auto df_with_probes = df_with_tags.Define("probe_mask", R"cpp(
    std::vector<bool> mask;
    for (bool is_tag : tag_mask) {
        mask.push_back(!is_tag);
    }
    return mask;
    )cpp");

//auto count_probe_quality = df_with_probes.Filter("Sum(probe_mask) > 0").Count();
//std::cout << "After probe_mask: " << count_probe_quality.GetValue() << std::endl;
auto nProbes = df_with_probes.Define("nProbes", "Sum(probe_mask)").Sum("nProbes");
std::cout << "Total probe muons: " << nProbes.GetValue() << std::endl;

// Any standalone muon updated at vertex with more than one hit and more than one station in the muon system (i.e. good track-hit χ2).
auto df_with_good_probes = df_with_probes.Define("probe_quality_mask", R"cpp(
    std::vector<bool> mask(Muon_pt.size(), false);
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!probe_mask[i]) continue;
        if (Muon_pt[i] < 20) continue;
        if (TMath::Abs(Muon_eta[i]) >= 2.4) continue;

        // Only Standalone muons updated at vertex
        //if (Muon_isStandalone[i] == false) continue; 
        if (Muon_isTracker[i] == false) continue; //If you're measuring general tightId efficiency, you'd normally want tracker muons or any reco muon.

        // Require at least two muon stations
        if (Muon_nStations[i] <= 1) continue;

        // Optional: require minimum tracker layers
        // if (Muon_nTrackerLayers[i] <= 5) continue;

        // Optional: use Muon_highPtId (0=fail, 1=tracker high pT (loose), 2=global high pT, which includes tracker high pT (tight))
        //if (Muon_highPtId[i] < 1) continue;
        
        mask[i] = true; // Only set to true if all cuts are passed
    }
    return mask;
)cpp");

//auto count_good_probe_quality = df_with_good_probes.Filter("Sum(probe_quality_mask) > 0").Count();
//std::cout << "After probe_quality_mask: " << count_good_probe_quality.GetValue() << std::endl;
auto nGoodProbes = df_with_good_probes.Define("nGoodProbes", "Sum(probe_quality_mask)").Sum("nGoodProbes");
std::cout << "Total good probe muons: " << nGoodProbes.GetValue() << std::endl;

auto df_with_passing_probes = df_with_good_probes.Define("passing_probe_mask", R"cpp(
    std::vector<bool> mask(Muon_pt.size(), false);
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!probe_quality_mask[i]) continue;
        if (!Muon_tightId[i]) continue;

        mask[i] = true; 
    }
    return mask;
)cpp");

//auto count_passing_probe_quality = df_with_passing_probes.Filter("Sum(passing_probe_mask) > 0").Count();
//std::cout << "After passing_probe_mask: " << count_passing_probe_quality.GetValue() << std::endl;
auto nPassingProbes = df_with_passing_probes.Define("nPassingProbes", "Sum(passing_probe_mask)").Sum("nPassingProbes");
std::cout << "Total passing probe muons: " << nPassingProbes.GetValue() << std::endl;

// Step 3: Build invariant mass of PASSING tag-probe pairs
// Passing pt
df_with_passing_probes = df_with_passing_probes.Define("probe_pt_pass", R"cpp(
    std::vector<double> pts;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
            }
        }
        if (local_max_pt > 0) {
            pts.push_back(local_max_pt);
        //break; // optional: break if you want only one tag-probe pair per event
        }
    }
    return pts;
)cpp")

// Passing eta
.Define("probe_eta_pass", R"cpp(
    std::vector<double> etas;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        int best_probe_idx = -1;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
                best_probe_idx = j;
            }
        }
        if (local_max_pt > 0) {
            etas.push_back(Muon_eta[best_probe_idx]);
        //break; // optional: break if you want only one tag-probe pair per event
        }
    }
    return etas;
)cpp")

// Passing phi
.Define("probe_phi_pass", R"cpp(
    std::vector<double> phis;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        int best_probe_idx = -1;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
                best_probe_idx = j;
            }
        }
        if (local_max_pt > 0) {
            phis.push_back(Muon_phi[best_probe_idx]);
        //break; // optional: break if you want only one tag-probe pair per event
        }
    }
    return phis;
)cpp")

// Passing dimuon mass
.Define("probe_mass_pass", R"cpp(
    std::vector<double> masses;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        int best_probe_idx = -1;
        double best_mass = -1.0;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j] || !passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
                best_probe_idx = j;
                best_mass = m;
            }
        }
        if (local_max_pt > 0 && best_mass > 0) {
            masses.push_back(best_mass);
        //break; // optional: break if you want only one tag-probe pair per event
        }
    }
    return masses;
)cpp");

auto df_flat_pass = df_with_passing_probes
    .Filter("probe_pt_pass.size() > 0 && probe_eta_pass.size() > 0 && probe_phi_pass.size() > 0 && probe_mass_pass.size() > 0")
    .Define("flat_probe_pt_pass", "probe_pt_pass[0]")
    .Define("flat_probe_eta_pass", "probe_eta_pass[0]")
    .Define("flat_probe_phi_pass", "probe_phi_pass[0]")
    .Define("flat_probe_mass_pass", "probe_mass_pass[0]");
/*
auto h2_pt_eta_pass = df_flat_pass
    .Filter("flat_probe_pt_pass >= 0 && flat_probe_eta_pass >= -2.4 && flat_probe_eta_pass <= 2.4")
    .Histo2D({"h2_pt_eta_pass", "Passing Probes;#eta; p_{T} [GeV]", 24, -2.4, 2.4, 30, 0, 150}, "flat_probe_eta_pass", "flat_probe_pt_pass");
*/
// Step 4: Build invariant mass of FAILING tag-probe pairs
// Failing pt
df_with_passing_probes = df_with_passing_probes.Define("probe_pt_fail", R"cpp(
    std::vector<double> pts;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
            }
        }
        if (local_max_pt > 0) {
            pts.push_back(local_max_pt);
        //break; // Optional: remove this if you want to allow multiple tags
        }
    }
    return pts;
)cpp")

// Failing eta
.Define("probe_eta_fail", R"cpp(
    std::vector<double> etas;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        int best_probe_idx = -1;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
                best_probe_idx = j;
            }
        }
        if (local_max_pt > 0) {
            etas.push_back(Muon_eta[best_probe_idx]);
        //break; // Optional: remove this if you want to allow multiple tags
        } 
    }
    return etas;
)cpp")

// Failing phi
.Define("probe_phi_fail", R"cpp(
    std::vector<double> phis;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        int best_probe_idx = -1;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
                best_probe_idx = j;
            }
        }
        if (local_max_pt > 0) {
            phis.push_back(Muon_phi[best_probe_idx]);
            //break; // Optional: remove this if you want to allow multiple tags
        }
    }
    return phis;
)cpp")

// Failing dimuon mass
.Define("probe_mass_fail", R"cpp(
    std::vector<double> masses;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        double local_max_pt = -1.0;
        int best_probe_idx = -1;
        double best_mass = -1.0;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            if (!probe_quality_mask[j]) continue;
            if (passing_probe_mask[j]) continue;

            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 70 && m < 115 && Muon_pt[j] > local_max_pt) {
                local_max_pt = Muon_pt[j];
                best_probe_idx = j;
                best_mass = m;
            }
        }
        if (local_max_pt > 0 && best_mass > 0) {
            masses.push_back(best_mass);
        //break; // Optional: remove this if you want to allow multiple tags
        }
    }
    return masses;
)cpp");

auto df_flat_fail = df_with_passing_probes
    .Filter("probe_pt_fail.size() > 0 && probe_eta_fail.size() > 0 && probe_phi_fail.size() > 0 && probe_mass_fail.size() > 0")
    .Define("flat_probe_pt_fail", "probe_pt_fail[0]")
    .Define("flat_probe_eta_fail", "probe_eta_fail[0]")
    .Define("flat_probe_phi_fail", "probe_phi_fail[0]")
    .Define("flat_probe_mass_fail", "probe_mass_fail[0]");
/*
auto h2_pt_eta_fail = df_flat_fail
    .Filter("flat_probe_pt_fail >= 0 && flat_probe_eta_fail >= -2.4 && flat_probe_eta_fail <= 2.4")
    .Histo2D({"h2_pt_eta_fail", "Failing Probes;#eta; p_{T} [GeV]", 24, -2.4, 2.4, 30, 0, 150}, "flat_probe_eta_fail", "flat_probe_pt_fail");
*/
//std::vector<float> eta_bins = {-2.4, -2.1, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.1, 2.4};
std::vector<float> eta_bins = {-2.4, -1.6, -0.9, 0.9, 1.6, 2.4};
std::vector<float> pt_bins  = {20, 25, 30, 35, 40, 50, 60, 80, 120};
//std::vector<float> phi_bins = {-3.14, -2.4, -1.6, -0.8, 0.0, 0.8, 1.6, 2.4, 3.14};
std::vector<float> phi_bins = {-3.1416, -1.571, 0.0, 1.571, 3.1416};

// Step 5: Plot histograms
TFile* outFile = new TFile(("TagAndProbe_DimuonInvMass_" + baseTag + ".root").c_str(), "RECREATE");

// Etas
for (size_t i = 0; i < eta_bins.size() - 1; ++i) {    
    float eta_min = eta_bins[i];
    float eta_max = eta_bins[i+1];
    std::string bin_label = "eta_" + std::to_string(i);
   
    auto h_eta_pass = df_flat_pass.Filter(Form("flat_probe_eta_pass > %.3f && flat_probe_eta_pass <= %.3f", eta_min, eta_max),bin_label + "_pass");
    auto h_eta_fail = df_flat_fail.Filter(Form("flat_probe_eta_fail > %.3f && flat_probe_eta_fail <= %.3f", eta_min, eta_max), bin_label + "_fail");

    auto h_mass_pass = h_eta_pass.Histo1D(TH1DModel(("h_mass_pass_" + bin_label).c_str(), Form("Passing Probes (%.1f < #eta < %.1f); m_{#mu#mu} [GeV]; Events", eta_min, eta_max), 45, 70, 115), "flat_probe_mass_pass");
    auto h_mass_fail = h_eta_fail.Histo1D(TH1DModel(("h_mass_fail_" + bin_label).c_str(), Form("Passing Probes (%.1f < #eta < %.1f); m_{#mu#mu} [GeV]; Events", eta_min, eta_max), 45, 70, 115), "flat_probe_mass_fail");

    // Canvas for passing probes
    TCanvas* c_eta_pass = new TCanvas(("c_eta_pass_" + bin_label).c_str(), ("Pass Mass " + bin_label).c_str(), 800, 600);
    h_mass_pass->SetMarkerStyle(20);      
    h_mass_pass->SetMarkerSize(1.2);     
    h_mass_pass->SetMarkerColor(kBlack);
    h_mass_pass->SetTitle(Form("Passing Probes (%.1f < #eta < %.1f);m_{#mu#mu} [GeV];Events", eta_min, eta_max));
    h_mass_pass->Draw("PE");

    // Canvas for failing probes
    TCanvas* c_eta_fail = new TCanvas(("c_eta_fail_" + bin_label).c_str(), ("Fail Mass " + bin_label).c_str(), 800, 600);
    h_mass_fail->SetMarkerStyle(20);      
    h_mass_fail->SetMarkerSize(1.2);     
    h_mass_fail->SetMarkerColor(kBlack);
    h_mass_fail->SetTitle(Form("Failing Probes (%.1f < #eta < %.1f);m_{#mu#mu} [GeV];Events", eta_min, eta_max));
    h_mass_fail->Draw("PE");
    
    h_mass_pass->Write();
    h_mass_fail->Write();
}

// Pts
for (size_t i = 0; i < pt_bins.size() - 1; ++i) {    
    float pt_min = pt_bins[i];
    float pt_max = pt_bins[i+1];
    std::string bin_label = "pt_" + std::to_string(i);
   
    auto h_pt_pass = df_flat_pass.Filter(Form("flat_probe_pt_pass > %.3f && flat_probe_pt_pass <= %.3f", pt_min, pt_max),bin_label + "_pass");
    auto h_pt_fail = df_flat_fail.Filter(Form("flat_probe_pt_fail > %.3f && flat_probe_pt_fail <= %.3f", pt_min, pt_max), bin_label + "_fail");

    auto h_mass_pass = h_pt_pass.Histo1D(TH1DModel(("h_mass_pass_" + bin_label).c_str(), Form("Passing Probes (%.1f < p_{T} < %.1f); m_{#mu#mu} [GeV]; Events", pt_min, pt_max), 45, 70, 115), "flat_probe_mass_pass");
    auto h_mass_fail = h_pt_fail.Histo1D(TH1DModel(("h_mass_fail_" + bin_label).c_str(), Form("Passing Probes (%.1f < p_{T} < %.1f); m_{#mu#mu} [GeV]; Events", pt_min, pt_max), 45, 70, 115), "flat_probe_mass_fail");

    // Canvas for passing probes
    TCanvas* c_pt_pass = new TCanvas(("c_pt_pass_" + bin_label).c_str(), ("Pass Mass " + bin_label).c_str(), 800, 600);
    h_mass_pass->SetMarkerStyle(20);      
    h_mass_pass->SetMarkerSize(1.2);     
    h_mass_pass->SetMarkerColor(kBlack);
    h_mass_pass->SetTitle(Form("Passing Probes (%.1f < p_{T} < %.1f);m_{#mu#mu} [GeV];Events", pt_min, pt_max));
    h_mass_pass->Draw("PE");

    // Canvas for failing probes
    TCanvas* c_pt_fail = new TCanvas(("c_pt_fail_" + bin_label).c_str(), ("Fail Mass " + bin_label).c_str(), 800, 600);
    h_mass_fail->SetMarkerStyle(20);      
    h_mass_fail->SetMarkerSize(1.2);     
    h_mass_fail->SetMarkerColor(kBlack);
    h_mass_fail->SetTitle(Form("Failing Probes (%.1f < p_{T} < %.1f);m_{#mu#mu} [GeV];Events", pt_min, pt_max));
    h_mass_fail->Draw("PE");

    h_mass_pass->Write();
    h_mass_fail->Write();
}

// Phis
for (size_t i = 0; i < phi_bins.size() - 1; ++i) {    
    float phi_min = phi_bins[i];
    float phi_max = phi_bins[i+1];
    std::string bin_label = "phi_" + std::to_string(i);
   
    auto h_phi_pass = df_flat_pass.Filter(Form("flat_probe_phi_pass > %.3f && flat_probe_phi_pass <= %.3f", phi_min, phi_max),bin_label + "_pass");
    auto h_phi_fail = df_flat_fail.Filter(Form("flat_probe_phi_fail > %.3f && flat_probe_phi_fail <= %.3f", phi_min, phi_max), bin_label + "_fail");

    auto h_mass_pass = h_phi_pass.Histo1D(TH1DModel(("h_mass_pass_" + bin_label).c_str(), Form("Passing Probes (%.1f < #phi < %.1f); m_{#mu#mu} [GeV]; Events", phi_min, phi_max), 45, 70, 115), "flat_probe_mass_pass");
    auto h_mass_fail = h_phi_fail.Histo1D(TH1DModel(("h_mass_fail_" + bin_label).c_str(), Form("Passing Probes (%.1f < #phi < %.1f); m_{#mu#mu} [GeV]; Events", phi_min, phi_max), 45, 70, 115), "flat_probe_mass_fail");

    // Canvas for passing probes
    TCanvas* c_phi_pass = new TCanvas(("c_phi_pass_" + bin_label).c_str(), ("Pass Mass " + bin_label).c_str(), 800, 600);
    h_mass_pass->SetMarkerStyle(20);      
    h_mass_pass->SetMarkerSize(1.2);     
    h_mass_pass->SetMarkerColor(kBlack);
    h_mass_pass->SetTitle(Form("Passing Probes (%.1f < #phi < %.1f);m_{#mu#mu} [GeV];Events", phi_min, phi_max));
    h_mass_pass->Draw("PE");

    TCanvas* c_phi_fail = new TCanvas(("c_phi_fail_" + bin_label).c_str(), ("Fail Mass " + bin_label).c_str(), 800, 600);
    h_mass_fail->SetMarkerStyle(20);      
    h_mass_fail->SetMarkerSize(1.2);     
    h_mass_fail->SetMarkerColor(kBlack);
    h_mass_fail->SetTitle(Form("Failing Probes (%.1f < #phi < %.1f);m_{#mu#mu} [GeV];Events", phi_min, phi_max));
    h_mass_fail->Draw("PE");

    h_mass_pass->Write();
    h_mass_fail->Write();
}

// Pts vs Etas
TFile* outfile = new TFile(("TagAndProbe_Heatmaps_" + baseTag + ".root").c_str(), "RECREATE");

// Define new binning
//const int nEtaBins = 12;
//const double etaBins[nEtaBins + 1] = {-2.4, -2.1, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.1, 2.4};
const int nEtaBins = 5;
const double etaBins[nEtaBins + 1] = {-2.4, -1.6, -0.9, 0.9, 1.6, 2.4};
const int nPtBins = 8;
const double ptBins[nPtBins + 1] = {20, 25, 30, 35, 40, 50, 60, 80, 120};

auto h2_pt_eta_pass = df_flat_pass.Histo2D({"h2_pt_eta_pass", "Passing Probes;#eta;p_{T} [GeV]", nEtaBins, etaBins, nPtBins, ptBins}, "flat_probe_eta_pass", "flat_probe_pt_pass");
auto h2_pt_eta_fail = df_flat_fail.Histo2D({"h2_pt_eta_fail", "Failing Probes;#eta;p_{T} [GeV]", nEtaBins, etaBins, nPtBins, ptBins}, "flat_probe_eta_fail", "flat_probe_pt_fail");

// Create re-binned versions
/*TH2D* h2_pt_eta_pass_rebinned = new TH2D("h2_pt_eta_pass_rebinned", "Passing Probes;#eta;p_{T} [GeV]", nEtaBins, etaBins, nPtBins, ptBins);
TH2D* h2_pt_eta_fail_rebinned = new TH2D("h2_pt_eta_fail_rebinned", "Failing Probes;#eta;p_{T} [GeV]", nEtaBins, etaBins, nPtBins, ptBins);

// Fill new histograms with existing contents
for (int x = 1; x <= h2_pt_eta_pass->GetNbinsX(); ++x) {
    for (int y = 1; y <= h2_pt_eta_pass->GetNbinsY(); ++y) {
        double content_pass = h2_pt_eta_pass->GetBinContent(x, y);
        double content_fail = h2_pt_eta_fail->GetBinContent(x, y);
        double eta = h2_pt_eta_pass->GetXaxis()->GetBinCenter(x);
        double pt = h2_pt_eta_pass->GetYaxis()->GetBinCenter(y);
        h2_pt_eta_pass_rebinned->Fill(eta, pt, content_pass);
        h2_pt_eta_fail_rebinned->Fill(eta, pt, content_fail);
    }
}
*/
// Make total and efficiency histograms
TH2D* h2_pt_eta_total = (TH2D*) h2_pt_eta_pass->Clone("h2_pt_eta_total");
h2_pt_eta_total->Add(h2_pt_eta_fail.GetPtr());

// Draw passing heatmap
TCanvas* c_pass = new TCanvas("c_pass", "Passing Probes pt vs eta", 800, 600);
h2_pt_eta_pass->SetTitle("Passing Probes;#eta;p_{T} [GeV]");
h2_pt_eta_pass->Draw("COLZ");

// Draw failing heatmap
TCanvas* c_fail = new TCanvas("c_fail", "Failing Probes pt vs eta", 800, 600);
h2_pt_eta_fail->SetTitle("Failing Probes;#eta;p_{T} [GeV]");
h2_pt_eta_fail->Draw("COLZ");

// Draw efficiency heatmap
TH2D* h2_pt_eta_eff = (TH2D*) h2_pt_eta_pass->Clone("h2_pt_eta_eff");
h2_pt_eta_eff->Divide(h2_pt_eta_pass.GetPtr(), h2_pt_eta_total);
h2_pt_eta_eff->SetTitle("DY-Run22 p_{T} vs. #eta Efficiency Map;#eta;p_{T} [GeV]");
h2_pt_eta_eff->SetMinimum(0.5);
h2_pt_eta_eff->SetMaximum(1.0);

TCanvas* c_eff = new TCanvas("c_eff", "Efficiency pt vs eta", 800, 600);
gStyle->SetOptStat(0);
h2_pt_eta_eff->Draw("COLZ");

TLatex latex;
latex.SetNDC();
latex.SetTextFont(42);
latex.SetTextSize(0.04);
latex.DrawLatex(0.13, 0.91, "CMS Preliminary");
latex.DrawLatex(0.62, 0.91, "Simulation (Tight ID)");

h2_pt_eta_pass->Write("h2_pt_eta_pass");
h2_pt_eta_fail->Write("h2_pt_eta_fail");
h2_pt_eta_eff->Write("h2_pt_eta_eff");
outfile->Close();
}
