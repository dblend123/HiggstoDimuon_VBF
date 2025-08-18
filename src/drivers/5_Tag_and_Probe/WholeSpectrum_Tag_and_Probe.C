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

void WholeSpectrum_Tag_and_Probe(const char* inputList, const char* tag) {
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
//auto df2mu = df.Range(10000);

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
             if (Muon_pfRelIso04_all[i] > 0.15) continue;
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

// Step 3: Build PASSING tag-probe pairs
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
)cpp");

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
)cpp");

// Step 5: Plot histograms
    TFile* outputFile = new TFile(("WholeSpectrum_" + baseTag + ".root").c_str(), "RECREATE");
    
    // Pts
    auto h_pt_pass  = df_with_passing_probes.Histo1D({"probe_pt_pass", "Passing Probes; Probe p_{T} [GeV]; Events", 25, 20, 120}, "probe_pt_pass");
    //df_with_passing_probes.Display("probe_pt_pass")->Print();
    auto h_pt_fail  = df_with_passing_probes.Histo1D({"probe_pt_fail", "Failing Probes; Probe p_{T} [GeV]; Events", 25, 20, 120}, "probe_pt_fail");
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
    auto h_phi_pass = df_with_passing_probes.Histo1D({"probe_phi_pass", "Passing Probes; Probe #phi; Events", 20, -3.14, 3.14}, "probe_phi_pass");
    auto h_phi_fail = df_with_passing_probes.Histo1D({"probe_phi_fail", "Failing Probes; Probe #phi; Events", 20, -3.14, 3.14}, "probe_phi_fail");
    // Clone passing phi histogram
    auto h_phi_total = (TH1D*) h_phi_pass.GetPtr()->Clone("h_phi_total");
    h_phi_total->Add(h_phi_fail.GetPtr());
    // Create efficiency graph
    auto eff_phi = new TGraphAsymmErrors();
    eff_phi->BayesDivide(h_phi_pass.GetPtr(), h_phi_total);

    // Step 6: Draw
    // pT
    TCanvas* c_pt = new TCanvas("c_pt", "Probe pT", 800, 600);
    h_pt_pass->SetLineColor(kGreen + 2);
    h_pt_pass->SetLineWidth(2);
    h_pt_fail->SetLineColor(kRed + 1);
    h_pt_fail->SetLineStyle(2);
    h_pt_fail->SetLineWidth(2);

    // Normalize histograms for shape comparison (optional)
    h_pt_pass->Scale(1.0 / h_pt_pass->Integral());
    h_pt_fail->Scale(1.0 / h_pt_fail->Integral());
    h_pt_pass->SetTitle("Normalized Probe p_{T};muon p_{T} (GeV);Normalized Events");

    // Adjust axis
    h_pt_pass->SetMaximum(std::max(h_pt_pass->GetMaximum(), h_pt_fail->GetMaximum()) * 1.2);
    gStyle->SetOptStat(0);

    // Draw 
    h_pt_pass->Draw("HIST");
    h_pt_fail->Draw("HIST SAME");
    
    auto leg_pt = new TLegend(0.65, 0.75, 0.88, 0.9);
    leg_pt->AddEntry(h_pt_pass.GetPtr(), "Passing probes", "l");
    leg_pt->AddEntry(h_pt_fail.GetPtr(), "Failing probes", "l");
    leg_pt->Draw();
    
    c_pt->SaveAs(("DY_pt_pass_fail_" + baseTag + ".pdf").c_str());

    auto c_eff_pt = new TCanvas("c_eff_pt", "Efficiency vs Probe p_{T}", 800, 600);
    eff_pt->SetMarkerStyle(20);
    eff_pt->SetMarkerColor(kBlue);
    eff_pt->SetTitle("Efficiency vs Probe p_{T};muon p_{T} (GeV);Efficiency");
    eff_pt->SetMinimum(0.85);
    eff_pt->SetMaximum(1.1);
    eff_pt->Draw("AP");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.18, 0.82, "#font[62]{CMS}");
    latex.DrawLatex(0.18, 0.78, "#it{Simulation}");
    latex.DrawLatex(0.18, 0.74, "#it{Preliminary}");
    latex.DrawLatex(0.65, 0.82, "Drell-Yan");
    latex.DrawLatex(0.55, 0.78, "20.0 GeV < muon p_{T} < 120.0 GeV");
    latex.DrawLatex(0.65, 0.74, "-2.4 < muon  #eta < 2.4");
    latex.DrawLatex(0.65, 0.70, "-#pi < muon  #phi < #pi");
    
    eff_pt->SetName(("eff_pt_" + baseTag).c_str());
    c_eff_pt->SaveAs(("DY_efficiency_pt_" + baseTag + ".pdf").c_str());
    
    h_pt_pass->Write();
    h_pt_fail->Write(); 
    h_pt_total->Write();
    leg_pt->Write();
    eff_pt->Write();

    // Eta
    TCanvas* c_eta = new TCanvas("c_eta", "Probe Eta", 800, 600);
    h_eta_pass->SetLineColor(kGreen + 2);
    h_eta_pass->SetLineWidth(2);
    h_eta_fail->SetLineColor(kRed + 1);
    h_eta_fail->SetLineStyle(2);
    h_eta_fail->SetLineWidth(2);
    
    // Normalize histograms for shape comparison (optional)
    h_eta_pass->Scale(1.0 / h_eta_pass->Integral());
    h_eta_fail->Scale(1.0 / h_eta_fail->Integral());
    h_eta_pass->SetTitle("Normalized Probe  #eta;muon  #eta;Normalized Events");

    // Adjust axis
    h_eta_pass->SetMaximum(std::max(h_eta_pass->GetMaximum(), h_eta_fail->GetMaximum()) * 1.2);
    gStyle->SetOptStat(0);

    // Draw
    h_eta_pass->Draw("HIST");
    h_eta_fail->Draw("HIST SAME");

    auto leg_eta = new TLegend(0.65, 0.75, 0.88, 0.9);
    leg_eta->AddEntry(h_eta_pass.GetPtr(), "Passing probes", "l");
    leg_eta->AddEntry(h_eta_fail.GetPtr(), "Failing probes", "l");
    leg_eta->Draw();
    c_eta->SaveAs(("DY_eta_pass_fail_" + baseTag + ".pdf").c_str());

    auto c_eff_eta = new TCanvas("c_eff_eta", "Efficiency vs Probe  #eta", 800, 600);
    eff_eta->SetMarkerStyle(20);
    eff_eta->SetMarkerColor(kBlue);
    eff_eta->SetTitle("Efficiency vs Probe  #eta;muon  #eta;Efficiency");
    eff_eta->SetMinimum(0.85);
    eff_eta->SetMaximum(1.1);
    eff_eta->Draw("AP");

    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.18, 0.82, "#font[62]{CMS}");
    latex.DrawLatex(0.18, 0.78, "#it{Simulation}");
    latex.DrawLatex(0.18, 0.74, "#it{Preliminary}");
    latex.DrawLatex(0.65, 0.82, "Drell-Yan");
    latex.DrawLatex(0.55, 0.78, "20.0 GeV < muon p_{T} < 120.0 GeV");
    latex.DrawLatex(0.65, 0.74, "-2.4 < muon  #eta < 2.4");
    latex.DrawLatex(0.65, 0.70, "-#pi < muon  #phi < #pi");
    
    
    eff_eta->SetName(("eff_eta_" + baseTag).c_str());
    c_eff_eta->SaveAs(("DY_efficiency_eta_" + baseTag + ".pdf").c_str());
    
    h_eta_pass->Write();
    h_eta_fail->Write(); 
    h_eta_total->Write();
    leg_eta->Write();
    eff_eta->Write();
    
    // Phi
    TCanvas* c_phi = new TCanvas("c_phi", "Probe Phi", 800, 600);
    h_phi_pass->SetLineColor(kGreen + 2);
    h_phi_pass->SetLineWidth(2);
    h_phi_fail->SetLineColor(kRed + 1);
    h_phi_fail->SetLineStyle(2);
    h_phi_fail->SetLineWidth(2);
    
    // Normalize histograms for shape comparison (optional)
    h_phi_pass->Scale(1.0 / h_phi_pass->Integral());
    h_phi_fail->Scale(1.0 / h_phi_fail->Integral());
    h_phi_pass->SetTitle("Normalized Probe  #phi;Probe  #phi;Normalized Events");

    // Adjust axis
    h_phi_pass->SetMaximum(std::max(h_phi_pass->GetMaximum(), h_phi_fail->GetMaximum()) * 1.2);
    gStyle->SetOptStat(0);
    
    // Draw
    h_phi_pass->Draw("HIST");
    h_phi_fail->Draw("HIST SAME");

    auto leg_phi = new TLegend(0.65, 0.75, 0.88, 0.9);
    leg_phi->AddEntry(h_phi_pass.GetPtr(), "Passing probes", "l");
    leg_phi->AddEntry(h_phi_fail.GetPtr(), "Failing probes", "l");
    leg_phi->Draw();
    c_phi->SaveAs(("DY_phi_pass_fail_" + baseTag + ".pdf").c_str());

    auto c_eff_phi = new TCanvas("c_eff_phi", "Efficiency vs Probe #phi", 800, 600);
    eff_phi->SetMarkerStyle(20);
    eff_phi->SetMarkerColor(kBlue);
    eff_phi->SetTitle("Efficiency vs Probe  #phi;muon  #phi;Efficiency");
    eff_phi->SetMinimum(0.85);
    eff_phi->SetMaximum(1.1);
    eff_phi->Draw("AP");

    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.18, 0.82, "#font[62]{CMS}");
    latex.DrawLatex(0.18, 0.78, "#it{Simulation}");
    latex.DrawLatex(0.18, 0.74, "#it{Preliminary}");
    latex.DrawLatex(0.65, 0.82, "Drell-Yan");
    latex.DrawLatex(0.55, 0.78, "20.0 GeV < muon p_{T} < 120.0 GeV");
    latex.DrawLatex(0.65, 0.74, "-2.4 < muon  #eta < 2.4");
    latex.DrawLatex(0.65, 0.70, "-#pi < muon  #phi < #pi");
    
    eff_phi->SetName(("eff_phi_" + baseTag).c_str());
    c_eff_phi->SaveAs(("DY_efficiency_phi_" + baseTag + ".pdf").c_str());
    
    h_phi_pass->Write();
    h_phi_fail->Write(); 
    h_phi_total->Write();
    leg_phi->Write();
    eff_phi->Write();

    outputFile->Close();
}
