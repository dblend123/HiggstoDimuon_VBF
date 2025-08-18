#include <fstream>  
#include <iostream>
#include <vector>
#include <string>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace ROOT;
using namespace ROOT::VecOps;

void Tag_and_Probe() {
    std::ifstream infile("test.list");
    std::vector<std::string> fileList;
    std::string line;

    while (std::getline(infile, line)) {
        if (!line.empty()) fileList.push_back(line);
    }

    ROOT::RDataFrame df("Events", fileList);

// Step 1: Only keep events with ≥2 muons
auto df2mu = df.Filter("nMuon >= 2");

// Step 2: Define tag and probe masks  
// Tags: muons that pass tight selection
auto df_with_tags = df2mu.Define("tag_mask", "Muon_pt > 29 && abs(Muon_eta) < 2.4 && Muon_isGlobal && Muon_isTracker && Muon_tightId");
// Probes: all other muons in event
auto df_with_probes = df_with_tags.Define("probe_mask", "!tag_mask");

// Step 3: Build invariant mass of passing tag-probe pairs
// Passing Probes
df_with_probes = df_with_probes.Define("probe_variables_pass", R"cpp(
    std::vector<double> pt, eta, phi, mass;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            
            if (m > 60 && m < 120 && Muon_tightId[j]) {
                pt.push_back(Muon_pt[j]);
                eta.push_back(Muon_eta[j]);
                phi.push_back(Muon_phi[j]);
                mass.push_back(m);
            }
        }
    }
    return std::make_tuple(pt, eta, phi, mass);;
)cpp");

// Step 4: Build invariant mass of failing tag-probe pairs
// Failing Probes
df_with_probes = df_with_probes.Define("probe_variables_fail", R"cpp(
    std::vector<double> pt, eta, phi, mass;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            
            if (m > 60 && m < 120 && !Muon_tightId[j]) {
                pt.push_back(Muon_pt[j]);
                eta.push_back(Muon_eta[j]);
                phi.push_back(Muon_phi[j]);
                mass.push_back(m);
            }
        }
    }
    return std::make_tuple(pt, eta, phi, mass);;
)cpp");

    df_with_probes = df_with_probes
    .Define("probe_pt_pass",  "std::get<0>(probe_variables_pass)")
    .Define("probe_eta_pass", "std::get<1>(probe_variables_pass)")
    .Define("probe_phi_pass", "std::get<2>(probe_variables_pass)")
    .Define("probe_mass_pass","std::get<3>(probe_variables_pass)");

// 5. Define pt bins
std::vector<float> pt_bins = {0, 10, 20, 30, 40, 50, 60, 80, 100};
const int n_pt_bins = pt_bins.size() - 1;

// 6. Create vectors to hold histograms
std::vector<ROOT::RDF::RResultPtr<TH1D>> h_mass_pass_bins;
std::vector<ROOT::RDF::RResultPtr<TH1D>> h_mass_fail_bins;

// 7. Fill histograms per bin
for (int i = 0; i < n_pt_bins; ++i) {
    std::string bin_sel_pass = Form("Sum((probe_pt_pass >= %.1f) && (probe_pt_pass < %.1f)) > 0", pt_bins[i], pt_bins[i+1]);
    std::string bin_sel_fail = Form("Sum((probe_pt_fail >= %.1f) && (probe_pt_fail < %.1f)) > 0", pt_bins[i], pt_bins[i+1]);
    
    h_mass_pass_bins.push_back(df_with_probes.Filter(bin_sel_pass).Histo1D({"h_pass", "", 60, 60, 120}, "probe_mass_pass"));
    h_mass_fail_bins.push_back(df_with_probes.Filter(bin_sel_fail).Histo1D({"h_fail", "", 60, 60, 120}, "probe_mass_fail"));
}

// 8. Create vectors to hold efficiencies
std::vector<double> efficiencies, eff_errors;

// 9. Loop over bins and fit
for (int i = 0; i < n_pt_bins; ++i) {
    // Passing probes
    auto h_pass = h_mass_pass_bins[i];
    TF1* fit_pass = new TF1(Form("fit_pass_bin%d", i), "gaus(0) + expo(3)", 60, 120);
    h_pass->Fit(fit_pass, "RQ0"); // Quiet, No draw
    double sig_amp_pass = fit_pass->GetParameter(0);
    double sig_sigma_pass = fit_pass->GetParameter(2);
    double sig_yield_pass = sig_amp_pass * sig_sigma_pass * std::sqrt(2 * M_PI);

    // Failing probes
    auto h_fail = h_mass_fail_bins[i];
    TF1* fit_fail = new TF1(Form("fit_fail_bin%d", i), "gaus(0) + expo(3)", 60, 120);
    h_fail->Fit(fit_fail, "RQ0");
    double sig_amp_fail = fit_fail->GetParameter(0);
    double sig_sigma_fail = fit_fail->GetParameter(2);
    double sig_yield_fail = sig_amp_fail * sig_sigma_fail * std::sqrt(2 * M_PI);

    // Compute efficiency
    double denom = sig_yield_pass + sig_yield_fail;
    double eff = (denom > 0) ? sig_yield_pass / denom : 0.0;
    efficiencies.push_back(eff);

    // Simple error (binomial)
    double err = (denom > 0) ? std::sqrt(eff * (1 - eff) / denom) : 0.0;
    eff_errors.push_back(err);
}

// 10. Create efficiency plot
auto g_eff = new TGraphErrors(n_pt_bins);
for (int i = 0; i < n_pt_bins; ++i) {
    double pt_center = 0.5 * (pt_bins[i] + pt_bins[i+1]);
    double pt_error = 0.5 * (pt_bins[i+1] - pt_bins[i]);
    g_eff->SetPoint(i, pt_center, efficiencies[i]);
    g_eff->SetPointError(i, pt_error, eff_errors[i]);
}

g_eff->SetTitle("Tag-and-Probe Efficiency vs p_{T}; Probe p_{T} [GeV]; Efficiency");
g_eff->SetMarkerStyle(20);
g_eff->SetMarkerColor(kBlue+1);
g_eff->SetLineColor(kBlue+1);

// 11. Draw and save
TCanvas* c_eff = new TCanvas("c_eff", "Efficiency vs pt", 800, 600);
g_eff->Draw("AP");
c_eff->SaveAs("efficiency_vs_pt.png");

}
