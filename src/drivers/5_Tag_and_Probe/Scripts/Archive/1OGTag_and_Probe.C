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
#include <ROOT/RDFHelpers.hxx>

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
// Passing pt
df_with_probes = df_with_probes.Define("probe_pt_pass", R"cpp(
    std::vector<double> pts;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && Muon_tightId[j]) {
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
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && Muon_tightId[j]) {
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
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && Muon_tightId[j]) {
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
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && Muon_tightId[j]) {
                masses.push_back(m);
            }
        }
    }
    return masses;
)cpp");

// Step 4: Build invariant mass of failing tag-probe pairs
// Failing pt
df_with_probes = df_with_probes.Define("probe_pt_fail", R"cpp(
    std::vector<double> pts;
    for (size_t i = 0; i < Muon_pt.size(); ++i) {
        if (!tag_mask[i]) continue;
        for (size_t j = 0; j < Muon_pt.size(); ++j) {
            if (tag_mask[j] || i == j || Muon_charge[i] * Muon_charge[j] >= 0) continue;
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && !Muon_tightId[j]) {
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
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && !Muon_tightId[j]) {
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
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && !Muon_tightId[j]) {
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
            TLorentzVector tag, probe;
            tag.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
            probe.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], 0.105);
            double m = (tag + probe).M();
            if (m > 60 && m < 120 && !Muon_tightId[j]) {
                masses.push_back(m);
            }
        }
    }
    return masses;
)cpp");

// 1. Define pt bins
std::vector<float> pt_bins = {0, 10, 20, 30, 40, 50, 60, 80, 100};
const int n_pt_bins = pt_bins.size() - 1;

// 2. Define per-bin flags for passing and failing probes
for (int i = 0; i < n_pt_bins; ++i) {
    double pt_low = pt_bins[i];
    double pt_high = pt_bins[i + 1];

    std::string bin_name_pass = Form("pass_bin_pt_%d", i);
    std::string bin_name_fail = Form("fail_bin_pt_%d", i);

    df_with_probes = df_with_probes.Define(bin_name_pass, [pt_low, pt_high](const std::vector<double>& pts) {
        for (const auto& pt : pts)
            if (pt >= pt_low && pt < pt_high)
                return true;
        return false;
    }, {"probe_pt_pass"});

    df_with_probes = df_with_probes.Define(bin_name_fail, [pt_low, pt_high](const std::vector<double>& pts) {
        for (const auto& pt : pts)
            if (pt >= pt_low && pt < pt_high)
                return true;
        return false;
    }, {"probe_pt_fail"});
}

// 3. Create vectors to hold histograms
std::vector<TH1D*> hists_pass, hists_fail;

// 4. Book histograms per pt bin
for (int i = 0; i < n_pt_bins; ++i) {
    std::string bin_filter_pass = Form("pass_bin_pt_%d", i);
    std::string bin_filter_fail = Form("fail_bin_pt_%d", i);

    auto df_pass_bin = df_with_probes.Filter(bin_filter_pass + " == true");
    auto df_fail_bin = df_with_probes.Filter(bin_filter_fail + " == true");

    ROOT::RDF::TH1DModel h_pass_model(("h_mass_pass_bin_" + std::to_string(i)).c_str(), "Passing; m_{#mu#mu} [GeV]", 60, 60, 120);
    ROOT::RDF::TH1DModel h_fail_model(("h_mass_fail_bin_" + std::to_string(i)).c_str(), "Failing; m_{#mu#mu} [GeV]", 60, 60, 120);

    auto hist_pass = df_pass_bin.Histo1D(h_pass_model, "probe_mass_pass");
    auto h_pass_ptr = hist_pass.GetPtr();
    
    if (h_pass_ptr && h_pass_ptr->GetEntries() > 0) {
        hists_pass.push_back(h_pass_ptr);
    } else {
        hists_pass.push_back(nullptr);  // placeholder for skipped bin
    }
    
    auto hist_fail = df_fail_bin.Histo1D(h_fail_model, "probe_mass_fail");
    auto h_fail_ptr = hist_fail.GetPtr();
    
    if (h_fail_ptr && h_fail_ptr->GetEntries() > 0) {
        hists_fail.push_back(h_fail_ptr);
    } else {
        hists_fail.push_back(nullptr);  // placeholder for skipped bin
    }

    std::cout << "Pass bin " << i << " entries: " << df_pass_bin.Count().GetValue() << std::endl;
    std::cout << "Pass bin " << i << " hist entries: " << hist_pass->GetEntries() << std::endl;
    std::cout << "Fail bin " << i << " entries: " << df_fail_bin.Count().GetValue() << std::endl;
    std::cout << "Fail bin " << i << " hist entries: " << hist_fail->GetEntries() << std::endl;
}

// 5. Compute efficiencies
std::vector<double> efficiencies, eff_errors;

if (hists_pass.size() < n_pt_bins || hists_fail.size() < n_pt_bins) {
    std::cerr << "Mismatch between number of pt bins and histograms!" << std::endl;
    std::cerr << "n_pt_bins = " << n_pt_bins 
              << ", hists_pass.size() = " << hists_pass.size()
              << ", hists_fail.size() = " << hists_fail.size() << std::endl;
    return;
}

for (int i = 0; i < n_pt_bins; ++i) {
    auto h_pass = hists_pass[i];
    auto h_fail = hists_fail[i];

     if (!h_pass || h_pass->GetEntries() == 0) {
        std::cerr << "Skipping fit for empty or invalid pass bin " << i << std::endl;
        continue;
    }

    if (!h_fail || h_fail->GetEntries() == 0) {
        std::cerr << "Skipping fit for empty or invalid fail bin " << i << std::endl;
        continue;
    }

    // Fit signal+background (Gaussian + exponential)
    TF1* fit_pass = new TF1(Form("fit_pass_bin%d", i), "gaus(0) + expo(3)", 60, 120);
    TF1* fit_fail = new TF1(Form("fit_fail_bin%d", i), "gaus(0) + expo(3)", 60, 120);

    h_pass->Fit(fit_pass, "RQ0");
    h_fail->Fit(fit_fail, "RQ0");

    double sig_yield_pass = fit_pass->GetParameter(0) * fit_pass->GetParameter(2) * std::sqrt(2 * M_PI);
    double sig_yield_fail = fit_fail->GetParameter(0) * fit_fail->GetParameter(2) * std::sqrt(2 * M_PI);

    double denom = sig_yield_pass + sig_yield_fail;
    double eff = (denom > 0) ? sig_yield_pass / denom : 0.0;
    double err = (denom > 0) ? std::sqrt(eff * (1 - eff) / denom) : 0.0;

    efficiencies.push_back(eff);
    eff_errors.push_back(err);
}

// 6. Draw efficiency vs pt
auto g_eff = new TGraphErrors(n_pt_bins);
for (int i = 0; i < n_pt_bins; ++i) {
    double pt_center = 0.5 * (pt_bins[i] + pt_bins[i + 1]);
    double pt_error = 0.5 * (pt_bins[i + 1] - pt_bins[i]);
    g_eff->SetPoint(i, pt_center, efficiencies[i]);
    g_eff->SetPointError(i, pt_error, eff_errors[i]);
}

g_eff->SetTitle("Tag-and-Probe Efficiency vs p_{T}; Probe p_{T} [GeV]; Efficiency");
g_eff->SetMarkerStyle(20);
g_eff->SetMarkerColor(kBlue + 1);
g_eff->SetLineColor(kBlue + 1);

TCanvas* c_eff = new TCanvas("c_eff", "Efficiency vs pt", 800, 600);
g_eff->Draw("AP");
c_eff->SaveAs("efficiency_vs_pt.png");

}
