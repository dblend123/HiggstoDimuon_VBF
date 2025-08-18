#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TKey.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>

void CombineEfficiencyWholeSpectrum(const char* listfile = "filelist.txt") {
    std::ifstream infile(listfile);
    if (!infile.is_open()) {
        std::cerr << "Cannot open input list file: " << listfile << std::endl;
        return;
    }

    std::vector<std::string> filenames;
    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty()) filenames.push_back(line);
    }
    infile.close();

    // Containers for graphs from all files
    std::vector<TGraphAsymmErrors*> eff_pt_graphs;
    std::vector<TGraphAsymmErrors*> eff_eta_graphs;
    std::vector<TGraphAsymmErrors*> eff_phi_graphs;

    // Read graphs from each file
    for (const auto& fname : filenames) {
        TFile* f = TFile::Open(fname.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Failed to open file: " << fname << std::endl;
            continue;
        }

        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom("TGraphAsymmErrors")) {
                TGraphAsymmErrors* g = (TGraphAsymmErrors*)obj;
                std::string title = g->GetTitle();

                if (title.find("Efficiency vs Probe p_{T}") != std::string::npos) {
                    eff_pt_graphs.push_back((TGraphAsymmErrors*)g->Clone());
                } else if (title.find("Efficiency vs Probe  #eta") != std::string::npos) {
                    eff_eta_graphs.push_back((TGraphAsymmErrors*)g->Clone());
                } else if (title.find("Efficiency vs Probe  #phi") != std::string::npos) {
                    eff_phi_graphs.push_back((TGraphAsymmErrors*)g->Clone());
                }
            }
        }
        f->Close();
    }

    if (eff_pt_graphs.empty() && eff_eta_graphs.empty() && eff_phi_graphs.empty()) {
        std::cerr << "No efficiency graphs found!" << std::endl;
        return;
    }

    // Function to average graphs (weighted average bin-by-bin)
    auto averageGraphs = [](const std::vector<TGraphAsymmErrors*>& graphs) -> TGraphAsymmErrors* {
        if (graphs.empty()) return nullptr;

        int npoints = graphs[0]->GetN();
        // Check all graphs have same number of points
        for (auto& g : graphs) {
            if (g->GetN() != npoints) {
                std::cerr << "Graphs have different number of points!" << std::endl;
                return nullptr;
            }
        }

        TGraphAsymmErrors* avg = (TGraphAsymmErrors*)graphs[0]->Clone("avgGraph");
        avg->SetTitle(graphs[0]->GetTitle());

        for (int i = 0; i < npoints; i++) {
            double x = 0, y_num = 0, y_den = 0;
            double ex_low = avg->GetErrorXlow(i);
            double ex_high = avg->GetErrorXhigh(i);

            // Weighted average: weights = 1/(error^2)
            // Use symmetric error as average of low/high error on y
            for (auto& g : graphs) {
                double yval, xval;
                g->GetPoint(i, xval, yval);
                double ey_low = g->GetErrorYlow(i);
                double ey_high = g->GetErrorYhigh(i);
                double eysym = 0.5*(ey_low + ey_high);

                if (eysym <= 0) {
                    // Avoid division by zero or nonsensical weight, skip this point
                    continue;
                }

                double w = 1.0 / (eysym * eysym);
                y_num += yval * w;
                y_den += w;
                x = xval;  // all xvals should be equal
            }

            if (y_den > 0) {
                double yavg = y_num / y_den;
                double eavg = 1.0 / sqrt(y_den);

                avg->SetPoint(i, x, yavg);
                avg->SetPointError(i, ex_low, ex_high, eavg, eavg);
            }
        }

        return avg;
    };

    // Average each set
    TGraphAsymmErrors* avg_pt = averageGraphs(eff_pt_graphs);
    TGraphAsymmErrors* avg_eta = averageGraphs(eff_eta_graphs);
    TGraphAsymmErrors* avg_phi = averageGraphs(eff_phi_graphs);

    // Set up TLatex for CMS-style labels
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextSize(0.04);

    // Lambda to draw and save a graph nicely
    auto drawAndSave = [&](TGraphAsymmErrors* g, const char* filename) {
        if (!g) return;

        TCanvas* c = new TCanvas("c", "Efficiency", 800, 600);
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlue);
        g->SetTitle(g->GetTitle());
        g->SetMinimum(0.75);
        g->SetMaximum(1.1);
        g->Draw("AP");

        latex.DrawLatex(0.18, 0.82, "#font[62]{CMS}");
        latex.DrawLatex(0.18, 0.78, "#it{Simulation}");
        latex.DrawLatex(0.18, 0.74, "#it{Preliminary}");
        latex.DrawLatex(0.65, 0.82, "Data Run 2022");
        latex.DrawLatex(0.55, 0.78, "20.0 GeV < muon p_{T} < 120.0 GeV");
        latex.DrawLatex(0.65, 0.74, "-2.4 < muon  #eta < 2.4");
        latex.DrawLatex(0.65, 0.70, "-#pi < muon  #phi < #pi");

        c->SaveAs(filename);
        delete c;
    };

    //drawAndSave(avg_pt, "Combined_efficiency_pt.pdf");
    //drawAndSave(avg_eta, "Combined_efficiency_eta.pdf");
    //drawAndSave(avg_phi, "Combined_efficiency_phi.pdf");
    drawAndSave(avg_pt, "Combined_efficiency_pt.png");
    drawAndSave(avg_eta, "Combined_efficiency_eta.png");
    drawAndSave(avg_phi, "Combined_efficiency_phi.png");

    std::cout << "Combined efficiency graphs saved as PDFs." << std::endl;
}

