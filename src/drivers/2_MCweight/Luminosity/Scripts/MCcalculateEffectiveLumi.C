#include <ROOT/RDataFrame.hxx>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <chrono>

void MCcalculateEffectiveLumi(
    const char* inputList = "test.list", double crossSection_pb = 1.0) {
    ROOT::EnableImplicitMT(); // optional multithreading

    // Start total timer
    auto start_total = std::chrono::high_resolution_clock::now();

    // Read the file list
    std::vector<std::string> files;
    std::ifstream infile(inputList);
    std::string line;
    
    int fileCounter = 0;
    
    while (std::getline(infile, line)) {
        if (!line.empty()) { 
           // Try to open the file quickly just to check
            TFile *f = TFile::Open(line.c_str());
            if (f && !f->IsZombie()) {
                files.push_back(line);
                std::cout << "File " << ++fileCounter << ": " << line << " [OK]" << std::endl;
            } else {
                std::cerr << "File " << line << " [FAILED TO OPEN]" << std::endl;
            }
            if (f) f->Close();
        }
    }

    if (files.empty()) {
        std::cerr << "No valid ROOT files found in " << inputList << std::endl;
        return;
    }

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Total valid files to process: " << files.size() << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    // Start processing timer
    auto start_processing = std::chrono::high_resolution_clock::now();// Start processing timer

    std::vector<double> cumulativeLumi_fb;
    std::vector<double> nFiles;
    double totalSumWeight = 0.0;

    for (size_t i = 0; i < files.size(); ++i) {
        ROOT::RDataFrame df("Events", files[i]);
        auto colNames = df.GetColumnNames();
        bool hasGenWeight = std::find(colNames.begin(), colNames.end(), "genWeight") != colNames.end();

        double fileWeight = 0.0;
        if (hasGenWeight) {
            fileWeight = *df.Sum("genWeight");
        } else {
            fileWeight = df.Count().GetValue(); // assume weight = 1
        }

        totalSumWeight += fileWeight;
        double lumi_fb = (totalSumWeight / crossSection_pb) / 1000.0; // to /fb

        cumulativeLumi_fb.push_back(lumi_fb);
        nFiles.push_back(i + 1);
    }

     // End processing timer
    auto end_processing = std::chrono::high_resolution_clock::now();

    // Compute effective luminosity
    double effectiveLumi_pb = totalSumWeight / crossSection_pb;
    double effectiveLumi_fb = effectiveLumi_pb / 1000.0;

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Input list: " << inputList << std::endl;
    std::cout << "Number of ROOT files = " << files.size() << std::endl;
    std::cout << "Cross section = " << crossSection_pb << " pb" << std::endl;
    std::cout << "Sum of event weights = " << totalSumWeight << std::endl;
    std::cout << "Effective Luminosity = " << effectiveLumi_pb << " /pb" << std::endl;
    std::cout << "Effective Luminosity = " << effectiveLumi_fb << " /fb" << std::endl;
    std::cout << "-------------------------------------" << std::endl;

    // Stop total timer
    auto end_total = std::chrono::high_resolution_clock::now();

    // Timing output
    std::chrono::duration<double> duration_total = end_total - start_total;
    std::chrono::duration<double> duration_processing = end_processing - start_processing;
    std::cout << "Time spent processing events: " << duration_processing.count() << " seconds" << std::endl;
    std::cout << "Total execution time: " << duration_total.count() << " seconds" << std::endl;

    // Create plot
    TCanvas* c = new TCanvas("c", "Cumulative Effective Lumi vs Files", 800, 600);
    c->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.04);

    TGraph* gr = new TGraph(nFiles.size(), &nFiles[0], &cumulativeLumi_fb[0]);
    gr->SetTitle("Cumulative Effective Luminosity;Number of Files;Effective Lumi (fb^{-1})");
    gr->GetYaxis()->SetTitleOffset(1.4);
    gr->SetLineColor(kAzure+7);
    gr->SetLineWidth(2);
    gr->Draw("AL");

    // CMS-style label
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextFont(62); // Bold for "CMS"
    latex.DrawLatex(0.15, 0.85, "CMS");

    latex.SetTextFont(52); // Plain for "Simulation Preliminary"
    latex.DrawLatex(0.24, 0.85, "Simulation Preliminary");

    latex.SetTextFont(42); // Plain text for lumi
    TString lumiLabel = Form("MC @ %.3f pb", crossSection_pb);
    latex.DrawLatex(0.65, 0.85, lumiLabel);

    // Save
    c->SaveAs("cumulative_lumi_vs_nfiles.pdf");
}

