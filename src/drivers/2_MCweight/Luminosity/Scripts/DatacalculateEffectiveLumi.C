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

void DatacalculateEffectiveLumi(const char* inputList = "1test.list") {
    ROOT::EnableImplicitMT();

    auto start_total = std::chrono::high_resolution_clock::now();

    std::vector<std::string> files;
    std::ifstream infile(inputList);
    std::string line;

    std::vector<double> nFiles;
    std::vector<double> cumulativeLumi_fb;

    double totalLumi_pb = 0.0;
    int fileCounter = 0;

    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        TFile* f = TFile::Open(line.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Failed to open: " << line << std::endl;
            continue;
        }

        std::cout << "File " << ++fileCounter << ": " << line << " [OK]" << std::endl;

        TTree* lumiTree = (TTree*) f->Get("LuminosityBlocks");
        if (!lumiTree) {
            std::cerr << "No LuminosityBlocks tree in: " << line << std::endl;
            f->Close();
            continue;
        }

        Float_t intLumi = 0;
        if (lumiTree->GetBranch("integratedLumi")) {
            lumiTree->SetBranchAddress("integratedLumi", &intLumi);
        } else {
            std::cerr << "No 'integratedLumi' branch in: " << line << std::endl;
            f->Close();
            continue;
        }

        for (Long64_t i = 0; i < lumiTree->GetEntries(); ++i) {
            lumiTree->GetEntry(i);
            totalLumi_pb += intLumi;
        }

        f->Close();

        // Save cumulative value for plotting
        nFiles.push_back(fileCounter);
        cumulativeLumi_fb.push_back(totalLumi_pb / 1000.0);
    }

    if (nFiles.empty()) {
        std::cerr << "No valid files processed." << std::endl;
        return;
    }

    auto end_processing = std::chrono::high_resolution_clock::now();
    auto end_total = std::chrono::high_resolution_clock::now();

    double totalLumi_fb = totalLumi_pb / 1000.0;

    std::cout << "-------------------------------------------" << std::endl;
    std::cout << "Input list: " << inputList << std::endl;
    std::cout << "Number of real data files processed: " << fileCounter << std::endl;
    std::cout << "Total integrated luminosity: " << totalLumi_pb << " /pb" << std::endl;
    std::cout << "                              = " << totalLumi_fb << " /fb" << std::endl;
    std::cout << "-------------------------------------------" << std::endl;

    std::chrono::duration<double> duration_processing = end_processing - end_total;
    std::chrono::duration<double> duration_total = end_total - start_total;

    std::cout << "Time spent processing: " << duration_processing.count() << " s" << std::endl;
    std::cout << "Total time: " << duration_total.count() << " s" << std::endl;

    // Plot
    TCanvas* c = new TCanvas("c", "Total Integrated Luminosity vs Files", 800, 600);
    c->SetLeftMargin(0.15);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.04);

    TGraph* gr = new TGraph(nFiles.size(), &nFiles[0], &cumulativeLumi_fb[0]);
    gr->SetTitle("Cumulative Integrated Luminosity (Real Data);Number of Files;Lumi (fb^{-1})");
    gr->GetYaxis()->SetTitleOffset(1.4);
    gr->SetLineColor(kBlue+2);
    gr->SetLineWidth(2);
    gr->Draw("AL");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextFont(62);
    latex.DrawLatex(0.15, 0.85, "CMS");

    latex.SetTextFont(52);
    latex.DrawLatex(0.24, 0.85, "Data Preliminary");

    c->SaveAs("realdata_cumulative_lumi_vs_nfiles.pdf");
}

