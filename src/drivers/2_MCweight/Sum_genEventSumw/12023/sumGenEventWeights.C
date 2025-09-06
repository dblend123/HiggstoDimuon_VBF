void sumGenEventWeights() {
    // Define all filelists (label, path)
    std::vector<std::pair<std::string, std::string>> filelists = {
        // Backgrounds
        {"DY",        "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/DYto2L-2Jets/fullpath.list"},
        {"EWK",       "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/EWK/fullpath.list"},
        {"TTto2L2Nu", "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/TTto2L2Nu/fullpath.list"},
        {"WWW",       "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/WWW/fullpath.list"},
        {"WWZ",       "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/WWZ/fullpath.list"},
        {"WZZ",       "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/WZZ/fullpath.list"},
        {"WZto2L2Q",  "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/WZto2L2Q/fullpath.list"},
        {"WZto3LNu",  "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/WZto3LNu/fullpath.list"},
        {"ZZZ",       "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/ZZZ/fullpath.list"},
        {"ZZto2L2Nu", "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/ZZto2L2Nu/fullpath.list"},
        {"ZZto2L2Q",  "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/ZZto2L2Q/fullpath.list"},
        {"ZZto4L",    "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Background_listfiles/ZZto4L/fullpath.list"},

        // Signals
        {"VBFHTo2Mu", "/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/12023/Signal_listfiles/VBFHto2Mu/fullpath.list"},
    };

    // Output file to store the results
    //std::ofstream outfile("justtheonefile_genEventSumw_output.txt");
    std::ofstream outfile("genEventSumw_output.txt");
    //std::ofstream outfile("Ext_genEventSumw_output.txt");
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return;
    }

    for (const auto& [label, filelist] : filelists) {
        std::ifstream infile(filelist);
        if (!infile.is_open()) {
            std::cerr << "Could not open list: " << filelist << std::endl;
            outfile << label << "\tERROR: Could not open list file.\n";
            continue;
        }

        std::string filename;
        double totalWeight = 0.0;
        int nFiles = 0;

        while (std::getline(infile, filename)) {
            if (filename.empty()) continue;

            TFile* f = TFile::Open(filename.c_str());
            if (!f || f->IsZombie()) {
                std::cerr << "Could not open file: " << filename << std::endl;
                continue;
            }

            TTree* runs = (TTree*)f->Get("Runs");
            if (!runs) {
                std::cerr << "Runs tree not found in: " << filename << std::endl;
                f->Close();
                continue;
            }

            Double_t genEventSumw = 0;
            runs->SetBranchAddress("genEventSumw", &genEventSumw);

            for (Long64_t i = 0; i < runs->GetEntries(); ++i) {
                runs->GetEntry(i);
                totalWeight += genEventSumw;
            }

            nFiles++;
            f->Close();
        }

        std::cout << "==== " << label << " ====" << std::endl;
        std::cout << "Processed " << nFiles << " files." << std::endl;
        std::cout << "Total genEventSumw: " << totalWeight << "\n" << std::endl;

        outfile << label << "\t" << nFiles << "\t" << totalWeight << "\n";
    }

    outfile.close();
    std::cout << "All results written to genEventSumw_output.txt" << std::endl;
    //std::cout << "All results written to Ext_genEventSumw_output.txt" << std::endl;
}

