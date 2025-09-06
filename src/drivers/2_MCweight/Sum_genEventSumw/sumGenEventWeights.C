void sumGenEventWeights() {
    // Define all filelists (label, path)
    std::vector<std::pair<std::string, std::string>> filelists = {
        //Data
        //{"Muon", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/DYBackground_listfiles/v5-v2/v5-v2_fullpaths.list"},

        // Backgrounds
        //{"DY", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/DYBackground_listfiles/v5-v2/v5-v2_fullpaths.list"},
        //{"DY_Ext", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/DYBackground_listfiles/v5_ext1-v1/v5_ext1-v1_fullpaths.list"},
        //{"EWK", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/EWKBackground_listfiles/EWK_2L2J_All.list"},
        //{"TTTo2L2Nu", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/TTBackground_listfiles/TTTo2L2Nu_All.list"},
        //{"TTTo2L2Nu_Ext", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/TTBackground_listfiles/TTTo2L2Nu_Ext_All.list"},
        //{"WZTo2L2Q", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/WZBackground_listfiles/WZto2L2Q_All.txt"},
        //{"WZTo2L2Q_Ext", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/WZBackground_listfiles/WZto2L2Q_Ext_All.txt"},
        //{"WZto3LNu", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/WZto3LNuBackground_listfiles/WZto3LNu_full.list"},
        //{"WZto3LNu_Ext", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/WZto3LNuBackground_listfiles/WZto3LNu_ext1_v3_Ext_full.list"},
        //{"ZZTo2L2Nu", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/ZZto2L2NuBackground_listfiles/ZZTo2L2Nu_All.list"},
        //{"ZZTo2L2Nu_Ext", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/ZZto2L2NuBackground_listfiles/ZZTo2L2Nu_Ext_All.list"},
        //{"ZZTo2L2Q", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/ZZto2L2QBackground_listfiles/ZZTo2L2Q_All.list"},
        //{"ZZTo2L2Q_Ext", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Background_listfiles/ZZto2L2QBackground_listfiles/ZZTo2L2Q_Ext_All.list"},

        // Signals
        //{"VBFHTo2Mu", "/afs/cern.ch/user/d/dblend/test-analysis/build/src/drivers/1_ListFiles/2022/Signal_listfiles/VBFSignal_listfiles/VBFSignal_FullRun22.lsit"},
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

