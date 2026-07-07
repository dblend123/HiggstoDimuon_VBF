// To run:
// 
// for era in 12022 22022 12023 22023 2024 12025
//do
//    echo
//    echo "=================================================="
//    echo "Running era: $era"
//    echo "=================================================="
//
//    root -l -b -q \
//    "PrintSystematics.C(\"/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/KIT/noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW\")"
//done
//
// root -l -b -q 'PrintSystematics.C("/eos/user/d/dblend/Analyzed_SkimmedFiles/12022/KIT/noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW")'

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>
#include <map>
#include <limits>

static TH1* getHist(TFile* f, const std::string& path) {
    if (!f) return nullptr;
    TH1* h = dynamic_cast<TH1*>(f->Get(path.c_str()));
    if (!h) {
        std::cerr << "[MISSING] " << path << "\n";
        return nullptr;
    }
    TH1* c = (TH1*)h->Clone();
    c->SetDirectory(nullptr);
    return c;
}

static void addHist(TH1* total, TFile* f, const std::string& path) {
    TH1* h = getHist(f, path);
    if (h) {
        total->Add(h);
        delete h;
    }
}

static double integralAll(TH1* h) {
    if (!h) return 0.0;
    return h->Integral(1, h->GetNbinsX());
}

static double integralWindow(TH1* h, double xmin, double xmax) {
    if (!h) return 0.0;
    int b1 = h->FindBin(xmin);
    int b2 = h->FindBin(xmax);
    return h->Integral(b1, b2);
}

static TH1* buildTotal(
    const std::string& name,
    const std::string& path,
    const std::vector<TFile*>& files,
    TH1* h_template
) {
    TH1* total = (TH1*)h_template->Clone(name.c_str());
    total->SetDirectory(nullptr);
    total->Reset("ICES");

    for (auto* f : files) addHist(total, f, path);
    return total;
}

static double envelopePct(TH1* h_nom, TH1* h_up, TH1* h_down, bool useWindow, double xmin, double xmax) {
    const double Nnom  = useWindow ? integralWindow(h_nom, xmin, xmax)  : integralAll(h_nom);
    const double Nup   = useWindow ? integralWindow(h_up, xmin, xmax)   : integralAll(h_up);
    const double Ndown = useWindow ? integralWindow(h_down, xmin, xmax) : integralAll(h_down);

    if (Nnom <= 0.0) return -1.0;

    const double upPct   = 100.0 * (Nup   - Nnom) / Nnom;
    const double downPct = 100.0 * (Ndown - Nnom) / Nnom;

    return std::max(std::abs(upPct), std::abs(downPct));
}

static void printSingle(
    const std::string& label,
    TH1* h_nom,
    TH1* h_up,
    TH1* h_down,
    bool useWindow,
    double xmin,
    double xmax
) {
    const double Nnom  = useWindow ? integralWindow(h_nom, xmin, xmax)  : integralAll(h_nom);
    const double Nup   = useWindow ? integralWindow(h_up, xmin, xmax)   : integralAll(h_up);
    const double Ndown = useWindow ? integralWindow(h_down, xmin, xmax) : integralAll(h_down);

    if (Nnom <= 0.0) {
        std::cout << std::left << std::setw(15) << label << "  Nnom <= 0\n";
        return;
    }

    const double upPct   = 100.0 * (Nup   - Nnom) / Nnom;
    const double downPct = 100.0 * (Ndown - Nnom) / Nnom;
    const double envPct  = std::max(std::abs(upPct), std::abs(downPct));

    std::cout << std::left << std::setw(15) << label
              << std::right
              << std::setw(12) << std::fixed << std::setprecision(4) << upPct
              << std::setw(12) << downPct
              << std::setw(12) << envPct
              << "\n";
}

void PrintSystematics(
    TString basePath,
    TString histName = "DiMuonMass_Higgs",
    TString nominalDir = "corr_DM_VBF_SR",
    bool useWindow = false,
    double xmin = 115.0,
    double xmax = 135.0
) {
    if (!basePath.EndsWith("/")) basePath += "/";

    std::unique_ptr<TFile> f_DY      (TFile::Open(basePath + "DY/DY_all.root"));
    std::unique_ptr<TFile> f_EWK     (TFile::Open(basePath + "EWK/EWK_all.root"));
    std::unique_ptr<TFile> f_TT      (TFile::Open(basePath + "TT/TT_all.root"));
    std::unique_ptr<TFile> f_WZ2L2Q  (TFile::Open(basePath + "WZto2L2Q/WZto2L2Q_all.root"));
    std::unique_ptr<TFile> f_WZ3LNu  (TFile::Open(basePath + "WZto3LNu/WZto3LNu_all.root"));
    std::unique_ptr<TFile> f_ZZ2L2Nu (TFile::Open(basePath + "ZZto2L2Nu/ZZto2L2Nu_all.root"));
    std::unique_ptr<TFile> f_ZZ2L2Q  (TFile::Open(basePath + "ZZto2L2Q/ZZto2L2Q_all.root"));
    std::unique_ptr<TFile> f_ZZ4L    (TFile::Open(basePath + "ZZto4L/ZZto4L_all.root"));
    std::unique_ptr<TFile> f_WWW     (TFile::Open(basePath + "WWW/WWW_all.root"));
    std::unique_ptr<TFile> f_WWZ     (TFile::Open(basePath + "WWZ/WWZ_all.root"));
    std::unique_ptr<TFile> f_WZZ     (TFile::Open(basePath + "WZZ/WZZ_all.root"));
    std::unique_ptr<TFile> f_ZZZ     (TFile::Open(basePath + "ZZZ/ZZZ_all.root"));

    std::vector<TFile*> files = {
        f_DY.get(), f_EWK.get(), f_TT.get(),
        f_WZ2L2Q.get(), f_WZ3LNu.get(),
        f_ZZ2L2Nu.get(), f_ZZ2L2Q.get(), f_ZZ4L.get(),
        f_WWW.get(), f_WWZ.get(), f_WZZ.get(), f_ZZZ.get()
    };

    struct SystDef {
        std::string label;
        std::string dir;
        std::string up;
        std::string down;
    };

    std::vector<SystDef> systs = {
        {"Scale",   "BDTTheorySyst_VBF",       "scaleUp",    "scaleDown"},
        {"PDF",     "BDTTheorySyst_VBF",       "pdfUp",      "pdfDown"},
        {"ISR",     "BDTTheorySyst_VBF",       "isrUp",      "isrDown"},
        {"FSR",     "BDTTheorySyst_VBF",       "fsrUp",      "fsrDown"},
        {"Pileup",  "BDTExperimentalSyst_VBF", "puUp",       "puDown"},
        {"MuonSF",  "BDTExperimentalSyst_VBF", "muonSFUp",   "muonSFDown"},
        {"Trigger", "BDTExperimentalSyst_VBF", "triggerUp",  "triggerDown"},
        {"BTag",    "BDTExperimentalSyst_VBF", "btagUp",     "btagDown"},
        {"ZPt",     "BDTExperimentalSyst_VBF", "zptUp",      "zptDown"}
    };

    std::vector<std::string> bdtVars = {
        "DiMuonMass_Higgs",
        "DiMuonpt",
        "DiMuonrapidity",
        "Muon1PTdDM",
        "Muon2PTdDM",
        "Muon1eta",
        "Muon2eta",
        "CollinsSoper",
        "cosCollinsSoper",
        "PtBalance",
        "Zeppenfeld"
    };

    // --------------------------
    // Single histogram printout
    // --------------------------
    const std::string nomPath = std::string(nominalDir.Data()) + "/" + std::string(histName.Data());
    TH1* h_template = getHist(f_DY.get(), nomPath);

    if (!h_template) {
        std::cerr << "[ERROR] Could not get template histogram: " << nomPath << "\n";
        return;
    }

    TH1* h_nom = buildTotal("h_nom_total", nomPath, files, h_template);

    std::cout << "\nSystematic yield variations for: " << histName << "\n";
    std::cout << "Nominal directory: " << nominalDir << "\n";
    if (useWindow) std::cout << "Integration window: [" << xmin << ", " << xmax << "]\n";
    else           std::cout << "Integration: full histogram\n";

    std::cout << "Nominal yield: "
              << (useWindow ? integralWindow(h_nom, xmin, xmax) : integralAll(h_nom))
              << "\n\n";

    std::cout << std::left << std::setw(15) << "Systematic"
              << std::right << std::setw(12) << "Up (%)"
              << std::setw(12) << "Down (%)"
              << std::setw(12) << "Envelope"
              << "\n";
    std::cout << "---------------------------------------------------\n";

    for (const auto& s : systs) {
        std::string upPath   = s.dir + "/" + std::string(histName.Data()) + "_" + s.up;
        std::string downPath = s.dir + "/" + std::string(histName.Data()) + "_" + s.down;

        TH1* h_up   = buildTotal("h_" + s.label + "_up",   upPath,   files, h_template);
        TH1* h_down = buildTotal("h_" + s.label + "_down", downPath, files, h_template);

        printSingle(s.label, h_nom, h_up, h_down, useWindow, xmin, xmax);

        delete h_up;
        delete h_down;
    }

    delete h_template;
    delete h_nom;

    // --------------------------
    // Range table across BDT vars
    // --------------------------
    std::map<std::string, double> minVal;
    std::map<std::string, double> maxVal;

    for (const auto& s : systs) {
        minVal[s.label] = std::numeric_limits<double>::infinity();
        maxVal[s.label] = -std::numeric_limits<double>::infinity();
    }

    for (const auto& var : bdtVars) {
        const std::string varNomPath = std::string(nominalDir.Data()) + "/" + var;

        TH1* h_temp = getHist(f_DY.get(), varNomPath);
        if (!h_temp) continue;

        TH1* h_var_nom = buildTotal("h_nom_" + var, varNomPath, files, h_temp);

        for (const auto& s : systs) {
            const std::string upPath   = s.dir + "/" + var + "_" + s.up;
            const std::string downPath = s.dir + "/" + var + "_" + s.down;

            TH1* h_up   = buildTotal("h_up_" + s.label + "_" + var, upPath, files, h_temp);
            TH1* h_down = buildTotal("h_down_" + s.label + "_" + var, downPath, files, h_temp);

            double env = envelopePct(h_var_nom, h_up, h_down, false, xmin, xmax);

            if (env >= 0.0 && std::isfinite(env)) {
                minVal[s.label] = std::min(minVal[s.label], env);
                maxVal[s.label] = std::max(maxVal[s.label], env);
            }

            delete h_up;
            delete h_down;
        }

        delete h_var_nom;
        delete h_temp;
    }

    std::cout << "\nRange across BDT variables\n";
    std::cout << "Nominal directory: " << nominalDir << "\n";
    std::cout << "Integration: full histogram for each BDT input\n\n";

    std::cout << std::left << std::setw(20) << "Uncertainty source"
              << std::right << std::setw(18) << "Range (%)"
              << "\n";
    std::cout << "--------------------------------------\n";

    for (const auto& s : systs) {
        if (!std::isfinite(minVal[s.label]) || !std::isfinite(maxVal[s.label])) {
            std::cout << std::left << std::setw(20) << s.label
                      << std::right << std::setw(18) << "n/a"
                      << "\n";
            continue;
        }

        std::ostringstream range;
        range << std::fixed << std::setprecision(2)
              << minVal[s.label] << " - " << maxVal[s.label];

        std::cout << std::left << std::setw(20) << s.label
                  << std::right << std::setw(18) << range.str()
                  << "\n";
    }
}
