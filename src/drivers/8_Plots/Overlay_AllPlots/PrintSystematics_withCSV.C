// To run one era and write CSV:
// root -l -b -q 'PrintSystematics_withCSV.C("/eos/user/d/dblend/Analyzed_SkimmedFiles/12022/KIT/noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW", "DiMuonMass_Higgs", "corr_DM_VBF_SR", false, 115.0, 135.0, "systematics_summary.csv", "12022")'
//
// To run all eras:
// rm -f systematics_summary.csv
// for era in 12022 22022 12023 22023 2024 2025; do
//   root -l -b -q "PrintSystematics_withCSV.C(\"/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/KIT/noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW\", \"DiMuonMass_Higgs\", \"corr_DM_VBF_SR\", false, 115.0, 135.0, \"systematics_summary.csv\", \"${era}\")"
// done

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
#include <sstream>
#include <fstream>
#include <string>

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

struct VariationResult {
    double nominal = 0.0;
    double upPct = 0.0;
    double downPct = 0.0;
    double envelopePct = -1.0;
    bool valid = false;
};

static VariationResult calcVariation(
    TH1* h_nom,
    TH1* h_up,
    TH1* h_down,
    bool useWindow,
    double xmin,
    double xmax
) {
    VariationResult r;
    r.nominal = useWindow ? integralWindow(h_nom, xmin, xmax) : integralAll(h_nom);
    const double Nup   = useWindow ? integralWindow(h_up, xmin, xmax) : integralAll(h_up);
    const double Ndown = useWindow ? integralWindow(h_down, xmin, xmax) : integralAll(h_down);

    if (r.nominal <= 0.0) return r;

    r.upPct       = 100.0 * (Nup   - r.nominal) / r.nominal;
    r.downPct     = 100.0 * (Ndown - r.nominal) / r.nominal;
    r.envelopePct = std::max(std::abs(r.upPct), std::abs(r.downPct));
    r.valid       = true;
    return r;
}

static std::string csvEscape(const std::string& in) {
    bool needsQuotes = false;
    for (char c : in) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            needsQuotes = true;
            break;
        }
    }
    if (!needsQuotes) return in;

    std::string out = "\"";
    for (char c : in) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += "\"";
    return out;
}

static bool fileIsEmptyOrMissing(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in.good()) return true;
    return in.peek() == std::ifstream::traits_type::eof();
}

static void writeCsvHeaderIfNeeded(std::ofstream& csv, const std::string& path) {
    if (!fileIsEmptyOrMissing(path)) return;
    csv << "row_type,era,hist_name,nominal_dir,integration,xmin,xmax,nominal_yield,"
        << "source,up_pct,down_pct,envelope_pct,min_range_pct,max_range_pct,range_text\n";
}

static void writeSingleCsvRow(
    std::ofstream& csv,
    const std::string& era,
    const std::string& histName,
    const std::string& nominalDir,
    const std::string& integration,
    double xmin,
    double xmax,
    const std::string& source,
    const VariationResult& r
) {
    if (!csv.good() || !r.valid) return;

    csv << "single,"
        << csvEscape(era) << ","
        << csvEscape(histName) << ","
        << csvEscape(nominalDir) << ","
        << csvEscape(integration) << ","
        << xmin << ","
        << xmax << ","
        << std::setprecision(10) << r.nominal << ","
        << csvEscape(source) << ","
        << std::fixed << std::setprecision(4) << r.upPct << ","
        << r.downPct << ","
        << r.envelopePct << ",,,\n";
}

static void writeRangeCsvRow(
    std::ofstream& csv,
    const std::string& era,
    const std::string& histName,
    const std::string& nominalDir,
    const std::string& source,
    double minVal,
    double maxVal
) {
    if (!csv.good() || !std::isfinite(minVal) || !std::isfinite(maxVal)) return;

    std::ostringstream range;
    range << std::fixed << std::setprecision(2) << minVal << " - " << maxVal;

    csv << "range,"
        << csvEscape(era) << ","
        << csvEscape(histName) << ","
        << csvEscape(nominalDir) << ",full_bdt_inputs,,,,"
        << csvEscape(source) << ",,,,"
        << std::fixed << std::setprecision(4) << minVal << ","
        << maxVal << ","
        << csvEscape(range.str()) << "\n";
}

void PrintSystematics(
    TString basePath,
    TString histName = "DiMuonMass_Higgs",
    TString nominalDir = "corr_DM_VBF_SR",
    bool useWindow = false,
    double xmin = 115.0,
    double xmax = 135.0,
    TString outCsv = "",
    TString eraLabel = ""
) {
    if (!basePath.EndsWith("/")) basePath += "/";

    const std::string era = eraLabel.Length() ? std::string(eraLabel.Data()) : "unknown";
    const std::string integration = useWindow ? "window" : "full_histogram";

    std::ofstream csv;
    std::string outCsvPath = std::string(outCsv.Data());
    if (outCsv.Length()) {
        const bool needHeader = fileIsEmptyOrMissing(outCsvPath);
        csv.open(outCsvPath.c_str(), std::ios::out | std::ios::app);
        if (!csv.good()) {
            std::cerr << "[WARNING] Could not open CSV for writing: " << outCsvPath << "\n";
        } else if (needHeader) {
            csv << "row_type,era,hist_name,nominal_dir,integration,xmin,xmax,nominal_yield,"
                << "source,up_pct,down_pct,envelope_pct,min_range_pct,max_range_pct,range_text\n";
        }
    }

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
        {"QCD",     "BDTTheorySyst_VBF",       "scaleUp",    "scaleDown"},
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

    const std::string nomPath = std::string(nominalDir.Data()) + "/" + std::string(histName.Data());
    TH1* h_template = getHist(f_DY.get(), nomPath);

    if (!h_template) {
        std::cerr << "[ERROR] Could not get template histogram: " << nomPath << "\n";
        return;
    }

    TH1* h_nom = buildTotal("h_nom_total", nomPath, files, h_template);
    const double nominalYield = useWindow ? integralWindow(h_nom, xmin, xmax) : integralAll(h_nom);

    std::cout << "\nSystematic yield variations for: " << histName << "\n";
    std::cout << "Nominal directory: " << nominalDir << "\n";
    if (useWindow) std::cout << "Integration window: [" << xmin << ", " << xmax << "]\n";
    else           std::cout << "Integration: full histogram\n";
    std::cout << "Nominal yield: " << nominalYield << "\n\n";

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

        VariationResult r = calcVariation(h_nom, h_up, h_down, useWindow, xmin, xmax);

        if (!r.valid) {
            std::cout << std::left << std::setw(15) << s.label << "  Nnom <= 0\n";
        } else {
            std::cout << std::left << std::setw(15) << s.label
                      << std::right
                      << std::setw(12) << std::fixed << std::setprecision(4) << r.upPct
                      << std::setw(12) << r.downPct
                      << std::setw(12) << r.envelopePct
                      << "\n";
            writeSingleCsvRow(csv, era, histName.Data(), nominalDir.Data(), integration, xmin, xmax, s.label, r);
        }

        delete h_up;
        delete h_down;
    }

    delete h_template;
    delete h_nom;

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

            VariationResult r = calcVariation(h_var_nom, h_up, h_down, false, xmin, xmax);
            if (r.valid && std::isfinite(r.envelopePct)) {
                minVal[s.label] = std::min(minVal[s.label], r.envelopePct);
                maxVal[s.label] = std::max(maxVal[s.label], r.envelopePct);
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

        writeRangeCsvRow(csv, era, histName.Data(), nominalDir.Data(), s.label, minVal[s.label], maxVal[s.label]);
    }

    if (csv.good()) {
        csv.close();
        std::cout << "\n[OK] Wrote CSV rows to: " << outCsvPath << "\n";
    }
}

void PrintSystematics_withCSV(
    TString basePath,
    TString histName = "DiMuonMass_Higgs",
    TString nominalDir = "corr_DM_VBF_SR",
    bool useWindow = false,
    double xmin = 115.0,
    double xmax = 135.0,
    TString outCsv = "",
    TString eraLabel = ""
) {
    PrintSystematics(basePath, histName, nominalDir, useWindow, xmin, xmax, outCsv, eraLabel);
}
