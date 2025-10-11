#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>

#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <algorithm>

static bool endsWith(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() &&
           std::equal(suf.rbegin(), suf.rend(), s.rbegin());
}

static void collectRootFiles(const std::string& folder,
                             std::vector<std::string>& out,
                             const std::string& patternSuffix = "_all.root")
{
    TSystemDirectory dir("dir", folder.c_str());
    TList* entries = dir.GetListOfFiles();
    if (!entries) return;

    entries->Sort();
    TIter it(entries);
    while (TObject* obj = it()) {
        auto f = dynamic_cast<TSystemFile*>(obj);
        if (!f) continue;
        std::string name = f->GetName();
        if (name == "." || name == "..") continue;

        std::string path = folder + "/" + name;
        if (f->IsDirectory()) {
            collectRootFiles(path, out, patternSuffix); // recurse
        } else {
            if (endsWith(name, patternSuffix)) out.push_back(path);
        }
    }
}

void CompareHists(const char* folder = "/eos/user/d/dblend/Analyzed_SkimmedFiles/12022/MC",
                  const char* subdir = "corr_DM_VBF",
                  const char* outReport = "compare_hists_report.txt")
{
    // 1) Recursively collect *_all.root files
    std::vector<std::string> filepaths;
    collectRootFiles(folder, filepaths, "_all.root");

    if (filepaths.empty()) {
        std::cerr << "No *_all.root files found under " << folder << "\n";
        return;
    }
    std::cout << "Found " << filepaths.size() << " files under " << folder << "\n";

    // 2) For each file, collect histogram names in the subdir
    std::map<std::string, std::set<std::string>> file2hists;
    std::set<std::string> allHists;

    for (const auto& path : filepaths) {
        TFile f(path.c_str(), "READ");
        if (f.IsZombie()) { std::cerr << "Could not open " << path << "\n"; continue; }

        TDirectory* d = (subdir && std::string(subdir).size())
                        ? dynamic_cast<TDirectory*>(f.Get(subdir))
                        : (TDirectory*)&f;
        if (!d) {
            std::cerr << "[MISS DIR] " << path << " has no subdir " << subdir << "\n";
            continue;
        }

        std::set<std::string> names;
        TIter nextkey(d->GetListOfKeys());
        while (TKey* key = (TKey*)nextkey()) {
            std::string cls = key->GetClassName();
            if (cls.rfind("TH", 0) == 0) {
                names.insert(key->GetName());
                allHists.insert(key->GetName());
            }
        }
        file2hists[path] = std::move(names);
    }

    // 3) Report
    std::ofstream rep(outReport);
    std::ostream& out = rep ? rep : std::cout;

    out << "Folder: " << folder << "\n";
    out << "Subdir: " << subdir << "\n";
    out << "Files scanned: " << filepaths.size() << "\n\n";

    out << "=== Per-file histogram counts in " << subdir << " ===\n";
    for (const auto& p : file2hists) {
        out << p.first << " : " << p.second.size() << " hists\n";
    }
    out << "\n";

    out << "=== Missing matrix (histogram -> files missing it) ===\n";
    size_t missingRows = 0;
    for (const auto& hname : allHists) {
        std::vector<std::string> missingIn;
        for (const auto& kv : file2hists) {
            if (!kv.second.count(hname)) missingIn.push_back(kv.first);
        }
        if (!missingIn.empty()) {
            out << hname << " :\n";
            for (const auto& fp : missingIn) out << "   - " << fp << "\n";
            ++missingRows;
        }
    }
    if (missingRows == 0) out << "(No missing histograms; all files consistent.)\n";

    if (rep) {
        rep.close();
        std::cout << "Wrote report: " << outReport << "\n";
    }
}

