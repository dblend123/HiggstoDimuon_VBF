// makeJSONfromList.C  — build a CMS-style (run,[ls_min,ls_max]) JSON
// from a list of skimmed NanoAOD ROOT files.
//
// Usage:
//   root -l -q 'makeJSONfromList.C("Muon_12022.list","seenLS_12022.json")'
// Then:
//   brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json \
//                 -u /fb -i seenLS_12022.json

#include "TChain.h"
#include "TBranch.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>

// collapse a sorted set of lumi numbers into [a,b] intervals
static std::vector<std::pair<unsigned,unsigned>>
toIntervals(const std::set<unsigned>& s) {
  std::vector<std::pair<unsigned,unsigned>> out;
  if (s.empty()) return out;
  auto it = s.begin();
  unsigned a = *it, b = *it;
  ++it;
  for (; it != s.end(); ++it) {
    if (*it == b + 1) { b = *it; continue; }
    out.emplace_back(a, b);
    a = b = *it;
  }
  out.emplace_back(a, b);
  return out;
}

// tiny trim
static inline void trim(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                  [](unsigned char c){ return !std::isspace(c); }));
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char c){ return !std::isspace(c); }).base(), s.end());
}

void makeJSONfromList(const char* listFile="files.list",
                      const char* outJSON ="seenLS.json",
                      const char* treeName="Events")
{
  // 1) Build a TChain from the file list
  TChain ch(treeName);
  std::ifstream fin(listFile);
  if (!fin) { std::cerr << "[ERROR] cannot open list: " << listFile << "\n"; return; }

  std::string line; int nAdded = 0;
  while (std::getline(fin, line)) {
    trim(line);
    if (line.empty() || line[0]=='#') continue;
    if (ch.Add(line.c_str()) == 0) {
      std::cerr << "[WARN] could not add: " << line << "\n";
    } else {
      ++nAdded;
    }
  }
  if (nAdded == 0) { std::cerr << "[ERROR] no files were added from " << listFile << "\n"; return; }

  // 2) Verify required branches exist
  if (!ch.GetBranch("run") || !ch.GetBranch("luminosityBlock")) {
    std::cerr << "[ERROR] Tree '" << treeName
              << "' must have branches 'run' and 'luminosityBlock'.\n";
    return;
  }

  // 3) Speed/IO setup (read only needed branches)
  ch.SetBranchStatus("*", 0);
  ch.SetBranchStatus("run", 1);
  ch.SetBranchStatus("luminosityBlock", 1);

  UInt_t run = 0, ls = 0;
  ch.SetBranchAddress("run", &run);
  ch.SetBranchAddress("luminosityBlock", &ls);

  // 4) Collect (run -> set of lumis)
  std::map<unsigned, std::set<unsigned>> seen;
  const Long64_t N = ch.GetEntries();
  std::cout << "[INFO] scanning " << N << " entries from " << nAdded << " files\n";

  for (Long64_t i = 0; i < N; ++i) {
    ch.GetEntry(i);
    seen[run].insert(ls);
    if ((i % 1000000) == 0) { std::cout << "  ... " << i << "/" << N << "\r" << std::flush; }
  }
  std::cout << "\n[INFO] collected " << seen.size() << " runs\n";

  // 5) Write CMS-style JSON
  std::ofstream out(outJSON);
  if (!out) { std::cerr << "[ERROR] cannot write " << outJSON << "\n"; return; }

  out << "{\n";
  bool firstRun = true;
  for (const auto& kv : seen) {
    if (!firstRun) out << ",\n";
    firstRun = false;
    out << "  \"" << kv.first << "\": [";
    auto ivals = toIntervals(kv.second);
    for (size_t i = 0; i < ivals.size(); ++i) {
      if (i) out << ", ";
      out << "[" << ivals[i].first << ", " << ivals[i].second << "]";
    }
    out << "]";
  }
  out << "\n}\n";
  out.close();

  std::cout << "[INFO] wrote " << outJSON << " (runs=" << seen.size() << ")\n";
}

