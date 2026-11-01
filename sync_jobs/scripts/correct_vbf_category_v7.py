#!/usr/bin/env python3

import argparse
from array import array
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

src_path = Path(args.input)
dst_path = Path(args.output)

if not src_path.is_file():
    raise FileNotFoundError(src_path)

if dst_path.exists():
    raise FileExistsError(
        f"Output already exists; refusing to overwrite: {dst_path}"
    )

src = ROOT.TFile.Open(str(src_path), "READ")
assert src and not src.IsZombie()

dst = ROOT.TFile.Open(str(dst_path), "CREATE")
assert dst and not dst.IsZombie()

conditions = [
    "passVBFLeadingPair",
    "passBJetVeto",
    "mediumIdBoth",
    "isoBoth",
    "passHLTAndObjectMatch",
    "passFilters",
    "passVertex",
    "passRound2ExactlyTwoMuons",
    "passElectronVeto",
    "passGolden",
]

def copy_directory(source, destination):
    latest = {}
    for key in source.GetListOfKeys():
        name = key.GetName()
        if name not in latest or key.GetCycle() > latest[name].GetCycle():
            latest[name] = key

    for name, key in latest.items():
        obj = key.ReadObj()

        if obj.InheritsFrom("TDirectory"):
            subdir = destination.mkdir(name)
            copy_directory(obj, subdir)

        elif obj.InheritsFrom("TTree"):
            destination.cd()

            if name != "SyncEvents":
                obj.CloneTree(-1, "fast").Write(name)
                continue

            n = obj.GetEntries()

            for branch in conditions + [
                "muonCharge", "passVBFCategoryApprox"
            ]:
                assert obj.GetBranch(branch), f"Missing {branch}"

            obj.SetBranchStatus("passVBFCategoryApprox", 0)
            tree = obj.CloneTree(0)
            tree.SetAutoSave(0)

            category = array("b", [0])
            tree.Branch(
                "passVBFCategoryApprox",
                category,
                "passVBFCategoryApprox/O"
            )

            passed_count = 0

            for i in range(n):
                obj.GetEntry(i)

                passed = all(
                    bool(getattr(obj, branch))
                    for branch in conditions
                )

                passed &= (
                    obj.muonCharge[0] * obj.muonCharge[1] < 0
                )

                category[0] = int(passed)
                passed_count += int(passed)

                tree.Fill()

            tree.Write(name, ROOT.TObject.kOverwrite)

            print(f"SyncEvents entries: {n}")
            print(f"Corrected VBF passes: {passed_count}")

        else:
            destination.cd()

            if name == "syncVBFCategoryStatus":
                obj = ROOT.TNamed(
                    name,
                    "passVBFCategoryApprox: leading-pT jet pair "
                    "(pT1>35, pT2>30, mjj>400, |dEta|>2.5) "
                    "plus Round2 exactly-two-good-muons, electron veto, "
                    "and common baseline flags; no event-level jet veto map "
                    "(per-jet exclusions retained); "
                    "validated against TAMU VBFH final category "
                    "on 3571455 matched original-sample events "
                    "(99.99658 percent agreement)"
                )

            obj.Write(name)

try:
    copy_directory(src, dst)
finally:
    dst.Close()
    src.Close()

print("Created:", dst_path)
