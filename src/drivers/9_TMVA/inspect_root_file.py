#!/usr/bin/env python3
"""Quick script to inspect ROOT file structure"""

import argparse
import ROOT

def main():
    p = argparse.ArgumentParser(description="Inspect ROOT file structure")
    p.add_argument("file", help="Path to ROOT file")
    args = p.parse_args()

    f = ROOT.TFile.Open(args.file)
    if not f or f.IsZombie():
        print(f"[ERROR] Could not open {args.file}")
        return

    print(f"\n=== File: {args.file} ===\n")
    
    # List all keys (trees, histograms, etc.)
    print("Available objects in file:")
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        print(f"  - {key.GetName():30s} (Type: {key.GetClassName()})")
        
        # If it's a tree, show branches
        if isinstance(obj, ROOT.TTree):
            tree = obj
            print(f"    Entries: {tree.GetEntries()}")
            print(f"    Branches:")
            for branch in tree.GetListOfBranches():
                print(f"      • {branch.GetName()}")

    f.Close()

if __name__ == "__main__":
    main()
