#!/usr/bin/env python3
"""
Diagnostic script to investigate why EWK BDT has poor performance on 2024.

Usage:
python3 diagnose_ewk_problem.py --era 2024 --dir DataPreparation_Files/2024
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from tmva100_DataPreparation import VARIABLES
from tmva101_Training import load_folds, TREENAME

def parse_args():
    p = argparse.ArgumentParser(description="Diagnose EWK model performance")
    p.add_argument("--era", required=True, help="Era (e.g., 2024)")
    p.add_argument("--dir", default="DataPreparation_Files", help="Directory with prepared files")
    p.add_argument("--model", help="Optional: path to trained EWK model for feature importance")
    return p.parse_args()

def load_category_data(dir_path, era, category):
    """Load training+test data for a specific background category"""
    sig_files = [os.path.join(dir_path, f"{era}_all_signal.root")]
    bkg_files = [os.path.join(dir_path, f"{era}_all_{category}.root")]
    
    for f in sig_files + bkg_files:
        if not os.path.exists(f):
            print(f"[WARN] Missing: {f}")
            return None
    
    # Load all data (not just test fold)
    try:
        rdf_sig = ROOT.RDataFrame(TREENAME, ROOT.std.vector("string")(sig_files))
        rdf_bkg = ROOT.RDataFrame(TREENAME, ROOT.std.vector("string")(bkg_files))
        
        # Get basic stats
        sig_data = rdf_sig.AsNumpy(VARIABLES + ["weight"] if "weight" in dir(rdf_sig) else VARIABLES)
        bkg_data = rdf_bkg.AsNumpy(VARIABLES + ["weight"] if "weight" in dir(rdf_bkg) else VARIABLES)
        
        return sig_data, bkg_data, rdf_sig.GetEntries(), rdf_bkg.GetEntries()
    except Exception as e:
        print(f"[ERROR] Failed to load {category}: {e}")
        return None

def compare_distributions(sig_data, bkg_data, category, outdir="Outputs"):
    """Compare feature distributions between signal and background"""
    os.makedirs(outdir, exist_ok=True)
    
    # Create comparison plots
    fig, axes = plt.subplots(4, 5, figsize=(20, 16))
    axes = axes.flatten()
    
    for i, var in enumerate(VARIABLES):
        ax = axes[i]
        
        sig_vals = sig_data[var]
        bkg_vals = bkg_data[var]
        
        # Remove NaN/Inf
        sig_vals = sig_vals[np.isfinite(sig_vals)]
        bkg_vals = bkg_vals[np.isfinite(bkg_vals)]
        
        # Normalize histograms
        ax.hist(sig_vals, bins=50, alpha=0.6, label="Signal", density=True, color="blue")
        ax.hist(bkg_vals, bins=50, alpha=0.6, label=f"{category}", density=True, color="red")
        
        # Calculate overlap metric (integral of min(hist1, hist2))
        sig_h, _ = np.histogram(sig_vals, bins=50, density=True)
        bkg_h, _ = np.histogram(bkg_vals, bins=50, density=True)
        overlap = np.sum(np.minimum(sig_h, bkg_h)) / len(sig_h)
        
        ax.set_title(f"{var}\nOverlap: {overlap:.3f}")
        ax.set_yscale("log")
        ax.legend()
    
    # Remove empty subplots
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])
    
    plt.tight_layout()
    outfile = os.path.join(outdir, f"compare_{category}_vs_signal.png")
    plt.savefig(outfile, dpi=100, bbox_inches="tight")
    print(f"[INFO] Saved: {outfile}")
    plt.close()

def compute_separability_metrics(sig_data, bkg_data, category):
    """Compute metrics to quantify signal/background separability"""
    print(f"\n{'='*70}")
    print(f"SEPARABILITY METRICS FOR {category} vs SIGNAL")
    print(f"{'='*70}")
    
    print(f"\nEvent counts:")
    print(f"  Signal:     {len(sig_data[VARIABLES[0]])}")
    print(f"  {category:8s}: {len(bkg_data[VARIABLES[0]])}")
    
    # Check for feature overlap via Kolmogorov-Smirnov test
    from scipy.stats import ks_2samp
    
    print(f"\nFeature separability (via KS test p-value):")
    print(f"  p-value > 0.05 suggests distributions are same (BAD for BDT)")
    print(f"  p-value < 0.05 suggests distributions are different (GOOD)\n")
    
    separability_scores = []
    for var in VARIABLES:
        sig_vals = sig_data[var]
        bkg_vals = bkg_data[var]
        
        # Remove NaN/Inf
        sig_vals = sig_vals[np.isfinite(sig_vals)]
        bkg_vals = bkg_vals[np.isfinite(bkg_vals)]
        
        if len(sig_vals) == 0 or len(bkg_vals) == 0:
            print(f"  {var:20s}: SKIPPED (no valid data)")
            continue
        
        stat, pval = ks_2samp(sig_vals, bkg_vals)
        separability_scores.append((var, stat, pval))
        
        # Color-code output
        if pval < 0.001:
            marker = "✓✓✓"  # Very good
        elif pval < 0.05:
            marker = "✓"    # Good
        else:
            marker = "✗"    # Poor
        
        print(f"  {var:20s}: stat={stat:.4f}  p={pval:.2e}  {marker}")
    
    # Summary
    poor_features = sum(1 for _, _, p in separability_scores if p > 0.05)
    print(f"\nSummary: {poor_features}/{len(VARIABLES)} features show poor separation (p > 0.05)")
    
    if poor_features > len(VARIABLES) / 2:
        print("⚠️  WARNING: More than half the features cannot separate Signal from Background!")
        print("   This explains the poor BDT performance.")
    
    return separability_scores

def check_data_integrity(sig_data, bkg_data, category):
    """Check for data quality issues"""
    print(f"\n{'='*70}")
    print(f"DATA QUALITY CHECK FOR {category}")
    print(f"{'='*70}\n")
    
    issues = []
    
    for var in VARIABLES:
        sig_vals = sig_data[var]
        bkg_vals = bkg_data[var]
        
        sig_nan = np.isnan(sig_vals).sum()
        sig_inf = np.isinf(sig_vals).sum()
        bkg_nan = np.isnan(bkg_vals).sum()
        bkg_inf = np.isinf(bkg_vals).sum()
        
        if sig_nan > 0 or sig_inf > 0 or bkg_nan > 0 or bkg_inf > 0:
            print(f"  {var:20s}: Signal NaN={sig_nan} Inf={sig_inf}  |  Bkg NaN={bkg_nan} Inf={bkg_inf}")
            issues.append(var)
    
    if issues:
        print(f"\n⚠️  WARNING: {len(issues)} variables have NaN/Inf values")
        print(f"   Variables: {', '.join(issues)}")
    else:
        print("✓ All variables appear valid (no NaN/Inf)")
    
    # Check for constant features
    print("\nConstant feature check:")
    for var in VARIABLES:
        sig_std = np.std(sig_data[var][np.isfinite(sig_data[var])])
        bkg_std = np.std(bkg_data[var][np.isfinite(bkg_data[var])])
        
        if sig_std < 1e-6 or bkg_std < 1e-6:
            print(f"  {var:20s}: CONSTANT or near-constant! (Signal std={sig_std:.2e}, Bkg std={bkg_std:.2e})")
            issues.append(var)
    
    if not issues:
        print("✓ No constant features detected")

def main():
    args = parse_args()
    
    print(f"[diagnose_ewk] Loading data for era {args.era} from {args.dir}")
    
    # Load EWK data
    result = load_category_data(args.dir, args.era, "EWK")
    if result is None:
        print(f"[ERROR] Could not load EWK data")
        return
    
    sig_data, bkg_data, n_sig, n_bkg = result
    
    print(f"[diagnose_ewk] Loaded {n_sig} signal and {n_bkg} EWK events")
    
    # Perform diagnostics
    check_data_integrity(sig_data, bkg_data, "EWK")
    separability_scores = compute_separability_metrics(sig_data, bkg_data, "EWK")
    compare_distributions(sig_data, bkg_data, "EWK")
    
    # Summary recommendations
    print(f"\n{'='*70}")
    print("RECOMMENDED NEXT STEPS:")
    print(f"{'='*70}")
    
    poor_features = sum(1 for _, _, p in separability_scores if p > 0.05)
    if poor_features > len(VARIABLES) / 2:
        print("1. EWK is physically too similar to Signal for this variable set")
        print("2. Consider:")
        print("   - Adding more discriminating variables (e.g., MET, collimation, etc.)")
        print("   - Training separate models for DY/EWK vs others")
        print("   - Using a different approach for EWK estimation")
    else:
        print("1. Variables seem separable; problem may be in training")
        print("2. Check:")
        print("   - Training set size and balance")
        print("   - Model hyperparameters and regularization")
        print("   - Whether EWK training data is complete/valid")

if __name__ == "__main__":
    main()
