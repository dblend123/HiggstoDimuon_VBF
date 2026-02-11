#!/usr/bin/env python3
# view_bdt_scores.py
# View and analyze BDT scores from the scored ROOT file

import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="View BDT scores from scored ROOT file")
    p.add_argument("--input", required=True, help="Scored ROOT file")
    p.add_argument("--treename", default="mva", help="Tree name")
    p.add_argument("--output", default="bdt_score_analysis.pdf", help="Output PDF")
    p.add_argument("--n-show", type=int, default=20, help="Number of events to print")
    return p.parse_args()

def main():
    args = parse_args()
    
    print("="*80)
    print("BDT SCORE ANALYSIS")
    print("="*80)
    print(f"\nInput file: {args.input}")
    
    # Open file
    f = ROOT.TFile.Open(args.input)
    if not f or f.IsZombie():
        print(f"ERROR: Cannot open {args.input}")
        return
    
    tree = f.Get(args.treename)
    if not tree:
        print(f"ERROR: Cannot find tree '{args.treename}'")
        f.Close()
        return
    
    n_entries = tree.GetEntries()
    print(f"Total entries: {n_entries:,}")
    
    # Check if bdt_score exists
    branches = [b.GetName() for b in tree.GetListOfBranches()]
    if 'bdt_score' not in branches:
        print("ERROR: No 'bdt_score' branch found!")
        print(f"Available branches: {branches[:10]}...")
        f.Close()
        return
    
    print("\n✓ Found 'bdt_score' branch")
    
    # Load data using RDataFrame for efficiency
    df = ROOT.RDataFrame(args.treename, args.input)
    
    # Get columns we need
    cols_to_get = ['bdt_score']
    
    # Add optional columns if they exist
    optional_cols = ['weight', 'isSignal', 'fold', 'm_mumu', 'pt_mumu', 'njets']
    for col in optional_cols:
        if col in branches:
            cols_to_get.append(col)
    
    print(f"Loading columns: {cols_to_get}")
    
    arrays = df.AsNumpy(cols_to_get)
    scores = np.asarray(arrays['bdt_score'], dtype=float)
    
    # Get weights if available
    if 'weight' in arrays:
        weights = np.asarray(arrays['weight'], dtype=float)
        has_weights = True
    else:
        weights = np.ones(len(scores))
        has_weights = False
    
    # Get signal/background labels if available
    if 'isSignal' in arrays:
        labels = np.asarray(arrays['isSignal'], dtype=int)
        has_labels = True
    else:
        labels = None
        has_labels = False
    
    # Print statistics
    print("\n" + "="*80)
    print("BDT SCORE STATISTICS")
    print("="*80)
    print(f"\nOverall statistics:")
    print(f"  Count:      {len(scores):,}")
    print(f"  Mean:       {scores.mean():.6f}")
    print(f"  Std Dev:    {scores.std():.6f}")
    print(f"  Min:        {scores.min():.6f}")
    print(f"  Max:        {scores.max():.6f}")
    print(f"  Median:     {np.median(scores):.6f}")
    print(f"\nPercentiles:")
    for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
        print(f"  {p:2d}th:       {np.percentile(scores, p):.6f}")
    
    if has_labels:
        sig_mask = labels == 1
        bkg_mask = labels == 0
        
        print(f"\n" + "-"*80)
        print("SIGNAL vs BACKGROUND")
        print("-"*80)
        print(f"\nSignal events: {sig_mask.sum():,}")
        print(f"  Mean score: {scores[sig_mask].mean():.6f}")
        print(f"  Std Dev:    {scores[sig_mask].std():.6f}")
        print(f"  Median:     {np.median(scores[sig_mask]):.6f}")
        
        print(f"\nBackground events: {bkg_mask.sum():,}")
        print(f"  Mean score: {scores[bkg_mask].mean():.6f}")
        print(f"  Std Dev:    {scores[bkg_mask].std():.6f}")
        print(f"  Median:     {np.median(scores[bkg_mask]):.6f}")
        
        # Calculate separation
        separation = abs(scores[sig_mask].mean() - scores[bkg_mask].mean())
        print(f"\nSeparation (|mean_sig - mean_bkg|): {separation:.6f}")
    
    # Show some example events
    print(f"\n" + "="*80)
    print(f"FIRST {args.n_show} EVENTS")
    print("="*80)
    
    tree.GetEntry(0)  # Reset to first entry
    
    print(f"\n{'Event':<8s} {'BDT Score':<12s}", end="")
    if has_labels:
        print(f" {'Label':<8s}", end="")
    if has_weights:
        print(f" {'Weight':<12s}", end="")
    if 'm_mumu' in arrays:
        print(f" {'m_mumu':<10s}", end="")
    if 'pt_mumu' in arrays:
        print(f" {'pt_mumu':<10s}", end="")
    if 'njets' in arrays:
        print(f" {'njets':<6s}", end="")
    print()
    print("-"*80)
    
    for i in range(min(args.n_show, n_entries)):
        tree.GetEntry(i)
        print(f"{i:<8d} {tree.bdt_score:<12.6f}", end="")
        if has_labels:
            label_str = "Signal" if labels[i] == 1 else "Bkg"
            print(f" {label_str:<8s}", end="")
        if has_weights:
            print(f" {weights[i]:<12.6e}", end="")
        if 'm_mumu' in arrays:
            print(f" {arrays['m_mumu'][i]:<10.3f}", end="")
        if 'pt_mumu' in arrays:
            print(f" {arrays['pt_mumu'][i]:<10.3f}", end="")
        if 'njets' in arrays:
            print(f" {int(arrays['njets'][i]):<6d}", end="")
        print()
    
    # Create plots
    print(f"\n" + "="*80)
    print(f"CREATING PLOTS -> {args.output}")
    print("="*80)
    
    with PdfPages(args.output) as pdf:
        
        # Page 1: Overall score distribution
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))
        
        # Linear scale
        ax = axes[0]
        if has_labels:
            ax.hist(scores[sig_mask], bins=50, range=(0, 1), 
                   weights=weights[sig_mask], alpha=0.5, 
                   label=f'Signal (n={sig_mask.sum():,})', color='blue', density=True)
            ax.hist(scores[bkg_mask], bins=50, range=(0, 1), 
                   weights=weights[bkg_mask], alpha=0.5, 
                   label=f'Background (n={bkg_mask.sum():,})', color='red', density=True)
        else:
            ax.hist(scores, bins=50, range=(0, 1), weights=weights, 
                   alpha=0.7, label=f'All events (n={len(scores):,})', color='green', density=True)
        
        ax.set_xlabel('BDT Score', fontsize=12)
        ax.set_ylabel('Normalized Events', fontsize=12)
        ax.set_title('BDT Score Distribution', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        # Log scale
        ax = axes[1]
        if has_labels:
            ax.hist(scores[sig_mask], bins=50, range=(0, 1), 
                   weights=weights[sig_mask], alpha=0.5, 
                   label=f'Signal', color='blue', density=True)
            ax.hist(scores[bkg_mask], bins=50, range=(0, 1), 
                   weights=weights[bkg_mask], alpha=0.5, 
                   label=f'Background', color='red', density=True)
        else:
            ax.hist(scores, bins=50, range=(0, 1), weights=weights, 
                   alpha=0.7, label=f'All events', color='green', density=True)
        
        ax.set_xlabel('BDT Score', fontsize=12)
        ax.set_ylabel('Normalized Events (log scale)', fontsize=12)
        ax.set_title('BDT Score Distribution (Log Scale)', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_yscale('log')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # Page 2: Score vs other variables
        if 'm_mumu' in arrays and 'pt_mumu' in arrays:
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            
            m_mumu = np.asarray(arrays['m_mumu'], dtype=float)
            pt_mumu = np.asarray(arrays['pt_mumu'], dtype=float)
            
            # Score vs m_mumu
            ax = axes[0, 0]
            if has_labels:
                ax.scatter(m_mumu[sig_mask], scores[sig_mask], 
                          alpha=0.3, s=1, label='Signal', c='blue')
                ax.scatter(m_mumu[bkg_mask], scores[bkg_mask], 
                          alpha=0.3, s=1, label='Background', c='red')
            else:
                ax.scatter(m_mumu, scores, alpha=0.3, s=1, c='green')
            ax.set_xlabel('m_μμ [GeV]', fontsize=10)
            ax.set_ylabel('BDT Score', fontsize=10)
            ax.set_title('BDT Score vs Dimuon Mass', fontsize=11)
            if has_labels:
                ax.legend(fontsize=8, markerscale=5)
            ax.grid(True, alpha=0.3)
            
            # Score vs pt_mumu
            ax = axes[0, 1]
            if has_labels:
                ax.scatter(pt_mumu[sig_mask], scores[sig_mask], 
                          alpha=0.3, s=1, label='Signal', c='blue')
                ax.scatter(pt_mumu[bkg_mask], scores[bkg_mask], 
                          alpha=0.3, s=1, label='Background', c='red')
            else:
                ax.scatter(pt_mumu, scores, alpha=0.3, s=1, c='green')
            ax.set_xlabel('pT_μμ [GeV]', fontsize=10)
            ax.set_ylabel('BDT Score', fontsize=10)
            ax.set_title('BDT Score vs Dimuon pT', fontsize=11)
            if has_labels:
                ax.legend(fontsize=8, markerscale=5)
            ax.grid(True, alpha=0.3)
            
            # m_mumu distribution colored by score
            ax = axes[1, 0]
            sc = ax.scatter(m_mumu, np.random.uniform(-0.1, 0.1, len(m_mumu)), 
                          c=scores, s=5, alpha=0.5, cmap='RdYlBu_r', vmin=0, vmax=1)
            ax.set_xlabel('m_μμ [GeV]', fontsize=10)
            ax.set_ylabel('Random (for visibility)', fontsize=10)
            ax.set_title('Dimuon Mass colored by BDT Score', fontsize=11)
            plt.colorbar(sc, ax=ax, label='BDT Score')
            ax.grid(True, alpha=0.3)
            
            # Score distribution in mass windows
            ax = axes[1, 1]
            mass_windows = [
                (110, 120, 'Low mass'),
                (120, 130, 'Signal region'),
                (130, 140, 'High mass')
            ]
            for m_low, m_high, label in mass_windows:
                mask = (m_mumu >= m_low) & (m_mumu < m_high)
                if mask.sum() > 0:
                    ax.hist(scores[mask], bins=30, range=(0, 1), 
                           alpha=0.5, label=f'{label} ({m_low}-{m_high} GeV)', 
                           density=True, histtype='step', linewidth=2)
            ax.set_xlabel('BDT Score', fontsize=10)
            ax.set_ylabel('Normalized Events', fontsize=10)
            ax.set_title('BDT Score in Different Mass Windows', fontsize=11)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
    
    print(f"✓ Saved plots to {args.output}")
    
    f.Close()
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"\n✓ BDT scores successfully loaded from: {args.input}")
    print(f"✓ Score range: [{scores.min():.6f}, {scores.max():.6f}]")
    print(f"✓ Mean score: {scores.mean():.6f}")
    if has_labels:
        print(f"✓ Signal mean: {scores[sig_mask].mean():.6f}")
        print(f"✓ Background mean: {scores[bkg_mask].mean():.6f}")
    print(f"✓ Plots saved to: {args.output}")
    print("\nNext steps:")
    print("  1. Check the plots to verify score distributions look reasonable")
    print("  2. Apply score cuts for your analysis (e.g., score > 0.8)")
    print("  3. Use scores for further analysis or plotting")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()