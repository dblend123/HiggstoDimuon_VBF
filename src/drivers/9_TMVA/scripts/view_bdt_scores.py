#!/usr/bin/env python3
# view_bdt_scores.py
# analysis of BDT scores

import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os

def parse_args():
    p = argparse.ArgumentParser(description="Analyze BDT scores with detailed diagnostics")
    p.add_argument("--input", required=True, help="Scored ROOT file (combined signal+background)")
    p.add_argument("--treename", default="mva", help="Tree name")
    p.add_argument("--output", default="bdt_score_analysis.pdf", help="Output PDF")
    p.add_argument("--signal-file", help="Optional: separate signal-only file to compare")
    p.add_argument("--bkg-file", help="Optional: separate background-only file to compare")
    return p.parse_args()

def load_data(filename, treename):
    """Load BDT scores and other variables from ROOT file"""
    df = ROOT.RDataFrame(treename, filename)
    
    # Get available columns
    cols = list(df.GetColumnNames())
    
    # Essential columns
    if 'bdt_score' not in cols:
        raise ValueError(f"No 'bdt_score' branch in {filename}")
    
    # Load columns that exist
    cols_to_load = ['bdt_score']
    optional = ['weight', 'fold', 'm_mumu', 'pt_mumu', 'njets', 'z', 'pt_balance', 
                'deta_jj', 'dphi_mumu_jj', 'mu1_pt_over_m', 'mu2_pt_over_m']
    
    for col in optional:
        if col in cols:
            cols_to_load.append(col)
    
    arrays = df.AsNumpy(cols_to_load)
    
    # Convert to dict of numpy arrays
    data = {}
    for col in cols_to_load:
        data[col] = np.asarray(arrays[col], dtype=float)
    
    # Add weights if not present
    if 'weight' not in data:
        data['weight'] = np.ones(len(data['bdt_score']))
    
    return data, cols

def print_statistics(data, label=""):
    """Print detailed statistics about BDT scores"""
    scores = data['bdt_score']
    weights = data.get('weight', np.ones(len(scores)))
    
    print(f"\n{'='*80}")
    print(f"BDT SCORE STATISTICS {label}")
    print(f"{'='*80}")
    
    print(f"\nBasic stats:")
    print(f"  Events:     {len(scores):,}")
    print(f"  Mean:       {scores.mean():.6f}")
    print(f"  Weighted:   {np.average(scores, weights=weights):.6f}")
    print(f"  Std Dev:    {scores.std():.6f}")
    print(f"  Median:     {np.median(scores):.6f}")
    print(f"  Min:        {scores.min():.6f}")
    print(f"  Max:        {scores.max():.6f}")
    
    print(f"\nPercentiles:")
    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    for p in percentiles:
        print(f"  {p:2d}th:       {np.percentile(scores, p):.6f}")
    
    print(f"\nScore distribution:")
    bins = [(0.0, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 0.95), (0.95, 1.0)]
    for low, high in bins:
        mask = (scores >= low) & (scores < high)
        frac = mask.sum() / len(scores) * 100
        print(f"  [{low:.2f}, {high:.2f}): {mask.sum():8,} events ({frac:5.2f}%)")
    
    # Weight statistics if available
    if 'weight' in data:
        w = data['weight']
        print(f"\nWeight statistics:")
        print(f"  Sum:      {w.sum():.4e}")
        print(f"  Mean:     {w.mean():.4e}")
        print(f"  Negative: {(w < 0).sum():,} ({100*(w < 0).sum()/len(w):.2f}%)")
        print(f"  Zero:     {(w == 0).sum():,} ({100*(w == 0).sum()/len(w):.2f}%)")

def create_diagnostic_plots(data_combined, data_sig, data_bkg, output_pdf):
    """Create comprehensive diagnostic plots"""
    
    with PdfPages(output_pdf) as pdf:
        
        # ========== PAGE 1: Score distributions ==========
        fig = plt.figure(figsize=(11, 14))
        
        # Plot 1: Overall distribution (linear)
        ax1 = plt.subplot(4, 1, 1)
        scores_comb = data_combined['bdt_score']
        weights_comb = data_combined.get('weight', np.ones(len(scores_comb)))
        
        if data_sig is not None and data_bkg is not None:
            scores_sig = data_sig['bdt_score']
            scores_bkg = data_bkg['bdt_score']
            weights_sig = data_sig.get('weight', np.ones(len(scores_sig)))
            weights_bkg = data_bkg.get('weight', np.ones(len(scores_bkg)))
            
            ax1.hist(scores_sig, bins=50, range=(0, 1), weights=weights_sig,
                    alpha=0.5, label=f'Signal (n={len(scores_sig):,})', 
                    color='blue', density=True, histtype='stepfilled')
            ax1.hist(scores_bkg, bins=50, range=(0, 1), weights=weights_bkg,
                    alpha=0.5, label=f'Background (n={len(scores_bkg):,})', 
                    color='red', density=True, histtype='stepfilled')
            
            # Add vertical lines for means
            ax1.axvline(np.average(scores_sig, weights=weights_sig), 
                       color='blue', linestyle='--', linewidth=2, alpha=0.7,
                       label=f'Signal mean: {np.average(scores_sig, weights=weights_sig):.3f}')
            ax1.axvline(np.average(scores_bkg, weights=weights_bkg), 
                       color='red', linestyle='--', linewidth=2, alpha=0.7,
                       label=f'Bkg mean: {np.average(scores_bkg, weights=weights_bkg):.3f}')
        else:
            ax1.hist(scores_comb, bins=50, range=(0, 1), weights=weights_comb,
                    alpha=0.7, label=f'All events (n={len(scores_comb):,})', 
                    color='green', density=True)
        
        ax1.set_xlabel('BDT Score', fontsize=12)
        ax1.set_ylabel('Normalized Density', fontsize=12)
        ax1.set_title('BDT Score Distribution (Linear Scale)', fontsize=13, fontweight='bold')
        ax1.legend(fontsize=9, loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, 1)
        
        # Plot 2: Log scale
        ax2 = plt.subplot(4, 1, 2)
        if data_sig is not None and data_bkg is not None:
            ax2.hist(scores_sig, bins=50, range=(0, 1), weights=weights_sig,
                    alpha=0.5, label='Signal', color='blue', density=True, histtype='stepfilled')
            ax2.hist(scores_bkg, bins=50, range=(0, 1), weights=weights_bkg,
                    alpha=0.5, label='Background', color='red', density=True, histtype='stepfilled')
        else:
            ax2.hist(scores_comb, bins=50, range=(0, 1), weights=weights_comb,
                    alpha=0.7, label='All events', color='green', density=True)
        
        ax2.set_xlabel('BDT Score', fontsize=12)
        ax2.set_ylabel('Normalized Density (log)', fontsize=12)
        ax2.set_title('BDT Score Distribution (Log Scale)', fontsize=13, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(0, 1)
        ax2.set_yscale('log')
        ax2.set_ylim(bottom=1e-3)
        
        # Plot 3: Cumulative distributions
        ax3 = plt.subplot(4, 1, 3)
        if data_sig is not None and data_bkg is not None:
            # Signal efficiency vs threshold
            thresholds = np.linspace(0, 1, 200)
            sig_eff = [np.sum(weights_sig[scores_sig >= t]) / np.sum(weights_sig) for t in thresholds]
            bkg_eff = [np.sum(weights_bkg[scores_bkg >= t]) / np.sum(weights_bkg) for t in thresholds]
            
            ax3.plot(thresholds, sig_eff, 'b-', linewidth=2, label='Signal efficiency')
            ax3.plot(thresholds, bkg_eff, 'r-', linewidth=2, label='Background efficiency')
            ax3.plot(thresholds, np.array(sig_eff) / (np.array(bkg_eff) + 1e-10), 
                    'g--', linewidth=2, label='S/B improvement')
            
            # Mark interesting points
            ax3.axvline(0.5, color='gray', linestyle=':', alpha=0.5)
            ax3.axvline(0.865, color='purple', linestyle='--', linewidth=2, 
                       alpha=0.7, label='Optimal cut (0.865)')
            
            ax3.set_xlabel('BDT Score Threshold', fontsize=12)
            ax3.set_ylabel('Efficiency / S/B', fontsize=12)
            ax3.set_title('Signal/Background Efficiency vs Threshold', fontsize=13, fontweight='bold')
            ax3.legend(fontsize=9)
            ax3.grid(True, alpha=0.3)
            ax3.set_xlim(0, 1)
        else:
            # Just show cumulative distribution
            sorted_scores = np.sort(scores_comb)
            cumsum = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
            ax3.plot(sorted_scores, cumsum, 'g-', linewidth=2)
            ax3.set_xlabel('BDT Score', fontsize=12)
            ax3.set_ylabel('Cumulative Fraction', fontsize=12)
            ax3.set_title('Cumulative Distribution', fontsize=13, fontweight='bold')
            ax3.grid(True, alpha=0.3)
        
        # Plot 4: Score bins summary
        ax4 = plt.subplot(4, 1, 4)
        bins = [(0.0, 0.2, 'Very Bkg-like'), (0.2, 0.5, 'Bkg-enriched'), 
                (0.5, 0.8, 'Mixed'), (0.8, 0.95, 'Sig-enriched'), 
                (0.95, 1.0, 'Very Sig-like')]
        
        if data_sig is not None and data_bkg is not None:
            sig_counts = []
            bkg_counts = []
            labels = []
            
            for low, high, label in bins:
                sig_mask = (scores_sig >= low) & (scores_sig < high)
                bkg_mask = (scores_bkg >= low) & (scores_bkg < high)
                sig_counts.append(np.sum(weights_sig[sig_mask]))
                bkg_counts.append(np.sum(weights_bkg[bkg_mask]))
                labels.append(f'{label}\n[{low:.2f},{high:.2f})')
            
            x = np.arange(len(labels))
            width = 0.35
            
            ax4.bar(x - width/2, sig_counts, width, label='Signal', color='blue', alpha=0.7)
            ax4.bar(x + width/2, bkg_counts, width, label='Background', color='red', alpha=0.7)
            
            ax4.set_xlabel('BDT Score Range', fontsize=12)
            ax4.set_ylabel('Weighted Events', fontsize=12)
            ax4.set_title('Event Distribution by Score Range', fontsize=13, fontweight='bold')
            ax4.set_xticks(x)
            ax4.set_xticklabels(labels, fontsize=9)
            ax4.legend(fontsize=10)
            ax4.grid(True, alpha=0.3, axis='y')
            ax4.set_yscale('log')
        else:
            counts = []
            labels = []
            for low, high, label in bins:
                mask = (scores_comb >= low) & (scores_comb < high)
                counts.append(mask.sum())
                labels.append(f'{label}\n[{low:.2f},{high:.2f})')
            
            ax4.bar(labels, counts, color='green', alpha=0.7)
            ax4.set_xlabel('BDT Score Range', fontsize=12)
            ax4.set_ylabel('Events', fontsize=12)
            ax4.set_title('Event Distribution by Score Range', fontsize=13, fontweight='bold')
            ax4.tick_params(axis='x', labelsize=9)
            ax4.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # ========== PAGE 2: Score vs Physics Variables ==========
        if 'm_mumu' in data_combined and 'pt_mumu' in data_combined:
            fig = plt.figure(figsize=(11, 14))
            
            m_mumu = data_combined['m_mumu']
            pt_mumu = data_combined['pt_mumu']
            scores = data_combined['bdt_score']
            
            # Plot 1: Score vs m_mumu (2D histogram)
            ax1 = plt.subplot(3, 2, 1)
            h = ax1.hist2d(m_mumu, scores, bins=[100, 50], 
                          range=[[100, 150], [0, 1]], cmap='YlOrRd', cmin=1)
            plt.colorbar(h[3], ax=ax1, label='Events')
            ax1.axvline(125, color='blue', linestyle='--', linewidth=2, label='Higgs mass')
            ax1.set_xlabel('m_μμ [GeV]', fontsize=11)
            ax1.set_ylabel('BDT Score', fontsize=11)
            ax1.set_title('BDT Score vs Dimuon Mass', fontsize=12, fontweight='bold')
            ax1.legend(fontsize=9)
            ax1.grid(True, alpha=0.3)
            
            # Plot 2: Score in mass windows
            ax2 = plt.subplot(3, 2, 2)
            mass_windows = [
                (110, 120, 'Low [110-120]', 'orange'),
                (120, 130, 'Signal [120-130]', 'green'),
                (130, 140, 'High [130-140]', 'purple')
            ]
            for m_low, m_high, label, color in mass_windows:
                mask = (m_mumu >= m_low) & (m_mumu < m_high)
                if mask.sum() > 10:
                    ax2.hist(scores[mask], bins=30, range=(0, 1), 
                           alpha=0.5, label=label, density=True, 
                           histtype='stepfilled', color=color)
            ax2.set_xlabel('BDT Score', fontsize=11)
            ax2.set_ylabel('Normalized Density', fontsize=11)
            ax2.set_title('BDT Score in Different Mass Windows', fontsize=12, fontweight='bold')
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)
            
            # Plot 3: Mean score vs mass
            ax3 = plt.subplot(3, 2, 3)
            mass_bins = np.linspace(100, 150, 26)
            mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
            mean_scores = []
            for i in range(len(mass_bins) - 1):
                mask = (m_mumu >= mass_bins[i]) & (m_mumu < mass_bins[i+1])
                if mask.sum() > 0:
                    mean_scores.append(scores[mask].mean())
                else:
                    mean_scores.append(np.nan)
            ax3.plot(mass_centers, mean_scores, 'o-', linewidth=2, markersize=5)
            ax3.axvline(125, color='blue', linestyle='--', linewidth=2, alpha=0.5)
            ax3.set_xlabel('m_μμ [GeV]', fontsize=11)
            ax3.set_ylabel('Mean BDT Score', fontsize=11)
            ax3.set_title('Mean BDT Score vs Mass', fontsize=12, fontweight='bold')
            ax3.grid(True, alpha=0.3)
            
            # Plot 4: Score vs pT
            ax4 = plt.subplot(3, 2, 4)
            h = ax4.hist2d(pt_mumu, scores, bins=[50, 50], 
                          range=[[0, 300], [0, 1]], cmap='YlOrRd', cmin=1)
            plt.colorbar(h[3], ax=ax4, label='Events')
            ax4.set_xlabel('pT_μμ [GeV]', fontsize=11)
            ax4.set_ylabel('BDT Score', fontsize=11)
            ax4.set_title('BDT Score vs Dimuon pT', fontsize=12, fontweight='bold')
            ax4.grid(True, alpha=0.3)
            
            # Plot 5: Score distribution by pT bins
            ax5 = plt.subplot(3, 2, 5)
            pt_windows = [
                (0, 50, 'Low pT [0-50]', 'red'),
                (50, 100, 'Med pT [50-100]', 'orange'),
                (100, 200, 'High pT [100-200]', 'green'),
                (200, 500, 'Very High pT [200+]', 'blue')
            ]
            for pt_low, pt_high, label, color in pt_windows:
                mask = (pt_mumu >= pt_low) & (pt_mumu < pt_high)
                if mask.sum() > 10:
                    ax5.hist(scores[mask], bins=30, range=(0, 1), 
                           alpha=0.5, label=label, density=True, 
                           histtype='stepfilled', color=color)
            ax5.set_xlabel('BDT Score', fontsize=11)
            ax5.set_ylabel('Normalized Density', fontsize=11)
            ax5.set_title('BDT Score by pT Range', fontsize=12, fontweight='bold')
            ax5.legend(fontsize=9)
            ax5.grid(True, alpha=0.3)
            
            # Plot 6: njets if available
            ax6 = plt.subplot(3, 2, 6)
            if 'njets' in data_combined:
                njets = data_combined['njets']
                for nj in [2, 3, 4, 5]:
                    mask = (njets == nj)
                    if mask.sum() > 10:
                        ax6.hist(scores[mask], bins=30, range=(0, 1), 
                               alpha=0.5, label=f'{nj} jets (n={mask.sum():,})', 
                               density=True, histtype='stepfilled')
                ax6.set_xlabel('BDT Score', fontsize=11)
                ax6.set_ylabel('Normalized Density', fontsize=11)
                ax6.set_title('BDT Score by Jet Multiplicity', fontsize=12, fontweight='bold')
                ax6.legend(fontsize=9)
                ax6.grid(True, alpha=0.3)
            else:
                ax6.text(0.5, 0.5, 'njets not available', 
                        ha='center', va='center', fontsize=14)
                ax6.axis('off')
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
        
        # ========== PAGE 3: Feature correlations with score ==========
        fig = plt.figure(figsize=(11, 14))
        
        feature_vars = ['pt_mumu', 'z', 'pt_balance', 'deta_jj', 'dphi_mumu_jj']
        available_features = [f for f in feature_vars if f in data_combined]
        
        if len(available_features) >= 4:
            for idx, feat in enumerate(available_features[:6]):
                ax = plt.subplot(3, 2, idx + 1)
                
                feat_vals = data_combined[feat]
                scores = data_combined['bdt_score']
                
                # 2D histogram
                h = ax.hist2d(feat_vals, scores, bins=[50, 50], cmap='YlOrRd', cmin=1)
                plt.colorbar(h[3], ax=ax, label='Events')
                ax.set_xlabel(feat, fontsize=11)
                ax.set_ylabel('BDT Score', fontsize=11)
                ax.set_title(f'BDT Score vs {feat}', fontsize=11, fontweight='bold')
                ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
        
        print(f"✓ Created diagnostic PDF: {output_pdf}")

def main():
    args = parse_args()
    
    print("="*80)
    print("ENHANCED BDT SCORE ANALYSIS")
    print("="*80)
    
    # Load combined file
    print(f"\nLoading combined file: {args.input}")
    data_combined, cols_combined = load_data(args.input, args.treename)
    print(f"  Loaded {len(data_combined['bdt_score']):,} events")
    print(f"  Available columns: {', '.join(sorted([str(c) for c in cols_combined]))}")
    
    # Print statistics for combined
    print_statistics(data_combined, "(Combined)")
    
    # Optionally load separate signal/background files
    data_sig = None
    data_bkg = None
    
    if args.signal_file and os.path.exists(args.signal_file):
        print(f"\nLoading signal file: {args.signal_file}")
        data_sig, _ = load_data(args.signal_file, args.treename)
        print_statistics(data_sig, "(Signal only)")
    
    if args.bkg_file and os.path.exists(args.bkg_file):
        print(f"\nLoading background file: {args.bkg_file}")
        data_bkg, _ = load_data(args.bkg_file, args.treename)
        print_statistics(data_bkg, "(Background only)")
    
    # Create plots
    print(f"\n{'='*80}")
    print("CREATING DIAGNOSTIC PLOTS")
    print("="*80)
    create_diagnostic_plots(data_combined, data_sig, data_bkg, args.output)
    
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("="*80)
    print(f"\n✓ Analysis complete!")
    print(f"✓ Plots saved to: {args.output}")
    print(f"\nNext steps:")
    print(f"  1. Review the multi-page PDF with comprehensive diagnostics")
    print(f"  2. Check if BDT score peaks at signal mass (125 GeV)")
    print(f"  3. Verify high-pT events have higher scores (VBF signature)")
    print(f"  4. Apply score cuts for your analysis")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()