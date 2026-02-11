#!/usr/bin/env python3
# comprehensive_hmumu_analysis.py
# Complete H→μμ VBF analysis with BDT categories, fits, and limits

# Basic run (without individual background components)
# python3 comprehensive_hmumu_analysis.py \
#   --data Outputs/2024/2024_data_scored.root \
#   --signal Outputs/2024/2024_scored_signal.root \
#   --background Outputs/2024/2024_scored_background.root \
#   --treename mva \
#   --output hmumu_full_analysis.pdf \
#   --lumi 104.80585 \
#   --signal-xsec 0.000909568

# # Advanced run (with individual background components for stacking)
# python3 comprehensive_hmumu_analysis.py \
#   --data Outputs/2024/2024_data_scored.root \
#   --signal Outputs/2024/2024_scored_signal.root \
#   --background Outputs/2024/2024_scored_background.root \
#   --bkg-components Outputs/2024/2024_scored_DY.root \
#                     Outputs/2024/2024_scored_TT.root \
#                     Outputs/2024/2024_scored_VV.root \
#                     Outputs/2024/2024_scored_VVV.root \
#   --treename mva \
#   --output hmumu_full_analysis.pdf \
#   --lumi 104.80585 \
#   --signal-xsec 0.000909568




import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from scipy.optimize import curve_fit
import argparse
import os

# Suppress ROOT warnings
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

def parse_args():
    p = argparse.ArgumentParser(description="Comprehensive H→μμ analysis")
    p.add_argument("--data", required=True, help="Scored data file")
    p.add_argument("--signal", required=True, help="Scored signal MC file")
    p.add_argument("--background", required=True, help="Scored background MC file")
    p.add_argument("--bkg-components", nargs="+", 
                   help="Individual background component files (DY, TT, VV, VVV, EWK)")
    p.add_argument("--treename", default="mva", help="Tree name")
    p.add_argument("--output", default="hmumu_full_analysis.pdf", help="Output PDF")
    p.add_argument("--lumi", type=float, default=1.0, help="Integrated luminosity in fb^-1")
    p.add_argument("--signal-xsec", type=float, default=1.0, 
                   help="Signal cross-section in pb")
    return p.parse_args()

def load_data(filename, treename):
    """Load scored events with all necessary branches"""
    df = ROOT.RDataFrame(treename, filename)
    cols = list(df.GetColumnNames())
    
    # Load essential branches
    branches = ['bdt_score', 'm_mumu']
    optional = ['weight', 'pt_mumu', 'njets', 'z', 'pt_balance', 'deta_jj']
    
    for b in optional:
        if b in cols:
            branches.append(b)
    
    arrays = df.AsNumpy(branches)
    
    data = {}
    for b in branches:
        data[b] = np.asarray(arrays[b], dtype=float)
    
    # Ensure weights exist
    if 'weight' not in data:
        data['weight'] = np.ones(len(data['bdt_score']))
    
    return data

def define_categories():
    """Define BDT score categories"""
    return {
        'tight': (0.95, 1.00, 'VBF Tight'),
        'loose': (0.80, 0.95, 'VBF Loose'),
        'fail': (0.00, 0.80, 'Fail (CR)')
    }

def signal_model(m, norm, mean, sigma):
    """Crystal Ball or Double-sided Crystal Ball for signal"""
    # Simplified Gaussian for now - can be upgraded to CB
    return norm * np.exp(-0.5 * ((m - mean) / sigma)**2)

def background_model(m, a0, a1, a2):
    """Exponential or polynomial background"""
    # Bernstein polynomial (more stable than exponential)
    return a0 + a1*m + a2*m**2

def combined_model(m, sig_norm, sig_mean, sig_sigma, bkg_a0, bkg_a1, bkg_a2):
    """Signal + Background model"""
    return signal_model(m, sig_norm, sig_mean, sig_sigma) + \
           background_model(m, bkg_a0, bkg_a1, bkg_a2)

def fit_mass_distribution(m_data, weights_data, m_sig, weights_sig, 
                          mass_range=(110, 150), signal_region=(120, 130)):
    """
    Fit signal+background to data
    Returns: fit parameters, uncertainties, chi2, significance
    """
    
    # Create histograms
    bins = np.linspace(mass_range[0], mass_range[1], 41)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    
    data_hist, _ = np.histogram(m_data, bins=bins, weights=weights_data)
    data_errors = np.sqrt(np.histogram(m_data, bins=bins, weights=weights_data**2)[0])
    
    # Estimate signal shape from MC
    sig_hist, _ = np.histogram(m_sig, bins=bins, weights=weights_sig)
    
    # Initial parameter guesses
    sig_mean_init = 125.0
    sig_sigma_init = 2.0
    sig_norm_init = np.sum(data_hist[(bin_centers >= signal_region[0]) & 
                                      (bin_centers < signal_region[1])]) * 0.1
    
    # Background estimated from sidebands
    sideband_mask = ((bin_centers < signal_region[0]) | 
                     (bin_centers > signal_region[1]))
    bkg_init = np.mean(data_hist[sideband_mask])
    
    p0 = [sig_norm_init, sig_mean_init, sig_sigma_init, bkg_init, 0, 0]
    
    # Fit
    try:
        # Mask for fit (use full range or exclude blinded region if needed)
        fit_mask = np.ones(len(bin_centers), dtype=bool)
        
        popt, pcov = curve_fit(
            combined_model, 
            bin_centers[fit_mask], 
            data_hist[fit_mask],
            p0=p0,
            sigma=data_errors[fit_mask] + 1e-10,
            absolute_sigma=True,
            maxfev=10000
        )
        
        perr = np.sqrt(np.diag(pcov))
        
        # Calculate chi2
        y_fit = combined_model(bin_centers, *popt)
        chi2 = np.sum(((data_hist - y_fit) / (data_errors + 1e-10))**2)
        ndf = len(bin_centers) - len(popt)
        
        # Significance (simple Z = S/sqrt(B) in signal region)
        sr_mask = (bin_centers >= signal_region[0]) & (bin_centers < signal_region[1])
        S = np.sum(signal_model(bin_centers[sr_mask], *popt[:3]))
        B = np.sum(background_model(bin_centers[sr_mask], *popt[3:]))
        
        if B > 0:
            significance = S / np.sqrt(B)
        else:
            significance = 0.0
        
        return {
            'params': popt,
            'errors': perr,
            'chi2': chi2,
            'ndf': ndf,
            'significance': significance,
            'signal_yield': S,
            'background_yield': B
        }
        
    except Exception as e:
        print(f"    Fit failed: {e}")
        return None

def calculate_systematics(data_sig, data_bkg, categories):
    """
    Calculate systematic uncertainties
    """
    systematics = {}
    
    # 1. BDT modeling uncertainty (vary training/test split)
    # Approximate as variation in yields between folds
    if 'fold' in data_sig:
        fold_variations = []
        for cat_name, (low, high, label) in categories.items():
            if cat_name == 'fail':
                continue
            
            yields_by_fold = []
            for fold in [0, 1, 2, 3]:
                mask = ((data_sig['fold'] == fold) & 
                       (data_sig['bdt_score'] >= low) & 
                       (data_sig['bdt_score'] < high))
                yields_by_fold.append(np.sum(data_sig['weight'][mask]))
            
            if len(yields_by_fold) > 1:
                variation = np.std(yields_by_fold) / np.mean(yields_by_fold)
                fold_variations.append(variation)
        
        systematics['bdt_modeling'] = np.mean(fold_variations) if fold_variations else 0.05
    else:
        systematics['bdt_modeling'] = 0.05  # 5% default
    
    # 2. MC statistical uncertainty
    systematics['mc_stats'] = 0.02  # 2% (compute properly from MC weights)
    
    # 3. Luminosity
    systematics['luminosity'] = 0.025  # 2.5% (typical for CMS)
    
    # 4. Lepton efficiency (trigger + ID + isolation)
    systematics['lepton_eff'] = 0.03  # 3%
    
    # 5. Background normalization (from control region)
    systematics['background_norm'] = 0.10  # 10%
    
    # 6. Signal acceptance (PDF, scale variations)
    systematics['signal_acceptance'] = 0.05  # 5%
    
    # Total systematic
    systematics['total'] = np.sqrt(sum(v**2 for k, v in systematics.items() 
                                       if k != 'total'))
    
    return systematics

def calculate_limits(observed, expected_sig, expected_bkg, systematics):
    """
    Calculate 95% CL upper limit on signal strength μ
    Using simplified CLs method
    """
    
    # Observed limit (Poisson with background)
    # For observed events N, background B, signal S:
    # Find μ such that P(N | μS + B) = 0.05
    
    # Simple approach: Profile likelihood ratio
    # -2ln(λ) = -2ln(L(μ)/L(μ_hat))
    
    # For now, use simple frequentist approach
    # μ_95 such that μ*S = 1.96 * sqrt(B + (σ_sys*B)²)
    
    total_bkg_unc = expected_bkg * systematics['total']
    
    if expected_sig > 0:
        # Observed limit
        obs_limit = 1.96 * np.sqrt(observed + total_bkg_unc**2) / expected_sig
        
        # Expected limit (median)
        exp_limit = 1.96 * np.sqrt(expected_bkg + total_bkg_unc**2) / expected_sig
        
        # Expected ±1σ bands
        exp_limit_p1 = 1.96 * np.sqrt(expected_bkg + 1.0*np.sqrt(expected_bkg) + total_bkg_unc**2) / expected_sig
        exp_limit_m1 = max(0, 1.96 * np.sqrt(expected_bkg - 1.0*np.sqrt(expected_bkg) + total_bkg_unc**2) / expected_sig)
        
        return {
            'observed': obs_limit,
            'expected': exp_limit,
            'expected_p1sigma': exp_limit_p1,
            'expected_m1sigma': exp_limit_m1
        }
    else:
        return None

def create_comprehensive_pdf(data_data, data_sig, data_bkg, bkg_components, 
                             categories, args, output_pdf):
    """
    Create comprehensive multi-page analysis PDF
    """
    
    with PdfPages(output_pdf) as pdf:
        
        # ========== PAGE 1: Title & Summary ==========
        fig = plt.figure(figsize=(11, 8.5))
        fig.text(0.5, 0.9, 'VBF H→μμ Analysis', 
                ha='center', fontsize=24, fontweight='bold')
        fig.text(0.5, 0.85, f'Luminosity: {args.lumi:.2f} fb⁻¹', 
                ha='center', fontsize=16)
        
        summary_text = f"""
Analysis Overview
{'='*60}

Dataset:
  • Data events: {len(data_data['bdt_score']):,}
  • Signal MC: {len(data_sig['bdt_score']):,}
  • Background MC: {len(data_bkg['bdt_score']):,}

BDT Categories:
"""
        for cat_name, (low, high, label) in categories.items():
            n_data = np.sum((data_data['bdt_score'] >= low) & 
                           (data_data['bdt_score'] < high))
            summary_text += f"  • {label:15s} [{low:.2f}, {high:.2f}): {n_data:,} events\n"
        
        fig.text(0.1, 0.65, summary_text, fontsize=11, family='monospace',
                verticalalignment='top')
        
        plt.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # ========== PAGE 2: BDT Score Distributions ==========
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        scores_data = data_data['bdt_score']
        scores_sig = data_sig['bdt_score']
        scores_bkg = data_bkg['bdt_score']
        weights_sig = data_sig['weight']
        weights_bkg = data_bkg['weight']
        
        # Normalize MC to data
        scale_sig = args.lumi * args.signal_xsec * 1000  # pb to fb
        scale_bkg = args.lumi  # Assume background already has proper weights
        
        # Plot 1: BDT comparison
        ax = axes[0, 0]
        ax.hist(scores_data, bins=50, range=(0, 1), alpha=0.7, 
               label=f'Data ({len(scores_data):,})', color='black', 
               histtype='step', linewidth=2)
        ax.hist(scores_sig, bins=50, range=(0, 1), weights=weights_sig*scale_sig,
               alpha=0.5, label='Signal (VBF H→μμ)', color='blue')
        ax.hist(scores_bkg, bins=50, range=(0, 1), weights=weights_bkg*scale_bkg,
               alpha=0.5, label='Background', color='red')
        ax.set_xlabel('BDT Score')
        ax.set_ylabel('Events')
        ax.set_title('BDT Score Distribution')
        ax.legend()
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Category definitions
        ax = axes[0, 1]
        for cat_name, (low, high, label) in categories.items():
            ax.axvspan(low, high, alpha=0.3, 
                      label=f'{label} [{low:.2f},{high:.2f})')
        ax.set_xlabel('BDT Score')
        ax.set_ylabel('Category')
        ax.set_title('Analysis Categories')
        ax.set_xlim(0, 1)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Data/MC comparison per category
        ax = axes[1, 0]
        cat_names = []
        data_yields = []
        mc_yields = []
        
        for cat_name, (low, high, label) in categories.items():
            mask_data = (scores_data >= low) & (scores_data < high)
            mask_sig = (scores_sig >= low) & (scores_sig < high)
            mask_bkg = (scores_bkg >= low) & (scores_bkg < high)
            
            cat_names.append(label)
            data_yields.append(np.sum(mask_data))
            mc_yields.append(np.sum(weights_sig[mask_sig])*scale_sig + 
                           np.sum(weights_bkg[mask_bkg])*scale_bkg)
        
        x = np.arange(len(cat_names))
        width = 0.35
        ax.bar(x - width/2, data_yields, width, label='Data', color='black', alpha=0.7)
        ax.bar(x + width/2, mc_yields, width, label='MC', color='red', alpha=0.7)
        ax.set_ylabel('Events')
        ax.set_title('Event Yields by Category')
        ax.set_xticks(x)
        ax.set_xticklabels(cat_names, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # Plot 4: S/B ratio vs threshold
        ax = axes[1, 1]
        thresholds = np.linspace(0, 1, 100)
        sb_ratios = []
        
        for thresh in thresholds:
            S = np.sum(weights_sig[scores_sig >= thresh]) * scale_sig
            B = np.sum(weights_bkg[scores_bkg >= thresh]) * scale_bkg
            if B > 0:
                sb_ratios.append(S / B)
            else:
                sb_ratios.append(0)
        
        ax.plot(thresholds, sb_ratios, 'b-', linewidth=2)
        ax.axvline(0.80, color='orange', linestyle='--', label='Loose cut')
        ax.axvline(0.95, color='red', linestyle='--', label='Tight cut')
        ax.set_xlabel('BDT Score Threshold')
        ax.set_ylabel('S/B Ratio')
        ax.set_title('Signal-to-Background Ratio')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # ========== PAGE 3-5: Mass Distributions per Category ==========
        mass_range = (110, 150)
        signal_region = (120, 130)
        
        fit_results = {}
        
        for cat_name, (bdt_low, bdt_high, cat_label) in categories.items():
            if cat_name == 'fail':
                continue  # Skip control region for signal extraction
            
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle(f'{cat_label} Category: BDT ∈ [{bdt_low:.2f}, {bdt_high:.2f})', 
                        fontsize=14, fontweight='bold')
            
            # Select events in this category
            mask_data = ((data_data['bdt_score'] >= bdt_low) & 
                        (data_data['bdt_score'] < bdt_high))
            mask_sig = ((data_sig['bdt_score'] >= bdt_low) & 
                       (data_sig['bdt_score'] < bdt_high))
            mask_bkg = ((data_bkg['bdt_score'] >= bdt_low) & 
                       (data_bkg['bdt_score'] < bdt_high))
            
            m_data_cat = data_data['m_mumu'][mask_data]
            m_sig_cat = data_sig['m_mumu'][mask_sig]
            m_bkg_cat = data_bkg['m_mumu'][mask_bkg]
            w_sig_cat = data_sig['weight'][mask_sig] * scale_sig
            w_bkg_cat = data_bkg['weight'][mask_bkg] * scale_bkg
            
            # Plot 1: Stacked mass distribution
            ax = axes[0, 0]
            bins = np.linspace(mass_range[0], mass_range[1], 41)
            
            # Stack backgrounds
            if bkg_components:
                bkg_hists = []
                bkg_labels = []
                for bkg_file in bkg_components:
                    try:
                        bkg_data = load_data(bkg_file, args.treename)
                        mask = ((bkg_data['bdt_score'] >= bdt_low) & 
                               (bkg_data['bdt_score'] < bdt_high))
                        if np.sum(mask) > 0:
                            bkg_name = os.path.basename(bkg_file).split('_')[2].split('.')[0]
                            hist, _ = np.histogram(bkg_data['m_mumu'][mask], 
                                                  bins=bins, 
                                                  weights=bkg_data['weight'][mask]*scale_bkg)
                            bkg_hists.append(hist)
                            bkg_labels.append(bkg_name)
                    except:
                        pass
                
                if bkg_hists:
                    ax.hist([bins[:-1] for _ in bkg_hists], bins=bins, 
                           weights=bkg_hists, stacked=True, 
                           label=bkg_labels, alpha=0.7)
            else:
                ax.hist(m_bkg_cat, bins=bins, weights=w_bkg_cat, 
                       alpha=0.5, label='Background', color='red')
            
            # Overlay signal
            ax.hist(m_sig_cat, bins=bins, weights=w_sig_cat, 
                   alpha=0.5, label='Signal (VBF H→μμ)', 
                   color='blue', histtype='step', linewidth=2)
            
            # Data points
            data_hist, _ = np.histogram(m_data_cat, bins=bins)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            ax.errorbar(bin_centers, data_hist, yerr=np.sqrt(data_hist), 
                       fmt='ko', label='Data', markersize=4)
            
            ax.axvline(125, color='black', linestyle='--', alpha=0.5, linewidth=1)
            ax.axvspan(signal_region[0], signal_region[1], alpha=0.1, color='green')
            ax.set_xlabel('m_μμ [GeV]')
            ax.set_ylabel('Events / 1 GeV')
            ax.set_title('Dimuon Mass Distribution')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            
            # Plot 2: Data with fit
            ax = axes[0, 1]
            
            # Perform fit
            print(f"\n  Fitting {cat_label} category...")
            fit_result = fit_mass_distribution(
                m_data_cat, np.ones(len(m_data_cat)),  # unit weights for data
                m_sig_cat, w_sig_cat,
                mass_range=mass_range,
                signal_region=signal_region
            )
            
            if fit_result:
                fit_results[cat_name] = fit_result
                
                # Plot data
                ax.errorbar(bin_centers, data_hist, yerr=np.sqrt(data_hist), 
                           fmt='ko', label='Data', markersize=4)
                
                # Plot fit
                m_fine = np.linspace(mass_range[0], mass_range[1], 200)
                y_fit = combined_model(m_fine, *fit_result['params'])
                y_sig = signal_model(m_fine, *fit_result['params'][:3])
                y_bkg = background_model(m_fine, *fit_result['params'][3:])
                
                ax.plot(m_fine, y_fit, 'r-', linewidth=2, label='Signal+Background Fit')
                ax.plot(m_fine, y_sig, 'b--', linewidth=2, label='Signal Component')
                ax.plot(m_fine, y_bkg, 'g--', linewidth=2, label='Background Component')
                
                ax.axvline(125, color='black', linestyle='--', alpha=0.5)
                ax.set_xlabel('m_μμ [GeV]')
                ax.set_ylabel('Events / 1 GeV')
                ax.set_title(f'Fit: S={fit_result["signal_yield"]:.1f}, ' + 
                           f'B={fit_result["background_yield"]:.1f}')
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                
                # Add fit parameters as text
                fit_text = (f'μ = {fit_result["params"][1]:.2f} ± {fit_result["errors"][1]:.2f} GeV\n'
                          f'σ = {fit_result["params"][2]:.2f} ± {fit_result["errors"][2]:.2f} GeV\n'
                          f'χ²/ndf = {fit_result["chi2"]:.1f}/{fit_result["ndf"]}\n'
                          f'Significance = {fit_result["significance"]:.2f}σ')
                ax.text(0.05, 0.95, fit_text, transform=ax.transAxes,
                       fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Plot 3: Pull distribution
            ax = axes[1, 0]
            if fit_result:
                y_fit_bins = combined_model(bin_centers, *fit_result['params'])
                pulls = (data_hist - y_fit_bins) / np.sqrt(data_hist + 1e-10)
                
                ax.axhline(0, color='black', linestyle='-', linewidth=1)
                ax.axhline(2, color='red', linestyle='--', linewidth=1, alpha=0.5)
                ax.axhline(-2, color='red', linestyle='--', linewidth=1, alpha=0.5)
                ax.errorbar(bin_centers, pulls, yerr=1, fmt='ko', markersize=4)
                ax.set_xlabel('m_μμ [GeV]')
                ax.set_ylabel('Pull (σ)')
                ax.set_title('Fit Residuals')
                ax.grid(True, alpha=0.3)
                ax.set_ylim(-4, 4)
            
            # Plot 4: Signal region zoom
            ax = axes[1, 1]
            sr_bins = np.linspace(signal_region[0], signal_region[1], 21)
            sr_mask_data = ((m_data_cat >= signal_region[0]) & 
                           (m_data_cat < signal_region[1]))
            sr_mask_sig = ((m_sig_cat >= signal_region[0]) & 
                          (m_sig_cat < signal_region[1]))
            sr_mask_bkg = ((m_bkg_cat >= signal_region[0]) & 
                          (m_bkg_cat < signal_region[1]))
            
            ax.hist(m_bkg_cat[sr_mask_bkg], bins=sr_bins, weights=w_bkg_cat[sr_mask_bkg],
                   alpha=0.5, label='Background', color='red')
            ax.hist(m_sig_cat[sr_mask_sig], bins=sr_bins, weights=w_sig_cat[sr_mask_sig],
                   alpha=0.5, label='Signal', color='blue', histtype='step', linewidth=2)
            
            sr_data_hist, _ = np.histogram(m_data_cat[sr_mask_data], bins=sr_bins)
            sr_centers = (sr_bins[:-1] + sr_bins[1:]) / 2
            ax.errorbar(sr_centers, sr_data_hist, yerr=np.sqrt(sr_data_hist),
                       fmt='ko', label='Data', markersize=4)
            
            ax.axvline(125, color='black', linestyle='--', alpha=0.5)
            ax.set_xlabel('m_μμ [GeV]')
            ax.set_ylabel('Events / 0.5 GeV')
            ax.set_title('Signal Region [120-130 GeV]')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
        
        # ========== PAGE 6: Systematic Uncertainties ==========
        print("\n  Calculating systematic uncertainties...")
        systematics = calculate_systematics(data_sig, data_bkg, categories)
        
        fig = plt.figure(figsize=(11, 8.5))
        fig.text(0.5, 0.95, 'Systematic Uncertainties', 
                ha='center', fontsize=18, fontweight='bold')
        
        # Table of systematics
        sys_text = """
Source                          Uncertainty
{'='*60}
"""
        for name, value in systematics.items():
            if name != 'total':
                sys_text += f"{name:30s}    {value*100:5.1f}%\n"
        sys_text += f"{'='*60}\n"
        sys_text += f"{'Total':30s}    {systematics['total']*100:5.1f}%\n"
        
        fig.text(0.2, 0.7, sys_text, fontsize=11, family='monospace',
                verticalalignment='top')
        
        # Bar chart
        ax = fig.add_subplot(212)
        sys_names = [k for k in systematics.keys() if k != 'total']
        sys_values = [systematics[k]*100 for k in sys_names]
        
        colors = ['lightcoral' if k == 'total' else 'skyblue' for k in sys_names]
        ax.barh(sys_names, sys_values, color=colors, alpha=0.7)
        ax.axvline(systematics['total']*100, color='red', linestyle='--', 
                  linewidth=2, label=f'Total: {systematics["total"]*100:.1f}%')
        ax.set_xlabel('Uncertainty (%)')
        ax.set_title('Systematic Uncertainty Breakdown')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # ========== PAGE 7: Combined Results & Limits ==========
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Combined Analysis Results', fontsize=18, fontweight='bold')
        
        # Combine categories
        total_obs = 0
        total_sig = 0
        total_bkg = 0
        
        results_text = """
Category-by-Category Results
{'='*70}
Category          Observed    Expected Sig    Expected Bkg    S/√B
{'='*70}
"""
        
        for cat_name in ['tight', 'loose']:
            if cat_name not in fit_results:
                continue
            
            bdt_low, bdt_high, cat_label = categories[cat_name]
            
            mask_data = ((data_data['bdt_score'] >= bdt_low) & 
                        (data_data['bdt_score'] < bdt_high) &
                        (data_data['m_mumu'] >= signal_region[0]) &
                        (data_data['m_mumu'] < signal_region[1]))
            
            mask_sig = ((data_sig['bdt_score'] >= bdt_low) & 
                       (data_sig['bdt_score'] < bdt_high) &
                       (data_sig['m_mumu'] >= signal_region[0]) &
                       (data_sig['m_mumu'] < signal_region[1]))
            
            mask_bkg = ((data_bkg['bdt_score'] >= bdt_low) & 
                       (data_bkg['bdt_score'] < bdt_high) &
                       (data_bkg['m_mumu'] >= signal_region[0]) &
                       (data_bkg['m_mumu'] < signal_region[1]))
            
            obs = np.sum(mask_data)
            sig = np.sum(data_sig['weight'][mask_sig]) * scale_sig
            bkg = np.sum(data_bkg['weight'][mask_bkg]) * scale_bkg
            
            significance = sig / np.sqrt(bkg) if bkg > 0 else 0
            
            results_text += f"{cat_label:15s}   {obs:8.0f}    {sig:12.2f}    {bkg:12.2f}    {significance:6.2f}\n"
            
            total_obs += obs
            total_sig += sig
            total_bkg += bkg
        
        results_text += f"{'='*70}\n"
        total_significance = total_sig / np.sqrt(total_bkg) if total_bkg > 0 else 0
        results_text += f"{'Combined':15s}   {total_obs:8.0f}    {total_sig:12.2f}    {total_bkg:12.2f}    {total_significance:6.2f}\n"
        
        fig.text(0.1, 0.85, results_text, fontsize=10, family='monospace',
                verticalalignment='top')
        
        # Calculate limits
        print("\n  Calculating limits...")
        limits = calculate_limits(total_obs, total_sig, total_bkg, systematics)
        
        if limits:
            limit_text = f"""
95% CL Upper Limits on Signal Strength μ
{'='*50}
Observed limit:           {limits['observed']:.2f}
Expected limit:           {limits['expected']:.2f}
Expected +1σ:             {limits['expected_p1sigma']:.2f}
Expected -1σ:             {limits['expected_m1sigma']:.2f}

Interpretation:
  μ = 1 corresponds to SM Higgs production
  σ(VBF) × BR(H→μμ) < {limits['observed'] * args.signal_xsec:.3f} pb
"""
            fig.text(0.1, 0.45, limit_text, fontsize=11, family='monospace',
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        
        plt.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # ========== PAGE 8: Validation Plots ==========
        fig, axes = plt.subplots(2, 3, figsize=(14, 10))
        fig.suptitle('Data/MC Validation in Control Regions', fontsize=14, fontweight='bold')
        
        # Use fail category as control region
        fail_low, fail_high, _ = categories['fail']
        mask_data_cr = ((data_data['bdt_score'] >= fail_low) & 
                       (data_data['bdt_score'] < fail_high))
        mask_sig_cr = ((data_sig['bdt_score'] >= fail_low) & 
                      (data_sig['bdt_score'] < fail_high))
        mask_bkg_cr = ((data_bkg['bdt_score'] >= fail_low) & 
                      (data_bkg['bdt_score'] < fail_high))
        
        # Plot various distributions
        variables = [
            ('m_mumu', 'm_μμ [GeV]', (110, 150), 40),
            ('pt_mumu', 'pT_μμ [GeV]', (0, 300), 30),
            ('njets', 'Number of jets', (0, 10), 10)
        ]
        
        for idx, (var, label, range_val, nbins) in enumerate(variables):
            if var not in data_data:
                continue
            
            ax = axes[idx // 3, idx % 3]
            
            bins = np.linspace(range_val[0], range_val[1], nbins+1)
            
            # MC
            if var in data_bkg:
                ax.hist(data_bkg[var][mask_bkg_cr], bins=bins, 
                       weights=data_bkg['weight'][mask_bkg_cr]*scale_bkg,
                       alpha=0.5, label='Background', color='red')
            if var in data_sig:
                ax.hist(data_sig[var][mask_sig_cr], bins=bins, 
                       weights=data_sig['weight'][mask_sig_cr]*scale_sig,
                       alpha=0.5, label='Signal', color='blue', 
                       histtype='step', linewidth=2)
            
            # Data
            data_hist, _ = np.histogram(data_data[var][mask_data_cr], bins=bins)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            ax.errorbar(bin_centers, data_hist, yerr=np.sqrt(data_hist),
                       fmt='ko', label='Data', markersize=4)
            
            ax.set_xlabel(label)
            ax.set_ylabel('Events')
            ax.set_title(f'{label} (Control Region)')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        print(f"\n✓ Created comprehensive PDF: {output_pdf}")

def main():
    args = parse_args()
    
    print("="*80)
    print("COMPREHENSIVE H→μμ VBF ANALYSIS")
    print("="*80)
    
    # Load data
    print("\nLoading data files...")
    print(f"  Data: {args.data}")
    data_data = load_data(args.data, args.treename)
    
    print(f"  Signal MC: {args.signal}")
    data_sig = load_data(args.signal, args.treename)
    
    print(f"  Background MC: {args.background}")
    data_bkg = load_data(args.background, args.treename)
    
    # Load individual background components if provided
    bkg_components = []
    if args.bkg_components:
        print(f"  Background components:")
        for bkg_file in args.bkg_components:
            if os.path.exists(bkg_file):
                print(f"    - {os.path.basename(bkg_file)}")
                bkg_components.append(bkg_file)
    
    # Define categories
    categories = define_categories()
    
    print(f"\nAnalysis configuration:")
    print(f"  Luminosity: {args.lumi:.2f} fb⁻¹")
    print(f"  Signal cross-section: {args.signal_xsec:.4f} pb")
    print(f"  Categories: {len(categories)}")
    
    # Create comprehensive PDF
    print(f"\nCreating comprehensive analysis PDF...")
    create_comprehensive_pdf(
        data_data, data_sig, data_bkg, bkg_components,
        categories, args, args.output
    )
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\n✓ Output saved to: {args.output}")
    print("\nThe PDF contains:")
    print("  • Page 1: Analysis summary")
    print("  • Page 2: BDT score distributions and S/B ratios")
    print("  • Pages 3-5: Mass distributions and fits per category")
    print("  • Page 6: Systematic uncertainties")
    print("  • Page 7: Combined results and 95% CL limits")
    print("  • Page 8: Data/MC validation plots")
    print("="*80 + "\n")

if __name__ == "__main__":
    main()