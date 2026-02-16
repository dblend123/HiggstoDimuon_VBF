#!/usr/bin/env python3
"""
Create combined mass plot for all eras using ONLY events at AMS-optimized BDT cuts
Each era uses its own optimal cut from the BDT analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ============================================
# DATA AT AMS-OPTIMIZED BDT CUTS
# ============================================

data_masses = {
    '12022 (BDT≥0.958)': [
        111.40, 113.13, 114.51, 116.32, 119.71, 121.35, 121.93, 122.88,
        132.67, 133.55, 133.92, 137.01, 142.71, 142.98, 143.07, 144.02,
        144.54, 144.94, 145.57, 146.36, 149.08
    ],
    
    '22022 (BDT≥0.976)': [
        112.00, 112.57, 113.14, 113.26, 114.16, 114.26, 114.94, 115.42,
        116.29, 116.54, 117.36, 117.94, 118.21, 118.59, 120.60, 121.89,
        122.90, 128.89, 130.57, 132.42, 132.66, 141.42, 142.79, 148.48
    ],
    
    '12023 (BDT≥0.990)': [
        112.90, 112.90, 115.15, 115.15, 116.84, 116.84, 118.90, 118.90,
        122.16, 122.16, 123.80, 123.80, 125.24, 125.24, 131.22, 131.22,
        134.03, 134.03, 144.22, 144.22, 146.81, 146.81, 147.20, 147.20
    ],
    
    '22023 (BDT≥0.930)': [
        110.34, 110.65, 110.68, 110.99, 111.11, 111.54, 111.63, 111.71,
        112.06, 112.17, 112.40, 112.52, 112.53, 112.83, 112.91, 113.07,
        113.09, 113.18, 113.43, 114.08, 114.19, 114.31, 114.36, 114.38,
        114.60, 114.61, 115.02, 115.39, 115.49, 115.64, 115.88, 116.02,
        116.12, 116.18, 116.23, 116.32, 116.38, 116.76, 117.00, 117.33,
        117.37, 117.38, 117.42, 117.58, 117.68, 117.97, 118.73, 118.73,
        119.04, 119.39, 119.94, 120.01, 120.07, 120.54, 120.64, 120.75,
        120.86, 120.95, 121.28, 121.37, 121.41, 121.62, 121.78, 121.89,
        122.08, 122.16, 122.25, 122.63, 122.66, 122.84, 122.89, 123.37,
        123.40, 123.60, 123.90, 124.40, 124.74, 125.16, 125.19, 125.51,
        125.88, 125.96, 125.98, 126.00, 126.09, 126.23, 126.47, 126.48,
        126.55, 126.72, 126.85, 127.18, 127.68, 128.11, 128.58, 128.83,
        129.14, 129.56, 129.66, 129.72, 130.20, 131.40, 132.95, 133.29,
        133.39, 133.60, 133.71, 133.91, 134.02, 134.82, 135.12, 135.29,
        135.33, 135.58, 135.82, 135.95, 135.96, 136.22, 136.54, 136.99,
        137.28, 137.44, 137.47, 137.53, 138.53, 139.08, 139.13, 139.42,
        139.91, 142.71, 143.18, 143.22, 143.96, 144.94, 145.06, 145.10,
        145.21, 146.40, 147.24, 147.99, 148.40, 149.31, 149.31, 149.33,
        149.65
    ],
    
    '2024 (BDT≥0.986)': [
        113.52, 117.65, 123.44, 125.66
    ]
}

# Combine all masses
all_masses = []
total_in_higgs = 0

print("="*60)
print("EVENTS AT AMS-OPTIMIZED BDT CUTS")
print("="*60)

for era, masses in data_masses.items():
    all_masses.extend(masses)
    in_higgs = sum(120 <= m <= 130 for m in masses)
    total_in_higgs += in_higgs
    print(f"\n{era}:")
    print(f"  Total events: {len(masses)}")
    print(f"  Events in Higgs window (120-130 GeV): {in_higgs}")
    if in_higgs > 0:
        higgs_events = [m for m in masses if 120 <= m <= 130]
        print(f"  Masses: {', '.join([f'{m:.2f}' for m in sorted(higgs_gs_events])}")

print("\n" + "="*60)
print(f"TOTAL EVENTS: {len(all_masses)}")
print(f"TOTAL in Higgs window (120-130 GeV): {total_in_higgs}")
print("="*60)

# Create the plot
fig, axes = plt.subplots(2, 1, figsize=(14, 10), 
                          gridspec_kw={'height_ratios': [3, 1]})

# Define mass bins (2 GeV bins from 110-150)
bins = np.linspace(110, 150, 21)
bin_centers = (bins[:-1] + bins[1:]) / 2

# Plot 1: Mass distribution
axes[0].hist(all_masses, bins=bins, histtype='stepfilled', 
             alpha=0.7, color='steelblue', edgecolor='black',
             label=f'Data (all eras, {len(all_masses)} events)')

# Add vertical line at Higgs mass
axes[0].axvline(x=125, color='red', linestyle='--', linewidth=2, 
                alpha=0.8, label='Higgs mass (125 GeV)')

# Shade Higgs window
axes[0].axvspan(120, 130, alpha=0.2, color='yellow', 
                label='Signal region (120-130 GeV)')

# Add text boxes for each era
text_y = 0.95
for era, masses in data_masses.items():
    n_events = len(masses)
    axes[0].text(0.02, text_y, f'{era}: {n_events} events', 
                 transform=axes[0].transAxes, fontsize=9,
                 verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    text_y -= 0.05

axes[0].set_xlabel('Dimuon Mass [GeV]', fontsize=12)
axes[0].set_ylabel('Events / 2 GeV', fontsize=12)
axes[0].set_title('Run3 Data at AMS-Optimized BDT Cuts', fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3)
axes[0].legend(loc='upper right')
axes[0].set_xlim(110, 150)
axes[0].set_ylim(0, max(30, np.histogram(all_masses, bins=bins)[0].max() + 2))

# Plot 2: Cumulative distribution
sorted_masses = sorted(all_masses)
cumulative = np.arange(1, len(sorted_masses) + 1) / len(sorted_masses) * 100

axes[1].plot(sorted_masses, cumulative, 'b-', linewidth=2, label='Cumulative distribution')
axes[1].axvline(x=125, color='red', linestyle='--', alpha=0.8)
axes[1].axvspan(120, 130, alpha=0.2, color='yellow')

# Find percentage in Higgs window
percent_in_higgs = total_in_higgs / len(all_masses) * 100
axes[1].axhline(y=percent_in_higgs, color='green', linestyle=':', 
                label=f'{percent_in_higgs:.1f}% in Higgs window')

axes[1].set_xlabel('Dimuon Mass [GeV]', fontsize=12)
axes[1].set_ylabel('Cumulative Percentage [%]', fontsize=12)
axes[1].set_title('Cumulative Distribution', fontsize=12)
axes[1].grid(True, alpha=0.3)
axes[1].legend(loc='lower right')
axes[1].set_xlim(110, 150)
axes[1].set_ylim(0, 100)

plt.tight_layout()

# Save the plot
plt.savefig('run3_data_ams_cuts_only.pdf', bbox_inches='tight')
plt.savefig('run3_data_ams_cuts_only.png', dpi=150, bbox_inches='tight')

print("\n✓ Plot saved as 'run3_data_ams_cuts_only.pdf'")
print("✓ Also saved as 'run3_data_ams_cuts_only.png'")

# Optional: Create a zoomed-in plot of the Higgs region
fig_zoom, ax_zoom = plt.subplots(figsize=(12, 6))

zoom_masses = [m for m in all_masses if 115 <= m <= 135]
zoom_bins = np.linspace(115, 135, 21)

ax_zoom.hist(zoom_masses, bins=zoom_bins, histtype='stepfilled',
             alpha=0.7, color='steelblue', edgecolor='black')
ax_zoom.axvline(x=125, color='red', linestyle='--', linewidth=2, alpha=0.8)
ax_zoom.axvspan(120, 130, alpha=0.2, color='yellow')

# Count events in each 1 GeV bin for labeling
for i in range(120, 131):
    count = sum(i <= m < i+1 for m in all_masses)
    if count > 0:
        ax_zoom.text(i+0.5, count+0.2, str(count), 
                    ha='center', va='bottom', fontweight='bold')

ax_zoom.set_xlabel('Dimuon Mass [GeV]', fontsize=12)
ax_zoom.set_ylabel('Events / 1 GeV', fontsize=12)
ax_zoom.set_title('Zoom: Higgs Region (115-135 GeV)', fontsize=14, fontweight='bold')
ax_zoom.grid(True, alpha=0.3)
ax_zoom.set_xlim(115, 135)

plt.tight_layout()
plt.savefig('run3_data_higgs_region_zoom.pdf', bbox_inches='tight')
print("✓ Zoomed plot saved as 'run3_data_higgs_region_zoom.pdf'")

plt.show()