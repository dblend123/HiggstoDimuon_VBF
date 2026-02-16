#!/bin/bash

# Script to run COMBINED BDT analysis on all Run3 eras
# This will produce a single PDF with all eras combined

# Combined luminosity for all eras (update these values!)
TOTAL_LUMI=172.0  # 41.5 + 41.5 + 44.3 + 44.3 = 171.6 fb^-1

echo "========================================="
echo "Running COMBINED BDT analysis for all eras"
echo "Total luminosity: ${TOTAL_LUMI} fb^-1"
echo "========================================="

# Build the command with all data and MC files
CMD="python3 scripts/bdt_analysis.py \
    --outputs-dir Outputs \
    --eras 12022 22022 12023 22023 \
    --treename mva \
    --bdt-branch bdt_score \
    --weight-branch weight \
    --mass-branch m_mumu \
    --mass-min 110 \
    --mass-max 150 \
    --alpha 1e-3 \
    --dthr-pct 0.20 \
    --abs-weights \
    --output bdt_analysis_all_eras_combined.pdf \
    --retrain-ks \
    --ks-dir DataPreparation_Files \
    --datamc-nbins 40 \
    --datamc-logy \
    --lumi ${TOTAL_LUMI} \
    --model-name myBDT"

# Add data files for each era
for ERA in 12022 22022 12023 22023; do
    CMD="${CMD} --data-era ${ERA}:Outputs/${ERA}/${ERA}_data_scored.root"
done

# Add MC files for each era and process
for ERA in 12022 22022 12023 22023; do
    for PROC in DY TT EWK VV VVV; do
        CMD="${CMD} --mc-era ${ERA}:${PROC}:Outputs/${ERA}/${ERA}_scored_${PROC}.root"
    done
    CMD="${CMD} --mc-era ${ERA}:VBF:Outputs/${ERA}/${ERA}_scored_signal.root"
done

# Add scan files (using wildcards - these will be expanded by the script)
CMD="${CMD} --scan-sig DataPreparation_Files/*/*_signal.root"
CMD="${CMD} --scan-bkg DataPreparation_Files/*/*_background.root"
CMD="${CMD} --scan-cuts 0.00:1.00:401"

# Execute the command
echo "Running command:"
echo "${CMD}"
echo ""
eval ${CMD}

if [ $? -eq 0 ]; then
    echo "✓ Combined analysis completed successfully!"
    echo "Output: bdt_analysis_all_eras_combined.pdf"
else
    echo "✗ ERROR: Combined analysis failed"
    exit 1
fi