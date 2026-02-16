#!/bin/bash

# Script to run BDT analysis on all Run3 eras with correct luminosities
# Eras and their luminosities (fb^-1):
# 12022 (2022C-D): 7.9805
# 22022 (2022E-G): 26.672
# 12023 (2023C): 17.794
# 22023 (2023D): 9.4510
# 2024: 104.8

# Define the eras and their corresponding luminosities
declare -A LUMI_MAP=(
    ["12022"]="7.9805"
    ["22022"]="26.672"
    ["12023"]="17.794"
    ["22023"]="9.4510"
)

# Loop over each era separately
for ERA in 12022 22022 12023 22023; do
    echo "========================================="
    echo "Running BDT analysis for era: ${ERA}"
    echo "Luminosity: ${LUMI_MAP[$ERA]} fb^-1"
    echo "========================================="
    
    # Get luminosity for this era
    LUMI=${LUMI_MAP[$ERA]}
    
    # Run the analysis for this era
    python3 scripts/bdt_analysis.py \
        --outputs-dir "Outputs" \
        --eras ${ERA} \
        --treename "mva" \
        --bdt-branch "bdt_score" \
        --weight-branch "weight" \
        --mass-branch "m_mumu" \
        --mass-min 110 \
        --mass-max 150 \
        --alpha 1e-3 \
        --dthr-pct 0.20 \
        --abs-weights \
        --output "bdt_analysis_${ERA}.pdf" \
        --retrain-ks \
        --ks-dir DataPreparation_Files \
        --data-era ${ERA}:Outputs/${ERA}/${ERA}_data_scored.root \
        --mc-era ${ERA}:DY:Outputs/${ERA}/${ERA}_scored_DY.root \
        --mc-era ${ERA}:Top:Outputs/${ERA}/${ERA}_scored_TT.root \
        --mc-era ${ERA}:EWK:Outputs/${ERA}/${ERA}_scored_EWK.root \
        --mc-era ${ERA}:Diboson:Outputs/${ERA}/${ERA}_scored_VV.root \
        --mc-era ${ERA}:Diboson:Outputs/${ERA}/${ERA}_scored_VVV.root \
        --mc-era ${ERA}:VBF:Outputs/${ERA}/${ERA}_scored_signal.root \
        --datamc-nbins 40 \
        --datamc-logy \
        --lumi ${LUMI}
    
    # Check if the command succeeded
    if [ $? -eq 0 ]; then
        echo "✓ Analysis completed successfully for era ${ERA}"
    else
        echo "✗ ERROR: Analysis failed for era ${ERA}"
        exit 1
    fi
    
    echo ""
done

echo "========================================="
echo "All individual era analyses completed!"
echo "========================================="

# Optional: Create a summary of high-BDT events across all eras
echo ""
echo "=== HIGH-BDT EVENT SUMMARY ==="
echo ""

for ERA in 12022 22022 12023 22023; do
    echo "Era ${ERA} (Lumi = ${LUMI_MAP[$ERA]} fb^-1):"
    echo "  Output PDF: bdt_analysis_${ERA}.pdf"
    echo "  Check the log above for high-BDT event masses"
    echo ""
done

