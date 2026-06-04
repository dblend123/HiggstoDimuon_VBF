#!/bin/bash

# Example:
# ./run_zmass_resolution_all_eras.sh KIT noJERC_noJetVeto_noBtagSF_noZPtRW sigma68

CORR_MODE="${1:-KIT}"
CONFIG_DIR="${2:-noScaleSmear_noJERC_noJER_BTagRW_ZPtRW}"
METRIC="${3:-sigma68}"   # sigma68, sigmaCB, sigmaFwhm

ERAS=(12022 22022 12023 22023 2024 2025)

# ---- SET YOUR MACRO DIRECTORY ----
MACRO_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/8_Plots/Z_Boson"

FIT_MACRO="${MACRO_DIR}/MassResolution_ZCR_AllRegions.C"
SUMMARY_MACRO="${MACRO_DIR}/ZPeakAndResolutionSummaryPlot.C"

for ERA in "${ERAS[@]}"; do

    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/${CORR_MODE}/${CONFIG_DIR}"

    echo "=================================================="
    echo " Running Z Mass Resolution Workflow for ERA: $ERA"
    echo " Path: $BASE_PATH"
    echo "=================================================="

    if [ ! -d "$BASE_PATH" ]; then
        echo "Directory does not exist, skipping: $BASE_PATH"
        echo
        continue
    fi

    # --------------------------------------------------
    # STEP 1: Fit all regions (creates summary txt file)
    # --------------------------------------------------
    echo "[STEP 1] Running ZCR fits"

    root -l -b <<EOF
.L ${FIT_MACRO}
MassResolution_ZCR_AllRegions("${BASE_PATH}", "${ERA}");
.q
EOF

    echo "Finished fits for ERA: $ERA"

    # --------------------------------------------------
    # STEP 2: Build summary plot (peak + resolution)
    # --------------------------------------------------

    SUMMARY_FILE="/eos/user/d/dblend/StackedPlots/${ERA}/Mass_Resolutions/ZCR/${CORR_MODE}/${CONFIG_DIR}/ZPeakSummary_${ERA}.txt"

    if [ -f "$SUMMARY_FILE" ]; then
        echo "[STEP 2] Making combined summary plot"

        root -l -b <<EOF
.L ${SUMMARY_MACRO}
ZPeakAndResolutionSummaryPlot("${SUMMARY_FILE}", "${ERA}", "${METRIC}", true);
.q
EOF

        echo "Finished summary plot for ERA: $ERA"
    else
        echo "Summary file not found, skipping plot:"
        echo "$SUMMARY_FILE"
    fi

    echo "=================================================="
    echo

done
