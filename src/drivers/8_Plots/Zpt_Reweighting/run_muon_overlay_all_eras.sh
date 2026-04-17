#!/bin/bash

# Example: ./run_muon_overlay_all_eras.sh KIT noScaleSmear_noBtagSF_noZPtRW

CORR_MODE="${1:-KIT}"
CONFIG_DIR="${2:-noScaleSmear_noBtagSF_noZPtRW}"

#ERAS=(12022 22022 12023 22023 2024)
ERAS=(2024)

# ---- SET YOUR MACRO DIRECTORY HERE ----
MACRO_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/8_Plots/Zpt_Reweighting"

MUON_MACRO="${MACRO_DIR}/Muon_Overlay_DiMuonpt.C"

for ERA in "${ERAS[@]}"; do

    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/${CORR_MODE}/${CONFIG_DIR}"

    echo "=================================================="
    echo " Running Muon DiMuonpt Overlay for ERA: $ERA"
    echo " Path: $BASE_PATH"
    echo "=================================================="

    if [ ! -d "$BASE_PATH" ]; then
        echo "Directory does not exist, skipping: $BASE_PATH"
        echo
        continue
    fi

    echo "[STEP 3] Muon DiMuonpt Overlay"

    root -l -b <<EOF
.L ${MUON_MACRO}
Muon_Overlay_DiMuonpt("${BASE_PATH}", "${ERA}");
.q
EOF

    echo "Finished Muon overlay for ERA: $ERA"
    echo "=================================================="
    echo

done
