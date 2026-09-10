#!/bin/bash

# Example:
# ./run_muon_and_jet_overlays_all_eras.sh KIT noScaleSmear_noBtagSF_ZPtRW

CORR_MODE="${1:-KIT}"
#CONFIG_DIR="${2:-noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW}"
#CONFIG_DIR="${2:-noScaleSmear_noJERC_noJER_BTagRW_ZPtRW}"
CONFIG_DIR="${2:-noScaleSmear_BTagRW_ZPtRW}"

#ERAS=(12022 22022 12023 22023 2024 2025)
ERAS=(2024)

# ---- SET YOUR MACRO DIRECTORY HERE ----
MACRO_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/8_Plots/Overlay_AllPlots"

MUON_MACRO="${MACRO_DIR}/Muon_Overlay_AllPlots.C"
JET_MACRO="${MACRO_DIR}/Jet_Overlay_AllPlots.C"

if [ ! -f "$MUON_MACRO" ]; then
    echo "[ERROR] Muon macro not found:"
    echo "  $MUON_MACRO"
    exit 1
fi

if [ ! -f "$JET_MACRO" ]; then
    echo "[ERROR] Jet macro not found:"
    echo "  $JET_MACRO"
    exit 1
fi

for ERA in "${ERAS[@]}"; do

    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/${CORR_MODE}/${CONFIG_DIR}"

    echo "=================================================="
    echo " Running Muon + Jet Overlay Workflow for ERA: $ERA"
    echo " Path: $BASE_PATH"
    echo "=================================================="

    if [ ! -d "$BASE_PATH" ]; then
        echo "[WARN] Directory does not exist, skipping:"
        echo "  $BASE_PATH"
        echo
        continue
    fi

    # ------------------- Muon Overlay -------------------
    echo "[STEP 1] Muon overlay plots"

    root -l -b <<EOF
.L ${MUON_MACRO}
Muon_Overlay_AllPlots("${BASE_PATH}", "${ERA}");
.q
EOF

    echo "Finished muon overlays for ERA: $ERA"

    # ------------------- Jet Overlay -------------------
    echo "[STEP 2] Jet overlay plots"

    root -l -b <<EOF
.L ${JET_MACRO}
Jet_Overlay_AllPlots("${BASE_PATH}", "${ERA}");
.q
EOF

    echo "Finished jet overlays for ERA: $ERA"

    echo "=================================================="
    echo

done
