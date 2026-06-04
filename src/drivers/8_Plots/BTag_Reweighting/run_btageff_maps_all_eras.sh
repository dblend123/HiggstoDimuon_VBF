#!/bin/bash

# Example:
#   ./run_btageff_maps_all_eras.sh KIT noScaleSmear_noBtagSF_noZPtRW

CORR_MODE="${1:-KIT}"
CONFIG_DIR="${2:-noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW}"

ERAS=(12022 22022 12023 22023 2024 2025)

# ---- SET YOUR MACRO DIRECTORY HERE ----
MACRO_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/8_Plots/BTag_Reweighting"

BTAG_MACRO="${MACRO_DIR}/BTagEff_Maps.C"

DEST_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/JME"

for ERA in "${ERAS[@]}"; do

    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/${CORR_MODE}/${CONFIG_DIR}"

    echo "=================================================="
    echo " Running BTag Efficiency Map Workflow for ERA: $ERA"
    echo " Path: $BASE_PATH"
    echo "=================================================="

    if [ ! -d "$BASE_PATH" ]; then
        echo "Directory does not exist, skipping: $BASE_PATH"
        echo
        continue
    fi

    echo "[STEP 1] Build BTag efficiency maps"

    root -l -b <<EOF
.L ${BTAG_MACRO}
BTagEff_Maps("${BASE_PATH}");
.q
EOF

    SRC_FILE="/eos/user/d/dblend/StackedPlots/${ERA}/BTagEfficiency/${CORR_MODE}/${CONFIG_DIR}/BTagEffMaps_${ERA}_${CORR_MODE}_${CONFIG_DIR}.root"

    if [ -f "$SRC_FILE" ]; then
        echo "[COPY] $SRC_FILE → $DEST_DIR"
        cp "$SRC_FILE" "${DEST_DIR}/BTagEffMaps_${ERA}.root"
    else
        echo "[WARNING] Missing output file!"
    fi

    echo "Finished BTag efficiency maps for ERA: $ERA"
    echo "=================================================="
    echo

done
