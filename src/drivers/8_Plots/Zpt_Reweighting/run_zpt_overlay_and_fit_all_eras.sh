#!/bin/bash

# Example:
# ./run_zpt_overlay_and_fit_all_eras.sh KIT noScaleSmear_noBtagSF_noZPtRW 2

CORR_MODE="${1:-KIT}"
CONFIG_DIR="${2:-noScaleSmear_noBtagSF_noZPtRW}"
REBIN="${3:-1}"

ERAS=(12022 22022 12023 22023)
# ERAS=(12022 22022 12023 22023 2024)

MACRO_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/8_Plots/Zpt_Reweighting"

OVERLAY_MACRO="${MACRO_DIR}/ZPt_Overlay_Plots.C"
FIT_MACRO="${MACRO_DIR}/FitZPtSF_Piecewise_StackedPlots.C"

DEST_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/7_Corrections/ZptReweight"

mkdir -p "$DEST_DIR"

for ERA in "${ERAS[@]}"; do

    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/${CORR_MODE}/${CONFIG_DIR}"

    echo "=================================================="
    echo " Running ZPt Overlay + Fit Workflow for ERA: $ERA"
    echo " Path: $BASE_PATH"
    echo "=================================================="

    if [ ! -d "$BASE_PATH" ]; then
        echo "[WARNING] Directory does not exist, skipping: $BASE_PATH"
        echo
        continue
    fi

    # ------------------- ZPt Overlay -------------------
    echo "[STEP 1] ZPt Overlay"

    root -l -b <<EOF
.L ${OVERLAY_MACRO}
ZPt_Overlay_Plots("${BASE_PATH}", "${ERA}", ${REBIN});
.q
EOF

    echo "[OK] Finished ZPt overlay for ERA: $ERA"

    # ------------------- Piecewise Fit -------------------
    SF_FILE="/eos/user/d/dblend/StackedPlots/${ERA}/ZPt_Corrections/${CORR_MODE}/${CONFIG_DIR}/ZPtSF_${ERA}_${CORR_MODE}_${CONFIG_DIR}.root"

    if [ -f "$SF_FILE" ]; then
        echo "[STEP 2] Piecewise Fit"

        root -l -b <<EOF
.L ${FIT_MACRO}
FitZPtSF_Piecewise_StackedPlots("${SF_FILE}");
.q
EOF

        echo "[OK] Finished fits for ERA: $ERA"
    else
        echo "[WARNING] SF file not found, skipping fit: $SF_FILE"
        echo
        continue
    fi

    # ------------------- Copy to AFS with short name -------------------
    FIT_DIR="/eos/user/d/dblend/StackedPlots/${ERA}/ZPt_Corrections/${CORR_MODE}/${CONFIG_DIR}/Fits"

    SRC_FILE="${FIT_DIR}/ZPtSF_Fits_${ERA}_${CORR_MODE}_${CONFIG_DIR}.root"
    DEST_FILE="${DEST_DIR}/ZPtSF_Fits_${ERA}.root"

    if [ -f "$SRC_FILE" ]; then
        echo "[STEP 3] Copying ZPt SF file to AFS"
        cp -f "$SRC_FILE" "$DEST_FILE"
        echo "[OK] Copied:"
        echo "     From: $SRC_FILE"
        echo "     To:   $DEST_FILE"
    else
        echo "[WARNING] Fit output not found, skipping copy: $SRC_FILE"
        echo "[DEBUG] Available files in Fits directory:"
        ls -lh "$FIT_DIR"
    fi

    echo "=================================================="
    echo

done
