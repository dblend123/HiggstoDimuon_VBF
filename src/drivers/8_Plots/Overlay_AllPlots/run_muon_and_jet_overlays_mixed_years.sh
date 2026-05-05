#!/bin/bash

# Example:
# ./run_muon_and_jet_overlays_mixed_years.sh KIT noScaleSmear_noBtagSF_ZPtRW 2025 2024
#
# Args:
#   $1 = correction mode
#   $2 = config dir
#   $3 = DATA era
#   $4 = MC era

CORR_MODE="${1:-KIT}"
CONFIG_DIR="${2:-noScaleSmear_noBtagSF_ZPtRW}"
DATA_ERA="${3:-2025}"
MC_ERA="${4:-2024}"

ERA_LABEL="${DATA_ERA}Data_${MC_ERA}MC"

# ---- SET YOUR MACRO DIRECTORY HERE ----
MACRO_DIR="/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/8_Plots/Overlay_AllPlots"

MUON_MACRO="${MACRO_DIR}/Muon_Overlay_MixedYears.C"
JET_MACRO="${MACRO_DIR}/Jet_Overlay_MixedYears.C"

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

DATA_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${DATA_ERA}/${CORR_MODE}/${CONFIG_DIR}"
MC_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${MC_ERA}/${CORR_MODE}/${CONFIG_DIR}"

echo "=================================================="
echo " Running Mixed-Year Muon + Jet Overlay Workflow"
echo " DATA era : ${DATA_ERA}"
echo " MC era   : ${MC_ERA}"
echo " Label    : ${ERA_LABEL}"
echo " DATA path: ${DATA_PATH}"
echo " MC path  : ${MC_PATH}"
echo "=================================================="

if [ ! -d "$DATA_PATH" ]; then
    echo "[ERROR] DATA directory does not exist:"
    echo "  $DATA_PATH"
    exit 1
fi

if [ ! -d "$MC_PATH" ]; then
    echo "[ERROR] MC directory does not exist:"
    echo "  $MC_PATH"
    exit 1
fi

# ------------------- Muon Overlay -------------------
echo "[STEP 1] Muon mixed-year overlay plots"

root -l -b <<EOF
.L ${MUON_MACRO}
Muon_Overlay_MixedYears("${DATA_PATH}", "${MC_PATH}", "${ERA_LABEL}");
.q
EOF

echo "Finished muon overlays: ${ERA_LABEL}"

# ------------------- Jet Overlay -------------------
echo "[STEP 2] Jet mixed-year overlay plots"

root -l -b <<EOF
.L ${JET_MACRO}
Jet_Overlay_MixedYears("${DATA_PATH}", "${MC_PATH}", "${ERA_LABEL}");
.q
EOF

echo "Finished jet overlays: ${ERA_LABEL}"

echo "=================================================="
echo "Done."
echo "=================================================="
