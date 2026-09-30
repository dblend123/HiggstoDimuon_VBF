#!/bin/bash
# run_full.sh — End-to-end BDT pipeline (prepare → train → test → score) for one era.
#
# Usage:
#   ./run_full.sh -e <era> [--base <path>] [--skip-training] [--skip-testing] [--skip-scoring]
#
# Defaults:
#   --base = /eos/user/d/dblend/Analyzed_SkimmedFiles/<era>/Rochester/noScaleSmear_noBtagSF

set -euo pipefail
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ERA=""; BASE_PATH=""
SKIP_TRAINING=false; SKIP_TESTING=false; SKIP_SCORING=false

usage() { sed -n '1,12p' "$0"; exit 1; }

while [[ $# -gt 0 ]]; do
    case $1 in
        -e|--era)         ERA="$2";   shift 2 ;;
        --base)           BASE_PATH="$2"; shift 2 ;;
        --skip-training)  SKIP_TRAINING=true; shift ;;
        --skip-testing)   SKIP_TESTING=true;  shift ;;
        --skip-scoring)   SKIP_SCORING=true;  shift ;;
        --help|-h)        usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done
[[ -z "$ERA" ]] && { echo "ERROR: --era is required"; usage; }
[[ -z "$BASE_PATH" ]] && \
    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/KIT/noScaleSmear_BTagRW_ZPtRW"

DPREP="DataPreparation_Files/${ERA}"
MODELDIR="Outputs/Training/${ERA}"
SCOREDIR="Outputs/${ERA}"
mkdir -p "$MODELDIR" "$SCOREDIR"

echo "========================================="
echo "TMVA BDT workflow — era ${ERA}"
echo "  base    = ${BASE_PATH}"
echo "========================================="

# Step 0: Build the per-era snapshots used by everything downstream
python3 tmva100_DataPreparation.py \
    --era "${ERA}" \
    --base "${BASE_PATH}" \
    --signal-subdir VBF \
    --treename-in mva \
    --make-all \
    --keep-branches "m_mumu,run,lumi,event,weight_puUp,weight_puDown,weight_muSFUp,weight_muSFDown,weight_triggerUp,weight_triggerDown,weight_btagSFUp,weight_btagSFDown,weight_btagEffStatUp,weight_btagEffStatDown,weight_zptUp,weight_zptDown,weight_scale0,weight_scale1,weight_scale3,weight_scale5,weight_scale7,weight_scale8,weight_pdfUp,weight_pdfDown,weight_isrUp,weight_isrDown,weight_fsrUp,weight_fsrDown"

# Step 1: Train one BDT against the merged background
if [[ "$SKIP_TRAINING" = false ]]; then
    echo ">>> [1/3] Training BDT"
    python3 tmva101_Training.py --era "${ERA}" --base-dir "${DPREP}" --bkg-cats background
    echo "✓ training done"
else
    echo ">>> [1/3] training skipped"
fi

MODEL_PATTERN="${MODELDIR}/${ERA}_tmva101_background_fold{fold}.root"

# Step 2: ROC / AUC on held-out fold==3
if [[ "$SKIP_TESTING" = false ]]; then
    echo ">>> [2/3] Testing"
    python3 tmva102_Testing.py --era "${ERA}" --dir "${DPREP}" --model "${MODEL}"
    echo "✓ testing done"
else
    echo ">>> [2/3] testing skipped"
fi

# Step 3: Score signal, per-cat backgrounds, merged background, and data
if [[ "$SKIP_SCORING" = false ]]; then
    echo ">>> [3/3] Scoring all samples"
    score() {
        local in="$1" out="$2"
        python3 tmva103_Application.py \
            --input "$in" \
            --output "$out" \
            --treename mva \
            --model-pattern "${MODEL_PATTERN}" \
            --model-name myBDT
    }
    score "${DPREP}/${ERA}_all_signal.root"     "${SCOREDIR}/${ERA}_scored_signal.root"
    score "${DPREP}/${ERA}_all_background.root" "${SCOREDIR}/${ERA}_scored_background.root"
    for cat in DY EWK TT VV VVV; do
        in="${DPREP}/${ERA}_all_${cat}.root"
        [[ -f "$in" ]] && score "$in" "${SCOREDIR}/${ERA}_scored_${cat}.root"
    done
    score "${BASE_PATH}/DATA/DATA_all.root" "${SCOREDIR}/${ERA}_data_scored.root"
    echo "✓ scoring done"
else
    echo ">>> [3/3] scoring skipped"
fi

if [[ "$SKIP_TESTING" = false ]]; then
    echo ">>> Testing combined out-of-fold scores"

    python3 tmva102_Testing.py \
        --era "${ERA}" \
        --signal "${SCOREDIR}/${ERA}_scored_signal.root" \
        --background "${SCOREDIR}/${ERA}_scored_background.root" \
        --treename mva \
        --output-dir "Outputs/Testing/${ERA}"

    echo "✓ testing done"
else
    echo ">>> testing skipped"
fi

echo "========================================="
echo "✓ done.  Models: ${MODELDIR}/   Scored: ${SCOREDIR}/"
echo "========================================="
