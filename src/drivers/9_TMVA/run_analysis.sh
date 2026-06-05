#!/bin/bash
# run_analysis.sh — Post-scoring analysis: per-era diagnostics + per-era
# categorization + combined multi-era categorization.
#
# Usage:
#   ./run_analysis.sh                 # all eras (default)
#   ./run_analysis.sh -e 12022        # single era (skips the combined step)

set -e
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
cd /eos/user/m/memmad/HiggstoDimuon_VBF
export MPLBACKEND=Agg

OUT=src/drivers/9_TMVA/Outputs
DIAG_SCRIPT=src/drivers/9_TMVA/vbf_bdt/bdt_validation.py
CAT_SCRIPT=src/drivers/9_TMVA/vbf_bdt/bdt_event_categorization.py
ALL_ERAS=(12022 22022 12023 22023 2024 2025)
ERAS=("${ALL_ERAS[@]}")

while [[ $# -gt 0 ]]; do
    case $1 in
        -e|--era) ERAS=("$2"); shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Step 1 — Per-era diagnostics (KS folds, score ROC, feature checks)
for era in "${ERAS[@]}"; do
    echo ">>> [Diagnostics] era=${era}"
    ODIR="$OUT/VBFValidation/$era"; mkdir -p "$ODIR"
    python3 "$DIAG_SCRIPT" --era "$era" --outputs-dir "$OUT" \
        --output-dir "$ODIR" --log-level INFO |& tee "$ODIR/run.log"
done

# Step 2 — Per-era categorization
for era in "${ERAS[@]}"; do
    echo ">>> [Categorization] era=${era}"
    CATDIR="$OUT/Categorization/$era"; mkdir -p "$CATDIR"
    python3 "$CAT_SCRIPT" --eras "$era" --outputs-dir "$OUT" \
        --outdir "$CATDIR" --log-level INFO |& tee "$CATDIR/run.log"
done

# Step 3 — Combined multi-era categorization (skipped for single-era runs)
if [[ ${#ERAS[@]} -gt 1 ]]; then
    echo ">>> [Categorization] combined (${ALL_ERAS[*]})"
    CATDIR="$OUT/Categorization/combined"; mkdir -p "$CATDIR"
    python3 "$CAT_SCRIPT" --eras "${ALL_ERAS[@]}" --outputs-dir "$OUT" \
        --outdir "$CATDIR" --log-level INFO |& tee "$CATDIR/run.log"
fi
