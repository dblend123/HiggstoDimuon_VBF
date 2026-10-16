#!/bin/bash
set -euo pipefail

SAMPLE="$1"
CHUNK="$2"

cd "$HOME/test-analysis"

EXE="$PWD/build/src/drivers/process_HiggsAnalysis_wCutsExtended_NoPairing_v0"
INPUT="$PWD/sync_jobs/chunks/$CHUNK"

if [[ "$SAMPLE" == "VBF" ]]; then
    ISMC=1
else
    ISMC=0
fi

OUTPUT="${CHUNK}_sync_v6.root"
LOCAL_OUTPUT="${TMPDIR:-/tmp}/$OUTPUT"

echo "Host: $(hostname)"
echo "Sample: $SAMPLE"
echo "Input: $INPUT"

"$EXE" \
    --input "$INPUT" \
    --isMC "$ISMC" \
    --sample "$SAMPLE" \
    --era 2024 \
    --muCorr kit \
    --noBtagSF \
    --noZPtRW \
    --sync \
    --output "$LOCAL_OUTPUT"

test -s "$LOCAL_OUTPUT"

DEST="root://eosuser.cern.ch//eos/user/d/dblend/Hmumu_Synchronization/2024/${SAMPLE}/${OUTPUT}"

xrdcp -f "$LOCAL_OUTPUT" "$DEST"

echo "Finished: $DEST"
