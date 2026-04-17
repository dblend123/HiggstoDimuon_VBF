#!/bin/bash

# Default settings (can override via arguments)
CORR_MODE="${1:-KIT}"
CONFIG_DIR="${2:-noPU_noScaleSmear_noMuonSF_noJERC_noJetVeto_noGolden_noFlags_noBtagSF_noJetCriteria}"

#ERAS=(12022 22022 12023 22023 2024)
ERAS=(2024)

for ERA in "${ERAS[@]}"; do

    BASE_PATH="/eos/user/d/dblend/Analyzed_SkimmedFiles/${ERA}/${CORR_MODE}/${CONFIG_DIR}"

    echo "=================================================="
    echo " Running Pileup Overlay for ERA: $ERA"
    echo " Path: $BASE_PATH"
    echo "=================================================="

    if [ ! -d "$BASE_PATH" ]; then
        echo "Directory does not exist, skipping: $BASE_PATH"
        echo
        continue
    fi

    root -l -b <<EOF
.L Pileup_Overlay_AllPlots.C
Pileup_Overlay_AllPlots("${BASE_PATH}", "${ERA}");
.q
EOF

    echo "Finished ERA: $ERA"
    echo
done
