#!/bin/bash

eras=(
    12022
    22022
    12023
    22023
    2024
    2025
)

for era in "${eras[@]}"; do
    echo
    echo "=================================================="
    echo "Running PrintSystematics for ERA: $era"
    echo "=================================================="

    root -l -b -q \
    "PrintSystematics.C(\"/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/KIT/noScaleSmear_noJERC_noJER_noBTagRW_noZPtRW\")"
done

echo
echo "=================================================="
echo "All eras complete!"
echo "=================================================="
