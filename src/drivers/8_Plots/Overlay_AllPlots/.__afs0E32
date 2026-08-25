#!/bin/bash

# Run all eras through the CSV-enabled systematics printout macro.
# This creates one combined CSV that can be read by make_muon_jet_slides_withSystematics.py.

eras=(
    12022
    22022
    12023
    22023
    2024
    2025
)

macro="PrintSystematics_withCSV.C"
csv_out="/eos/user/d/dblend/PowerPoint/systematics_summary.csv"

# Remove old CSV so the first era writes a fresh header.
rm -f "${csv_out}"

for era in "${eras[@]}"; do
    echo
    echo "=================================================="
    echo "Running ${macro} for ERA: ${era}"
    echo "=================================================="

    root -l -b -q \
    "${macro}(\"/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/KIT/noScaleSmear_noBTagRW_noZPtRW\", \"DiMuonMass_Higgs\", \"corr_DM_VBF_SR\", false, 115.0, 135.0, \"${csv_out}\", \"${era}\")"
done

echo
echo "=================================================="
echo "All eras complete!"
echo "CSV written to: ${csv_out}"
echo "=================================================="
