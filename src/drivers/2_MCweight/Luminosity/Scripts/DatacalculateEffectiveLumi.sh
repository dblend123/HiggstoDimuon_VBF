#!/bin/bash

# EOS directory with JSONs
EOS_DIR="/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/13p6TeV"

# Output file
OUTFILE="lumis_2022.csv"
echo "era,lumi_fb" > "$OUTFILE"

# List of eras to check (manually controlled to avoid duplicates)
for ERA in B C D E F G; do
    JSON="$EOS_DIR/Cert_Collisions2022_era${ERA}_*_Golden.json"

    # Expand only the first match (avoid duplicates)
    MATCHED_JSON=$(ls $JSON 2>/dev/null | head -n 1)

    if [[ ! -f "$MATCHED_JSON" ]]; then
        echo "Skipping 2022$ERA — JSON not found"
        continue
    fi

    echo "Processing Run 2022$ERA..."

    # Run brilcalc and capture recorded lumi
    RECORDED=$(brilcalc lumi \
        --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json \
        -u /fb -i "$MATCHED_JSON" \
        | awk '/^#Tot recorded/ {print $4}')

    if [[ -z "$RECORDED" ]]; then
        echo "Warning: no recorded lumi found for 2022$ERA"
        RECORDED="0.0"
    fi

    echo "2022$ERA,$RECORDED" >> "$OUTFILE"
done

echo "✅ Saved luminosities to $OUTFILE"

