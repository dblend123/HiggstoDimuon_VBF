#!/bin/bash

# Usage:
#   ./output_lumi_JSONfiles.sh 2025/*Golden.json
# or:
#   ./output_lumi_JSONfiles.sh 2025/Cert_*.json

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 path/to/*.json"
    exit 1
fi

mkdir -p lumi_outputs

for json_file in "$@"; do
    echo ""
    echo "📁 Processing: $json_file"

    if [[ ! -f "$json_file" ]]; then
        echo "❌ Error: File '$json_file' not found! Skipping..."
        continue
    fi

    base_name=$(basename "$json_file" .json)
    output_file="lumi_outputs/${base_name}_lumi.txt"

    brilcalc lumi \
      --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json \
      -u /fb \
      -i "$json_file" | tee "$output_file"

    total_lumi=$(awk -F'|' '
        /totrecorded/ {getline; getline; print $6}
    ' "$output_file" | xargs)

    echo "✅ Output saved to: $output_file"
    echo "➡️  Total recorded luminosity: ${total_lumi} /fb"
    echo "---------------------------------------------"
done
