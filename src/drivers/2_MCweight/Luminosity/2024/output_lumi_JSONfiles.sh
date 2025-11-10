#!/bin/bash

# Prompt for multiple JSON files
read -p "Enter the path(s) to your JSON file(s), separated by space: " -a json_files

# Loop over each JSON file
for json_file in "${json_files[@]}"; do

    echo ""
    echo "📁 Processing: $json_file"

    # Check if file exists
    if [[ ! -f "$json_file" ]]; then
        echo "❌ Error: File '$json_file' not found! Skipping..."
        continue
    fi

    # Derive output filename
    base_name=$(basename "$json_file" .json)
    output_file="${base_name}_lumi.txt"

    # Run brilcalc and save output
    brilcalc lumi \
      --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json \
      -u /fb \
      -i "$json_file" | tee "$output_file"

    # Extract total lumi
    total_lumi=$(grep -vE '^\s*#' "$output_file" | tail -n 1 | awk '{print $10}')

    # Display result
    echo "✅ Output saved to: $output_file"
    echo "➡️  Total recorded luminosity: ${total_lumi} /fb"
    echo "---------------------------------------------"

done

