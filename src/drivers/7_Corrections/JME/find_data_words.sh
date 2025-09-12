#!/bin/bash

# List of JSON files
files=(
  "jet_jerc12022.json"
  "jet_jerc12023.json"
  "jet_jerc22022.json"
  "jet_jerc22023.json"
)

# Loop through each file
for f in "${files[@]}"; do
  echo "==== $f ===="
  # Extract all words containing DATA
  grep -o '\w*DATA\w*' "$f" | sort -u
done

