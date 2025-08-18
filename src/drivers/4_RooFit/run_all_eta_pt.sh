#!/bin/bash

# List of eta tags
eta_tags=("eta_0" "eta_1" "eta_2" "eta_3" "eta_4")

# List of pt tags
pt_tags=("pt_0" "pt_1" "pt_2" "pt_3" "pt_4" "pt_5" "pt_6" "pt_7")

# Path to your ROOT macro
# macro="Voigtian_ZMassRegion.C"
macro="CBG_ZMassRegion.C"
# name="Voigitian"
name="CBG"

# Output directory for logs
mkdir -p logs

# Run for each eta tag
for tag in "${eta_tags[@]}"; do
    echo "Running for $tag..."
    root -l -b -q "${macro}+(\"$tag\")" > "logs/log_${name}_${tag}.txt"
done

# Run for each pt tag
for tag in "${pt_tags[@]}"; do
    echo "Running for $tag..."
    root -l -b -q "${macro}+(\"$tag\")" > "logs/log_${name}_${tag}.txt"

done
