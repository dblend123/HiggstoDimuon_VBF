#!/bin/bash

output_file="combined_output.txt"
> "$output_file"  # clear previous content

for txtfile in *.txt; do
    grep -E '^(Processing|nSigP|nSigF|Tag-and-Probe|Manual|Chi|NLL)' "$txtfile" >> "$output_file"
done

