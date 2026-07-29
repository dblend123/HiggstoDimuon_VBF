#!/bin/bash

files=(
  jetvetomaps12022.json
  jetvetomaps22022.json
  jetvetomaps12023.json
  jetvetomaps22023.json
  jetvetomaps2024.json
  jetvetomaps2025.json
)

for f in "${files[@]}"; do
    echo
    echo "==== $f ===="

    jq -r '.. | .name? // empty' "$f" | sort -u

    echo
    jq -r '.. | .key? // empty' "$f" | sort -u
done
