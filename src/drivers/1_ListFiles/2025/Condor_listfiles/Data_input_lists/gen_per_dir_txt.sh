#!/usr/bin/env bash

#Example: ./gen_per_dir_txt.sh /afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/2025/Data_listfiles

set -euo pipefail

BASEDIR="${1:-/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/2025/Data_listfiles}"
SUFFIX="${SUFFIX:-DataRun2025}"

declare -A ALIAS=(
  ["Muon0"]="Muon0"
  ["Muon1"]="Muon1"
)

if [[ ! -d "$BASEDIR" ]]; then
  echo "ERROR: directory not found: $BASEDIR" >&2
  exit 2
fi

echo "Scanning process folders in: $BASEDIR"
shopt -s nullglob

for procdir in "$BASEDIR"/*/; do
  [[ -d "$procdir" ]] || continue

  procbase="$(basename "$procdir")"

  # Skip helper dirs like latest_only if present
  if [[ "$procbase" != "Muon0" && "$procbase" != "Muon1" ]]; then
    echo "Skipping non-data folder: $procbase"
    continue
  fi

  alias="${ALIAS[$procbase]:-$procbase}"
  outfile="${BASEDIR}/${procbase}.txt"

  echo "→ Building ${outfile}"
  : > "$outfile"

  find "$procdir" -maxdepth 1 -type f -name '*.list' -print0 \
    | sort -zV \
    | while IFS= read -r -d '' f; do
        base="$(basename "$f" .list)"
        # Example base: Muon0_Cpart017

        run_letter="$(sed -E 's/^.*_([A-Z])part[0-9]+.*$/\1/' <<< "$base")"
        idx="$(sed -E 's/^.*part([0-9]+).*$/\1/' <<< "$base")"

        if [[ -z "$run_letter" || "$run_letter" == "$base" ]]; then
          echo "WARNING: could not extract run letter from $base; using X" >&2
          run_letter="X"
        fi

        if [[ -z "$idx" || "$idx" == "$base" ]]; then
          echo "WARNING: could not extract part index from $base; using 000" >&2
          idx="000"
        fi

        ext=""
        [[ "$base" =~ [Ee][Xx][Tt] ]] && ext="_Ext"

        # Old: Muon0_017DataRun2025
        # New: Muon0_C017DataRun2025
        tag="${alias}_${run_letter}${idx}${ext}${SUFFIX:+$SUFFIX}"

        printf "%s %s\n" "$f" "$tag" >> "$outfile"
      done

  if [[ ! -s "$outfile" ]]; then
    echo "  (no .list files found in $procdir — skipping)"
    rm -f "$outfile"
  else
    echo "  Wrote $(wc -l < "$outfile") lines."
  fi
done

echo "Done."
