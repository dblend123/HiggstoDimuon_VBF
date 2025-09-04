#!/usr/bin/env bash
set -euo pipefail

# Base dir that contains per-process folders (each with part*.list)
# Pass as arg 1, or it’ll use this default.
BASEDIR="${1:-/afs/cern.ch/user/d/dblend/test-analysis/src/drivers/1_ListFiles/22022/Data_listfiles}"

# Optional: suffix appended to each tag (leave empty for none)
# Example values: MCRun2022EE, MCRun23C, MCRun24v15, etc.
SUFFIX="${SUFFIX:-MCRun22022}"

# Optional: alias map (folder name -> short tag prefix).
# If a folder isn’t listed here, its folder name is used as-is.
declare -A ALIAS=(
  ["E"]="E"
  ["F"]="F"
  ["G"]="G"
)

if [[ ! -d "$BASEDIR" ]]; then
  echo "ERROR: directory not found: $BASEDIR" >&2
  exit 2
fi

echo "Scanning process folders in: $BASEDIR"
shopt -s nullglob

# Iterate immediate subdirectories (process folders)
for procdir in "$BASEDIR"/*/; do
  [[ -d "$procdir" ]] || continue
  procbase="$(basename "$procdir")"
  alias="${ALIAS[$procbase]:-$procbase}"
  outfile="${BASEDIR}/${procbase}.txt"

  echo "→ Building ${outfile}"
  : > "$outfile"

  # Find *.list files and write lines: "<full_path> <TAG>"
  # Natural sort so part000, part001, ..., part010
  find "$procdir" -maxdepth 1 -type f -name '*.list' -print0 \
    | sort -zV \
    | while IFS= read -r -d '' f; do
        base="$(basename "$f" .list)"

        # Extract trailing number for index (defaults to 0 if missing)
        idx="$(grep -oE '[0-9]+' <<< "$base" | tail -n1 || true)"
        [[ -z "$idx" ]] && idx="0"

        # Optional "_Ext" marker if filename mentions Ext (case-insensitive)
        ext=""
        [[ "$base" =~ [Ee][Xx][Tt] ]] && ext="_Ext"

        # Compose tag. Matches your style: e.g. DY_0MCRun22CD (no extra underscore)
        tag="${alias}_${idx}${ext}${SUFFIX:+$SUFFIX}"

        printf "%s %s\n" "$f" "$tag" >> "$outfile"
      done

  # If no lines were written, remove the empty file
  if [[ ! -s "$outfile" ]]; then
    echo "  (no .list files found in $procdir — skipping)"
    rm -f "$outfile"
  else
    echo "  Wrote $(wc -l < "$outfile") lines."
  fi
done

echo "Done."

