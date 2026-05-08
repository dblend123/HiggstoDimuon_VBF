#!/usr/bin/env bash
set -euo pipefail

unset PYTHONHOME
unset PYTHONPATH

# ---- config ----
REDIR="root://cms-xrd-global.cern.ch"
LINES_PER_LIST="${LINES_PER_LIST:-100}"   # override: LINES_PER_LIST=150 ./make_lists_data.sh
OUTROOT="lists"                           # base output dir
SHUFFLE="${SHUFFLE:-yes}"                 # set SHUFFLE=no to keep DAS order

# ---- DATASETS: data only ----
DATASETS=$(cat <<'EOF'
/Muon0/Run2025B-PromptReco-v1/NANOAOD
/Muon0/Run2025C-PromptReco-v1/NANOAOD
/Muon0/Run2025C-PromptReco-v2/NANOAOD
/Muon0/Run2025D-PromptReco-v1/NANOAOD
/Muon0/Run2025E-PromptReco-v1/NANOAOD
/Muon0/Run2025F-PromptReco-v1/NANOAOD
/Muon0/Run2025F-PromptReco-v2/NANOAOD
/Muon0/Run2025G-PromptReco-v1/NANOAOD

/Muon1/Run2025B-PromptReco-v1/NANOAOD
/Muon1/Run2025C-PromptReco-v1/NANOAOD
/Muon1/Run2025C-PromptReco-v2/NANOAOD
/Muon1/Run2025D-PromptReco-v1/NANOAOD
/Muon1/Run2025E-PromptReco-v1/NANOAOD
/Muon1/Run2025F-PromptReco-v1/NANOAOD
/Muon1/Run2025F-PromptReco-v2/NANOAOD
/Muon1/Run2025G-PromptReco-v1/NANOAOD
EOF
)

# ---- helpers ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH"; exit 2; }; }

# Extract leading tag: Muon0 / Muon1
extract_muon_tag() {
  awk -F'/' '{print $2}' <<< "$1"
}

# Extract run block: Run2025B (take part before first '-')
extract_run_block() {
  awk -F'/' '{split($3,a,"-"); print a[1]}' <<< "$1"
}

echo "Checking dependencies..."
need dasgoclient
need split
if [[ "$SHUFFLE" == "yes" ]]; then need shuf; fi

mkdir -p "$OUTROOT/tmp"

echo "Starting… (LINES_PER_LIST=$LINES_PER_LIST, OUTROOT=$OUTROOT, SHUFFLE=$SHUFFLE)"

# 1) Collect files per (MuonX, run-letter) into tmp text files
while IFS= read -r ds; do
  # trim whitespace
  ds="${ds#"${ds%%[![:space:]]*}"}"
  ds="${ds%"${ds##*[![:space:]]}"}"

  # skip blanks and comments
  [[ -z "$ds" ]] && continue
  [[ "${ds:0:1}" == "#" ]] && continue

  # ensure leading slash
  [[ "${ds:0:1}" != "/" ]] && ds="/$ds"

  muon_tag=$(extract_muon_tag "$ds")      # Muon0 or Muon1
  run_block=$(extract_run_block "$ds")    # Run2025B, Run2025C, ...
  run_letter="${run_block: -1}"           # B, C, D, ...

  tmp="$OUTROOT/tmp/${muon_tag}_${run_letter}_files.txt"

  echo "Querying DAS for: $ds"
  files=$(dasgoclient -query "file dataset=$ds" --limit=0 || true)
  if [[ -z "$files" ]]; then
    echo "WARNING: no files for $ds"
    continue
  fi

  # prepend redirector and (optionally) shuffle, then append to tmp list
  if [[ "$SHUFFLE" == "yes" ]]; then
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" | shuf >> "$tmp"
  else
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" >> "$tmp"
  fi

done <<< "$DATASETS"

# 2) For each tmp file, split into MuonX_<letter>partNNN.list under lists/MuonX/
for tmpfile in "$OUTROOT"/tmp/*_files.txt; do
  [[ ! -e "$tmpfile" ]] && continue  # handle case of no tmp files

  fname=$(basename "$tmpfile")       # e.g. Muon0_C_files.txt
  muon_tag="${fname%%_*}"            # Muon0
  rest="${fname#*_}"                 # C_files.txt
  run_letter="${rest%%_*}"           # C

  outdir="$OUTROOT/$muon_tag"
  mkdir -p "$outdir"

  prefix="$outdir/${muon_tag}_${run_letter}part"
  echo "Splitting $tmpfile -> ${prefix}NNN.list"

  split -d -a 3 -l "$LINES_PER_LIST" --additional-suffix=.list "$tmpfile" "$prefix"

  echo "→ Wrote $(ls -1 "$outdir"/${muon_tag}_${run_letter}part*.list | wc -l) list(s) to $outdir/"
  rm -f "$tmpfile"
done

rmdir "$OUTROOT/tmp" 2>/dev/null || true

echo "All done."
