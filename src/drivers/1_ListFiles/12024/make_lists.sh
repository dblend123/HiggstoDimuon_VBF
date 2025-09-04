#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Config
# -----------------------------
REDIR="root://cms-xrd-global.cern.ch"
LINES_PER_LIST="${LINES_PER_LIST:-200}"   # override: LINES_PER_LIST=150 ./make_lists.sh
OUTROOT="lists"                           # base output dir
SHUFFLE="${SHUFFLE:-yes}"                 # set SHUFFLE=no to keep DAS order

# -----------------------------
# Datasets
# -----------------------------
DATASETS=$(cat <<'EOF'
/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v6/NANOAODSIM
/DYto2Tau-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v5/NANOAODSIM
EOF
)

# -----------------------------
# Helpers
# -----------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH"; exit 2; }; }

# Directory name = first token (before '_') of the primary dataset component
# e.g. /VBFH-Hto2Mu_Par-... -> VBFH-Hto2Mu
process_name() {
  awk -F'/' '{split($2,a,"_"); print a[1]}' <<<"$1"
}

# -----------------------------
# Run
# -----------------------------
echo "Checking dependencies..."
need dasgoclient
need split
if [[ "$SHUFFLE" == "yes" ]]; then need shuf; fi

mkdir -p "$OUTROOT"

echo "Starting… (LINES_PER_LIST=$LINES_PER_LIST, OUTROOT=$OUTROOT, SHUFFLE=$SHUFFLE)"
# shellcheck disable=SC2086
for ds in $DATASETS; do
  [[ -z "${ds// }" ]] && continue
  [[ "${ds:0:1}" != "/" ]] && ds="/$ds"

  pname=$(process_name "$ds")
  outdir="$OUTROOT/$pname"
  tmp="$outdir/files.txt"
  mkdir -p "$outdir"

  echo "Querying DAS for: $ds"
  # Guard against DAS hiccups/typos so the loop continues
  if ! files=$(dasgoclient -query "file dataset=$ds"); then
    echo "ERROR querying DAS for: $ds (skipping)"
    continue
  fi
  if [[ -z "$files" ]]; then
    echo "WARNING: no files for $ds"
    continue
  fi

  # Prepend redirector and optionally shuffle
  if [[ "$SHUFFLE" == "yes" ]]; then
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" | shuf > "$tmp"
  else
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" > "$tmp"
  fi

  # Split into chunks (part000.list, part001.list, ...)
  split -d -a 3 -l "$LINES_PER_LIST" --additional-suffix=.list "$tmp" "$outdir/part"
  rm -f "$tmp"

  echo "→ Wrote $(ls -1 "$outdir"/*.list | wc -l) list(s) to $outdir/"
done

echo "All done."

