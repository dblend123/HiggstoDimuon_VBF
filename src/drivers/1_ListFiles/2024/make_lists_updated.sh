#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
REDIR="root://cms-xrd-global.cern.ch"
LINES_PER_LIST="${LINES_PER_LIST:-40}"   # override: LINES_PER_LIST=150 ./make_lists_updated.sh
OUTROOT="lists"                           # base output dir
SHUFFLE="${SHUFFLE:-yes}"                 # set SHUFFLE=no to keep DAS order

# ---- datasets ----
DATASETS=$(cat <<'EOF'
# --- MC datasets ---
/VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v6/NANOAODSIM
/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM
/EWK-2Mu2J_Bin-MLL-105to160_TuneCH3_13p6TeV_madgraph-herwig7/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/EWK-2L2J_Bin-MLL-50-MJJ-120_TuneCH3_13p6TeV_madgraph-herwig7/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM

/WZZ-5F_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/ZZZ-5F_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/WWW-4F_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
/WWZ-4F_TuneCP5_13p6TeV_amcatnlo-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM

# --- Data datasets ---
/Muon0/Run2024C-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024D-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024E-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024F-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024G-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024H-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024I-MINIv6NANOv15-v1/NANOAOD
/Muon0/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD

/Muon1/Run2024C-MINIv6NANOv15-v1/NANOAOD
/Muon1/Run2024D-MINIv6NANOv15-v1/NANOAOD
/Muon1/Run2024E-MINIv6NANOv15-v1/NANOAOD
/Muon1/Run2024F-MINIv6NANOv15-v1/NANOAOD
/Muon1/Run2024G-MINIv6NANOv15-v2/NANOAOD
/Muon1/Run2024H-MINIv6NANOv15-v2/NANOAOD
/Muon1/Run2024I-MINIv6NANOv15-v1/NANOAOD
/Muon1/Run2024I-MINIv6NANOv15_v2-v1/NANOAOD
EOF
)

# ---- helpers ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH"; exit 2; }; }
process_name() {  # use first token before '_' in primary name (e.g. VBFH-Hto2Mu_Par-M-125 -> VBFH-Hto2Mu; Muon0 stays Muon0)
  awk -F'/' '{split($2,a,"_"); print a[1]}' <<<"$1"
}

echo "Checking dependencies..."
need dasgoclient
need split
if [[ "$SHUFFLE" == "yes" ]]; then need shuf; fi

mkdir -p "$OUTROOT"

echo "Starting… (LINES_PER_LIST=$LINES_PER_LIST, OUTROOT=$OUTROOT, SHUFFLE=$SHUFFLE)"

# iterate line-by-line to preserve datasets and skip comments/blank lines
while IFS= read -r ds; do
  # trim leading/trailing whitespace
  ds="${ds#"${ds%%[![:space:]]*}"}"
  ds="${ds%"${ds##*[![:space:]]}"}"

  # skip blanks and comments
  [[ -z "$ds" ]] && continue
  [[ "${ds:0:1}" == "#" ]] && continue

  # ensure leading slash
  [[ "${ds:0:1}" != "/" ]] && ds="/$ds"

  pname=$(process_name "$ds")
  outdir="$OUTROOT/$pname"
  tmp="$outdir/files.txt"
  mkdir -p "$outdir"

  echo "Querying DAS for: $ds"
  files=$(dasgoclient -query "file dataset=$ds" || true)
  if [[ -z "$files" ]]; then
    echo "WARNING: no files for $ds"
    continue
  fi

  if [[ "$SHUFFLE" == "yes" ]]; then
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" | shuf > "$tmp"
  else
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" > "$tmp"
  fi

  split -d -a 3 -l "$LINES_PER_LIST" --additional-suffix=.list "$tmp" "$outdir/part"
  rm -f "$tmp"

  echo "→ Wrote $(ls -1 "$outdir"/*.list | wc -l) list(s) to $outdir/"
done <<< "$DATASETS"

echo "All done."
