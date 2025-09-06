#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
REDIR="root://cms-xrd-global.cern.ch"
LINES_PER_LIST="${LINES_PER_LIST:-40}"   # override: LINES_PER_LIST=150 ./make_lists.sh
OUTROOT="lists"                           # base output dir
SHUFFLE="${SHUFFLE:-yes}"                 # set SHUFFLE=no to keep DAS order

# ---- datasets ----
DATASETS=$(cat <<'EOF'
/VBFHto2Mu_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM

/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM

/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM

/EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM

/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM

/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM

/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM

/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM
/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6_ext1-v2/NANOAODSIM

/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM

/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM

/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM

/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM

/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM

/Muon0/Run2023D-22Sep2023_v1-v1/NANOAOD
/Muon0/Run2023D-22Sep2023_v2-v1/NANOAOD
/Muon1/Run2023D-22Sep2023_v1-v1/NANOAOD
/Muon1/Run2023D-22Sep2023_v2-v1/NANOAOD
EOF
)

# ---- helpers ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH"; exit 2; }; }
process_name() {  # use first token before '_' in primary name (e.g. VBFHto2Mu_M-125 -> VBFHto2Mu; Muon0 stays Muon0)
  awk -F'/' '{split($2,a,"_"); print a[1]}' <<<"$1"
}

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
  files=$(dasgoclient -query "file dataset=$ds")
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
done

echo "All done."

