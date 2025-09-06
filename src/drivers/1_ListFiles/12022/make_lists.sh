#!/usr/bin/env bash
set -euo pipefail

# ---- config ----
REDIR="root://cms-xrd-global.cern.ch"
LINES_PER_LIST="${LINES_PER_LIST:-40}"   # override: LINES_PER_LIST=150 ./make_lists.sh
OUTROOT="lists"                           # base output dir
SHUFFLE="${SHUFFLE:-yes}"                 # set SHUFFLE=no to keep DAS order
SUFFIX_LEN=3                              # part000.list, part001.list, ...
DEDUP="${DEDUP:-yes}"                     # set DEDUP=no to skip de-dup across reruns

# ---- datasets ----
DATASETS=$(cat <<'EOF'
/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v1/NANOAODSIM

/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM

/EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM

/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM

/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM

/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM

/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v3/NANOAODSIM

/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
/ZZto4L_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM

/WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM

/WWZ_4F_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM

/WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM

/ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM
EOF
)

# ---- helpers ----
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH"; exit 2; }; }
process_name() {  # use first token before '_' in primary name
  awk -F'/' '{split($2,a,"_"); print a[1]}' <<<"$1"
}

echo "Checking dependencies..."
need dasgoclient
need split
if [[ "$SHUFFLE" == "yes" ]]; then need shuf; fi

mkdir -p "$OUTROOT"

echo "Starting… (LINES_PER_LIST=$LINES_PER_LIST, OUTROOT=$OUTROOT, SHUFFLE=$SHUFFLE, DEDUP=$DEDUP)"
# shellcheck disable=SC2086
for ds in $DATASETS; do
  [[ -z "${ds// }" ]] && continue
  [[ "${ds:0:1}" != "/" ]] && ds="/$ds"

  pname=$(process_name "$ds")
  outdir="$OUTROOT/$pname"
  tmp="$outdir/files.txt"
  mkdir -p "$outdir"

  echo "Querying DAS for: $ds"

  # Make dasgoclient failures non-fatal with set -e
  if ! files=$(dasgoclient -query "file dataset=$ds" 2>/dev/null); then
    echo "WARNING: dasgoclient failed for $ds; skipping."
    continue
  fi
  if [[ -z "$files" ]]; then
    echo "WARNING: no files for $ds"
    continue
  fi

  # Build list of full xrootd paths
  if [[ "$SHUFFLE" == "yes" ]]; then
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" | shuf > "$tmp"
  else
    printf "%s\n" "$files" | sed "s#^#${REDIR}/#" > "$tmp"
  fi

  # Optional: de-duplicate vs anything already present in outdir (part*.list)
  if [[ "$DEDUP" == "yes" ]]; then
    shopt -s nullglob
    existing_lists=( "$outdir"/part*.list )
    if ((${#existing_lists[@]})); then
      already="$outdir/.already.txt"
      cat "${existing_lists[@]}" | sort -u > "$already"
      # Filter out lines already used (ignore exit code for no matches)
      grep -vxFf "$already" "$tmp" > "$tmp.new" || true
      mv -f "$tmp.new" "$tmp"
      rm -f "$already"
    fi
    shopt -u nullglob
  fi

  # Nothing new? skip
  if [[ ! -s "$tmp" ]]; then
    echo "→ No new files to add for $pname"
    rm -f "$tmp"
    continue
  fi

  # --- find next numeric suffix safely ---
  shopt -s nullglob
  existing_lists=( "$outdir"/part*.list )
  if ((${#existing_lists[@]})); then
    last=$(
      printf '%s\n' "${existing_lists[@]}" \
        | sed -E 's/.*part([0-9]+)\.list/\1/' \
        | sort -n | tail -1
    )
    start=$((10#$last + 1))
  else
    start=0
  fi
  shopt -u nullglob

  # --- split into a temp dir first, then rename without overwriting ---
  tmpdir="$(mktemp -d)"
  # Make chunks: chunk000.list, chunk001.list, ...
  split -d -a "$SUFFIX_LEN" -l "$LINES_PER_LIST" \
        --additional-suffix=.list \
        "$tmp" "$tmpdir/chunk"

  # Rename sequentially to part%03d.list starting at $start
  idx="$start"
  for f in "$tmpdir"/chunk*.list; do
    [[ -e "$f" ]] || continue
    printf -v out "%s/part%0*d.list" "$outdir" "$SUFFIX_LEN" "$idx"
    while [[ -e "$out" ]]; do
      idx=$((idx+1))
      printf -v out "%s/part%0*d.list" "$outdir" "$SUFFIX_LEN" "$idx"
    done
    mv "$f" "$out"
    idx=$((idx+1))
  done
  rmdir "$tmpdir"
  rm -f "$tmp"

  count_now=$(ls -1 "$outdir"/part*.list 2>/dev/null | wc -l | awk '{print $1}')
  echo "→ Wrote/now have $count_now list(s) in $outdir/"
done

echo "All done."

