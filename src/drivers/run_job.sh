#!/bin/bash
set -euo pipefail

# Prefer Condor-provided proxy; fall back to shipped file
if [[ -z "${X509_USER_PROXY:-}" || ! -f "${X509_USER_PROXY}" ]]; then
  export X509_USER_PROXY=x509up.cert
fi

# (Optional) CMS env for xrdcp/root tools (harmless if absent)
if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
  # shellcheck disable=SC1091
  source /cvmfs/cms.cern.ch/cmsset_default.sh
fi

# --------------------
# Args from JDL
# --------------------
inputList_remote="${1:?missing inputList}"
tag="${2:?missing tag}"
sample="${3:-DATA}"          # e.g. DY, TT, WZ, ... or DATA
weights_remote="${4:-None}"  # transferred by Condor (we'll use its basename if present)
maxEvents="${5:--1}"         # -1 = all
era="${6:-22022}"
muCorr="${7:-roch}"

# Booleans from JDL (8–16)
noPU="${8:-false}"
noScaleSmear="${9:-false}"
noMuonSF="${10:-false}"
noJERC="${11:-false}"
noJER="${12:-false}"
noJetVeto="${13:-false}"
noGolden="${14:-false}"
noFlags="${15:-false}"
noFSR="${16:-false}"
noBtagSF="${17:-false}"
noJetCriteria="${18:-false}"

noZPtRW="${19:-true}"
zptRWFile_remote="${20:-None}"
zptRWFile="$(basename -- "$zptRWFile_remote")"

btagEffFile_remote="${21:-None}"
btagEffFile="$(basename -- "$btagEffFile_remote")"

noBeamConstraint="${22:-false}"

echo "=== Job $(date) ==="
echo "tag=$tag  sample=$sample  maxEvents=$maxEvents"
echo "inputList (remote) = $inputList_remote"
echo "weights (remote)   = $weights_remote"
echo "Toggles: noPU=$noPU noScaleSmear=$noScaleSmear noMuonSF=$noMuonSF noJERC=$noJERC noJER=$noJER noJetVeto=$noJetVeto noGolden=$noGolden noFlags=$noFlags noFSR=$noFSR noBtagSF=$noBtagSF noJetCriteria=$noJetCriteria noZPtRW=$noZPtRW noBeamConstraint=$noBeamConstraint"

# Resolve to local filenames (Condor puts transferred files in CWD under their basename)
inputList="$(basename -- "$inputList_remote")"
echo "inputList (local)  = $inputList"
[[ -s "$inputList" ]] || { echo "ERROR: list '$inputList' missing in $(pwd)"; ls -la; exit 66; }

# ------------------------------------------------------------------
# Redirect skimmed-file inputs to the shared working area.
# This changes ONLY where the analyzer reads skimmed ROOT files from.
# Stage-out below still writes to /eos/user/d/dblend/Analyzed_SkimmedFiles.
#
# Override at submit time only if needed, e.g.
#   environment = "SKIMMED_BASE=/some/other/base"
# ------------------------------------------------------------------
SKIMMED_BASE="${SKIMMED_BASE:-/eos/experiment/higgstodimuon-vbf}"
inputList_original="${inputList}.original_paths"
cp -p "$inputList" "$inputList_original"

python3 - "$inputList" "$SKIMMED_BASE" <<'PY'
import re
import sys
from pathlib import Path

list_path = Path(sys.argv[1])
new_base = sys.argv[2].rstrip('/')

# Known old skimmed-file bases that may appear in the transferred list files.
old_bases = [
    "/eos/user/d/dblend/SkimmedFiles",
    "/eos/user/o/okoseyan/Personal/SkimmedFiles",
]

text = list_path.read_text()
for old in old_bases:
    text = text.replace(old, new_base)
    # Also handle paths already written with an EOS xrootd prefix.
    text = text.replace(f"root://eosuser.cern.ch//{old}", f"root://eosuser.cern.ch//{new_base}")
    text = text.replace(f"root://eoscms.cern.ch//{old}", f"root://eoscms.cern.ch//{new_base}")

list_path.write_text(text)
PY

echo "Skimmed input base redirected to: $SKIMMED_BASE"
echo "First usable input after redirect:"
awk 'NF && $1 !~ /^#/{print "  "$1; exit}' "$inputList" || true

# Count usable lines (non-empty, not comments)
nfiles=$(awk 'NF && $1 !~ /^#/{c++} END{print c+0}' "$inputList")
echo "Usable files in list: $nfiles"
if [[ "$nfiles" -eq 0 ]]; then
  echo "ERROR: no usable entries in '$inputList'"; exit 65
fi

# DATA vs MC + weights
isMC=true
[[ "$sample" == "DATA" ]] && isMC=false

weights_arg=()
if $isMC; then
  weights_base="$(basename -- "${weights_remote:-}" || true)"
  if [[ -n "${weights_remote:-}" && "$weights_remote" != "None" && -s "$weights_base" ]]; then
    echo "Using weights file : $weights_base"
    weights_arg=( --weights "$weights_base" )
  else
    echo "No usable weights file found; using built-in defaults."
  fi
fi

# Build correction toggle flags for the C++ binary
extra_args=()
[[ "$noPU"         == "true" ]] && extra_args+=( --noPU )
[[ "$noScaleSmear" == "true" ]] && extra_args+=( --noScaleSmear )
[[ "$noMuonSF"     == "true" ]] && extra_args+=( --noMuonSF )
[[ "$noJERC"       == "true" ]] && extra_args+=( --noJERC )
[[ "$noJER"        == "true" ]] && extra_args+=( --noJER )
[[ "$noJetVeto"    == "true" ]] && extra_args+=( --noJetVeto )
[[ "$noGolden"     == "true" ]] && extra_args+=( --noGolden )
[[ "$noFlags"      == "true" ]] && extra_args+=( --noFlags )
[[ "$noFSR"        == "true" ]] && extra_args+=( --noFSR )
[[ "$noBtagSF"     == "true" ]] && extra_args+=( --noBtagSF )
[[ "$noJetCriteria" == "true" ]] && extra_args+=( --noJetCriteria )
[[ "$noZPtRW"       == "true" ]] && extra_args+=( --noZPtRW )
[[ "$noZPtRW"       != "true" ]] && extra_args+=( --zptRWFile "$zptRWFile" )
[[ "$noBeamConstraint" == "true" ]] && extra_args+=( --disableBeamConstraint )
[[ "$isMC" == "true" && "$noBtagSF" != "true" && "$btagEffFile_remote" != "None" ]] && \
  extra_args+=( --btagEffFile "$btagEffFile" )

# Build args; use the full (already-edited) list as-is
bin="./process_HiggsAnalysis_wCutsExtended_NoPairing_v0"
args=( --input "$inputList"
       --output "${tag}.root"
       --isMC "$isMC"
       --sample "$sample"
       --era "$era"
       --muCorr "$muCorr" )

# Only pass --maxEvents if binary supports it
if "$bin" --help 2>&1 | grep -q -- '--maxEvents'; then
  echo "Binary supports --maxEvents; using $maxEvents"
  args+=( --maxEvents "$maxEvents" )
else
  echo "Binary does not support --maxEvents; not passing it."
fi

# Append weights if any
(( ${#weights_arg[@]} )) && args+=( "${weights_arg[@]}" )

# Append extra correction toggles if any
(( ${#extra_args[@]} )) && args+=( "${extra_args[@]}" )

echo
echo "Running:"
printf '%q ' "$bin" "${args[@]}"; echo
echo

set +e
"$bin" "${args[@]}"
rc=$?
set -e
echo "[WRAPPER] binary exit code = $rc"

echo
echo "Produced files:"
ls -lh || true

# Normalize and validate correction mode
muCorr=$(echo "$muCorr" | tr '[:upper:]' '[:lower:]')
case "$muCorr" in
  none) CORR_DIR="NoMuCorr" ;;
  kit)  CORR_DIR="KIT" ;;
  *)    CORR_DIR="Rochester" ;;  # default roch
esac

# Build a label of which corrections are OFF for directory naming
cfg_parts=()

[[ "$noPU"          == "true" ]] && cfg_parts+=(noPU)
[[ "$noScaleSmear"  == "true" ]] && cfg_parts+=(noScaleSmear)
[[ "$noMuonSF"      == "true" ]] && cfg_parts+=(noMuonSF)
[[ "$noJERC"        == "true" ]] && cfg_parts+=(noJERC)
[[ "$noJER"         == "true" ]] && cfg_parts+=(noJER)
[[ "$noJetVeto"     == "true" ]] && cfg_parts+=(noJetVeto)
[[ "$noGolden"      == "true" ]] && cfg_parts+=(noGolden)
[[ "$noFlags"       == "true" ]] && cfg_parts+=(noFlags)
[[ "$noFSR"         == "true" ]] && cfg_parts+=(noFSR)
[[ "$noJetCriteria" == "true" ]] && cfg_parts+=(noJetCriteria)
[[ "$noBeamConstraint" == "true" ]] && cfg_parts+=(noBeamConstraint)

# --- BTag Reweighting flag ---
if [[ "$noBtagSF" == "true" ]]; then
  cfg_parts+=(noBTagRW)
else
  cfg_parts+=(BTagRW)
fi

# --- ZPt Reweighting flag ---
if [[ "$noZPtRW" == "true" ]]; then
  cfg_parts+=(noZPtRW)
else
  cfg_parts+=(ZPtRW)
fi

if (( ${#cfg_parts[@]} == 0 )); then
  CFG_DIR="AllCorrOn"
else
  CFG_DIR="$(IFS=_; echo "${cfg_parts[*]}")"
fi

# Build EOS base: include CORR_DIR and CFG_DIR so runs don't overwrite each other
if [[ "$isMC" == "false" ]]; then
  EOS_BASE="/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/${CORR_DIR}/${CFG_DIR}/DATA"
else
  EOS_BASE="/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/${CORR_DIR}/${CFG_DIR}/${sample}"
fi

EOS_URL="root://eosuser.cern.ch//${EOS_BASE}"

echo "Staging ${tag}.root  ->  ${EOS_URL}/${tag}.root"
xrdfs root://eosuser.cern.ch mkdir -p "${EOS_BASE}" || true
xrdcp -f "${tag}.root" "${EOS_URL}/${tag}.root"
xrdfs root://eosuser.cern.ch stat "${EOS_BASE}/${tag}.root" \
  || { echo "ERROR: xrdcp verification failed"; exit 73; }

# Optional: keep sandbox light
rm -f "${tag}.root" || true

echo "Stage-out done to ${EOS_BASE}/"
echo "=== Done $(date) ==="

