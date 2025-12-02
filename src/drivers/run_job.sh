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
noJetVeto="${12:-false}"
noGolden="${13:-false}"
noFlags="${14:-false}"
noFSR="${15:-false}"
noBtagSF="${16:-false}"
noJetCriteria="${17:-false}"

echo "=== Job $(date) ==="
echo "tag=$tag  sample=$sample  maxEvents=$maxEvents"
echo "inputList (remote) = $inputList_remote"
echo "weights (remote)   = $weights_remote"
echo "Toggles: noPU=$noPU noScaleSmear=$noScaleSmear noMuonSF=$noMuonSF noJERC=$noJERC noJetVeto=$noJetVeto noGolden=$noGolden noFlags=$noFlags noFSR=$noFSR noBtagSF=$noBtagSF noJetCriteria=$noJetCriteria"

# Resolve to local filenames (Condor puts transferred files in CWD under their basename)
inputList="$(basename -- "$inputList_remote")"
echo "inputList (local)  = $inputList"
[[ -s "$inputList" ]] || { echo "ERROR: list '$inputList' missing in $(pwd)"; ls -la; exit 66; }

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
[[ "$noJetVeto"    == "true" ]] && extra_args+=( --noJetVeto )
[[ "$noGolden"     == "true" ]] && extra_args+=( --noGolden )
[[ "$noFlags"      == "true" ]] && extra_args+=( --noFlags )
[[ "$noFSR"        == "true" ]] && extra_args+=( --noFSR )
[[ "$noBtagSF"     == "true" ]] && extra_args+=( --noBtagSF )
[[ "$noJetCriteria" == "true" ]] && extra_args+=( --noJetCriteria )

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

"$bin" "${args[@]}"

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
[[ "$noPU"         == "true" ]] && cfg_parts+=(noPU)
[[ "$noScaleSmear" == "true" ]] && cfg_parts+=(noScaleSmear)
[[ "$noMuonSF"     == "true" ]] && cfg_parts+=(noMuonSF)
[[ "$noJERC"       == "true" ]] && cfg_parts+=(noJERC)
[[ "$noJetVeto"    == "true" ]] && cfg_parts+=(noJetVeto)
[[ "$noGolden"     == "true" ]] && cfg_parts+=(noGolden)
[[ "$noFlags"      == "true" ]] && cfg_parts+=(noFlags)
[[ "$noFSR"        == "true" ]] && cfg_parts+=(noFSR)
[[ "$noBtagSF"     == "true" ]] && cfg_parts+=(noBtagSF)
[[ "$noJetCriteria" == "true" ]] && cfg_parts+=(noJetCriteria)

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

