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

echo "=== Job $(date) ==="
echo "tag=$tag  sample=$sample  maxEvents=$maxEvents"
echo "inputList (remote) = $inputList_remote"
echo "weights (remote)   = $weights_remote"

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

# Build args; use the full (already-edited) list as-is
bin="./process_HiggsAnalysis_wCutsExtended_NoPairing_v0"
args=( --input "$inputList" --output "${tag}.root" --isMC "$isMC" --sample "$sample" --era "$era" --muCorr "$muCorr" )

# Only pass --maxEvents if binary supports it
if "$bin" --help 2>&1 | grep -q -- '--maxEvents'; then
  echo "Binary supports --maxEvents; using $maxEvents"
  args+=( --maxEvents "$maxEvents" )
else
  echo "Binary does not support --maxEvents; not passing it."
fi

# Append weights if any
(( ${#weights_arg[@]} )) && args+=( "${weights_arg[@]}" )

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

# Build EOS base: include CORR_DIR even for DATA so runs don't overwrite each other
if [[ "$isMC" == "false" ]]; then
  EOS_BASE="/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/${CORR_DIR}/DATA"
else
  EOS_BASE="/eos/user/d/dblend/Analyzed_SkimmedFiles/${era}/${CORR_DIR}/${sample}"
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

