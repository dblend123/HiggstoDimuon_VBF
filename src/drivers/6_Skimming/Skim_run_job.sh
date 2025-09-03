#!/usr/bin/env bash
# Wrapper for Condor with JDL:
#   arguments = $(inputList) $(tag) $(sample)
#   transfer_output_files   = $(tag).root
#   transfer_output_remaps  = "$(tag).root=/eos/.../$(sample)/$(tag).root"

set -euo pipefail

INPUT_LIST_RAW="${1:?need inputList}"
TAG="${2:?need tag}"
SAMPLE="${3:?need sample}"

# Default max events (since JDL doesn't pass it). You can override via:
#   condor_submit ... -append 'environment = "MAX_EVENTS=50000"'
MAX_EVENTS="${MAX_EVENTS:-${maxEvents:- -1}}"

echo "== Job started: $(date)"
echo "Host: $(hostname -f)"
echo "Args: inputList='${INPUT_LIST_RAW}'  tag='${TAG}'  sample='${SAMPLE}'"
echo "MAX_EVENTS: ${MAX_EVENTS}"

# ---------------------------
# Proxy
# ---------------------------
if [[ -n "${X509_USER_PROXY:-}" && -f "${X509_USER_PROXY}" ]]; then
  echo "Using Condor proxy: ${X509_USER_PROXY}"
elif [[ -f ./x509up.cert ]]; then
  export X509_USER_PROXY="$PWD/x509up.cert"
  echo "Using shipped proxy: ${X509_USER_PROXY}"
else
  echo "WARNING: No X509 proxy found; remote reads may fail."
fi

# ---------------------------
# Environment (ROOT from CVMFS) — relax nounset while sourcing
# ---------------------------
set +u
if [[ -f /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh ]]; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
elif [[ -f /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh ]]; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
elif [[ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]]; then
  export SCRAM_ARCH="${SCRAM_ARCH:-el9_amd64_gcc12}"
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  scramv1 project CMSSW CMSSW_13_0_14
  pushd CMSSW_13_0_14/src >/dev/null
  eval "$(scramv1 runtime -sh)"
  popd >/dev/null
else
  echo "FATAL: No ROOT/CMSSW environment available on worker." >&2
  # create placeholder so transfer doesn't fail
  root -l -q -b -e "TFile f(\"${TAG}.root\",\"RECREATE\"); f.Close();" || true
  exit 99
fi
set -u

echo "ROOT: $(which root)"
echo "ROOT version: $(root-config --version)"
echo "GCC: $(c++ --version | head -1)"


# ---------------------------
# Resolve input list (Condor usually transfers by basename)
# ---------------------------
INPUT_LIST="${INPUT_LIST_RAW}"
if [[ ! -f "${INPUT_LIST}" ]]; then
  BASE="$(basename "${INPUT_LIST_RAW}")"
  if [[ -f "${BASE}" ]]; then
    INPUT_LIST="${BASE}"
  else
    echo "FATAL: input list not found: '${INPUT_LIST_RAW}' or '${BASE}'" >&2
    # ensure condor can still transfer something
    root -l -q -b -e "TFile f(\"${TAG}.root\",\"RECREATE\"); f.Close();" || true
    exit 2
  fi
fi
echo "Using input list: ${INPUT_LIST}"


echo "Input list has $(grep -cv '^[[:space:]]*$' "${INPUT_LIST}") lines"
head -n 3 "${INPUT_LIST}" || true


# ---------------------------
# Binary presence & perms
# ---------------------------
if [[ ! -f ./my_skim ]]; then
  echo "FATAL: ./my_skim missing." >&2
  root -l -q -b -e "TFile f(\"${TAG}.root\",\"RECREATE\"); f.Close();" || true
  exit 3
fi
chmod +x ./my_skim || true

# ---------------------------
# Run the skimmer
# ---------------------------
echo "Running: ./my_skim '${INPUT_LIST}' '${MAX_EVENTS}'"
set +e
./my_skim "${INPUT_LIST}" "${MAX_EVENTS}"
SKIM_RC=$?
set -e
echo "my_skim exit code: ${SKIM_RC}"

# ---------------------------
# Stage out: copy each skim file to EOS, then delete local copies
# ---------------------------
EOS_BASE="${EOS_BASE:-root://eosuser.cern.ch//eos/user/d/dblend/Outputs/${SAMPLE}/${TAG}}"
echo "Staging to: ${EOS_BASE}"
xrdfs root://eosuser.cern.ch mkdir -p "/eos/user/d/dblend/Outputs/${SAMPLE}/${TAG}" || true

mapfile -t OUTS < <(find ./skims -type f -name 'Skimmed_*.root' 2>/dev/null | sort || true)
echo "Found ${#OUTS[@]} skim files."

FAILS=0
if (( ${#OUTS[@]} > 0 )); then
  for f in "${OUTS[@]}"; do
    base=$(basename "$f")
    echo "xrdcp -> ${EOS_BASE}/${base}"
    if xrdcp -f "$f" "${EOS_BASE}/${base}"; then
      rm -f "$f" || true     # remove local file after successful copy
    else
      echo "ERROR: xrdcp failed for ${f}" >&2
      ((FAILS++))
    fi
  done
else
  echo "WARNING: no skim outputs found"
fi

# Remove (now empty) directory to avoid leaving anything behind
rm -rf ./skims || true

# If any copies failed, exit non-zero so you notice in logs/monitoring
if (( FAILS > 0 )); then
  echo "FATAL: ${FAILS} file(s) failed to stage to EOS." >&2
  exit 4
fi

echo "Done staging."
echo "== Job finished: $(date)"

# Prefer success so Condor transfers the file; inspect logs for SKIM_RC
exit 0

