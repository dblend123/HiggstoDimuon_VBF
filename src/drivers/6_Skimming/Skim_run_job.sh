#!/usr/bin/env bash
# Skim_run_job.sh — HTCondor wrapper for MyAnaSkim
# Expects arguments: <inputList> <tag> <sample>
# Uses env vars: MAX_EVENTS (default: -1), SKIM_BIN (default: my_skim), EOS_BASE (optional)

set -euo pipefail

# ---------------------------
# Parse args / defaults
# ---------------------------
INPUT_LIST="${1:-}"
TAG="${2:-}"
SAMPLE="${3:-}"
if [[ -z "${INPUT_LIST}" || -z "${TAG}" || -z "${SAMPLE}" ]]; then
  echo "Usage: $0 <inputList> <tag> <sample>" >&2
  exit 2
fi
ERA="${4:-}"

# my_skim takes: <filelist> <max_events>
MAX_EVENTS="${MAX_EVENTS:- -1}"
# normalize MAX_EVENTS to avoid stray leading spaces
MAX_EVENTS="${MAX_EVENTS//[[:space:]]/}"
if [[ -z "${MAX_EVENTS}" ]]; then MAX_EVENTS="-1"; fi

SKIM_BIN="${SKIM_BIN:-my_skim}"

start_ts="$(date)"
host="$(hostname -f || hostname)"
echo "== Job started: ${start_ts}"
echo "Host: ${host}"
echo "Args: inputList='${INPUT_LIST}'  tag='${TAG}'  sample='${SAMPLE}'"
echo "ERA (4th arg): '${ERA}'"
echo "MAX_EVENTS:  ${MAX_EVENTS}"

# Normalize the list file to the transferred local copy
LIST_BASE="$(basename "${INPUT_LIST}")"
LIST_LOCAL="${PWD}/${LIST_BASE}"
echo "PWD: ${PWD}"
echo "Input list arg: '${INPUT_LIST}' -> using local: '${LIST_LOCAL}'"

if [[ -f "${LIST_LOCAL}" ]]; then
  : # OK, use the transferred list
elif [[ -f "${INPUT_LIST}" ]]; then
  # Rare case: worker can see AFS; copy for consistency
  cp -v "${INPUT_LIST}" "${LIST_LOCAL}"
else
  echo "FATAL: neither '${LIST_LOCAL}' (transferred) nor '${INPUT_LIST}' (original) exists on the worker." >&2
  echo "Listing working dir:"; ls -l
  exit 2
fi

# ---------------------------
# Auth (optional): print proxy path if present
# ---------------------------
if [[ -n "${X509_USER_PROXY:-}" && -f "${X509_USER_PROXY}" ]]; then
  echo "Using Condor proxy: ${X509_USER_PROXY}"
fi

# ---------------------------
# Environment: ensure ROOT from LCG view
# ---------------------------
LCG_VIEW="${LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh}"
set +u
if [[ -f "${LCG_VIEW}" ]]; then
  # shellcheck disable=SC1090
  source "${LCG_VIEW}"
elif [[ -f /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh ]]; then
  # Fallback (older) — not ideal but better than system ROOT
  # shellcheck disable=SC1090
  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
elif [[ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]]; then
  # Last-resort CMSSW runtime (only if no LCG view available)
  export SCRAM_ARCH="${SCRAM_ARCH:-el9_amd64_gcc12}"
  # shellcheck disable=SC1091
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  scramv1 project CMSSW CMSSW_13_0_14
  pushd CMSSW_13_0_14/src >/dev/null
  eval "$(scramv1 runtime -sh)"
  popd >/dev/null
else
  echo "FATAL: No ROOT/LCG environment available on worker." >&2
  exit 99
fi
set -u

echo "ROOT bin:   $(command -v root || echo 'not found')"
echo "ROOT ver:   $(root-config --version || echo 'unknown')"
echo "libdir:     $(root-config --libdir || echo 'unknown')"

# ---------------------------
# Input list sanity
# ---------------------------
echo "Using input list: ${LIST_LOCAL}"
if [[ ! -s "${LIST_LOCAL}" ]]; then
  echo "FATAL: Input list missing or empty: ${LIST_LOCAL}" >&2  
  exit 3
fi
echo "Input list has $(wc -l < "${LIST_LOCAL}") lines"
head -n 3 "${LIST_LOCAL}" | sed 's|$| ...|' || true

# ---------------------------
# Prepare skims dir
# ---------------------------
mkdir -p ./skims

# ---------------------------
# Ensure correct binary & show fingerprint
# ---------------------------
if [[ ! -x "./${SKIM_BIN}" ]]; then
  echo "FATAL: Skimmer binary './${SKIM_BIN}' not found or not executable." >&2
  exit 10
fi
echo "Using skimmer: ${SKIM_BIN}"
command -v md5sum >/dev/null && md5sum "./${SKIM_BIN}" || stat "./${SKIM_BIN}"

# Capture ldd once, print it, then test it
LDD_OUT="$(ldd "./${SKIM_BIN}" 2>&1)"
if grep -q '/usr/lib64/root' <<<"${LDD_OUT}"; then
  echo "FATAL: ${SKIM_BIN} is linked against system ROOT (/usr/lib64/root)." >&2
  echo "${LDD_OUT}" >&2
  exit 88
fi

# Optional: show RPATH/RUNPATH for extra confidence
if command -v readelf >/dev/null; then
  readelf -d "./${SKIM_BIN}" | grep -E 'RPATH|RUNPATH' || true
fi

# ---------------------------
# Clean TAG to avoid odd EOS dirs
# ---------------------------
TAG_CLEAN="${TAG//[[:space:]]/_}"
TAG_CLEAN="${TAG_CLEAN//[^A-Za-z0-9._-]/_}"
echo "Using cleaned TAG: ${TAG_CLEAN}"

# ---------------------------
# Run the skimmer
# ---------------------------
set +e
echo "Running: ./${SKIM_BIN} '${LIST_LOCAL}' '${MAX_EVENTS}'"

"./${SKIM_BIN}" "${LIST_LOCAL}" "${MAX_EVENTS}" > skim.stdout.log 2> skim.stderr.log
SKIM_RC=$?

# Keep condor logs small:
if [[ ${SKIM_RC} -ne 0 ]]; then
  echo "my_skim FAILED (exit ${SKIM_RC})" >&2
  echo "---- last 200 lines of skim.stderr.log ----" >&2
  tail -n 200 skim.stderr.log >&2 || true
else
  # optional: show a tiny amount of normal output
  echo "---- last 50 lines of skim.stdout.log ----"
  tail -n 50 skim.stdout.log || true
fi
set -e

# ---------------------------
# Stage-out to EOS (and remove local files)
# ---------------------------
# Allow override via EOS_BASE; otherwise default to your user path
DEFAULT_EOS_BASE="root://eosuser.cern.ch//eos/experiment/higgstodimuon-vbf/${ERA}/${SAMPLE}/${TAG_CLEAN}"
EOS_BASE="${EOS_BASE:-${DEFAULT_EOS_BASE}}"

echo "Staging to: ${EOS_BASE}"
# Make dir on EOS (ignore if exists)
xrdfs root://eosuser.cern.ch mkdir -p "/eos/experiment/higgstodimuon-vbf/${ERA}/${SAMPLE}/${TAG_CLEAN}" || true

# Gather outputs
mapfile -t OUTS < <(find ./skims -type f -name 'Skimmed_*.root' 2>/dev/null | LC_ALL=C sort || true)
echo "Found ${#OUTS[@]} skim files."

FAILS=0
if (( ${#OUTS[@]} > 0 )); then
  for f in "${OUTS[@]}"; do
    base="$(basename "$f")"
    echo "xrdcp -> ${EOS_BASE}/${base}"
    if xrdcp -f "$f" "${EOS_BASE}/${base}"; then
      rm -f "$f" || true
    else
      echo "ERROR: xrdcp failed for ${f}" >&2
      ((FAILS++))
    fi
  done
  # remove empty dir
  rm -rf ./skims || true
else
  echo "WARNING: no skim outputs found"
fi

# ---------------------------
# Exit status logic
# ---------------------------
if (( SKIM_RC != 0 )); then
  echo "FATAL: my_skim failed with code ${SKIM_RC}" >&2
  exit "${SKIM_RC}"
fi
if (( ${#OUTS[@]} == 0 )); then
  echo "FATAL: skimmer produced no outputs." >&2
  exit 5
fi
if (( FAILS > 0 )); then
  echo "FATAL: ${FAILS} file(s) failed to stage to EOS." >&2
  exit 4
fi

echo "Done staging."
echo "== Job finished: $(date)"
exit 0

