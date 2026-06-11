#!/usr/bin/env bash
# Example usage:
#   ./batch_condor_submit.sh 12022,22022,12023,22023
#   ./batch_condor_submit.sh 12022,22022 DY EWK DATA
#
# Optional environment overrides before running:
#   NOZPTRW=true ./batch_condor_submit.sh 12022,22022
#   NOZPTRW=false NOJERC=true NOBTAGSF=true ./batch_condor_submit.sh 2024 DY
#
# Notes:
# - The script edits/submits JDL_FILE, restoring from BASE_FILE before each job.
# - Default behavior below sets noZPtRW=true unless overridden by env.
# - Adjust JDL_FILE if your real submit file is not condor_submit.jdl.

set -euo pipefail

[[ $# -lt 1 ]] && {
  echo "Usage: $0 <ERA1,ERA2,...> [SAMPLE1 SAMPLE2 ...]"
  exit 2
}

# -----------------------------
# Files
# -----------------------------
JDL_FILE="condor_submit.jdl"
BASE_FILE="${JDL_FILE}.base"

# If you actually use condor_submit.txt instead, change the line above to:
# JDL_FILE="condor_submit.txt"

# -----------------------------
# Parse eras
# -----------------------------
IFS=',' read -ra ERAS <<< "$1"
shift

# -----------------------------
# Samples
# -----------------------------
DEFAULT_SAMPLES=(
  DY EWK TT WWW WWZ WZZ WZto2L2Q WZto3LNu
  ZZZ ZZto2L2Nu ZZto2L2Q ZZto4L VBF DATA
)

SAMPLES=("$@")
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  SAMPLES=("${DEFAULT_SAMPLES[@]}")
fi

# -----------------------------
# Toggle defaults
# Override from environment if desired
# -----------------------------
NOJERC="${NOJERC:-true}"
NOJER="${NOJER:-true}"
NOJETCRITERIA="${NOJETCRITERIA:-false}"
NOJETVETO="${NOJETVETO:-false}"
NOBTAGSF="${NOBTAGSF:-true}"
NOSCALESMEAR="${NOSCALESMEAR:-true}"
NOBEAMCONSTRAINT="${NOBEAMCONSTRAINT:-false}"
NOZPTRW="${NOZPTRW:-true}"       
NOMUONSF="${NOMUONSF:-false}"
NOFSR="${NOFSR:-false}"
NOPU="${NOPU:-false}"
NOGOLDEN="${NOGOLDEN:-false}"
NOFLAGS="${NOFLAGS:-false}"
MUCORR="${MUCORR:-kit}"

PRINTED_CONFIG=false

# -----------------------------
# Base-file handling
# -----------------------------
if [[ ! -f "$JDL_FILE" ]]; then
  echo "[ERROR] JDL file not found: $JDL_FILE"
  exit 1
fi

if [[ ! -f "$BASE_FILE" ]]; then
  cp -p "$JDL_FILE" "$BASE_FILE"
  echo "[INFO] Created base file: $BASE_FILE"
fi

restore_base() {
  cp -p "$BASE_FILE" "$JDL_FILE"
}

# -----------------------------
# Generic setter for "key = value"
# -----------------------------
set_key_value() {
  local key="$1"
  local val="$2"

  if grep -qE "^[[:space:]]*#?[[:space:]]*${key}[[:space:]]*=" "$JDL_FILE"; then
    sed -i -E "s|^[[:space:]]*#?[[:space:]]*${key}[[:space:]]*=.*|${key}             = ${val}|" "$JDL_FILE"
  else
    printf '\n%s             = %s\n' "$key" "$val" >> "$JDL_FILE"
  fi
}

set_era() {
  local era_val="$1"
  set_key_value "era" "$era_val"
}

set_muCorr() {
  local mode="$1"
  set_key_value "muCorr" "$mode"
}

set_toggle() {
  local key="$1"
  local val="$2"
  set_key_value "$key" "$val"
}

# -----------------------------
# Disable all sample/queue blocks
# -----------------------------
comment_all_queues() {
  sed -i -E 's/^([[:space:]]*)(queue[[:space:]].*)/\1# \2/g' "$JDL_FILE"
}

comment_all_samples() {
  sed -i -E 's/^([[:space:]]*)(sample[[:space:]]*=[[:space:]]*.+)/\1# \2/g' "$JDL_FILE"
}

# -----------------------------
# Activate one sample block
# Assumes each block has:
#   # sample = NAME
#   # queue ...
# or
#   sample = NAME
#   queue ...
# -----------------------------
activate_block() {
  local name="$1"

  awk -v name="$name" '
  BEGIN { inblk=0; donequeue=0 }
  {
    if ($0 ~ "^[[:space:]]*#?[[:space:]]*sample[[:space:]]*=[[:space:]]*" name "[[:space:]]*$") {
      line=$0
      sub(/^[ \t]*#?[ \t]*/, "", line)
      print line
      inblk=1
      donequeue=0
      next
    }

    if (inblk==1) {
      if ($0 ~ /^[[:space:]]*$/) {
        inblk=0
      }
      else if (donequeue==0 && $0 ~ "^[[:space:]]*#?[[:space:]]*queue[[:space:]]+") {
        line=$0
        sub(/^[ \t]*#?[ \t]*/, "", line)
        print line
        donequeue=1
        next
      }
    }

    print
  }
  ' "$JDL_FILE" > "${JDL_FILE}.tmp" && mv "${JDL_FILE}.tmp" "$JDL_FILE"
}

# -----------------------------
# Sanity check sample exists in JDL/base
# -----------------------------
sample_exists() {
  local name="$1"
  grep -qE "^[[:space:]]*#?[[:space:]]*sample[[:space:]]*=[[:space:]]*${name}[[:space:]]*$" "$BASE_FILE"
}

# -----------------------------
# Submit one era/sample combo
# -----------------------------
submit_one() {
  local era="$1"
  local sample="$2"

  restore_base

  set_era "$era"
  set_muCorr "$MUCORR"

  # Toggles
  set_toggle "noJERC" "$NOJERC"
  set_toggle "noJER" "$NOJER"
  set_toggle "noJetCriteria" "$NOJETCRITERIA"
  set_toggle "noJetVeto" "$NOJETVETO"
  set_toggle "noBtagSF" "$NOBTAGSF"
  set_toggle "noScaleSmear" "$NOSCALESMEAR"
  set_toggle "noBeamConstraint" "$NOBEAMCONSTRAINT"
  set_toggle "noZPtRW" "$NOZPTRW"
  set_toggle "noMuonSF" "$NOMUONSF"
  set_toggle "noFSR" "$NOFSR"
  set_toggle "noPU" "$NOPU"
  set_toggle "noGolden" "$NOGOLDEN"
  set_toggle "noFlags" "$NOFLAGS"

  comment_all_queues
  comment_all_samples
  activate_block "$sample"

  if [[ "$PRINTED_CONFIG" = false ]]; then
    echo "[INFO] Active config:"
    grep -nE '^[[:space:]]*(muCorr|noJERC|noJER|noJetCriteria|noJetVeto|noBtagSF|noScaleSmear|noBeamConstraint|noZPtRW|noMuonSF|noFSR|noPU|noGolden|noFlags)[[:space:]]*=' "$JDL_FILE"
    PRINTED_CONFIG=true
  fi

  condor_submit "$JDL_FILE"
}

# -----------------------------
# Main loop
# -----------------------------
main() {
  echo "[INFO] JDL_FILE   = $JDL_FILE"
  echo "[INFO] BASE_FILE  = $BASE_FILE"
  echo "[INFO] ERAs       = ${ERAS[*]}"
  echo "[INFO] Samples    = ${SAMPLES[*]}"
  echo "[INFO] muCorr     = $MUCORR"
  echo "[INFO] noZPtRW    = $NOZPTRW"
  echo "[INFO] noJERC     = $NOJERC"
  echo "[INFO] noJER      = $NOJER"
  echo "[INFO] noJetCriteria = $NOJETCRITERIA"
  echo "[INFO] noJetVeto  = $NOJETVETO"
  echo "[INFO] noBtagSF   = $NOBTAGSF"
  echo "[INFO] noScaleSmear = $NOSCALESMEAR"
  echo "[INFO] noBeamConstraint = $NOBEAMCONSTRAINT"
  echo "[INFO] noMuonSF   = $NOMUONSF"
  echo "[INFO] noFSR      = $NOFSR"
  echo "[INFO] noPU       = $NOPU"
  echo "[INFO] noGolden   = $NOGOLDEN"
  echo "[INFO] noFlags    = $NOFLAGS"

  for era in "${ERAS[@]}"; do
    for sample in "${SAMPLES[@]}"; do
      submit_one "$era" "$sample"
    done
  done

  echo
  echo "[DONE] All submissions attempted."
}

main
