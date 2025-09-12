#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./batch_condor_submit_v3.sh <ERA> [SAMPLE1 SAMPLE2 ...]
#
# Prefers .jdl; falls back to .txt
JDL_FILE="Skim_condor_submit.jdl"
if [[ ! -f "$JDL_FILE" ]]; then
  JDL_FILE="Skim_condor_submit.txt"
fi
BASE_FILE="${JDL_FILE}.base"

DEFAULT_SAMPLES=(
  DY
  EWK
  TT
  WWW
  WWZ
  WZZ
  WZto2L2Q
  WZto3LNu
  ZZZ
  ZZto2L2Nu
  ZZto2L2Q
  ZZto4L
  VBF
)

usage(){ echo "Usage: $0 <ERA> [SAMPLE1 SAMPLE2 ...]"; exit 2; }
[[ $# -lt 1 ]] && usage
ERA="$1"; shift
SAMPLES=("$@"); [[ ${#SAMPLES[@]} -eq 0 ]] && SAMPLES=("${DEFAULT_SAMPLES[@]}")

# Create a pristine base if missing
if [[ ! -f "$BASE_FILE" ]]; then
  cp -p "$JDL_FILE" "$BASE_FILE"
fi

restore_base(){ cp -p "$BASE_FILE" "$JDL_FILE"; }

set_era(){
  local era_val="$1"
  sed -i -E 's/^([[:space:]]*)#?[[:space:]]*era[[:space:]]*=[[:space:]]*[0-9]+([[:space:]]*)$/\1# era       = X\2/g' "$JDL_FILE"
  if grep -qE '^\s*#\s*era\s*=\s*X\s*$' "$JDL_FILE"; then
    sed -i -E "0,/^([[:space:]]*)#\s*era\s*=\s*X\s*$/s//\1era       = ${era_val}/" "$JDL_FILE"
  else
    printf '\nera       = %s\n' "$era_val" >> "$JDL_FILE"
  fi
}

fix_arguments_line(){
  if grep -qE '^\s*arguments\s*=' "$JDL_FILE"; then
    sed -i -E 's/^\s*arguments\s*=.*/arguments = $(inputList) $(tag) $(sample) $(era)/' "$JDL_FILE"
  else
    printf '\narguments = $(inputList) $(tag) $(sample) $(era)\n' >> "$JDL_FILE"
  fi
}

comment_all_queues(){
  sed -i -E 's/^([[:space:]]*)(queue[[:space:]].*)/\1# \2/g' "$JDL_FILE"
}

comment_all_samples(){
  sed -i -E 's/^([[:space:]]*)(sample[[:space:]]*=[[:space:]]*.+)/\1# \2/g' "$JDL_FILE"
}

activate_block(){
  local name="$1"
  awk -v name="$name" '
    BEGIN { inblk=0; donequeue=0 }
    {
      # Un-comment the exact "sample = NAME" line (even if already commented by sed)
      if ($0 ~ "^[[:space:]]*#?[[:space:]]*sample[[:space:]]*=[[:space:]]*" name "[[:space:]]*$") {
        line = $0
        sub(/^[ \t]*#?[ \t]*/, "", line)
        print line
        inblk=1; donequeue=0
        next
      }
      if (inblk==1) {
        if ($0 ~ /^[[:space:]]*$/) { inblk=0 }  # blank line ends block
        else if (donequeue==0 && $0 ~ "^[[:space:]]*#?[[:space:]]*queue[[:space:]]+") {
          line = $0
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

submit_one(){
  local name="$1"
  echo "==> Submitting: ${name} (ERA=${ERA})"
  restore_base
  set_era "$ERA"
  fix_arguments_line
  comment_all_queues
  comment_all_samples
  activate_block "$name"
  condor_submit "$JDL_FILE"
}

main(){
  for s in "${SAMPLES[@]}"; do
    submit_one "$s"
  done
  echo "All submissions attempted."
}

main "$@"
