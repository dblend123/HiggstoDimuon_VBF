#!/usr/bin/env bash
set -euo pipefail

JDL_FILE="condor_submit.jdl"

# Samples to iterate over (must match names used in your JDL + *_input_lists.txt)
SAMPLES=(
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

# Make a one-time backup (won't overwrite if already exists)
[[ -f "${JDL_FILE}.bak" ]] || cp -p "${JDL_FILE}" "${JDL_FILE}.bak"

uncomment_block() {
  local name="$1"

  # Uncomment the exact "sample = NAME" line (allowing arbitrary spaces and an optional leading '#')
  sed -i -E \
    "s|^([[:space:]]*)#?[[:space:]]*sample[[:space:]]*=[[:space:]]*${name}[[:space:]]*$|\\1sample = ${name}|g" \
    "${JDL_FILE}"

  # Uncomment the exact "queue ... from $(path)/NAME_input_lists.txt" line
  # (allow arbitrary spaces and an optional leading '#')
  sed -i -E \
    "s|^([[:space:]]*)#?[[:space:]]*queue[[:space:]]+inputList, tag from[[:space:]]*\\$\\(path\\)/${name}_input_lists\\.txt[[:space:]]*$|\\1queue inputList, tag from \\\$(path)/${name}_input_lists.txt|g" \
    "${JDL_FILE}"
}

recomment_block() {
  local name="$1"

  # Re-comment the now-active sample line
  sed -i -E \
    "s|^([[:space:]]*)sample[[:space:]]*=[[:space:]]*${name}[[:space:]]*$|\\1#sample    = ${name}|g" \
    "${JDL_FILE}"

  # Re-comment the now-active queue line
  sed -i -E \
    "s|^([[:space:]]*)queue[[:space:]]+inputList, tag from[[:space:]]*\\$\\(path\\)/${name}_input_lists\\.txt[[:space:]]*$|\\1#queue inputList, tag from \\\$(path)/${name}_input_lists.txt|g" \
    "${JDL_FILE}"
}

submit_one() {
  local name="$1"
  echo "==> Submitting: ${name}"

  uncomment_block "${name}"

  # Run condor_submit; if it fails, restore comments then propagate the error
  if ! condor_submit "${JDL_FILE}"; then
    echo "!! condor_submit failed for ${name}; re-commenting and aborting." >&2
    recomment_block "${name}"
    return 1
  fi

  recomment_block "${name}"
}

main() {
  # Quick sanity check: ensure the file appears to be the right one
  if ! grep -qE '^\s*#\s*sample\s*=\s*DY\b' "${JDL_FILE}"; then
    echo "Warning: ${JDL_FILE} doesn't look like the expected file (DY block not found)." >&2
  fi

  for s in "${SAMPLES[@]}"; do
    submit_one "${s}"
  done

  echo "All submissions attempted."
}

main "$@"

