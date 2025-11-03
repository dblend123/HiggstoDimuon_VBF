    #!/usr/bin/env bash
    # Robust, Bash-3-compatible generator for input .txt mapping files from .list files.
    # No path args needed (auto-detects list folder).
    #
    # Usage:
    #   ./make_input_txts_from_lists_auto_v2.sh            # outputs next to the script
    #   ./make_input_txts_from_lists_auto_v2.sh <OUT_DIR>  # optional output directory
    #
    # Behavior:
    #   - Looks for .list files (excluding 'all_*.list') in:
    #       1) The script directory
    #       2) A sibling '../Skimmed_listfiles' (if script is in 'Skimmed_input_lists')
    #       3) $LISTS_BASE_DIR (if set and valid)
    #   - Groups by top-level dataset folder under BASE_DIR (e.g., DY, DATA, EWK, ...).
    #   - Writes:
    #       <OUT_DIR>/<DATASET>_input_lists.txt
    #       <OUT_DIR>/all_input_lists.txt
    #     Each line:  <absolute_path_to_list> <basename_without_.list>
    #
    set -euo pipefail

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    OUT_DIR="${1:-${SCRIPT_DIR}}"
    mkdir -p "${OUT_DIR}"

    # Build candidate bases to search for .list files
    CANDIDATES=()
    CANDIDATES+=("${SCRIPT_DIR}")
    if [[ -d "${SCRIPT_DIR}/../Skimmed_listfiles" ]]; then
      CANDIDATES+=("$(cd "${SCRIPT_DIR}/../Skimmed_listfiles" && pwd)")
    fi
    if [[ "${LISTS_BASE_DIR-}" != "" && -d "${LISTS_BASE_DIR}" ]]; then
      CANDIDATES+=("$(cd "${LISTS_BASE_DIR}" && pwd)")
    fi

    BASE_DIR=""
    for cand in "${CANDIDATES[@]}"; do
      if find "$cand" -type f -name '*.list' ! -name 'all_*.list' -print -quit | grep -q . ; then
        BASE_DIR="$cand"
        break
      fi
    done

    if [[ -z "$BASE_DIR" ]]; then
      echo "No .list files found." >&2
      echo "Checked candidates:" >&2
      for cand in "${CANDIDATES[@]}"; do
        echo "  - $cand" >&2
      done
      exit 1
    fi

    echo "[Using LISTS BASE] $BASE_DIR"
    echo "[Writing OUTPUTS ] $OUT_DIR"

    # Collect all .list files (excluding any master 'all_*.list')
    ALL_LISTS=()
    while IFS= read -r -d '' f; do
      ALL_LISTS+=("$f")
    done < <(find "$BASE_DIR" -type f -name '*.list' ! -name 'all_*.list' -print0 | sort -z)

    if [[ ${#ALL_LISTS[@]} -eq 0 ]]; then
      echo "No .list files found under: $BASE_DIR"
      exit 0
    fi

    # Build a unique list of dataset names (top-level dirs under BASE_DIR)
    groups_tmp="$(mktemp)"
    for f in "${ALL_LISTS[@]}"; do
      rel="${f#${BASE_DIR}/}"
      first="${rel%%/*}"
      if [[ "$rel" == "$first" ]]; then
        dataset="_BASE_"
      else
        dataset="$first"
      fi
      echo "$dataset" >> "$groups_tmp"
    done
    sort -u "$groups_tmp" -o "$groups_tmp"

    # Helper: label is the basename without .list
    label() {
      local p="$1"
      local b
      b="$(basename "$p")"
      printf "%s" "${b%.list}"
    }

    # Helper: absolute path with fallbacks
    abspath() {
      local p="$1"
      if command -v readlink >/dev/null 2>&1; then
        if readlink -f / >/dev/null 2>&1; then
          readlink -f "$p" && return 0
        fi
      fi
      if command -v realpath >/dev/null 2>&1; then
        realpath "$p" && return 0
      fi
      if command -v python3 >/dev/null 2>&1; then
        python3 - "$p" <<'PY'
import os, sys
print(os.path.abspath(sys.argv[1]))
PY
        return 0
      fi
      # Fallback: echo as-is
      printf "%s
" "$p"
    }

    master="${OUT_DIR}/all_input_lists.txt"
    : > "$master"

    while IFS= read -r dataset; do
      out="${OUT_DIR}/${dataset}_input_lists.txt"
      : > "$out"
      for f in "${ALL_LISTS[@]}"; do
        rel="${f#${BASE_DIR}/}"
        first="${rel%%/*}"
        belongs="$([[ "$rel" == "$first" ]] && echo "_BASE_" || echo "$first")"
        if [[ "$belongs" == "$dataset" ]]; then
          # Skip dataset aggregate file if present (e.g., <BASE_DIR>/DY.list)
          if [[ "$f" == "${BASE_DIR}/${dataset}.list" ]]; then
            continue
          fi
          abs="$(abspath "$f")"
          printf "%s %s\n" "$abs" "$(label "$f")" | tee -a "$out" >> "$master"
        fi
      done
      echo "[Wrote] $out ($(wc -l < "$out") lines)"
    done < "$groups_tmp"

    echo "[Wrote] $master ($(wc -l < "$master") lines)"
    rm -f "$groups_tmp"
