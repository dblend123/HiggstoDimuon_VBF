#!/bin/bash
# Create .list files for all subdirectories of a given base directory.
# Each .list will contain full paths to all .root files in that subdir.
#
# Usage:
#   ./make_eos_lists.sh /eos/user/d/dblend/SkimmedFiles/2024/DATA_2024_Muon0
#   ./make_eos_lists.sh /eos/user/d/dblend/SkimmedFiles/2024/DATA_2024_Muon1
#
# Optional: specify an output directory for the .list files:
#   ./make_eos_lists.sh <base_dir> <output_dir>
#
# Then you'll get, for example:
#   <output_dir>/Muon0C_000DataRun2024.list
#   <output_dir>/Muon0C_001DataRun2024.list
#   ...

set -euo pipefail

BASE_DIR="${1:-}"
OUT_DIR="${2:-}"

if [[ -z "${BASE_DIR}" ]]; then
    echo "Usage: $0 <base_dir> [output_dir]" >&2
    exit 1
fi

if [[ ! -d "${BASE_DIR}" ]]; then
    echo "Error: base_dir '${BASE_DIR}' does not exist or is not a directory." >&2
    exit 1
fi

# If no output directory is given, put the .list files in BASE_DIR
if [[ -z "${OUT_DIR}" ]]; then
    OUT_DIR="${BASE_DIR}"
fi

mkdir -p "${OUT_DIR}"

echo "Base directory : ${BASE_DIR}"
echo "Output dir     : ${OUT_DIR}"
echo

shopt -s nullglob
subdirs=( "${BASE_DIR}"/*/ )
shopt -u nullglob

if (( ${#subdirs[@]} == 0 )); then
    echo "No subdirectories found in ${BASE_DIR}."
    exit 0
fi

for subdir in "${subdirs[@]}"; do
    # Strip trailing slash
    subdir="${subdir%/}"
    sample_name="$(basename "${subdir}")"
    list_file="${OUT_DIR}/${sample_name}.list"

    echo ">>> Building list for ${sample_name}"
    echo "    Subdir : ${subdir}"
    echo "    Output : ${list_file}"

    # Find all .root files under this subdir, sort, and write to list
    # (If you want non-recursive, add: -maxdepth 1)
    find "${subdir}" -type f -name "*.root" | sort > "${list_file}"

    n_files=$(wc -l < "${list_file}" || echo 0)
    echo "    Added ${n_files} .root files"
    echo
done

echo "Done."

