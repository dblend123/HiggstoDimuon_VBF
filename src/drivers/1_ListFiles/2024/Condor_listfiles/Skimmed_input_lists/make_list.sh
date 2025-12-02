#!/bin/bash
# Make a .list file that contains the full path + short name
# for every .list file inside a given directory.
#
# Usage:
#   ./make_list_of_lists.sh /path/to/Skimmed_listfiles/DATA output_MASTER.list

LIST_DIR="$1"
OUTPUT_FILE="$2"

if [[ -z "$LIST_DIR" || -z "$OUTPUT_FILE" ]]; then
    echo "Usage: $0 <directory_with_lists> <output_list>"
    exit 1
fi

if [[ ! -d "$LIST_DIR" ]]; then
    echo "Error: directory '$LIST_DIR' not found."
    exit 1
fi

echo "Creating master list: $OUTPUT_FILE"
> "$OUTPUT_FILE"

shopt -s nullglob
for f in "$LIST_DIR"/*.list; do
    base="$(basename "$f" .list)"
    fullpath="$(readlink -f "$f")"
    echo "$fullpath $base" >> "$OUTPUT_FILE"
done
shopt -u nullglob

echo "Done. Wrote $(wc -l < "$OUTPUT_FILE") entries."

