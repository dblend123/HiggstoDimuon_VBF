    #!/usr/bin/env bash
    # Create .list files of ROOT files per dataset directory AND per immediate subdirectory,
    # plus a master list.
    #
    # Usage:
    #   ./make_root_lists_per_subdir.sh [/eos/user/d/dblend/SkimmedFiles/22022] [OUTPUT_DIR]
    #
    # Output layout (example):
    #   OUTPUT_DIR/
    #     all_22022.list
    #     DY.list                         # all ROOTs under DY/*
    #     DY/                             # per-subdirectory lists
    #       DY_000MCRun22022.list
    #       DY_001MCRun22022.list
    #       ...
    #     DATA.list
    #     DATA/
    #       DataRun2022C.list
    #       ...
    #
    set -euo pipefail

    BASE_DIR="${1:-/eos/user/d/dblend/SkimmedFiles/22022}"
    OUT_DIR="${2:-.}"

    if [[ ! -d "${BASE_DIR}" ]]; then
      echo "Base directory not found: ${BASE_DIR}" >&2
      exit 1
    fi

    mkdir -p "${OUT_DIR}"

    base_name="$(basename "${BASE_DIR}")"
    master_list="${OUT_DIR}/all_${base_name}.list"
    : > "${master_list}"

    echo "Scanning base: ${BASE_DIR}"

    shopt -s nullglob dotglob

    # Iterate top-level dataset directories
    for dataset_dir in "${BASE_DIR}"/*/ ; do
      [[ -d "$dataset_dir" ]] || continue
      dataset="$(basename "$dataset_dir")"
      dataset_out_dir="${OUT_DIR}/${dataset}"
      mkdir -p "${dataset_out_dir}"

      dataset_list="${OUT_DIR}/${dataset}.list"
      : > "${dataset_list}"
      echo "[Dataset] ${dataset}"

      # Include any ROOT files directly under dataset_dir (rare but possible)
      tmp_dataset="$(mktemp)"
      find "$dataset_dir" -maxdepth 1 -type f -name '*.root' -print0         | xargs -0 -I{} printf '%s
' "{}" > "${tmp_dataset}" || true
      if [[ -s "${tmp_dataset}" ]]; then
        sort -u "${tmp_dataset}" -o "${tmp_dataset}"
        cat "${tmp_dataset}" >> "${dataset_list}"
        cat "${tmp_dataset}" >> "${master_list}"
      fi
      rm -f "${tmp_dataset}"

      # Per-immediate-subdirectory lists under this dataset
      for subdir in "${dataset_dir}"*/ ; do
        [[ -d "$subdir" ]] || continue
        subname="$(basename "$subdir")"
        sub_list="${dataset_out_dir}/${subname}.list"
        : > "${sub_list}"

        tmp_sub="$(mktemp)"
        # Recurse inside each immediate child subdir
        find "$subdir" -type f -name '*.root' -print0           | xargs -0 -I{} printf '%s
' "{}" > "${tmp_sub}" || true

        if [[ -s "${tmp_sub}" ]]; then
          sort -u "${tmp_sub}" -o "${tmp_sub}"
          cat "${tmp_sub}" >> "${sub_list}"
          cat "${tmp_sub}" >> "${dataset_list}"
          cat "${tmp_sub}" >> "${master_list}"
          count="$(wc -l < "${tmp_sub}")"
          echo "  - ${dataset}/${subname}: ${count} file(s)"
        else
          echo "  - ${dataset}/${subname}: 0 file(s)"
        fi
        rm -f "${tmp_sub}"
      done

      # Final count for dataset
      if [[ -s "${dataset_list}" ]]; then
        dcount="$(wc -l < "${dataset_list}")"
        echo "    => ${dataset}.list: ${dcount} file(s)"
      else
        echo "    => ${dataset}.list: 0 file(s)"
      fi
    done

    total="$(wc -l < "${master_list}" || echo 0)"
    echo "Done. Total ROOT files across all datasets: ${total}"
    echo "Master list: ${master_list}"

    # Tips:
    # To use EOS client instead of plain POSIX paths, you can replace 'find' with:
    #   eos find --name '*.root' <DIR>
    # To convert local EOS paths to xrootd URLs:
    #   sed -E 's|^/eos/user/([a-z])/([^/]+)/|root://eosuser.cern.ch//eos/user/\1/\2/|' INPUT.list > OUTPUT.xrd.list
