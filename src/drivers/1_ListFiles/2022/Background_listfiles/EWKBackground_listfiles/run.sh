#!/bin/bash

# Dataset and redirector
DATASET="/EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count lines
TOTAL=$(wc -l < all_files.txt)

# Step 4: Compute lines per file (roughly equal distribution across 6 files)
LINES_PER_FILE=$(( (TOTAL + 1) / 2 ))  # ensures rounding up

# Step 5: Split into 6 list files
split -l $LINES_PER_FILE -d --additional-suffix=.list all_files.txt EWK_2L2J_

# Step 6: Rename split files to use 1-based numbering
for i in {0..5}; do
  mv EWK_2L2J_$(printf "%02d" $i).list EWK_2L2J_$((i + 1)).list
done

# Final output
echo "Done. Files split:"
wc -l EWK_2L2J_*.list

