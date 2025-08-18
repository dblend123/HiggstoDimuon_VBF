#!/bin/bash

# Dataset and redirector
DATASET="/SingleMuon/Run2022C-22Sep2023-v1/NANOAOD"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count total lines
TOTAL=$(wc -l < all_files.txt)
THIRD=$(( (TOTAL + 2) / 3 ))  # round up for 3 parts

# Step 4: Split into 3 files
head -n $THIRD all_files.txt > SingleMuon_Run2022C_1.list
tail -n +$((THIRD + 1)) all_files.txt | head -n $THIRD > SingleMuon_Run2022C_2.list
tail -n +$((2 * THIRD + 1)) all_files.txt > SingleMuon_Run2022C_3.list

# Final output
echo "Done. Files split:"
wc -l SingleMuon_Run2022C_*.list

