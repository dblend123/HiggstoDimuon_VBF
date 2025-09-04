#!/bin/bash

# Dataset and redirector
DATASET="/DoubleMuon/Run2022C-22Sep2023-v1/NANOAOD"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count total lines
TOTAL=$(wc -l < all_files.txt)
HALF=$(( (TOTAL + 1) / 2 ))  # round up to avoid missing any file

# Step 4: Split into 2 files
head -n $HALF all_files.txt > DoubleMuon_Run2022C_1.list
tail -n +$((HALF + 1)) all_files.txt > DoubleMuon_Run2022C_2.list

# Final output
echo "Done. Files split:"
wc -l DoubleMuon_Run2022C_*.list

