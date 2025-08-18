#!/bin/bash

# Dataset and redirector
DATASET="/WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v3/NANOAODSIM"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count total lines
TOTAL=$(wc -l < all_files.txt)
HALF=$(( (TOTAL + 1) / 2 ))  # round up for 2 parts

# Step 4: Split into 2 files with 'Ext' tag and unique filenames
head -n $HALF all_files.txt > WZto3LNu_ext1_v3_Ext_1.list
tail -n +$((HALF + 1)) all_files.txt > WZto3LNu_ext1_v3_Ext_2.list

# Final output
echo "Done. Files split:"
wc -l WZto3LNu_ext1_v3_Ext_*.list

