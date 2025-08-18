#!/bin/bash

# Dataset and redirector
DATASET="/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM"
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

# Step 4: Split into 2 files with 'Ext' in the name
head -n $HALF all_files.txt > WZTo2L2Q_Ext_1.list
tail -n +$((HALF + 1)) all_files.txt > WZTo2L2Q_Ext_2.list

# Final output
echo "Done. Files split:"
wc -l WZTo2L2Q_Ext_*.list

