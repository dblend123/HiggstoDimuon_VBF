#!/bin/bash

# Dataset and redirector
DATASET="/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v1/NANOAODSIM"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count total lines
TOTAL=$(wc -l < all_files.txt)
THIRD=$(( (TOTAL + 2) / 3 ))  # round up to ensure all lines are included

# Step 4: Split into 3 files with 'Ext' in the name
head -n $THIRD all_files.txt > DYto2L2Jets_Ext_1.list
tail -n +$((THIRD + 1)) all_files.txt | head -n $THIRD > DYto2L2Jets_Ext_2.list
tail -n +$((2 * THIRD + 1)) all_files.txt > DYto2L2Jets_Ext_3.list

# Final output
echo "Done. Files split:"
wc -l DYto2L2Jets_Ext_*.list

