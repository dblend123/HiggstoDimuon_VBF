#!/bin/bash

# Dataset and redirector
DATASET="/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count total lines
TOTAL=$(wc -l < all_files.txt)
SIXTH=$(( (TOTAL + 5) / 6 ))  # round up for 6 parts

# Step 4: Split into 6 files
head -n $SIXTH all_files.txt > ZZTo2L2Q_Ext_1.list
tail -n +$((SIXTH + 1)) all_files.txt | head -n $SIXTH > ZZTo2L2Q_Ext_2.list
tail -n +$((2 * SIXTH + 1)) all_files.txt | head -n $SIXTH > ZZTo2L2Q_Ext_3.list
tail -n +$((3 * SIXTH + 1)) all_files.txt | head -n $SIXTH > ZZTo2L2Q_Ext_4.list
tail -n +$((4 * SIXTH + 1)) all_files.txt | head -n $SIXTH > ZZTo2L2Q_Ext_5.list
tail -n +$((5 * SIXTH + 1)) all_files.txt > ZZTo2L2Q_Ext_6.list

# Final output
echo "Done. Files split:"
wc -l ZZTo2L2Q_Ext_*.list

