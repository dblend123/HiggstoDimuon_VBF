#!/bin/bash

# Dataset and redirector
DATASET="/TTH_Hto2Mu_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM"
REDIRECTOR="root://cmsxrootd.fnal.gov/"

# Step 1: Query DAS and save file list
echo "Querying DAS for dataset..."
dasgoclient --query="file dataset=${DATASET}" --limit=0 > all_files.txt

# Step 2: Prepend redirector
echo "Adding redirector..."
sed -i "s|^|${REDIRECTOR}|" all_files.txt

# Step 3: Count lines
TOTAL=$(wc -l < all_files.txt)

# Step 4: Compute line distribution
LINES1=$(( (TOTAL + 2) / 3 ))
LINES2=$LINES1
LINES3=$(( TOTAL - LINES1 - LINES2 ))

# Step 5: Split into 3 list files
head -n $LINES1 all_files.txt > TTHTo2Mu_1.list
tail -n +$((LINES1 + 1)) all_files.txt | head -n $LINES2 > TTHTo2Mu_2.list
tail -n $LINES3 all_files.txt > TTHTo2Mu_3.list

# Final output
echo "Done. Files split:"
wc -l TTHTo2Mu_*.list

