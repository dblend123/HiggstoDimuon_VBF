#!/bin/bash
# run_full.sh - Complete H→μμ analysis workflow

set -e  # Exit on error

echo "Starting full H→μμ VBF analysis..."

# Configuration
MODEL="Outputs/Training/2024/2024_tmva101_background.root"
TREENAME="mva"
OUTDIR="Outputs/2024"

# 1. Score signal
echo "[1/6] Scoring signal MC..."
python3 tmva103_Application.py \
  --input DataPreparation_Files/2024/2024_all_signal.root \
  --output ${OUTDIR}/2024_scored_signal.root \
  --treename ${TREENAME} \
  --model-file ${MODEL} \
  --model-name myBDT

# 2. Score combined background
echo "[2/6] Scoring combined background MC..."
python3 tmva103_Application.py \
  --input DataPreparation_Files/2024/2024_all_background.root \
  --output ${OUTDIR}/2024_scored_background.root \
  --treename ${TREENAME} \
  --model-file ${MODEL} \
  --model-name myBDT

# 3. Score individual background components
echo "[3/6] Scoring individual background components..."
for cat in DY EWK TT VV VVV; do
  echo "  - Scoring ${cat}..."
  python3 tmva103_Application.py \
    --input DataPreparation_Files/2024/2024_all_${cat}.root \
    --output ${OUTDIR}/2024_scored_${cat}.root \
    --treename ${TREENAME} \
    --model-file ${MODEL} \
    --model-name myBDT
done

# 4. Score real data
echo "[4/6] Scoring collision data..."
python3 tmva103_Application.py \
  --input /eos/user/m/memmad/Analyzed_SkimmedFiles/2024/Rochester/noScaleSmear_noBtagSF/DATA/DATA_all.root \
  --output ${OUTDIR}/2024_data_scored.root \
  --treename ${TREENAME} \
  --model-file ${MODEL} \
  --model-name myBDT

# 5. BDT diagnostic plots
echo "[5/6] Creating BDT diagnostic plots..."
python3 scripts/view_bdt_scores.py \
  --input ${OUTDIR}/2024_scored_test.root \
  --signal-file ${OUTDIR}/2024_scored_signal.root \
  --bkg-file ${OUTDIR}/2024_scored_background.root \
  --treename ${TREENAME} \
  --output ${OUTDIR}/bdt_analysis_enhanced.pdf

# 6. Full comprehensive analysis
echo "[6/6] Running comprehensive analysis..."
python3 comprehensive_hmumu_analysis.py \
  --data ${OUTDIR}/2024_data_scored.root \
  --signal ${OUTDIR}/2024_scored_signal.root \
  --background ${OUTDIR}/2024_scored_background.root \
  --bkg-components ${OUTDIR}/2024_scored_DY.root \
                    ${OUTDIR}/2024_scored_EWK.root \
                    ${OUTDIR}/2024_scored_TT.root \
                    ${OUTDIR}/2024_scored_VV.root \
                    ${OUTDIR}/2024_scored_VVV.root \
  --treename ${TREENAME} \
  --output hmumu_full_analysis.pdf \
  --lumi 0.13015 \
  --signal-xsec 0.000909568

echo ""
echo "===================="
echo "ANALYSIS COMPLETE!"
echo "===================="
echo "Output files:"
echo "  - BDT diagnostics: ${OUTDIR}/bdt_analysis_enhanced.pdf"
echo "  - Full analysis:   hmumu_full_analysis.pdf"
echo "===================="