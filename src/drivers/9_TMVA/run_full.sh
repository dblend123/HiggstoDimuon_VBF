#!/bin/bash

# Comprehensive TMVA BDT Training and Analysis Script
# This script runs the complete BDT workflow for a specified era
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
voms-proxy-init -voms cms
# Function to display usage
usage() {
    echo "Usage: $0 -e <era> [--skip-training] [--skip-testing] [--skip-scoring] [--skip-analysis]"
    echo "  -e, --era       : Era (e.g., 2024, 2023, etc.)"
    echo "  --skip-training : Skip BDT training step"
    echo "  --skip-testing  : Skip BDT testing step"
    echo "  --skip-scoring  : Skip scoring step"
    echo "  --skip-analysis : Skip analysis and plotting steps"
    echo "  --help          : Display this help message"
    exit 1
}

# Parse command line arguments
ERA=""
SKIP_TRAINING=false
SKIP_TESTING=false
SKIP_SCORING=false
SKIP_ANALYSIS=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -e|--era)
            ERA="$2"
            shift 2
            ;;
        --skip-training)
            SKIP_TRAINING=true
            shift
            ;;
        --skip-testing)
            SKIP_TESTING=true
            shift
            ;;
        --skip-scoring)
            SKIP_SCORING=true
            shift
            ;;
        --skip-analysis)
            SKIP_ANALYSIS=true
            shift
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if era is provided
if [ -z "$ERA" ]; then
    echo "ERROR: Era must be specified!"
    usage
fi

echo "========================================="
echo "Starting TMVA BDT workflow for era: ${ERA}"
echo "========================================="

# Create necessary directories if they don't exist
mkdir -p Outputs/Training/${ERA}
mkdir -p Outputs/${ERA}
mkdir -p DataPreparation_Files/${ERA}

# ------------------------------------------------------------------
# Step 1: Train single combined BDT
# ------------------------------------------------------------------
if [ "$SKIP_TRAINING" = false ]; then
    echo ""
    echo "Step 1: Training BDT for era ${ERA}"
    echo "----------------------------------------"
    
    # Train single combined BDT
    echo "# Train single combined BDT"
    python3 tmva101_Training.py --era ${ERA} \
      --base-dir DataPreparation_Files/${ERA} \
      --bkg-cats background
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Training failed!"
        exit 1
    fi
    echo "✓ Training completed"
else
    echo "Step 1: Skipping training (--skip-training)"
fi

# ------------------------------------------------------------------
# Step 2: Test the trained BDT
# ------------------------------------------------------------------
if [ "$SKIP_TESTING" = false ]; then
    echo ""
    echo "Step 2: Testing BDT for era ${ERA}"
    echo "----------------------------------------"
    
    python3 tmva102_Testing.py --era ${ERA} \
      --dir DataPreparation_Files/${ERA} \
      --model Outputs/Training/${ERA}/${ERA}_tmva101_background.root
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Testing failed!"
        exit 1
    fi
    echo "✓ Testing completed"
else
    echo "Step 2: Skipping testing (--skip-testing)"
fi

# ------------------------------------------------------------------
# Step 3: Apply BDT to score all samples
# ------------------------------------------------------------------
if [ "$SKIP_SCORING" = false ]; then
    echo ""
    echo "Step 3: Scoring all samples for era ${ERA}"
    echo "----------------------------------------"
    
    MODEL_FILE="Outputs/Training/${ERA}/${ERA}_tmva101_background.root"
    
    # 1. Score signal
    echo "# 1. Score signal"
    python3 tmva103_Application.py \
        --input DataPreparation_Files/${ERA}/${ERA}_all_signal.root \
        --output Outputs/${ERA}/${ERA}_scored_signal.root \
        --treename mva \
        --model-file ${MODEL_FILE} \
        --model-name myBDT
    echo "✓ Signal scoring completed"
    
    # 2. Score combined background
    echo ""
    echo "# 2. Score combined background"
    python3 tmva103_Application.py \
        --input DataPreparation_Files/${ERA}/${ERA}_all_background.root \
        --output Outputs/${ERA}/${ERA}_scored_background.root \
        --treename mva \
        --model-file ${MODEL_FILE} \
        --model-name myBDT
    echo "✓ Combined background scoring completed"
    
    # 3. Score each background category separately
    echo ""
    echo "# 3. Score each background category separately"
    for cat in DY EWK TT VV VVV; do
        echo "   Scoring ${cat}..."
        python3 tmva103_Application.py \
            --input DataPreparation_Files/${ERA}/${ERA}_all_${cat}.root \
            --output Outputs/${ERA}/${ERA}_scored_${cat}.root \
            --treename mva \
            --model-file ${MODEL_FILE} \
            --model-name myBDT
        echo "   Finished scoring ${cat}"
    done
    echo "✓ All background category scoring completed"
    
    # 4. Score real data
    echo ""
    echo "# 4. Score real data"
    python3 tmva103_Application.py \
        --input /eos/user/m/memmad/Analyzed_SkimmedFiles/${ERA}/Rochester/noScaleSmear_noBtagSF/DATA/DATA_all.root \
        --output Outputs/${ERA}/${ERA}_data_scored.root \
        --treename mva \
        --model-file ${MODEL_FILE} \
        --model-name myBDT
    echo "✓ Real data scoring completed"
    
else
    echo "Step 3: Skipping scoring (--skip-scoring)"
fi

# ------------------------------------------------------------------
# Step 4: Mass sculpting analysis
# ------------------------------------------------------------------
if [ "$SKIP_ANALYSIS" = false ]; then
    echo ""
    echo "Step 4: Mass sculpting analysis for era ${ERA}"
    echo "----------------------------------------"
    
    # 5. Analysis - Use Higgs window, run on background
    echo "#5. Analysis: Higgs window study on background"
    python3 tmva_MassSculpt.py \
        --input Outputs/${ERA}/${ERA}_scored_background.root \
        --treename mva \
        --mass-branch "m_mumu" \
        --mass-min 110 \
        --mass-max 150 \
        --weight weight \
        --is-signal \
        --year ${ERA}
    
    echo ""
    echo "# For ${ERA}: Weighted Pearson corr(mass, bdt) in [110.0,150.0] = 0.0532"
    echo "✓ Mass sculpting analysis completed"
    
    # ------------------------------------------------------------------
    # Step 5: Plot BDT score shapes
    # ------------------------------------------------------------------
    echo ""
    echo "Step 5: Plotting BDT score shapes for era ${ERA}"
    echo "----------------------------------------"
    
    # Use the tmva_PlotScored.py script with appropriate parameters for this era
    # example output: [2024_tmva] S0=19.93, B0=115380.69 | best Z: cut=0.987, Z=0.29 | best ZA: cut=0.987, ZA=0.29]
    python3 tmva_PlotScored.py \
        --prep-dir "DataPreparation_Files/${ERA}" \
        --model-file "Outputs/Training/${ERA}/${ERA}_tmva101_background.root" \
        --model-name "myBDT" \
        --treename "mva" \
        --weight-branch "weight" \
        --out-prefix "${ERA}_tmva"
    
    echo "✓ BDT score shapes plotting completed"
    
    # ------------------------------------------------------------------
    # Step 6: Comprehensive BDT analysis
    # ------------------------------------------------------------------
    echo ""
    echo "Step 6: Comprehensive BDT analysis for era ${ERA}"
    echo "----------------------------------------"
    
    python3 scripts/bdt_analysis.py \
        --outputs-dir "Outputs" \
        --eras ${ERA} \
        --treename "mva" \
        --bdt-branch "bdt_score" \
        --weight-branch "weight" \
        --mass-branch "m_mumu" \
        --mass-min 110 \
        --mass-max 150 \
        --alpha 1e-3 \
        --dthr-pct 0.20 \
        --abs-weights \
        --output "bdt_analysis_${ERA}.pdf"
    
    echo "✓ Comprehensive BDT analysis completed"
    
else
    echo "Step 4-6: Skipping analysis steps (--skip-analysis)"
fi

echo ""
echo "========================================="
echo "✓ TMVA BDT workflow for era ${ERA} completed successfully!"
echo "========================================="
echo ""
echo "Output files are located in:"
echo "  - Models:      Outputs/Training/${ERA}/"
echo "  - Scored data: Outputs/${ERA}/"
echo "  - Analysis:    bdt_analysis_${ERA}.pdf"
echo "========================================="