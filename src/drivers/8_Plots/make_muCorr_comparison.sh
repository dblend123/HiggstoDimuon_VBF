#!/bin/bash
# -------------------------------------------------------------------
# Run the Muon MuCorr comparison plotting macro for a given era
# Usage:
#   ./make_muCorr_comparison.sh 12022
# -------------------------------------------------------------------

ERA=${1:-12022}  # Default to 12022 if not provided

echo "============================================================"
echo "  Running Muon MuCorr Comparison for Era: ${ERA}"
echo "============================================================"

# Load CMS environment if needed
if [ -z "$ROOTSYS" ]; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
fi

# Run ROOT macro
root -l -b -q "Muon_MuCorrComparison.C(\"${ERA}\")"

echo "============================================================"
echo "  Done! Plots are in:"
echo "  /eos/user/d/dblend/StackedPlots/${ERA}/MuCorr_Comparison/"
echo "============================================================"

