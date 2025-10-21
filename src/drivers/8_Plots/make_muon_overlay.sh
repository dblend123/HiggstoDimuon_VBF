#!/bin/bash
# -------------------------------------------------------------------
# Run the Muon Overlay plotting macro for a given era and correction mode
# Usage:
#   ./make_muon_overlay.sh <ERA> [MODE]
#
#   MODE options:
#     Rochester   (default)
#     KIT
#     NoMuCorr
#
# Examples:
#   ./make_muon_overlay.sh 12022
#   ./make_muon_overlay.sh 12022 KIT
#   ./make_muon_overlay.sh 12022 NoMuCorr
# -------------------------------------------------------------------

ERA=${1:-12022}
MODE=${2:-Rochester}

echo "============================================================"
echo "  Running Muon Overlay for Era: ${ERA}"
echo "  Correction Mode: ${MODE}"
echo "============================================================"

# Load ROOT (only if not already loaded)
if [ -z "$ROOTSYS" ]; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-centos7-gcc11-opt/setup.sh
fi

# Run ROOT macro
root -l -b -q "Muon_Overlay.C(\"Muon_Overlay_Z(\\\"${ERA}\\\",\\\"${MODE}\\\",true,false)\")"

echo "============================================================"
echo "  Done! Plots are in:"
echo "  /eos/user/d/dblend/StackedPlots/${ERA}/Muons/${MODE}/"
echo "============================================================"

