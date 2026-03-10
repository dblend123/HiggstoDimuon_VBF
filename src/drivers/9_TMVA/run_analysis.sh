source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
cd /eos/user/m/memmad/HiggstoDimuon_VBF
export MPLBACKEND=Agg

OUT=src/drivers/9_TMVA/Outputs
SCRIPT=src/drivers/9_TMVA/scripts/vbf_bdt_validation.py

for era in 12022 12023 22022 22023 2024; do
  ODIR="$OUT/VBFValidation/$era"
  mkdir -p "$ODIR"
  python3 "$SCRIPT" \
    --era "$era" \
    --outputs-dir "$OUT" \
    --output-dir "$ODIR" \
    --log-level INFO |& tee "$ODIR/run.log"
done