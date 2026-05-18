source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
cd /eos/user/m/memmad/HiggstoDimuon_VBF
export MPLBACKEND=Agg

OUT=src/drivers/9_TMVA/Outputs
DIAG_SCRIPT=src/drivers/9_TMVA/vbf_bdt/bdt_validation.py
CAT_SCRIPT=src/drivers/9_TMVA/vbf_bdt/bdt_event_categorization.py

# ------------------------------------------------------------------
# Step 1: Per-era diagnostic reports (KS folds, ROC, feature checks)
# ------------------------------------------------------------------
for era in 12022 22022 12023 22023 2024 2025; do
  ODIR="$OUT/VBFValidation/$era"
  mkdir -p "$ODIR"
  python3 "$DIAG_SCRIPT" \
    --era "$era" \
    --outputs-dir "$OUT" \
    --output-dir "$ODIR" \
    --log-level INFO |& tee "$ODIR/run.log"
done

# ------------------------------------------------------------------
# Step 2: Multi-era categorization, category optimisation, blinded
#         sideband validation report (the primary analysis output)
# ------------------------------------------------------------------
CATDIR="$OUT/Categorization"
mkdir -p "$CATDIR"
python3 "$CAT_SCRIPT" \
  --eras 12022 22022 12023 22023 2024 2025 \
  --outputs-dir "$OUT" \
  --outdir "$CATDIR" \
  --log-level INFO |& tee "$CATDIR/run.log"