# `9_TMVA/` — VBF H→μμ BDT Workflow

Complete training, scoring, and statistical-analysis pipeline for the Run-3
VBF H(125)→μμ search. Follows the strategy of the analysis note
[`AN2019_090_UCSD.pdf`](../../../AN2019_090_UCSD.pdf) and extends it to all
Run-3 eras (12022, 22022, 12023, 22023, 2024, 2025).

The pipeline is split into two driver scripts:

| Driver              | Stage                                  | Per-era? |
|---------------------|----------------------------------------|----------|
| `run_full.sh`       | Prepare → Train → Test → Score         | Yes (one era at a time) |
| `run_analysis.sh`   | Diagnostics → Categorization → Combined| Loops eras |

---

## Directory layout

```
9_TMVA/
├── run_full.sh                 # one-era training + scoring pipeline
├── run_analysis.sh             # post-scoring diagnostics + categorization
│
├── tmva100_DataPreparation.py  # define 18 BDT inputs, 4-fold split, snapshot ROOT files
├── tmva101_Training.py         # train one XGBoost BDT per (era, bkg category)
├── tmva102_Testing.py          # weighted ROC/AUC on the held-out fold
├── tmva103_Application.py      # score any *_all.root file → add `bdt_score` branch
│
├── vbf_bdt/                    # post-scoring analysis library
│   ├── bdt_validation.py       # per-era KS, score Data/MC, feature checks (Step 1 of analysis)
│   ├── bdt_event_categorization.py  # per-era + combined BDT categorization (Steps 2–3)
│   ├── bdt_analysis.py         # shared utilities (imported by validation/categorization)
│   ├── significance_inference.py    # Asimov Z, χ², discrete-profiling helpers
│   └── bdt_figures.py          # optional publication-style figure exporter
│
├── DataPreparation_Files/      # per-era snapshots (output of tmva100, input to tmva101–103)
└── Outputs/                    # everything produced by the workflow
    ├── Training/<era>/         # tmva101 model files + feature-importance txt
    ├── <era>/                  # tmva103 scored ROOT files (signal/bkg cats/data)
    ├── VBFValidation/<era>/    # bdt_validation reports
    ├── Categorization/<era>/   # per-era categorization reports
    ├── Categorization/combined/# multi-era combined categorization
    └── _archive/               # legacy PDFs/PNGs (kept for reference; not regenerated)
```

---

## `run_full.sh` — Training + scoring for one era

```bash
./run_full.sh -e <era> [--base <path>] \
              [--skip-training] [--skip-testing] [--skip-scoring]
```

**Defaults**

- `--base` → `/eos/user/d/dblend/Analyzed_SkimmedFiles/<era>/Rochester/noScaleSmear_noBtagSF`
- Model name inside TMVA file: `myBDT`
- Model file written to `Outputs/Training/<era>/<era>_tmva101_background.root`

### Step 0 — Data preparation (`tmva100_DataPreparation.py`)
- Reads the per-process NanoAOD-derived `mva` trees from `--base`.
- Defines the **18 BDT input variables** (see AN §6, Table 22 for the Run-2
  baseline; the Run-3 set here adds Collins–Soper angles φ_CS, cosθ_CS):

  `z` (Zeppenfeld), `pt_balance`, `mu1_pt_over_m`, `mu2_pt_over_m`,
  `mu1_eta`, `mu2_eta`, `y_mumu`, `pt_mumu`, `jet1_eta`, `deta_jj`,
  `dphi_jj`, `njets`, `pt_j1`, `pt_j2`, `dphi_mumu_jj`, `deta_mumu_jj`,
  `phi_CS`, `costheta_CS`.

  `m_μμ` is **deliberately excluded** to keep the score uncorrelated with the
  fit observable (AN §6.4).
- Assigns each event to one of 4 folds via a deterministic
  `(run, luminosityBlock, event)` hash. Folds **2** and **3** are reserved
  for validation and testing; folds **0–1** are used for training. This is
  the closest practical analog to the 4-fold scheme in AN §6.5.
- Writes:
  - `<era>_train_<cat>.root`, `<era>_test_<cat>.root`,
    `<era>_all_<cat>.root` for each of `signal, DY, EWK, TT, VV, VVV`
  - chained `<era>_train_background.root` / `<era>_test_background.root` /
    `<era>_all_background.root`

### Step 1 — BDT training (`tmva101_Training.py`)
- `run_full.sh` invokes the script with `--bkg-cats background` so a **single
  classifier per era** is trained against the merged background. The script
  also supports per-category training (`DY EWK TT VV VVV`) for diagnostics.
- Model: **`xgboost.XGBClassifier`**
  - `tree_method=hist, max_bin=30, max_depth=4, learning_rate=0.1,
     n_estimators=1000, subsample=0.5, min_child_weight=0.03,
     reg_lambda=1, scale_pos_weight=0.7, eval_metric=logloss`
- Weight handling (essential for weighted MC, AN §6.3):
  1. Sanitize (NaN→0, |w|, floor at 1e-9)
  2. **Negative-weight annihilation** — cancel equal sums of positive and
     negative MC weights and keep the surplus.
  3. **Class re-balance** so `Σw_S = Σw_B = √(Σw_S · Σw_B)`.
- Folds: train = (0, 1), valid = (2), test = (3) → ~50 / 25 / 25.
- Outputs (per era):
  - `Outputs/Training/<era>/<era>_xgb_<cat>.json` — native XGBoost
  - `Outputs/Training/<era>/<era>_tmva101_<cat>.root` — TMVA `RBDT` model
  - `Outputs/Training/<era>/<era>_featureimportance_<cat>.txt`

### Step 2 — Testing (`tmva102_Testing.py`)
- Loads the TMVA model and evaluates on fold==3 only.
- Reports **weighted AUC**, **TPR@1% FPR**, and the **Youden J** working
  point. Writes ROC PDF/PNG plus `tmva102_metrics.txt` and
  `tmva102_roc_points.csv` under `Outputs/`.

### Step 3 — Scoring (`tmva103_Application.py`)
- Clones any `_all.root` file and appends a single `bdt_score / F` branch
  using the trained model. Run for: signal, each background category, the
  merged background, and `DATA_all.root`. Output goes to `Outputs/<era>/`.

---

## `run_analysis.sh` — Diagnostics + categorization

```bash
./run_analysis.sh                # all eras + combined
./run_analysis.sh -e <era>       # single era (no combined step)
```

`MPLBACKEND=Agg` is forced so matplotlib does not need a display. The script
assumes the working directory is the analysis-note repo root.

### Step 1 — Per-era diagnostics (`vbf_bdt/bdt_validation.py`)
- KS test of the BDT score between folds (overfitting check, AN §6.6)
- Score distributions for data vs. stacked MC in the m_μμ sidebands
- Input-feature data/MC closure
- Output: `Outputs/VBFValidation/<era>/`

### Step 2 — Per-era categorization (`vbf_bdt/bdt_event_categorization.py`)
- Iterative quantile scan over the signal-score CDF to find category
  boundaries that maximize the per-era Asimov significance summed in
  quadrature over categories (AN §7.2). Background shape in each category
  is fit with a family of analytic PDFs (`BWZRedux`, `BWZGamma`,
  `BWZ×Bernstein`, `Landau×Bernstein`, `S-PowerLaw`, `S-Exponential`,
  `Bernstein`) and the best **Baker–Cousins χ²** (AN Ref [67], Table 27) is
  retained for each category.
- Output: `Outputs/Categorization/<era>/`

### Step 3 — Combined multi-era categorization
- Same script, called with all six eras simultaneously. Re-derives a single
  set of category boundaries that maximize the **combined-era Asimov Z** in
  the ±HWHM mass window of a double-Crystal-Ball signal model (AN §7.3).
  Skipped when `run_analysis.sh` is invoked with a single `-e <era>`.
- Output: `Outputs/Categorization/combined/`

---

## End-to-end recipe

```bash
cd src/drivers/9_TMVA

# 1) Per-era training + scoring (run once per era)
for era in 12022 22022 12023 22023 2024 2025; do
    ./run_full.sh -e $era
done

# 2) Diagnostics + per-era + combined categorization
./run_analysis.sh
```

The final expected significance numbers, χ²/ndf tables, and category
boundaries land in `Outputs/Categorization/combined/run.log`.

---

## Notes & references

- All references to the analysis note (`AN-2019/090`) above point to the
  PDF stored at the repository root: [`AN2019_090_UCSD.pdf`](../../../AN2019_090_UCSD.pdf).
- `Outputs/_archive/` holds legacy plots (`bdt_analysis_*.pdf`,
  `bdt_combined_result.pdf`, `2024_tmva_*.{pdf,png}`) that were produced by
  earlier iterations of the workflow and are no longer regenerated.
- Removed in the May 2026 cleanup (now superseded by `vbf_bdt/`):
  `plots_bdt_cuts.py`, `tmva_DataMC_BDT.py`, `tmva_KS_perfold.py`,
  `tmva_MassSculpt.py`, `tmva_PlotScored.py`, `tmva_ScanYields.py`.
