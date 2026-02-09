# VBF H→μμ Analysis Framework

**Analysis targeting VBF Higgs boson production with H→μμ decay in Run-3 pp collisions at √s = 13.6 TeV**

This framework processes NanoAOD datasets from 2022, 2023, and 2024 data-taking periods, applying VBF topology selections to maximize signal efficiency while suppressing Drell–Yan, top, and electroweak backgrounds.

---

## Table of Contents

- [Quick Start](#quick-start)
- [Analysis Workflow](#analysis-workflow)
- [Usage Guide](#usage-guide)
  - [0. Prerequisites & Build](#0-prerequisites--build)
  - [1. File Lists Preparation](#1-file-lists-preparation)
  - [2. Event Skimming](#2-event-skimming)
  - [3. Main Analysis](#3-main-analysis)
  - [4. Multivariate Analysis (TMVA)](#4-multivariate-analysis-tmva)
  - [5. Validation & Plotting](#5-validation--plotting)
- [Datasets & Event Selection](#datasets--event-selection)
- [Validation Results](#validation-results)
- [Next Steps](#next-steps)

---

## Quick Start

**Prerequisites:** CERN lxplus access with AFS workspace, EOS storage, HTCondor batch system

**Required Software:** CMake 3.15+, GCC 11+, ROOT 6.26+, Python 3.8+, Boost (via CVMFS)

```bash
# 1. Login and clone
ssh -Y username@lxplus.cern.ch
cd /afs/cern.ch/work/u/username/
git clone --recursive <repository_url>
cd HiggstoDimuon_VBF

# 2. Build
mkdir -p build && cd build
cmake ..
make -j4

# 3. Verify executables
ls src/drivers/6_Skimming/my_skim
ls src/drivers/process_HiggsAnalysis_wCutsExtended_NoPairing_v0
```

---

## Analysis Workflow

```
┌──────────────────────────────────────────────────────────────┐
│              CMS Data Acquisition System (DAS)               │
│         Raw collision data & Monte Carlo samples             │
│                  (NANOAOD/NANOAODSIM format)                 │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 1: File Lists (1_ListFiles)                           │
│  • Query DAS for dataset file paths per era/process          │
│  • Organize by data-taking period (2022/2023/2024)           │
│  • Generate input lists for batch job submission             │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 2: MC Weights (2_MCweight)                            │
│  • Compute Sum_genEventSumw for MC normalization             │
│  • Generate pileup reweighting histograms                    │
│  • Create weights.ini with cross-sections                    │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 3: RooFit Analysis (4_RooFit)                         │
│  • Fit Z→μμ peak to extract mass scale/resolution            │
│  • Model signal/background shapes for systematics            │
│  • Validate background modeling in control regions           │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 4: Tag-and-Probe (5_Tag_and_Probe)                    │
│  • Measure muon efficiencies using Z→μμ events               │
│  • Extract scale factors binned in pT and η:                 │
│    - ID: NUM_MediumID_DEN_TrackerMuons                       │
│    - Isolation: NUM_LoosePFIso_DEN_MediumID                  │
│    - Trigger: NUM_IsoMu24_DEN_CutBasedIdMedium               │
│  • Output: SF ROOT files → 7_Corrections/MUO/                │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 5: Corrections Repository (7_Corrections)             │
│  • Muon: Rochester corrections, efficiency SFs, FSR recovery │
│  • Jets: JEC/JER (era-specific), veto maps                   │
│  • Pileup: Reweighting histograms per Golden JSON            │
│  • Z-pT: Reweighting in bins of N_jets (0, 1, ≥2)            │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 6: Event Skimming (6_Skimming)                        │
│  Soft preselection to reduce dataset size (~90% reduction)   │
│  • Exactly 2 opposite-sign muons                             │
│  • Leading muon pT > 26 GeV (trigger-matched)                │
│  • 70 < m_μμ < 200 GeV                                        │
│  • Apply corrections: muon SF, JEC/JER, pileup, Z-pT         │
│  → Output: Skimmed ntuples to EOS                            │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 7: Main Analysis (process_HiggsAnalysis_*)            │
│  VBF topology selection and high-level variable computation  │
│  • Muon: leading pT >= 26 GeV, subleading pT >= 20 GeV       |
|  • |η| < 2.4 for both muons                                  |
|  • Muon isolation (pfRelIso04_all) <= 0.25 for both muons    |
|  • BDT mass window for filling: 110 < m_mumu < 150 GeV       |
|  • Z-region for Z-pT reweighting: 70 < m_mumu < 110 GeV      |
|  • Jets (VBF): lead pT >= 35 GeV, sublead pT >= 25 GeV       │
│  • Dijet: m_jj >= 400 GeV, |Δη_jj| >= 2.5                    │
│  • Note: b-jet veto and a central-jet veto are computed but  |
|    the veto lines are commented out.                         │
│  • Calculate: Zeppenfeld (z), pT balance, Δφ(μμ,jj), etc.    │
│  → Output: Analyzed ntuples with MVA variables to EOS        │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 8: Multivariate Analysis (9_TMVA)                     │
│  Train XGBoost BDT for signal/background separation          │
│  1. tmva100: Prepare training/test samples (5-fold split)    │
│  2. tmva101: Train model with hyperparameter tuning          │
│  3. tmva102: Validate performance (ROC, feature importance)  │
│  4. tmva103: Apply scores to all events                      │
│  → Output: Ntuples with BDT scores (0=bkg, 1=signal)         │
│                                                              │
│  ⚠ RETRAIN after changes to cuts, corrections, or variables │
└──────────────────────────┬───────────────────────────────────┘
                           │
                           ▼
┌──────────────────────────────────────────────────────────────┐
│  STAGE 9: Plotting & Results (8_Plots)                       │
│  • Control distributions (Z region validation)               │
│  • VBF topology: m_jj, Δη_jj, centrality                     │
│  • MVA performance and cut optimization                      │
│  • Final physics results: signal extraction, significance    │
└──────────────────────────────────────────────────────────────┘
```

---

## Usage Guide

### 0. Prerequisites & Build

Ensure valid proxy before submitting HTCondor jobs:

```bash
voms-proxy-init --voms cms -valid 192:00
```

Build the analysis framework:

```bash
cd /afs/cern.ch/work/u/username/HiggstoDimuon_VBF
mkdir -p build && cd build
cmake ..
make -j4

# Verify binaries
ls ../src/drivers/6_Skimming/my_skim
ls ../src/drivers/process_HiggsAnalysis_wCutsExtended_NoPairing_v0
```

---

### 1. File Lists Preparation

**Directory Structure:**
```
src/drivers/1_ListFiles/
├── 12022/          # 2022CD era
├── 22022/          # 2022EFG era (2022EE)
├── 12023/          # 2023C era
├── 22023/          # 2023D era (BPix)
└── 2024/           # 2024 era
```

**Query DAS and generate file lists:**

```bash
cd src/drivers/1_ListFiles/

# Data: Muon primary dataset (Run 2024 example)
dasgoclient --query "file dataset=/Muon/Run2024A-MINIv6NANOv15-v1/NANOAOD" > 2024/Muon_2024A.txt
dasgoclient --query "file dataset=/Muon/Run2024B-MINIv6NANOv15-v1/NANOAOD" > 2024/Muon_2024B.txt

# Signal MC: VBF H→μμ
dasgoclient --query "file dataset=/VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM" > 2024/VBF_H125.txt

# Background MC: Drell–Yan
dasgoclient --query "file dataset=/DYto2Mu-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v1/NANOAODSIM" > 2024/DYto2Mu-2Jets.txt

# Similar queries for: TTto2L2Nu, WZto*, ZZto*, WWW, WWZ, WZZ, ZZZ, etc.
```

**Verify:**
```bash
wc -l 2024/*.txt
head 2024/Muon_2024A.txt  # Check ROOT xrootd format
```

---

### 2. Event Skimming

**Purpose:** Apply soft preselection, reduce output by ~90%, store at CERN EOS

**Selection Criteria (`MyAnaSkim.C`):**
- Exactly 2 opposite-sign muons
- Leading muon pT > 26 GeV (trigger-matched)
- 70 < m_μμ < 200 GeV
- Muon SF, JEC/JER, pileup weights applied
- **Output tree:** `Events` with kinematic variables and per-event weights

**Submit jobs:**

```bash
cd src/drivers/6_Skimming/

# Edit Skim_condor_submit.jdl: set era and populate queue section
voms-proxy-init --voms cms -valid 192:00
condor_submit Skim_condor_submit.jdl

# Monitor
condor_q username
tail -f logs/*.log
```

**Job configuration:**
- **Executable:** `Skim_run_job.sh`
- **Binary:** `my_skim <inputList> <maxEvents>`
- **Output location:** `/eos/user/USERNAME/Skimmed_SkimmedFiles/[Era]/[Process]/`

**Validate outputs:**
```bash
xrdfs eosuser ls /eos/user/d/dblend/Skimmed_SkimmedFiles/2024/
root -l -q 'skim_file.root->Events->GetEntries()'
```

---

### 3. Main Analysis

**Purpose:** Apply VBF topology selection, calculate discriminants for TMVA

**Selection (`process_HiggsAnalysis_wCutsExtended_NoPairing_v0.cc`):**

**Dimuon:**
- Opposite-sign, isolated muons
- Leading muon pT > 26 GeV
- **110 < m_μμ < 160 GeV** (Higgs region)

**VBF Jets:**
- Two leading jets: pT > 35, 30 GeV
- Δη_jj > 3.5
- m_jj > 400 GeV

**Vetoes:**
- b-jet veto: ≤0 medium b-tags
- Central jet veto: N_jets(|η| < 2.5, pT > 25) ≤ 1

**Analysis variables:** Zeppenfeld (z), pT balance, Δφ(μμ,jj), plus 15+ MVA inputs

**Output tree:** `mva` with selection variables and per-event weights

**Build analyzer:**
```bash
cd src/drivers
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j4
```

**Submit jobs:**

```bash
cd src/drivers/

# Edit condor_submit.jdl: set era and populate queue section
voms-proxy-init --voms cms -valid 192:00
condor_submit condor_submit.jdl

# Monitor
condor_q username
tail -f HTCondor/logs/*.log
```

**Job arguments:**
```bash
--input <listfile>       # Skim file list
--output <tag.root>      # Output analyzed ntuple
--isMC <true|false>      # Is Monte Carlo?
--sample <sampleType>    # DY, TT, VV, VVV, VBF, DATA
--era <eraCode>          # 12022, 22022, 12023, 22023, 2024
--muCorr <type>          # roch (Rochester), none
--weights <filename>     # weights.ini with cross-sections
```

**Optional toggles** (append `--noXYZ` to disable):
- `--noPU`: Disable pileup reweighting
- `--noScaleSmear`: Disable muon scale/smear
- `--noMuonSF`: Disable muon efficiency SFs
- `--noJERC`: Disable jet energy corrections
- `--noJetVeto`: Disable jet veto map
- `--noGolden`: Disable Golden JSON (data only)

**Output location:** `/eos/user/USERNAME/Analyzed_SkimmedFiles/[Era]/[Process]/`

---

### 4. Multivariate Analysis (TMVA)

**CRITICAL:** Retrain TMVA whenever you:
- Change VBF selection cuts or mass windows
- Modify correction procedures (muon SF, JEC/JER, Z-pT reweighting)
- Add/remove analysis variables
- Update MC samples

**Step 1: Data Preparation**

```bash
cd src/drivers/9_TMVA/
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh

python3 tmva100_DataPreparation.py \
  --era 12022 \
  --base /eos/user/d/dblend/Analyzed_SkimmedFiles/{era}/ \
  --signal-subdir VBF \
  --treename-in mva \
  --treename-out mva \
  --weight-branch weight \
  --make-all \
  --outdir DataPreparation_Files/12022
```

**Outputs:**
- `12022_signal_fold0.root` ... `fold4.root` (5 folds)
- `12022_background_DY_fold0.root` ... (backgrounds by type)
- `12022_signal_all.root`, `12022_background_all.root` (combined)
- Tree: `mva` with variables + `weight` branch

**Step 2: Model Training**

```bash
python3 tmva101_Training.py \
  --era 12022 \
  --base-dir DataPreparation_Files/12022 \
  --bkg-cats background \
  --outdir Outputs/12022
```

**Hyperparameters:**
- n_estimators: 100
- max_depth: 5–6
- learning_rate: 0.1
- subsample: 0.8

**Outputs:**
- `12022_tmva101.root` (trained XGBoost model: "myBDT")
- `importance_*.png` (feature importance plots)
- `training_*.log` (convergence metrics)

**Step 3: Model Testing**

```bash
python3 tmva102_Testing.py \
  --era 12022 \
  --dir DataPreparation_Files/12022 \
  --model Outputs/12022/12022_tmva101.root \
  --name myBDT
```

**Outputs:** ROC curves, correlation matrices, score distributions

**Step 4: Apply Scores**

```bash
# Score signal + backgrounds
python3 tmva103_Application.py \
  --input DataPreparation_Files/12022_signal_all.root \
           DataPreparation_Files/12022_background_all.root \
  --output Outputs/12022/scored_all.root \
  --treename mva \
  --model-file Outputs/12022/12022_tmva101.root \
  --model-name myBDT

# Score data
python3 tmva103_Application.py \
  --input /eos/user/d/dblend/Analyzed_SkimmedFiles/12022/MC/DATA/DATA_all.root \
  --output Outputs/12022/data_scored_all.root \
  --treename mva \
  --model-file Outputs/12022/12022_tmva101.root \
  --model-name myBDT
```

**Output:** Tree `mva` + new branch `myBDT_score` (0=background-like, 1=signal-like)

**Top Discriminating Variables by Era:**

| Era | Top 3 Features |
|-----|----------------|
| **2022** | 1. pt_balance (13.46%)<br>2. z (11.82%)<br>3. njets (11.26%) |
| **2022EE** | 1. pt_balance (23.00%)<br>2. pt_mumu (10.02%)<br>3. z (9.42%) |
| **2023** | 1. pt_balance (17.46%)<br>2. pt_j2 (9.47%)<br>3. z (8.69%) |
| **2023BPix** | 1. pt_balance (15.47%)<br>2. z (9.33%)<br>3. dphi_mumu_jj (8.72%) |

---

### 5. Validation & Plotting

**Validation plots (Z control region):**

```bash
cd src/drivers/8_Plots/
root -l

# Muon momentum corrections
.L Checking_Corrections/Muon_Overlay_AllPlots.C+
.x Checking_Corrections/Muon_Overlay_AllPlots.C()

# Z-pT reweighting (critical for jet multiplicity)
.L ZPt_Corrections/ZPt_Overlay_Plots.C+
.x ZPt_Corrections/ZPt_Overlay_Plots.C()
```

**Final result plots:**
- Invariant mass (m_μμ) in signal region (110–160 GeV)
- MVA score distributions and cut optimization
- VBF topology: m_jj, Δη_jj, centrality, pT balance
- Trigger/ID/Isolation efficiencies

---

## Complete Workflow Example (2024 Era)

```bash
# === Setup ===
ssh -Y username@lxplus.cern.ch
cd /afs/cern.ch/work/u/username/
git clone --recursive <repo_url>
cd HiggstoDimuon_VBF
mkdir -p build && cd build
cmake .. && make -j4
cd ..
voms-proxy-init --voms cms

# === File Lists ===
cd src/drivers/1_ListFiles/2024/
dasgoclient --query "file dataset=/Muon/Run2024A-MINIv6NANOv15-v1/NANOAOD" > Muon_2024A.txt
# ... (repeat for all samples)
cd ../..

# === Skimming ===
cd 6_Skimming/
# Edit Skim_condor_submit.jdl
condor_submit Skim_condor_submit.jdl
cd ..

# === Main Analysis ===
# Edit condor_submit.jdl
condor_submit condor_submit.jdl
cd ..

# === TMVA Training ===
cd 9_TMVA/
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt/setup.sh
python3 tmva100_DataPreparation.py --era 2024 --signal-subdir VBF --make-all --outdir DataPreparation_Files/2024
python3 tmva101_Training.py --era 2024 --base-dir DataPreparation_Files/2024 --bkg-cats background --outdir Outputs/2024
python3 tmva102_Testing.py --era 2024 --dir DataPreparation_Files/2024 --model Outputs/2024/2024_tmva101.root
python3 tmva103_Application.py --input DataPreparation_Files/2024_signal_all.root DataPreparation_Files/2024_background_all.root --output Outputs/2024/scored_all.root --model-file Outputs/2024/2024_tmva101.root
cd ..

# === Validation ===
cd 8_Plots/
root -l 'Checking_Corrections/Muon_Overlay_AllPlots.C()'
root -l 'ZPt_Corrections/ZPt_Overlay_Plots.C()'
```

---

## Datasets & Event Selection

### 1. Data Samples

**Muon primary datasets** are used for all data-taking years:

**Run 2022:**
```
/Muon/Run2022*-22Sep2023-v*/NANOAOD
```

**Run 2023:**
```
/Muon[0,1]/Run2023*-22Sep2023_v*-v*/NANOAOD
```

**Run 2024:**
```
/Muon[0,1]/Run2024*-MINIv6NANOv15-v*/NANOAOD
```

Only certified good runs and lumisections are considered. Events require at least one good reconstructed primary vertex (`npvsGood ≥ 1`).

---

### 2. Monte Carlo Samples

#### Signal

**Vector Boson Fusion Higgs production with decay to dimuons:**

**2022–2023:**
```
/VBFHto2Mu_M-125_TuneCP5_withDipoleRecoil_13p6TeV_powheg-pythia8
/Run3*NanoAODv12-130X_mcRun3_202*_realistic_v*/NANOAODSIM
```

**2024:**
```
/VBFH-Hto2Mu_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8
/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM
```

**Production cross section × branching ratio:**
```
σ × BR (VBF H → μμ) = 9.09568 × 10⁻⁴ pb
```

#### Backgrounds

- **Drell–Yan + jets:**
  ```
  /DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8
  ```
  (2024: `DYto2Mu-2Jets`, binned in MLL)

- **Electroweak Z + 2 jets:**
  ```
  /EWK_2L2J_TuneCH3_13p6TeV_madgraph-herwig7
  ```

- **Top pair (dileptonic):**
  ```
  /TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8
  ```

- **Diboson:**
  ```
  /WZto2L2Q, /ZZto2L2Q, /ZZto2L2Nu, /WZto3LNu, /ZZto4L
  ```

- **Triboson:**
  ```
  /WWW_4F, /WWZ_4F, /WZZ, /ZZZ
  ```

All MC samples generated at √s = 13.6 TeV, processed with NanoAODv12 (2022–2023) or NanoAODv15 (2024).

---

### 3. Event Preselection

#### Muons

- Exactly two reconstructed muons
- Opposite electric charge
- At least one muon matched to `HLT_IsoMu24`
- **Leading trigger-matched muon:** pT > 26 GeV
- **Additional lepton veto:**
  - No extra muons with pT > 20 GeV and |η| < 2.4
  - No electrons with pT > 20 GeV and |η| < 2.5

#### Jets

- Reconstructed with AK4 PUPPI
- **Jet ID:** `jetId ≥ 4` (Tight or TightLeptonVeto)
- **Kinematic requirements:**
  - pT > 25 GeV
  - |η| < 4.7
- **Muon–jet overlap removal:** ΔR(μ, jet) > 0.4

#### Vertices

- At least one good primary vertex: `npvsGood ≥ 1`

---

### 4. Jet Energy Corrections and Resolution

Jet energy corrections (JEC) and jet energy resolution (JER) applied using year- and era-specific calibrations:

#### 2022 (CD)

**JEC:**
- `Summer22_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi`
- `Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi`
- `Summer22_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi`

**JER:**
- `Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi`
- `Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi`

#### 2022 (EE)

**JEC:**
- `Summer22EE_22Sep2023_V3_MC_L1FastJet_AK4PFPuppi`
- `Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi`
- `Summer22EE_22Sep2023_V3_MC_L3Absolute_AK4PFPuppi`

**JER:**
- `Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi`
- `Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi`

#### 2023 (C)

**JEC:**
- `Summer23Prompt23_V2_MC_L1FastJet_AK4PFPuppi`
- `Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi`
- `Summer23Prompt23_V2_MC_L3Absolute_AK4PFPuppi`

**JER:**
- `Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi`
- `Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi`

#### 2023 (D, BPix)

**JEC:**
- `Summer23BPixPrompt23_V3_MC_L1FastJet_AK4PFPuppi`
- `Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi`
- `Summer23BPixPrompt23_V3_MC_L3Absolute_AK4PFPuppi`

**JER:**
- `Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi`
- `Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi`

---

### 5. Jet Veto Maps

Applied to suppress detector noise and problematic regions:

- **2022 CD:** `Summer22_23Sep2023_RunCD_V1`
- **2022 EE:** `Summer22EE_23Sep2023_RunEFG_V1`
- **2023:** `Summer23Prompt23_RunC_V1`
- **2023 BPix:** `Summer22_23Sep2023_RunCD_V1`
- **2024:** `Summer24Prompt24_RunBCDEFGHI_V1`

---

### 6. b-Tagging

No b-tagging scale factors currently applied. Available configurations (not in use):

- `particleNet_comb`
- `particleNet_light`
- `particleNet_wp_values`

---

### 7. Cross Sections

Background cross sections from official CMS MATRIX recommendations at √s = 13.6 TeV:

https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV

---

### 8. Applied Corrections

#### Muons

**Efficiency Scale Factors (from tag-and-probe):**
- ID: `NUM_MediumID_DEN_TrackerMuons`
- Isolation: `NUM_LoosePFIso_DEN_MediumID`
- Trigger: `NUM_IsoMu24_DEN_CutBasedIdMedium`

**Rochester corrections:** under re-validation

**FSR photon recovery:** `Muon_pt` → `Muon_bsConstrainedPt`

#### Jets (Era-Specific)

See [Section 4: Jet Energy Corrections](#4-jet-energy-corrections-and-resolution) for full JEC/JER details.

#### Pileup Reweighting

Era-specific Golden JSON profiles:
- **2022CD:** `Collisions2022_355100_357900_eraBCD_GoldenJson`
- **2022EFG:** `Collisions2022_359022_362760_eraEFG_GoldenJson`
- **2023C:** `Collisions2023_366403_369802_eraBC_GoldenJson`
- **2023D:** `Collisions2023_369803_370790_eraD_GoldenJson`

#### Z-pT Reweighting

Corrects Drell-Yan pT mismodeling in bins of:
- N_jets: 0, 1, ≥2
- pT ranges: 0–60, 60–120, 120–200 GeV

**Critical for:** Jet multiplicity agreement and MVA training accuracy

---

## Validation Results

Validation performed in both **Z control region** (muon-dominated) and **VBF-enriched region** (jet-dominated). Results show strong muon modeling across all eras, with discrepancies in VBF region driven by jet-related observables.

### Z Control Region

**Selection:**
- Opposite-sign dimuons
- 70 < m_μμ < 200 GeV
- No VBF jet requirements

**Purpose:** Validate muon reconstruction, momentum scale, trigger/ID efficiencies, baseline background normalization

| Era | Data/MC | Events | Assessment |
|-----|---------|--------|------------|
| **2022** | ~1.01 | 6.3M | ✓ Very good |
| **2022EE** | ~0.99 | 15.8M | ✓ Very good |
| **2023** | ~0.98 | 13.3M | ✓ Very good |
| **2023BPix** | ~1.01 | 6.9M | ✓ Very good |
| **2024** | ~0.99 | 74.2M | ✓ Very good |

**Interpretation:**
- Overall muon Data/MC agreement is good across all eras
- Small residual discrepancies in dimuon pT distributions consistent with recoil modeling
- Effects subdominant before jet-based selections

---

### Jet Validation (Before VBF Criteria)

Jets validated using same events as Z control region, prior to VBF topology cuts.

| Era | Jet Data/MC (pre-cuts) | Assessment |
|-----|------------------------|------------|
| **2022** | ~100% | ✓ Good |
| **2022EE** | ~100% | ✓ Good |
| **2023** | ~100% | ✓ Good |
| **2023BPix** | ~100% | ✓ Good |
| **2024** | ~95–100% | ✓ Acceptable |

**Interpretation:**
- Jet reconstruction and baseline JEC/JER consistent with muons before tight selections
- No major reconstruction-level issues confirmed

---

### VBF-Enriched Region

**Selection:**
- m_jj > 400 GeV
- Δη_jj > 3.5
- Full event preselection applied

**Purpose:** Sensitive to jet energy scale, resolution, and correlations

| Era | Data/MC | Events | Assessment |
|-----|---------|--------|------------|
| **2022** | ~0.85 | 725K | ⚠ Moderate deficit |
| **2022EE** | ~0.85 | 1.7M | ⚠ Moderate deficit |
| **2023** | ~0.85 | 1.3M | ⚠ Moderate deficit |
| **2023BPix** | ~0.87 | 679K | ⚠ Moderate deficit |
| **2024** | ~0.50 | 8.4M | ❗ Significant deficit |

**Interpretation:**

After VBF cuts, Data/MC agreement degrades significantly:
- **~15% deficit in 2022–2023**
- **~50% deficit in 2024**

Discrepancy not present before jet criteria, indicating issues with:
- JEC/JER modeling at high m_jj
- Jet correlations (balance, topology)
- Possibly incomplete/preliminary 2024 calibrations

---

### MVA Feature Sensitivity

Feature importance confirms dominant discriminating power is jet-driven, making observed discrepancies physically consistent.

| Era | Top Features |
|-----|--------------|
| **2022** | pt_balance (13.46%), z (11.82%), njets (11.26%) |
| **2022EE** | pt_balance (23.00%), pt_mumu (10.02%), z (9.42%) |
| **2023** | pt_balance (17.46%), pt_j2 (9.47%), z (8.69%) |
| **2023BPix** | pt_balance (15.47%), z (9.33%), dphi_mumu_jj (8.72%) |

**Implication:**
- Leading feature (pt_balance) directly couples muons and jets
- Amplifies jet mis-modeling after VBF cuts
- Good Z-region agreement does not guarantee VBF-region agreement

---

### Summary

✅ **Muon reconstruction and calibration validated across Run 3**

⚠ **Jet modeling is dominant source of discrepancy, especially after VBF selections**

❗ **2024 requires correction revalidation before inclusion in combined results**

**The observed behavior is consistent with MVA feature sensitivities, not a pathological failure**

---

## Next Steps

### Immediate

- Re-implement and validate JEC/JER, especially 2024
- Re-check jet veto maps and era matching
- Repeat VBF validation before and after each correction step

### Near-Term

- Re-validate Rochester corrections (cross-check vs Z peak)
- Introduce b-tagging scale factors (ParticleNet) and assess impact on jet balance
- Monitor MVA mass sculpting after corrections

### Final

- Freeze per-era corrections
- Combine eras with systematic uncertainties
- Proceed to statistical inference

---