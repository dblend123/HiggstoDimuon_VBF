# CMS H→μμ Run 3 — 2024 Synchronization

## Status

**Iowa synchronization ntuples published and validated.**

The 2024 synchronization work covers Run2024H Muon0/Muon1 data
and VBF H→μμ Monte Carlo, including the ext1 extension sample.

The final VBF category agreement between Iowa and TAMU is
**99.99658% across 3,571,455 matched original VBFH events**.

EOS directory:

    /eos/cms/store/group/phys_higgs/cmshmm/run3_sync/Iowa/update_after_Hmm_15Oct/

| Dataset | SyncEvents entries | Final VBF passes |
|---|---:|---:|
| Muon0 Run2024H | 2,128,353 | 10,603 |
| Muon1 Run2024H | 2,125,575 | 10,626 |
| VBFH 2024 + ext1 | 7,319,472 | 2,243,346 |

## Final VBF category

The final synchronization category is stored in:

    SyncEvents/passVBFCategoryApprox

Its requirements are:

- Leading-pT selected jet pair satisfying the VBF thresholds:
  pT1 > 35 GeV, pT2 > 30 GeV,
  mjj > 400 GeV, and |DeltaEta_jj| > 2.5.
- Opposite-sign muons.
- Both muons passing medium ID and isolation requirements.
- HLT event acceptance and trigger-object matching.
- Event filters and primary-vertex requirements.
- Exactly two Round 2 good muons.
- Electron veto.
- b-jet veto.
- Golden JSON requirement where applicable.

The Round 2 good-muon definition uses raw pT > 20 GeV,
medium ID, and PF relative isolation < 0.25.

Jet veto map exclusions remain applied at the individual-jet
selection level. The separate event-level veto-map
requirement is not part of the final category.

The original diagnostic branches remain available.

## Iowa–TAMU comparison

Matched original VBFH events: 3,571,455.

| Selection | Agreement | Disagreements |
|---|---:|---:|
| Previous Iowa category | 98.46337% | 54,880 |
| Round 2 muon multiplicity | 99.28449% | 25,554 |
| Final v7 category | 99.99658% | 122 |

Final v7 disagreement breakdown:

- Both pass: 1,109,804
- Iowa only: 93
- TAMU only: 29
- Both fail: 2,461,529

Of the 122 remaining disagreements, 100 coincide with
leading-jet VBF selection differences and 22 with
trigger-object matching differences.

These residual differences are deferred pending further
synchronization feedback.

## Implementation

Main analyzer:

    src/drivers/process_HiggsAnalysis_wCutsExtended_NoPairing_v0.cc

Input-file and event-count safeguards:

    src/defs/Streamer.h

ROOT category correction utility:

    sync_jobs/scripts/correct_vbf_category_v7.py

The correction utility recalculates the
passVBFCategoryApprox branch in existing ROOT ntuples.

It does not regenerate the underlying analyzer histograms,
other physics branches, or systematic variation trees.

## Validation

The published v7 files were checked for:

- Expected event counts.
- Correct default SyncEvents tree and ROOT key cycles.
- Exact agreement with the independently reconstructed v7
  category for every stored event.
- Preservation of all other SyncEvents branches.
- Preservation of the latest-cycle ROOT object inventory
  and histogram summary values.
- Correct published event counts and final category yields.

The finalized synchronization implementation is maintained
on the Git branch:

    hmumu-sync
