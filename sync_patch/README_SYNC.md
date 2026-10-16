# H→μμ synchronization: initial implementation (NOT validated production output)

The diff `hmumu_sync_initial.patch` modifies `src/drivers/process_HiggsAnalysis_wCutsExtended_NoPairing_v0.cc` in the submitted repository. Apply from repository root with `patch -p1 < hmumu_sync_initial.patch` (on a clean copy matching the uploaded ZIP).

## Implemented in this initial patch
- `--sync` option creates a `SyncEvents` ROOT TTree, one row per event meeting raw two-muon pT >20 GeV and raw dimuon mass 50–200 GeV, using the first pair by NanoAOD index that satisfies those conditions. No opposite-sign, medium ID, iso, event filter, Golden JSON, vertex, extra-lepton or trigger cut is applied to output. Such decisions are recorded as flags (the standard analysis path is bypassed when --sync is set).
- Stores run/lumi/event/sample, selected pair indices, raw and correction-base pT, beam-constrained pT, corrected pT, raw/corrected-before-FSR/corrected-after-FSR dimuon mass and pT, charge, eta, phi, isolation, raw jet quantities, MET input, component generator/PU and approximate muon SF/trigger factors. Disables Z-pT reweighting in sync mode.
- Corrects FSR candidate association for the normal analysis: use the actual selected muon indices rather than 0/1. Binds `FsrPhoton_relIso03` if available and maps the unsigned electron-index sentinel 65535 to -1.

## Important remaining work before submitting to CMS synchronization
- Needs ROOT/CMS compilation and runtime tests against the exact Run2024H Muon0, Muon1 and VBF Summer24 NanoAODv15 files. ROOT was not available in the editing environment. No events or ROOT output have been validated.
- Raw jet quantities currently include *input* `Jet_pt` and `Jet_pt*(1-Jet_rawFactor)`, not the analysis's fully recalculated JERC pT. Record final corrected jets, cleaning, veto, jet pair indices and dijet kinematics in the next iteration.
- Store separate muon ID and isolation SF factors (currently combined as `weightMuonIDISO`), and clarify/verify trigger SF combination. Store b-tag SF and other standard event weights separately; current fields are an initial subset.
- Expose region/category and electron/b-veto flags. Current `passAdditionalLeptonVeto` is broader than those separately requested decisions.
- Confirm appropriate event HLT paths for Run2024H; initial event branch uses `HLT_IsoMu24[0]`. Object matching is separately computed by the existing passMuonHLT implementation. Audit whether HLT path OR and filter bits match the approved 2024 menu.
- Confirm Golden JSON selection and list of calibration files and their applicability for 2024H (currently inherited). Ensure all object arrays used in 2024 NanoAOD have correct types/branches.
- Audit KIT `kScaleOnly` (already configured by default), beam-constrained input base pT and both corrections in current analysis with data/MC reference values. `muonBasePt` and `muonBsConstrainedPt` are saved to support that comparison.
- The pair definition currently chooses first passing pair by original index; synchronize candidate arbitration with other groups and clarify how to handle 3+ muon events.
- The broader-MPV Landau bias investigation remains paused.

Do not upload or represent this prototype as the final synchronized ntuples.
