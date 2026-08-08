# Experiment 04 (Phase 1) — Results

**The ESM-2 within-assay depth×scaling interaction does not depend on assays
exposed to long-sequence windowing. It is fully present, and marginally stronger,
on the 201 assays ProteinGym never windows. Threat 1 (scoring-protocol fairness,
windowing mechanism) is retired on published scores. Stage 3 re-scoring is not
triggered — no GPU.**

Phase 1 was metadata-only and reused the Experiment 01 estimator verbatim
(`rho ~ lp + lp:ld + C(assay)`, cluster-robust SE by `UniProt_ID`, upper segment
650M/3B/15B). Only the assay set changed. Preregistered in PROTOCOL.md, frozen
before the confirmatory coefficient was computed; benchmark composition was
inspected first and recorded as such.

## Confirmatory result

| Set | beta (lp:ld) | p | 95% CI | assays / clusters |
|---|---|---|---|---|
| full (217) — sanity, recovers 01 | -0.0153 | 0.011 | [-0.0270, -0.0036] | 217 / 186 |
| context-safe (<= 1022 residues) | **-0.0180** | 0.009 | [-0.0314, -0.0045] | 201 / 173 |
| windowed-only (16) — exploratory | +0.0033 | 0.836 | [-0.0276, +0.0342] | 16 / 13 |

Retention |beta_safe / beta_full| = **1.18**, same sign. Preregistered rubric
(pass >= 0.75): **PASS**. The full-set fit reproduces the published -0.0153 /
p = 0.011 exactly, confirming the estimator is faithful before the subset test is
read.

Confound term (`lp:ld:windowed`): +0.021, p = 0.20 — the windowed shift points
toward zero and is not significant. The 16 windowed assays carry no interaction
on their own (beta ~ 0, p = 0.84), as expected from 8 Low + 8 Medium + 0 High:
no depth contrast to support one.

## What this establishes

The context-windowing / long-sequence-handling mechanism cannot be the source of
the Experiment 01 result. The interaction is identified entirely on the 92.6% of
assays where every checkpoint receives the complete target sequence, so no
size-differential windowing operates. The composition predicted this: exposure is
zero for all 72 High-depth assays and anti-correlated with depth, so the windowed
assays are shallow and, if anything, dilute the effect. Both the confound term and
the near-zero windowed-only fit are consistent with that direction.

Persistence where the mechanism cannot operate is direct falsification, not a
plausibility argument. The claim the arc can now carry: *the interaction persists
when every assay requiring context windowing is removed.*

## What this does not establish

- The strengthening to -0.0180 is the direction the secondary prediction named,
  but it is within noise — CIs overlap and dropping 16 shallow assays shifts
  leverage. The defensible statement is unchanged-to-slightly-stronger, passing
  retention; not that windowing was masking the effect. Significance was reported,
  not used as the threshold.
- Scope is the windowing mechanism only. Other size-dependent scoring effects
  (numerical precision, checkpoint-specific implementation, score extraction) are
  untested and out of scope by design. Assay fixed effects remove assay-constant
  artefacts, not every model-size-dependent implementation effect.
- Nothing here addresses threat 2 (ProGen short-ladder power), which passes to
  Experiment 05.

## Open item (now non-blocking)

Which strategy produced the published ESM-2 long-assay scores (`wt-marginals` +
`overlapping` vs `masked-marginals` + `optimal`) was not pinned. It is moot for
this conclusion because Stage 3 is not triggered; pin it only if the 16 windowed
assays are ever characterised for their own sake.

## Outputs

- `notebooks/04_scoring_window_fairness.py`
- `results/04_exposure_manifest.csv`
- `results/04_interaction_by_exposure.csv`
- `results/provenance_04_scoring_window.json`

## Next

Experiment 05 — attack the remaining cross-architecture / ladder-power question
(ProGen2/3), threat 2. Experiment 04 Stage 2/3 (controlled re-scoring) remains
available but is not warranted on this result.
