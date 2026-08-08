# Experiment 04 — Does the measurement protocol create the ESM-2 within-assay depth×scaling result? (PREREGISTERED)

**Status: protocol only. No Experiment 04 confirmatory analysis has been run.
This file is committed before any 04 regression.** The benchmark composition
below — sequence lengths, windowing exposure, and its association with depth —
was inspected first. That is the cheap metadata diagnostic STATUS.md authorized;
it is descriptive, not the confirmatory statistic. The confirmatory quantity, the
fixed-effects `lp:ld` coefficient on the context-safe subset, has not been
computed. Composition was seen; the coefficient was not. This prediction is made
with the composition known and says so, rather than pretending to be blind.

## Estimand

Do model-context / windowing rules induce a depth-dependent distortion in the
published ProteinGym scaling estimates — enough to create the Experiment 01 ESM-2
within-assay depth×scaling interaction?

Not "does truncation hurt the 15B model on long proteins." Experiment 01's
anatomy rules that framing out before we start: the interaction is carried by the
650M->3B contrast (beta = -0.0205, excludes zero), while 3B->15B is weaker and
inconclusive (beta = -0.0142, includes zero). A 15B-specific truncation theory is
already misaligned with where the effect lives. The estimand is the general one:
any size-differential scoring distortion that tracks depth.

**Primary hypothesis.** The ESM-2 within-assay depth×scaling interaction does not
depend materially on the assays exposed to long-sequence windowing.

## Why this experiment, now

The whole Experiment 01 arc rests on one number: the ESM-2 depth×scaling
interaction survives assay fixed effects (beta = -0.0153, p = 0.011). It is a
within-assay result — identified purely from how each assay's Spearman moves
across the ladder, with everything constant within an assay absorbed. One class of
confound is not absorbed: a scoring-protocol effect that varies across model size
*within* an assay. STATUS.md ranks this threat 1, above ProGen ladder power
(threat 2), because it sits under the result the arc depends on. ESM-2 is the
primary target for the same reason. ProGen is a secondary cross-architecture
diagnostic: its result is already unresolved on power grounds, and its signal is
between-assay, which a within-assay scoring confound would not produce.

## Stage 1 — code audit (locked, no outcome analysis)

Pin the exact scoring implementation at ProteinGym commit `144fe22b`, for every
ESM-2 checkpoint used in Experiment 01 and, secondarily, ProGen2/3. Record
effective context length, window rule, overlap, edge treatment, per-position
aggregation, and — the load-bearing question — whether any of these vary by model
size. Derive every threshold and exposure variable from the implementation, not
from memory. No assay is classified and no coefficient is fit in this stage.

Findings already pinned for ESM-2, `proteingym/baselines/esm/compute_fitness.py`
with `utils/scoring_utils.get_optimal_window`:

- `model_window = 1024` tokens (about 1022 residues plus BOS/EOS), **identical
  across 650M, 3B, and 15B.** Windowing is therefore an assay-level constant
  across the ladder; it cannot enter `lp:ld` through different windows per size.
  The only admissible mechanism is size-differential distortion on the same input.
- protein <= 1022 residues -> window `[0, seq_len]`. **No windowing; the complete
  protein is scored.** Because the threshold is size-invariant, "context-safe" is
  checkpoint-invariant for ESM-2 — an assay safe for one size is safe for all.
- longer protein, interior mutation -> 1024-token window centered on the mutation
  `[mut-512, mut+512]`; near a terminus -> the terminal 1024 block.
- `--scoring-window` default `optimal` (the centered rule above); the
  `wt-marginals` long-sequence branch only executes the overlapping tiling when
  `--scoring-window overlapping`, while `masked-marginals` supports the centered
  `optimal` window. The default pair therefore does not by itself identify the
  published long-sequence configuration. Which strategy produced the published
  ESM-2 scores must be read from the scoring config and recorded verbatim (open
  item). It does not affect the context-safe subset: every strategy scores a
  <=1022 protein identically.

ProGen2/3 use a different, autoregressive scoring path and their own context
limits; audit and record separately, deferred as secondary.

## Stage 2 — metadata-only falsification (CPU, cheap)

**Exposure variable, derived from the pinned rule.** Per assay,
`exposure = max(0, L - 1022) / L` — the fraction of the target protein that cannot
fit a single window. Continuous, not a long/short label. A per-mutation refinement
(fraction of the protein outside each mutation's centered window) is available if
the DMS files are pulled; not required for the falsification, which turns only on
the binary context-safe split.

**Context-safe subset.** Assays whose target protein is <= 1022 residues, so the
window rule never fires at any checkpoint. On these assays every checkpoint
receives the complete target sequence, so the estimated coefficient is free of the
context-windowing / truncation mechanism under test. It is not thereby free of
every conceivable size-dependent scoring effect — numerical precision,
checkpoint-specific implementation, score extraction — which assay fixed effects
do not rule out and this experiment does not claim to rule out. The scope is
deliberately narrow: this stage tests the long-sequence-handling mechanism
STATUS.md named, not an unbounded ProteinGym fairness audit.

**Composition (already inspected, descriptive).**

| Set | Assays |
|---|---|
| context-safe (<= 1022 residues, exposure = 0) | **201 / 217 (92.6%)** |
| exposed (> 1022 residues, exposure > 0) | 16 / 217 (7.4%) |

The 16 exposed assays are 8 Low + 8 Medium depth, **0 High**, and collapse to
**13 distinct proteins** (CAR11, NPC1, SPIKE each carry two assays), so the exposed
side is small and not independent. Exposure is *identically zero* for every one of
the 72 High-depth assays and anti-correlated with depth (length vs Neff/L Spearman
-0.62). Long proteins are the shallow ones.

**Direction of the composition argument — stated carefully.** Under the mechanism
originally proposed, context restriction preferentially degrading larger models,
this composition predicts a bias toward zero or the opposite sign in the full-set
estimate, not inflation of the observed negative interaction: the distortion would
land on shallow assays and flatten the shallow-versus-deep scaling gap. We do not
assume windowing's effect on a Spearman must be degradation, or monotonic in model
size — a restricted context could in principle improve one checkpoint's ranking
more than another's, which on shallow assays could move `lp:ld` in either
direction. Phase 1 tests the effect empirically rather than assuming its
direction.

**Tests.**
1. Association: is exposure related to MSA depth? Report it. Direction informs
   interpretation; the association itself is only opportunity for confounding, not
   evidence of it.
2. Falsification (decisive): fit `rho ~ lp + lp:ld + C(assay)`, cluster SE by
   `UniProt_ID`, upper segment (650M/3B/15B), on (a) all 217 and (b) the 201
   context-safe assays — exactly the Experiment 01 estimator and clustering, only
   the assay set changes.
3. Confound term: add `exposure:lp:ld` on the full set; test whether the carrier
   is in the exposed assays.
4. Exploratory diagnostic (not confirmatory): the same estimator on the 16 exposed
   assays alone. With 8 Low + 8 Medium + 0 High it has restricted depth support;
   its job is to flag whether those assays behave strangely, not to establish a
   depth-scaling law independently.

## Primary prediction and decision (stated to fail)

Reference: the Experiment 01 upper-segment fixed-effects estimate,
`beta_full = -0.0153`.

| Outcome | Conclusion |
|---|---|
| Context-safe retains \|beta_unwindowed / beta_full\| >= 0.75, same sign | scoring-protocol explanation strongly weakened; stop before GPU unless another preregistered trigger fires. (>= 0.75 = <= 25% attenuation = the repo's existing "strong survival" rubric, 01b.) |
| Retains 0.50–0.75 | possible partial contribution. (New tier, a convention introduced for this experiment, not a prior repo rule.) |
| Retains < 0.50, or beta ~ 0, or sign reversal | protocol artefact remains live; controlled re-scoring required (Stage 3). |
| Harmonized re-scoring restores >= 0.50 of the lost coefficient | scoring protocol materially explains Experiment 01. |
| Re-scoring moves assay scores but attenuates beta < 25% | protocol affects scores but does not explain depth×scaling. |
| Too few independent exposed / context-safe assays | metadata stage underpowered; no causal claim. Fires on the *exposed* side (13 proteins), not the safe side (201 assays): the falsification stays well-powered while windowed-specific characterization does not. |

Significance is reported but is not the principal threshold. Dropping 16 assays
changes the standard error mechanically, so a p that crosses 0.05 is not itself
failure — the estimand is attenuation of the coefficient, read as the ratio above.

**Secondary directional prediction (not pass/fail).** Under the originally
proposed mechanism the context-safe beta would be equal to or more negative than
`beta_full`. Report the direction; do not gate the decision on it. Removing
observations changes leverage and sampling variance, so a modestly weaker estimate
does not by itself falsify.

## Stage 3 — controlled re-scoring (GPU/weights, conditional)

Triggered only if Stage 2 leaves the mechanism alive (ratio < 0.75, or the
confound term carries the effect). ESM-2 cannot be run at literal full length:
proteins beyond the architecture's context cannot be passed whole, and the
published overlapping path itself constructs multiple valid 1024-token evaluations
rather than extending the context. So Stage 3 compares windowing schemes, not
window versus impossible full context.

- **A (reproduce):** the exact published ProteinGym scoring protocol.
- **B (harmonized alternative):** a single scheme applied identically to all three
  checkpoints — either overlapping tiled scoring with fixed weighting, or
  mutation-centered scoring. Preregister one as primary; run both if compute
  permits.

Re-score a preregistered informative subset deliberately broader than the 16
exposed assays (too few and non-independent — 13 distinct proteins — to support a
windowed-specific coefficient). **Endpoint: attenuation of the same `lp:ld`
coefficient, A versus B**, on identical assays and specification as Experiment 01.
Not whether individual assay Spearman values move — re-scoring can improve many
assays without touching beta; if scores move and beta does not, scoring quality
changed and the Experiment 01 explanation still failed.

## Kill criteria

- Reference `seq_len` or `target_seq` malformed for the exposed set -> cannot build
  exposure reliably; run Stage 2 on the clean subset and mark the rest "not run."
- The exact ESM-2 long-sequence strategy cannot be recovered from the pinned
  config -> Stage 2 falsification (the 201) runs unaffected; Stage 3 is blocked
  until the rule is pinned, and stops rather than assuming a window.

## Outputs

- `notebooks/04_scoring_window_fairness.py`
- `results/04_exposure_manifest.csv` — per-assay `seq_len`, `exposure`,
  context-safe flag, depth category, selection type, `UniProt_ID`
- `results/04_interaction_by_exposure.csv` — full / context-safe estimates and the
  `exposure:lp:ld` confound term
- `results/provenance_04_scoring_window.json` — pinned commit, scoring file and
  line refs, per-checkpoint window rule, strategy used for the reference scores
- `figures/04_exposure_vs_depth.png`
- `figures/04_interaction_by_stratum.png`

## Open item to pin before Stage 3

Which strategy produced the published ESM-2 substitution scores: `wt-marginals` +
`overlapping` (edge-tapered tiling) or `masked-marginals` + `optimal`
(mutation-centered window). Read it from the scoring launcher at the pinned commit;
record verbatim. Stage 2 does not depend on it; Stage 3 cannot proceed without it.

## What this cannot establish

Stage 2 clears the ESM-2 within-assay result of the context-windowing mechanism on
the 92.6% of assays that fit in context — which includes every deep family, where
the effect lives — but says nothing about the 16 exposed assays beyond removing
them, and does not rule out size-dependent scoring effects other than windowing.
It does not address threat 2, short-ladder power for ProGen, which passes to
Experiment 05. And removing the windowing mechanism is not removing model
misspecification: the interaction on the context-safe set is a real property of
published ESM-2 scores, still on one benchmark, still not a general scaling law.
