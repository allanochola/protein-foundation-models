# Experiment 05 — Is the ProGen within-assay null architecture or power? (PREREGISTERED)

**Status: protocol only. No Experiment 05 confirmatory statistic has been run.
This file is committed before any 05 regression, injection simulation, or
sub-ladder fit.** The ladder geometry below (checkpoint counts, log-decade
spans, which sizes ProteinGym already scored) is metadata inspected first and
recorded as descriptive, exactly as Experiment 04 recorded benchmark
composition before its confirmatory coefficient. Composition is seen; no
within-assay coefficient on any new assay set or ladder is computed here.

## The question

Experiment 01 left one live question, and 02–04 closed everything around it.
The ESM-2 depth×scaling interaction is within-assay: it survives assay fixed
effects (beta = -0.0153, p = 0.011). In both ProGen ladders the same
interaction replicates on the confirmatory model and **collapses under fixed
effects** — the signal is between-assay only. Two hypotheses remain, and they
are not artefact-rejection candidates. They are live scientific alternatives:

- **H_arch — architecture.** ProGen genuinely lacks the within-assay
  depth×scaling response. The between-assay-only signal is a real, different
  behaviour from ESM-2's.
- **H_power — power.** ProGen has the within-assay effect too, but the short
  upper-segment ladders cannot identify it under fixed effects, where
  within-assay leverage is scarcest.

This is the first experiment in the series whose job is positive
discrimination between two hypotheses, not removal of a dead one. A null result
does not settle it; the design must make each hypothesis win in a way the other
cannot fake.

## The design premise STATUS.md got wrong, and why it makes 05 cheap

STATUS.md asserts "ProteinGym cannot lengthen the ProGen ladders." That is true
only for *upward* extension of the *upper segment* — there is no released ProGen2
checkpoint above 6.4B, none for ProGen3 above 3B. It is false for the ladders as
a whole. ProteinGym at pinned commit `144fe22b` already scored every released
size of both models on the same 217 assays, and Experiment 01c verified this at
its gate. The confirmatory fits then **discarded the lower checkpoints** by
restricting to the upper segment (764M+ for ProGen2, 762m+ for ProGen3), by
analogy to ESM-2's 650M breakpoint.

Those discarded checkpoints are unused identification leverage that costs no GPU
to recover:

| Ladder | Upper segment (used in 01c) | Full ladder (already scored) |
|---|---|---|
| ESM-2 | 650M/3B/15B — 3 pts, 1.36 dec | 8M/35M/150M/650M/3B/15B — 6 pts, 3.27 dec |
| ProGen2 | 764M/2.7B/6.4B — 3 pts, 0.92 dec | 151M/764M/2.7B/6.4B — 4 pts, 1.63 dec |
| ProGen3 | 762m/1b/3b — 3 pts, 0.60 dec | 112m/219m/339m/762m/1b/3b — 6 pts, 1.43 dec |

The number the fixed-effects collapse was blamed on — short ladder — is partly a
choice, not a benchmark limit. ProGen3's full ladder carries six checkpoints
over 1.43 decades: more points and nearly the span of the ESM-2 upper segment
(1.36 dec) where the within-assay effect *survives*. If power is the whole story,
that leverage should be enough to recover it. The heavy GPU step STATUS.md
sequenced last may not be needed at all.

## Stage 1 — audit (locked, no outcome analysis)

Confirm before any fit, derive from the data not from memory:

1. **Full-ladder scores load at the pinned commit.** Verify ProteinGym
   `144fe22b` supplies substitution scores for ProGen2 {151M, 764M, 2.7B, 6.4B},
   ProGen3 {112m, 219m, 339m, 762m, 1b, 3b}, and ESM-2 {8M, 35M, 150M, 650M, 3B,
   15B} across the same 217 assays. Gate: any missing size drops from its
   ladder and is recorded "not scored," not imputed.
2. **The existing pairwise contrasts are confirmatory, not fixed-effects.**
   `notebooks/01c_replication_analysis.py:149` fits the pairwise decomposition
   with `PRIMARY` (covariate model), not `FIXED_EFFECTS`. So the Stage 2
   Analysis C sub-ladder FE fits are genuinely new evidence, not a relabel of
   numbers already in `results/`. Recorded so the novelty claim is auditable.
3. **ProGen3 provenance and weight availability.** Record where the existing
   ProGen3 scores came from (ProteinGym column vs Allan-scored) and whether
   Profluent-AI/progen3 weights are public, because only the *conditional*
   Stage 3 needs weights; Stage 2 needs none. Non-blocking for Stage 2.

No coefficient is fit in Stage 1.

## Stage 2 — the leverage already in the benchmark (CPU, decisive)

Three analyses, each attacking the architecture/power split from a different
side. They are built so the two hypotheses win by different routes: **H_power
wins on a positive** (a real within-assay effect appears once leverage is
restored), **H_arch wins on an informed null** (ProGen stays null on a ladder
demonstrably powerful enough, while a matched real effect survives the same
handicap). All three reuse the Experiment 01 estimator verbatim —
`rho ~ lp + lp:ld + C(assay)`, cluster-robust SE by `UniProt_ID` — changing only
the ladder or the data source.

### Analysis A — restore the discarded leverage (already-scored lower checkpoints)

Two fits per model, both on data already on disk.

**A1 — upper-segment extent sensitivity.** The 762m ProGen3 breakpoint was set by
analogy to ESM-2, never fitted. Lower it using the scored sub-762m checkpoints
and re-fit the FE interaction as the upper segment grows:

- ProGen3: {762m/1b/3b, 0.60 dec, 3 pts} → {339m+, 0.95 dec, 4 pts} →
  {219m+, 1.13 dec, 5 pts} → {112m+, 1.43 dec, 6 pts}.
- ProGen2: {764M+, 0.92 dec, 3 pts} → {151M+, 1.63 dec, 4 pts}.
- ESM-2 anchor: {650M+, 1.36 dec, 3 pts} → {150M+, 2.0 dec, 4 pts} →
  {35M+, 2.6 dec, 5 pts} → {8M+, 3.27 dec, 6 pts}.

**A2 — full-ladder global FE fit, plus a segmented low/high decomposition.** Fit
the global `lp:ld` on the complete ladder for each model, and separately fit a
piecewise specification that estimates the depth×scaling slope below and above
the breakpoint at once: `rho ~ lp_lo + lp_hi + lp_lo:ld + lp_hi:ld + C(assay)`,
where `lp_lo`, `lp_hi` are the hinge terms at the model's breakpoint. The
coefficient on `lp_hi:ld` isolates the upper-regime interaction while the lower
checkpoints anchor the assay means and the low-regime slope.

The asymmetric reading is preregistered and load-bearing, but it is stated as an
inferential rule, not a mathematical guarantee:

- A **significantly negative** ProGen `lp:ld` on an extended or full ladder is
  strong evidence for **H_power** — *provided the segmented fit shows the negative
  slope is carried by the upper regime, not the lower one.*
- A **null** on the extended/full ladder is **ambiguous** — genuine absence
  (H_arch) or low-regime dilution — and does not on its own support H_arch.

I am not claiming extension can only attenuate. Adding lower checkpoints shifts
the within-assay lp variance, the lp centering, and the lp:ld covariance under
fixed effects, and those can move the coefficient in either direction; the
convex-combination-of-segment-slopes picture holds only under a flat, zero-
interaction low regime, which is assumed, not established. The lower segment's
interaction sign was never characterised — 01a defined the effect above the
breakpoint and said nothing below it. So a recovered negative global slope could
in principle originate entirely in the low regime and have nothing to do with the
upper-regime null. The segmented decomposition is the guard: a global negative
counts as H_power evidence only when `lp_hi:ld` carries it and `lp_lo:ld` does
not. If the low regime turns out to carry its own negative interaction, that is a
separate, reportable finding, not confirmation of the upper-regime effect.

The ESM-2 anchor calibrates whether extension preserves a real effect at all. If
ESM-2's known within-assay interaction survives its own downward extension —
global or, more likely, in `lp_hi:ld` — extension preserves real effects and a
ProGen null is informative. If ESM-2's effect itself washes to null once the low
regime is added, extension destroys real effects and neither ProGen null means
anything; A is uninformative and the verdict rests on B and C.

### Analysis B — did the upper design ever have the power to see it? (injection)

Quantify detectability directly. On each ProGen ladder's *actual* geometry —
upper segment and full — inject a known within-assay `lp:ld` slope, regenerate
synthetic Spearman values, refit the FE estimator, and count detections.

Construction, per replicate:
1. Fit a **residual-source model** on the real ProGen data; keep its fitted
   values and residuals.
2. Synthetic outcome = fitted value + `beta_inject * (lp_c * ld_c)` +
   residual, where `lp_c`, `ld_c` are the centred regressors on the segment being
   fitted and the residual is **block-resampled by assay** (all sizes of a
   sampled assay move together) to inherit ProGen's real per-observation noise,
   heteroskedasticity, and the 3-plus-per-assay correlation. Semi-parametric, not
   Gaussian — the noise model is empirical.
3. Refit `rho ~ lp + lp:ld + C(assay)`, cluster by `UniProt_ID`. Record the
   recovered `lp:ld`, its CI, and whether it is negative with p < 0.05.

**Residual-source model — the contamination guard, run as a bracket.** Residuals
carry whatever the source model omits. Residualising on `rho ~ C(assay)` leaves
ProGen's *real* within-assay structure — the lp main effect and, critically, any
real `lp:ld` — inside the residuals, so block-resampling reinjects it and the
synthetic slope is `beta_inject` plus ProGen's real (unknown) slope. That is
exactly the effect being measured, so the power estimate is contaminated by the
answer. Note that `rho ~ assay + ld` does not help — `ld` is assay-constant and
fully absorbed by `C(assay)` — and `rho ~ assay + lp` removes only the main
effect, leaving the interaction. Removing the interaction structure requires a
source model that contains `lp:ld`. So the bracket is three source models:

- `rho ~ C(assay)` — **primary (preregistered).** Retains real structure as
  noise, so it *over*-states residual variance and *under*-states power. It is
  the conservative choice for the risky conclusion: it makes an informed null
  (H_arch) harder to reach, not easier.
- `rho ~ lp + C(assay)` — removes the lp main effect.
- `rho ~ lp + lp:ld + C(assay)` — removes the interaction too; least
  contaminated, so it *over*-states power. The ceiling of the bracket.

The primary and the interaction-removed model bound the truth from opposite
sides. If the verdict (which side of the MDE thresholds below the recovered power
lands) is stable across the bracket, contamination is immaterial. If it flips,
the power estimate is fragile and is reported as a bracket, not a point.

`beta_inject` grid, preregistered: **0** (calibration), **-0.0153** (the ESM-2
within-assay magnitude), each ProGen's **own confirmatory magnitude** (-0.0245
ProGen2, -0.0283 ProGen3, as a plausible upper bound), and a sweep -0.004 to
-0.030 in 0.002 steps for a power curve dense enough to read a minimum detectable
effect. 2000 reps per cell, seed 0.

Calibration kill criterion: injecting 0 must yield a false-positive rate in
[0.03, 0.07] under the primary source model. Outside that, the simulator is
misspecified and B is discarded rather than reported.

**Reading — bounded, not binary.** The output is not "architecture or power." It
is the **minimum detectable effect (MDE)**: the smallest `|beta_inject|` the
ProGen design detects at 80% power on its actual geometry. Power >= 0.8 at
-0.0153 on the **upper segment**, with the real upper-segment fit null, licenses
only *"ProGen has no within-assay effect as large as ESM-2's"* — it does not
exclude a smaller real effect, e.g. -0.007, which is a genuine version of
H_power the -0.0153 injection cannot see. The honest informed-null claim is
therefore quantitative: *the within-assay effect, if present, is smaller than the
MDE.* The full-ladder power cell gives the MDE the recovered leverage achieves,
and predicts whether A2 should have been able to resolve anything.

### Analysis C — does a real effect survive a ProGen-sized handicap? (ESM-2)

The mirror of B, run on real data instead of synthetic. Take ESM-2, where the
within-assay effect is confirmed, and cripple its ladder to ProGen leverage by
re-fitting the FE estimator on sub-ladders:

- 650M/3B/15B — 1.36 dec, 3 pts (reference, must recover -0.0153).
- 650M→3B — 0.66 dec, 2 pts (**overshoots the handicap**: shorter than ProGen3
  upper's 0.60 dec is not reachable with ESM-2 checkpoints; 0.66 is the closest
  and slightly harder than ProGen2 upper).
- 3B→15B — 0.70 dec, 2 pts.
- 650M→15B — 1.36 dec, 2 pts (full span, minimum points: isolates the
  point-count handicap from the span handicap).

These FE 2-point fits are new (Stage 1 item 2). If ESM-2's real within-assay
effect **survives** a 0.66-decade 2-point FE fit — leverage at or below ProGen3
upper's — then short span alone does not kill a real within-assay effect, and the
ProGen null is harder to blame on span: leans **H_arch**. If ESM-2's effect
**collapses** at that handicap, the ProGen upper null is exactly what an
underpowered real effect looks like: leans **H_power**.

**C is corroboration, not a co-equal pillar, and the protocol treats it that
way.** A 2-point FE interaction is identified from a single within-assay contrast
per assay — the endpoint difference — so it is high-variance by construction, and
ESM-2's checkpoint spacing forbids a short-span *3-point* sub-ladder (every
consecutive ESM-2 triple spans ~1.3 decades), so C cannot be made quantitatively
rigorous no matter how it is cut. It is a demonstration argument — *a real effect
can or cannot survive this handicap* — not a precise estimate. To keep the
instability visible rather than asserted, each 2-point fit is bootstrapped by
assay cluster (2000 reps, seed 0) and reported with its bootstrap interval, so
"survives" or "collapses" is read off a distribution, not a single fragile point
estimate. C moves the verdict only at the margin, breaking ties between A and B;
it never overrides them.

## Primary prediction and decision (stated to fail)

Reference: ESM-2 upper-segment fixed-effects `beta = -0.0153`, and each ProGen's
upper-segment FE null from 01c (ProGen2 +0.0064 p=0.60; ProGen3 -0.0031 p=0.77).

The pre-stated expectation is **H_power** — that restored leverage recovers a
negative ProGen within-assay slope. Stating the bet first, as with 03's
composition prediction, so a null is a real falsification and not a quiet
reinterpretation.

Verdicts are stated as bounded claims about effect *magnitude*, not as an
unqualified architecture/power dichotomy. "H_arch" everywhere below means the
specific, falsifiable statement *ProGen has no within-assay effect as large as
the MDE* — never "no effect at all."

| Analysis A (extended, segmented) | Analysis B (MDE on upper seg) | Analysis C (ESM-2 handicap) | Verdict |
|---|---|---|---|
| ProGen `lp_hi:ld` sig. negative, carried by upper regime, ESM-2 anchor survives | any | any | **H_power** — effect real, upper ladder hid it. Confirmed, no GPU. |
| ProGen null, ESM-2 anchor survives | MDE <= -0.0153 (>= 0.8 power at ESM magnitude) | ESM-2 survives 0.66-dec handicap | **H_arch, bounded** — no within-assay effect as large as ESM-2's; a smaller effect is not excluded. State the MDE as the bound. No GPU. |
| ProGen null, ESM-2 anchor survives | MDE >= -0.0153 (<= 0.5 power at ESM magnitude) | ESM-2 collapses at handicap | **H_power live, unresolved** — upper design cannot see an ESM-sized effect; Stage 3 warranted. |
| ProGen global negative but carried by `lp_lo:ld` | any | any | **Low-regime finding, not H_power** — the upper-regime null stands; report separately. |
| ESM-2 anchor itself washes to null under extension | (B, C stand alone) | (B, C stand alone) | A uninformative; decide on B+C; if they split, Stage 3. |
| A, B, C point different directions | — | — | Unresolved; Stage 3, and report the split rather than force a call. |

Significance is reported throughout but is not the sole threshold; changing the
ladder changes the standard error mechanically. Alongside every fit, report the
effect-size ratios `beta_recovered / beta_ESM` (-0.0153) and
`beta_recovered / beta_confirmatory` (each ProGen's 01c between-assay estimate),
each with a bootstrap CI from the assay-cluster resample, mirroring the retention
ratios the repo already uses in 01b and 04. A ratio with a CI is the primary
read; recovery-or-not against the 01c null is secondary; the MDE bound is what a
null actually buys.

## Stage 3 — fresh inference (GPU/weights, conditional, possibly never)

Triggered only by the last three rows above. Defined *after* Stage 2 names the
specific residual ambiguity, to keep scope bounded — as Experiment 04 gated its
Stage 3 on the metadata outcome. Because the audit shows ProteinGym already
scored every released ProGen size, "lengthen the ladder with more of the same
model" adds nothing; Stage 3 is therefore not ladder extension but confound
breaking:

- The architecture signal is perfectly confounded with ladder length across
  01/01c — ESM-2 (within-assay) is the long masked-LM ladder; ProGen
  (between-assay) is the short autoregressive one. Break it by scoring a
  **masked-LM on a deliberately short ladder** or an **autoregressive model on a
  long, dense ladder** on the same 217 assays, fixed scoring protocol across
  checkpoints, and testing whether the within-vs-between split tracks
  architecture family or ladder geometry.
- Preregister the specific model and checkpoint set at trigger time, not now.
  Do not pre-commit compute.

## Kill criteria

- Full-ladder scores absent for a size at the pinned commit → that size drops,
  its ladder shortens, recorded "not scored." Analysis proceeds on what loads.
- Injection calibration outside [0.03, 0.07] false-positive at `beta_inject = 0`
  → Analysis B discarded, not patched; verdict rests on A and C.
- ProGen3 provenance unrecoverable or weights unavailable → affects only a
  possible Stage 3 arm; Stage 2 unaffected.

## Outputs

- `notebooks/05_architecture_vs_power.py` — Stage 2 A/B/C, parameterised by ladder
- `results/05_ladder_extent_sensitivity.csv` — A1 extended-segment fits and A2 global + segmented (`lp_lo:ld`, `lp_hi:ld`) fits per model, each with pts, log-span, and the ratios `beta/beta_ESM` and `beta/beta_confirmatory` with bootstrap CIs
- `results/05_injection_power.csv` — B detection rates per ladder × `beta_inject` × residual-source model, the 0-injection calibration row, and the derived MDE per ladder
- `results/05_esm2_handicap.csv` — C sub-ladder FE fits with span, point count, and bootstrap interval
- `results/provenance_05_architecture_vs_power.json` — pinned commit, sizes loaded per model, seed, residual-resampling scheme and source-model bracket
- `figures/05_lp_ld_by_ladder_length.png` — central figure: within-assay `lp:ld` vs ladder span/points, all three models, ProGen extension trajectory against the ESM-2 anchor
- `figures/05_injection_power_curve.png` — power vs `beta_inject` per ladder, MDE marked, bracket band across residual-source models

## What this cannot establish

Even a clean H_arch verdict is a claim about ProteinGym-scored zero-shot
variant-effect prediction on 217 assays, not a general statement about
autoregressive protein language models. The injection power calculation assumes
the only ESM-2/ProGen difference relevant to detectability is ladder geometry and
noise; if ProGen's scoring is biased in some size-dependent way not captured by
the resampled residuals, B's power estimate is optimistic. Analysis A trades
estimand cleanliness for leverage by construction — its null is deliberately weak
evidence and is treated as such. And distinguishing architecture from power does
not make the surviving effect a scaling law: it stays a benchmark-specific,
depth-dependent property of published scores, now attributed to one model family
or to measurement resolution, on one benchmark.
