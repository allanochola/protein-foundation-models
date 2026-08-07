# Experiment 01c — ProGen2 replication (PREREGISTERED, NOT YET RUN)

**Status: protocol only. No 01c analysis has been run. Commit this file
before writing any 01c code.**

01b was run in the same session it was designed. The kill criterion was
preregistered, but two of the four specifications were added after seeing the
first two results. They moved the estimate against the hypothesis, so nothing
was gained by it — but that is unverifiable from the outside, which is the
whole reason preregistration exists. 01c fixes the process.

## Question

Does the negative model-size × MSA-depth interaction found in ESM-2 above
650M reproduce on the ProGen2 ladder — a different architecture, a different
training objective, scored on the same 217 assays?

## Why this is the decisive test

01b established the ESM-2 effect is not driven by benchmark composition. It
did not establish that the effect is a property of scaling rather than of
ESM-2. ProGen2 is autoregressive where ESM-2 is masked, Salesforce-trained on
UniRef/BFD where ESM-2 used UniRef, and shares no weights or code. If the
sign flip reproduces there, the claim generalizes from "ESM-2 does this" to
"scaling pLMs past a certain size does this on well-sampled families."

## Data

Same two files, same pinned commit
`144fe22b07dfaeec2b366f2346203a9838a55b4c`. No new downloads.

ProGen2 columns in the scores file: S, M, Base, L, XL.

**Ladder construction rule, fixed before fitting.** ~~Base and M differ by
training corpus, not parameter count~~ — **verified, and the original wording
was imprecise.** Nijkamp et al. Table 1 gives small 151M, medium 764M, base
764M, large 2.7B, xlarge 6.4B. Medium and base share parameter count, layers,
heads and head dimension; they differ in context length (1,024 vs 2,048),
learning rate, warm-up, and total steps (350k vs 400k). Two points at the same
x with different y would bias the slope regardless of why they differ, so the
rule stands: drop `Progen2 Base`, keep `Progen2 M`.

Verification gate **passed**. Ladder: 151M / 764M / 2.7B / 6.4B.

## Amendment 1 — ProGen3 added as a second replication ladder

**Written and committed before any 01c outcome was inspected.** Recorded
here rather than silently folded in.

While confirming column availability, ProteinGym was found to also carry a
ProGen3 ladder: 112m, 219m, 339m, 762m, 1B, 3B, all scored on the same 217
assays. This was not known when the protocol was written.

ProGen3 is added as a **second preregistered replication ladder**, analysed
identically and reported alongside ProGen2. Two independent replications are
strictly more informative than one, and adding it now — before any result is
seen — costs nothing in inferential validity.

Parameter counts are taken from ProteinGym's own column labels and must be
verified against the ProGen3 paper before fitting, under the same gate that
applied to ProGen2. If they cannot be verified, ProGen3 is dropped and ProGen2
proceeds alone.

Upper segment for ProGen3: 762m and above, by the same analogy to ESM-2's
650M breakpoint that fixed ProGen2's at 764M. Not chosen by inspecting
results.

**This amendment does not change the ProGen2 analysis, the prediction, the
decision rule, or the kill criteria.** Those remain as originally committed.

## Amendment 2 — ProGen3 verification gate, and two limits recorded in advance

**Written before any 01c outcome was inspected.**

Verified against the official Profluent-AI/progen3 release: 112m, 219m, 339m,
762m, 1b, 3b. These match ProteinGym's column labels. **Gate passed.**

Two properties of ProGen3 were not known when Amendment 1 was written and
must be recorded before results are seen.

**Sparse mixture of experts.** ProGen3 activates roughly 27% of parameters
per forward pass. Because the sparsity fraction is constant across the
ladder, total and active parameters stay proportional and log-scaling within
the ladder is unaffected up to an additive constant — the within-ladder
interaction is well defined. Cross-model *magnitude* comparison against dense
ESM-2 and ProGen2 is not well defined, and no such comparison will be made on
magnitude alone.

**Truncated upper segment.** ProteinGym scores ProGen3 only to 3B; the 46B
model is absent. Upper-segment scaling range by ladder:

| Ladder | Upper segment | log10 span |
|---|---|---|
| ESM-2 | 650M – 15B | 1.36 |
| ProGen2 | 764M – 6.4B | 0.92 |
| ProGen3 | 762m – 3B | 0.60 |

Both replications have less leverage than the original. Wider intervals are
expected by construction. A wide interval in ProGen3 is **inconclusive, not
null**, and the decision table already says so.

## Amendment 3 — cross-model comparison is paired, not independent

**Written before any 01c outcome was inspected.**

All three ladders are scored on the same 217 assays. The estimates are
therefore **correlated**, and this is replication across architectures, not
across benchmarks. A shared benchmark artefact — assay composition, scoring
protocol, context truncation — would reproduce in all three and look exactly
like independent confirmation.

Consequences, fixed now:

- Random-effects pooling across models is **not** applied. DerSimonian-Laird
  assumes independent estimates; these are not independent.
- The cross-model comparison is a **paired** contrast on shared assays:
  per-assay model-specific slopes, differenced within assay.
- `cross_model_summary.csv` reports each ladder's beta, CI, and bootstrap
  P(beta < 0) side by side, with the shared-assay caveat stated in the file
  itself, not only in prose.

## Specification

Primary, on the upper segment (models at or above the ESM-2 breakpoint
equivalent — see breakpoint rule below):

```
rho ~ lp + lp:ld + C(assay)          cluster SE by UniProt_ID
```

Assay fixed effects, because they absorb taxon, length, selection type and
any unmeasured assay-level property in one step. Coefficient of interest:
`lp:ld`.

Secondary, reported alongside and not instead:

```
rho ~ lp*ld + lp*taxon + lp*llen + lp*sel
rho ~ lp*ld                          full ladder, all sizes
```

**Breakpoint rule.** ESM-2's flip occurs above 650M. ProGen2's ladder is
151M / 764M / 2.7B / 6.4B. The upper segment is defined as 764M and above —
chosen by analogy to ESM-2's breakpoint, not by inspecting ProGen2 results.
If a different breakpoint later looks better, that is exploratory and
labelled as such.

## Prediction

`lp:ld` is negative in the primary specification. Stated plainly so it can
fail: I expect this to reproduce, based on 01b.

## Decision rule, fixed in advance

Reported for every specification: β, 95% CI, p, and attenuation relative to
the ESM-2 estimate of −0.0153.

| Outcome | Reading |
|---|---|
| Negative, CI excludes 0, \|β\| within 50% of ESM-2 | Replicates — architecture-independent |
| Negative, CI includes 0 | Directionally consistent, underpowered — report as such, do not claim replication |
| Near zero | ESM-2-specific, not a scaling property |
| Positive, CI excludes 0 | Contradicts 01a/01b; the ESM-2 result needs re-examination |

Survival is not defined as p < 0.05. A wide interval is inconclusive, not
null.

## Robustness battery, specified now

Run all of these, report all of them, regardless of what the primary shows:

- HC3 robust SE
- cluster by assay as well as by UniProt
- leave-one-taxon-out (4 fits)
- leave-one-selection-type-out (5 fits)
- drop observations with Cook's D > 4/n
- full ladder as well as upper segment

01b's lesson: R5 removed 6.6% of observations and cost all the significance.
Any 01c result that fails R5 is reported as fragile in the headline, not in a
footnote.

## Kill criteria

- `lp:ld` near zero or positive with a tight CI → the ESM-2 finding does not
  generalize; amend `models/esm2.md` to scope the claim to ESM-2 and stop
  this line.
- Parameter counts cannot be verified → do not fit; the ladder is undefined.

## Multiplicity

Two ladders means two primary tests. No correction is applied, because the
decision rule is not a significance threshold — it reads sign, magnitude and
interval together. Both ladders are reported in full whatever they show, so
there is no selective-reporting channel to correct for. A reader who prefers
a corrected threshold has both results and can apply one.

## What this cannot establish

ProGen2 and ESM-2 were scored by ProteinGym under its own protocol. A shared
artefact of that protocol — context truncation affecting large models,
batching differences — would produce correlated results in both ladders and
look like replication. Ruling that out requires scoring models directly,
which is a separate experiment.
