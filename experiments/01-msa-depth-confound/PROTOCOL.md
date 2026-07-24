# Experiment 01 — MSA depth as a confound in pLM scaling

Written before the analysis, amended after. Amendments are marked and the
original wording is left intact.

## Question

Do reported ESM-2 scaling gains on zero-shot variant effect prediction
survive a control for how well sampled the protein's family is?

## Original prediction

Scaling gains concentrate in high-MSA-depth families. Under depth control the
scaling curve flattens.

## Status: prediction failed, different effect found

**Amendment, pilot run.** The prediction is not supported. Across the full
8M-to-15B ladder the gain is roughly constant across depth strata
(Low +0.197, Medium +0.242, High +0.197 mean Spearman) and the size x depth
interaction is null (beta = -0.006, p = 0.30).

The effect is real but lives in the upper half of the ladder. Restricting to
650M and above, the interaction is significant (beta = -0.015, p = 0.002)
while the main effect of size is not distinguishable from zero
(beta = -0.010, p = 0.11). Paired per assay:

| MSA depth | mean 15B − 650M | assays where 15B wins |
|---|---|---|
| Low | +0.009 | 61% |
| Medium | +0.007 | 49% |
| High | −0.054 | 35% |

**Revised claim.** Past 650M, scaling ESM-2 degrades performance on
well-sampled families and helps marginally on poorly-sampled ones. The
aggregate number hides a sign flip. Reporting a single mean Spearman across
ProteinGym conceals that the largest model is the worse choice for most
assays in the benchmark.

The SaProt paper noted in passing that ESM-2 15B fails to beat 650M. This
locates where and by how much.

## Data

- ProteinGym per-assay zero-shot Spearman, DMS substitutions, 217 assays.
- ProteinGym reference file for `MSA_Neff_L` (continuous) and
  `MSA_Neff_L_category` (Low 36 / Medium 109 / High 72).
- Both pulled from the ProteinGym repo at run time. Nothing is vendored.

## Method

OLS of Spearman on centred log10 parameters, centred log10 Neff/L, and their
interaction. Standard errors clustered by `UniProt_ID`, since several assays
share a protein. Paired per-assay differences reported separately because
means hide direction.

## Kill criteria

- Interaction null in both the full ladder and the upper segment → confound
  story is wrong for this task, stop. **Not triggered.**
- Effect disappears once `taxon`, `seq_len`, and `selection_type` are added →
  the pattern is a composition artefact, not a depth effect. **Run (01b), not
  triggered.** See below.
- Pattern absent on the ProGen2 ladder → ESM-specific, much weaker claim.
  **Not yet run.**

## Amendment, 01b — covariate controls

The confound is real: prokaryotic assays average roughly 0.5 log units more
Neff/L than human ones, and sequence length varies by taxon too. So the
control was necessary, not a formality.

The interaction survives it. Upper segment, 650M and above:

| Specification | beta | p | 95% CI |
|---|---|---|---|
| 1. baseline (01a) | −0.0153 | 0.002 | [−0.025, −0.006] |
| 2. covariates as main effects only | −0.0153 | 0.002 | [−0.025, −0.006] |
| 3. size × covariate interactions | −0.0172 | 0.010 | [−0.030, −0.004] |
| 4. assay fixed effects | −0.0153 | 0.011 | [−0.027, −0.004] |

Across the full 8M–15B ladder it stays null under every specification
(p = 0.14 to 0.34), as in 01a.

**Two things worth recording.**

Specification 2 changed nothing — beta identical to four decimal places.
Assay-level covariates shift levels; they cannot absorb an interaction with a
variable that varies within assay. "We controlled for confounds" by adding
main effects is a control that does no work, and it is the standard move.

Specification 4 is the real test. Assay fixed effects absorb every
assay-level property, measured or not, and identify the interaction purely
from how each assay's score moves across model sizes. `ld` drops out as
collinear. Nothing about which proteins are in the benchmark can produce this
result.

**Cost of the honest version.** p rises from 0.002 to 0.011. Significant at
0.05, not at 0.01. Reporting the baseline alone would have overstated the
evidence.

### Robustness battery (R1-R6)

Full table in `results/01b_robustness_checks.csv`. Sign is negative in all 14
fits. Significance is not.

| Check | beta | attenuation | p |
|---|---|---|---|
| confirmatory (cluster UniProt) | −0.0172 | −13% | 0.010 |
| R1 HC3 | −0.0172 | −13% | **0.285** |
| R2 cluster by assay | −0.0172 | −13% | 0.004 |
| R3 drop Human (−44% of data) | −0.0166 | −8% | **0.087** |
| R4 drop OrganismalFitness | −0.0107 | **+30%** | **0.128** |
| R4 drop Stability | −0.0202 | −32% | 0.003 |
| R5 drop Cook's D > 4/n (43 obs) | −0.0133 | +13% | **0.163** |
| R6 assay fixed effects | −0.0153 | 0% | 0.011 |

**Verdict: partial survival.** Not a composition artefact — the sign never
flips and the fixed-effects estimate is unattenuated. But R5 is the finding
that matters: removing 6.6% of observations by influence costs 13% of the
magnitude and all of the significance, so a small number of assays carry the
effect. R4 agrees from another angle. Any writeup leads with this, not with
the confirmatory p-value.

### Process failure, recorded

01b was run in the same session it was designed. The kill criterion was
preregistered; specifications 3 and 4 were added after seeing 1 and 2. They
moved the estimate against the hypothesis, so nothing was gained — but that
is unverifiable from outside, which is the point. 01c is preregistered and
committed before any 01c code is written.

## What is not yet established

The pilot is a pattern, not a result. MSA depth is not randomly assigned:
deep families skew toward well-studied human and viral proteins with
particular assay types. Before this is publishable it needs the covariate
controls, the ProGen2 replication, and a check on whether ProteinGym's own
scoring protocol treats the large models fairly (context truncation, batching
differences).

ProteinGym already reports mean performance binned by MSA depth. The novel
part here is treating depth as continuous, fitting the interaction with size,
and finding the sign flip in the upper segment. Do not re-report their table
and call it a finding.

## Next

1. ~~Add covariates.~~ **Done (01b).** Interaction survives assay fixed
   effects at p = 0.011.
2. Replicate on ProGen2 S/M/L/XL — independent architecture, same assays.
   Confirm parameter counts first: Base and M differ by training corpus, not
   size, so one must be dropped.
3. Check whether ProteinGym's scoring protocol handles the large models
   fairly (context truncation, batching).
4. Only then run models yourself, to extend the ladder beyond what
   ProteinGym has scored.
