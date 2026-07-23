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
  the pattern is a composition artefact, not a depth effect. **Not yet run.**
- Pattern absent on the ProGen2 ladder → ESM-specific, much weaker claim.
  **Not yet run.**

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

1. Add covariates. Does the upper-segment interaction survive?
2. Replicate on ProGen2 S/M/L/XL. Confirm parameter counts first — Base and M
   differ by training corpus, not size, so one must be dropped.
3. Only then run models yourself (notebook 01b), to extend the ladder to
   models ProteinGym has not scored.
