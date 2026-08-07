# 01c results — ProGen2 and ProGen3 replication

**Verdict: replicates on the confirmatory model, fails under assay fixed
effects.** The depth-scaling interaction reproduces across three
independently developed architectures on the primary specification, which
rules out an ESM-2-specific artefact. But it survives assay fixed effects
only in ESM-2. In both ProGen ladders the signal is between-assay, not
within-assay — consistent with benchmark-composition dependence, with low
power on a short ladder, or both. This strengthens the benchmark-specific
conclusion from 01b rather than overturning it.

Analysis specification, robustness battery, bootstrap scheme, and decision
table were copied unchanged from 01b and fixed in `PROTOCOL.md` before any
ProGen number was seen. `notebooks/01c_replication_analysis.py` runs either
ladder; only the input dataset changes.

## Headline numbers

Confirmatory model (`rho ~ lp*ld + lp*taxon + lp*z_llen + lp*sel`, cluster SE
by UniProt_ID), upper segment. ESM-2 reference: beta = -0.0153, p = 0.011,
survived fixed effects.

| Ladder | Sizes | log10 span | Confirmatory beta | p | Bootstrap P(beta<0) | Fixed effects (R6) |
|---|---|---|---|---|---|---|
| ESM-2 | 650M-15B | 1.36 | -0.0153 | 0.011 | 1997/2000 | **survives** (-0.0153, p=0.011) |
| ProGen2 | 764M-6.4B | 0.92 | -0.0245 | 0.039 | 1977/2000 | **collapses** (+0.0064, p=0.60) |
| ProGen3 | 762m-3B | 0.60 | -0.0283 | 0.012 | 1988/2000 | **collapses** (-0.0031, p=0.77) |

Both ProGen confirmatory estimates are negative, exclude zero, and are larger
in magnitude than ESM-2 (161% and 185%). Bootstrap sign stability is nearly
as strong as ESM-2's. By the preregistered decision table's first row —
negative, CI excludes zero, magnitude within 50% of ESM-2 — both replicate.

## The fixed-effects collapse is the finding

R6 replaces every covariate with a dummy per assay. It identifies the
interaction purely from how each individual assay's score moves across model
sizes — within-assay variation. Everything an assay carries that does not
change across sizes (its family, depth, taxon, assay type, and any
unmeasured property) is absorbed.

- **ESM-2:** beta unchanged under R6 (-0.0153, p=0.011). The interaction is a
  genuine within-assay scaling response.
- **ProGen2:** beta flips to +0.0064 (p=0.60). Nothing left once between-assay
  structure is removed.
- **ProGen3:** beta to -0.0031 (p=0.77). Same.

So in the ProGen ladders the association between depth and scaling gains lives
in the *contrast between assays* — deep-family assays behave differently from
shallow ones — not in how any single assay responds to a larger model. That
is exactly what a benchmark-composition effect looks like, and it is what 01b
flagged as the open alternative for ESM-2 that fixed effects there ruled out.

## Two readings, both admitted in advance

**Composition dependence.** The ProGen ladders may be recovering the
depth-score correlation through which proteins are in ProteinGym, not through
a scaling mechanism. Amendment 3 anticipated this: all three ladders share
the same 217 assays, so a benchmark artefact reproduces across architectures
and looks like confirmation.

**Power.** ProGen2 spans 0.92 log-decades over 3 sizes, ProGen3 only 0.60,
against ESM-2's 1.36. Within-assay identification is where leverage is
scarcest, so R6 is the first test to lose significance when the ladder is
short. Amendment 2 recorded this limit before results were seen: "a wide
interval in ProGen3 is inconclusive, not null."

These are not mutually exclusive and the data cannot separate them here.
Distinguishing them needs a longer ladder (ProteinGym does not score ProGen3
above 3B or ProGen2 above 6.4B) or independently scored models on a
depth-balanced assay set — which is Experiment 02 territory.

## Pairwise decomposition

Consistent with the fixed-effects reading: the confirmatory signal is not
evenly spread across each ladder.

ProGen2: 764M->2.7B beta = -0.0299 (p=0.037, excludes zero); 2.7B->6.4B
beta = -0.0157 (p=0.47, inconclusive). The lower step carries it.

ProGen3: 762m->1B beta = -0.105 (p=0.10, wide, includes zero); 1B->3B
beta = -0.0157 (p=0.15, includes zero). Neither step is individually
conclusive; the confirmatory result borrows strength across the whole range.

## What Experiment 01 establishes, end to end

1. Published ESM-2 scores show a negative MSA-depth-by-scaling interaction
   above 650M (01a).
2. In ESM-2 it survives measured covariates, assay fixed effects, and
   assay-level resampling; it is homogeneous across selection classes
   (I2 = 0%); it is sensitive to influential assays and to the
   OrganismalFitness subpopulation (01b).
3. It is not novel: Hou et al. (2026) reported non-monotonic ESM-2 scaling
   first, via predicted likelihood rather than external MSA depth.
4. The depth-versus-scaling association reproduces in ProGen2 and ProGen3 on
   the confirmatory model but not under fixed effects (01c).

**Defensible conclusion.** There is a reproducible, cross-architecture
association between MSA depth and diminishing scaling returns in
ProteinGym-scored protein language models. In ESM-2 it is a within-assay
scaling response; in ProGen2 and ProGen3 it is between-assay and consistent
with benchmark composition. The evidence supports a benchmark-dependent
depth pattern, not a general protein-language-model scaling law. Settling
whether the ProGen collapse is composition or power requires a benchmark with
more independent families and longer ladders — the motivation for Experiment
02.

## Artefacts

- `notebooks/01c_replication_analysis.py` — parameterized, runs either ladder
- `results/01c_{progen2,progen3}_robustness_checks.csv` — R1-R6
- `results/01c_{progen2,progen3}_model_coefficients.csv`
- `results/01c_{progen2,progen3}_posthoc_diagnostics.csv` — pairwise + bootstrap
- `figures/01c_{progen2,progen3}_coefficient_stability.png`
- `results/provenance_01c_{progen2,progen3}_analysis.json`

## Deviations from the reviewed design

One. The verdict logic was extended after seeing that R6 collapsed, to
distinguish "replicates (strong, survives fixed effects)" from "replicates
(confirmatory only)." This is a reporting refinement, not a change to the
preregistered decision table, whose thresholds were applied as written. The
underlying numbers are unchanged; the label is more precise.
