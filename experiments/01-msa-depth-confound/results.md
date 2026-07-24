# 01b results

**Verdict: partial survival.**

Experiment 01b partially supports the interaction observed in 01a. The
coefficient is negative in all 14 specifications, is not attenuated by
measured covariates or by assay fixed effects, and an assay-level cluster
bootstrap places 1,997 of 2,000 replicates below zero with a 95% percentile
interval of [−0.0300, −0.0060]. Pairwise decomposition shows the pattern is
not created solely by the 15B endpoint.

However, significance is lost under HC3 uncertainty, under removal of
OrganismalFitness assays, under removal of Human assays, and under exclusion
of the 6.6% of observations flagged by Cook's distance. High-influence
observations are disproportionately low-depth assays, including cases where
all three models score near chance.

**The measured composition-artefact explanation is not supported.** That is
narrower than saying composition is ruled out. Assay fixed effects protect
against stable assay-level properties; they do not protect against an
interaction between model size and an unmeasured assay property. Still open:
endpoint-specific measurement noise, model-specific failure modes correlated
with depth, benchmark construction effects, scoring-protocol effects such as
context truncation, and nonlinear size-response relationships.

The leave-category-out sensitivity limits the scope of the claim and may
reflect genuine effect heterogeneity, benchmark composition dependence, or
both. The class-specific analysis below argues against heterogeneity and for
composition, but does not settle it.

**Defensible conclusion.** The published ESM-2 ladder contains a reproducible
negative association between MSA depth and scaling gains above 650M. The sign
is stable under assay-level resampling and survives measured covariate
adjustment and assay fixed effects. It is not solely a 15B endpoint artefact.
Its magnitude and statistical certainty depend on influential assays and on
the largest benchmark subpopulation. The evidence supports a
benchmark-specific depth-dependent pattern in published ESM-2 scores — not a
general protein-language-model scaling law. Establishing that requires an
independent architecture (01c) and, ultimately, independently scored
models.

## Applying the repository's attenuation rubric

Confirmatory model: β01b = −0.0172 against β01a = −0.0153. Attenuation
**−12.7%** — the estimate grows slightly rather than shrinking. Same sign,
under 25% attenuation, CI [−0.0304, −0.0040] excluding zero. On the
confirmatory model alone that reads as **strong survival**.

The robustness battery contradicts that reading. Four of fourteen fits have
CIs including zero, and the two that also attenuate are the informative ones:
dropping OrganismalFitness (+30% attenuation) and dropping high-influence
observations (+13%). The honest classification across the whole battery is
**partial survival**.

Reporting only the confirmatory model would have earned "strong survival" on
a rubric designed to prevent exactly that.

## Where the effect lives

The raw picture, upper segment, depth terciles:

| MSA depth | 650M | 3B | 15B |
|---|---|---|---|
| shallow | 0.348 | 0.353 | 0.366 |
| mid | 0.450 | 0.457 | 0.447 |
| deep | 0.518 | 0.487 | 0.463 |

Deep families decline monotonically. Shallow ones rise slightly. The
interaction is that divergence, and it is visible without any model.

## Influence diagnostics

43 of 651 observations (6.6%) exceed Cook's D > 4/n, spanning 32 distinct
assays. Their composition against the overall benchmark:

| Stratum | over/under-represented |
|---|---|
| Virus | 2.1× |
| Binding | 2.3× |
| Activity | 1.9× |
| Low MSA depth | 1.8× |
| Human | 0.68× |
| OrganismalFitness | 0.52× |

**The 3B model contributes zero high-influence observations.** All 43 come
from 650M (53%) and 15B (47%). Influence is therefore concentrated at the
endpoints, so the estimated interaction is disproportionately — not
exclusively — determined by endpoint contrasts rather than evenly supported
across the ladder. The 3B rows still constrain the fitted slope; zero
individually flagged 3B observations does not mean the middle model carries
no information. The pairwise decomposition below shows it carries a great
deal.

**Low-depth assays are 1.8× overrepresented.** The top influential cases
include two GFP-homolog assays (Somermeyer 2022) where ESM-2's Spearman is
approximately −0.02 — models are performing at chance. The "15B wins on
shallow families" half of the finding partly rests on assays where every
model fails and small differences are noise.

That is the single most important caveat 01b produced, and it was not
visible until the influence diagnostics were run.

## Post hoc diagnostics

Motivated by the robustness results, not preregistered, and not confirmatory.

**Pairwise endpoint decomposition.**

| Contrast | beta | p | 95% CI |
|---|---|---|---|
| 650M → 3B | −0.0205 | 0.029 | [−0.0388, −0.0021] |
| 3B → 15B | −0.0142 | 0.071 | [−0.0296, +0.0012] |
| 650M → 15B | −0.0172 | 0.011 | [−0.0305, −0.0040] |

The interaction is not created solely by the 15B endpoint. The 650M-to-3B
contrast is negative and excludes zero (beta = −0.0205, 95% CI [−0.0388,
−0.0021]), while the 3B-to-15B contrast has the same estimated direction but
remains inconclusive (beta = −0.0142, 95% CI [−0.0296, +0.0012]). The middle
model therefore materially constrains the result, but the data do not
establish the interaction separately in both scaling steps.

**Assay-level cluster bootstrap.**

| Parameter | Value |
|---|---|
| Sampling unit | assay |
| Sampling | assays drawn with replacement, n = 217 per replicate |
| Rows retained | all three model-size observations per sampled assay |
| Refit model | primary covariate-adjusted specification |
| Replicates | 2,000 |
| Failed fits | 0 |
| Median beta | −0.0174 |
| 95% percentile interval | [−0.0300, −0.0060] |
| Replicates with beta < 0 | 1,997 of 2,000 |

Resampling complete assay clusters, not individual rows, is required by the
repeated-measures structure — each assay contributes three correlated
observations.

**Interaction by selection class.** Deleting a class tests only whether it
was load-bearing. Estimating within each class tests whether the effect
differs across them.

| Class | beta | 95% CI | assays |
|---|---|---|---|
| Activity | −0.0096 | [−0.0261, +0.0069] | 43 |
| Binding | −0.0161 | [−0.0625, +0.0304] | 13 |
| Expression | −0.0183 | [−0.0504, +0.0138] | 18 |
| OrganismalFitness | −0.0305 | [−0.0585, −0.0024] | 77 |
| Stability | −0.0051 | [−0.0320, +0.0219] | 66 |

Random-effects pooling (DerSimonian-Laird): beta = −0.0137, 95% CI
[−0.0250, −0.0023]. Q = 2.10 on 4 df, **I-squared = 0%**.

All five classes point the same direction with no detectable heterogeneity.
The OrganismalFitness sensitivity therefore reflects loss of precision — it
is the largest class at 77 assays — rather than a genuinely heterogeneous
effect or a misspecification averaging incompatible processes. With five
classes and wide individual intervals, the heterogeneity test is weakly
powered, so this narrows that question rather than settling it.

The bootstrap and the leave-out checks disagree, and the disagreement is
informative rather than contradictory. The bootstrap asks whether the
estimate is stable under sampling variability — it is, strongly. Leave-one-
category-out asks whether it depends on a particular subpopulation — it does.
Both hold simultaneously.

## Artefacts

- `results/01b_merged_assay_data.csv` — assay-level merged data with
  influence columns
- `results/01b_model_coefficients.csv` — confirmatory and secondary specs
- `results/01b_robustness_checks.csv` — R1–R6, 14 fits
- `results/01b_influence_diagnostics.csv` — Cook's D and leverage per
  observation, sorted
- `results/01b_posthoc_diagnostics.csv` — pairwise contrasts and bootstrap
- `results/01b_by_selection_class.csv` — class-specific interaction estimates
- `figures/01b_coefficient_stability.png` — forest plot, all 14 fits
- `figures/01b_adjusted_interaction.png` — terciles with influential points
  overlaid

## Deviations from the reviewed design

Three, all recorded rather than silently taken.

**Confirmatory specification.** The reviewed design put the coefficient of
interest on `upper_segment:log_msa_depth`. That tests whether the *level*
association between depth and score shifts above 650M, not whether the
*scaling slope* varies with depth. The quantity of interest is
`log_size:log_depth` estimated within the upper segment, which is what was
fitted.

**Assay fixed effects added.** Not in the reviewed design. They absorb taxon,
length, selection type and every unmeasured assay-level property in one step,
which makes the covariate-selection argument moot. Reported as R6.

**Directory.** 01b lives under `experiments/01-msa-depth-confound/` rather
than a separate `01b-msa-depth-covariates/`. It is the same question with
the same data; splitting it would separate the finding from the protocol
that predicted it.

## Process failure

01b was run in the same session it was designed. The kill criterion was
preregistered; specifications 3 and 4 were not. They moved the estimate
against the hypothesis, so nothing was gained — but that is unverifiable from
outside, which is the point of preregistering. 01c is preregistered and
committed before any 01c code exists.
