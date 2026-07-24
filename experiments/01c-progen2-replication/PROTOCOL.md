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

**Ladder construction rule, fixed before fitting.** Base and M differ by
training corpus, not parameter count — both are 764M. Including both would
put two points at the same x with different y and bias the slope. Drop
`Progen2 Base`, keep `Progen2 M`. Parameter counts to be verified against
Nijkamp et al. 2023 before the fit; if verification contradicts these
figures, the corrected values are used and this paragraph is amended with the
original left struck.

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

## What this cannot establish

ProGen2 and ESM-2 were scored by ProteinGym under its own protocol. A shared
artefact of that protocol — context truncation affecting large models,
batching differences — would produce correlated results in both ladders and
look like replication. Ruling that out requires scoring models directly,
which is a separate experiment.
