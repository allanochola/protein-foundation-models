# STATUS — where things stand

Last updated: end of Experiment 02.

## Done

- **Experiment 01** (complete, on GitHub): depth × scaling interaction in
  ProteinGym zero-shot variant-effect prediction. Within-assay in ESM-2
  (survives assay fixed effects, beta = -0.0153, p = 0.011); between-assay in
  ProGen2 and ProGen3 (both replicate on the confirmatory model but collapse
  under fixed effects). Reconciled with Hou et al. 2026 — the novelty is the
  fixed-effects + cross-architecture angle, not the existence of non-monotonic
  scaling, which Hou reported first via predicted likelihood.
  See `experiments/01-msa-depth-confound/results.md` and
  `experiments/01c-progen2-replication/results.md`.

- **Experiment 02** (complete): ProteinGym family-independence audit.
  217 assays collapse to only 178 independent sequences at 50% identity
  (172 at 30%); the preregistered prediction of "<100" failed. ProteinGym is
  NOT a few dominant families — largest cluster is 1.8%. Coverage-sensitivity
  checked (floor ~160). All four decision rules pass.
  See `experiments/02-proteingym-family-audit/results.md`.

## What 02 changes for 01

Removes one of three explanations for the ProGen between-assay signal: it is
NOT an artefact of a few over-represented families. Does NOT resolve between
the remaining two — genuine composition effect across many independent
proteins, or low power from short ProGen ladders. Keep that restraint; do not
over-read benchmark diversity as evidence the depth effect is mechanistic.

## Next: Experiment 03

Cross-family study is warranted (all 02 rules pass). Sample one assay per
50%-identity cluster → ~178 units, well-powered, no hybrid case-study fallback
needed. PROTOCOL FIRST — committed before any sampling or analysis, as with
01 and 02. Alternative slot: structure-aware models under pLDDT stratification.

## Workflow notes for the next session

- Work in the repo directly (Colab clone or local git), not via file deltas —
  file transport was the main friction this cycle.
- Verify every push: `git rev-parse HEAD` must match `git ls-remote origin main`.
- Analysis inputs pinned to ProteinGym commit 144fe22b07dfaeec2b366f2346203a9838a55b4c.
- mmseqs2 static binary (for Tier B clustering) pulls from the MMseqs2 GitHub
  releases page — no bioconda needed.

## How to resume in a new chat

Open a fresh conversation and say: "Continuing work on
github.com/allanochola/protein-foundation-models — read STATUS.md, the
results.md files, and the Experiment 03 protocol if it exists. We're designing
Experiment 03." The repo carries all context; the chat does not.
