# STATUS — where things stand

Last updated: end of Experiment 03.

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

- **Experiment 03** (complete): cross-family robustness of the depth–scaling
  association. Preregistered prediction FAILED in the informative direction.
  Balancing cross-family depth composition did NOT preferentially attenuate
  ProGen relative to ESM-2 — ESM-2 was balance-invariant (-0.0130 → -0.0124),
  both ProGen ladders grew more negative under balancing, and the differential
  D = A_ProGen − A_ESM2 leaned negative across both ProGen ladders and both
  chance-floor settings (P(D>0) 0.33–0.40). Arms 0-2 weakened three candidate
  artefacts: homologous-assay multiplicity (random-representative bootstrap
  strengthens every coefficient, so the ~39 redundant assays are not carrying
  the effect), the deterministic representative rule (conservative, not
  inflationary), and weak-signal assays (chance floor barely moves the estimate;
  the "holds only with the floor" row is not triggered). Direction is stable
  across every robustness cut; precision is what varies.
  See `experiments/03-cross-family-independence/results.md`.

## What 03 settles, and what it leaves open

Settled: the depth–scaling association is not an artefact of assay multiplicity,
representative choice, weak-signal assays, or simple Low/Medium/High benchmark
composition. The composition-sensitivity hypothesis for ProGen — the one I would
have bet on — is contradicted.

Still open, and now the whole remaining question:
1. **Scoring-protocol fairness.** Every 01/03 number is a *published* ProteinGym
   score. If context truncation/windowing disadvantages the 15B model on long
   (often deep-family) proteins, it could manufacture the ESM-2 within-assay
   result that the entire arc rests on. Untested.
2. **Composition vs power for ProGen.** ProteinGym cannot lengthen the ProGen
   ladders, so 03 could not separate a genuine cross-protein effect from
   short-ladder low power. ProGen stays negative but never clears 95% on any cut.

## Next: Experiment 04

Target fairness before power. Rationale: threat 1 sits under the ESM-2
within-assay result everything depends on, so it is higher-leverage than
threat 2, which only limits what ProGen can tell us. Threat 1 is also partly
checkable cheaply — from metadata already in the reference file (sequence length
vs each model's truncation threshold, correlated with MSA depth) — before
committing to any model re-scoring. Sequence: metadata diagnostic first; if it
finds the interaction concentrated in truncated / long-sequence assays, then
controlled re-scoring under a fixed protocol becomes mandatory (that is the
heavy, GPU/weights step 01's protocol always sequenced last). PROTOCOL FIRST —
committed before any diagnostic or scoring, as with 01/02/03. Power (threat 2)
is the natural Experiment 05: longer ladders or independently scored models,
not more benchmark reweighting.

## Workflow notes for the next session

- Work in the repo directly (Colab clone or local git), not via file deltas —
  file transport was the main friction; Colab sessions also recycle, wiping the
  clone, the git identity, and the PAT, so expect to re-clone, re-set
  `git config user.name/email`, and re-enter the token when resuming.
- Push from Python to avoid shell-interpolation of the token:
  `subprocess.run(["git","push", f"https://{token}@github.com/{REPO}.git","HEAD:main"])`,
  and scrub the token from any printed error.
- Verify every push by the DELTA, not just agreement: the new `git rev-parse HEAD`
  must differ from the pre-push hash AND match `git ls-remote origin main`.
  Matching hashes alone can mean "nothing pushed."
- Analysis inputs pinned to ProteinGym commit 144fe22b07dfaeec2b366f2346203a9838a55b4c.
  For 04, pin the exact truncation/windowing rule from the ProteinGym scoring
  code at that commit and record it in provenance — do not assume a context length.

## How to resume in a new chat

Open a fresh conversation and say: "Continuing work on
github.com/allanochola/protein-foundation-models — read STATUS.md, the
results.md files, and the Experiment 04 protocol if it exists. We're designing
Experiment 04." The repo carries all context; the chat does not.
