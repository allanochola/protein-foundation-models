# STATUS — where things stand

Last updated: end of Experiment 04 (Phase 1).

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

- **Experiment 04 — Phase 1** (complete): scoring-protocol fairness, windowing
  mechanism. Preregistered; PROTOCOL.md frozen before analysis. Reused the
  Experiment 01 estimator verbatim (`rho ~ lp + lp:ld + C(assay)`, cluster SE by
  UniProt_ID, upper segment) and changed only the assay set. Only 16/217 assays
  exceed ProteinGym's 1024-token window (all Low/Medium depth, 0 High); the
  effect lives on the 92.6% that fit fully in context. On the 201 never-windowed
  assays beta = -0.0180 vs the full-set -0.0153, retention 1.18, same sign —
  PASS. The 16 windowed assays carry none of the interaction (windowed-only
  beta ~ 0, p = 0.84; confound term lp:ld:windowed = +0.021, p = 0.20). Threat 1
  retired on published scores; Stage 3 re-scoring not triggered. Scope: clears
  the windowing mechanism, not every size-dependent scoring effect.
  See `experiments/04-scoring-window-fairness/results.md`.

## What is settled, and what it leaves open

Settled through 03: the depth–scaling association is not an artefact of assay
multiplicity, representative choice, weak-signal assays, or simple Low/Medium/High
benchmark composition. The composition-sensitivity hypothesis for ProGen — the one
I would have bet on — is contradicted.

Settled by 04: threat 1 — scoring-protocol fairness — is retired. The ESM-2
within-assay result holds on the 201 assays the windowing rule never touches, so
context windowing cannot have manufactured the number the arc rests on. (Caveat:
this clears the windowing mechanism specifically, not every conceivable
size-dependent scoring effect.)

Still open, and now the whole remaining question:
1. **Composition vs power for ProGen.** ProteinGym cannot lengthen the ProGen
   ladders, so 03 could not separate a genuine cross-protein effect from
   short-ladder low power. ProGen stays negative but never clears 95% on any cut.

## Next: Experiment 05

Threat 1 is closed, so power (threat 2) is the whole remaining question. ProGen's
depth–scaling signal is between-assay and never clears 95% on any 03 cut, and
ProteinGym cannot lengthen the ProGen ladders — so the open question is whether a
genuine cross-protein effect is present but under-powered, or absent. Attack it
with longer ladders or independently scored models (fresh inference), not more
benchmark reweighting. PROTOCOL FIRST — committed before any scoring, as with
01/02/03/04. The heavy GPU/weights step (running ProGen checkpoints for fresh
scores) is sequenced last, after any cheap metadata / power diagnostic.

## Workflow notes for the next session

- Work in the repo directly (Colab clone or local git), not via file deltas —
  file transport was the main friction; Colab sessions also recycle, wiping the
  clone, the git identity, and the PAT, so expect to re-clone, re-set
  `git config user.name/email`, and re-enter the token when resuming. Push with a
  `getpass`-typed token so it stays out of cell output and the notebook file.
- Verify every push by the DELTA, not just agreement: the new `git rev-parse HEAD`
  must differ from the pre-push hash AND match `git ls-remote origin main`.
  Matching hashes alone can mean "nothing pushed."
- Analysis inputs pinned to ProteinGym commit 144fe22b07dfaeec2b366f2346203a9838a55b4c.
  04's windowing rule is pinned: ESM-2 `model_window = 1024` (~1022 residues),
  size-invariant across the ladder; full sequence scored when len <= 1022. Still
  unpinned (non-blocking — only needed if the 16 windowed assays are ever
  re-scored): which long-sequence strategy produced the published ESM-2 scores
  (`wt-marginals`+`overlapping` vs `masked-marginals`+`optimal`).

## How to resume in a new chat

Open a fresh conversation and say: "Continuing work on
github.com/allanochola/protein-foundation-models — read STATUS.md and the
results.md files (through Experiment 04). We're designing Experiment 05." The repo
carries all context; the chat does not.
