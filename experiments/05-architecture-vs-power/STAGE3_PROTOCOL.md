# Experiment 05 — Stage 3 protocol (ProtGPT3 confound break)

Frozen before any ProtGPT3 sequence is scored. Triggered by the Stage 2 result
(`results.md`): the preregistered H_power prediction failed on ProGen, but the
ProGen upper ladders were underpowered for an ESM-2-sized within-assay effect, so
the null could not bound absence. Architecture and ladder geometry stayed
confounded. Stage 3 breaks that confound with the object the frozen Stage 3 plan
called for and Stage 2 could not supply: a single-sequence autoregressive protein
LM whose ladder spans wider than ProGen's.

## The question

Does the ESM-2 within-assay depth×scaling interaction (`lp:ld < 0`) appear in an
autoregressive protein LM once the AR ladder carries substantially more
parameter-span leverage than ProGen? ESM-2's effect is established on a 1.36-decade
upper ladder. ProGen2 (1.63) and ProGen3 (1.43) showed nothing, but underpowered.
ProtGPT3 spans 1.96 decades (total) / 1.91 (active). If the interaction is real and
merely hidden in ProGen by short-ladder identification, it should surface here. If
it stays absent at this span, "ProGen was merely too short" becomes hard to
maintain.

## What Stage 3 can and cannot separate

ProtGPT3 is an autoregressive **Mixtral-style sparse MoE** family, not a dense
transformer. So it breaks the ladder-length confound but not the architecture
confound cleanly: masked-vs-causal and dense-vs-MoE remain partially entangled
across the arc (ESM-2 is masked-dense; ProGen3 and ProtGPT3 are causal-MoE). A null
here therefore supports the *narrow* claim — the ESM-2 interaction fails to
reproduce in a second autoregressive family, one with far greater span than ProGen
— not the broad claim that causal LMs lack it. Stated before results so the
inference is not widened afterward.

## Model (locked)

Three single-sequence base checkpoints, `AI4PD/ProtGPT3-{112M,1.3B,10B}`,
`MixtralForCausalLM`, 8 experts, top-2 routing, 31-token residue-level vocabulary,
1025-position context, identical across scales. Exact counts audited from
checkpoint metadata, no reliance on the name tags:

| checkpoint | total params | active params (top-2) | active frac |
|---|---|---|---|
| 112M | 0.109B | 0.034B | 31% |
| 1.3B | 1.328B | 0.366B | 28% |
| 10B  | 10.000B | 2.752B | 28% |

Span: total 1.96 decades, active 1.91. Excluded and why: `-MSA` (MSA-promptable —
reintroduces the depth axis into the model, which is the independent variable);
`-dpo` (DPO alignment changes the likelihood surface scoring depends on). The base,
unaligned, single-sequence checkpoints only.

`lp` primary axis = `log10(total params)`, for consistency with 01–05 and with the
ProGen3 MoE ladder already handled that way. **Active parameters are a locked
robustness axis**: if the sign of `lp:ld` or the verdict changes between the total-
and active-parameter fits, Stage 3 is inconclusive — not a post-hoc choice.

## Scoring convention (frozen across all three checkpoints)

ProtGPT3 was not benchmarked on ProteinGym (absent from the v1.3 baseline set; the
paper is generation-focused), so there is no authors' rule to reproduce and we
define one, fixed identically across scales:

- **Full-sequence log-likelihood ratio.** score(variant) = log P(x_mut) −
  log P(x_wt), summed over positions autoregressively (N→C, the only direction the
  family supports). Matches the convention ProGen was scored under, so the
  within-AR comparison is fair by the Experiment 04 standard.
- **Prepend `<|bos|>`, append `<|eos|>`.** The tokenizer does not add `<|bos|>`
  even with `add_special_tokens=True`, but the model was trained BOS-led; scoring
  without it puts position-0 likelihoods off-distribution. The bos/eos terms are
  included in the sum identically for WT and mutant.
- Multi-substitution variants: one forward per mutant sequence (no additive
  approximation), same as WT.
- **WT-marginal is a locked robustness pass**, not the primary: prefix-conditional
  log-odds at the mutated position from a single WT forward. Cheap (~201 forwards
  per checkpoint) and prespecified. If the sign of `lp:ld` disagrees between
  full-LLR and WT-marginal, Stage 3 is inconclusive.

## Assay set (locked)

The Experiment 04 **context-safe subset**: 201 DMS substitution assays with
reference sequences ≤ 1022 residues, so `<|bos|>` + sequence + `<|eos|>` fits the
1025 context with no windowing. This reuses the exact subset where windowing was
retired in 04, and it comes with a matched, same-assay ESM-2 reference already on
record: **`lp:ld` = −0.0180, p = 0.009** on these 201 assays. The full-217 fit is a
secondary, flagged with the windowing caveat; it is not the primary comparison.

## Estimator (unchanged from 01c/05)

`rho ~ lp + lp:ld + C(assay)`, cluster-robust SE by `UniProt_ID`, `lp` and `ld`
centered on the fitted ladder. `rho` is per-assay Spearman between ProtGPT3 scores
and DMS measurements. The Stage 2 fast within-estimator and its statsmodels gate
carry over unchanged.

## Analysis

**A — the interaction.** Single FE fit on the 3-point ProtGPT3 ladder. Report
`lp:ld` with cluster-robust p and an assay-cluster bootstrap CI, plus the effect-
size ratios `beta/beta_ESM` against the 04 context-safe reference (−0.0180) and
against the full-217 ESM reference (−0.0153), each with a bootstrap CI. Corroborate
with the three 2-point sub-ladders (112M→1.3B, 1.3B→10B, 112M→10B), bootstrap
intervals only, exactly as Stage 2 Analysis C handicapped ESM-2 — instability
expected and shown, not asserted away.

**B — MDE before any null is read.** Reuse the Stage 2 injection machinery verbatim
(no-interaction base + two-source noise bracket + dual-SE calibration) on ProtGPT3's
actual 3-point geometry and fresh-score noise. Report the 80% MDE and the power at
−0.0180 and −0.0153. A ProtGPT3 null is only informative if the MDE clears the
reference magnitude; otherwise the outcome is "directional, unresolved," same
discipline that governed Stage 2.

## Primary prediction (stated to fail) and decision

Pre-registered expectation: **H_power** — a negative `lp:ld` appears in ProtGPT3,
because the AR ladder now exceeds the span that carries the ESM-2 effect. Stating
the bet first, as in Stage 2, so a null is a real falsification and not a quiet
reinterpretation. Note honestly that Stage 2's evidence leaned the other way (A and
C toward architecture); this bet is the falsifiable one, not the likely one.

| ProtGPT3 result | Verdict |
|---|---|
| `lp:ld` significantly negative, magnitude broadly compatible with ESM-2 (ratio CI overlapping ~1), stable across total/active and LLR/marginal | **H_power supported** — an AR family with sufficient span shows the within-assay response; the ProGen null plausibly reflected identification/power. |
| Null, and MDE ≤ reference magnitude at ≥ 0.8 power | **Bounded H_arch (narrow)** — the ESM-2 behavior does not reproduce in a long-span AR (MoE) family; a smaller effect is not excluded. State the MDE bound. |
| Null, MDE > reference magnitude | **Directional, unresolved** — more AR span did not surface it, but the design still cannot exclude an ESM-sized effect. |
| Significantly positive | **Architecture-specific divergence** — report as its own finding, not folded into the binary. |
| Sign flips across total/active params or LLR/marginal | **Inconclusive** — robustness disagreement, reported as such. |

## Kill criteria (scoring sanity, before the FE fit)

- Scoring pipeline validated before any coefficient: mean per-assay ρ in a sane
  range (roughly 0.3–0.5 for a family of this scale) and non-degenerate per-assay ρ
  distribution; larger ProtGPT3 should not score *worse* on aggregate. If the
  pipeline fails these, stop and fix before fitting — a broken scorer manufactures
  or hides interactions.
- Any assay whose sequence still exceeds context after bos/eos drops and is
  recorded "not scored," never truncated silently.
- Calibration of Analysis B (DGP-validity rate) outside [0.015, 0.040] → B
  discarded, verdict rests on A and the MDE-free reading, as in Stage 2.

## What this cannot establish

- Masked-vs-causal and dense-vs-MoE stay entangled: a null is evidence about
  autoregressive-MoE families, not causal LMs in general. Dense-AR rests on ProGen.
- Three checkpoints. The interaction test and its MDE are the whole story; there is
  no segmented within-ladder decomposition to lean on.
- Zero-shot ProteinGym substitution prediction on ≤1022-residue assays only. Not a
  scaling law; a benchmark-specific, depth-dependent property, now tested for
  reproduction in one more family.
- Fresh scores are ours, not ProteinGym's; the ESM-2 comparison is a magnitude
  reference across differing (necessarily architecture-specific) scoring functions,
  not a pooled fit.

## Outputs

- `notebooks/05_stage3_protgpt3.py` — scoring + the Stage 2 estimator/analyses
- `results/05_stage3_protgpt3_scores.csv` — per-assay ρ per checkpoint (both conventions)
- `results/05_stage3_interaction.csv` — A: `lp:ld`, ratios, bootstrap CIs, total/active
- `results/05_stage3_injection_power.csv`, `results/05_stage3_calibration.csv` — B
- `results/provenance_05_stage3.json` — ProtGPT3 revisions, exact/active params, scoring convention, seed, reps, context-safe assay list

## Compute

Full-LLR is one forward per variant: ~2.5M variants × 3 checkpoints on the
context-safe subset. The 10B is MoE (2.75B active), ~20 GB bf16 weights — fits a
single A100/H100 40 GB; vLLM/SGLang `prompt_logprobs` scores prompts efficiently.
Estimate ~1 A100-day for the 10B, far less for the two smaller checkpoints. This is
the GPU arm the frozen protocol reserved as "conditional"; it is now warranted and
tractable. WT-marginal robustness is ~201 forwards per checkpoint — negligible.
