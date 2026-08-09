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

ProtGPT3 was never benchmarked on ProteinGym (absent from the v1.3 baseline set;
the paper is generation-only), so there is no authors' rule to reproduce. We fix
one rule, identical across the three scales.

**Primary score — autoregressive wild-type marginal.** For a substitution at
residue position `p` (1-indexed) changing amino acid `a` → `b`, score it by the
model's next-token distribution at that position conditioned on the wild-type
prefix, read from a single forward pass over the wild-type sequence:

    x   = [<|bos|>, r_1, ..., r_L, <|eos|>]     # residue r_i at input index i
    lo  = model(x).logits                        # lo[t] predicts the token at index t+1
    z   = log_softmax(lo[p-1])                   # distribution over residue p | WT prefix r_1..r_{p-1}
    s_p = z[b] - z[a]                            # = lo[p-1][b] - lo[p-1][a]  (normalizer cancels)
    score(variant) = sum over substituted positions p of s_p     # additive; ONE WT forward per assay

Per-assay `rho` = Spearman(score, DMS_score) over the assay's variants.

This is a **left-context-only** estimand: each mutated residue is scored given the
N-terminal prefix, never the C-terminal suffix. That is the precise sense in which
an AR wild-type marginal differs from an ESM masked marginal, which conditions on
both flanks. State the consequence plainly: the ProtGPT3 primary score and the
ESM-2 reference (`-0.0180`, masked-marginal) are **different estimands**, so the
cross-family comparison tests the *presence, sign, and rough magnitude* of the
within-assay interaction, not a coefficient match. The within-ProtGPT3 `lp:ld`
interaction — one estimand applied identically across the ladder — is the primary
result; ESM-2 is a directional reference only.

WT-marginal is primary, not full-LLR, for two reasons: it is the only convention
that makes 217 assays × 3 checkpoints tractable (one forward per assay vs one per
variant — minutes vs A100-days on the 10B, see the runtime gate), and for the
mutated-position *difference* the softmax normalizer cancels, so the score is
invariant to the 20-vs-full-vocabulary choice that full-LLR is sensitive to.

**Mutation universe.** Restrict to the 20 standard amino acids
(`ACDEFGHIKLMNPQRSTVWY`). A variant is **dropped and recorded "not scored"**, never
silently mis-tokenized, if any substituted `a` or `b` is non-standard
(`X/B/Z/U/O/*` or gap), the parsed position lies outside `1..L`, or the stated
wild-type `a` disagrees with the reference residue at `p`. Assays reduced below a
usable variant count are recorded and excluded. Per-assay drop counts are reported.

**`<|bos|>` / `<|eos|>` handling.** Prepend `<|bos|>`, append `<|eos|>` — the
tokenizer omits `<|bos|>` even with `add_special_tokens=True`, but the model is
BOS-led, and position-0 conditionals are off-distribution without it.

**Full-sequence LLR is the locked validation convention — not primary, not
confirmatory.** score = logP(x_mut), one forward per variant, N→C, `<|bos|>`-led /
`<|eos|>`-terminated: the exact ProGen2 ProteinGym convention (mean-reduced over
residue tokens; within an assay this ranks variants identically to the sum or to
the WT-subtracted ratio). It is computed only on the predetermined validation
subset (runtime gate), only to test whether WT-marginal faithfully proxies it,
against a threshold fixed before the panel is run.

**Pre-protocol feasibility runs.** The three assays already scored with full-LLR at
112M (`A0A192B1T2_9HIV1_Haddox_2018`, `A0A1I9GEU1_NEIME_Kennouche_2019`,
`A0A247D711_LISMN_Stadelmann_2021`) were pre-freeze feasibility scouting. They are
recorded as such and are **not** confirmatory evidence; the confirmatory run
recomputes every assay fresh under the frozen convention.

## Assay set (locked)

The Experiment 04 **context-safe subset**: 201 DMS substitution assays with
reference sequences ≤ 1022 residues, so `<|bos|>` + sequence + `<|eos|>` fits the
1025 context with no windowing. This reuses the exact subset where windowing was
retired in 04, and it comes with a matched, same-assay ESM-2 reference already on
record: **`lp:ld` = −0.0180, p = 0.009** on these 201 assays. The full-217 fit is a
secondary, flagged with the windowing caveat; it is not the primary comparison.

## Runtime-feasibility gate (pre-inference, outcome-blind)

A frozen gate, run before any confirmatory scoring. It may inspect wall time,
throughput, peak memory, failures, and WT-marginal-vs-full-LLR scoring agreement.
It may **not** compute or inspect any cross-checkpoint `lp:ld` coefficient or MDE.
Runtime and scoring agreement decide whether to proceed — never the outcome.

**Panel (frozen here by geometry, before benchmarking, never by rho).** Five
context-safe assays selected deterministically from the reference by sequence
length and variant count (the selection rule and the resulting IDs are both frozen;
no post-hoc substitution):

| role | DMS_id | len | variants |
|---|---|---|---|
| short / low-variant | `OTU7A_HUMAN_Tsuboyama_2023_2L2D` | 42 | 635 |
| medium | `ESTA_BACSU_Nutschel_2020` | 212 | 2,172 |
| long / low-variant | `RAF1_HUMAN_Zinkus-Boltz_2019` | 648 | 297 |
| high-variant | `SPG1_STRSG_Olson_2014` | 448 | 536,962 |
| long + high-variant stress | `A0A192B1T2_9HIV1_Haddox_2018` | 852 | 12,577 |

Rule: short/low-variant = cheapest assay with length and variants both ≤ Q25;
medium = nearest to the median (length, variants) in rank space; long/low-variant =
fewest variants among length ≥ Q75; high-variant = most variants; stress = the fixed
852-aa Haddox case.

**Full-LLR is subsampled** to a fixed 2,000 variants per panel assay (seed 0) for
the agreement check — otherwise `SPG1`'s 537k variants alone would blow the budget.
WT-marginal scores every variant (one forward), so the WT-marginal side is never
subsampled. `SPG1` is a double-mutant landscape, so its agreement additionally
probes WT-marginal's additive multi-substitution handling; per-assay agreement is
reported so this is visible and not buried in the pooled statistic.

**Metrics recorded**, per checkpoint × convention: wall time, variants/sec, peak
GPU memory, OOM/failure count, projected GPU-hours to score all 217 assays.

**Numeric budget (fixed now, no post-hoc discretion).**
- 10B peak memory for a single WT-marginal forward (length ≤ 1024) ≤ **40 GB**,
  decided by the *measured* peak — not by the ~20 GB raw-weight estimate, which
  ignores kernels, router/expert buffers, logits, attention state, and framework
  overhead. If the measured peak exceeds 40 GB, the memory gate fails.
- WT-marginal production projected ≤ **3 A100-hours**. Projection uses 217 assays as
  a deliberately conservative bound though the confirmatory set is the 201
  context-safe assays.
- Full-LLR validation subset projected ≤ **8 A100-hours** (112M and 1.3B over the
  subsampled panel; 10B over only panel assays with ≤ 5,000 variants).

**Agreement criterion for WT-marginal (fixed before the panel).** On the validation
subset, WT-marginal is retained as primary only if BOTH hold:
- variant-level: median per-assay Spearman(WT-marginal scores, full-LLR scores)
  ≥ **0.90**; and
- assay-level: per-assay `rho`-vs-DMS under the two conventions agrees in sign on
  **every** subset assay, and Spearman(rho_marg, rho_llr) across the subset ≥ **0.80**.

**Failure behavior (transparent, no silent narrowing).**
- Any budget exceeded, or the 10B OOMs at 40 GB → record the feasibility failure and
  **stop**. Stage 3 does not silently drop to a 112M→1.3B ladder — that discards the
  long-span confound break that is its entire point. Any redesign (80 GB / offload /
  a different feasible long-span AR family) is a new, transparent decision.
- Agreement criterion failed → WT-marginal is not a faithful proxy and full-LLR is
  infeasible on the 10B, so Stage 3 as designed cannot proceed. Record and redesign;
  do not quietly switch primary conventions.

Only after the gate PASSES are the context-safe per-assay `rho` computed and the
cross-checkpoint FE `lp:ld` / MDE exposed.

## Estimator (unchanged from 01c/05)

`rho ~ lp + lp:ld + C(assay)`, cluster-robust SE by `UniProt_ID`, `lp` and `ld`
centered on the fitted ladder. `rho` is the per-assay Spearman between the primary
(WT-marginal) ProtGPT3 scores and DMS measurements. The Stage 2 fast within-estimator
and its statsmodels gate carry over unchanged.

## Analysis

**A — the interaction.** Single FE fit on the 3-point ProtGPT3 ladder, on the
primary (WT-marginal) `rho`. Report `lp:ld` with cluster-robust p and an
assay-cluster bootstrap CI, plus the effect-size ratios `beta/beta_ESM` against the
04 context-safe reference (−0.0180) and the full-217 ESM reference (−0.0153), each
with a bootstrap CI. Refit on the active-parameter axis as the locked robustness
(a sign change → inconclusive). Corroborate with the three 2-point sub-ladders
(112M→1.3B, 1.3B→10B, 112M→10B), bootstrap intervals only, as Stage 2 Analysis C
handicapped ESM-2 — instability expected and shown, not asserted away. The
WT-marginal-vs-full-LLR question is settled earlier, at the scoring level, by the
runtime gate; it is not re-litigated as a coefficient comparison here.

**B — MDE before any null is read.** Reuse the Stage 2 injection machinery verbatim
(no-interaction base + two-source noise bracket + dual-SE calibration) on ProtGPT3's
actual 3-point geometry and fresh-score noise. The primary bounded-`H_arch`
threshold is **MDE ≤ 0.0180**, the context-safe ESM-2 reference on the same
201-assay set that Stage 3 uses; **0.0153** (full-217 ESM-2) is a stricter secondary
benchmark. Fixing the primary reference to 0.0180 now removes any post-hoc choice of
whichever reference flips the verdict. Report the 80% MDE and power at both. A
ProtGPT3 null is only informative if the MDE clears 0.0180; otherwise the outcome is
"directional, unresolved," the Stage 2 discipline.

## Primary prediction (stated to fail) and decision

Pre-registered expectation: **H_power** — a negative `lp:ld` appears in ProtGPT3,
because the AR ladder now exceeds the span that carries the ESM-2 effect. Stating
the bet first, as in Stage 2, so a null is a real falsification and not a quiet
reinterpretation. Note honestly that Stage 2's evidence leaned the other way (A and
C toward architecture); this bet is the falsifiable one, not the likely one.

| ProtGPT3 result | Verdict |
|---|---|
| `lp:ld` significantly negative, magnitude broadly compatible with ESM-2 (ratio CI overlapping ~1), stable across the total- and active-parameter axes | **H_power supported** — an AR family with sufficient span shows the within-assay response; the ProGen null plausibly reflected identification/power. |
| Null, and MDE ≤ 0.0180 (primary context-safe ESM-2 reference) at ≥ 0.8 power | **Bounded H_arch (narrow)** — the ESM-2 behavior does not reproduce in a long-span AR (MoE) family; a smaller effect is not excluded. State the MDE bound; 0.0153 is the stricter secondary check. |
| Null, MDE > 0.0180 | **Directional, unresolved** — more AR span did not surface it, but the design still cannot exclude an ESM-sized effect. |
| Significantly positive | **Architecture-specific divergence** — report as its own finding, not folded into the binary. |
| Sign flips between the total- and active-parameter axes | **Inconclusive** — parameter-axis robustness disagreement, reported as such. (WT-marginal-vs-full-LLR agreement is settled in the gate, not here.) |

## Kill criteria (technical scoring sanity, before the FE fit)

These identify a **broken scorer**, never an unfavorable biological result. None
inspects cross-checkpoint ordering or the sign/size of any interaction.

- Scores finite and non-constant within each assay (a constant or all-NaN column is
  a scorer failure, not a result).
- Correct wild-type token lookup: the reference residue at each scored position
  matches the tokenized wild-type; mismatches are counted and, above a small
  tolerance, halt.
- Synthetic unit test passes: on a constructed sequence with a known injected
  next-token bias, the scorer recovers the expected mutation sign and the
  batched score equals the per-sequence score (the offline-validated identities).
- Sufficient usable variants per assay after the mutation-universe filter; assays
  below the floor are recorded "not scored," not silently fit.
- No systematic parsing/tokenization failure: the fraction of dropped/unparsable
  variants is within a fixed tolerance, else halt and inspect.
- Any assay exceeding context after bos/eos is recorded "not scored," never
  truncated silently.

Descriptive diagnostics (mean per-assay `rho`, its spread, per-checkpoint means)
are recorded for transparency but are **not** kill criteria — a low or size-ordered
`rho` may be the real Stage 3 result and must not gate the fit.

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

- `notebooks/05_stage3_protgpt3_scoring.py` — WT-marginal (primary) + full-LLR scoring
- `notebooks/05_stage3_runtime_gate.py` — the outcome-blind feasibility + agreement gate
- `notebooks/05_stage3_analysis.py` — Stage 2 estimator + injection MDE, run only post-gate
- `results/05_stage3_runtime_gate.csv` — per checkpoint × convention: wall time, variants/sec, peak GPU memory, failures, projected GPU-hours
- `results/05_stage3_agreement.csv` — validation-subset variant- and assay-level WT-marginal-vs-full-LLR agreement against the fixed thresholds
- `results/05_stage3_protgpt3_scores.csv` — per-assay WT-marginal `rho` per checkpoint (all context-safe assays), plus full-LLR `rho` on the validation subset, with drop counts
- `results/05_stage3_interaction.csv` — A: `lp:ld`, ratios, bootstrap CIs, total- and active-parameter axes
- `results/05_stage3_injection_power.csv`, `results/05_stage3_calibration.csv` — B
- `results/provenance_05_stage3.json` — ProtGPT3 revisions, exact/active params, scoring convention, gate thresholds and outcome, seed, reps, context-safe assay list

## Sequence of work (frozen)

Freeze this protocol → write the gate notebook → run the gate → commit gate outputs
→ **only if the gate PASSES** on budget and agreement, run the full WT-marginal
scoring of 112M/1.3B/10B → compute all context-safe `rho` → expose the
cross-checkpoint FE `lp:ld` / MDE last. No cross-checkpoint coefficient is computed
before the gate clears.

## Compute

WT-marginal (primary) is one forward per assay: ~201 forwards × 3 checkpoints —
minutes. The 10B is 2.75B active, ~20 GB raw bf16 weights, but whether it fits a
single 40 GB A100 in practice is a *hypothesis the runtime gate tests* by measured
peak memory, not a guarantee — kernels, expert/router buffers, logits, and
attention state add to raw weights. Full-LLR (validation only) is one forward per
variant, bounded to the subsampled panel and to ≤ 5,000-variant assays on the 10B,
under the gate's ≤ 8 A100-hour ceiling. The full-benchmark full-LLR run the earlier
draft assumed (~A100-days, 10B-dominated) is deliberately avoided — that cost is
what demoted full-LLR from primary to a bounded validation convention.
