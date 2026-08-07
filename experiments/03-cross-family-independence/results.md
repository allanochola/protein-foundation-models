# 03 results — cross-family robustness of the depth–scaling association

**Verdict: the preregistered primary prediction failed, in the informative
direction. Equalizing cross-family depth composition did not preferentially
attenuate ProGen relative to ESM-2 — under balancing ESM-2 was nearly
unchanged and the ProGen coefficients grew more negative. Three candidate
artefacts (assay multiplicity, the representative rule, weak-signal assays) are
weakened; the simple Low/Medium/High composition explanation for ProGen is
weakened rather than confirmed. The composition-versus-power question for the
short ProGen ladders is not resolved and passes to Experiment 04.**

All numbers below are on published ProteinGym scores, pinned commit
`144fe22b`, seed 0, 2000 reps per bootstrap. Every specification was fixed in
`PROTOCOL.md` before any 03 number was seen.

## What the prediction was

The differential `D = A_ProGen − A_ESM2` was predicted positive with its
interval excluding zero: ProGen, carrying a between-assay signal, should lose
more of its interaction than within-assay ESM-2 when the cross-family depth
mix is balanced. It did not.

## Arms 0-2 — the association is robust to independence, representative choice, and the floor

Restricting ProteinGym to one representative per 50%-identity cluster (178
units) left the interaction essentially where 01 had it. The 217→178
reweighting moved nothing, as expected — 02 had already shown the benchmark was
near-independent, so there was little multiplicity to remove.

The informative arm is the representative-sensitivity check, not the
reweighting. Randomizing which member represents each multi-assay cluster made
every coefficient **more** negative, on both floor settings and all three
ladders:

| Ladder | deterministic (floor off) | random-rep (floor off) |
|---|---|---|
| ESM-2 | −0.0132 | −0.0142 |
| ProGen2 | −0.0171 | −0.0205 |
| ProGen3 | −0.0190 | −0.0208 |

If the ~39 redundant assays were carrying the association, randomizing their
representatives would have pulled the estimate toward zero. It moved the other
way. So multiplicity of sequence-related assays does not explain the effect —
this goes past 02, which established diversity; arm 2 shows the estimate does
not ride on which redundant system gets the weight.

The chance-floor comparison came out cleanly, and against the intuition the
floor was built on. Removing the 8 weak-signal ESM-2 assays did not reveal the
interaction — it slightly reduced precision:

| ESM-2 | β (boot median) | 95% interval | P(neg) |
|---|---|---|---|
| floor off (178) | −0.0132 | [−0.0270, −0.0005] | 0.978 |
| floor on (170) | −0.0138 | [−0.0281, +0.0009] | 0.969 |

The point estimate barely moves (0.0006); the interval widens and crosses zero
because 8 systems of power are gone. The decision-table row "holds only with
the floor, collapses without" is **not triggered** — the opposite holds, the
effect is marginally stronger without the floor. Weak-signal assays are not
carrying the ESM-2 result.

Across every arms-0-2 cut — independent families, floor on, floor off,
deterministic and random representatives — ESM-2 stayed in −0.013 to −0.017 and
the ProGen ladders in −0.017 to −0.024. **The coefficient is stable; the
precision is what varies.** The data support robustness of direction more
strongly than robustness of conventional significance. ESM-2 clears the 95%
line on the full independent set (p=0.045) but not after the floor (p=0.076);
ProGen never clears it on either — consistent with, though not proof of, the
longer-ladder / within-assay reading for ESM-2 versus short-ladder leverage for
ProGen.

## Arms 3-4 — depth balancing contradicts the prediction

Balancing the cross-family depth composition (equal Low/Medium/High weight)
left ESM-2 almost exactly where it was and pushed both ProGen ladders further
from zero. β_original is the native-depth confirmatory estimate on the common
set; β_balanced is the depth-balanced median.

| Floor | Ladder | β_original | β_balanced |
|---|---|---|---|
| on (167, n_balance 25) | ESM-2 | −0.0118 | −0.0116 |
| | ProGen2 | −0.0154 | −0.0247 |
| | ProGen3 | −0.0176 | −0.0264 |
| off (178, n_balance 29) | ESM-2 | −0.0130 | −0.0124 |
| | ProGen2 | −0.0171 | −0.0198 |
| | ProGen3 | −0.0184 | −0.0238 |

ESM-2 is balance-invariant; ProGen amplifies under balancing. That is a
directional fact about the raw coefficient, not about the ratio statistic, so
the heavy-tailed A and D intervals below do not threaten it. If ProGen's signal
were a composition artefact of ProteinGym's native depth mix, equal depth
weight should have pulled it toward zero. It went the other way.

The differential confirms the prediction failed, read by sign, not magnitude:

| Contrast | Floor | D median | P(D>0) |
|---|---|---|---|
| ProGen2 − ESM-2 | on | −0.560 | 0.329 |
| ProGen3 − ESM-2 | on | −0.438 | 0.338 |
| ProGen2 − ESM-2 | off | −0.249 | 0.401 |
| ProGen3 − ESM-2 | off | −0.323 | 0.367 |

D was predicted positive; all four are negative at the median with P(D>0)
between 0.33 and 0.40, consistent across both ProGen ladders and both floor
settings. The bootstrap intervals are enormous — A_esm2 spans [−1.90, +0.95] —
because A divides by β_original ≈ −0.012 on 75-family draws; **D's magnitude is
not interpretable and is not interpreted here.** P(D>0) is a sign count,
immune to the tail, and it is a clean, consistent signal that the differential
does not lean the predicted way. No positive differential attenuation; equally,
the intervals forbid claiming a reliably negative one.

Note on the `material_attenuation` column in `03_depth_balanced.csv`: it is an
implementation flag that fires whenever the balanced interval includes zero,
regardless of A's sign. All six rows trip it, while meaning opposite things
(ESM-2 unchanged, ProGen amplified). It is not a substantive result; the
components — β_original, β_balanced, A, P(D>0) — are reported separately above
and should be read instead.

## Against the decision table

- *ProGen attenuates strongly; ESM-2 stable; D excludes zero* — **not met.**
  ProGen did not attenuate; D does not exclude zero and does not lean positive.
- *All three stable / composition explanation weakened* — **this is the cell.**
  No interaction collapsed toward zero under balancing; ESM-2 held, ProGen
  strengthened. The simple Low/Medium/High composition account for ProGen loses
  ground.
- *Holds only with floor, collapses without* — **not triggered** (arms 0-2).
- *ProGen negative but bootstrap highly uncertain → power unresolved* — **still
  live.** ProGen stays negative and never clears 95% on any cut; ProteinGym
  cannot lengthen those ladders, so composition-versus-power is not settled by
  03.

## What this establishes and what it does not

Weakened: homologous-assay multiplicity, deterministic representative
selection, weak-signal assays, and simple benchmark depth composition as
explanations for the depth–scaling association. The predicted
composition-sensitivity of ProGen relative to ESM-2 is contradicted — gently,
given the power, but consistently across four differentials.

Not resolved: whether ProGen's between-assay signal is a genuine cross-protein
effect or short-ladder low power (03 cannot separate these on ProteinGym), and
whether any of it is a scoring-protocol artefact of context truncation on
published scores. Both pass to **Experiment 04**, which is now the natural next
step: score models under a controlled protocol, or use a longer ladder, rather
than reweight the benchmark further.

This is a genuine falsification. The composition explanation for ProGen was
preregistered and failed in the informative direction — the prediction I would
have bet on did not survive contact with the balanced design, which is exactly
why it was frozen first.
