# Experiment 05 — Results

**Experiment 05 is unresolved and triggers Stage 3. The preregistered H_power
prediction failed: restoring the already-scored lower ProGen checkpoints recovered
no negative within-assay depth×scaling slope on either ladder. But the ProGen
upper ladders are underpowered for an ESM-2-sized effect — the 80% minimum
detectable effect is −0.0248 (ProGen2) and −0.0266 (ProGen3), both larger in
magnitude than ESM-2's −0.0153, giving only 38.6% and 31.0% power at that
magnitude — so the nulls cannot establish absence. ESM-2 keeps essentially the
same coefficient under both downward ladder extension and severe short-ladder
handicaps. Analyses A and C lean architecture; B keeps power live. By the frozen
decision table, that split is unresolved and warrants Stage 3.**

Stage 2 is all-CPU on scores already in the repo. ProteinGym scored every released
ProGen size on the same 217 assays; 01c kept only the upper segment, so the lower
checkpoints were unused within-assay leverage, restored here at zero GPU cost. The
estimator is the 01/01c fixed-effects model verbatim (`rho ~ lp + lp:ld +
C(assay)`, cluster-robust SE by `UniProt_ID`), applied to progressively longer
ladders. Preregistered in PROTOCOL.md and frozen (commit 4f70a03, corrected 78dde69
after a pre-freeze dry run) before any Stage 2 coefficient was read; run at 2000
replicates, seed 0, ProteinGym pin 144fe22b. The fast within-transform estimator
used in the injection loop matches statsmodels 0.14.6 to machine precision, gated
with a hard stop on mismatch.

## Analysis A — restore the discarded leverage

The ESM-2 anchor calibrates whether extension preserves a real effect; the ProGen
ladders are the test. Effect-size ratios `beta / beta_ESM` (−0.0153) are the
primary read; the 95% CI is the assay-cluster bootstrap on `beta`.

| Ladder | segment | pts | log-span | beta (lp:ld) | p | 95% CI (boot) |
|---|---|---|---|---|---|---|
| ESM-2 | 650M+ (upper) | 3 | 1.36 | −0.0153 | 0.011 | [−0.0247, −0.0074] |
| ESM-2 | 35M+ | 5 | 2.63 | −0.0156 | 0.007 | [−0.0241, −0.0075] |
| ESM-2 | 8M+ (full) | 6 | 3.27 | −0.0055 | 0.34 | [−0.0133, +0.0029] |
| ESM-2 | segmented `lp_hi:ld` | 6 | 3.27 | **−0.0213** | 0.0002 | [−0.0326, −0.0099] |
| ProGen2 | 764M+ (upper) | 3 | 0.92 | +0.0064 | 0.60 | [−0.0116, +0.0246] |
| ProGen2 | 151M+ (full) | 4 | 1.63 | +0.0010 | 0.89 | [−0.0115, +0.0140] |
| ProGen2 | segmented `lp_hi:ld` | 4 | 1.63 | +0.0064 | 0.58 | [−0.0163, +0.0291] |
| ProGen3 | 762M+ (upper) | 3 | 0.59 | −0.0031 | 0.77 | [−0.0179, +0.0128] |
| ProGen3 | 339M+ | 4 | 0.95 | −0.0126 | 0.44 | [−0.0332, +0.0075] |
| ProGen3 | 112M+ (full) | 6 | 1.43 | +0.0031 | 0.82 | [−0.0153, +0.0230] |
| ProGen3 | segmented `lp_hi:ld` | 6 | 1.43 | −0.0104 | 0.38 | [−0.0335, +0.0128] |

The ESM-2 anchor survives extension. The coefficient holds from the 3-point upper
ladder (−0.0153) down through 35M+ (−0.0156, 5 points, 2.63 decades), then washes
only when the two smallest checkpoints enter at 8M+ (−0.0055, ns) — expected, since
those models sit in the low-capability regime. The segmented fit localises it: the
upper-regime slope `lp_hi:ld` is −0.0213 (ratio 1.39, p = 0.0002), the low-regime
`lp_lo:ld` +0.0046 (ns). Extension preserves a real effect and the segmentation
credits it to the upper regime. The anchor is informative.

ProGen recovers nothing. Every ProGen rung is null. ProGen2 stays positive-to-zero
across its ladder (+0.0064 → +0.0010); ProGen3 wanders without significance
(−0.0031, −0.0126, −0.0038, +0.0031); both segmented upper-regime slopes are
non-significant (ProGen2 +0.0064, p = 0.58; ProGen3 −0.0104, p = 0.38). Every
bootstrap CI spans zero. The ratio to ESM-2's effect is near zero or wrong-signed
with intervals crossing zero — no recovery, on either ladder, at any length.

## Analysis B — injection MDE (upper segment is the decision input)

| Ladder | segment | 80% MDE | power at −0.0153 |
|---|---|---|---|
| ProGen2 | upper | −0.0248 | 0.386 |
| ProGen2 | full | −0.0170 | 0.709 |
| ProGen3 | upper | −0.0266 | 0.310 |
| ProGen3 | full | −0.0250 | 0.377 |

Conservative noise source; the clean/ceiling source agrees within 0.001.

The upper ladders — the segments 01c actually used — cannot see an ESM-sized
effect. Both need a slope near −0.025 to −0.027 for 80% detection, larger in
magnitude than the −0.0153 in question, and both sit near one-third power at that
magnitude. Even the full ProGen2 ladder reaches only 70.9%. A null from a design
with this little power does not distinguish "no effect" from "an effect it cannot
see."

Calibration holds, so B stands. On the clean-null construction the DGP
false-positive rate is 0.0195–0.028 across cells (one-sided, target 0.025) — a
valid null. The same draws scored with the small-sample-corrected SE the analysis
uses give 0.007–0.010: the clustered test is conservative on these short
within-assay ladders, so the reported power is a floor and the MDE a mild upper
bound. Neither is a kill trigger; both are recorded.

## Analysis C — ESM-2 handicap (corroboration only)

| Sub-ladder | pts | log-span | beta | p | 95% CI (boot) |
|---|---|---|---|---|---|
| 650M/3B/15B | 3 | 1.36 | −0.0153 | 0.011 | [−0.0247, −0.0074] |
| 650M→3B | 2 | 0.66 | −0.0162 | 0.078 | [−0.0286, −0.0051] |
| 3B→15B | 2 | 0.70 | −0.0144 | 0.090 | [−0.0260, −0.0047] |
| 650M→15B | 2 | 1.36 | −0.0153 | 0.027 | [−0.0246, −0.0074] |

ESM-2 keeps its coefficient when handicapped to ProGen-sized leverage. Every
two-point sub-ladder returns essentially the full −0.0153 (ratios 0.94–1.06), and
every bootstrap CI excludes zero. The 0.66-decade handicap — shorter than
ProGen2's upper span — still returns −0.0162. The analytic p-values weaken to
0.078–0.090 on the two-point ladders, but that is the mechanical SE inflation of
two points, not a vanishing effect; the bootstrap intervals, which the protocol
reads as primary, exclude zero throughout. C is corroboration, not a co-equal test
— ESM-2's spacing forbids the 3-point short-span ladders ProGen has — but it points
one way: the effect does not need ProGen's leverage to appear.

## Decision (frozen table)

Row 1 (H_power confirmed) required a significant negative ProGen `lp_hi:ld`
carried by the upper regime. It did not appear. **The preregistered prediction
failed.**

Row 2 (H_arch, bounded) required the ProGen upper MDE to reach −0.0153 at ≥ 0.8
power. It does not — power is 0.31–0.39. So the null cannot be read as "no
within-assay effect as large as ESM-2's"; the design never had the resolution to
make that claim.

What remains is the split the final table rows name: A and C show ESM-2's effect
surviving both extension and handicap while ProGen shows nothing, but B shows the
ProGen nulls are underpowered for the effect in question. Architecture is the
leaning explanation; power is not excluded. **Unresolved — Stage 3 is warranted.**

## What this establishes

- The H_power mechanism as originally bet — that the upper ProGen ladders hid a
  real effect more checkpoints would reveal — is not supported. Adding every scored
  ProGen checkpoint moved neither ladder off null.
- The ESM-2 within-assay interaction is not an artefact of long, wide-spanned
  ladders. It survives compression to 2–3 points and spans as short as 0.66
  decades, and segmentation carries it in the upper regime.
- The estimator is faithful: the full-set ESM-2 fit reproduces −0.0153 / p = 0.011
  exactly, and the fast injection estimator matches statsmodels to machine
  precision.

## What this does not establish

- **A ProGen null is not ProGen absence.** The upper-ladder MDEs (0.025–0.027)
  exceed the effect they would test for (0.0153), so the nulls bound nothing
  useful. This is the reason Stage 3 exists.
- No claim about a ProGen effect *smaller* than the MDE, either direction. The
  conservative and clean noise sources agree within 0.001, which does not change
  this.
- C cannot carry an architecture verdict alone. It is corroboration built on
  ESM-2's spacing, not a matched test on ProGen's geometry.
- The conservative real-test type-I (0.007–0.010) means the whole B analysis leans
  against detection; if anything it understates ProGen power. That does not rescue
  the nulls, but it is stated for honesty.
- Scope is the within-assay depth×scaling slope only. ProGen's between-assay
  effects (the 01c confirmatory −0.0245 / −0.0283) are not in question here and are
  not evidence about the within-assay slope.

## Outputs

- `notebooks/05_architecture_vs_power.py`
- `results/05_ladder_extent_sensitivity.csv`
- `results/05_injection_power.csv`
- `results/05_calibration.csv`
- `results/05_esm2_handicap.csv`
- `results/provenance_05_architecture_vs_power.json`
- `figures/05_lp_ld_by_ladder_length.png`, `figures/05_injection_power_curve.png`

## Next

Stage 3 — fresh inference on ProGen with controlled leverage, the only thing that
can separate architecture from power now that the scored ladders are exhausted. The
model and weights choice is deliberately not made in this document: the
interpretation above is frozen against commit 292de3a first, so the Stage 3 design
cannot back-fit it.
