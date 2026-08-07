# 02 results — ProteinGym family-independence audit

**Verdict: ProteinGym is far more sequence-independent than expected. All four
decision rules pass; Experiment 03 gets a genuine cross-family design.**

The preregistered prediction failed, and in the informative direction. I
predicted independent families at 50% identity would fall "plausibly below
100." The actual number is 178. Homology clustering barely collapses the
benchmark.

## The prediction was wrong

| Threshold | Independent units | Largest cluster |
|---|---|---|
| 217 assays (nominal) | 217 | — |
| Protein identity (UniProt) | 186 | 1.8% |
| 90% identity | 184 | 1.8% |
| 70% identity | 179 | 1.8% |
| 50% identity | 178 | 1.8% |
| 30% identity | 172 | 1.8% |

The only real collapse is 217 → 186, and it comes from a handful of proteins
carrying multiple assays (BLAT_ECOLX and P53_HUMAN with four each). Below
that, sequence homology removes almost nothing: 178 clusters at 50%, still 172
at 30%. The largest cluster is 1.8% of the benchmark at every threshold.

ProteinGym's assays are drawn from proteins that are, by sequence, nearly all
distinct systems.

## Coverage sensitivity

mmseqs `easy-cluster` defaults to `-c 0.8` (80% length overlap required),
which could inflate the count by refusing to cluster proteins of very
different lengths even at high identity. Tested:

| 30% identity setting | Clusters |
|---|---|
| c=0.8 (default) | 172 |
| c=0.5 | 170 |
| c=0.0 (no coverage requirement) | 160 |
| c=0.8, cov-mode 1 (target) | 167 |

Even with the coverage requirement removed entirely, 160 clusters remain. The
mild collapse is a property of the benchmark, not of the clustering
parameters. The headline uses the default `-c 0.8`; the floor across settings
is ~160.

## Decision rules — all pass at the 50% primary threshold

1. At least 30 independent clusters: **178** ✓
2. No cluster exceeds 15% of the benchmark: largest is **1.8%** ✓
3. At least two selection classes with five or more clusters: **all five**
   coarse classes qualify ✓
4. Meaningful MSA-depth spread survives: **all three** depth categories
   represented ✓

**Experiment 03 cross-family design is warranted.**

## What this changes for Experiment 01

The 01c finding — that the depth-scaling association is between-assay in the
ProGen ladders — was flagged as possibly a benchmark-composition artefact, on
the worry that a few homologous families might dominate. This audit shows that
worry does not hold: ProteinGym is not composed of a few large families.
Between-assay structure exists, but it is spread across ~178 essentially
independent sequences, not concentrated.

That makes the 01c between-assay signal harder to dismiss as an artefact of a
handful of over-represented families, and correspondingly more likely to be
either a genuine composition effect across many independent proteins, or the
low-power result the short ProGen ladders predict. It does not resolve which,
but it removes one of the three candidate explanations.

## What this cannot establish

Sequence identity is a proxy for independence, not proof of it. Two sequences
below 30% identity can share a fold and a mechanism, so the true number of
mechanistically independent systems may be lower than 178. A structural tier
(Foldseek over predicted or experimental structures) would tighten this and is
the natural extension — recorded as a gap, not run in this pass. The audit
establishes homology-based independence, which is necessary but not sufficient
for statistical independence.

## Artefacts

- `notebooks/02_proteingym_family_audit.py` — Tier A (no tool) + Tier B (mmseqs2)
- `results/02_proteingym_assay_manifest.csv` — per-assay, with cluster IDs at each threshold
- `results/02_independence_summary.csv` — the headline table
- `results/provenance_02_audit.json`
- `figures/02_assays_per_cluster.png`, `figures/02_depth_within_independent_set.png`

## Deviation from the reviewed design

None to the analysis. One addition: a coverage-sensitivity check
(`-c` swept 0.0–0.8) was run post hoc to confirm the mild collapse is not a
clustering-parameter artefact. It is reported above and does not change the
headline, which uses the preregistered default settings.
