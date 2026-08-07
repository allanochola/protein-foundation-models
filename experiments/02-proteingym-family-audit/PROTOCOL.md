# Experiment 02 — ProteinGym family-independence audit (PREREGISTERED)

**Status: protocol only. No 02 analysis has been run. This file is committed
before any audit code.**

## Why this experiment, now

Experiment 01 found a depth-scaling association that is within-assay in ESM-2
but between-assay in ProGen2 and ProGen3 — consistent with benchmark
composition rather than a scaling mechanism. That interpretation rests on an
untested assumption: that ProteinGym's 217 assays are enough independent
biological systems to support cross-family inference at all.

This experiment tests that assumption directly. It is a feasibility audit, not
a model benchmark. It measures how many statistically independent systems
ProteinGym actually contains, so that Experiment 03 (and any re-reading of 01)
knows whether cross-family claims are supported or whether a mechanistic
case-study design is required instead.

## Question

How many statistically independent biological systems does ProteinGym's
substitution benchmark contain, under a range of defensible homology
definitions?

## Primary prediction

The nominal count of 217 assays and ~2.5M mutations substantially overstates
effective sample size. Assays cluster within proteins and within homologous
families, so the number of independent units is materially smaller than 217 —
and smaller still than the mutation count implies.

Stated to fail: I expect independent families under a 50% homology threshold
to number well below 217, plausibly below 100.

## Data

ProteinGym reference file `reference_files/DMS_substitutions.csv` at the pinned
commit `144fe22b07dfaeec2b366f2346203a9838a55b4c` — the same commit used
throughout Experiment 01. 217 assays, 46 columns, all `target_seq` present.
No new downloads for the identity and metadata tiers.

## Two tiers, by what the data supports

**Tier A — from the reference file directly, no external tool.**
Already computable: assays per UniProt ID, taxon distribution, MSA depth
distribution, selection-type distribution. The protein-identity row of the
headline table comes from here.

**Tier B — sequence homology clustering, requires a tool.**
The 90/70/50/30% identity thresholds and any structural grouping need a
clustering step over the 217 `target_seq` values (86k residues total — a
trivial job). Tool: `mmseqs2 easy-cluster` at each `--min-seq-id`. If mmseqs
is unavailable in the run environment, CD-HIT is an acceptable substitute and
the substitution is recorded. Structural family (Foldseek/SCOP) is **optional
and out of scope for the first pass** — recorded as a gap, not blocked on.

Missingness is explicit. A threshold that cannot be computed is reported as
"not run," never silently omitted or inferred.

## Manifest — one row per assay

`results/02_proteingym_assay_manifest.csv`, columns:

- assay_id, uniprot_id, protein_name, organism, taxon
- assay_type (coarse_selection_type), selection_type
- seq_len, n_single_mutants, n_total_mutants
- msa_neff_l, msa_neff_l_category
- homology_cluster_90, _70, _50, _30   (Tier B; "NA" if not run)
- structural_family                     (optional; "NA" if not run)
- in_experiment_01                      (all True — same 217)

## Headline table

`results/02_independence_summary.csv`:

| Threshold | Assays | Unique proteins | Independent clusters | Median assays/cluster | Largest cluster % |
|---|---|---|---|---|---|
| Protein identity (UniProt) | 217 | (Tier A) | (Tier A) | | |
| 90% identity | 217 | | (Tier B) | | |
| 70% identity | 217 | | (Tier B) | | |
| 50% identity | 217 | | (Tier B) | | |
| 30% identity | 217 | | (Tier B) | | |

The load-bearing number is the independent-cluster count, not the mutation
count. Mutation totals are reported once, for contrast, and never presented as
statistical power.

## Decision rules — fixed before running

A cross-family statistical study (Experiment 03) is warranted only if, at the
50% identity threshold (the primary threshold):

1. At least 30 independent clusters remain.
2. No single cluster exceeds 15% of clusters' total effective weight.
3. At least two coarse selection-type classes each have five or more
   independent clusters.
4. Meaningful spread in MSA depth survives within the independent set (not all
   clusters in one depth category).

If any fails: do not present mutation-level or assay-level counts as
statistical power. Switch Experiment 03 to a hybrid design — broad behavioural
reproduction plus mechanistic case studies on a small, deliberately chosen set
of proteins.

Thresholds are conventions for this repository, not universal laws, and may be
amended — but only in a commit made before the audit runs.

## Kill criteria

- `target_seq` incomplete or malformed → Tier B cannot run; report Tier A only
  and mark Tier B "not run."
- No clustering tool installable in the environment → Tier B "not run," and the
  50% decision rule cannot be evaluated; the experiment reports Tier A and
  stops short of the Experiment 03 go/no-go.

## Outputs

- `notebooks/02_proteingym_family_audit.py`
- `results/02_proteingym_assay_manifest.csv`
- `results/02_independence_summary.csv`
- `results/provenance_02_audit.json`
- `figures/02_assays_per_cluster.png`
- `figures/02_depth_within_independent_set.png`

## What this cannot establish

Sequence-identity clustering is a proxy for independence, not a guarantee of
it. Two proteins below 30% identity can share a fold and a mechanism; two
above 90% can behave differently under a given assay. A structural tier
(Foldseek) would tighten the lower thresholds and is the natural extension.
The audit measures homology-based independence, which is necessary but not
sufficient for statistical independence.
