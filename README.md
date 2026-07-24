# Protein Foundation Models — Reading Notes

Structured notes on protein and biological sequence foundation models, and the experiments they generate. Every entry answers the same questions so models stay comparable across architectures, modalities, and training scales.

Model files use these headings, in this order:

1. **Trained on** — corpus, size, objective, parameter count
2. **Does well** — tasks with reported wins, and what the win actually shows
3. **Limitations** — what the paper does not establish
4. **When not to use it** — concrete negative constraints
5. **Evaluated by** — benchmarks, metrics, splits, baselines
6. **Ideas** — follow-on experiments, each with a status

Heading 6 is the point of the exercise, and heading 4 is the part readers use most. Notes that stop at heading 5 are summaries; notes that reach heading 6 and link results back are a research agenda. See [experiment 01](experiments/01-msa-depth-confound/PROTOCOL.md) for the loop working: an idea from `esm2.md` was run, the prediction failed, and the note was corrected.

## Index

| Model | Year | Modality | Params | Pretraining corpus | Generative |
|---|---|---|---|---|---|
| [ProtT5](models/prott5.md) | 2021 | Sequence | 3B (XL-U50) | BFD (2.1B seqs) → UniRef50 | Weak (span infill) |
| [ProGen](models/progen.md) | 2020 / 2023 | Sequence + tags | 1.2B | ~280M seqs, >19k families | Yes |
| [ESM-2](models/esm2.md) | 2022 | Sequence | 8M – 15B | UniRef50/90, ~65M seqs | No (MLM) |
| [GearNet](models/gearnet.md) | 2023 | Structure graph | ~20M | 805K AlphaFold DB structures | No |
| [SaProt](models/saprot.md) | 2024 | Sequence + 3Di tokens | 650M | ~40M seqs + AFDB structures | No (MLM) |
| [ESM-3](models/esm3.md) | 2024 / 2025 | Seq + struct + function | 1.4B / 7B / 98B | 2.78B proteins, 771B tokens | Yes |
| [Evo 2](models/evo2.md) | 2025 / 2026 | DNA (nucleotide) | 7B / 40B | 9.3T nt, >128k genomes | Yes |

Machine-readable version: [models.csv](models.csv).

## Cross-cutting notes

- [Lineage of the reading set](docs/lineage.md) — who inherited what, and which comparisons are controlled.
- [Evaluation weaknesses shared across the field](evaluation-notes.md) — homology leakage, MSA-depth confounds, and why benchmark wins keep failing to transfer.

## Scope and safety

This repository evaluates general capabilities of public protein and genomic
foundation models using published benchmark data. It contains no experimental
protocols, no wet-lab methods, and no work directed at optimizing biological
agents. Where a model's dual-use profile is discussed — see
[evo2.md](models/evo2.md) — the aim is testing whether stated safety
mitigations do what they claim.

Model weights carry their own licenses, some of which propagate into
derivative work. See [docs/third-party-licenses.md](docs/third-party-licenses.md)
before adding a model.

## Comparability

These models do not sit on one axis, and the index table above should not be
read as a leaderboard. ESM-2, ProtT5 and SaProt are representation models;
ESM-3 and ProGen also generate; GearNet is a graph encoder requiring
structures; Evo 2 models DNA at nucleotide resolution.

They do overlap on one task. Zero-shot variant effect prediction is scored
natively for ESM-2, SaProt, ProGen2 and Evo 2 on the same assays, which makes
it the only place a cross-modality number means anything. Everything else is
model-specific and reported as such.

## Project phases

1. **Model survey** — what each model was trained on, what it does, where it breaks. This is the current phase.
2. **Unified benchmarking** — run the models on shared tasks with corrected splits. ProteinBench and PFMBench are the starting points; neither controls for homology leakage, which is the gap.
3. **Research contribution** — an evaluation method that separates capability from coverage. See [evaluation-notes.md](evaluation-notes.md).

## Conventions

- **Filenames are slugs, not numbers.** Ordering lives in the index table above. Numbered prefixes force a renumber every time a model is inserted, which breaks links and shared permalinks.
- **No ordinal performance labels.** "Excellent" and "moderate" cannot be contested. Every performance claim carries a metric, a benchmark, a split, and a baseline, or it is marked "not reported."
- **Limitations separate acknowledged from unaddressed.** The unaddressed ones are where the research questions live.
- **Every entry ends in Ideas, and every idea carries a status** — `[open]`, `[running]`, `[done]` with a link. Ideas that never get linked to whether they ran are how notes go stale.
- **"When not to use it" is required.** Negative constraints are more decision-useful than strengths and easier to falsify.
- **Tag claim provenance** — `[author]`, `[repro]`, `[community]`. An unreproduced result must not read like a confirmed one.
- Full citations in [references.bib](references.bib). Entry format in [template.md](template.md).

## Reading order

Start with ProtT5 and ESM-2 to fix the sequence-only baseline. Add GearNet next as the control: structure supervision with no language modeling, 100× less data, comparable function prediction. Then read SaProt and ESM-3 as two different answers to "add structure." Finish with ProGen and Evo 2, which move from representation to generation and change what evaluation has to prove.
