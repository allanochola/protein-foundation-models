# Protein Foundation Models

A living, experiment-driven knowledge base for protein and biological foundation models.

> Literature should not be static. Experiments should continuously revise what the literature notes claim.

Model summaries and experiments are not separate projects here. Every experiment feeds back into the documentation: a successful replication raises confidence in a claim, a failed one corrects or narrows it. [Experiment 01](experiments/01-msa-depth-confound/PROTOCOL.md) is the loop working — an idea in `models/esm2.md` was run, the prediction failed, and the note carries the correction with the original wording struck rather than deleted.

## Goals

1. Build structured notes on major protein foundation models.
2. Compare models using reproducible experiments.
3. Separate published claims from independently reproduced findings.
4. Produce reusable Colab experiments.
5. Keep an auditable history showing why every documentation change happened.

## Index

| Model | Year | Modality | Params | Pretraining corpus | Generative |
|---|---|---|---|---|---|
| [ProtT5](models/prott5.md) | 2021 | Sequence | 3B (XL-U50) | BFD (2.1B seqs) → UniRef50 | Weak (span infill) |
| [ProGen](models/progen.md) | 2020 / 2023 | Sequence + tags | 1.2B | ~280M seqs, >19k families | Yes |
| [ESM-2](models/esm2.md) | 2022 | Sequence | 8M – 15B | UniRef50/90, ~65M seqs | No (MLM) |
| [GearNet](models/gearnet.md) | 2023 | Structure graph | ~20M | 805K AlphaFold DB structures | No |
| [SaProt](models/saprot.md) | 2024 | Sequence + 3Di tokens | 650M / 1.3B | ~40M seqs + AFDB structures | No (MLM) |
| [ESM-3](models/esm3.md) | 2024 / 2025 | Seq + struct + function | 1.4B / 7B / 98B | 2.78B proteins, 771B tokens | Yes |
| [Evo 2](models/evo2.md) | 2025 / 2026 | DNA (nucleotide) | 7B / 40B | 9.3T nt, >128k genomes | Yes |

Machine-readable: [models.csv](models.csv). Ancestry: [docs/lineage.md](docs/lineage.md).

## This table is not a leaderboard

These models do not sit on one axis. ESM-2, ProtT5 and SaProt are representation models; ESM-3 and ProGen also generate; GearNet is a graph encoder requiring structures; Evo 2 models DNA at nucleotide resolution.

They overlap on exactly one task. Zero-shot variant effect prediction is scored natively for ESM-2, SaProt, ProGen2 and Evo 2 on the same assays, which makes it the only place a cross-modality number means anything. Everything else is model-specific and reported as such.

## Repository structure

```
├── README.md              this file
├── SETUP.md               Colab ↔ GitHub workflow
├── template.md            entry format for new model notes
├── models.csv             machine-readable comparison matrix
├── references.bib         bibliography
├── evaluation-notes.md    weaknesses shared across the field
│
├── models/                one note per model
├── experiments/           one protocol per experiment
├── notebooks/             canonical .py analysis scripts
├── results/               CSVs and provenance records
├── figures/
├── src/                   shared helpers
└── docs/                  lineage, third-party licenses
```

## Model note format

Six headings, in this order:

1. **Trained on** — corpus, size, objective, parameter count
2. **Does well** — reported wins, and what each actually demonstrates
3. **Limitations** — what the paper does not establish
4. **When not to use it** — concrete negative constraints
5. **Evaluated by** — benchmarks, metrics, splits, baselines
6. **Ideas** — follow-on experiments, each with a status

Heading 6 is the point of the exercise; heading 4 is what readers use most. Format details in [template.md](template.md).

## Provenance tags

Every factual claim carries one.

| Tag | Meaning |
|---|---|
| `[author]` | Reported by the original authors |
| `[repro]` | Independently reproduced |
| `[community]` | Community observation, not formally published |

Filling these in honestly is uncomfortable, and that is the point. Most claims in these notes are `[author]`. A result nobody has reproduced must not read like one three groups have confirmed.

## Idea status

Ideas are tracked as experiments, not brainstorming.

| Status | Meaning |
|---|---|
| `[open]` | Not started |
| `[running]` | Active experiment |
| `[done]` | Resolved — links to the experiment that resolved it |

## Experimental philosophy

Every experiment starts with a protocol containing the question, prediction, analysis plan, and kill criteria — written and committed **before** any analysis runs. Results are interpreted only after the protocol is fixed. When a prediction fails, the protocol is amended with the failure marked and the original wording intact.

## Current experiments

### [Experiment 01 — MSA depth as a confound in pLM scaling](experiments/01-msa-depth-confound/PROTOCOL.md)

**Question.** Do reported ESM-2 scaling gains on zero-shot variant effect prediction survive a control for how well sampled the protein's family is?

**Status: pilot complete, kill criteria outstanding.**

The original prediction failed. Across the full 8M–15B ladder the size × depth interaction is null (p = 0.30) — scaling helps roughly equally regardless of family depth. The effect lives in the upper segment: from 650M up the interaction is significant (β = −0.015, p = 0.002) while the main effect of size is not (p = 0.11). Paired per assay, 15B beats 650M on 61% of low-depth assays but only 35% of high-depth ones.

Past 650M, scaling ESM-2 degrades performance on well-sampled families and helps marginally on poorly-sampled ones. A single mean Spearman across ProteinGym hides that sign flip.

Not yet a result. The covariate controls and the ProGen2 replication are still to run.

## Reproducibility

```
GitHub  →  Colab  →  experiment  →  results  →  documentation update  →  commit
```

GitHub is the permanent record. Colab is temporary compute. Full workflow in [SETUP.md](SETUP.md).

Canonical experiment files are Python scripts in jupytext percent format. `.ipynb` is a generated artefact and is gitignored — edit the `.py` and regenerate, not the reverse.

Analysis inputs are pinned by commit SHA in the code that fetches them, and the same variable is written into the provenance record. Recording a version while fetching from a branch documents intent, not what the code consumed.

Excluded from the repository: model checkpoints, cached embeddings, generated notebooks, downloaded datasets. Experiments should be reproducible from scripts, not from stored artefacts.

## Reading order

Start with ProtT5 and ESM-2 to fix the sequence-only baseline. Add GearNet next as the control: structure supervision, no language modeling, 100× less data, comparable function prediction. Then SaProt and ESM-3 as two different answers to "add structure." Finish with ProGen and Evo 2, which move from representation to generation and change what evaluation has to prove.

Cross-cutting: [evaluation weaknesses shared across the field](evaluation-notes.md).

## Scope and safety

This repository evaluates general capabilities of public protein and genomic foundation models using published benchmark data. It contains no experimental protocols, no wet-lab methods, and no work directed at optimizing biological agents. Where a model's dual-use profile is discussed, the aim is testing whether stated safety mitigations do what they claim.

## Roadmap

- [x] Structured model notes
- [x] Experiment-driven documentation loop
- [x] Experiment 01a — published-score analysis
- [ ] Experiment 01b — covariate controls (taxon, length, selection type)
- [ ] Experiment 01c — ProGen2 ladder replication
- [ ] Experiment 02 — structure-aware models under pLDDT stratification
- [ ] Experiment 03 — generative model scoring
- [ ] Experiment 04 — robustness and out-of-distribution evaluation

## Citation

Cite the specific GitHub release corresponding to the experiments relied upon.

## License

Repository code and notes: MIT ([LICENSE](LICENSE)).

Model weights carry their own terms, some of which propagate into derivative work — ESM-3's non-commercial license is the one to plan around. See [docs/third-party-licenses.md](docs/third-party-licenses.md) before adding a model.
