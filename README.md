# Protein Foundation Models

A living, experiment-driven knowledge base for protein and biological foundation models.

This repository is built around a simple principle:

> Literature should not be static. Experiments should continuously revise what the literature notes claim.

Rather than treating model summaries and experiments as separate projects, every experiment feeds back into the documentation. A successful replication increases confidence in a claim. A failed replication corrects or narrows the claim. Documentation and experiments therefore evolve together.

---

# Goals

1. Build structured notes on major protein foundation models.
2. Compare models using reproducible experiments.
3. Distinguish published claims from independently reproduced findings.
4. Produce reusable Colab experiments.
5. Maintain an auditable history showing why every documentation change occurred.

---

# Repository Structure

```
protein-foundation-models/
│
├── README.md
├── SETUP.md
├── PROTOCOL.md
├── requirements.txt
├── models.csv
├── third-party-licenses.md
│
├── models/
│
├── experiments/
│
├── notebooks/
│
├── results/
│
├── figures/
│
└── src/
```

---

# Models Covered

Current scope:

- ESM-2
- ESM-3
- ProtT5
- ProGen
- Evo 2
- SaProt
- GearNet

Each model note follows the same template.

- Overview
- What it was trained on
- Tasks it performs well
- Limitations
- Evaluation
- When not to use it
- Ideas
- References

---

# Provenance Tags

Every factual claim should be tagged.

| Tag | Meaning |
|------|----------|
| **[author]** | Reported by the original authors |
| **[repro]** | Independently reproduced |
| **[community]** | Community observation |

Example

```
Zero-shot mutation prediction performs well. [author]

Replicated on ProteinGym. [repro]

Several users report instability on long proteins. [community]
```

---

# Idea Status

Ideas are tracked as experiments rather than brainstorming.

| Status | Meaning |
|---------|---------|
| **[open]** | Not started |
| **[running]** | Active experiment |
| **[done]** | Completed |

Completed ideas always link to the experiment that resolved them.

Example

```
[done]

Experiment 01

Prediction contradicted.

See:
experiments/01-msa-depth-confound/
```

---

# Experimental Philosophy

Every experiment begins with a protocol.

The protocol contains

- question
- hypothesis
- prediction
- stopping criterion
- analysis plan

before any analysis is run.

Results are interpreted only after the protocol has been fixed.

---

# Current Experiments

## Experiment 01

Question

Does MSA depth explain the reported scaling behaviour of ESM-2?

Outcome

The interaction was statistically significant but opposite to the original prediction.

This resulted in corrections to the ESM-2 documentation.

---

# Reproducibility

The repository follows one workflow.

```
GitHub
    ↓

Google Colab

    ↓

Experiment

    ↓

Results

    ↓

Documentation update

    ↓

Git commit
```

GitHub is the permanent record.

Colab is temporary compute.

---

# Canonical Files

The canonical experiment files are Python scripts.

Interactive notebooks are generated from them using Jupytext when needed.

```
.py
↓

.ipynb
```

The Python files are the source of truth.

---

# Large Files

The repository intentionally excludes

- model checkpoints
- cached embeddings
- temporary notebooks
- downloaded datasets

Experiments should be reproducible from scripts rather than stored artifacts.

---

# Citation

If you use this repository in academic work, please cite the specific GitHub release corresponding to the experiments you relied upon.

---

# License

Repository code is released under the MIT License.

Individual models retain their own licenses.

See `third-party-licenses.md`.

---

# Roadmap

- [x] Structured model notes
- [x] Experiment-driven documentation
- [x] Experiment 01
- [ ] Experiment 02 – Frozen embedding comparison
- [ ] Experiment 03 – Structure-aware models
- [ ] Experiment 04 – Generative protein models
- [ ] Experiment 05 – Robustness and out-of-distribution evaluation
