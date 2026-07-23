# PROTOCOL

Every experiment in this repository begins with a written protocol.

The protocol is completed **before any analysis is run**.

Its purpose is to separate prediction from interpretation.

If an experiment changes direction after results are known, the protocol should be updated only in a new commit with an explicit explanation.

---

# Scientific Principles

Experiments should be

- reproducible
- falsifiable
- minimal
- auditable

The goal is not to confirm expectations.

The goal is to reduce uncertainty.

---

# Experiment Lifecycle

Every experiment follows the same sequence.

```
Question
    ↓
Background
    ↓
Hypothesis
    ↓
Prediction
    ↓
Experimental Design
    ↓
Stopping Criterion
    ↓
Run Analysis
    ↓
Interpret Results
    ↓
Update Documentation
```

Interpretation never occurs before the analysis has completed.

---

# Experiment Template

Every experiment directory contains a `PROTOCOL.md`.

---

## Title

A concise description of the experiment.

Example

```
Experiment 01

MSA Depth and Reported ESM-2 Scaling
```

---

## Motivation

Why is the question interesting?

What published claim is being examined?

What uncertainty exists?

---

## Research Question

A single sentence.

Example

```
Does MSA depth explain the reported scaling behaviour of ESM-2?
```

---

## Background

Summarize the relevant literature.

Each statement should carry a provenance tag.

Example

```
Scaling improves with model size. [author]

MSA depth may influence evaluation. [community]
```

---

## Hypothesis

A testable statement.

Example

```
MSA depth explains a substantial fraction of the apparent scaling trend.
```

---

## Prediction

A concrete observable outcome.

Example

```
After controlling for MSA depth,
the model-size effect should shrink substantially.
```

Predictions should be quantitative whenever possible.

---

## Null Hypothesis

Explicitly state the alternative.

Example

```
MSA depth does not explain the reported scaling trend.
```

---

## Variables

Independent variables

Dependent variables

Control variables

Confounders

Example

```
Independent

Model size

MSA depth

Dependent

Reported performance

Controls

Assay

Protein family
```

---

## Dataset

Describe

- source
- version
- commit SHA
- download date

Never rely on

```
main
```

for external repositories.

Pin everything.

---

## Experimental Design

Describe

- preprocessing
- filtering
- exclusions
- statistical model
- evaluation procedure

Someone unfamiliar with the project should be able to reproduce the experiment.

---

## Expected Failure Modes

Examples

- missing data

- duplicated proteins

- inconsistent annotations

- dataset leakage

- confounding

Write these before running the experiment.

---

## Stopping Criterion

The experiment ends when

- hypothesis supported

or

- hypothesis contradicted

or

- assumptions violated

Avoid "collect more data until something becomes significant."

---

## Output Files

Specify expected outputs.

Example

```
results/

figures/

logs/

provenance.json
```

---

## Documentation Impact

List every document that may require updating.

Example

```
models/esm2.md

models.csv

README.md
```

---

## Reproducibility Checklist

Before analysis begins

- [ ] Protocol committed

- [ ] Dataset pinned

- [ ] Dependencies installed

- [ ] Random seed fixed

- [ ] Output directory created

- [ ] Provenance enabled

---

# After the Experiment

The following sections remain empty until analysis is complete.

---

## Results

Summarize

- numerical findings

- figures

- statistical significance

No discussion yet.

---

## Interpretation

Did the prediction survive?

Was the hypothesis supported?

Were assumptions violated?

---

## Documentation Changes

Record every documentation update.

Example

```
esm2.md

Changed

"Scaling improves uniformly."

↓

"Scaling improvement weakens after controlling for MSA depth."
```

---

## Lessons Learned

Document

- mistakes

- unexpected observations

- future experiments

---

# Provenance

Every experiment records

- Git commit

- experiment ID

- timestamp

- dataset commit SHA

- software versions

- operating system

- Python version

- package versions

- random seed

The provenance record should be machine-readable.

---

# Repository Philosophy

Experiments are primary.

Documentation is secondary.

Documentation changes because experiments produce evidence.

Evidence never changes because documentation expected it.
