# <Model name>

## Overview

| | |
|---|---|
| Organization | |
| Publication | <venue, year, arXiv/DOI> |
| Sizes | |
| Architecture | |
| Modality | <sequence / sequence+structure / multimodal / DNA> |
| Objective | <MLM / autoregressive / contrastive / masked-multimodal> |
| Open weights | <yes / partial — which sizes / no> |
| Code + checkpoints | <repo, HuggingFace> |

Open weights is a governance field, not a convenience field. Record which
sizes are public and whether the headline results come from a public model.

## Trained on

<Corpus, sequence count, structural and functional annotations, synthetic
data, curation and exclusions. Note exclusions explicitly — they are usually
safety claims and are usually untested.>

## Does well

- <Task — the reported number, the benchmark it came from, and what the
  number actually demonstrates.>

## Limitations

**Acknowledged by the authors**
- <>

**Not addressed**
- <What the paper does not establish. This is where the research questions
  come from.>

## When not to use it

- <Concrete negative constraint. More decision-useful than the strengths
  list, and easier to falsify.>

## Downstream performance

| Task | Metric | Benchmark | Split construction | Baseline | Result |
|---|---|---|---|---|---|
| | | | | | |

Tag the provenance of every row: `[author]` reported by the authors,
`[repro]` independently reproduced, `[community]` observed but not formally
published. A claim nobody has reproduced must not read like one three groups
have confirmed.

No ordinal labels. "Excellent" is not contestable; 0.68 mean TM-score on
CAMEO against AlphaFold2's 0.88 is. If a number is unavailable, write
"not reported" — that absence is itself a finding.

## Evaluated by

<Benchmarks, metrics, splits, baselines, zero-shot vs fine-tuned. Name the
split construction. Most claims break here.>

## Ideas

Each entry carries a status. The failure mode is not that ideas are
unsorted — it is that they never get linked back to whether anyone ran them.

- `[open]` <Experiment you could run this month, with your prediction.>
- `[running]` <In progress — link the experiment directory.>
- `[done]` <Result, and a link. If the prediction failed, say so here and
  strike the original wording rather than deleting it.>

## Links

<Original paper · follow-ups · benchmark papers that include this model ·
official repo · checkpoints. Full citations go in references.bib.>
