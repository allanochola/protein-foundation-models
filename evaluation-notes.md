# Cross-cutting: what the evaluations share, and where they break

Seven models, four modalities, five orders of magnitude in parameter count — and the same five evaluation weaknesses appear in all of them.

## 1. Homology leakage survives sequence-identity splits

Every model in this set uses sequence-identity clustering (UniRef50, 30% identity thresholds) to build held-out sets. Sequence identity is the wrong axis. Two proteins at 15% identity can share a fold, a function, and a binding site. Fold-level or Foldseek-based structural holdout is the correct split and almost nobody uses it.

Consequence: reported generalization is partly memorized homology, and the size of that fraction is unmeasured across the board.

## 2. MSA depth is an unmeasured confound

Sequence-only models are marketed as MSA-free. Their performance still tracks how many homologs exist for the target family. The MSA has been amortized into the weights, not eliminated. No paper in this set reports performance stratified by family depth, which means "works without an MSA" and "works on proteins with no homologs" are being conflated.

Fix: regress task performance on log family size. Report the residual. This is a one-figure addition to any paper and it would change how several headline claims read.

## 3. Wet-lab validation has no denominator

ProGen reports functional lysozymes. ESM-3 reports esmGFP. Neither reports the funnel: candidates generated, candidates selected, candidates synthesized, candidates that worked, and the same numbers for a baseline method. A success rate without a denominator is an anecdote, and capability claims built on anecdotes cannot be compared across models or tracked over time.

This matters beyond publication hygiene. Governance frameworks that key on demonstrated design capability need a measurable quantity, and the field is not currently producing one.

## 4. Benchmarks are loss-of-function biased

Deep mutational scanning assays overwhelmingly measure disruption. Models are therefore evaluated on their ability to predict what breaks a protein. Gain-of-function prediction — the direction that matters most for both engineering value and misuse risk — is barely represented in ProteinGym or its relatives. A model could be excellent on every published variant benchmark and uninformative about the cases anyone actually worries about.

## 5. Predicted structures contaminate structure-aware models

SaProt, GearNet, and ESM-3 all train on AlphaFold predictions. SaProt's own loss curves show validation loss on AF2 structures diverging from loss on real PDB structures. These models learn AlphaFold's error distribution as if it were signal, then get evaluated on tasks whose labels partly derive from the same predictions. The circularity is unquantified.

## What the unified suites fix, and what they don't

Two suites now attempt cross-model comparison, and they are complementary rather than interchangeable.

**ProteinBench** (Ye et al., ICLR 2025) covers generative challenges — protein design, structure prediction, conformational dynamics — and scores on four axes: quality, novelty, diversity, robustness. Separating novelty from quality is the right instinct and directly addresses the "is this model inventing or recovering?" problem.

**PFMBench** (Gao et al., 2025) runs 38 tasks across 8 areas on 17 models. Its most useful design decision is normalizing to roughly 1B parameters where multiple sizes exist, which controls the one confound the field most often leaves free. It also reports inter-task correlations, so you can see which of the 38 tasks are measuring the same thing.

Neither fixes items 1 through 5 above. Both inherit the underlying datasets' splits, so homology leakage passes straight through. Neither stratifies by MSA depth or by AF2 confidence. Neither addresses the loss-of-function bias, because they draw from the same DMS assays. A leaderboard built on contaminated splits ranks models confidently and wrongly, and it does so at greater scale than the individual papers did.

Worth tracking: at least one contamination-free benchmark effort has appeared (LiveProteinBench, arXiv 2512.22257) built around holding out post-cutoff depositions. Verify the details before citing — this note is from a passing reference, not a full read.

## What to build

A single evaluation harness with:

- structural rather than sequence-identity holdout
- performance stratified by MSA depth and by AF2 pLDDT
- reported denominators for any wet-lab claim
- explicit gain-of-function coverage
- the same tasks run across sequence, structure, and DNA models so modality comparisons are meaningful

Every model in this reading set could be re-scored under it. The most likely outcome is that several reported scaling advantages shrink, and GearNet — the smallest model here — looks better than the ranking currently suggests.

That result would matter for oversight, not just benchmarking. If the field cannot tell which reported capabilities are real, then neither developers nor regulators can contest a capability claim, and contestability is the precondition for the oversight to mean anything.
