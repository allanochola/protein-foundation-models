# ESM-2

Lin et al., *Science* (2023). "Evolutionary-scale prediction of atomic-level protein structure with a language model."

## Trained on

UniRef50 cluster representatives, with training sequences sampled from the corresponding UniRef90 members — roughly 65M unique sequences. Objective is masked language modeling over amino acids, no structure and no MSA at training time. Released at 8M, 35M, 150M, 650M, 3B, and 15B parameters. The 650M model is the default workhorse; most downstream papers benchmark against it.

ESMFold sits on top: the frozen ESM-2 trunk feeds a folding head trained on PDB structures.

## Does well

- **Contact and structure prediction from a single sequence.** ESMFold reaches roughly 0.68 mean TM-score on CAMEO targets, an order of magnitude faster than AlphaFold2 because it skips MSA search.
- **Zero-shot variant effect.** Masked-token pseudo-likelihood correlates with deep mutational scanning fitness without any supervised training.
- **General-purpose embeddings.** Per-residue and per-protein representations transfer to secondary structure, localization, stability, and binding tasks with a linear head.
- **Scale improves structure.** Contact precision and perplexity improve monotonically from 8M to 15B, which is the paper's central scaling claim.

## Limitations

- **The "no MSA" framing oversells it.** Performance tracks the depth of the sequence family the protein belongs to. Orphan proteins and shallow families degrade sharply. The model has not escaped MSA dependence — it has amortized the MSA into weights during pretraining.
- **15B is not uniformly better than 3B** on downstream tasks. Scaling gains show up in perplexity and contacts more reliably than in fine-tuned task accuracy.
- **Below AlphaFold2 on hard targets.** ESMFold trades accuracy for speed; the gap widens on low-homology CASP targets.
- **No function conditioning and no controllable generation.** MLM gives you a scoring function, not a design tool.
- **UniRef50 clustering does not prevent evaluation leakage.** Held-out sequences frequently share fold and remote homology with training data. The paper's splits do not isolate this.

## When not to use it

- **Orphan proteins and shallow families.** Performance tracks family depth despite the MSA-free framing.
- **The 15B model.** On ProteinGym substitutions it loses to 650M on 65% of high-depth assays (mean −0.054 Spearman). Use 650M unless you have a specific reason. See [experiment 01](../experiments/01-msa-depth-confound/PROTOCOL.md).
- **Accuracy-critical structure prediction on low-homology targets**, if you can afford AlphaFold2. ESMFold buys speed, not accuracy.
- **Controllable generation.** It is a scorer, not a design tool.

## Evaluated by

Perplexity on held-out UniRef; long-range contact precision at P@L, P@L/2, P@L/5; TM-score and pLDDT on CAMEO and CASP14; Spearman correlation on deep mutational scanning sets; TAPE-style downstream benchmarks. ESM Atlas — 617M predicted metagenomic structures — is presented as a scale demonstration rather than a benchmark.

## Ideas

- ~~Regress every reported ESM-2 win against MSA depth. Prediction: the scaling curve flattens.~~ **Run, prediction wrong.** Gains are constant across depth strata over the full ladder. The real effect is in the upper segment: past 650M, larger models lose to 650M on 65% of high-depth assays while marginally winning on low-depth ones. See [experiment 01](../experiments/01-msa-depth-confound/PROTOCOL.md).
- Build a strict-holdout benchmark using structural rather than sequence-identity splits (fold-level holdout via Foldseek clustering) and re-rank the ladder.
- Test whether pseudo-likelihood variant scores degrade differently for gain-of-function versus loss-of-function mutations. Most DMS assays are loss-biased, which would make the benchmark blind to exactly the direction that matters for risk assessment.
