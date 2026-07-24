# Evo 2

Brixi et al., Arc Institute / NVIDIA / Stanford / UC Berkeley / UCSF. bioRxiv (February 2025), published in *Nature* (March 2026). "Genome modeling and design across all domains of life with Evo 2." Successor to Evo (Nguyen et al., *Science*, November 2024).

## Trained on

OpenGenome2: 9.3 trillion nucleotides from more than 128,000 genomes across bacteria, archaea, and eukaryotes. Two model sizes, 7B and 40B parameters, built on StripedHyena 2 rather than a transformer — the earlier Evo paper showed transformers scale poorly at byte-level nucleotide resolution.

Context window of 1 million tokens at single-nucleotide resolution, reached through two-phase training: short-context pretraining, then long-context midtraining with data weighting that prioritizes functional elements. Training code, weights, and the full dataset are open — currently the largest fully open model release in the field.

Note the modality shift. Evo 2 is a DNA model. Protein capability is downstream of genome modeling, not the training target.

## Does well

- **Zero-shot variant effect prediction**, including clinically relevant human variants such as BRCA1 missense and splice variants, without task-specific supervision.
- **Cross-domain generality.** One model handles prokaryotic, eukaryotic, coding, and non-coding sequence, where earlier genomic models were domain-specific.
- **Long-range structure.** Effective retrieval across the full 1M-token context lets it model regulatory context that protein-level models cannot see at all.
- **Generation at genome scale**, including sequences comparable in length to minimal bacterial genomes.
- **Interpretability.** Sparse autoencoders trained on layer 26 of the 7B model recover biologically meaningful features — exon–intron boundaries, protein secondary structure elements, phage regions. This is the most interpretability-legible biological foundation model available, and the SAE work was published alongside the model rather than years later.

## Limitations

- **Worse than dedicated protein language models on protein tasks.** Modeling proteins through their DNA encoding is indirect; codon degeneracy and untranslated context add noise that pLMs never see.
- **Eukaryotic generation lags prokaryotic.** Quality metrics on generated eukaryotic sequence are noticeably weaker, and the paper is candid about this.
- **Generated genomes are evaluated in silico.** No generated genome has been synthesized and shown to function. "Genome design" is currently a claim about sequence plausibility scores.
- **Compute is prohibitive for replication.** Roughly 2,000 H100s for months. Independent verification of the headline results is effectively impossible for academic groups, which is a structural problem for a fully open release.
- **Dual-use is sharpest here.** The training set excludes pathogens of concern. That exclusion is a safety claim, and it is testable — nobody has published an evaluation of whether the exclusion actually removes the capability it is meant to remove, or whether the model reconstructs it from homologous non-excluded sequence.

## When not to use it

- **Protein-only representation learning.** Modeling proteins through their DNA encoding is indirect and worse; use a pLM.
- **On Colab, at any tier.** StripedHyena 2 needs custom CUDA kernels and the 7B does not fit the free or Pro runtimes.
- **Eukaryotic generation where quality matters.** The authors are candid that it lags prokaryotic.
- **As evidence that genome design works.** No generated genome has been synthesized and shown to function.

## Evaluated by

Zero-shot AUROC on ClinVar and BRCA1 saturation genome editing data; Spearman correlation on deep mutational scanning sets; gene essentiality prediction from mutation effects; perplexity and naturalness scores on generated sequences; SAE feature interpretability with expert annotation; the public Evo Designer and mechanistic interpretability visualizer as qualitative tools.

## Ideas

- **Test the exclusion.** Take a held-out excluded pathogen family, measure whether Evo 2's likelihoods on those sequences are distinguishable from a matched non-excluded control. If they are not, training-set exclusion does not do what its proponents claim, and a widely cited biosecurity mitigation needs rethinking. This is the single highest-value experiment on this list and the open weights make it feasible.
- Use the released SAEs as a contestability testbed: can a domain expert, given only feature activations, correctly predict and challenge a variant call? That is the concrete version of the interpretability-as-contestability thesis, on a model where the tooling already exists.
- Quantify the DNA-versus-protein penalty precisely. Same variant effect benchmark, Evo 2 versus ESM-2 versus SaProt. The size of the gap tells you whether genome-scale modeling is a different capability or a worse route to the same one.
- Check whether the 1M-token context is doing work on protein-level tasks or only on regulatory ones. If long context does not help protein prediction, the architectural argument for DNA models in protein applications weakens considerably.
