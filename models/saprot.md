# SaProt

Su et al., ICLR 2024 (spotlight). "SaProt: Protein Language Modeling with Structure-aware Vocabulary."

## Trained on

Roughly 40 million protein sequences paired with AlphaFold DB predicted structures. The paper model is 650M parameters, masked language modeling, trained on 64 A100s for about three months — deliberately matched to ESM-2 650M so the comparison isolates the vocabulary change. Released at 35M, 650M, and 1.3B; the 1.3B version (December 2024) comes in AF2-only and AF2+OMG+NCBI variants, and the authors report it beats 650M specifically on sequence-only tasks. A 35M residue-sequence-only model ships alongside as an explicit ablation control.

The idea is small and effective. Foldseek encodes local 3D structure into 20 discrete 3Di states, one per residue. Cross the 20 amino acids with the 20 structure states plus a mask and you get a 441-token structure-aware vocabulary. Each residue becomes a single token carrying both identity and local fold. No architectural change is needed — any residue-based pLM can consume the new alphabet.

## Does well

- **Beats ESM-2 650M across ten downstream tasks** at matched parameter count: clinical variant prediction, fitness landscape prediction, protein–protein interaction, thermostability, metal-ion binding, EC and GO prediction.
- **Cheap structure integration.** No geometric attention, no diffusion, no coordinate handling. The structure signal arrives through the tokenizer.
- **Drop-in compatibility.** Existing pLM training and fine-tuning code works unchanged.

## Limitations

- **Entirely dependent on Foldseek.** 3Di was designed for fast structure search, not for representation learning. Whatever Foldseek discards, SaProt cannot recover. The paper acknowledges this as its first limitation.
- **Best performance requires a structure at inference.** For proteins without experimental or confidently predicted structures — orphans, intrinsically disordered regions, designed sequences — the structure track is unreliable or absent, and those are precisely the cases where a model would add most value. The 1.3B release narrows this: it is reported stronger than 650M on sequence-only tasks, so the structure dependence is a scale-dependent weakness rather than an architectural one. Whether the sequence-only 1.3B still beats a matched ESM-2 is the number to look for, and it is not in the ICLR paper.
- **Training on AF2 predictions creates distribution shift.** The paper's own loss curves separate: validation loss on AF2 structures diverges from loss on real PDB structures. The model is partly learning AlphaFold's biases.
- **20 structure states is coarse.** Local geometry is compressed to a 20-way categorical, losing side-chain orientation and fine backbone detail.
- **No generation.** Encoder only.

## When not to use it

- **Proteins with no structure and no confident prediction.** Orphans, disordered regions, designed sequences — the advantage over ESM-2 comes from the structure track, and stratifying by pLDDT is the first check anyone should run.
- **De novo and designed proteins.** AF2 is confident there but off-distribution, which is the failure mode most likely to bite a generative filter.
- **The 650M in sequence-only mode.** Use 1.3B or ESM-2 instead.
- **Generation of any kind.** Encoder only.

## Evaluated by

Spearman correlation on ProteinGym-style deep mutational scanning; AUC on ClinVar clinical variant classification; Fmax on EC number and GO term prediction; accuracy on thermostability, metal-ion binding, and PPI benchmarks. Baselines are ESM-2, ESM-1b, ProtBert, and MIF-ST, with size held constant where possible.

## Ideas

- SaProt is the cleanest natural experiment in this reading set: identical size, identical objective, one variable changed. Use it to ask what fraction of any pLM's downstream performance is attributable to structure information versus scale. Nobody has run that decomposition properly.
- Stratify every SaProt-versus-ESM-2 win by AlphaFold pLDDT of the input structure. If the advantage disappears below pLDDT 70, the reported gains come from proteins that were already well-characterized — which is a confound, not a capability.
- Swap the tokenizer. Replace 3Di with a learned structure vocabulary at varying codebook sizes and measure the performance-versus-vocabulary-size curve. The 20-state choice is inherited from a search tool, not optimized for this purpose.
- Test on designed and de novo proteins, where AF2 confidence is high but the sequences are off-distribution. This is the failure mode most likely to matter for evaluating generative design pipelines that use SaProt as a filter.
