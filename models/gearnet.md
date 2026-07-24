# GearNet

Zhang et al., ICLR 2023. "Protein Representation Learning by Geometric Structure Pretraining." Included here as the control condition: structure supervision without language modeling.

## Trained on

805,000 predicted structures from the AlphaFold DB — 365K proteome-wide predictions plus Swiss-Prot. Roughly 20M parameters, three orders of magnitude smaller than the language models in this set and with about 1/80th the training data.

Architecture is a relational graph neural network over residue graphs. Edges are typed — sequential, radius, K-nearest-neighbor — and message passing is edge-type-specific. GearNet-Edge adds message passing between edges, which captures angular relationships. Everything is built from distances and angles, so the encoder is E(3)-invariant by construction rather than by augmentation.

Five self-supervised objectives are compared: multiview contrast, residue type prediction, distance prediction, angle prediction, dihedral prediction. Multiview contrast wins.

## Does well

- **Function prediction at a fraction of the data.** Matches or beats sequence-based models on EC number and GO term prediction while pretraining on 805K structures rather than tens of millions of sequences.
- **Fold classification**, including the hard superfamily and fold-level splits where sequence identity gives no signal.
- **Reaction classification.**
- **Clean ablations.** The paper isolates the contribution of relational convolution, edge message passing, and each pretraining objective — rarer than it should be.

## Limitations

- **Narrow task coverage.** Function and fold. No variant effect prediction, no structure prediction, no generation, no zero-shot anything.
- **Requires structure at inference,** with all the AF2-dependence problems that implies.
- **Residue-level coarse-graining** discards side-chain atoms, which limits binding-site and catalytic-mechanism tasks.
- **Small scale, and untested at large scale.** The authors say explicitly that they used only 805K of the 100M+ available structures. Nobody has run the scaling study, so it remains unknown whether the data-efficiency advantage survives.
- **Does not generalize to sequence-only settings** — the entire input assumption differs from the pLM family.

## When not to use it

- **Without structures.** There is no sequence-only fallback.
- **Variant effect prediction or generation.** It does function and fold classification, nothing else.
- **As a general-purpose embedding model.** Task coverage is narrow by design.
- **As evidence about scale.** It was trained on 805K of 100M+ available structures and the scaling study has never been run.

## Evaluated by

Fmax on EC and GO prediction (biological process, molecular function, cellular component); accuracy on fold classification under fold, superfamily, and family splits; reaction classification accuracy; AUPR on EC and GO. Baselines include CNN, ResNet, LSTM, transformer, GCN, GAT, GVP, 3DCNN-MQA, and IEConv.

## Ideas

- **This is the load-bearing comparison in the reading set.** GearNet matching sequence models on function prediction with 1% of the data is evidence that scale is not what produces function knowledge — structure is. Every scaling claim in the pLM literature should be checked against this baseline, and mostly is not.
- Run the scaling study the authors declined to run: GearNet on 100M AlphaFold structures. Either the curve flattens, in which case structure-based learning saturates early and the pLM scaling story is safe, or it does not, in which case the field has been spending compute in the wrong place. Either result is publishable.
- Use GearNet as the confound control for structure-aware pLMs. SaProt and ESM-3 both claim gains from adding structure. GearNet tells you what pure structure supervision buys. The difference is the interaction term, and nobody has measured it.
- Compare the geometric self-supervised objectives against masked language modeling on a shared evaluation. The five-objective ablation exists on the structure side and has no counterpart on the sequence side.
