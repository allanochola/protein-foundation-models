# ProGen

Madani et al., arXiv (2020) and *Nature Biotechnology* 41:1099–1106 (2023). "Large language models generate functional protein sequences across diverse families." See also ProGen2 (Nijkamp et al., *Cell Systems*, 2023), 151M–6.4B parameters.

## Trained on

Roughly 280 million sequences drawn from UniParc, UniProtKB, Pfam, and NCBI taxonomy, spanning more than 19,000 protein families. 1.2B parameters, CTRL-style autoregressive decoder. The distinguishing feature is conditioning: each sequence is prefixed with control tags encoding taxonomy, Pfam family, molecular function, cellular component, and biological process. Generation is steered by choosing tags.

For the lysozyme work the base model is fine-tuned on five curated lysozyme families.

## Does well

- **Controllable family-conditional generation.** Set the tag, get sequences from that family. This is the first protein model where the control interface is explicit rather than emergent from prompting.
- **Functional wet-lab validation.** Artificial lysozymes showed catalytic efficiencies comparable to natural enzymes at sequence identities as low as 31.4% to any natural protein. A designed variant was crystallized at 2.5 Å (PDB 7RGR) and matched the predicted fold.
- **Transfers across enzyme classes.** The same recipe worked for chorismate mutase and malate dehydrogenase.

## Limitations

- **Needs deep families.** Controllable generation degrades where homologous samples are scarce, which is where design would be most valuable. The model interpolates within families rather than producing new function.
- **Function is recapitulated, not invented.** Every validated design does the job its family already does. No result shows novel catalytic activity.
- **The wet-lab set is small and assay-friendly.** Lysozyme has a cheap fluorescence readout. Three enzyme families with easy assays is not evidence of general design capability, and the paper does not report the full synthesis-to-success funnel.
- **No structure modality and no structural constraint at generation time.** Filtering relies on downstream predictors.
- **Control tags are a misuse surface.** The interface that makes ProGen useful — specify function, receive sequences — is the same interface that makes screening necessary. The paper does not discuss it.

## When not to use it

- **Mutation effect prediction.** It is an autoregressive generator. ProGen2 is scored on ProteinGym, but ESM-2 and SaProt outperform it per unit of compute.
- **Families with few homologs** — exactly where design would be most valuable and where controllable generation degrades.
- **When you need novel function.** Every validated design recapitulates its family's existing activity.
- **Without a screening layer.** The control-tag interface that makes it useful is also the misuse surface, and the paper provides no filter.

## Evaluated by

Perplexity conditioned on tag; predicted secondary structure accuracy against family consensus; Rosetta and trRosetta conformational energy of generated sequences; sequence identity distribution to the nearest natural protein as a novelty measure; and, for lysozymes, measured kcat/KM against hen egg-white lysozyme plus X-ray crystallography.

## Ideas

- ProGen is the best available testbed for conditional-generation screening research: the control tags give a clean, labeled dimension along which to test whether a filter blocks hazardous conditioning without blocking benign conditioning. Nobody has published that experiment.
- Measure the identity-versus-function frontier systematically. At what sequence identity to the nearest natural protein does measured activity collapse? The lysozyme data hints at 31%; a proper curve would tell you whether generative models can leave the natural distribution at all.
- Compare ProGen's tag conditioning against ESM-3's function-track conditioning on matched prompts. Explicit tags versus learned multimodal tokens is a control-mechanism comparison with direct oversight implications — an explicit tag is auditable, a function token is not.
