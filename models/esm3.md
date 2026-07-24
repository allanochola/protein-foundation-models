# ESM-3

Hayes et al., EvolutionaryScale (2024 preprint; *Science*, 2025). "Simulating 500 million years of evolution with a language model."

## Trained on

2.78 billion proteins and 771 billion unique tokens, at 1.07 × 10²⁴ FLOPs for the largest run. Three scales: 1.4B, 7B, and 98B parameters. Sequence data from UniRef, MGnify, and JGI metagenomes; structure from PDB plus predicted structures; function annotations from InterPro keywords.

The architectural move is tokenization of all three modalities. Local atomic structure around each residue is encoded into discrete structure tokens by a learned VQ encoder, function annotations become keyword tokens, and all tracks are fused in one latent space. Geometric attention in the first block allows conditioning on raw coordinates. Training is masked-token prediction across tracks with a variable noise schedule, so any track can prompt any other. Because experimental structure and function labels are scarce, the training set is augmented with hundreds of millions of synthetic annotations.

## Does well

- **Any-to-any conditional generation.** Prompt with partial sequence, partial structure, function keywords, or any combination. This is the capability no earlier pLM had.
- **Single-sequence structure prediction.** ESM3-98B reaches 0.895 mean LDDT on CAMEO against 0.865 for ESMFold.
- **Scaling across all three tracks.** Test loss improves from 1.4B to 98B on sequence, structure, and function, with the largest gain on sequence.
- **Alignment-style post-training works.** Preference optimization raised the success rate on hard design prompts from 26.8% to 65.5%.
- **Wet-lab novelty.** esmGFP fluoresces at levels comparable to natural GFPs at 58% identity to its nearest known relative.

## Limitations

- **The headline claim is rhetoric, not a measurement.** "500 million years of evolution" is an inference from sequence distance to evolutionary time, calibrated against GFP family divergence rates. It is a single-protein extrapolation presented as a general property of the model.
- **esmGFP is one success.** The paper does not report the design funnel — how many candidates were synthesized, how many failed, what the hit rate was against a matched baseline. Without the denominator the result is not a capability estimate.
- **Structure tokenization is lossy** and the reconstruction error is not propagated into downstream confidence.
- **Synthetic training data inherits AlphaFold2's error distribution.** Regions where AF2 is systematically wrong become regions where ESM-3 is confidently wrong.
- **The 98B weights are closed.** Only the 1.4B open model is independently auditable, and it is not the model the headline results come from. Every safety claim about the frontier model rests on internal evaluation.
- **No published screening evaluation.** A model that generates functional sequences conditioned on function keywords has an obvious misuse surface, and the release materials describe policy commitments rather than measured refusal or filtering performance.

## When not to use it

- **Anything that must ship under a permissive license.** Cambrian Non-Commercial propagates to derivative models and methods. See [third-party licenses](../docs/third-party-licenses.md).
- **As a stand-in for the headline results.** The open 1.4B is not the 98B the paper reports. Conclusions from one do not transfer to the other.
- **Routine variant effect prediction.** ESM-2 650M or SaProt does it cheaper, and neither carries the license constraint.
- **Any claim about design capability without a denominator.** esmGFP is one success with an unreported funnel.

## Evaluated by

Negative log-likelihood per track as a function of training FLOPs; prompt faithfulness via backbone cRMSD, three-state secondary structure accuracy, SASA Spearman correlation, and keyword recovery; pTM and pLDDT for generation confidence; CAMEO LDDT for structure prediction; sequence identity and TM-score to nearest training neighbor as a novelty proxy; wet-lab fluorescence for esmGFP.

## Ideas

- Novelty is measured against the *training set*. Re-measure against the full sequence universe including post-cutoff depositions. If generations are novel relative to training but not relative to nature, the model is recovering rather than inventing.
- Probe the function track directly: do function keyword tokens have linearly decodable representations, and do they generalize to held-out InterPro families? This is a clean interpretability target and it maps onto contestability — if the function track is not readable, no one can contest a generation's claimed function.
- Design a capability evaluation with a reported denominator: fixed prompt set, fixed synthesis budget, hit rate against a ProGen and an ESM-2 baseline. The absence of this format across the field is why capability claims are not comparable.
- Compare the 1.4B open model to the 98B on any evaluation available for both. The gap sets a floor on how much open-model auditing can tell you about the closed frontier model — a governance-relevant number nobody has published.
