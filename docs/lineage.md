# Lineage

Who inherited what. Solid arrows are direct architectural or weight descent;
dashed arrows are dependencies on a component rather than a predecessor.

```
              BERT / T5 (NLP)
                    │
        ┌───────────┴────────────┐
        │                        │
   ProtBERT ── ProtT5       ESM-1b (2021, Meta)
   (2021, same paper)            │
                                 ▼
                            ESM-2 (2022, Meta)
                             │        │
                             │        └──▶ ESMFold (folding head on frozen trunk)
                             │
             ┌───────────────┴─────────────────┐
             │ backbone reused                  │ team departs Meta
             ▼                                  ▼
        SaProt (2024, Westlake)          ESM-3 (2024, EvolutionaryScale)
             ▲                                  │
             ┆ structure alphabet               ├──▶ ESM C (2024) — representation branch
             ┆                                  │
        Foldseek 3Di (2023)              multimodal masked generation


   CTRL (NLP)                    Hyena / StripedHyena
       │                                │
       ▼                                ▼
   ProGen (2020, Salesforce)       Evo (2024, Arc)
       │                                │
       ▼                                ▼
   ProGen2 (2023)                  Evo 2 (2025, Arc + NVIDIA)


   Relational GNNs ──▶ GearNet (2023, MILA)     [no shared ancestry with the above]
```

## What the graph is for

**SaProt inherits the ESM-2 backbone.** This is the most useful edge here.
Same architecture, same objective, matched parameter count, one variable
changed — the vocabulary. That makes the SaProt/ESM-2 comparison a controlled
experiment rather than a benchmark result, and it is why SaProt is the right
place to ask how much of a pLM's performance comes from structure versus
scale.

**ESM-2 to ESM-3 is not a version bump.** Different organization, different
architecture class, different licensing regime. Treating them as a single
scaling series — as several surveys do — mixes an MLM encoder with a
multimodal masked-generation model and reads the difference as scale.

**ProtBERT and ProtT5 are siblings, not successors.** Both come from the
ProtTrans paper. The paper is an architecture ablation, which is what makes it
worth revisiting.

**GearNet has no ancestry in common with anything else here.** That is the
point of including it. It is the only model in the set whose performance
cannot be attributed to language-model pretraining, which makes it the control
condition for every claim that scale produces function knowledge.

**Two generative lines, both borrowed from NLP.** ProGen from CTRL's control
codes, Evo from Hyena's subquadratic sequence mixing. Neither innovation
originated in biology, which is worth remembering when a paper frames an
architectural choice as biologically motivated.
