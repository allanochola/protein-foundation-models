# Third-party model licenses

This repository's MIT license covers our code and analysis only. Model
weights and upstream code carry their own terms, and some of them propagate
into derivative work.

Two entries below are verified against primary sources. The rest are marked
**unverified** — check before you rely on them. An unverified license in a
table is worse than no table, because it looks like it was checked.

| Model | Weights license | Propagates? | Status |
|---|---|---|---|
| ESM-3 open (1.4B) | Cambrian Non-Commercial | **Yes** | Verified |
| ESM C 300M + `esm` repo code | Cambrian Open (permissive) | No | Verified |
| ESM C 600M | Cambrian Non-Commercial | **Yes** | Verified |
| ESM C 6B, ESM-3 7B/98B | API only, no weights | n/a | Verified |
| ESM-2 (`fair-esm`) | MIT (believed) | No | Unverified |
| ProtT5 / ProtTrans | Check repo and HF model card separately | ? | Unverified |
| ProGen / ProGen2 | BSD-3-Clause (claimed) | No | Unverified |
| Evo 2 | Apache 2.0 (believed) | No | Unverified |
| SaProt | Check repo | ? | Unverified |
| GearNet | Check repo | ? | Unverified |
| ProteinGym data and scores | Check repo | ? | Unverified |

## The constraint that actually bites

The Cambrian Non-Commercial License requires that derivative models and
methods be released under the same terms, and that derivative work carry
"Built with ESM" attribution. Fine-tuning ESM-3 open, or building a method on
top of it, pulls non-commercial terms into that work — which conflicts with
publishing it under MIT here.

Practical consequence for this project: **keep ESM-3 work in a separate
directory with its own license notice, or restrict ESM-3 to scoring and
analysis that produces no derivative model.** The experiment 01 line of work
is safe on this front because it uses published benchmark scores rather than
weights.

Evo 2 is the opposite case. Fully open weights, code, and training data are
what make the training-set exclusion test in [evo2.md](../models/evo2.md)
feasible at all. That contrast — the model with the sharpest dual-use profile
being the most auditable, and the frontier protein model being the least — is
worth stating in any writeup.

## Before adding a model to this repo

1. Read the actual LICENSE file in the upstream repo, not the README summary.
2. Check the HuggingFace model card separately — weights and code frequently
   differ.
3. Record it in the table above with a link and the date checked.
4. If it propagates, decide where the work lives before writing any code.
