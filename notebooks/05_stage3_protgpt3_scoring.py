# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3
#     name: python3
# ---

# %% [markdown]
# # 05 Stage 3 — ProtGPT3 scoring (produces per-assay rho, no coefficient fit)
#
# Executes the scoring stage of `experiments/05-architecture-vs-power/STAGE3_PROTOCOL.md`.
# GPU. This notebook scores the three single-sequence ProtGPT3 base checkpoints on
# the context-safe 201-assay subset and writes per-assay Spearman rho for two frozen
# conventions (full-sequence LLR primary, WT-marginal robustness). It STOPS at the
# scores CSV and the kill-criteria sanity block. The fixed-effects interaction (A)
# and injection MDE (B) are a separate notebook, run only after scoring passes
# sanity — no `lp:ld` coefficient is computed here.
#
# The log-probability indexing below is a line-by-line mirror of a brute-force
# reference validated offline (full-LLR, WT-marginal, batched==per-sequence,
# LLR(wt,wt)=0).

# %%
import os, re, json, gc
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from transformers import AutoModelForCausalLM, AutoTokenizer

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "results").mkdir(exist_ok=True)

PROTEINGYM_COMMIT = os.environ.get("PROTEINGYM_COMMIT",
                                   "144fe22b07dfaeec2b366f2346203a9838a55b4c")
CTX_MAX_RESIDUES = 1022                 # + <|bos|> + <|eos|> = 1024 <= 1025 context
CHECKPOINTS = {                          # total params from audited checkpoint metadata
    "AI4PD/ProtGPT3-112M": dict(total=0.109e9, active=0.034e9),
    "AI4PD/ProtGPT3-1.3B": dict(total=1.328e9, active=0.366e9),
    "AI4PD/ProtGPT3-10B":  dict(total=10.000e9, active=2.752e9),
}
CHUNK = int(os.environ.get("CHUNK", "64"))          # mutants per forward; lower for the 10B
_sel = os.environ.get("CHECKPT")                     # e.g. CHECKPT=112M to validate on one checkpoint
if _sel:
    CHECKPOINTS = {k: v for k, v in CHECKPOINTS.items() if _sel in k}
    assert CHECKPOINTS, f"CHECKPT={_sel} matched no checkpoint"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DTYPE = torch.bfloat16

# DMS substitution assay CSVs (mutant, mutated_sequence, DMS_score), acquired exactly
# as in 01a/01c at the pinned commit. Set to wherever they live in this runtime.
DMS_DIR = Path(os.environ.get("DMS_DIR", "/content/ProteinGym/DMS_ProteinGym_substitutions"))

# %% [markdown]
# ## Assay metadata and context-safe subset
#
# Reference file gives per-assay `target_seq` (WT), `seq_len`, `MSA_Neff_L`,
# `UniProt_ID`. Keep assays with `seq_len <= 1022`. This is the Experiment 04
# context-safe subset; its matched ESM-2 reference (`lp:ld = -0.0180`) is the
# comparison, carried into Stage 3 analysis, not recomputed here.

# %%
REF_URL = (f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/"
           f"{PROTEINGYM_COMMIT}/reference_files/DMS_substitutions.csv")
ref = pd.read_csv(REF_URL)
ref = ref.dropna(subset=["target_seq", "MSA_Neff_L", "UniProt_ID"])
ref = ref[ref["seq_len"] <= CTX_MAX_RESIDUES].reset_index(drop=True)
print(f"context-safe assays (seq_len<={CTX_MAX_RESIDUES}): {len(ref)}")   # expect 201
assert 195 <= len(ref) <= 205, "context-safe subset size off; check reference/commit"

def assay_frame(dms_id):
    df = pd.read_csv(DMS_DIR / f"{dms_id}.csv")
    return df[["mutant", "mutated_sequence", "DMS_score"]].dropna()

# %% [markdown]
# ## Scorer (torch mirror of the offline-validated numpy reference)

# %%
def encode(tok, seq):
    ids = tok(seq, add_special_tokens=False)["input_ids"]
    return [tok.bos_token_id] + ids + [tok.eos_token_id]   # tokenizer does NOT add bos itself

@torch.no_grad()
def seq_logprobs(model, batch_ids, pad_id):
    """Sum log p(token_t | token_<t) over predicted positions 1..N-1 (bos is context
    only; eos included). Mirrors validated seq_logprob with a padding mask."""
    ids = torch.tensor(batch_ids, device=DEVICE)
    attn = (ids != pad_id).long()
    logits = model(input_ids=ids, attention_mask=attn).logits[:, :-1, :].float()
    logp = torch.log_softmax(logits, dim=-1)
    tgt = ids[:, 1:]
    tok_logp = logp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1)
    mask = attn[:, 1:].float()
    return (tok_logp * mask).sum(-1).cpu().numpy()

def pad_batch(seqs, pad_id):
    m = max(len(s) for s in seqs)
    return [s + [pad_id] * (m - len(s)) for s in seqs]

def score_llr(model, tok, wt_seq, mutated_seqs):
    """Full-sequence LLR = logP(mut) - logP(wt). Substitutions preserve length, so a
    per-assay batch is uniform-length; chunk for memory."""
    pad = tok.pad_token_id
    wt_lp = seq_logprobs(model, [encode(tok, wt_seq)], pad)[0]
    out = np.empty(len(mutated_seqs))
    enc = [encode(tok, s) for s in mutated_seqs]
    for i in range(0, len(enc), CHUNK):
        chunk = enc[i:i + CHUNK]
        out[i:i + len(chunk)] = seq_logprobs(model, pad_batch(chunk, pad), pad)
    return out - wt_lp

_MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")
@torch.no_grad()
def score_marginal(model, tok, wt_seq, mutant_strings):
    """WT-marginal: one WT forward; for each substitution at residue position p,
    logp(mut|prefix) - logp(wt|prefix) from logits[p-1] (input index p = residue p,
    since index 0 is bos). Multi-substitution = additive sum. Mirrors validated T4."""
    pad = tok.pad_token_id
    ids = torch.tensor([encode(tok, wt_seq)], device=DEVICE)
    logits = model(input_ids=ids, attention_mask=torch.ones_like(ids)).logits[0].float()
    logp = torch.log_softmax(logits, dim=-1).cpu().numpy()      # [N, V]
    aa_id = {a: tok(a, add_special_tokens=False)["input_ids"][0]
             for a in "ACDEFGHIKLMNPQRSTVWY"}
    scores = np.zeros(len(mutant_strings))
    for k, mstr in enumerate(mutant_strings):
        s = 0.0
        for sub in str(mstr).split(":"):
            m = _MUT.match(sub.strip())
            if not m:
                s = np.nan; break
            frm, pos, to = m.group(1), int(m.group(2)), m.group(3)
            row = logp[pos - 1]                                  # predicts residue at input index pos
            s += row[aa_id[to]] - row[aa_id[frm]]
        scores[k] = s
    return scores

# %% [markdown]
# ## Score all checkpoints

# %%
rows = []
for repo, meta in CHECKPOINTS.items():
    print(f"\nloading {repo} ...")
    tok = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        repo, trust_remote_code=True, torch_dtype=DTYPE).to(DEVICE).eval()
    for _, r in ref.iterrows():
        try:
            df = assay_frame(r["DMS_id"])
        except FileNotFoundError:
            print(f"  MISSING dms file: {r['DMS_id']} — recorded not scored"); continue
        wt = r["target_seq"]
        llr = score_llr(model, tok, wt, df["mutated_sequence"].tolist())
        marg = score_marginal(model, tok, wt, df["mutant"].tolist())
        y = df["DMS_score"].values
        rho_llr = spearmanr(llr, y).correlation
        rho_marg = spearmanr(marg, y).correlation
        rows.append(dict(assay=r["DMS_id"], repo=repo,
                         params=meta["total"], active_params=meta["active"],
                         MSA_Neff_L=r["MSA_Neff_L"], UniProt_ID=r["UniProt_ID"],
                         seq_len=r["seq_len"], n_variants=len(df),
                         rho_llr=rho_llr, rho_marg=rho_marg))
    del model; gc.collect(); torch.cuda.empty_cache()

scores = pd.DataFrame(rows)
scores.to_csv(ROOT / "results/05_stage3_protgpt3_scores.csv", index=False)
print(f"\nwrote {len(scores)} assay-checkpoint rows")

# %% [markdown]
# ## Kill-criteria sanity (before any coefficient fit)
#
# A broken scorer manufactures or hides interactions. Require: per-assay rho in a
# sane range on aggregate, non-degenerate spread, and larger ProtGPT3 not worse on
# average. If these fail, STOP and fix scoring — do not proceed to the FE fit.

# %%
print("mean per-assay rho (full-LLR) by checkpoint:")
agg = (scores.groupby("params")
       .agg(mean_llr=("rho_llr", "mean"), med_llr=("rho_llr", "median"),
            mean_marg=("rho_marg", "mean"), n=("rho_llr", "size"))
       .reset_index().sort_values("params"))
print(agg.to_string(index=False))

ok = True
if not (0.20 <= agg["mean_llr"].min() and agg["mean_llr"].max() <= 0.60):
    print("WARN: mean rho outside the sane ~0.2-0.6 band for a family of this scale"); ok = False
if agg["mean_llr"].iloc[-1] < agg["mean_llr"].iloc[0] - 0.05:
    print("WARN: largest checkpoint scores materially worse than smallest"); ok = False
if scores["rho_llr"].std() < 0.05:
    print("WARN: degenerate rho spread across assays"); ok = False
print("\nSANITY:", "PASS — proceed to the Stage 3 FE notebook" if ok
      else "FAIL — fix scoring before any lp:ld fit")

# %%
from src.provenance import stamp  # noqa: E402
import sys; sys.path.insert(0, str(ROOT))
stamp("05_stage3_protgpt3_scoring", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, checkpoints=list(CHECKPOINTS),
      context_max_residues=CTX_MAX_RESIDUES, n_assays=int(len(ref)),
      conventions=["full_sequence_llr", "wt_marginal"], device=str(DEVICE))
