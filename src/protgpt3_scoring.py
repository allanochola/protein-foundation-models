"""Stage 3 ProtGPT3 scoring — shared by the runtime gate and the production scorer
so the two cannot drift on the convention.

Conventions are frozen in experiments/05-architecture-vs-power/STAGE3_PROTOCOL.md:

- Primary: autoregressive WT-marginal. One forward over the wild-type sequence
  ([<|bos|>] + residues + [<|eos|>]); a substitution a->b at residue position p is
  scored by z[b] - z[a] where z = log_softmax(logits[p-1]) is the model's
  distribution over residue p given the WT prefix. Additive across positions for
  multi-substitution variants. Left-context-only estimand.
- Validation: full-sequence LLR, logP(mutant) - logP(wt), one forward per variant.

The log-probability indexing here mirrors an offline brute-force-validated
reference (see selftest() at the bottom): seq-logprob == sum of per-position
log-softmax, WT-marginal position offset, batched == per-sequence, LLR(wt,wt)==0.
"""
import re
import numpy as np

AA20 = "ACDEFGHIKLMNPQRSTVWY"
_MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def cap_pair(n, n_cap=2000, seed=0):
    """Deterministic (cap2, cap4) for a size-n assay. cap2 is the production cap
    (size min(n_cap, n)); cap4 is the up-to-2x superset used by the gate's
    cap-stability check. cap2 is nested in cap4 by construction, so production scores
    exactly the set the gate validated."""
    rng = np.random.default_rng(seed)
    big = np.sort(rng.choice(n, size=min(2 * n_cap, n), replace=False))
    return big[:n_cap], big


def cap_indices(n, n_cap=2000, seed=0):
    """The production cap: first n_cap of the seed-0 2x draw (all of it if n<=n_cap)."""
    return cap_pair(n, n_cap, seed)[0]


def load(repo, dtype=None, device="cuda"):
    import torch
    from transformers import AutoModelForCausalLM, AutoTokenizer
    dtype = dtype or torch.bfloat16
    tok = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
    model = AutoModelForCausalLM.from_pretrained(
        repo, trust_remote_code=True, dtype=dtype).to(device).eval()
    return model, tok


def aa_ids(tok):
    return {a: tok(a, add_special_tokens=False)["input_ids"][0] for a in AA20}


def encode_wt(tok, seq):
    # tokenizer omits <|bos|> even with add_special_tokens=True — add both ends here
    return [tok.bos_token_id] + tok(seq, add_special_tokens=False)["input_ids"] + [tok.eos_token_id]


def parse_variant(mstr, wt_seq):
    """Return [(pos0, from_aa, to_aa), ...] under the 20-AA universe, or None if the
    variant must be dropped: non-standard AA, out-of-range position, or a stated
    wild-type residue that disagrees with the reference."""
    subs = []
    for sub in str(mstr).split(":"):
        m = _MUT.match(sub.strip())
        if not m:
            return None
        frm, pos, to = m.group(1), int(m.group(2)), m.group(3)
        if frm not in AA20 or to not in AA20:
            return None
        if pos < 1 or pos > len(wt_seq):
            return None
        if wt_seq[pos - 1] != frm:
            return None
        subs.append((pos - 1, frm, to))
    return subs


def _marginal_from_logp(logp, subs, aid):
    """Pure math: sum over substitutions of z[to]-z[from] at the mutated position.
    logp[p0] is log_softmax of the logits that predict residue at 1-indexed position
    p0+1 (input index p0+1, predicted by logits[p0])."""
    s = 0.0
    for pos0, frm, to in subs:
        row = logp[pos0]
        s += row[aid[to]] - row[aid[frm]]
    return s


def wt_marginal(model, tok, wt_seq, mutant_strings, device="cuda"):
    """One WT forward; per-variant additive WT-marginal score. Returns
    (scores[np.nan for dropped], dropped_indices)."""
    import torch
    ids = torch.tensor([encode_wt(tok, wt_seq)], device=device)
    with torch.no_grad():
        logits = model(input_ids=ids, attention_mask=torch.ones_like(ids)).logits[0].float()
    logp = torch.log_softmax(logits, dim=-1).cpu().numpy()      # [N, V]
    aid = aa_ids(tok)
    scores = np.full(len(mutant_strings), np.nan)
    dropped = []
    for k, mstr in enumerate(mutant_strings):
        subs = parse_variant(mstr, wt_seq)
        if subs is None:
            dropped.append(k)
            continue
        scores[k] = _marginal_from_logp(logp, subs, aid)
    return scores, dropped


def _seq_logprobs_batch(model, tok, seqs, device="cuda"):
    import torch
    pad = tok.pad_token_id
    enc = [encode_wt(tok, s) for s in seqs]
    mx = max(len(e) for e in enc)
    ids = torch.tensor([e + [pad] * (mx - len(e)) for e in enc], device=device)
    attn = (ids != pad).long()
    with torch.no_grad():
        logits = model(input_ids=ids, attention_mask=attn).logits[:, :-1, :].float()
    logp = torch.log_softmax(logits, dim=-1)
    tgt = ids[:, 1:]
    lp = logp.gather(-1, tgt.unsqueeze(-1)).squeeze(-1) * attn[:, 1:].float()
    return lp.sum(-1).cpu().numpy()


def full_llr(model, tok, wt_seq, mutated_seqs, device="cuda", chunk=32):
    """Full-sequence LLR = logP(mut) - logP(wt), one forward per variant (batched).
    WT offset is Spearman-irrelevant; kept for interpretability."""
    wt_lp = _seq_logprobs_batch(model, tok, [wt_seq], device)[0]
    out = np.empty(len(mutated_seqs))
    for i in range(0, len(mutated_seqs), chunk):
        out[i:i + chunk] = _seq_logprobs_batch(model, tok, list(mutated_seqs[i:i + chunk]), device)
    return out - wt_lp


def selftest():
    """Offline validation of the pure math with numpy (no model, no torch)."""
    rng = np.random.default_rng(0)

    def log_softmax(x):
        x = x - x.max(-1, keepdims=True)
        return x - np.log(np.exp(x).sum(-1, keepdims=True))

    # parse_variant
    assert parse_variant("A2C", "MACDE") == [(1, "A", "C")]
    assert parse_variant("A2C:D4E", "MACDE") == [(1, "A", "C"), (3, "D", "E")]
    assert parse_variant("A2C", "MMCDE") is None      # WT mismatch (pos2 is M, not A)
    assert parse_variant("A9C", "MACDE") is None      # out of range
    assert parse_variant("A2X", "MACDE") is None      # non-standard AA
    assert parse_variant("bad", "MACDE") is None

    # WT-marginal sum matches brute force over log-softmax rows
    V = 25
    aid = {a: i for i, a in enumerate(AA20)}
    logp = log_softmax(rng.normal(size=(10, V)))
    subs = [(2, "A", "C"), (5, "D", "E")]
    got = _marginal_from_logp(logp, subs, aid)
    ref = sum(logp[p0][aid[to]] - logp[p0][aid[frm]] for p0, frm, to in subs)
    assert abs(got - ref) < 1e-12

    # WT-marginal difference is invariant to the softmax normalizer (20-vs-full vocab)
    raw = rng.normal(size=(10, V))
    d_norm = _marginal_from_logp(log_softmax(raw), subs, aid)
    d_raw = sum(raw[p0][aid[to]] - raw[p0][aid[frm]] for p0, frm, to in subs)
    assert abs(d_norm - d_raw) < 1e-10

    # cap: cap2 nested in cap4, deterministic, correct sizes
    c2, c4 = cap_pair(536962, n_cap=2000, seed=0)
    assert len(c2) == 2000 and len(c4) == 4000 and set(c2).issubset(set(c4))
    assert np.array_equal(cap_indices(536962), c2)
    c2b, c4b = cap_pair(536962, n_cap=2000, seed=0)
    assert np.array_equal(c2, c2b) and np.array_equal(c4, c4b)          # reproducible
    small2, small4 = cap_pair(1500, n_cap=2000)
    assert len(small2) == 1500 and len(small4) == 1500                  # n<=n_cap -> all
    return "protgpt3_scoring selftest: OK"


if __name__ == "__main__":
    print(selftest())
