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
# # 05 Stage 3 — ProtGPT3 production scoring (full-LLR primary, capped)
#
# Runs ONLY after the runtime + cap-stability gate PASSES. Scores the three
# single-sequence ProtGPT3 base checkpoints with the frozen primary convention —
# full-sequence LLR — over all context-safe assays, capped at N_CAP=2000 variants
# per assay (seed 0, same variants across checkpoints), and writes per-assay Spearman
# `rho`. Computes **no** `lp:ld`; the FE interaction and injection MDE are the
# separate analysis step. Scoring math is the offline-validated `src.protgpt3_scoring`
# module, shared with the gate so the convention cannot drift.
#
# The WT-marginal convention was withdrawn after failing its gate (commit 55b0c22);
# it is not used here.

# %%
import os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
import sys; sys.path.insert(0, str(ROOT))
from src import protgpt3_scoring as pg
print(pg.selftest())                                  # synthetic unit test before scoring

(ROOT / "results").mkdir(exist_ok=True)
PROTEINGYM_COMMIT = os.environ.get("PROTEINGYM_COMMIT",
                                   "144fe22b07dfaeec2b366f2346203a9838a55b4c")
CTX_MAX_RESIDUES = 1022
DMS_DIR = Path(os.environ.get("DMS_DIR", "/content/ProteinGym/DMS_ProteinGym_substitutions"))
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
CHUNK = int(os.environ.get("CHUNK", "16"))
SEED, N_CAP = 0, 2000
MIN_VARIANTS = 10
MAX_DROP_FRAC = 0.20
CHECKPOINTS = {
    "AI4PD/ProtGPT3-112M": dict(total=0.109e9, active=0.034e9),
    "AI4PD/ProtGPT3-1.3B": dict(total=1.328e9, active=0.366e9),
    "AI4PD/ProtGPT3-10B":  dict(total=10.000e9, active=2.752e9),
}
_sel = os.environ.get("CHECKPT")
if _sel:
    CHECKPOINTS = {k: v for k, v in CHECKPOINTS.items() if _sel in k}

# %%
REF_URL = (f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/"
           f"{PROTEINGYM_COMMIT}/reference_files/DMS_substitutions.csv")
ref = pd.read_csv(REF_URL).dropna(subset=["target_seq", "MSA_Neff_L", "UniProt_ID"])
ref = ref[ref["seq_len"] <= CTX_MAX_RESIDUES].reset_index(drop=True)
print(f"context-safe assays: {len(ref)}")             # expect 201

def prep(dms_id, wt):
    """20-AA-universe filter, then the seed-0 cap. Returns (mutated_seqs, DMS, n_usable,
    n_dropped, cap_applied)."""
    df = pd.read_csv(DMS_DIR / f"{dms_id}.csv")
    keep = np.array([pg.parse_variant(m, wt) is not None for m in df["mutant"].astype(str)])
    u = df[keep].reset_index(drop=True)
    n_drop = int((~keep).sum())
    cap = len(u) > N_CAP
    if cap:
        u = u.iloc[pg.cap_indices(len(u), N_CAP, SEED)].reset_index(drop=True)
    return u["mutated_sequence"].values, u["DMS_score"].values, int(keep.sum()), n_drop, cap

# %%
rows = []
for repo, meta in CHECKPOINTS.items():
    print(f"\nloading {repo} ...")
    model, tok = pg.load(repo, device=DEVICE)
    for _, r in ref.iterrows():
        f = DMS_DIR / f"{r['DMS_id']}.csv"
        if not f.exists():
            rows.append(dict(assay=r["DMS_id"], repo=repo, status="dms_missing")); continue
        seqs, y, n_usable, n_drop, cap = prep(r["DMS_id"], r["target_seq"])
        drop_frac = n_drop / (n_usable + n_drop) if (n_usable + n_drop) else 1.0
        if n_usable < MIN_VARIANTS:
            rows.append(dict(assay=r["DMS_id"], repo=repo, status="below_variant_floor",
                             n_usable=n_usable, n_dropped=n_drop)); continue
        if drop_frac > MAX_DROP_FRAC:
            rows.append(dict(assay=r["DMS_id"], repo=repo, status="high_drop_frac",
                             n_usable=n_usable, n_dropped=n_drop, drop_frac=round(drop_frac, 4))); continue
        scores = pg.full_llr(model, tok, r["target_seq"], seqs, DEVICE, CHUNK)
        ok = np.isfinite(scores) & np.isfinite(y)
        status = "ok" if np.unique(scores[ok]).size >= 2 else "constant_scores"
        rho = spearmanr(scores[ok], y[ok]).correlation if status == "ok" else np.nan
        rows.append(dict(assay=r["DMS_id"], repo=repo, status=status,
                         params=meta["total"], active_params=meta["active"],
                         MSA_Neff_L=r["MSA_Neff_L"], UniProt_ID=r["UniProt_ID"],
                         seq_len=r["seq_len"], n_scored=len(seqs), n_dropped=n_drop,
                         drop_frac=round(drop_frac, 4), cap_applied=cap, rho_llr=rho))
    del model; torch.cuda.empty_cache() if DEVICE == "cuda" else None

scores = pd.DataFrame(rows)
scores.to_csv(ROOT / "results/05_stage3_protgpt3_scores.csv", index=False)
print(f"\nwrote {len(scores)} assay-checkpoint rows")

# %% [markdown]
# ## Technical sanity summary (descriptive rho is NOT a kill criterion)

# %%
ok = scores[scores.status == "ok"]
print("status counts:\n", scores.status.value_counts().to_string())
if len(ok):
    print("\ndescriptive per-checkpoint mean rho (diagnostic only, does not gate the fit):")
    print(ok.groupby("params")["rho_llr"].agg(["mean", "median", "size"]).to_string())
broken = scores[scores.status.isin(["constant_scores", "dms_missing"])]
if len(broken):
    print(f"\nWARNING: {len(broken)} rows flagged (constant/missing) — inspect before the FE step.")

# %%
from src.provenance import stamp  # noqa: E402
stamp("05_stage3_protgpt3_scoring", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, convention="full_llr_capped",
      n_cap=N_CAP, seed=SEED, context_max_residues=CTX_MAX_RESIDUES,
      n_assays=int(len(ref)), checkpoints=list(CHECKPOINTS), device=str(DEVICE),
      note="production scoring; no lp:ld computed here")
