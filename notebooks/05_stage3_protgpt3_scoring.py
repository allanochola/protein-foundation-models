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
# # 05 Stage 3 — ProtGPT3 production scoring (WT-marginal primary)
#
# Runs ONLY after the runtime gate PASSES. Scores the three single-sequence ProtGPT3
# base checkpoints with the frozen primary convention — autoregressive WT-marginal —
# over all context-safe assays, and writes per-assay Spearman `rho`. Full-LLR is not
# run here (it is the gate's validation convention). This notebook computes **no**
# `lp:ld`; the FE interaction and injection MDE are the separate analysis step.
#
# Scoring math is the offline-validated `src.protgpt3_scoring` module, shared with
# the gate so the convention cannot drift.

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
MIN_VARIANTS = 10                                     # usable-variant floor per assay
MAX_DROP_FRAC = 0.20                                  # per-assay parse/universe drop tolerance
CHECKPOINTS = {
    "AI4PD/ProtGPT3-112M": dict(total=0.109e9, active=0.034e9),
    "AI4PD/ProtGPT3-1.3B": dict(total=1.328e9, active=0.366e9),
    "AI4PD/ProtGPT3-10B":  dict(total=10.000e9, active=2.752e9),
}
_sel = os.environ.get("CHECKPT")
if _sel:
    CHECKPOINTS = {k: v for k, v in CHECKPOINTS.items() if _sel in k}

# %% [markdown]
# ## Context-safe assays

# %%
REF_URL = (f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/"
           f"{PROTEINGYM_COMMIT}/reference_files/DMS_substitutions.csv")
ref = pd.read_csv(REF_URL).dropna(subset=["target_seq", "MSA_Neff_L", "UniProt_ID"])
ref = ref[ref["seq_len"] <= CTX_MAX_RESIDUES].reset_index(drop=True)
print(f"context-safe assays: {len(ref)}")             # expect 201

# %% [markdown]
# ## Score (WT-marginal) with technical sanity per assay

# %%
rows = []
for repo, meta in CHECKPOINTS.items():
    print(f"\nloading {repo} ...")
    model, tok = pg.load(repo, device=DEVICE)
    for _, r in ref.iterrows():
        f = DMS_DIR / f"{r['DMS_id']}.csv"
        if not f.exists():
            rows.append(dict(assay=r["DMS_id"], repo=repo, status="dms_missing")); continue
        df = pd.read_csv(f)
        scores, dropped = pg.wt_marginal(model, tok, r["target_seq"], df["mutant"].tolist(), DEVICE)
        y = df["DMS_score"].values
        ok = np.isfinite(scores) & np.isfinite(y)
        drop_frac = 1 - ok.sum() / len(df)
        status = "ok"
        if ok.sum() < MIN_VARIANTS:
            status = "below_variant_floor"
        elif drop_frac > MAX_DROP_FRAC:
            status = "high_drop_frac"
        elif np.unique(scores[ok]).size < 2:
            status = "constant_scores"                # broken scorer, not a result
        rho = spearmanr(scores[ok], y[ok]).correlation if status == "ok" else np.nan
        rows.append(dict(assay=r["DMS_id"], repo=repo, status=status,
                         params=meta["total"], active_params=meta["active"],
                         MSA_Neff_L=r["MSA_Neff_L"], UniProt_ID=r["UniProt_ID"],
                         seq_len=r["seq_len"], n_variants=len(df),
                         n_dropped=len(dropped), drop_frac=round(drop_frac, 4),
                         rho_marg=rho))
    del model; torch.cuda.empty_cache() if DEVICE == "cuda" else None

scores = pd.DataFrame(rows)
scores.to_csv(ROOT / "results/05_stage3_protgpt3_scores.csv", index=False)
print(f"\nwrote {len(scores)} assay-checkpoint rows")

# %% [markdown]
# ## Technical sanity summary (descriptive rho is NOT a kill criterion)

# %%
ok = scores[scores.status == "ok"]
print("status counts:\n", scores.status.value_counts().to_string())
print("\nglobal parse/universe drop fraction:",
      round(scores["n_dropped"].sum() / scores["n_variants"].sum(), 4)
      if "n_variants" in scores else "n/a")
if len(ok):
    print("\ndescriptive per-checkpoint mean rho (diagnostic only, does not gate the fit):")
    print(ok.groupby("params")["rho_marg"].agg(["mean", "median", "size"]).to_string())
broken = scores[scores.status.isin(["constant_scores", "dms_missing"])]
if len(broken):
    print(f"\nWARNING: {len(broken)} rows flagged (constant/missing) — inspect before the FE step.")

# %%
from src.provenance import stamp  # noqa: E402
stamp("05_stage3_protgpt3_scoring", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, convention="wt_marginal_primary",
      context_max_residues=CTX_MAX_RESIDUES, n_assays=int(len(ref)),
      checkpoints=list(CHECKPOINTS), device=str(DEVICE),
      note="production scoring; no lp:ld computed here")
