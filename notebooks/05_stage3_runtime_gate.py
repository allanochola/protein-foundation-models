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
# # 05 Stage 3 — runtime + agreement gate (outcome-blind)
#
# Executes the runtime-feasibility gate in `STAGE3_PROTOCOL.md`. It measures wall
# time, throughput, peak GPU memory, failures, and WT-marginal-vs-full-LLR scoring
# agreement on the **five frozen panel assays**, then reports PASS/FAIL against the
# numeric budget and the locked agreement thresholds. It computes **no**
# cross-checkpoint `lp:ld` — it never assembles the three-point ladder, and it
# touches only five assays chosen by geometry.
#
# Run each checkpoint with `CHECKPT=112M|1.3B|10B`. The 10B step needs an A100-40GB
# (that is the memory gate). For the GPU-hour budget to be in A100-hours, run the
# gate (at least the 10B) on an A100.

# %%
import os, time, json
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
print(pg.selftest())                                  # synthetic unit test before anything runs

(ROOT / "results").mkdir(exist_ok=True)
PROTEINGYM_COMMIT = os.environ.get("PROTEINGYM_COMMIT",
                                   "144fe22b07dfaeec2b366f2346203a9838a55b4c")
DMS_DIR = Path(os.environ.get("DMS_DIR", "/content/ProteinGym/DMS_ProteinGym_substitutions"))
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
CHUNK = int(os.environ.get("CHUNK", "16"))
SEED, N_SUB = 0, 2000                                 # full-LLR subsample per panel assay

CHECKPOINTS = {
    "AI4PD/ProtGPT3-112M": dict(total=0.109e9, active=0.034e9),
    "AI4PD/ProtGPT3-1.3B": dict(total=1.328e9, active=0.366e9),
    "AI4PD/ProtGPT3-10B":  dict(total=10.000e9, active=2.752e9),
}
_sel = os.environ.get("CHECKPT")
if _sel:
    CHECKPOINTS = {k: v for k, v in CHECKPOINTS.items() if _sel in k}
    assert CHECKPOINTS, f"CHECKPT={_sel} matched no checkpoint"

# Frozen panel (STAGE3_PROTOCOL.md), selected by geometry — NOT by rho.
PANEL = ["OTU7A_HUMAN_Tsuboyama_2023_2L2D",   # short / low-variant
         "ESTA_BACSU_Nutschel_2020",           # medium
         "RAF1_HUMAN_Zinkus-Boltz_2019",       # long / low-variant
         "SPG1_STRSG_Olson_2014",              # high-variant (double-mutant landscape)
         "A0A192B1T2_9HIV1_Haddox_2018"]       # long + high-variant stress

# Locked thresholds.
GATE = dict(reserved_gb_max=40.0, wt_prod_hours_max=3.0, llr_val_hours_max=8.0,
            var_spearman_min=0.90, assay_spearman_min=0.80)
LLR_10B_MAX_VARIANTS = 5000                            # 10B full-LLR only on small panel assays

# %%
REF_URL = (f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/"
           f"{PROTEINGYM_COMMIT}/reference_files/DMS_substitutions.csv")
ref = pd.read_csv(REF_URL).set_index("DMS_id")
print("device:", torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu")

def load_assay(dms_id):
    df = pd.read_csv(DMS_DIR / f"{dms_id}.csv")
    return df, ref.loc[dms_id, "target_seq"]

# %%
gate_rows, agree_rows = [], []
for repo, meta in CHECKPOINTS.items():
    is10b = "10B" in repo
    print(f"\n=== {repo} ===")
    if DEVICE == "cuda":
        torch.cuda.reset_peak_memory_stats(); torch.cuda.empty_cache()
    model, tok = pg.load(repo, device=DEVICE)
    marg_time = llr_time = 0.0; n_marg = n_llr = 0
    for dms in PANEL:
        df, wt = load_assay(dms)
        n = len(df)
        rng = np.random.default_rng(SEED)
        sub = np.sort(rng.choice(n, size=min(n, N_SUB), replace=False))

        t = time.time()
        marg_all, dropped = pg.wt_marginal(model, tok, wt, df["mutant"].tolist(), DEVICE)
        marg_time += time.time() - t; n_marg += n

        run_llr = not (is10b and n > LLR_10B_MAX_VARIANTS)   # skip 10B full-LLR on huge assays
        if run_llr:
            t = time.time()
            llr_sub = pg.full_llr(model, tok, wt, df["mutated_sequence"].values[sub], DEVICE, CHUNK)
            llr_time += time.time() - t; n_llr += len(sub)
        else:
            llr_sub = np.full(len(sub), np.nan)

        y = df["DMS_score"].values
        marg_sub = marg_all[sub]
        fin = np.isfinite(marg_sub) & np.isfinite(llr_sub)
        var_sp = spearmanr(marg_sub[fin], llr_sub[fin]).correlation if fin.sum() > 10 else np.nan
        fm = np.isfinite(marg_all) & np.isfinite(y)
        rho_marg = spearmanr(marg_all[fm], y[fm]).correlation if fm.sum() > 10 else np.nan
        fl = np.isfinite(llr_sub) & np.isfinite(y[sub])
        rho_llr = spearmanr(llr_sub[fl], y[sub][fl]).correlation if fl.sum() > 10 else np.nan
        agree_rows.append(dict(repo=repo, assay=dms, n=n, n_sub=int(len(sub)),
                               dropped=len(dropped), var_spearman=var_sp,
                               rho_marg=rho_marg, rho_llr=rho_llr,
                               sign_agree=(np.sign(rho_marg) == np.sign(rho_llr))
                               if np.isfinite(rho_marg) and np.isfinite(rho_llr) else np.nan,
                               llr_run=run_llr))
        print(f"  {dms:38s} n={n:>6d} var_sp={var_sp if var_sp==var_sp else float('nan'):.3f}"
              f"  marg={marg_time:.0f}s llr={llr_time:.0f}s")

    peak_res = torch.cuda.max_memory_reserved() / 1e9 if DEVICE == "cuda" else np.nan
    peak_alloc = torch.cuda.max_memory_allocated() / 1e9 if DEVICE == "cuda" else np.nan
    # projections: WT-marginal is one forward/assay -> mean panel marg time x 217 (conservative)
    wt_prod_h = (marg_time / len(PANEL)) * 217 / 3600
    llr_val_h = llr_time / 3600                          # panel IS the validation set
    gate_rows.append(dict(repo=repo, total_params=meta["total"], active_params=meta["active"],
                          marg_wall_s=marg_time, llr_wall_s=llr_time,
                          marg_vars_per_s=n_marg / marg_time if marg_time else np.nan,
                          llr_vars_per_s=n_llr / llr_time if llr_time else np.nan,
                          peak_reserved_gb=peak_res, peak_allocated_gb=peak_alloc,
                          proj_wt_prod_hours=wt_prod_h, proj_llr_val_hours=llr_val_h))
    print(f"  peak reserved={peak_res:.1f}GB allocated={peak_alloc:.1f}GB "
          f"| proj WT-prod={wt_prod_h:.2f}h  full-LLR-val={llr_val_h:.2f}h")
    del model; torch.cuda.empty_cache() if DEVICE == "cuda" else None

G = pd.DataFrame(gate_rows); A = pd.DataFrame(agree_rows)
G.to_csv(ROOT / "results/05_stage3_runtime_gate.csv", index=False)
A.to_csv(ROOT / "results/05_stage3_agreement.csv", index=False)

# %% [markdown]
# ## PASS / FAIL against the frozen budget and agreement thresholds
#
# Only meaningful once all three checkpoints have been run (the 10B memory gate is
# the binding one). Per-checkpoint runs accumulate into the two CSVs; this cell
# evaluates whatever is present and says what is still missing.

# %%
have = set(G["repo"].str.extract(r"(112M|1\.3B|10B)")[0])
median_var = np.nanmedian(A["var_spearman"].values)
assay_sp = (spearmanr(A["rho_marg"], A["rho_llr"], nan_policy="omit").correlation
            if A["rho_llr"].notna().sum() > 2 else np.nan)
sign_all = bool(A["sign_agree"].dropna().all()) if A["sign_agree"].notna().any() else False
worst_reserved = np.nanmax(G["peak_reserved_gb"].values)

checks = {
    "10B reserved <= 40GB": (G.loc[G.repo.str.contains("10B"), "peak_reserved_gb"].max() <= GATE["reserved_gb_max"])
                            if G.repo.str.contains("10B").any() else None,
    "WT-prod proj <= 3h":   bool((G["proj_wt_prod_hours"] <= GATE["wt_prod_hours_max"]).all()),
    "full-LLR val proj <= 8h": bool(G["proj_llr_val_hours"].sum() <= GATE["llr_val_hours_max"]),
    "variant median Spearman >= 0.90": bool(median_var >= GATE["var_spearman_min"]),
    "assay Spearman >= 0.80": bool(assay_sp >= GATE["assay_spearman_min"]) if assay_sp == assay_sp else None,
    "assay sign agreement (all)": sign_all,
}
print("checkpoints present:", sorted(have))
print(f"median variant Spearman={median_var:.3f} | assay Spearman={assay_sp:.3f} "
      f"| sign-agree-all={sign_all} | worst reserved={worst_reserved:.1f}GB\n")
for k, v in checks.items():
    print(f"  [{'PASS' if v is True else 'FAIL' if v is False else '—pending—'}] {k}")

missing = {"112M", "1.3B", "10B"} - have
decided = [v for v in checks.values() if v is not None]
verdict = ("PASS" if (not missing and all(decided)) else
           "FAIL" if any(v is False for v in checks.values()) else "INCOMPLETE")
print(f"\nGATE: {verdict}" + (f"  (missing checkpoints: {sorted(missing)})" if missing else ""))
if verdict == "FAIL":
    print("Per protocol: record the feasibility failure and STOP/redesign — do NOT "
          "silently narrow to 112M->1.3B, and do NOT proceed to production scoring.")

# %%
from src.provenance import stamp  # noqa: E402
stamp("05_stage3_runtime_gate", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, panel=PANEL, thresholds=GATE,
      subsample=N_SUB, seed=SEED, device=str(DEVICE),
      checkpoints_run=sorted(have), verdict=verdict,
      note="outcome-blind: no cross-checkpoint lp:ld computed")
