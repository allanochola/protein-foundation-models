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
# # 05 Stage 3 — runtime + cap-stability gate (outcome-blind)
#
# Executes the runtime-feasibility gate in `STAGE3_PROTOCOL.md` **after the
# convention was redesigned to full-LLR primary** (the WT-marginal proxy failed its
# agreement gate, commit `55b0c22`). It measures full-LLR wall time, throughput,
# peak GPU memory, failures, and the **cap-stability** check (does per-assay `rho`
# at N_CAP=2000 reproduce `rho` at 4000?), then reports PASS/FAIL against the numeric
# budget. It computes **no** cross-checkpoint `lp:ld`.
#
# Run each checkpoint with `CHECKPT=112M|1.3B|10B`. The 10B needs an A100-40GB — that
# is the binding memory + GPU-hour gate under full-LLR primary.

# %%
import os, time
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
print(pg.selftest())

(ROOT / "results").mkdir(exist_ok=True)
PROTEINGYM_COMMIT = os.environ.get("PROTEINGYM_COMMIT",
                                   "144fe22b07dfaeec2b366f2346203a9838a55b4c")
DMS_DIR = Path(os.environ.get("DMS_DIR", "/content/ProteinGym/DMS_ProteinGym_substitutions"))
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
CHUNK = int(os.environ.get("CHUNK", "16"))
SEED, N_CAP = 0, 2000

CHECKPOINTS = {
    "AI4PD/ProtGPT3-112M": dict(total=0.109e9, active=0.034e9),
    "AI4PD/ProtGPT3-1.3B": dict(total=1.328e9, active=0.366e9),
    "AI4PD/ProtGPT3-10B":  dict(total=10.000e9, active=2.752e9),
}
_sel = os.environ.get("CHECKPT")
if _sel:
    CHECKPOINTS = {k: v for k, v in CHECKPOINTS.items() if _sel in k}
    assert CHECKPOINTS, f"CHECKPT={_sel} matched no checkpoint"

PANEL = ["OTU7A_HUMAN_Tsuboyama_2023_2L2D", "ESTA_BACSU_Nutschel_2020",
         "RAF1_HUMAN_Zinkus-Boltz_2019", "SPG1_STRSG_Olson_2014",
         "A0A192B1T2_9HIV1_Haddox_2018"]

GATE = dict(reserved_gb_max=40.0, prod_hours_max=12.0, drho_max=0.03)

# %%
REF_URL = (f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/"
           f"{PROTEINGYM_COMMIT}/reference_files/DMS_substitutions.csv")
ref = pd.read_csv(REF_URL).set_index("DMS_id")
print("device:", torch.cuda.get_device_name(0) if DEVICE == "cuda" else "cpu")

def usable(df, wt):
    """Rows whose mutant parses under the 20-AA universe (matches production filter)."""
    keep = [pg.parse_variant(m, wt) is not None for m in df["mutant"].astype(str)]
    return df[keep].reset_index(drop=True)

# %%
gate_rows, cap_rows = [], []
for repo, meta in CHECKPOINTS.items():
    print(f"\n=== {repo} ===")
    if DEVICE == "cuda":
        torch.cuda.reset_peak_memory_stats(); torch.cuda.empty_cache()
    model, tok = pg.load(repo, device=DEVICE)
    wall = 0.0; n_scored = 0
    for dms in PANEL:
        df = usable(pd.read_csv(DMS_DIR / f"{dms}.csv"), ref.loc[dms, "target_seq"])
        wt = ref.loc[dms, "target_seq"]; n = len(df)
        cap2, idxB = pg.cap_pair(n, N_CAP, SEED)          # cap2 (production) nested in idxB (<=4000)
        seqsB = df["mutated_sequence"].values[idxB]
        yB = df["DMS_score"].values[idxB]
        t = time.time()
        sB = pg.full_llr(model, tok, wt, seqsB, DEVICE, CHUNK)
        wall += time.time() - t; n_scored += len(idxB)
        # cap-stability: rho at first 2000 (subset of idxB) vs rho at all of idxB (<=4000)
        if len(idxB) > N_CAP:
            a = slice(0, N_CAP)
            rho2 = spearmanr(sB[a], yB[a]).correlation
            rho4 = spearmanr(sB, yB).correlation
            cap_rows.append(dict(repo=repo, assay=dms, n_usable=n, n2=N_CAP, n4=len(idxB),
                                 rho2=rho2, rho4=rho4, drho=abs(rho2 - rho4),
                                 sign_ok=(np.sign(rho2) == np.sign(rho4))))
            print(f"  {dms:34s} n={n:>7d}  |rho2-rho4|={abs(rho2-rho4):.3f}  ({wall:.0f}s)")
        else:
            print(f"  {dms:34s} n={n:>7d}  (<=2000: scored in full, no cap)  ({wall:.0f}s)")
    peak_res = torch.cuda.max_memory_reserved() / 1e9 if DEVICE == "cuda" else np.nan
    peak_alloc = torch.cuda.max_memory_allocated() / 1e9 if DEVICE == "cuda" else np.nan
    vps = n_scored / wall if wall else np.nan
    proj_h = (217 * N_CAP) / vps / 3600 if vps == vps and vps > 0 else np.nan   # capped production
    gate_rows.append(dict(repo=repo, total_params=meta["total"], active_params=meta["active"],
                          full_llr_wall_s=wall, variants_per_s=vps,
                          peak_reserved_gb=peak_res, peak_allocated_gb=peak_alloc,
                          proj_prod_hours=proj_h))
    print(f"  peak reserved={peak_res:.1f}GB allocated={peak_alloc:.1f}GB "
          f"| {vps:.0f} var/s | proj capped-production={proj_h:.2f}h")
    del model; torch.cuda.empty_cache() if DEVICE == "cuda" else None

G = pd.DataFrame(gate_rows); C = pd.DataFrame(cap_rows)
G.to_csv(ROOT / "results/05_stage3_runtime_gate.csv", index=False)
C.to_csv(ROOT / "results/05_stage3_capstability.csv", index=False)

# %% [markdown]
# ## PASS / FAIL against the frozen budget and cap-stability
#
# The 10B is the binding memory + GPU-hour check; per-checkpoint runs accumulate.

# %%
have = set(G["repo"].str.extract(r"(112M|1\.3B|10B)")[0])
worst_drho = C["drho"].max() if len(C) else np.nan
sign_ok_all = bool(C["sign_ok"].all()) if len(C) else True
total_proj = G["proj_prod_hours"].sum()
ten = G[G.repo.str.contains("10B")]

checks = {
    "10B reserved <= 40GB": (ten["peak_reserved_gb"].max() <= GATE["reserved_gb_max"]) if len(ten) else None,
    "capped production <= 12 A100-h (sum)": bool(total_proj <= GATE["prod_hours_max"]),
    "cap-stability |drho| <= 0.03 (all)": bool(worst_drho <= GATE["drho_max"]) if len(C) else None,
    "cap-stability sign agreement (all)": sign_ok_all,
}
print("checkpoints present:", sorted(have))
print(f"worst |rho2-rho4|={worst_drho:.3f} | sign-ok-all={sign_ok_all} "
      f"| summed proj={total_proj:.2f}h | worst reserved={np.nanmax(G['peak_reserved_gb']):.1f}GB\n")
for k, v in checks.items():
    print(f"  [{'PASS' if v is True else 'FAIL' if v is False else '—pending—'}] {k}")

missing = {"112M", "1.3B", "10B"} - have
verdict = ("PASS" if (not missing and all(v is True for v in checks.values())) else
           "FAIL" if any(v is False for v in checks.values()) else "INCOMPLETE")
print(f"\nGATE: {verdict}" + (f"  (missing: {sorted(missing)})" if missing else ""))
if verdict == "FAIL":
    print("Per protocol: record and STOP/redesign — raise N_CAP and re-gate, or shrink "
          "the assay set transparently. Do NOT silently narrow the ladder or proceed.")

# %%
from src.provenance import stamp  # noqa: E402
stamp("05_stage3_runtime_gate", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, panel=PANEL, thresholds=GATE,
      n_cap=N_CAP, seed=SEED, device=str(DEVICE),
      checkpoints_run=sorted(have), verdict=verdict,
      convention="full_llr_capped", note="outcome-blind: no cross-checkpoint lp:ld")
