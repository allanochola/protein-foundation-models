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
# # 03 — Cross-family robustness: independent set, chance floor, arms 1-2
#
# CPU only, a few minutes (the two 2000-rep bootstraps dominate).
#
# This notebook **executes** the frozen protocol in
# `experiments/03-cross-family-independence/PROTOCOL.md`. It designs nothing.
# Every rule below — the deterministic representative, the 0.10 chance floor,
# the cluster-balanced bootstrap, the random-representative sensitivity arm —
# was committed before this file was written. Nothing here inspects a result
# before its specification is applied.
#
# Covers: the independent set (arm 0), the per-ladder chance floor, arm 1
# (cluster-balanced bootstrap), and arm 2 (random-representative sensitivity).
# Depth balancing and the differential (arms 3-4) are in
# `03_depth_balanced_resampling.py`.
#
# **Implementation decisions not fixed by the protocol** (surfaced, not buried;
# flip them here before running if the reviewer disagrees):
#   - On the independent set each 50%-cluster contributes exactly one assay, so
#     clustering SEs by `homology_cluster_50` equals clustering by assay — it
#     groups the 3-4 size rows of each representative. That is the intended
#     repeated-measures correction, not a coarser one.
#   - Centering is redone on each analysis subset. The interaction `lp:ld` is
#     invariant to centering, so this changes only main-effect readings.

# %%
try:
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf
except ImportError as exc:
    raise RuntimeError("Missing dependencies. Run: pip install -r requirements.txt") from exc

from pathlib import Path

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "results").mkdir(exist_ok=True)

SEED = 0
REPS = 2000
FLOOR = 0.10          # |Spearman| chance floor, fixed in PROTOCOL.md
LADDERS = ["esm2", "progen2", "progen3"]

# Copied unchanged from 01b/01c. Not re-derived.
PRIMARY = "rho ~ lp*ld + lp*taxon + lp*z_llen + lp*sel"

MANIFEST = ROOT / "results/02_proteingym_assay_manifest.csv"
assert MANIFEST.exists(), f"missing {MANIFEST} — Experiment 02 must be committed first"


# %% [markdown]
# ## Loaders
#
# One loader for all three ladders. ESM-2 is stored already-upper-segment in
# `01b_merged_assay_data.csv`; the ProGen ladders carry an `upper_segment`
# flag and a `selection_type` column that 01b's formula calls `sel`. After
# alignment the three are treated identically.

# %%
def load_ladder(name):
    if name == "esm2":
        d = pd.read_csv(ROOT / "results/01b_merged_assay_data.csv")
    else:
        d = pd.read_csv(ROOT / f"results/01c_{name}_assay_data.csv")
        d = d.rename(columns={"selection_type": "sel"})
        d = d[d.upper_segment].copy()
    return recentre(d)


def recentre(d):
    d = d.copy()
    for col, src in [("lp", "params"), ("ld", "MSA_Neff_L"), ("llen", "seq_len")]:
        v = np.log10(d[src].clip(lower=0.01))
        d[col] = v - v.mean()
    d["z_llen"] = d["llen"] / d["llen"].std()
    return d


def beta_ld(data, **fit_kw):
    """Interaction coefficient from one fit, or NaN if the fit fails."""
    try:
        return smf.ols(PRIMARY, data).fit(**fit_kw).params["lp:ld"]
    except Exception:
        return np.nan


def weak_assays(d, thr=FLOOR):
    """Assays where the whole ladder has weak predictive signal (|rho| < thr on
    every model size). Observable weakness — not a claim the assay is noise."""
    g = d.groupby("assay")["rho"].apply(lambda s: (s.abs() < thr).all())
    return set(g[g].index)


# %% [markdown]
# ## Independent set — deterministic representative per 50%-identity cluster
#
# One assay per `homology_cluster_50`: largest `n_single_mutants`, ties broken
# lexicographically by `assay_id`. Model-independent, zero researcher degrees
# of freedom. The rule can shift the depth marginal (it favours well-studied,
# often deeper proteins), so the independent-set depth mix is reported against
# the full 217 here and the random-representative arm below measures the
# sensitivity.

# %%
man = pd.read_csv(MANIFEST)
reps = (man.sort_values(["homology_cluster_50", "n_single_mutants", "assay_id"],
                        ascending=[True, False, True])
           .groupby("homology_cluster_50", as_index=False).first())
INDEP = set(reps.assay_id)
CLUSTER_OF = dict(zip(man.assay_id, man.homology_cluster_50))
MEMBERS = man.groupby("homology_cluster_50")["assay_id"].apply(list).to_dict()
CAT = dict(zip(man.assay_id, man.msa_neff_l_category))
print(f"independent set: {len(INDEP)} representatives from {man.homology_cluster_50.nunique()} clusters")

full_mix = man.msa_neff_l_category.value_counts().to_dict()
indep_mix = reps.msa_neff_l_category.value_counts().to_dict()
order = ["Low", "Medium", "High"]
print("depth mix   full 217 :", {k: full_mix.get(k, 0) for k in order})
print("depth mix   indep set:", {k: indep_mix.get(k, 0) for k in order})

reps_out = reps[["homology_cluster_50", "assay_id", "uniprot_id", "taxon",
                 "selection_type", "msa_neff_l", "msa_neff_l_category",
                 "n_single_mutants"]].copy()
reps_out["rule"] = "max n_single_mutants; ties lexicographic by assay_id"
reps_out.to_csv(ROOT / "results/03_independent_set_manifest.csv", index=False)


# %% [markdown]
# ## Arm 1 — cluster-balanced bootstrap on the independent set
#
# Point estimate: the confirmatory model on the ~178 representatives that
# survive that ladder's chance floor, SEs clustered by `homology_cluster_50`.
# Primary uncertainty: resample the representatives with replacement, carrying
# each one's complete model ladder, 2000 reps, fixed seed. Equal
# biological-system weight, not equal assay weight.

# %%
def independent_frame(d, weak):
    keep = d[d.assay.isin(INDEP) & ~d.assay.isin(weak)].copy()
    return recentre(keep)


ladder_cache = {L: load_ladder(L) for L in LADDERS}
arm1_rows, dropped = [], {}
# The decision table is run BOTH with the chance floor (floor=on) and without
# it (floor=off), per PROTOCOL.md — the with/without contrast is what tells a
# weak-signal dependence apart from a family-composition one.
for apply_floor in (True, False):
    tag = "on" if apply_floor else "off"
    for L in LADDERS:
        d = ladder_cache[L]
        weak = weak_assays(d) if apply_floor else set()
        if apply_floor:
            dropped[L] = sorted(INDEP & weak & set(d.assay.unique()))
        di = independent_frame(d, weak)
        di["cl"] = di.assay.map(CLUSTER_OF)   # one representative per cluster => clusters == assays here
        n_assay = di.assay.nunique()

        m = smf.ols(PRIMARY, di).fit(cov_type="cluster", cov_kwds={"groups": di.cl})
        b, p = m.params["lp:ld"], m.pvalues["lp:ld"]
        lo, hi = m.conf_int().loc["lp:ld"]

        rng = np.random.default_rng(SEED)
        assays = di.assay.unique()
        groups = {a: g for a, g in di.groupby("assay")}
        betas = []
        for _ in range(REPS):
            pick = rng.choice(assays, size=len(assays), replace=True)
            bs = recentre(pd.concat([groups[a] for a in pick], ignore_index=True))
            b_ = beta_ld(bs)
            if not np.isnan(b_):
                betas.append(b_)
        betas = np.array(betas)
        blo, bmed, bhi = np.percentile(betas, [2.5, 50, 97.5])
        pneg = (betas < 0).mean()
        n_drop = len(dropped[L]) if apply_floor else 0
        print(f"floor={tag:3s} {L:8s} n={n_assay:3d} dropped={n_drop:2d} | point {b:+.4f} "
              f"CI[{lo:+.4f},{hi:+.4f}] p={p:.4f} | boot med {bmed:+.4f} "
              f"[{blo:+.4f},{bhi:+.4f}] P(neg)={pneg:.3f} conv={len(betas)}/{REPS}")
        arm1_rows.append(dict(ladder=L, floor=tag, n_assays=n_assay, n_dropped_floor=n_drop,
                              point_beta=b, point_p=p, point_ci_lo=lo, point_ci_hi=hi,
                              boot_median=bmed, boot_ci_lo=blo, boot_ci_hi=bhi,
                              boot_p_negative=pneg, boot_converged=len(betas)))

pd.DataFrame(arm1_rows).to_csv(ROOT / "results/03_cluster_balanced_bootstrap.csv", index=False)


# %% [markdown]
# ## Arm 2 — random-representative sensitivity
#
# The same cluster-balanced bootstrap, but the within-cluster representative is
# drawn uniformly at random each replicate (from members present in the ladder
# that pass the floor) instead of the deterministic rule. Single-assay clusters
# are unaffected; only multi-assay clusters move. If the distribution shifts
# materially against arm 1, the deterministic rule's depth bias is doing real
# work and the writeup must say so.

# %%
arm2_rows = []
for apply_floor in (True, False):
    tag = "on" if apply_floor else "off"
    for L in LADDERS:
        d = ladder_cache[L]
        present = set(d.assay.unique())
        weak = weak_assays(d) if apply_floor else set()
        groups = {a: g for a, g in d.groupby("assay")}
        clusters = sorted({CLUSTER_OF[a] for a in INDEP})
        elig = {c: [a for a in MEMBERS[c] if a in present and a not in weak] for c in clusters}
        elig = {c: v for c, v in elig.items() if v}          # keep clusters with an eligible member
        cl_ids = list(elig)

        rng = np.random.default_rng(SEED)
        betas = []
        for _ in range(REPS):
            pick = rng.choice(cl_ids, size=len(cl_ids), replace=True)
            rows = [groups[elig[c][rng.integers(len(elig[c]))]] for c in pick]
            bs = recentre(pd.concat(rows, ignore_index=True))
            b_ = beta_ld(bs)
            if not np.isnan(b_):
                betas.append(b_)
        betas = np.array(betas)
        blo, bmed, bhi = np.percentile(betas, [2.5, 50, 97.5])
        pneg = (betas < 0).mean()
        print(f"floor={tag:3s} {L:8s} random-rep boot med {bmed:+.4f} [{blo:+.4f},{bhi:+.4f}] "
              f"P(neg)={pneg:.3f} conv={len(betas)}/{REPS}")
        arm2_rows.append(dict(ladder=L, floor=tag, boot_median=bmed, boot_ci_lo=blo,
                              boot_ci_hi=bhi, boot_p_negative=pneg, boot_converged=len(betas)))

pd.DataFrame(arm2_rows).to_csv(ROOT / "results/03_representative_sensitivity.csv", index=False)

# %% [markdown]
# ## Provenance

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("03_independence", out_dir=ROOT / "results",
      independent_set=len(INDEP), floor=FLOOR, bootstrap_reps=REPS, seed=SEED,
      dropped_by_floor={L: dropped[L] for L in LADDERS})
print("\ndone — arms 0-2 written to results/03_*.csv")
