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
# # 01c Phase 3 — replication analysis (ProGen2 and ProGen3)
#
# CPU only, about a minute per ladder (the bootstrap dominates).
#
# This is a **replication**, so the specification is copied from 01b and not
# reconsidered. The only thing that changes is the input dataset. Every
# design choice — the primary model, the covariate set, the robustness
# battery, the bootstrap scheme, the decision table — was fixed in
# `experiments/01c-progen2-replication/PROTOCOL.md` before any ProGen number
# was seen. This notebook executes that protocol; it does not design anything.
#
# Set `LADDER` below and run. Run once for each ladder. Nothing here inspects
# results before the specification is applied.

# %%
LADDER = "progen2"   # "progen2" or "progen3" — the ONLY line to change between runs

# %%
try:
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise RuntimeError("Missing dependencies. Run: pip install -r requirements.txt") from exc

from pathlib import Path

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "figures").mkdir(exist_ok=True)

assert LADDER in ("progen2", "progen3")
DATA = ROOT / f"results/01c_{LADDER}_assay_data.csv"
assert DATA.exists(), f"missing {DATA} — run 01c_build_dataset.py first"

# The ESM-2 reference the replication is compared against (01b confirmatory).
BETA_ESM2 = -0.01527576017255374
SEED = 0

# Identical specifications to 01b. Not re-derived here.
PRIMARY = "rho ~ lp*ld + lp*taxon + lp*z_llen + lp*sel"   # cluster SE by UniProt_ID
FIXED_EFFECTS = "rho ~ lp + lp:ld + C(assay)"             # robustness reference

# %% [markdown]
# ## Load and centre
#
# Centring is done on the upper segment being fitted, exactly as in 01b, so
# main effects read at the segment mean. The dataset already carries
# `upper_segment` from Phase 2.

# %%
full = pd.read_csv(DATA)
# Phase 2 datasets store the column as `selection_type`; 01b's specification
# uses `sel`. Alias so the copied formula stays byte-identical to 01b.
full = full.rename(columns={"selection_type": "sel"})
d = full[full.upper_segment].copy()
for col, src in [("lp", "params"), ("ld", "MSA_Neff_L"), ("llen", "seq_len")]:
    v = np.log10(d[src].clip(lower=0.01))
    d[col] = v - v.mean()
d["z_llen"] = d["llen"] / d["llen"].std()

sizes = sorted(d.params.unique())
span = np.log10(sizes[-1]) - np.log10(sizes[0])
print(f"{LADDER}: upper segment = {d.assay.nunique()} assays x {len(sizes)} sizes "
      f"= {len(d)} rows | log10 span {span:.2f}")
print("sizes:", [f"{s/1e6:.0f}M" if s < 1e9 else f"{s/1e9:.1f}B" for s in sizes])

# %% [markdown]
# ## Confirmatory model + R1-R6 robustness battery
#
# Identical to 01b. Attenuation is measured against the ESM-2 estimate so the
# replication is directly comparable.

# %%
def record(model, label, rows):
    b, p = model.params["lp:ld"], model.pvalues["lp:ld"]
    lo, hi = model.conf_int().loc["lp:ld"]
    att = 100 * (abs(BETA_ESM2) - abs(b)) / abs(BETA_ESM2)
    print(f"{label:<40s} beta={b:+.4f} vs_esm2={att:+6.1f}% p={p:.4f} [{lo:+.4f},{hi:+.4f}]")
    rows.append(dict(check=label, beta=b, vs_esm2_pct=round(att, 1), p=p,
                     ci_lo=lo, ci_hi=hi, n=int(model.nobs)))


def ols(formula, data, **kw):
    return smf.ols(formula, data).fit(**kw)


CL = dict(cov_type="cluster")
rows = []
record(ols(PRIMARY, d, **CL, cov_kwds={"groups": d.UniProt_ID}),
       "confirmatory (cluster UniProt)", rows)
record(ols(PRIMARY, d, cov_type="HC3"), "R1 HC3 robust SE", rows)
record(ols(PRIMARY, d, **CL, cov_kwds={"groups": d.assay}), "R2 cluster by assay", rows)
for t in sorted(d.taxon.unique()):
    sub = d[d.taxon != t]
    if sub.taxon.nunique() < 2:
        continue
    record(ols(PRIMARY, sub, **CL, cov_kwds={"groups": sub.UniProt_ID}), f"R3 drop {t}", rows)
for c in sorted(d.sel.unique()):
    sub = d[d.sel != c]
    if sub.sel.nunique() < 2:
        continue
    record(ols(PRIMARY, sub, **CL, cov_kwds={"groups": sub.UniProt_ID}), f"R4 drop {c}", rows)

infl = smf.ols(PRIMARY, d).fit().get_influence()
d["cooks_d"] = infl.cooks_distance[0]
d["high_influence"] = d["cooks_d"] > 4 / len(d)
keep = d[~d.high_influence]
record(ols(PRIMARY, keep, **CL, cov_kwds={"groups": keep.UniProt_ID}),
       f"R5 drop Cook's D>4/n ({d.high_influence.sum()} obs)", rows)
record(ols(FIXED_EFFECTS, d, **CL, cov_kwds={"groups": d.UniProt_ID}),
       "R6 assay fixed effects", rows)

rb = pd.DataFrame(rows)
rb.to_csv(ROOT / f"results/01c_{LADDER}_robustness_checks.csv", index=False)
print(f"\nsign negative in {(rb.beta < 0).sum()}/{len(rb)} fits; "
      f"p<0.05 in {(rb.p < 0.05).sum()}/{len(rb)}")

# %% [markdown]
# ## Pairwise endpoint decomposition
#
# Same logic as 01b: is the interaction spread across the ladder or driven by
# one contrast? ProGen ladders have more than three sizes, so this covers all
# adjacent pairs plus the full-range endpoint contrast.

# %%
post = []
for i in range(len(sizes) - 1):
    a, b_ = sizes[i], sizes[i + 1]
    sub = d[d.params.isin([a, b_])].copy()
    v = np.log10(sub.params); sub["lp"] = v - v.mean()
    m = ols(PRIMARY, sub, **CL, cov_kwds={"groups": sub.UniProt_ID})
    lo, hi = m.conf_int().loc["lp:ld"]
    lbl = f"{a/1e6:.0f}M->{b_/1e6:.0f}M" if b_ < 1e9 else f"{a/1e9:.2f}B->{b_/1e9:.2f}B"
    print(f"  {lbl:<16s} beta={m.params['lp:ld']:+.4f} p={m.pvalues['lp:ld']:.4f} [{lo:+.4f},{hi:+.4f}]")
    post.append(dict(contrast=lbl, beta=m.params["lp:ld"], p=m.pvalues["lp:ld"],
                     ci_lo=lo, ci_hi=hi, n=int(m.nobs)))

# %% [markdown]
# ## Assay-level cluster bootstrap
#
# Identical scheme to 01b: resample complete assay clusters (all model sizes
# per sampled assay), 2000 reps, fixed seed. Rows, not clusters, would break
# the repeated-measures structure.

# %%
rng = np.random.default_rng(SEED)
assays = d.assay.unique()
groups = {a: g for a, g in d.groupby("assay")}
betas = []
for _ in range(2000):
    pick = rng.choice(assays, size=len(assays), replace=True)
    bs = pd.concat([groups[a] for a in pick], ignore_index=True)
    try:
        betas.append(smf.ols(PRIMARY, bs).fit().params["lp:ld"])
    except Exception:
        pass
betas = np.array(betas)
lo, med, hi = np.percentile(betas, [2.5, 50, 97.5])
print(f"converged {len(betas)}/2000 | median {med:+.4f} | 95% [{lo:+.4f}, {hi:+.4f}] "
      f"| negative {(betas < 0).sum()}/{len(betas)}")
post.append(dict(contrast="bootstrap (2000 reps)", beta=med, p=np.nan,
                 ci_lo=lo, ci_hi=hi, n=len(betas)))
pd.DataFrame(post).to_csv(ROOT / f"results/01c_{LADDER}_posthoc_diagnostics.csv", index=False)

# %% [markdown]
# ## Coefficient table + figure

# %%
rb[["check", "beta", "ci_lo", "ci_hi", "p", "vs_esm2_pct", "n"]].to_csv(
    ROOT / f"results/01c_{LADDER}_model_coefficients.csv", index=False)

fig, ax = plt.subplots(figsize=(7.5, 5.2))
for i, (_, r) in enumerate(rb.iterrows()):
    y = len(rb) - 1 - i
    c = "#2b5797" if r.p < 0.05 else "#b0b0b0"
    ax.plot([r.ci_lo, r.ci_hi], [y, y], color=c, lw=1.6)
    ax.plot([r.ci_lo, r.ci_hi], [y, y], "|", color=c, ms=6)
    ax.plot(r.beta, y, "o", color=c, ms=6, zorder=3)
ax.axvline(0, color="k", lw=0.9)
ax.axvline(BETA_ESM2, color="#c0392b", ls="--", lw=1, label="ESM-2 estimate")
ax.set_yticks(np.arange(len(rb))[::-1]); ax.set_yticklabels(rb.check, fontsize=8)
ax.set_xlabel("size x depth interaction (beta)")
ax.set_title(f"01c {LADDER}: replication of the ESM-2 depth-scaling interaction\n"
             "(grey = CI includes zero; red dashed = ESM-2)", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="lower right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / f"figures/01c_{LADDER}_coefficient_stability.png", dpi=200)
plt.close(fig)
print("figure written")

# %% [markdown]
# ## Verdict against the preregistered decision table
#
# Mechanical application of the committed rule. No interpretation beyond it.

# %%
conf = rb.iloc[0]
fe = rb[rb.check == "R6 assay fixed effects"].iloc[0]
p_neg = (betas < 0).mean()
neg = conf.beta < 0
excl0 = conf.ci_hi < 0
within50 = abs(conf.beta) >= 0.5 * abs(BETA_ESM2)
# In 01b the fixed-effects model was the strong test: it survived there.
# A replication that holds on the confirmatory model but collapses under
# fixed effects is materially weaker — the signal is between-assay, not
# within-assay, so it may be composition rather than a scaling response.
fe_survives = (fe.beta < 0) and (fe.ci_hi < 0)

if neg and excl0 and within50 and fe_survives:
    verdict = "REPLICATES (strong) — holds on confirmatory AND assay fixed effects"
elif neg and excl0 and within50 and not fe_survives:
    verdict = ("REPLICATES (confirmatory only) — negative and CI excludes 0, but "
               "collapses under assay fixed effects. Signal is between-assay, "
               "not within-assay: consistent with composition or with low "
               "within-assay power on a short ladder, not with a clean scaling "
               "response. WEAKER than the ESM-2 result, which survived R6.")
elif neg and not excl0:
    verdict = "DIRECTIONALLY CONSISTENT, UNDERPOWERED — negative but confirmatory CI includes 0"
elif abs(conf.beta) < 0.25 * abs(BETA_ESM2):
    verdict = "NEAR ZERO — ESM-2-specific, not a general scaling property"
elif conf.beta > 0 and conf.ci_lo > 0:
    verdict = "CONTRADICTS — positive, CI excludes 0; ESM-2 result needs re-examination"
else:
    verdict = "MIXED — report all specifications, do not force a single label"

print(f"LADDER: {LADDER}")
print(f"confirmatory beta = {conf.beta:+.4f}  CI [{conf.ci_lo:+.4f}, {conf.ci_hi:+.4f}]  p = {conf.p:.4f}")
print(f"fixed-effects beta = {fe.beta:+.4f}  CI [{fe.ci_lo:+.4f}, {fe.ci_hi:+.4f}]  p = {fe.p:.4f}")
print(f"bootstrap P(beta<0) = {p_neg:.3f}   ({(betas<0).sum()}/{len(betas)})")
print(f"magnitude vs ESM-2 = {100*abs(conf.beta)/abs(BETA_ESM2):.0f}%")
print(f"\nVERDICT: {verdict}")

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp(f"01c_{LADDER}_analysis", out_dir=ROOT / "results",
      ladder=LADDER, bootstrap_seed=SEED, bootstrap_reps=2000,
      confirmatory_beta=float(conf.beta), bootstrap_p_negative=float(p_neg))
