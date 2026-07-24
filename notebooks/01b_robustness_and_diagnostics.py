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
# # 01b diagnostics — robustness, influence, bootstrap, class heterogeneity
#
# CPU only. About two minutes, almost all of it the bootstrap.
#
# Run `01b_covariate_controls.py` first; this reproduces everything after the
# confirmatory model. All of it is **post hoc** — motivated by the robustness
# results, not preregistered, and not confirmatory. See
# `experiments/01-msa-depth-confound/results.md`.

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

import os
from pathlib import Path

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
for sub in ("results", "figures"):
    (ROOT / sub).mkdir(exist_ok=True)

PROTEINGYM_COMMIT = os.environ.get(
    "PROTEINGYM_COMMIT", "144fe22b07dfaeec2b366f2346203a9838a55b4c"
)
assert PROTEINGYM_COMMIT != "main" and len(PROTEINGYM_COMMIT) == 40
RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"

BETA_01A = -0.01527576017255374   # reference for attenuation
SEED = 0
UPPER = {"ESM2 (650M)": 650e6, "ESM2 (3B)": 3e9, "ESM2 (15B)": 15e9}
PRIMARY = "rho ~ lp*ld + lp*taxon + lp*z_llen + lp*sel"

# %% [markdown]
# ## Build the assay-level frame
#
# Merge at assay level. Never merge covariates onto depth-bin averages — that
# would be an ecological analysis in which category means hide composition.

# %%
scores = pd.read_csv(
    f"{RAW}/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
)
ref = pd.read_csv(f"{RAW}/reference_files/DMS_substitutions.csv")

d = (
    scores[["DMS ID"] + list(UPPER)]
    .melt("DMS ID", var_name="model", value_name="rho")
    .assign(params=lambda x: x["model"].map(UPPER))
    .merge(
        ref[["DMS_id", "MSA_Neff_L", "MSA_Neff_L_category", "UniProt_ID",
             "taxon", "seq_len", "coarse_selection_type"]],
        left_on="DMS ID", right_on="DMS_id",
    )
    .dropna(subset=["rho", "MSA_Neff_L"])
    .rename(columns={"DMS ID": "assay", "coarse_selection_type": "sel"})
)
for col, src in [("lp", "params"), ("ld", "MSA_Neff_L"), ("llen", "seq_len")]:
    v = np.log10(d[src].clip(lower=0.01))
    d[col] = v - v.mean()
d["z_llen"] = d["llen"] / d["llen"].std()

# %% [markdown]
# ## R1-R6 robustness battery

# %%
def record(model, label, rows):
    b, p = model.params["lp:ld"], model.pvalues["lp:ld"]
    lo, hi = model.conf_int().loc["lp:ld"]
    att = 100 * (abs(BETA_01A) - abs(b)) / abs(BETA_01A)
    print(f"{label:<40s} beta={b:+.4f} att={att:+6.1f}% p={p:.4f} [{lo:+.4f},{hi:+.4f}]")
    rows.append(dict(check=label, beta=b, attenuation_pct=round(att, 1), p=p,
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
    record(ols(PRIMARY, sub, **CL, cov_kwds={"groups": sub.UniProt_ID}), f"R3 drop {t}", rows)
for c in sorted(d.sel.unique()):
    sub = d[d.sel != c]
    record(ols(PRIMARY, sub, **CL, cov_kwds={"groups": sub.UniProt_ID}), f"R4 drop {c}", rows)

infl = smf.ols(PRIMARY, d).fit().get_influence()
d["cooks_d"] = infl.cooks_distance[0]
d["leverage"] = infl.hat_matrix_diag
d["high_influence"] = d["cooks_d"] > 4 / len(d)
keep = d[~d.high_influence]
record(ols(PRIMARY, keep, **CL, cov_kwds={"groups": keep.UniProt_ID}),
       f"R5 drop Cook's D>4/n ({d.high_influence.sum()} obs)", rows)
record(ols("rho ~ lp + lp:ld + C(assay)", d, **CL, cov_kwds={"groups": d.UniProt_ID}),
       "R6 assay fixed effects", rows)

rb = pd.DataFrame(rows)
rb.to_csv(ROOT / "results/01b_robustness_checks.csv", index=False)
rb[["check", "beta", "ci_lo", "ci_hi", "p", "attenuation_pct", "n"]].to_csv(
    ROOT / "results/01b_model_coefficients.csv", index=False)
d.to_csv(ROOT / "results/01b_merged_assay_data.csv", index=False)
print(f"\nsign negative in {(rb.beta < 0).sum()}/{len(rb)} fits; p<0.05 in {(rb.p < 0.05).sum()}/{len(rb)}")

# %% [markdown]
# ## Influence diagnostics
#
# Which observations carry the result, and are they systematically anything?

# %%
d.sort_values("cooks_d", ascending=False)[
    ["assay", "model", "taxon", "sel", "MSA_Neff_L_category", "seq_len",
     "rho", "cooks_d", "leverage", "high_influence"]
].to_csv(ROOT / "results/01b_influence_diagnostics.csv", index=False)

hi = d[d.high_influence]
print(f"{len(hi)} of {len(d)} obs ({len(hi)/len(d):.1%}) across {hi.assay.nunique()} assays\n")
for col in ["taxon", "sel", "MSA_Neff_L_category", "model"]:
    comp = pd.DataFrame({"high_infl_%": 100 * hi[col].value_counts(normalize=True),
                         "overall_%": 100 * d[col].value_counts(normalize=True)})
    comp["ratio"] = (comp["high_infl_%"] / comp["overall_%"]).round(2)
    print(col.upper()); print(comp.dropna().round(1).sort_values("ratio", ascending=False)); print()

# %% [markdown]
# ## Pairwise endpoint decomposition
#
# Is the interaction an endpoint contrast, or present across the ladder?

# %%
post = []
S = {"650M": 650e6, "3B": 3e9, "15B": 15e9}
for a, b in [("650M", "3B"), ("3B", "15B"), ("650M", "15B")]:
    sub = d[d.params.isin([S[a], S[b]])].copy()
    v = np.log10(sub["params"]); sub["lp"] = v - v.mean()
    m = ols(PRIMARY, sub, **CL, cov_kwds={"groups": sub.UniProt_ID})
    be, p = m.params["lp:ld"], m.pvalues["lp:ld"]
    lo, hi_ = m.conf_int().loc["lp:ld"]
    print(f"  {a:>5s} -> {b:<5s} beta={be:+.4f} p={p:.4f} [{lo:+.4f},{hi_:+.4f}]")
    post.append(dict(contrast=f"{a}->{b}", beta=be, p=p, ci_lo=lo, ci_hi=hi_, n=int(m.nobs)))

# %% [markdown]
# ## Assay-level cluster bootstrap
#
# Resample **complete assay clusters**, not individual rows. Each assay
# contributes three correlated observations; resampling rows would break the
# repeated-measures structure and understate uncertainty.

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
lo, med, hi_ = np.percentile(betas, [2.5, 50, 97.5])
print(f"converged {len(betas)}/2000 | median {med:+.4f} | 95% [{lo:+.4f}, {hi_:+.4f}] "
      f"| negative {(betas < 0).sum()}/{len(betas)}")
post.append(dict(contrast="bootstrap (2000 reps)", beta=med, p=np.nan,
                 ci_lo=lo, ci_hi=hi_, n=len(betas)))
pd.DataFrame(post).to_csv(ROOT / "results/01b_posthoc_diagnostics.csv", index=False)

# %% [markdown]
# ## Interaction by selection class
#
# Deleting a class tests whether it was load-bearing. Estimating within each
# class tests whether the effect differs across them. These answer different
# questions and here they point opposite ways.

# %%
cls = []
FC = "rho ~ lp*ld + lp*taxon + lp*z_llen"
for c in sorted(d.sel.unique()):
    sub = d[d.sel == c]
    if sub.assay.nunique() < 8:
        continue
    m = ols(FC, sub, **CL, cov_kwds={"groups": sub.UniProt_ID})
    lo_, hi2 = m.conf_int().loc["lp:ld"]
    cls.append(dict(sel=c, beta=m.params["lp:ld"], se=m.bse["lp:ld"],
                    ci_lo=lo_, ci_hi=hi2, assays=sub.assay.nunique()))
cdf = pd.DataFrame(cls)
print(cdf.round(4).to_string(index=False))

w = 1 / cdf.se ** 2
fe = (w * cdf.beta).sum() / w.sum()
Q = (w * (cdf.beta - fe) ** 2).sum()
k = len(cdf)
I2 = max(0.0, 100 * (Q - (k - 1)) / Q) if Q > 0 else 0.0
tau2 = max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
wr = 1 / (cdf.se ** 2 + tau2)
re, se_re = (wr * cdf.beta).sum() / wr.sum(), np.sqrt(1 / wr.sum())
print(f"\nrandom-effects pooled beta={re:+.4f} 95% CI "
      f"[{re - 1.96 * se_re:+.4f}, {re + 1.96 * se_re:+.4f}]")
print(f"Q={Q:.2f} on {k-1} df, I2={I2:.1f}%")
cdf.to_csv(ROOT / "results/01b_by_selection_class.csv", index=False)

# %% [markdown]
# ## Figures

# %%
fig, ax = plt.subplots(figsize=(7.5, 5.2))
for i, (_, r) in enumerate(rb.iterrows()):
    y = len(rb) - 1 - i
    c = "#2b5797" if r.p < 0.05 else "#b0b0b0"
    ax.plot([r.ci_lo, r.ci_hi], [y, y], color=c, lw=1.6)
    ax.plot([r.ci_lo, r.ci_hi], [y, y], "|", color=c, ms=6)
    ax.plot(r.beta, y, "o", color=c, ms=6, zorder=3)
ax.axvline(0, color="k", lw=0.9)
ax.axvline(BETA_01A, color="#c0392b", ls="--", lw=1, label="01a estimate")
ax.set_yticks(np.arange(len(rb))[::-1]); ax.set_yticklabels(rb.check, fontsize=8)
ax.set_xlabel("size x depth interaction (beta)")
ax.set_title("01b coefficient stability: sign robust, significance fragile\n(grey = CI includes zero)",
             fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="lower right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT / "figures/01b_coefficient_stability.png", dpi=200)
plt.close(fig)

d["tercile"] = pd.qcut(d.MSA_Neff_L, 3, labels=["shallow", "mid", "deep"])
g = d.groupby(["tercile", "params"], observed=True).rho.mean().unstack()
fig, ax = plt.subplots(figsize=(6.4, 4.2))
for t, st, c in [("shallow", "o--", "#2b5797"), ("mid", "s-", "#7f8c8d"), ("deep", "^-", "#c0392b")]:
    ax.plot(np.log10(g.columns), g.loc[t], st, color=c, label=f"{t} MSA")
jit = np.random.default_rng(SEED).uniform(-0.03, 0.03, len(hi))
ax.scatter(np.log10(hi.params) + jit, hi.rho, s=9, c="k", alpha=0.3, zorder=1,
           label=f"high-influence (n={len(hi)})")
ax.set_xlabel("log10 parameters"); ax.set_ylabel("Spearman")
ax.set_title("Above 650M: deep families decline, shallow ones do not", fontsize=10)
ax.legend(frameon=False, fontsize=8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT / "figures/01b_adjusted_interaction.png", dpi=200)
plt.close(fig)
print("figures written")

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("01b_diagnostics", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, bootstrap_seed=SEED, bootstrap_reps=2000)
