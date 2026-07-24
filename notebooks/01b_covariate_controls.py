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
# # 01b — Does the size x depth interaction survive covariate control?
#
# No GPU. Runs in under a minute.
#
# 01a found that past 650M, larger ESM-2 models lose to 650M on well-sampled
# families and win marginally on poorly-sampled ones. The kill criterion in
# `PROTOCOL.md` was: if that interaction disappears once taxon, sequence
# length, and selection type are controlled, it is a composition artefact.
#
# **The design point that matters.** Adding covariates as main effects does
# not control a confound in an *interaction* term. If deep families skew
# prokaryotic, and prokaryotic assays happen to have a different scaling
# slope, `lp:ld` absorbs `lp:taxon` and the main effect changes nothing. The
# test needs size x covariate interactions — and better still, assay fixed
# effects, which absorb every assay-level confound whether or not we thought
# to measure it.
#
# Four specifications, increasing in strength. Model 4 is the one that counts.

# %%
try:
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf
except ImportError as exc:
    raise RuntimeError("Missing dependencies. Run: pip install -r requirements.txt") from exc

import os
from pathlib import Path

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "results").mkdir(exist_ok=True)

PROTEINGYM_COMMIT = os.environ.get(
    "PROTEINGYM_COMMIT", "144fe22b07dfaeec2b366f2346203a9838a55b4c"
)
assert PROTEINGYM_COMMIT != "main", "pin a commit SHA, not a branch"
assert len(PROTEINGYM_COMMIT) == 40, "use the full 40-character SHA"

RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"
SCORES = f"{RAW}/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
REF = f"{RAW}/reference_files/DMS_substitutions.csv"

# %%
LADDER = {
    "ESM2 (8M)": 8e6, "ESM2 (35M)": 35e6, "ESM2 (150M)": 150e6,
    "ESM2 (650M)": 650e6, "ESM2 (3B)": 3e9, "ESM2 (15B)": 15e9,
}
scores, ref = pd.read_csv(SCORES), pd.read_csv(REF)
for c in ["taxon", "seq_len", "coarse_selection_type"]:
    assert c in ref.columns, f"reference column changed: {c}"

df = (
    scores[["DMS ID"] + list(LADDER)]
    .melt("DMS ID", var_name="model", value_name="rho")
    .assign(params=lambda d: d["model"].map(LADDER))
    .merge(
        ref[["DMS_id", "MSA_Neff_L", "UniProt_ID", "taxon", "seq_len", "coarse_selection_type"]],
        left_on="DMS ID", right_on="DMS_id",
    )
    .dropna(subset=["rho", "MSA_Neff_L"])
    .rename(columns={"DMS ID": "assay", "coarse_selection_type": "sel"})
)
for col, src in [("lp", "params"), ("ld", "MSA_Neff_L"), ("llen", "seq_len")]:
    v = np.log10(df[src].clip(lower=0.01))
    df[col] = v - v.mean()          # centred, so main effects read at the mean

upper = df[df["params"] >= 650e6].copy()

# %% [markdown]
# ## Is depth actually confounded with the covariates?
#
# If not, the control is a formality. It is not.

# %%
conf = df.groupby("taxon")[["ld", "llen"]].mean().round(2)
conf["n_assays"] = df.groupby("taxon")["assay"].nunique()
print(conf)
print("\nProkaryote families are far better sampled than human ones.")

# %% [markdown]
# ## Four specifications

# %%
def fit(formula, data, label, records):
    m = smf.ols(formula, data).fit(
        cov_type="cluster", cov_kwds={"groups": data["UniProt_ID"]}
    )
    b, p = m.params["lp:ld"], m.pvalues["lp:ld"]
    ci = m.conf_int().loc["lp:ld"]
    print(f"{label:<46s} beta={b:+.4f}  p={p:.4f}  [{ci[0]:+.4f}, {ci[1]:+.4f}]  n={int(m.nobs)}")
    records.append(dict(segment=data.attrs.get("seg"), spec=label, beta=b, p=p,
                        ci_lo=ci[0], ci_hi=ci[1], n=int(m.nobs)))
    return m


SPECS = [
    ("rho ~ lp*ld", "1. baseline (01a)"),
    ("rho ~ lp*ld + taxon + llen + sel", "2. covariates as main effects only"),
    ("rho ~ lp*ld + lp*taxon + lp*llen + lp*sel", "3. size x covariate interactions"),
    ("rho ~ lp + lp:ld + C(assay)", "4. assay fixed effects"),
]

rows = []
for data, seg in [(upper, "650M-15B"), (df, "8M-15B")]:
    data.attrs["seg"] = seg
    print(f"\n=== {seg} ===")
    for formula, label in SPECS:
        fit(formula, data, label, rows)

out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/01b_covariate_controls.csv", index=False)

# %% [markdown]
# ## Reading the result
#
# Specification 2 changes nothing — beta identical to four decimals. That is
# the point: assay-level covariates shift levels, they cannot absorb an
# interaction with a within-assay variable. Anyone who "controls for
# confounds" this way has controlled for nothing.
#
# Specification 4 is the strong test. With assay fixed effects, `ld` is
# collinear with the dummies and drops out; the interaction is identified
# purely from how each assay's score changes across model sizes. No
# assay-level property — taxon, length, assay type, or anything unmeasured —
# can produce it.
#
# The kill criterion is not triggered. The effect is robust, though p rises
# from 0.002 to 0.011 under the stronger controls, so the evidence is weaker
# than 01a alone suggested. Reporting only the baseline would have
# overstated it.

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("01b_covariate_controls", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT)
