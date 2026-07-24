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
# # 01a — Does ESM-2 scaling survive a control for MSA depth?
#
# No GPU. Runs in about a minute on Colab CPU.
#
# ProteinGym publishes per-assay zero-shot Spearman for the full ESM-2 size
# ladder (8M through 15B) and, separately, MSA depth for every assay. Joining
# them answers the question without running a single model.
#
# **Hypothesis.** Reported ESM-2 scaling gains are partly a coverage effect:
# bigger models help most where the protein family is already well sampled.
#
# **Kill criterion.** If the size x depth interaction is null across the whole
# ladder and in the upper segment, the confound story is wrong for this task
# and the experiment stops here.

# %%
# Fail loudly rather than repairing the environment with whatever PyPI
# serves today. Run `pip install -r requirements.txt` first.
try:
    import numpy as np
    import pandas as pd
    import statsmodels.formula.api as smf
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise RuntimeError(
        "Missing dependencies. Run: pip install -r requirements.txt"
    ) from exc

import os
from pathlib import Path

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:          # notebook cell, no __file__
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "figures").mkdir(exist_ok=True)

# The input is pinned, not just recorded. ProteinGym's main branch moves;
# fetching from it and stamping a SHA afterwards documents intent, not what
# the code consumed.
PROTEINGYM_COMMIT = os.environ.get(
    "PROTEINGYM_COMMIT", "144fe22b07dfaeec2b366f2346203a9838a55b4c"
)
assert PROTEINGYM_COMMIT != "main", "pin a commit SHA, not a branch"
assert len(PROTEINGYM_COMMIT) == 40, "use the full 40-character SHA"

RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"
SCORES = f"{RAW}/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
REF = f"{RAW}/reference_files/DMS_substitutions.csv"

scores = pd.read_csv(SCORES)
ref = pd.read_csv(REF)
print(scores.shape, ref.shape)

# %% [markdown]
# Schemas drift between ProteinGym releases. Assert rather than assume.

# %%
LADDER = {
    "ESM2 (8M)": 8e6,
    "ESM2 (35M)": 35e6,
    "ESM2 (150M)": 150e6,
    "ESM2 (650M)": 650e6,
    "ESM2 (3B)": 3e9,
    "ESM2 (15B)": 15e9,
}
missing = [c for c in LADDER if c not in scores.columns]
assert not missing, f"score columns changed: {missing}"
for c in ["DMS_id", "MSA_Neff_L", "MSA_Neff_L_category", "UniProt_ID", "taxon", "seq_len"]:
    assert c in ref.columns, f"reference column changed: {c}"

# %%
df = (
    scores[["DMS ID"] + list(LADDER)]
    .melt("DMS ID", var_name="model", value_name="rho")
    .assign(params=lambda d: d["model"].map(LADDER))
    .merge(
        ref[["DMS_id", "MSA_Neff_L", "MSA_Neff_L_category", "UniProt_ID", "taxon", "seq_len"]],
        left_on="DMS ID", right_on="DMS_id",
    )
    .dropna(subset=["rho", "MSA_Neff_L"])
)
df["log_params"] = np.log10(df["params"])
df["log_depth"] = np.log10(df["MSA_Neff_L"].clip(lower=0.01))
# Centre so the main effects read as "at average depth" / "at average size".
df["lp"] = df["log_params"] - df["log_params"].mean()
df["ld"] = df["log_depth"] - df["log_depth"].mean()

print(f"{df['DMS ID'].nunique()} assays x {df['model'].nunique()} model sizes")

# %% [markdown]
# ## The scaling curve, split by depth

# %%
tab = df.pivot_table(index="model", columns="MSA_Neff_L_category", values="rho", aggfunc="mean")
tab["All"] = df.groupby("model")["rho"].mean()
tab = tab.reindex(list(LADDER))[["Low", "Medium", "High", "All"]]
print(tab.round(3))
tab.to_csv(ROOT / "results/esm2_ladder_by_depth.csv")

# %% [markdown]
# ## Regression
#
# Two fits. The full ladder answers "does scale help at all". The 650M-and-up
# segment answers the question people actually argue about.

# %%
def fit(data, label):
    m = smf.ols("rho ~ lp*ld", data).fit(
        cov_type="cluster", cov_kwds={"groups": data["UniProt_ID"]}
    )
    print(f"\n=== {label}  (n={len(data)}, R2={m.rsquared:.3f}) ===")
    print(m.summary().tables[1])
    return m

full = fit(df, "full ladder, 8M to 15B")
upper = fit(df[df["params"] >= 650e6], "upper segment, 650M to 15B")

# %% [markdown]
# ## Paired within-assay comparison
#
# Averages hide direction. Ask instead: on what fraction of assays does 15B
# actually beat 650M?

# %%
wide = df.pivot_table(index="DMS ID", columns="model", values="rho").join(
    ref.set_index("DMS_id")[["MSA_Neff_L_category"]]
)
wide["delta"] = wide["ESM2 (15B)"] - wide["ESM2 (650M)"]
paired = (
    wide.groupby("MSA_Neff_L_category")["delta"]
    .agg(mean="mean", median="median", n="count", frac_15B_wins=lambda x: (x > 0).mean())
    .reindex(["Low", "Medium", "High"])
)
print(paired.round(3))
paired.to_csv(ROOT / "results/esm2_15B_vs_650M_paired.csv")

# %% [markdown]
# ## Figure

# %%
fig, ax = plt.subplots(figsize=(6, 4))
for stratum, style in [("Low", "o--"), ("Medium", "s-"), ("High", "^-")]:
    ax.plot(np.log10(list(LADDER.values())), tab[stratum], style, label=f"{stratum} MSA depth")
ax.set_xlabel("log10 parameters")
ax.set_ylabel("mean Spearman (ProteinGym substitutions)")
ax.set_title("ESM-2 scaling by MSA depth")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(ROOT / "figures/esm2_scaling_by_depth.png", dpi=200)

# %% [markdown]
# ## Provenance
#
# Reads the same variable the fetch used, so the record cannot drift from
# the input.

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("01a_msa_depth", out_dir=ROOT / "results", proteingym_commit=PROTEINGYM_COMMIT)

# %% [markdown]
# ## Replication on an independent ladder
#
# The same file carries ProGen2 S/M/Base/L/XL — different architecture,
# different training objective, same assays. If the pattern reproduces there,
# it is not an ESM artefact. Run this before believing anything above.

# %%
PROGEN = {"Progen2 S": 151e6, "Progen2 M": 764e6, "Progen2 Base": 764e6,
          "Progen2 L": 2.7e9, "Progen2 XL": 6.4e9}
# TODO: confirm ProGen2 parameter counts against Nijkamp et al. 2023 before
# fitting -- the Base/M distinction is a training-corpus difference, not a
# size difference, so one of them must be dropped from the ladder.

# %% [markdown]
# ## Confounds still open
#
# MSA depth is not randomly assigned. Deep families skew toward well-studied
# human and viral proteins with particular assay types. Before this is a
# result rather than a pattern, add `taxon`, `seq_len`, and `selection_type`
# as covariates and check the interaction survives.
