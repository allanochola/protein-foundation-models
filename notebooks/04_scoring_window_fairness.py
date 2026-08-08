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
# # 04 (Phase 1) — Does the ESM-2 depth x scaling interaction depend on
# # assays exposed to long-sequence windowing?
#
# No GPU. Metadata-only falsification. See PROTOCOL.md (frozen before this ran).
#
# The confirmatory test reuses the Experiment 01 estimator verbatim
# (`rho ~ lp + lp:ld + C(assay)`, cluster-robust SE by UniProt_ID, upper
# segment 650M/3B/15B) and changes only the assay set: full 217 vs the
# context-safe subset (target protein <= 1022 residues, so ProteinGym's
# 1024-token window never fires at any checkpoint).
#
# Preregistered decision (PROTOCOL.md):
#   reference beta_full = -0.0153
#   Pass:    same sign and |beta_safe / beta_full| >= 0.75
#   Partial: 0.50 - 0.75
#   Fail:    < 0.50, ~0, or sign reversal
#   Significance reported, not the principal threshold.

# %%
import os, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "results").mkdir(exist_ok=True)

PROTEINGYM_COMMIT = os.environ.get(
    "PROTEINGYM_COMMIT", "144fe22b07dfaeec2b366f2346203a9838a55b4c"
)
assert len(PROTEINGYM_COMMIT) == 40, "use the full 40-character SHA"
RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"
SCORES = f"{RAW}/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
REF = f"{RAW}/reference_files/DMS_substitutions.csv"

WINDOW_RESIDUES = 1022  # get_optimal_window: full sequence kept when len <= 1022

# %% [markdown]
# ## Build the Experiment 01 dataset (identical construction to 01b)

# %%
LADDER = {
    "ESM2 (8M)": 8e6, "ESM2 (35M)": 35e6, "ESM2 (150M)": 150e6,
    "ESM2 (650M)": 650e6, "ESM2 (3B)": 3e9, "ESM2 (15B)": 15e9,
}
scores, ref = pd.read_csv(SCORES), pd.read_csv(REF)

df = (
    scores[["DMS ID"] + list(LADDER)]
    .melt("DMS ID", var_name="model", value_name="rho")
    .assign(params=lambda d: d["model"].map(LADDER))
    .merge(
        ref[["DMS_id", "MSA_Neff_L", "MSA_Neff_L_category", "UniProt_ID",
             "taxon", "seq_len", "coarse_selection_type"]],
        left_on="DMS ID", right_on="DMS_id",
    )
    .dropna(subset=["rho", "MSA_Neff_L"])
    .rename(columns={"DMS ID": "assay", "coarse_selection_type": "sel"})
)
for col, src in [("lp", "params"), ("ld", "MSA_Neff_L"), ("llen", "seq_len")]:
    v = np.log10(df[src].clip(lower=0.01))
    df[col] = v - v.mean()

# exposure + context-safe flag, derived from the pinned rule
df["exposure"] = np.maximum(0, df["seq_len"] - WINDOW_RESIDUES) / df["seq_len"]
df["windowed"] = df["seq_len"] > WINDOW_RESIDUES
df["context_safe"] = ~df["windowed"]

upper = df[df["params"] >= 650e6].copy()

n_assays = upper["assay"].nunique()
n_safe = upper.loc[upper.context_safe, "assay"].nunique()
n_win = upper.loc[upper.windowed, "assay"].nunique()
print(f"upper-segment assays: {n_assays}  | context-safe: {n_safe}  | windowed: {n_win}")
print(f"windowed distinct proteins: {upper.loc[upper.windowed,'UniProt_ID'].nunique()}")

# %% [markdown]
# ## Estimator (verbatim from 01b, Model 4)

# %%
def fit_model4(data, label):
    m = smf.ols("rho ~ lp + lp:ld + C(assay)", data).fit(
        cov_type="cluster", cov_kwds={"groups": data["UniProt_ID"]}
    )
    b, p = m.params["lp:ld"], m.pvalues["lp:ld"]
    ci = m.conf_int().loc["lp:ld"]
    print(f"{label:<34s} beta={b:+.4f}  p={p:.4f}  "
          f"[{ci[0]:+.4f}, {ci[1]:+.4f}]  n={int(m.nobs)}  "
          f"assays={data['assay'].nunique()}  clusters={data['UniProt_ID'].nunique()}")
    return dict(spec=label, beta=float(b), p=float(p),
                ci_lo=float(ci[0]), ci_hi=float(ci[1]),
                n=int(m.nobs), assays=int(data["assay"].nunique()),
                clusters=int(data["UniProt_ID"].nunique()))

rows = []
print("\n=== confirmatory (Model 4, upper segment) ===")
rows.append(fit_model4(upper, "full (217)"))                              # sanity: recover -0.0153
rows.append(fit_model4(upper[upper.context_safe], "context-safe (<=1022)"))

beta_full = rows[0]["beta"]
beta_safe = rows[1]["beta"]
ratio = abs(beta_safe / beta_full)
same_sign = np.sign(beta_safe) == np.sign(beta_full)

# %% [markdown]
# ## Confound term and exploratory windowed-only fit

# %%
print("\n=== confound term: does the interaction differ in windowed assays? ===")
cf = smf.ols("rho ~ lp + lp:ld + lp:windowed + lp:ld:windowed + C(assay)", upper).fit(
    cov_type="cluster", cov_kwds={"groups": upper["UniProt_ID"]}
)
term = "lp:ld:windowed[T.True]" if "lp:ld:windowed[T.True]" in cf.params.index else \
       [k for k in cf.params.index if "ld:windowed" in k][0]
print(f"lp:ld (context-safe baseline)  beta={cf.params['lp:ld']:+.4f}  p={cf.pvalues['lp:ld']:.4f}")
print(f"{term} (windowed shift)  beta={cf.params[term]:+.4f}  p={cf.pvalues[term]:.4f}")

print("\n=== exploratory (windowed assays only; restricted depth support, diagnostic) ===")
expl = fit_model4(upper[upper.windowed], "windowed-only (16)")

# %% [markdown]
# ## Decision

# %%
if same_sign and ratio >= 0.75:
    verdict = "PASS - scoring-protocol (windowing) explanation strongly weakened"
elif same_sign and ratio >= 0.50:
    verdict = "PARTIAL - possible partial contribution"
else:
    verdict = "FAIL - windowing artefact remains live; Stage 3 required"

print(f"\nbeta_full={beta_full:+.4f}  beta_safe={beta_safe:+.4f}  "
      f"ratio={ratio:.3f}  same_sign={same_sign}")
print("VERDICT:", verdict)

# %% [markdown]
# ## Outputs

# %%
man_cols = ["assay", "UniProt_ID", "seq_len", "exposure", "context_safe",
            "MSA_Neff_L", "MSA_Neff_L_category", "sel"]
manifest = (upper[man_cols].drop_duplicates("assay").sort_values("seq_len", ascending=False))
manifest.to_csv(ROOT / "results/04_exposure_manifest.csv", index=False)

interaction = pd.DataFrame(rows + [dict(
    spec="confound: lp:ld:windowed", beta=float(cf.params[term]),
    p=float(cf.pvalues[term]), ci_lo=np.nan, ci_hi=np.nan,
    n=int(cf.nobs), assays=n_assays, clusters=upper["UniProt_ID"].nunique()),
    dict(spec="exploratory: windowed-only", **{k: expl[k] for k in
         ["beta", "p", "ci_lo", "ci_hi", "n", "assays", "clusters"]})])
interaction.to_csv(ROOT / "results/04_interaction_by_exposure.csv", index=False)

json.dump(
    dict(proteingym_commit=PROTEINGYM_COMMIT, window_residues=WINDOW_RESIDUES,
         estimator="rho ~ lp + lp:ld + C(assay); cluster SE by UniProt_ID; upper 650M/3B/15B",
         beta_full=beta_full, beta_safe=beta_safe, ratio=ratio,
         same_sign=bool(same_sign), verdict=verdict,
         n_assays=int(n_assays), n_context_safe=int(n_safe), n_windowed=int(n_win),
         scoring_note="ESM-2 model_window=1024, size-invariant; long-sequence strategy "
                      "(wt-marginals+overlapping vs masked-marginals+optimal) still to be "
                      "pinned from launcher config before Stage 3"),
    open(ROOT / "results/provenance_04_scoring_window.json", "w"), indent=2)
print("\nwrote results/04_exposure_manifest.csv, results/04_interaction_by_exposure.csv, "
      "results/provenance_04_scoring_window.json")
