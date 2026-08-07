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
# # 01c Phase 2 — build the assay-level datasets
#
# **This notebook fits nothing.** It builds and commits the ProGen2 and
# ProGen3 datasets so that the ladder-construction decisions are locked into
# their own commit, before any analysis exists that could be tuned around
# them.
#
# The decision with real degrees of freedom here is dropping `Progen2 Base`.
# Separating construction from analysis makes that choice verifiable from the
# git history rather than asserted.
#
# Runs in seconds on CPU.

# %%
try:
    import numpy as np
    import pandas as pd
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
assert PROTEINGYM_COMMIT != "main" and len(PROTEINGYM_COMMIT) == 40
RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"

# %% [markdown]
# ## Ladders, fixed by protocol
#
# **ProGen2.** Nijkamp et al. Table 1: small 151M, medium 764M, base 764M,
# large 2.7B, xlarge 6.4B. Medium and base share parameter count, layers,
# heads and head dimension, differing in context length (1,024 vs 2,048) and
# training schedule. Two points at the same x with different y would bias the
# slope regardless of why they differ, so `Progen2 Base` is dropped.
#
# **ProGen3.** Verified against Profluent-AI/progen3. Sparse mixture of
# experts, ~27% of parameters active per forward pass; sparsity is constant
# across the ladder, so within-ladder log-scaling is well defined. ProteinGym
# does not score the 46B model, so the ladder is truncated at 3B.

# %%
LADDERS = {
    "progen2": {
        "Progen2 S": 151e6,
        "Progen2 M": 764e6,
        # "Progen2 Base": 764e6,   # dropped by protocol: duplicate x
        "Progen2 L": 2.7e9,
        "Progen2 XL": 6.4e9,
    },
    "progen3": {
        "Progen3 112m": 112e6,
        "Progen3 219m": 219e6,
        "Progen3 339m": 339e6,
        "Progen3 762m": 762e6,
        "Progen3 1B": 1.0e9,
        "Progen3 3B": 3.0e9,
    },
}
BREAKPOINT = {"progen2": 764e6, "progen3": 762e6}   # by analogy to ESM-2's 650M

# %%
scores = pd.read_csv(
    f"{RAW}/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
)
ref = pd.read_csv(f"{RAW}/reference_files/DMS_substitutions.csv")

META = ["DMS_id", "MSA_Neff_L", "MSA_Neff_L_category", "UniProt_ID",
        "taxon", "seq_len", "coarse_selection_type"]
for c in META:
    assert c in ref.columns, f"reference column changed: {c}"

# %%
for name, ladder in LADDERS.items():
    missing = [c for c in ladder if c not in scores.columns]
    assert not missing, f"{name}: score columns changed: {missing}"

    d = (
        scores[["DMS ID"] + list(ladder)]
        .melt("DMS ID", var_name="model", value_name="rho")
        .assign(params=lambda x: x["model"].map(ladder))
        .merge(ref[META], left_on="DMS ID", right_on="DMS_id")
        .dropna(subset=["rho", "MSA_Neff_L"])
        .rename(columns={"DMS ID": "assay", "coarse_selection_type": "selection_type"})
    )
    # Centring is done inside the analysis notebook, on the segment being
    # fitted. Storing raw values here keeps this file a dataset, not a
    # half-finished model matrix.
    d["log_params"] = np.log10(d["params"])
    d["log_msa_depth"] = np.log10(d["MSA_Neff_L"].clip(lower=0.01))
    d["log_seq_len"] = np.log10(d["seq_len"])
    d["upper_segment"] = d["params"] >= BREAKPOINT[name]

    cols = ["assay", "UniProt_ID", "model", "params", "log_params",
            "taxon", "selection_type", "seq_len", "log_seq_len",
            "MSA_Neff_L", "log_msa_depth", "MSA_Neff_L_category",
            "upper_segment", "rho"]
    out = ROOT / f"results/01c_{name}_assay_data.csv"
    d[cols].to_csv(out, index=False)

    span = np.log10(d[d.upper_segment].params.max()) - np.log10(d[d.upper_segment].params.min())
    print(f"{name:8s} {d.assay.nunique():3d} assays x {d.model.nunique()} sizes "
          f"= {len(d):4d} rows | upper segment {d.upper_segment.sum():3d} rows, "
          f"log10 span {span:.2f}")

# %% [markdown]
# ## Sanity checks
#
# Structural only. Nothing here inspects `rho` against `log_msa_depth` —
# that is the analysis, and it belongs in the next notebook and the next
# commit.

# %%
for name in LADDERS:
    d = pd.read_csv(ROOT / f"results/01c_{name}_assay_data.csv")
    assert d.groupby("assay").size().nunique() == 1, "unbalanced panel"
    assert d.params.nunique() == len(LADDERS[name]), "ladder size mismatch"
    assert d.groupby("params").assay.nunique().nunique() == 1, "assays differ across sizes"
    assert d.rho.notna().all() and d.log_msa_depth.notna().all()
print("all structural checks passed")

# %% [markdown]
# ## Shared-assay check
#
# Confirms the correlated-estimates problem recorded in Amendment 3: the
# ladders are scored on the same assays, so these are replications across
# architectures, not across benchmarks.

# %%
a2 = set(pd.read_csv(ROOT / "results/01c_progen2_assay_data.csv").assay)
a3 = set(pd.read_csv(ROOT / "results/01c_progen3_assay_data.csv").assay)
print(f"ProGen2 assays: {len(a2)} | ProGen3 assays: {len(a3)} | shared: {len(a2 & a3)}")

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("01c_dataset", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT,
      ladders={k: sorted(v.values()) for k, v in LADDERS.items()},
      dropped=["Progen2 Base (duplicate parameter count with Progen2 M)"])
