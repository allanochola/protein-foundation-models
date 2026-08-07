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
# # 02 — ProteinGym family-independence audit
#
# CPU only. Tier A (identity, taxon, depth) runs in seconds with no external
# tool. Tier B (homology clustering) needs mmseqs2 or CD-HIT; if neither is
# present it is skipped and marked "not run", per the protocol's kill
# criterion.
#
# This is inventory, not inference. It counts independent systems; it fits no
# model. See `experiments/02-proteingym-family-audit/PROTOCOL.md`.

# %%
try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise RuntimeError("Missing dependencies. Run: pip install -r requirements.txt") from exc

import os
import shutil
import subprocess
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

ref = pd.read_csv(f"{RAW}/reference_files/DMS_substitutions.csv")
print(f"{len(ref)} assays, {ref.UniProt_ID.nunique()} unique UniProt IDs")

# %% [markdown]
# ## Tier A — independence from the reference file alone
#
# Protein identity, taxon, depth, selection type. No tool required.

# %%
manifest = pd.DataFrame({
    "assay_id": ref.DMS_id,
    "uniprot_id": ref.UniProt_ID,
    "organism": ref.source_organism,
    "taxon": ref.taxon,
    "assay_type": ref.coarse_selection_type,
    "selection_type": ref.selection_type,
    "seq_len": ref.seq_len,
    "n_single_mutants": ref.DMS_number_single_mutants,
    "n_total_mutants": ref.DMS_total_number_mutants,
    "msa_neff_l": ref.MSA_Neff_L,
    "msa_neff_l_category": ref.MSA_Neff_L_category,
    "in_experiment_01": True,
})
# Tier B columns initialised NA; filled below if a tool is available.
for c in ["homology_cluster_90", "homology_cluster_70",
          "homology_cluster_50", "homology_cluster_30", "structural_family"]:
    manifest[c] = "NA"

print(f"217 assays -> {manifest.uniprot_id.nunique()} unique proteins "
      f"(protein-identity independence)")
print("\nassays per protein, distribution:")
print(manifest.uniprot_id.value_counts().value_counts().sort_index()
      .rename_axis("assays_per_protein").rename("n_proteins"))

# %% [markdown]
# ## Tier B — homology clustering
#
# Cluster the 217 target sequences at 90/70/50/30% identity. Prefer mmseqs2;
# fall back to CD-HIT; if neither is installed, mark "not run" and continue.

# %%
def find_tool():
    for t in ("mmseqs", "cd-hit"):
        if shutil.which(t):
            return t
    return None


def write_fasta(path):
    with open(path, "w") as fh:
        for _, r in ref.iterrows():
            fh.write(f">{r.DMS_id}\n{r.target_seq}\n")


def cluster_mmseqs(fasta, min_id, tmp):
    out = f"{tmp}/clu_{int(min_id*100)}"
    subprocess.run(
        ["mmseqs", "easy-cluster", fasta, out, f"{tmp}/tmp",
         "--min-seq-id", str(min_id), "-c", "0.8", "--cov-mode", "0", "-v", "0"],
        check=True, capture_output=True,
    )
    tsv = pd.read_csv(f"{out}_cluster.tsv", sep="\t", header=None,
                      names=["rep", "member"])
    return dict(zip(tsv.member, tsv.rep))


TOOL = find_tool()
THRESHOLDS = [0.90, 0.70, 0.50, 0.30]
tier_b_ran = False

if TOOL == "mmseqs":
    import tempfile
    tmp = tempfile.mkdtemp()
    fasta = f"{tmp}/pg.fasta"
    write_fasta(fasta)
    for thr in THRESHOLDS:
        m = cluster_mmseqs(fasta, thr, tmp)
        col = f"homology_cluster_{int(thr*100)}"
        manifest[col] = manifest.assay_id.map(m)
    tier_b_ran = True
    print("Tier B complete via mmseqs2")
elif TOOL == "cd-hit":
    print("cd-hit present; wire up per its output format (not implemented here) "
          "-- marking Tier B not run for this pass")
else:
    print("No clustering tool found. Install mmseqs2 in Colab with:")
    print("  !apt-get -qq install -y mmseqs2   ||   !pip install mmseqs2")
    print("Tier B marked 'not run'. Tier A results below are complete.")

# %% [markdown]
# ## Independence summary

# %%
rows = [{
    "threshold": "protein_identity_uniprot",
    "assays": len(manifest),
    "independent_clusters": manifest.uniprot_id.nunique(),
    "median_assays_per_cluster": manifest.uniprot_id.value_counts().median(),
    "largest_cluster_pct": round(100 * manifest.uniprot_id.value_counts().max()
                                 / len(manifest), 1),
}]
if tier_b_ran:
    for thr in THRESHOLDS:
        col = f"homology_cluster_{int(thr*100)}"
        vc = manifest[col].value_counts()
        rows.append({
            "threshold": f"{int(thr*100)}pct_identity",
            "assays": len(manifest),
            "independent_clusters": manifest[col].nunique(),
            "median_assays_per_cluster": vc.median(),
            "largest_cluster_pct": round(100 * vc.max() / len(manifest), 1),
        })

summary = pd.DataFrame(rows)
print(summary.to_string(index=False))

manifest.to_csv(ROOT / "results/02_proteingym_assay_manifest.csv", index=False)
summary.to_csv(ROOT / "results/02_independence_summary.csv", index=False)

# %% [markdown]
# ## Decision rules (evaluated only if Tier B ran)

# %%
if tier_b_ran:
    col50 = "homology_cluster_50"
    n_clusters = manifest[col50].nunique()
    vc50 = manifest[col50].value_counts()
    largest_pct = 100 * vc50.max() / len(manifest)
    # clusters with >=5 members per coarse selection class
    per_class = (manifest.groupby("assay_type")[col50].nunique())
    classes_ge5 = (per_class >= 5).sum()
    depth_spread = manifest.groupby(col50).msa_neff_l_category.first().nunique()

    print(f"independent clusters @50%: {n_clusters}  (rule 1: >=30 -> {n_clusters>=30})")
    print(f"largest cluster: {largest_pct:.1f}%  (rule 2: <=15% -> {largest_pct<=15})")
    print(f"selection classes with >=5 clusters: {classes_ge5}  (rule 3: >=2 -> {classes_ge5>=2})")
    print(f"depth categories represented: {depth_spread}  (rule 4: >1 -> {depth_spread>1})")
    go = (n_clusters >= 30) and (largest_pct <= 15) and (classes_ge5 >= 2) and (depth_spread > 1)
    print(f"\nExperiment 03 cross-family design warranted: {go}")
else:
    print("Tier B not run -- 50% decision rule cannot be evaluated.")
    print("Tier A stands: 217 assays -> "
          f"{manifest.uniprot_id.nunique()} unique proteins.")

# %% [markdown]
# ## Figures

# %%
fig, ax = plt.subplots(figsize=(6, 4))
vc = manifest.uniprot_id.value_counts().value_counts().sort_index()
ax.bar(vc.index, vc.values, color="#2b5797")
ax.set_xlabel("assays sharing one protein")
ax.set_ylabel("number of proteins")
ax.set_title(f"217 assays collapse to {manifest.uniprot_id.nunique()} proteins", fontsize=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "figures/02_assays_per_cluster.png", dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 4))
order = ["Low", "Medium", "High"]
counts = manifest.msa_neff_l_category.value_counts().reindex(order)
ax.bar(order, counts.values, color="#7f8c8d")
ax.set_xlabel("MSA depth category")
ax.set_ylabel("assays")
ax.set_title("Depth spread across the benchmark", fontsize=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "figures/02_depth_within_independent_set.png", dpi=200)
plt.close(fig)
print("figures written")

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("02_audit", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT,
      tier_b_tool=TOOL or "none", tier_b_ran=tier_b_ran)
