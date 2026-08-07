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
# # 03 — Depth-balanced resampling (arm 3) and differential attenuation (arm 4)
#
# CPU only, a few minutes. The sharp test: does the depth-scaling interaction
# survive when the cross-family depth composition is balanced rather than taken
# as ProteinGym gives it? And is ProGen more composition-sensitive than ESM-2?
#
# Executes the frozen protocol; designs nothing. Run
# `03_cross_family_independence.py` first — it writes the independent-set
# manifest this notebook reads.
#
# **Implementation decisions not fixed by the protocol** (surfaced for review;
# flip here before running if the reviewer disagrees):
#   - *Common set* = deterministic representatives present in all three ladders
#     AND passing the 0.10 floor in all three. The protocol fixes a per-ladder
#     floor and requires presence in all three; requiring the floor to pass in
#     all three is the strict reading taken here, so every unit is a genuine
#     test assay for every model family.
#   - *Fixed denominator* `beta_original` is the confirmatory point estimate on
#     the common set with native depth weighting, per ladder — same base
#     population as the balanced draws, differing only in depth weighting, so
#     the ratio isolates the balancing manipulation. It is computed once and
#     held constant across all balanced replicates (never recomputed inside a
#     draw), which is what removes the `beta_original ~ 0` blow-up.
#   - *Balanced draw* is a stratified bootstrap: n_balance clusters drawn WITH
#     replacement from each depth stratum per replicate, so the min stratum
#     also carries resampling variance.

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

SEED = 0
REPS = 2000
FLOOR = 0.10
MIN_STRATUM = 10          # kill criterion: any stratum below this => depth-balance not run
LADDERS = ["esm2", "progen2", "progen3"]
STRATA = ["Low", "Medium", "High"]
PRIMARY = "rho ~ lp*ld + lp*taxon + lp*z_llen + lp*sel"

MANIFEST = ROOT / "results/02_proteingym_assay_manifest.csv"
INDEP_FILE = ROOT / "results/03_independent_set_manifest.csv"
assert INDEP_FILE.exists(), "run 03_cross_family_independence.py first"


# %%
def load_ladder(name):
    if name == "esm2":
        d = pd.read_csv(ROOT / "results/01b_merged_assay_data.csv")
    else:
        d = pd.read_csv(ROOT / f"results/01c_{name}_assay_data.csv")
        d = d.rename(columns={"selection_type": "sel"})
        d = d[d.upper_segment].copy()
    return d


def recentre(d):
    d = d.copy()
    for col, src in [("lp", "params"), ("ld", "MSA_Neff_L"), ("llen", "seq_len")]:
        v = np.log10(d[src].clip(lower=0.01))
        d[col] = v - v.mean()
    d["z_llen"] = d["llen"] / d["llen"].std()
    return d


def beta_ld(data):
    try:
        return smf.ols(PRIMARY, recentre(data)).fit().params["lp:ld"]
    except Exception:
        return np.nan


def weak_assays(d, thr=FLOOR):
    g = d.groupby("assay")["rho"].apply(lambda s: (s.abs() < thr).all())
    return set(g[g].index)


# %% [markdown]
# ## Common set — present in all three ladders and passing the floor in all three

# %%
man = pd.read_csv(MANIFEST)
reps = pd.read_csv(INDEP_FILE)
INDEP = set(reps.assay_id)
CAT = dict(zip(man.assay_id, man.msa_neff_l_category))

frames = {L: load_ladder(L) for L in LADDERS}
present = {L: set(frames[L].assay.unique()) for L in LADDERS}
weak = {L: weak_assays(frames[L]) for L in LADDERS}

common_on = {a for a in INDEP
             if all(a in present[L] and a not in weak[L] for L in LADDERS)}
print(f"common set (floor on): {len(common_on)} assays")

CL = dict(zip(man.assay_id, man.homology_cluster_50))
groups = {L: {a: g for a, g in frames[L].groupby("assay")} for L in LADDERS}


def pctl(x):
    x = np.asarray(x)
    return np.percentile(x, [2.5, 50, 97.5])


def build_common(apply_floor):
    ok = (lambda a, L: a not in weak[L]) if apply_floor else (lambda a, L: True)
    common = {a for a in INDEP if all(a in present[L] and ok(a, L) for L in LADDERS)}
    by_stratum = {s: sorted(a for a in common if CAT.get(a) == s) for s in STRATA}
    counts = {s: len(by_stratum[s]) for s in STRATA}
    return common, by_stratum, counts, min(counts.values())


# %% [markdown]
# ## Arms 3-4, run with and without the chance floor
#
# For each floor setting: build the common set, fix per-ladder native-depth
# denominators (`beta_original`, computed once and held constant), then draw
# n_balance clusters per depth stratum with replacement and refit under all
# three ladders on the same draw. `A = 1 - |beta_balanced| / |beta_original|`;
# `D = A_ProGen - A_ESM2`. The floor=off pass is what makes the "holds only
# with the floor" decision-table row evaluable.

# %%
def run_floor(apply_floor):
    tag = "on" if apply_floor else "off"
    common, by_stratum, counts, n_balance = build_common(apply_floor)
    print(f"\nfloor={tag}: common {len(common)} | strata {counts} | n_balance = {n_balance}")
    if n_balance < MIN_STRATUM:
        print(f"  KILL CRITERION: smallest stratum {n_balance} < {MIN_STRATUM}; "
              f"depth-balance NOT run for floor={tag}.")
        return [], [], None

    beta_orig = {}
    for L in LADDERS:
        di = recentre(frames[L][frames[L].assay.isin(common)].copy())
        di["cl"] = di.assay.map(CL)
        beta_orig[L] = smf.ols(PRIMARY, di).fit(
            cov_type="cluster", cov_kwds={"groups": di.cl}).params["lp:ld"]
        print(f"  {L:8s} beta_original = {beta_orig[L]:+.4f}")

    rng = np.random.default_rng(SEED)
    bal_ = {L: [] for L in LADDERS}
    Arec = {L: [] for L in LADDERS}
    Drec_ = {"progen2": [], "progen3": []}
    converged = 0
    for _ in range(REPS):
        drawn = []
        for s in STRATA:
            pool = by_stratum[s]
            drawn.extend(pool[i] for i in rng.integers(0, len(pool), size=n_balance))
        betas = {L: beta_ld(pd.concat([groups[L][a] for a in drawn], ignore_index=True))
                 for L in LADDERS}
        if any(np.isnan(v) for v in betas.values()):
            continue
        converged += 1
        A = {}
        for L in LADDERS:
            bal_[L].append(betas[L])
            A[L] = 1 - abs(betas[L]) / abs(beta_orig[L])
            Arec[L].append(A[L])
        Drec_["progen2"].append(A["progen2"] - A["esm2"])
        Drec_["progen3"].append(A["progen3"] - A["esm2"])
    print(f"  converged {converged}/{REPS}")

    d_rows, f_rows = [], []
    for L in LADDERS:
        blo, bmed, bhi = pctl(bal_[L]); alo, amed, ahi = pctl(Arec[L])
        includes_zero = blo <= 0 <= bhi
        sign_flip = np.sign(bmed) != np.sign(beta_orig[L])
        material = (amed >= 0.50) or sign_flip or includes_zero
        print(f"  {L:8s} beta_bal {bmed:+.4f} [{blo:+.4f},{bhi:+.4f}] | "
              f"A {amed:+.3f} [{alo:+.3f},{ahi:+.3f}] | material={material}")
        d_rows.append(dict(floor=tag, ladder=L, beta_original=beta_orig[L],
                           beta_balanced_median=bmed, beta_balanced_ci_lo=blo,
                           beta_balanced_ci_hi=bhi, A_median=amed, A_ci_lo=alo,
                           A_ci_hi=ahi, balanced_includes_zero=includes_zero,
                           sign_flip=sign_flip, material_attenuation=material))
    for L in ["progen2", "progen3"]:
        dlo, dmed, dhi = pctl(Drec_[L]); p_pos = float(np.mean(np.asarray(Drec_[L]) > 0))
        print(f"  D=A_{L}-A_esm2 {dmed:+.3f} [{dlo:+.3f},{dhi:+.3f}] P(D>0)={p_pos:.3f}")
        f_rows.append(dict(floor=tag, contrast=f"A_{L} - A_esm2", D_median=dmed,
                           D_ci_lo=dlo, D_ci_hi=dhi, p_D_positive=p_pos,
                           interval_excludes_zero=(dlo > 0 or dhi < 0)))
    return d_rows, f_rows, (bal_, beta_orig, Drec_)


depth_rows, diff_rows, plot_state = [], [], None
for apply_floor in (True, False):
    d_rows, f_rows, state = run_floor(apply_floor)
    depth_rows += d_rows
    diff_rows += f_rows
    if apply_floor:
        plot_state = state          # figures use the floor=on (primary) run

pd.DataFrame(depth_rows).to_csv(ROOT / "results/03_depth_balanced.csv", index=False)
pd.DataFrame(diff_rows).to_csv(ROOT / "results/03_differential_attenuation.csv", index=False)

if plot_state is None:
    raise SystemExit("floor=on run hit the kill criterion; no figures to draw.")
bal, beta_original, Drec = plot_state


# %% [markdown]
# ## Figures

# %%
fig, ax = plt.subplots(figsize=(7.5, 4.6))
colors = {"esm2": "#2b5797", "progen2": "#c0392b", "progen3": "#8e44ad"}
for L in LADDERS:
    ax.hist(bal[L], bins=40, alpha=0.5, color=colors[L], label=f"{L} balanced")
    ax.axvline(beta_original[L], color=colors[L], ls="--", lw=1.2)
ax.axvline(0, color="k", lw=0.9)
ax.set_xlabel("depth-balanced size x depth interaction (beta)")
ax.set_ylabel("replicates")
ax.set_title("03 arm 3: depth-balanced interaction by ladder\n"
             "(dashed = native-depth beta_original; solid black = 0)", fontsize=10)
ax.legend(frameon=False, fontsize=8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "figures/03_depth_balanced_distribution.png", dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.5, 4.0))
for L, c in [("progen2", "#c0392b"), ("progen3", "#8e44ad")]:
    ax.hist(Drec[L], bins=40, alpha=0.5, color=c, label=f"A_{L} - A_esm2")
ax.axvline(0, color="k", lw=0.9, label="D = 0 (no differential)")
ax.set_xlabel("differential attenuation D  (predicted > 0)")
ax.set_ylabel("replicates")
ax.set_title("03 arm 4: differential attenuation, ProGen vs ESM-2", fontsize=10)
ax.legend(frameon=False, fontsize=8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "figures/03_differential_attenuation.png", dpi=200)
plt.close(fig)
print("figures written")

# %% [markdown]
# ## Provenance

# %%
import sys

sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("03_depth_balanced", out_dir=ROOT / "results",
      common_set=len(common), strata=counts, n_balance=n_balance,
      floor=FLOOR, reps=REPS, seed=SEED, converged=converged,
      beta_original={L: float(beta_original[L]) for L in LADDERS})
print("\ndone — arms 3-4 written to results/03_*.csv and figures/03_*.png")
