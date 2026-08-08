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
# # 05 Stage 2 — architecture vs power for the ProGen within-assay null
#
# CPU only. Executes `experiments/05-architecture-vs-power/PROTOCOL.md`, frozen
# at commit 4f70a03. This notebook designs nothing; every specification below —
# the estimator, the ladders, the injection grid, the residual-source bracket,
# the decision rule — was fixed in the protocol before any 05 number was seen.
#
# The reframe the protocol rests on: ProteinGym already scored every released
# ProGen size on the same 217 assays. 01c kept only the upper segment. The lower
# checkpoints are unused within-assay leverage, recovered here at zero GPU cost.
#
# Three analyses:
#   A  restore the discarded leverage (extended + segmented FE fits)
#   B  injection MDE — did the upper design ever have power to see an ESM-sized effect
#   C  ESM-2 handicap — does a real effect survive ProGen-sized leverage (corroboration)

# %%
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
(ROOT / "figures").mkdir(exist_ok=True)
(ROOT / "results").mkdir(exist_ok=True)

SEED = 0
REPS = int(os.environ.get("REPS", "2000"))          # protocol: 2000; set REPS=200 for a smoke run
BETA_ESM2 = -0.01527576017255374                    # 01b confirmatory reference

# ProteinGym pin — identical to 01a. Only the ESM-2 full-ladder anchor (Analysis A)
# fetches from it; everything else runs on checked-in CSVs.
PROTEINGYM_COMMIT = os.environ.get(
    "PROTEINGYM_COMMIT", "144fe22b07dfaeec2b366f2346203a9838a55b4c"
)
assert len(PROTEINGYM_COMMIT) == 40 and PROTEINGYM_COMMIT != "main"

# Breakpoints (params) defining each model's upper segment, from 01c.
BREAKPOINT = {"esm2": 650e6, "progen2": 764e6, "progen3": 762e6}

# Injection magnitudes (Analysis B). Own-confirmatory values loaded below.
INJECT_GRID = np.round(np.arange(-0.040, -0.0039, 0.002), 4)   # sweep for the power curve
for extra in (0.0, BETA_ESM2):
    if not np.any(np.isclose(INJECT_GRID, extra)):
        INJECT_GRID = np.sort(np.append(INJECT_GRID, round(extra, 4)))[::-1]

# %% [markdown]
# ## Confirmatory reference magnitudes (ratio denominators)
#
# The ratio `beta_recovered / beta_confirmatory` uses each ProGen's own 01c
# between-assay estimate. Read it, do not hardcode.

# %%
def confirmatory_beta(ladder):
    c = pd.read_csv(ROOT / f"results/01c_{ladder}_model_coefficients.csv")
    row = c[c.check.str.startswith("confirmatory")]
    return float(row.beta.iloc[0])

BETA_CONF = {L: confirmatory_beta(L) for L in ("progen2", "progen3")}
print("reference betas:")
print(f"  ESM-2 within-assay (01b)   {BETA_ESM2:+.4f}")
for L in ("progen2", "progen3"):
    print(f"  {L} confirmatory (01c)     {BETA_CONF[L]:+.4f}")

# %% [markdown]
# ## Estimator
#
# The 01b/01c fixed-effects estimator, unchanged: `rho ~ lp + lp:ld + C(assay)`,
# cluster-robust SE by UniProt_ID. `fe_sm` is the statsmodels reference. `fe_fast`
# is an algebraically identical within-transform used only in the injection loop,
# where thousands of fits on a fixed design make dummy expansion wasteful. The
# gate cell asserts the two agree to machine precision before `fe_fast` is trusted.

# %%
def _center(frame, subset):
    d = subset.copy()
    for col, src in [("lp", "params"), ("ld", "MSA_Neff_L")]:
        v = np.log10(d[src].clip(lower=0.01))
        d[col] = v - v.mean()
    return d

def fe_sm(d):
    """statsmodels reference; returns (beta, se, p, n) for lp:ld."""
    m = smf.ols("rho ~ lp + lp:ld + C(assay)", d).fit(
        cov_type="cluster", cov_kwds={"groups": d.UniProt_ID})
    lo, hi = m.conf_int().loc["lp:ld"]
    return dict(beta=m.params["lp:ld"], se=m.bse["lp:ld"], p=m.pvalues["lp:ld"],
                ci_lo=lo, ci_hi=hi, n=int(m.nobs))

class FEDesign:
    """Precompute the fixed within-design for a ladder subset so the injection
    loop only swaps y. statsmodels cluster cov uses the normal dist (use_t=False,
    verified), so p comes from 2*norm.sf(|z|)."""
    def __init__(self, d):
        self._sort = np.argsort(d.assay.values, kind="stable")   # apply to y in fit()
        d = d.iloc[self._sort].reset_index(drop=True)
        self.assay = d.assay.values
        self.acode = pd.factorize(self.assay)[0]
        self.gcode = pd.factorize(d.UniProt_ID.values)[0]
        self.N = len(d); self.n_assay = d.assay.nunique(); self.G = d.UniProt_ID.nunique()
        self.K = 2 + self.n_assay                       # absorbed params: lp, lp:ld, assay dummies (incl intercept)
        X = np.column_stack([d.lp.values, (d.lp * d.ld).values])
        self.Xd = self._demean(X)
        self.XtXinv = np.linalg.inv(self.Xd.T @ self.Xd)
        self.c = (self.G / (self.G - 1)) * ((self.N - 1) / (self.N - self.K))
        self.lp_ld = (d.lp * d.ld).values              # centred interaction regressor (for injection)
    def _demean(self, M):
        M = np.atleast_2d(M.T).T if M.ndim == 1 else M
        out = np.empty_like(M, dtype=float)
        for j in range(M.shape[1]):
            s = pd.Series(M[:, j]).groupby(self.acode).transform("mean").values
            out[:, j] = M[:, j] - s
        return out
    def fit(self, y, use_c=True):
        y = np.asarray(y)[self._sort]
        yd = self._demean(y.reshape(-1, 1))[:, 0]
        beta = self.XtXinv @ (self.Xd.T @ yd)
        u = yd - self.Xd @ beta
        contrib = self.Xd * u[:, None]
        S = np.zeros((self.G, 2))
        np.add.at(S, self.gcode, contrib)
        meat = S.T @ S
        c = self.c if use_c else 1.0        # use_c=True matches the real test; False isolates DGP validity
        V = c * self.XtXinv @ meat @ self.XtXinv
        se = np.sqrt(V[1, 1]); b = beta[1]; z = b / se
        return b, se, 2 * stats.norm.sf(abs(z))

# %% [markdown]
# ## Gate — fe_fast must equal statsmodels on the real ladders
#
# Hard stop if it does not. A power sim on a mis-specified SE is worse than none.

# %%
def _prog(L):
    return pd.read_csv(ROOT / f"results/01c_{L}_assay_data.csv").rename(
        columns={"selection_type": "sel"})

for L in ("progen2", "progen3"):
    for up in (True, False):
        raw = _prog(L)
        sub = raw[raw.upper_segment] if up else raw
        d = _center(raw, sub)
        ref = fe_sm(d)
        b, se, p = FEDesign(d).fit(d.rho.values)
        assert abs(b - ref["beta"]) < 1e-10, (L, up, b, ref["beta"])
        assert abs(se - ref["se"]) / ref["se"] < 1e-9, (L, up, se, ref["se"])
        assert abs(p - ref["p"]) < 1e-9, (L, up, p, ref["p"])
print("gate passed: fe_fast == statsmodels to machine precision on all ProGen ladders")

# %% [markdown]
# ## Data loaders
#
# ProGen full ladders are checked in. ESM-2 upper segment (650M/3B/15B) is in the
# 01b merged frame. The ESM-2 **full** ladder (for the Analysis A anchor) is not
# stored per-assay — reconstruct it from ProteinGym at the pin, exactly as 01a
# built it. If the fetch fails (offline), the anchor is recorded deferred and the
# rest of Stage 2 proceeds, per the protocol's kill-criterion discipline.

# %%
def esm2_upper():
    d = pd.read_csv(ROOT / "results/01b_merged_assay_data.csv")
    d = d.rename(columns={"MSA_Neff_L": "MSA_Neff_L"})
    return d  # already has params, MSA_Neff_L, rho, assay, UniProt_ID (650M/3B/15B)

def esm2_full():
    RAW = f"https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/{PROTEINGYM_COMMIT}"
    SCORES = f"{RAW}/benchmarks/DMS_zero_shot/substitutions/Spearman/DMS_substitutions_Spearman_DMS_level.csv"
    REF = f"{RAW}/reference_files/DMS_substitutions.csv"
    LADDER = {"ESM2 (8M)": 8e6, "ESM2 (35M)": 35e6, "ESM2 (150M)": 150e6,
              "ESM2 (650M)": 650e6, "ESM2 (3B)": 3e9, "ESM2 (15B)": 15e9}
    scores = pd.read_csv(SCORES); ref = pd.read_csv(REF)
    df = (scores[["DMS ID"] + list(LADDER)]
          .melt("DMS ID", var_name="model", value_name="rho")
          .assign(params=lambda d: d["model"].map(LADDER))
          .merge(ref[["DMS_id", "MSA_Neff_L", "UniProt_ID", "taxon", "seq_len"]],
                 left_on="DMS ID", right_on="DMS_id")
          .dropna(subset=["rho", "MSA_Neff_L"]))
    df["assay"] = df["DMS_id"]
    return df

ESM2_FULL_OK = True
try:
    _esm_full = esm2_full()
    print(f"ESM-2 full ladder reconstructed: {_esm_full.assay.nunique()} assays x "
          f"{_esm_full.params.nunique()} sizes")
except Exception as exc:
    ESM2_FULL_OK = False
    _esm_full = None
    print(f"ESM-2 full ladder NOT reconstructed offline ({type(exc).__name__}); "
          "Analysis A ESM-2 anchor deferred, recorded not-available. Rest of Stage 2 runs.")

# %% [markdown]
# ## Ladder geometry helper

# %%
def geometry(params):
    s = np.sort(np.unique(params))
    return len(s), float(np.log10(s[-1]) - np.log10(s[0]))

def fmt(p):
    return f"{p/1e6:.0f}M" if p < 1e9 else f"{p/1e9:.2f}B"

# %% [markdown]
# # Analysis A — restore the discarded leverage
#
# A1: lower the breakpoint using already-scored checkpoints; refit the single-slope
# FE interaction as the segment grows. A2: full-ladder global fit, plus the
# segmented `lp_lo:ld / lp_hi:ld` decomposition that localises the interaction so a
# recovered negative is credited to the regime that carries it — the guard against
# a low-regime signal masquerading as recovery.

# %%
def cluster_bootstrap_beta(frame_subset, statistic, reps=REPS, seed=SEED):
    """Resample whole assay clusters; return percentile CI of `statistic(centered_df)`."""
    rng = np.random.default_rng(seed)
    d = frame_subset
    assays = d.assay.unique()
    groups = {a: g for a, g in d.groupby("assay")}
    out = []
    for _ in range(reps):
        pick = rng.choice(assays, size=len(assays), replace=True)
        bs = pd.concat([groups[a] for a in pick], ignore_index=True)
        try:
            out.append(statistic(_center(bs, bs)))
        except Exception:
            pass
    out = np.array(out)
    return np.percentile(out, [2.5, 50, 97.5]), out

def global_beta(centered):
    return FEDesign(centered).fit(centered.rho.values)[0]

def segmented_fit(raw_subset, breakpoint):
    """FE fit with hinge terms at the breakpoint; returns upper- and lower-regime
    depth-scaling slopes. Uses statsmodels (4 regressors)."""
    d = _center(raw_subset, raw_subset)
    lp_raw = np.log10(d.params.clip(lower=0.01)); bp = np.log10(breakpoint)
    d["lp_lo"] = np.minimum(lp_raw, bp) - np.minimum(lp_raw, bp).mean()
    d["lp_hi"] = np.maximum(lp_raw - bp, 0.0); d["lp_hi"] = d["lp_hi"] - d["lp_hi"].mean()
    m = smf.ols("rho ~ lp_lo + lp_hi + lp_lo:ld + lp_hi:ld + C(assay)", d).fit(
        cov_type="cluster", cov_kwds={"groups": d.UniProt_ID})
    def grab(term):
        lo, hi = m.conf_int().loc[term]
        return dict(beta=m.params[term], p=m.pvalues[term], ci_lo=lo, ci_hi=hi)
    return {"lp_hi:ld": grab("lp_hi:ld"), "lp_lo:ld": grab("lp_lo:ld"), "n": int(m.nobs)}

def analysis_A_rows():
    rows = []
    ladders = [("progen2", _prog("progen2"), BREAKPOINT["progen2"]),
               ("progen3", _prog("progen3"), BREAKPOINT["progen3"])]
    if ESM2_FULL_OK:
        ladders.append(("esm2", _esm_full, BREAKPOINT["esm2"]))

    for name, raw, bp in ladders:
        sizes = np.sort(raw.params.unique())
        # A1 rungs: progressively include lower checkpoints, from upper segment down.
        upper = sizes[sizes >= bp]
        lowers = sizes[sizes < bp][::-1]                 # add nearest-below first
        rung_sets, cur = [upper], list(upper)
        for s in lowers:
            cur = sorted(set(cur) | {s})
            rung_sets.append(np.array(cur))
        for rung in rung_sets:
            sub = raw[raw.params.isin(rung)]
            cen = _center(sub, sub)
            npts, span = geometry(rung)
            b, se, p = FEDesign(cen).fit(cen.rho.values)
            (blo, bmed, bhi), _ = cluster_bootstrap_beta(sub, global_beta)
            rows.append(dict(
                model=name, analysis="A1_extent", segment=f"{fmt(rung.min())}+",
                term="lp:ld", pts=npts, log_span=round(span, 3),
                beta=b, se=se, p=p, boot_lo=blo, boot_med=bmed, boot_hi=bhi,
                ratio_esm2=b / BETA_ESM2,
                ratio_conf=(b / BETA_CONF[name]) if name in BETA_CONF else np.nan))
        # A2 segmented on the full ladder.
        seg = segmented_fit(raw, bp)
        npts, span = geometry(sizes)
        for term in ("lp_hi:ld", "lp_lo:ld"):
            g = seg[term]
            rows.append(dict(
                model=name, analysis="A2_segmented", segment="full",
                term=term, pts=npts, log_span=round(span, 3),
                beta=g["beta"], se=np.nan, p=g["p"], boot_lo=g["ci_lo"],
                boot_med=g["beta"], boot_hi=g["ci_hi"],
                ratio_esm2=g["beta"] / BETA_ESM2,
                ratio_conf=(g["beta"] / BETA_CONF[name]) if name in BETA_CONF else np.nan))
    return pd.DataFrame(rows)

A = analysis_A_rows()
A.to_csv(ROOT / "results/05_ladder_extent_sensitivity.csv", index=False)
print("\nAnalysis A — within-assay lp:ld as leverage is restored")
for name in A.model.unique():
    sub = A[A.model == name]
    print(f"\n  {name}:")
    for _, r in sub.iterrows():
        star = "*" if r.p < 0.05 else " "
        print(f"    {r.analysis:<12s} {r.segment:<6s} {r.term:<8s} pts={r.pts} "
              f"span={r.log_span:.2f}  beta={r.beta:+.4f}{star} p={r.p:.3f} "
              f"boot[{r.boot_lo:+.4f},{r.boot_hi:+.4f}]  b/ESM={r.ratio_esm2:+.2f}")

# %% [markdown]
# # Analysis B — injection MDE with a residual-source bracket
#
# On each ProGen ladder's actual geometry, inject a known within-assay slope onto a
# no-interaction base (assay means + real lp main effect), draw block residuals by
# assay, refit, count detections. The noise source is bracketed: CONTAM
# (`rho~lp+C(assay)` residuals, real lp:ld left in the noise → conservative) and
# CLEAN (`rho~lp+lp:ld+C(assay)` residuals, pure noise → power ceiling and the
# true-null calibration). Output is the minimum detectable effect, not a verdict.

# %%
def base_and_noise(centered, noise_spec):
    """Injection scaffold. `base` = fitted from the NO-interaction model
    (assay means + real lp main effect), so beta_inject alone sets the total
    within-assay interaction and inject=0 is a true null. `noise` = residuals
    from `noise_spec`, block-resampled downstream. Both models fit lp, so neither
    noise source reintroduces the lp main effect that base already carries."""
    base = smf.ols("rho ~ lp + C(assay)", centered).fit().fittedvalues.values
    noise = smf.ols(noise_spec, centered).fit().resid.values
    return base, noise

def assay_blocks(centered):
    """Return list of row-index arrays, one per assay, each ordered by params."""
    d = centered.sort_values(["assay", "params"])
    order = d.index.values
    blocks = [g.index.values for _, g in d.groupby("assay")]
    return blocks

def injection_power(name, segment, spec, beta_inject, design, fitted, resid,
                    blocks, reps=REPS, seed=SEED):
    rng = np.random.default_rng(seed)
    inj = beta_inject * design.lp_ld
    hits = 0
    nblk = len(blocks)
    resid_by_block = [resid[b] for b in blocks]
    idx_by_block = blocks
    for _ in range(reps):
        pick = rng.integers(0, nblk, size=nblk)
        ysyn = np.empty_like(fitted)
        for tgt, src in zip(range(nblk), pick):
            ysyn[idx_by_block[tgt]] = (fitted[idx_by_block[tgt]]
                                       + inj[idx_by_block[tgt]]
                                       + resid_by_block[src])
        b, se, p = design.fit(ysyn)
        if b < 0 and p < 0.05:
            hits += 1
    return hits / reps

def mde(grid, powers, target=0.8):
    """Smallest |beta| reaching target power, linear-interpolated on the sweep."""
    g = np.abs(grid); order = np.argsort(g)
    g, pw = g[order], np.array(powers)[order]
    above = np.where(pw >= target)[0]
    if len(above) == 0:
        return np.nan
    i = above[0]
    if i == 0:
        return -g[0]
    x0, x1, y0, y1 = g[i - 1], g[i], pw[i - 1], pw[i]
    x = x0 + (target - y0) * (x1 - x0) / (y1 - y0) if y1 != y0 else g[i]
    return -x

def calibrate(design, base, noise, blocks, reps=REPS, seed=SEED):
    """Inject 0 on the clean-null construction and count one-sided detections two
    ways: use_c=False isolates DGP validity (expect ~0.025); use_c=True is the real
    test's own type-I on this design (conservative, because the cluster small-sample
    factor over-corrects inside a resampling DGP)."""
    rng = np.random.default_rng(seed)
    nblk = len(blocks)
    noise_by_block = [noise[b] for b in blocks]
    hd = ht = 0
    for _ in range(reps):
        pick = rng.integers(0, nblk, size=nblk)
        ysyn = np.empty_like(base)
        for tgt, src in zip(range(nblk), pick):
            ysyn[blocks[tgt]] = base[blocks[tgt]] + noise_by_block[src]
        b, _, p = design.fit(ysyn, use_c=False)
        if b < 0 and p < 0.05:
            hd += 1
        b2, _, p2 = design.fit(ysyn, use_c=True)
        if b2 < 0 and p2 < 0.05:
            ht += 1
    return hd / reps, ht / reps

# Noise bracket over a shared no-interaction base:
#   CONTAM retains ProGen's real lp:ld in the resampled noise (conservative — real
#          interaction variance inflates the estimator, lowers power);
#   CLEAN  removes it (pure noise → power ceiling, and the true-null calibration).
CONTAM, CLEAN = "rho ~ lp + C(assay)", "rho ~ lp + lp:ld + C(assay)"
SOURCES = [CONTAM, CLEAN]
B_rows, calib = [], []
for name in ("progen2", "progen3"):
    raw = _prog(name)
    for segment in ("upper", "full"):
        sub = raw[raw.upper_segment] if segment == "upper" else raw
        cen = _center(raw, sub)
        cen_sorted = cen.sort_values("assay").reset_index(drop=True)
        design = FEDesign(cen_sorted)
        blocks = assay_blocks(cen_sorted)
        grid = np.unique(np.append(INJECT_GRID, round(BETA_CONF[name], 4)))[::-1]
        for spec in SOURCES:
            base, noise = base_and_noise(cen_sorted, spec)
            powers = []
            for bi in grid:
                pw = injection_power(name, segment, spec, bi, design, base, noise, blocks)
                powers.append(pw)
                B_rows.append(dict(model=name, segment=segment, source=spec,
                                   beta_inject=round(bi, 4), power=pw))
            this_mde = mde(grid, powers)
            B_rows.append(dict(model=name, segment=segment, source=spec,
                               beta_inject=np.nan, power=np.nan, mde=this_mde))
            tag = "conservative" if spec == CONTAM else "clean/ceiling"
            print(f"  {name} {segment:<5s} [{tag:<12s}] MDE(80%)={this_mde:+.4f}")
            if spec == CLEAN:
                fp_dgp, fp_test = calibrate(design, base, noise, blocks)
                calib.append(dict(model=name, segment=segment, fp_dgp=fp_dgp, fp_test=fp_test))

Bdf = pd.DataFrame(B_rows)
Bdf.to_csv(ROOT / "results/05_injection_power.csv", index=False)

# Calibration kill criterion. Detection is one-sided (beta<0 AND p<0.05), so the
# expected null rate is 0.025. The DGP-validity check uses the un-inflated cluster
# SE (use_c=False): the small-sample factor is correct for one-shot inference on
# data with the assay means absorbed, but inside the resampling DGP the empirical
# spread already carries that finite-sample variability, so applying it double-
# counts and makes the *real test* conservative. fp_dgp validates the null (band
# [0.015, 0.040]); fp_test is the real test's own type-I on this design — expected
# below 0.025 and reported, not a failure. A liberal fp_dgp (>0.040) discards B.
bad = [c for c in calib if not (0.015 <= c["fp_dgp"] <= 0.040)]
print("\nCalibration (inject=0, clean-null construction):")
print("  DGP validity (un-inflated SE, target ~0.025) | real-test type-I (matches analysis)")
for c in calib:
    flag = "" if 0.015 <= c["fp_dgp"] <= 0.040 else "  <-- DGP OUT OF RANGE"
    print(f"  {c['model']} {c['segment']:<5s}  fp_dgp={c['fp_dgp']:.3f}   "
          f"fp_test={c['fp_test']:.3f}{flag}")
if bad:
    print("WARNING: DGP validity outside band at REPS>=2000 signals a mis-specified "
          "simulator; discard B there (protocol kill criterion).")
pd.DataFrame(calib).to_csv(ROOT / "results/05_calibration.csv", index=False)

# %% [markdown]
# # Analysis C — ESM-2 handicap (corroboration only)
#
# Does a confirmed real effect survive ProGen-sized leverage? Two- and three-point
# ESM-2 sub-ladders on the upper segment, each bootstrapped so instability is shown,
# not asserted. C breaks ties between A and B; it never overrides them.

# %%
E = esm2_upper()
sub_ladders = [
    ("650M/3B/15B", [650e6, 3e9, 15e9]),
    ("650M->3B",    [650e6, 3e9]),
    ("3B->15B",     [3e9, 15e9]),
    ("650M->15B",   [650e6, 15e9]),
]
C_rows = []
for label, sizes in sub_ladders:
    sub = E[E.params.isin(sizes)]
    cen = _center(sub, sub)
    b, se, p = FEDesign(cen).fit(cen.rho.values)
    npts, span = geometry(np.array(sizes))
    (blo, bmed, bhi), _ = cluster_bootstrap_beta(sub, global_beta)
    C_rows.append(dict(subladder=label, pts=npts, log_span=round(span, 3),
                       beta=b, se=se, p=p, boot_lo=blo, boot_med=bmed, boot_hi=bhi,
                       ratio_esm2=b / BETA_ESM2))
    star = "*" if p < 0.05 else " "
    print(f"  {label:<14s} pts={npts} span={span:.2f}  beta={b:+.4f}{star} "
          f"p={p:.3f}  boot[{blo:+.4f},{bhi:+.4f}]")
pd.DataFrame(C_rows).to_csv(ROOT / "results/05_esm2_handicap.csv", index=False)

# %% [markdown]
# # Figures

# %%
fig, ax = plt.subplots(figsize=(7.6, 5.2))
colors = {"progen2": "#2b5797", "progen3": "#c0392b", "esm2": "#2e8b57"}
for name in A.model.unique():
    sub = A[(A.model == name) & (A.analysis == "A1_extent")].sort_values("log_span")
    ax.plot(sub.log_span, sub.beta, "-o", color=colors.get(name, "#555"), label=name, ms=5)
    for _, r in sub.iterrows():
        ax.plot([r.log_span, r.log_span], [r.boot_lo, r.boot_hi],
                color=colors.get(name, "#555"), lw=1, alpha=0.6)
ax.axhline(0, color="k", lw=0.9)
ax.axhline(BETA_ESM2, color="#888", ls="--", lw=1, label="ESM-2 within-assay")
ax.set_xlabel("ladder log10 span (decades)"); ax.set_ylabel("within-assay lp:ld (beta)")
ax.set_title("05A: within-assay interaction as discarded leverage is restored", fontsize=10)
ax.legend(frameon=False, fontsize=8)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT / "figures/05_lp_ld_by_ladder_length.png", dpi=200)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.6, 5.2))
prim = Bdf[(Bdf.source == CONTAM) & Bdf.beta_inject.notna()]
for (name, seg), g in prim.groupby(["model", "segment"]):
    g = g.sort_values("beta_inject")
    ls = "-" if seg == "full" else "--"
    ax.plot(-g.beta_inject, g.power, ls, label=f"{name} {seg}")
ax.axhline(0.8, color="#888", ls=":", lw=1); ax.axvline(-BETA_ESM2, color="#c0392b", ls="--", lw=1, label="|ESM-2 beta|")
ax.set_xlabel("injected |beta|"); ax.set_ylabel("detection power")
ax.set_title("05B: power to detect an injected within-assay effect (conservative source)", fontsize=10)
ax.legend(frameon=False, fontsize=8, loc="lower right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(ROOT / "figures/05_injection_power_curve.png", dpi=200)
plt.close(fig)
print("figures written")

# %% [markdown]
# ## Provenance

# %%
import sys
sys.path.insert(0, str(ROOT))
from src.provenance import stamp  # noqa: E402

stamp("05_architecture_vs_power", out_dir=ROOT / "results",
      proteingym_commit=PROTEINGYM_COMMIT, seed=SEED, reps=REPS,
      esm2_full_anchor="reconstructed" if ESM2_FULL_OK else "deferred_offline",
      noise_sources={"conservative": CONTAM, "clean": CLEAN},
      calibration="dual: un-inflated SE for DGP validity, corrected SE for real-test type-I",
      beta_esm2=BETA_ESM2, beta_conf=BETA_CONF)
print("done")
