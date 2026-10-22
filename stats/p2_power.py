"""Paper 2: what the design could have detected, and how much the result depends on analysis choices.

1. Injection-recovery (detectability). Real case and control intervals are kept as measured. A synthetic tumour-driven
   change of rate delta (in multiples of the normative annual change) is added to every case interval for the part of
   the interval that falls within L years before diagnosis; intervals are re-standardised against the unchanged
   control reference and the pre-specified window comparison is repeated. Power = share of 1,000 matched-set bootstrap
   replicates, centred on the null by removing the observed estimate, in which the <1-year window shows significant
   acceleration (one-sided 2.5%); delta = 0 therefore returns the nominal 2.5%. This is a simulation for power; it adds known signal to real data and is labelled as such.
2. Specification curve. The <1-year (or 0.5/1.5-year) window estimate recomputed under 96 combinations of analysis
   choices: segmentation model x noise adjustment x covariate-adjusted reference x kernel filter x late-window cut-off
   x exclusion of scans within 90 days of diagnosis. 95% CIs from 500 matched-set bootstrap replicates.
3. Interval table for the raster figure (every case and control interval with its z).
Outputs: results/p2/power_grid.csv, power_type1.csv, spec_curve.csv, raster_intervals.csv
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "p2")
rng = np.random.default_rng(20261009)
long = pd.read_csv(os.path.join(RES, "partA_long.csv"), dtype={"pid": str, "set": str})
norm = pd.read_csv(os.path.join(RES, "partA_normative_model.csv"))
SIGN = {"psm_fat_pct": 1, "aut_fat_pct": 1, "psm_lvol": -1, "aut_lvol": -1, "psm_hu": -1, "aut_hu": -1}


def gamma_within(v):
    """noise coefficient from within-person variation in controls"""
    c = long[(long.group == "control") & long[v].notna()]
    y = c[v] - c.groupby("pid")[v].transform("mean")
    x = c.noise_c - c.groupby("pid").noise_c.transform("mean")
    return float((x * y).sum() / (x * x).sum())


def intervals(v, noise=True, gam=None, drop_gap90=False, kernel_consistent=False):
    d = long[long[v].notna()].copy()
    if kernel_consistent:
        d = d[d.kernel_consistent == 1]
    if drop_gap90:
        d = d[(d.group == "control") | (d.tau <= -90 / 365.25)]
    g = gam if gam is not None else gamma_within(v)
    d["y"] = d[v] - (g * d.noise_c if noise else 0)
    d = d.sort_values(["pid", "t"])
    nxt = d.groupby("pid").shift(-1)
    iv = d.assign(r1=nxt["round"], y1=nxt.y, t1=nxt.t, tau1=nxt.tau).dropna(subset=["r1"])
    iv = iv.assign(r0=iv["round"].astype(int), r1=iv.r1.astype(int), dt=iv.t1 - iv.t,
                   rate=(iv.y1 - iv.y) / (iv.t1 - iv.t), tau0=iv.tau, end_tau=iv.tau1)
    keep = d.groupby("pid").size()
    iv = iv[iv.pid.map(keep) >= 2]
    return iv[["pid", "group", "set", "age0", "sex", "smoke", "vendor0", "r0", "r1", "dt", "rate", "tau0", "end_tau"]]


def reference(iv, covars=True):
    """control reference per pair of rounds: (prediction function, residual SD)"""
    refs = {}
    for (a, b), c in iv[iv.group == "control"].groupby(["r0", "r1"]):
        if covars:
            X = pd.get_dummies(c[["age0", "sex", "smoke", "vendor0"]], columns=["sex", "smoke", "vendor0"],
                               drop_first=True, dtype=float)
            X.insert(0, "const", 1.0)
            beta, *_ = np.linalg.lstsq(X.values, c.rate.values, rcond=None)
            sd = np.std(c.rate.values - X.values @ beta, ddof=X.shape[1])
            refs[(a, b)] = (list(X.columns), beta, sd)
        else:
            refs[(a, b)] = (None, np.array([c.rate.mean()]), c.rate.std())
    return refs


def predict(iv, refs):
    pred, sd = np.full(len(iv), np.nan), np.full(len(iv), np.nan)
    for (a, b), (cols, beta, s) in refs.items():
        m = ((iv.r0 == a) & (iv.r1 == b)).values
        if not m.any():
            continue
        if cols is None:
            pred[m] = beta[0]
        else:
            X = pd.get_dummies(iv.loc[m, ["age0", "sex", "smoke", "vendor0"]], columns=["sex", "smoke", "vendor0"],
                               drop_first=True, dtype=float)
            X.insert(0, "const", 1.0)
            X = X.reindex(columns=cols, fill_value=0.0)
            pred[m] = X.values @ beta
        sd[m] = s
    return pred, sd


def window_stats(z, late, is_case, set_idx, nset, B, boot_sets=None):
    """late-window excess z (cases' late intervals minus controls) with matched-set bootstrap replicates"""
    est = z[is_case & late].mean() - z[~is_case].mean()
    sets = boot_sets if boot_sets is not None else rng.integers(0, nset, size=(B, nset))
    # per-set sums for fast resampling
    sc = np.bincount(set_idx[is_case & late], weights=z[is_case & late], minlength=nset)
    nc = np.bincount(set_idx[is_case & late], minlength=nset)
    sk = np.bincount(set_idx[~is_case], weights=z[~is_case], minlength=nset)
    nk = np.bincount(set_idx[~is_case], minlength=nset)
    bs = sc[sets].sum(1) / nc[sets].sum(1) - sk[sets].sum(1) / nk[sets].sum(1)
    return est, bs


# ---------------------------------------------------------------- 1. injection-recovery
out = []
type1 = []
for v, base_term in [("psm_fat_pct", "psm_fat_pct"), ("psm_lvol", "psm_lvol")]:
    gam = norm[(norm.measure == base_term) & (norm.term == "noise_c")].estimate.iloc[0]
    iv = intervals(v, gam=gam)
    refs = reference(iv)
    pred, sd = predict(iv, refs)
    sgn = SIGN[v]
    normal_rate = abs(norm[(norm.measure == base_term) & (norm.term == "t")].estimate.iloc[0])
    is_case = (iv.group == "case").values
    late = (iv.end_tau > -1).values
    sets = pd.Categorical(iv.set)
    set_idx, nset = sets.codes, len(sets.categories)
    B = 1000
    boot_sets = rng.integers(0, nset, size=(B, nset))
    a0, a1 = iv.tau0.values, iv.end_tau.values
    z0 = sgn * (iv.rate.values - pred) / sd
    est_obs, _ = window_stats(z0, late, is_case, set_idx, nset, B, boot_sets)   # observed (no injection)
    for L in [0.5, 1.0, 1.5, 2.0]:
        overlap = np.clip(np.minimum(a1, 0) - np.maximum(a0, -L), 0, None)
        for k in [0, 0.5, 1, 2, 3, 4, 6, 8, 12]:
            delta = k * normal_rate
            rate = iv.rate.values + np.where(is_case, sgn * delta * overlap / iv.dt.values, 0)
            z = sgn * (rate - pred) / sd
            est, bs = window_stats(z, late, is_case, set_idx, nset, B, boot_sets)
            # each bootstrap replicate is a replicate study: reject if its own one-sided 97.5% bound excludes 0;
            # SE of a replicate approximated by the bootstrap SD (the same design)
            se = bs.std()
            # replicates centred on the null (observed estimate removed), so only the injected signal remains
            power = float(np.mean(bs - est_obs - 1.96 * se > 0))
            out.append(dict(measure=v, onset_years=L, multiple_of_normal=k, delta_per_year=delta, estimate=est,
                            injected_shift=est - est_obs, se=se, power=power))
            if k == 0:
                type1.append(dict(measure=v, onset_years=L, rejection_rate=power, estimate=est))
pg = pd.DataFrame(out)
pg.to_csv(os.path.join(RES, "power_grid.csv"), index=False)
pd.DataFrame(type1).to_csv(os.path.join(RES, "power_type1.csv"), index=False)
print(pg.pivot_table(index=["measure", "multiple_of_normal"], columns="onset_years", values="power").round(2))

# ---------------------------------------------------------------- 2. specification curve
specs = []
for seg in ["psm", "aut"]:
    v = f"{seg}_fat_pct"
    for noise in [True, False]:
        for kc in [False, True]:
            for gap in [False, True]:
                iv = intervals(v, noise=noise, drop_gap90=gap, kernel_consistent=kc)
                for cov in [True, False]:
                    refs = reference(iv, covars=cov)
                    pred, sd = predict(iv, refs)
                    z = (iv.rate.values - pred) / sd
                    ok = ~np.isnan(z)
                    is_case = (iv.group == "case").values
                    sets = pd.Categorical(iv.set)
                    for cut in [0.5, 1.0, 1.5]:
                        late = (iv.end_tau > -cut).values
                        est, bs = window_stats(np.nan_to_num(z), late & ok, is_case & ok, sets.codes,
                                               len(sets.categories), 500)
                        specs.append(dict(segmentation=seg, noise_adjusted=noise, kernel_consistent_only=kc,
                                          exclude_scans_90d=gap, covariate_reference=cov, window_cut_years=cut,
                                          n_case_intervals=int((is_case & late & ok).sum()), estimate=est,
                                          lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
sc = pd.DataFrame(specs).sort_values("estimate").reset_index(drop=True)
sc.to_csv(os.path.join(RES, "spec_curve.csv"), index=False)
print(f"specification curve: {len(sc)} specifications; estimates {sc.estimate.min():.3f} to {sc.estimate.max():.3f}; "
      f"CI excludes 0 in {((sc.lo > 0) | (sc.hi < 0)).sum()}")

# ---------------------------------------------------------------- 3. raster intervals (primary measure, main settings)
gam = norm[(norm.measure == "psm_fat_pct") & (norm.term == "noise_c")].estimate.iloc[0]
iv = intervals("psm_fat_pct", gam=gam)
pred, sd = predict(iv, reference(iv))
iv["z"] = (iv.rate - pred) / sd
iv.to_csv(os.path.join(RES, "raster_intervals.csv"), index=False)
print("raster intervals:", len(iv))
