"""Paper 2: descriptive quantities for the additional figure panels (no new hypothesis tests).

Reads results/p2 outputs of stats/p2_partA.R and stats/p2_phenotype.R. All quantities are descriptive summaries of the
same analysis population; confidence intervals are matched-set bootstrap (1,000 replicates, seed 20261012) where given.
  1. extra_score_by_time.csv    reserve score at every QC-passing scan, cases and controls, by time to the case's
                                diagnosis (control reference models from controls' first scans, as in the score plan)
  2. extra_component_corr.csv   Spearman correlations of the five z-scores in controls (first scan)
  3. extra_examples.csv         rule-chosen illustration: the control men aged 60-64 at entry whose first-scan score
                                is nearest the 10th and the 90th centile of all controls' scores
  4. extra_shift_function.csv   case minus control difference at each decile of each z-score (shift function)
  5. extra_cross_vs_long.csv    cross-sectional age gradient (controls' first scans) versus longitudinal annual change
  6. extra_variance_explained.csv  share of between-person variance (controls' first scans) explained by age, sex,
                                smoking status and scanner (vendor + noise), drop-one R2
  7. extra_projected_reliability.csv  slope reliability for annual scans over 1-10 years, from the variance components
  8. extra_kernel_change.csv    individual radiodensity slopes of controls, kernel kept versus switched
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "p2")
rng = np.random.default_rng(20261012)


class OLS:
    """least squares on numeric and categorical columns (treatment coding, first level as reference)"""
    def __init__(self, df, y, num=(), cat=()):
        self.num, self.cat = list(num), list(cat)
        self.levels = {c: sorted(df[c].dropna().unique()) for c in self.cat}
        X = self.design(df)
        yy = df[y].values.astype(float)
        self.beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
        res = yy - X @ self.beta
        self.dof = len(yy) - X.shape[1]
        self.scale = (res ** 2).sum() / self.dof
        self.rsquared = 1 - (res ** 2).sum() / ((yy - yy.mean()) ** 2).sum()
        self.cov = self.scale * np.linalg.pinv(X.T @ X)
        self.names = ["const"] + self.num + [f"{c}[{l}]" for c in self.cat for l in self.levels[c][1:]]

    def design(self, df):
        cols = [np.ones(len(df))] + [df[c].values.astype(float) for c in self.num]
        for c in self.cat:
            for l in self.levels[c][1:]:
                cols.append((df[c].values == l).astype(float))
        return np.column_stack(cols)

    def predict(self, df):
        return self.design(df) @ self.beta

    def coef(self, name):
        i = self.names.index(name)
        se = np.sqrt(self.cov[i, i])
        return self.beta[i], self.beta[i] - 1.96 * se, self.beta[i] + 1.96 * se

long = pd.read_csv(os.path.join(RES, "partA_long.csv"), dtype={"pid": str, "set": str})
norm = pd.read_csv(os.path.join(RES, "partA_normative_model.csv"))
rel = pd.read_csv(os.path.join(RES, "partA_reliability.csv"))
sc = pd.read_csv(os.path.join(RES, "pheno_scores.csv"), dtype={"pid": str, "set": str})
COMP = {"T10_trab_hu": 1, "psm_lvol": 1, "lung_perc15_hu": 1, "psm_fat_pct": -1, "psm_hu": 1}
gam = {v: norm[(norm.measure == v) & (norm.term == "noise_c")].estimate.iloc[0] for v in COMP}
long["age0c"] = long.age0 - 62
for v in COMP:
    long[f"y_{v}"] = long[v] - gam[v] * long.noise_c
long = long.sort_values(["pid", "t"])
first = long.groupby("pid").head(1).copy()

# ---------------------------------------------------------------- 1. score at every scan
fits = {}
cf = first[first.group == "control"]
for v in COMP:
    d = cf[cf[f"y_{v}"].notna()]
    fits[v] = OLS(d, f"y_{v}", num=["age0c"], cat=["sex", "smoke", "vendor0"])


def zscores(df):
    out = pd.DataFrame(index=df.index)
    for v, sgn in COMP.items():
        dd = df.copy()
        known = set(fits[v].levels["vendor0"])
        dd.loc[~dd.vendor0.isin(known), "vendor0"] = cf.vendor0.mode().iloc[0]
        out[f"z_{v}"] = sgn * (dd[f"y_{v}"] - fits[v].predict(dd)) / np.sqrt(fits[v].scale)
    return out


zf = zscores(first)
sd3 = zf.loc[first.group == "control", ["z_T10_trab_hu", "z_psm_lvol", "z_lung_perc15_hu"]].mean(axis=1).std()
za = zscores(long)
long["S_scan"] = za[["z_T10_trab_hu", "z_psm_lvol", "z_lung_perc15_hu"]].mean(axis=1) / sd3
analysed_sets = set(sc[(sc.group == "case") & sc.S.notna()].set)
L = long[long.set.isin(analysed_sets) & long.S_scan.notna()].copy()
bins = [-7.5, -5, -4, -3, -2, -1, 0.01]
L["bin"] = pd.cut(L.tau, bins, right=False, labels=["<−5", "−5 to −4", "−4 to −3", "−3 to −2", "−2 to −1", "−1 to 0"])
rows = []
for b, g in L.groupby("bin", observed=True):
    r = dict(bin=b, mean_tau=g[g.group == "case"].tau.mean() if (g.group == "case").any() else g.tau.mean())
    for grp in ("case", "control"):
        x = g[g.group == grp].S_scan
        r.update({f"{grp}_n": len(x), f"{grp}_mean": x.mean(), f"{grp}_lo": x.mean() - 1.96 * x.std() / np.sqrt(len(x)),
                  f"{grp}_hi": x.mean() + 1.96 * x.std() / np.sqrt(len(x))})
    # case minus mean of the same sets' controls in the bin, matched-set bootstrap
    dif = g.groupby(["set", "group"]).S_scan.mean().unstack()
    dif = (dif["case"] - dif["control"]).dropna() if {"case", "control"} <= set(dif.columns) else pd.Series(dtype=float)
    bs = [rng.choice(dif.values, len(dif)).mean() for _ in range(1000)] if len(dif) > 10 else [np.nan]
    r.update(diff=dif.mean(), diff_lo=np.percentile(bs, 2.5), diff_hi=np.percentile(bs, 97.5), n_sets=len(dif))
    rows.append(r)
pd.DataFrame(rows).to_csv(os.path.join(RES, "extra_score_by_time.csv"), index=False)
print(pd.DataFrame(rows)[["bin", "case_n", "control_n", "case_mean", "control_mean", "diff", "diff_lo", "diff_hi"]].round(3))

# ---------------------------------------------------------------- 2. component correlations (controls)
zc = sc[sc.group == "control"][["z_lung_perc15_hu", "z_T10_trab_hu", "z_psm_lvol", "z_psm_hu", "z_psm_fat_pct"]]
zc.corr(method="spearman").to_csv(os.path.join(RES, "extra_component_corr.csv"))

# ---------------------------------------------------------------- 3. rule-chosen examples
c = sc[(sc.group == "control") & sc.S.notna()].merge(first[["pid", "case_id", "age0", "sex"]], on="pid")
pool = c[(c.sex == "male") & c.age0.between(60, 64)]
ex = []
for q in (0.1, 0.9):
    target = c.S.quantile(q)
    r = pool.iloc[(pool.S - target).abs().argmin()]
    ex.append(dict(centile=int(q * 100), target=target, pid=r.pid, case_id=r.case_id, age0=r.age0, S=r.S,
                   z_bone=r.z_T10_trab_hu, z_muscle=r.z_psm_lvol, z_lung=r.z_lung_perc15_hu))
pd.DataFrame(ex).to_csv(os.path.join(RES, "extra_examples.csv"), index=False)
print(pd.DataFrame(ex).round(2))

# ---------------------------------------------------------------- 4. shift function
s = sc[sc.S.notna() & sc.set.isin(analysed_sets)]
sets = np.array(sorted(analysed_sets))
by = {k: g for k, g in s.groupby("set")}
qs = np.arange(1, 10) / 10
rows = []
for v in ["S", "z_lung_perc15_hu", "z_T10_trab_hu", "z_psm_lvol"]:
    obs = np.quantile(s[s.group == "case"][v], qs) - np.quantile(s[s.group == "control"][v], qs)
    bs = []
    for _ in range(1000):
        pick = pd.concat([by[k] for k in rng.choice(sets, len(sets))])
        bs.append(np.quantile(pick[pick.group == "case"][v], qs) - np.quantile(pick[pick.group == "control"][v], qs))
    bs = np.array(bs)
    for i, q in enumerate(qs):
        rows.append(dict(measure=v, decile=q, diff=obs[i], lo=np.percentile(bs[:, i], 2.5), hi=np.percentile(bs[:, i], 97.5)))
pd.DataFrame(rows).to_csv(os.path.join(RES, "extra_shift_function.csv"), index=False)
print(pd.DataFrame(rows).pivot(index="decile", columns="measure", values="diff").round(2))

# ---------------------------------------------------------------- 5. cross-sectional vs longitudinal
rows = []
for v in ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]:
    d = cf[cf[v].notna()]
    m = OLS(d, v, num=["age0", "noise_c"], cat=["sex", "smoke", "vendor0"])
    est, lo, hi = m.coef("age0")
    t = norm[(norm.measure == v) & (norm.term == "t")].iloc[0]
    sd = d[v].std()
    rows.append(dict(measure=v, cross=est, cross_lo=lo, cross_hi=hi, long=t.estimate, long_lo=t.lo,
                     long_hi=t.hi, control_sd=sd))
pd.DataFrame(rows).to_csv(os.path.join(RES, "extra_cross_vs_long.csv"), index=False)
print(pd.DataFrame(rows).round(3))

# ---------------------------------------------------------------- 6. variance explained (drop-one R2)
rows = []
for v in ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]:
    d = cf[cf[v].notna()]
    terms = {"Age": (["age0"], []), "Sex": ([], ["sex"]), "Smoking": ([], ["smoke"]),
             "Scanner": (["noise_c"], ["vendor0"])}
    def r2(keep):
        num = sum((terms[k][0] for k in keep), []); cat = sum((terms[k][1] for k in keep), [])
        return OLS(d, v, num=num, cat=cat).rsquared
    full = r2(list(terms))
    r = dict(measure=v, full=full)
    for k in terms:
        r[k] = full - r2([kk for kk in terms if kk != k])
    r["Shared"] = full - sum(r[k] for k in terms)
    r["Unexplained"] = 1 - full
    rows.append(r)
pd.DataFrame(rows).to_csv(os.path.join(RES, "extra_variance_explained.csv"), index=False)
print(pd.DataFrame(rows).round(3))

# ---------------------------------------------------------------- 7. projected slope reliability
rows = []
for _, r in rel[rel.subset == "all controls"].iterrows():
    for T in range(1, 11):
        n = T + 1
        sxx = n * (n ** 2 - 1) / 12
        rr = r.slope_true_sd ** 2 / (r.slope_true_sd ** 2 + r.resid_sd ** 2 / sxx)
        rows.append(dict(measure=r.measure, years=T, n_scans=n, reliability=rr))
pd.DataFrame(rows).to_csv(os.path.join(RES, "extra_projected_reliability.csv"), index=False)

# ---------------------------------------------------------------- 8. kernel kept vs switched: individual slopes
maa = pd.read_csv(os.path.join(RES, "p2_maa.csv"), dtype={"pid": str})
kc = long.groupby("pid").kernel_consistent.first()
k = maa[maa.measure == "psm_hu"].merge(kc.rename("kernel_consistent"), left_on="pid", right_index=True)
k = k[k.pid.isin(set(long[long.group == "control"].pid))]
k[["pid", "kernel_consistent", "slope", "n_scans", "span"]].to_csv(os.path.join(RES, "extra_kernel_change.csv"),
                                                                   index=False)
print(k.groupby("kernel_consistent").slope.agg(["count", "mean", "std"]).round(3))
