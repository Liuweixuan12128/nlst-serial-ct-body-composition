"""Paper 2: statistics on the shape-normalised intramuscular-fat maps (results/p2/maps/*.npz from pipeline/p2_maps.py).

Person-level maps use each participant's first QC-passing scan (cases: first pre-diagnosis scan). Cell values are the
percentage of compartment voxels in the fat window; cells with <20 voxels are missing and filled with the person's
whole-compartment value. Analysis cells = cells with >=20 voxels in >=80% of controls (the "anatomical support").
  1. controls: mean fat map, mean HU map
  2. age map (controls): per-cell OLS of fat % on age and sex; slope per 10 years; max-T permutation (1,000
     permutations of age within sex) for family-wise error
  3. case-control map: per matched set, case minus mean of its controls; one-sample t; 1,000 sign-flip permutations,
     max-T family-wise error
  4. change map: per participant (last - first scan) / years, same case-control sign-flip test
  5. craniocaudal profiles (6 bins T11 -> T6): controls' mean and case-control difference (matched-set bootstrap CI)
Outputs: results/p2/map_stats.npz, results/p2/map_level_profiles.csv, results/p2/map_summary.csv
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "p2")
rng = np.random.default_rng(20261009)
long = pd.read_csv(os.path.join(RES, "partA_long.csv"), dtype={"pid": str, "set": str})
long = long[[os.path.exists(os.path.join(RES, "maps", f"{c}.npz")) for c in long.case_id]]


def load(cid):
    z = np.load(os.path.join(RES, "maps", f"{cid}.npz"))
    return z["cnt"], z["fat"], z["hu"], z["lev_cnt"], z["lev_fat"]


def person_map(cid, min_cnt=20):
    cnt, fat, hu, lc, lf = load(cid)
    whole = 100 * fat.sum() / max(cnt.sum(), 1)
    f = np.where(cnt >= min_cnt, 100 * fat / np.maximum(cnt, 1), np.nan)
    h = np.where(cnt >= min_cnt, hu / np.maximum(cnt, 1), np.nan)
    return f, h, cnt >= min_cnt, whole, 100 * lf / np.maximum(lc, 1)


srt = long.sort_values("t")
first = srt.groupby("pid").head(1).set_index("pid")
last = srt.groupby("pid").tail(1).set_index("pid")
pids = first.index.tolist()
F, H, S, W, LEV = [], [], [], [], []
for p in pids:
    f, h, s, w, lev = person_map(first.loc[p, "case_id"])
    F.append(f); H.append(h); S.append(s); W.append(w); LEV.append(lev)
F, H, S, W, LEV = map(np.array, (F, H, S, W, LEV))
grp = first.group.values
ctrl = grp == "control"
support = S[ctrl].mean(0) >= 0.8
Ff = np.where(np.isnan(F), W[:, None, None], F)
print(f"participants {len(pids)}; support cells {support.sum()} of {support.size}")

# 1. control means
mean_fat = np.nanmean(F[ctrl], 0)
mean_hu = np.nanmean(H[ctrl], 0)


def maxT_perm_regression(Y, x, Z, nperm=1000, strata=None):
    """per-column t for x in Y ~ x + Z; family-wise p by max |t| over permutations of x (within strata)"""
    X = np.column_stack([np.ones(len(x)), x, Z])
    def tstats(xv):
        Xv = X.copy(); Xv[:, 1] = xv
        beta, *_ = np.linalg.lstsq(Xv, Y, rcond=None)
        r = Y - Xv @ beta
        s2 = (r ** 2).sum(0) / (len(Y) - Xv.shape[1])
        cinv = np.linalg.inv(Xv.T @ Xv)[1, 1]
        return beta[1], beta[1] / np.sqrt(s2 * cinv)
    b, t = tstats(x)
    mx = np.empty(nperm)
    for i in range(nperm):
        xp = x.copy()
        if strata is None:
            xp = rng.permutation(x)
        else:
            for s in np.unique(strata):
                m = strata == s
                xp[m] = rng.permutation(x[m])
        mx[i] = np.abs(tstats(xp)[1]).max()
    p_fwe = (1 + (mx[None, :] >= np.abs(t)[:, None]).sum(1)) / (nperm + 1)
    return b, t, p_fwe


def signflip_maxT(D, nperm=1000):
    n = len(D)
    t = D.mean(0) / (D.std(0, ddof=1) / np.sqrt(n))
    mx = np.empty(nperm)
    for i in range(nperm):
        s = rng.choice([-1, 1], size=n)[:, None]
        Dp = D * s
        mx[i] = np.abs(Dp.mean(0) / (Dp.std(0, ddof=1) / np.sqrt(n))).max()
    p = (1 + (mx[None, :] >= np.abs(t)[:, None]).sum(1)) / (nperm + 1)
    return D.mean(0), t, p


cells = support.ravel()
# 2. age map (controls), age per 10 years, adjusted for sex, permuted within sex
Yc = Ff[ctrl].reshape(ctrl.sum(), -1)[:, cells]
age = first.age0.values[ctrl].astype(float) / 10
sex = (first.sex.values[ctrl] == "female").astype(float)
b_age, t_age, p_age = maxT_perm_regression(Yc, age - age.mean(), (sex - sex.mean())[:, None], strata=sex)
# 3. case-control map, first scan
sets = first.set.values
case_idx = {s: i for i, (s, g) in enumerate(zip(sets, grp)) if g == "case"}
Dl = []
Fv = Ff.reshape(len(pids), -1)[:, cells]
for s, i in case_idx.items():
    cm = np.where((sets == s) & ctrl)[0]
    if len(cm):
        Dl.append(Fv[i] - Fv[cm].mean(0))
d_cc, t_cc, p_cc = signflip_maxT(np.array(Dl))
# 4. change map (participants with >=2 scans and >=0.8 y apart)
Cm, keep = [], []
for k, p in enumerate(pids):
    c0, c1 = first.loc[p, "case_id"], last.loc[p, "case_id"]
    dt = last.loc[p, "t"] - first.loc[p, "t"]
    if c0 == c1 or dt < 0.8:
        Cm.append(np.full(cells.sum(), np.nan)); continue
    f1, _, _, w1, _ = person_map(c1)
    f1 = np.where(np.isnan(f1), w1, f1).ravel()[cells]
    Cm.append((f1 - Fv[k]) / dt)
Cm = np.array(Cm)
Dc = []
for s, i in case_idx.items():
    cm = np.where((sets == s) & ctrl & ~np.isnan(Cm[:, 0]))[0]
    if len(cm) and not np.isnan(Cm[i, 0]):
        Dc.append(Cm[i] - Cm[cm].mean(0))
d_ch, t_ch, p_ch = signflip_maxT(np.array(Dc))
ctrl_change = np.nanmean(Cm[ctrl], 0)


def unflat(v):
    m = np.full(support.size, np.nan)
    m[cells] = v
    return m.reshape(support.shape)


np.savez_compressed(os.path.join(RES, "map_stats.npz"), support=support, mean_fat=mean_fat, mean_hu=mean_hu,
                    coverage=S[ctrl].mean(0), age_b=unflat(b_age), age_t=unflat(t_age), age_p=unflat(p_age),
                    cc_d=unflat(d_cc), cc_t=unflat(t_cc), cc_p=unflat(p_cc), ch_d=unflat(d_ch), ch_t=unflat(t_ch),
                    ch_p=unflat(p_ch), ctrl_change=unflat(ctrl_change))
# 5. craniocaudal profiles
rows = []
for b in range(LEV.shape[1]):
    cm = LEV[ctrl, b]
    diffs = []
    for s, i in case_idx.items():
        cc = np.where((sets == s) & ctrl)[0]
        if len(cc):
            diffs.append(LEV[i, b] - LEV[cc, b].mean())
    diffs = np.array(diffs)
    bs = [rng.choice(diffs, len(diffs)).mean() for _ in range(1000)]
    rows.append(dict(level_bin=b, ctrl_mean=cm.mean(), ctrl_lo=np.percentile(cm, 25), ctrl_hi=np.percentile(cm, 75),
                     case_mean=LEV[grp == "case", b].mean(), diff=diffs.mean(), diff_lo=np.percentile(bs, 2.5),
                     diff_hi=np.percentile(bs, 97.5)))
pd.DataFrame(rows).to_csv(os.path.join(RES, "map_level_profiles.csv"), index=False)
summ = dict(n_participants=len(pids), n_sets_level=len(Dl), n_sets_change=len(Dc), support_cells=int(support.sum()),
            age_cells_fwe05=int((p_age < 0.05).sum()), age_slope_median=float(np.median(b_age)),
            cc_cells_fwe05=int((p_cc < 0.05).sum()), cc_min_p=float(p_cc.min()),
            change_cells_fwe05=int((p_ch < 0.05).sum()), change_min_p=float(p_ch.min()))
pd.Series(summ).to_csv(os.path.join(RES, "map_summary.csv"))
print(summ)
