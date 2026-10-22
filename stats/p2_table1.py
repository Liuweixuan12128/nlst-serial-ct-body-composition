"""Paper 2 Table 1 and key derived numbers (analysed participants only), from the selection files and Part A outputs.
Usage: python stats/p2_table1.py   ->  results/p2/p2_table1.csv, results/p2/p2_key_numbers.csv
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = lambda f: pd.read_csv(os.path.join(ROOT, "results", "p2", f))
D = lambda f: pd.read_csv(os.path.join(ROOT, "data", "nlst_serial", f), dtype={"pid": str, "matched_case": str})
cas = pd.concat([D("p2_cases.csv").assign(resected=True), D("p2b_cases.csv").assign(resected=False)])
ctl = pd.concat([D("p2_controls.csv"), D("p2b_controls.csv")])
ser = pd.concat([D("p2_series.csv"), D("p2b_series.csv")])
maa = R("p2_maa.csv").astype({"pid": str})
prim = open(os.path.join(ROOT, "results", "p2", "partA_primary_exposure.txt")).read().split()[1]
an = maa[maa.measure == prim]
ids = set(an.pid)
cas, ctl = cas[cas.pid.isin(ids)].copy(), ctl[ctl.pid.isin(ids)].copy()
prs = pd.read_csv(os.path.join(ROOT, "data", "nlst", "nlst_prsn_public.csv"), usecols=["pid", "race", "can_scr"],
                  dtype={"pid": str})
cas = cas.drop(columns=[c for c in ["race"] if c in cas]).merge(prs, on="pid", how="left")
ctl = ctl.drop(columns=[c for c in ["race"] if c in ctl]).merge(prs[["pid", "race"]], on="pid", how="left")
nsc = an.set_index("pid").n_scans
v0 = ser.sort_values("round").groupby("pid").vendor.first()
kc = ser.groupby("pid").kernel_consistent.first()


def block(d, name):
    n = len(d)
    out = {"n": f"{n}", "Age at entry, mean (SD), y": f"{d.age.mean():.1f} ({d.age.std():.1f})",
           "Men, n (%)": f"{(d.gender == 1).sum()} ({100 * (d.gender == 1).mean():.1f})",
           "Current smoker at entry, n (%)": f"{(d.cigsmok == 1).sum()} ({100 * (d.cigsmok == 1).mean():.1f})",
           "White race, n (%)": f"{(d.race == 1).sum()} ({100 * (d.race == 1).mean():.1f})",
           "Three analysed scans, n (%)": f"{(d.pid.map(nsc) == 3).sum()} ({100 * (d.pid.map(nsc) == 3).mean():.1f})",
           "Same reconstruction kernel in all rounds, n (%)": f"{d.pid.map(kc).sum()} ({100 * d.pid.map(kc).mean():.1f})"}
    for v in ["GE", "SE", "PH", "TO"]:
        k = (d.pid.map(v0) == v)
        out[f"Scanner vendor at first scan: {v}, n (%)"] = f"{k.sum()} ({100 * k.mean():.1f})"
    return pd.Series(out, name=name)


tab = pd.concat([block(cas, "Cases"), block(ctl, "Controls")], axis=1)
st = cas.de_stag_7thed.map({110: "IA", 120: "IB", 210: "IIA", 220: "IIB", 310: "IIIA", 320: "IIIB", 400: "IV"}).fillna("Unknown")
extra = {}
for g in ["IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IV", "Unknown"]:
    extra[f"Stage (7th ed.) {g}, n (%)"] = f"{(st == g).sum()} ({100 * (st == g).mean():.1f})"
def hist(m):
    m = int(m) if pd.notna(m) else -1
    if m in (8041, 8042, 8043, 8044, 8045): return "Small cell"
    if m in (8070, 8071, 8072, 8073, 8074): return "Squamous"
    if m in (8140, 8250, 8251, 8252, 8253, 8254, 8255, 8260, 8310, 8323, 8480, 8481, 8490, 8550, 8551, 8560): return "Adenocarcinoma"
    return "Other/unspecified"
h = cas.lc_morph.map(hist)
for g in ["Adenocarcinoma", "Squamous", "Small cell", "Other/unspecified"]:
    extra[f"Histology {g}, n (%)"] = f"{(h == g).sum()} ({100 * (h == g).mean():.1f})"
extra["Resected (pathological T recorded), n (%)"] = f"{cas.resected.sum()} ({100 * cas.resected.mean():.1f})"
extra["Screen-detected, n (%)"] = f"{(cas.can_scr == 1).sum()} ({100 * (cas.can_scr == 1).mean():.1f})"
last = cas.apply(lambda r: r.candx_days - max(r[f"scr_days{k}"] for k in (0, 1, 2)
                                               if pd.notna(r[f"scr_days{k}"]) and r[f"scr_days{k}"] <= r.candx_days), axis=1)
extra["Last screening CT to diagnosis, median (IQR), d"] = f"{last.median():.0f} ({last.quantile(.25):.0f}–{last.quantile(.75):.0f})"
extra["Age at diagnosis, mean (SD), y"] = f"{(cas.age + cas.candx_days / 365.25).mean():.1f} ({(cas.age + cas.candx_days / 365.25).std():.1f})"
tab = pd.concat([tab, pd.DataFrame({"Cases": extra, "Controls": "—"})])
tab.index.name = "Characteristic"
tab.to_csv(os.path.join(ROOT, "results", "p2", "p2_table1.csv"))
print(tab.to_string())

# key numbers: normative 2-year change vs MDC, standardised slope SD, upper equivalence bounds
norm, rel, div = R("partA_normative_model.csv"), R("partA_reliability.csv"), R("partA_divergence.csv")
k = []
for m in ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]:
    a = norm[(norm.measure == m) & (norm.term == "t")].iloc[0]
    r = rel[(rel.measure == m) & (rel.subset == "all controls")].iloc[0]
    w = div[(div.measure == m) & (div.analysis == "interval z: case window vs controls")]
    k.append(dict(measure=m, annual_change=a.estimate, annual_lo=a.lo, annual_hi=a.hi, two_year_change=2 * a.estimate,
                  mdc95_2y=r.mdc95_2y, years_to_exceed_mdc=r.mdc95_2y / abs(a.estimate),
                  slope_reliability=r.slope_reliability, rel_lo=r.rel_lo, rel_hi=r.rel_hi, level_icc=r.level_icc,
                  max_upper_bound_z=w.hi.max(), min_lower_bound_z=w.lo.min(),
                  max_upper_bound_true_sd=w.hi.max() / np.sqrt(max(r.slope_reliability, 1e-6))))
kk = pd.DataFrame(k)
kk.to_csv(os.path.join(ROOT, "results", "p2", "p2_key_numbers.csv"), index=False)
print(kk.round(3).to_string())
