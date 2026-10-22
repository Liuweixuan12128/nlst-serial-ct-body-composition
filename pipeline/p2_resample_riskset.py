"""Paper 2: re-draw the matched controls so that eligibility is judged at the case's diagnosis date (risk set).

The original selection (p2_select_serial.py) required controls to have images in all three screening rounds, while a
case needed only >=2 rounds before diagnosis. For a case diagnosed before T2, that made attendance after the case's
diagnosis a condition for being a control. Here a control is eligible if, like the original pool, never diagnosed with
lung cancer and still free of it on the case's diagnosis day, and imaged at the same rounds the case had before
diagnosis (other rounds may be missing). The matching factors are unchanged.

Every case's new risk set R contains its original eligible pool A. The original controls are a random sample of A, so a
random sample of 3 from R is obtained by keeping a hypergeometric number of them, n_A ~ Hypergeometric(|R|, |A|, 3), and
drawing the other 3 - n_A from R minus A. This reuses the processed scans; only the new controls need downloading.
For cases with three pre-diagnostic rounds R = A, so their controls are unchanged.

Usage: python pipeline/p2_resample_riskset.py [--seed 20261015]
Writes data/nlst_serial/{p2,p2b}_controls_riskset.csv (same columns as the original controls files, plus kept=0/1)
and {p2,p2b}_series_riskset_new.csv (series to fetch for the new controls).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nlst_select_download import table  # noqa: E402
from p2_select_serial import pick_series, series_index  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SER = os.path.join(ROOT, "data", "nlst_serial")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20261015)
    ap.add_argument("--ratio", type=int, default=3)
    a = ap.parse_args()
    from idc_index import IDCClient
    ser = series_index(IDCClient())
    prsn, canc = table("nlst_prsn"), table("nlst_canc")
    for t in (prsn, canc):
        t["pid"] = t["pid"].astype(str)
        for col in ["age", "candx_days", "scr_days0", "scr_days1", "scr_days2", "canc_free_days"]:
            if col in t:
                t[col] = pd.to_numeric(t[col], errors="coerce")
    rounds_img = ser.groupby("PatientID")["round"].apply(lambda r: frozenset(r))

    # never-lung-cancer pool, as in the original selection, without the all-three-rounds condition
    pool = prsn[~prsn.pid.isin(set(canc.pid)) & (prsn.cancyr.astype(str) == "")].copy()
    pool["img"] = pool.pid.map(lambda p: rounds_img.get(p, frozenset()))
    pool = pool[pool.img.map(len) >= 2].copy()
    pool_ser = ser[ser.PatientID.isin(set(pool.pid))]
    vend = pool_ser.groupby(["PatientID", "round"]).vendor.agg(lambda v: v.mode().iloc[0]).unstack()
    pool["all3"] = pool.img.map(lambda s: s >= {0, 1, 2})
    pool = pool.set_index("pid", drop=False)
    by_pid = dict(tuple(pool_ser.groupby("PatientID")))
    pick_cache = {}

    def can_pick(p, rounds):
        k = (p, rounds)
        if k not in pick_cache:
            pick_cache[k] = pick_series(by_pid[p], list(rounds))
        return pick_cache[k]

    orig = {s: pd.read_csv(os.path.join(SER, f"{s}_controls.csv"), dtype={"pid": str, "matched_case": str})
            for s in ("p2", "p2b")}
    used = set(pd.concat(orig.values()).pid)          # original controls stay unavailable as new draws
    rng = np.random.default_rng(a.seed)
    log = []
    for s in ("p2", "p2b"):
        cases = pd.read_csv(os.path.join(SER, f"{s}_cases.csv"), dtype={"pid": str, "rounds": str})
        oc = orig[s]
        rows, new_ser = [], []
        for _, r in cases.sort_values(["age", "gender"], ascending=[False, False]).iterrows():
            rounds = tuple(int(c) for c in r.rounds)
            mine = oc[oc.matched_case == r.pid]
            if len(rounds) == 3 or not len(mine):
                rows.append(mine.assign(kept=1))
                log.append(dict(series=s, case=r.pid, rounds=r.rounds, R=np.nan, A=np.nan, n_new=0))
                continue
            tol = int((mine.age - r.age).abs().max())
            v_case = r.vendor_first
            first = rounds[0]
            m = pool[(pool.gender == r.gender) & (pool.cigsmok == r.cigsmok) & ((pool.age - r.age).abs() <= tol)
                     & (pool.canc_free_days >= r.candx_days) & pool.img.map(lambda x: x >= set(rounds))]
            # vendor at the case's first round, as for the case itself
            m = m[m.pid.map(lambda p: vend.at[p, first] if first in vend.columns and p in vend.index else None)
                  == v_case]
            m = m[~m.pid.isin(used - set(mine.pid))]
            m = m[[can_pick(p, rounds) is not None for p in m.pid]]
            A = m[m.all3 & m.pid.map(lambda p: can_pick(p, (0, 1, 2)) is not None)]
            B = m[~m.pid.isin(set(A.pid))]
            nA = rng.hypergeometric(len(A), len(B), a.ratio) if len(B) else a.ratio
            nA = min(max(nA, a.ratio - len(B)), len(mine))
            keep_idx = rng.permutation(len(mine))[:nA]
            kept = mine.iloc[sorted(keep_idx)].assign(kept=1)
            n_new = a.ratio - nA
            take = B.pid.values[rng.permutation(len(B))[:n_new]] if n_new else []
            new = []
            for p in take:
                pk = can_pick(p, rounds)
                new_ser.append(pk.assign(pid=p, group="control"))
                q = pool.loc[p]
                new.append(dict(pid=p, matched_case=r.pid, index_days=r.candx_days, age=q.age, gender=q.gender,
                                cigsmok=q.cigsmok, race=q.race, canc_free_days=q.canc_free_days,
                                vendor_first=v_case, scr_days0=q.scr_days0, scr_days1=q.scr_days1,
                                scr_days2=q.scr_days2, kept=0))
                used.add(p)
            rows.append(pd.concat([kept, pd.DataFrame(new)], ignore_index=True))
            log.append(dict(series=s, case=r.pid, rounds=r.rounds, R=len(m), A=len(A), n_new=n_new))
        ctl = pd.concat(rows, ignore_index=True)
        ctl.to_csv(os.path.join(SER, f"{s}_controls_riskset.csv"), index=False)
        keep = ["pid", "group", "round", "SeriesInstanceUID", "SeriesDescription", "kernel", "soft", "thickness",
                "vendor", "instanceCount", "series_size_MB", "kernel_consistent"]
        ns = (pd.concat(new_ser)[keep].rename(columns={"instanceCount": "n_images"}) if new_ser
              else pd.DataFrame(columns=keep))
        ns["case_id"] = ns.pid + "_T" + ns["round"].astype(str)
        ns.to_csv(os.path.join(SER, f"{s}_series_riskset_new.csv"), index=False)
        print(f"[{s}] controls {len(ctl)} (kept {int(ctl.kept.sum())}, new {int((ctl.kept == 0).sum())}); "
              f"new series {len(ns)} ({ns.series_size_MB.sum() / 1e3:.0f} GB)", file=sys.stderr)
    lg = pd.DataFrame(log)
    lg.to_csv(os.path.join(SER, "riskset_draw_log.csv"), index=False)
    two = lg[lg.rounds.str.len() == 2]
    print(f"two-round cases: {len(two)}; median |R| {two.R.median():.0f}, median share outside the original pool "
          f"{((two.R - two.A) / two.R).median():.1%}; new controls {int(two.n_new.sum())}", file=sys.stderr)


if __name__ == "__main__":
    main()
