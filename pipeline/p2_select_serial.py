"""Paper 2 (serial muscle ageing): cases with >=2 pre-diagnosis screening CTs and matched cancer-free controls.

Public IDC tables only (nlst_prsn, nlst_canc, nlst_screen); no CDAS variables are used here.

Cases: CT-arm first primary lung cancer with a pathological T category (resection proxy, confirmed later against the
CDAS treatment table), images in IDC for >=2 screening rounds performed on or before the diagnosis date.
Controls: never diagnosed with lung cancer in the public cancer table, lung-cancer-free at the case's diagnosis day
(canc_free_days >= candx_days), images in IDC for all three rounds; matched 1:3 without replacement on sex, smoking
status at entry (cigsmok), age at entry (exact, widened to +-1 then +-2 years) and T0 scanner vendor.
Series per round: original primary axial, >=40 images; a soft kernel shared by every round is preferred so that
within-person change is not driven by a kernel switch; ties go to the thinnest slice <= 3 mm.

Usage: python p2_select_serial.py --out data/nlst_serial [--ratio 3] [--seed 20261009]
       python p2_select_serial.py --out data/nlst_serial --scope unresected --prefix p2b \
              --exclude-controls data/nlst_serial/p2_controls.csv      (adds the non-resected cancers)
Writes: p2_cases.csv, p2_controls.csv (one row per participant) and p2_series.csv (one row per series to process).
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nlst_select_download import SOFT, YEAR_OF_ROUND, parse_desc, table  # noqa: E402


def series_index(client):
    nl = client.index[(client.index.collection_id == "nlst") & (client.index.Modality == "CT")].copy()
    nl["PatientID"] = nl["PatientID"].astype(str)
    parsed = nl.SeriesDescription.map(parse_desc)
    nl = nl[parsed.notna()].copy()
    nl = pd.concat([nl, pd.DataFrame(list(parsed.dropna()), index=nl.index)], axis=1)
    nl = nl[(nl.img_type == "OPA") & (nl.instanceCount >= 40)].copy()
    nl["round"] = pd.to_datetime(nl.StudyDate).dt.year.map({v: k for k, v in YEAR_OF_ROUND.items()})
    nl = nl[nl["round"].notna()].copy()
    nl["round"] = nl["round"].astype(int)
    nl["soft"] = nl.kernel.isin(SOFT)
    nl["thin"] = nl.thickness <= 3
    return nl


def pick_series(s, rounds):
    """Choose one series per round, preferring a kernel available in every round."""
    s = s[s["round"].isin(rounds)]
    if set(s["round"]) != set(rounds):
        return None
    common = set.intersection(*[set(s[s["round"] == k].kernel) for k in rounds])
    def score(d):
        return (d.kernel.isin(common).astype(int) * 100 + d.soft.astype(int) * 10 + d.thin.astype(int) * 5
                - d.thickness.fillna(9) / 10)
    s = s.assign(score=score(s))
    if common:  # restrict to the best common kernel so all rounds share it
        best = (s[s.kernel.isin(common)].groupby("kernel").score.max().sort_values(ascending=False).index[0])
        s = s[s.kernel == best]
    out = s.sort_values("score", ascending=False).groupby("round").head(1).sort_values("round")
    return out.assign(kernel_consistent=int(out.kernel.nunique() == 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ratio", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20261009)
    ap.add_argument("--scope", choices=["resected", "unresected"], default="resected",
                    help="resected: pathological T recorded (original set); unresected: first primaries without one")
    ap.add_argument("--exclude-controls", default="", help="controls already used by an earlier selection")
    ap.add_argument("--prefix", default="p2")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from idc_index import IDCClient
    ser = series_index(IDCClient())
    prsn, canc = table("nlst_prsn"), table("nlst_canc")
    for t in (prsn, canc):
        t["pid"] = t["pid"].astype(str)
        for col in ["age", "candx_days", "scr_days0", "scr_days1", "scr_days2", "canc_free_days"]:
            if col in t:
                t[col] = pd.to_numeric(t[col], errors="coerce")
    rounds_img = ser.groupby("PatientID")["round"].apply(lambda r: frozenset(r))

    has_pt = canc.path_t.astype(str).isin(["100", "200", "300", "400"])
    first = canc[(canc.first_lc == 1) & (has_pt if a.scope == "resected" else ~has_pt)]
    cases = first.merge(prsn[["pid", "age", "gender", "cigsmok", "race", "scr_days0", "scr_days1", "scr_days2"]],
                        on="pid")
    cases = cases[cases.pid.isin(rounds_img.index)].copy()

    def prior_rounds(r):
        have = rounds_img.get(r.pid, frozenset())
        return tuple(k for k in (0, 1, 2)
                     if pd.notna(r[f"scr_days{k}"]) and r[f"scr_days{k}"] <= r.candx_days and k in have)
    cases["rounds"] = cases.apply(prior_rounds, axis=1)
    n_all = len(cases)
    cases = cases[cases.rounds.map(len) >= 2].copy()
    cases["age_at_dx"] = (cases.age + cases.candx_days / 365.25).round(1)

    sel_rows, case_rows = [], []
    for _, r in cases.iterrows():
        pk = pick_series(ser[ser.PatientID == r.pid], list(r.rounds))
        if pk is None:
            continue
        sel_rows.append(pk.assign(pid=r.pid, group="case"))
        case_rows.append(r)
    cases = pd.DataFrame(case_rows)
    vendor_t0 = {p: d.sort_values("round").vendor.iloc[0] for p, d in pd.concat(sel_rows).groupby("pid")}
    cases["vendor_first"] = cases.pid.map(vendor_t0)

    # control pool: never lung cancer, all three rounds imaged
    pool = prsn[~prsn.pid.isin(set(canc.pid)) & (prsn.cancyr.astype(str) == "")].copy()
    if a.exclude_controls:
        for f in a.exclude_controls.split(","):
            pool = pool[~pool.pid.isin(set(pd.read_csv(f, dtype={"pid": str}).pid))]
    pool = pool[pool.pid.map(lambda p: rounds_img.get(p, frozenset()) >= {0, 1, 2})].copy()
    pool_ser = ser[ser.PatientID.isin(set(pool.pid))]
    v0 = pool_ser[pool_ser["round"] == 0].groupby("PatientID").vendor.agg(lambda v: v.mode().iloc[0])
    pool["vendor_first"] = pool.pid.map(v0)
    rng = np.random.default_rng(a.seed)
    used, ctrl_rows, ctrl_ser = set(), [], []
    # hardest-to-match cases first (oldest, then women) so the pool is not exhausted by easy strata
    order = cases.sort_values(["age", "gender"], ascending=[False, False])
    for _, r in order.iterrows():
        chosen = []
        for tol in (0, 1, 2):
            m = pool[(pool.gender == r.gender) & (pool.cigsmok == r.cigsmok) & (pool.vendor_first == r.vendor_first)
                     & ((pool.age - r.age).abs() <= tol) & (pool.canc_free_days >= r.candx_days)
                     & ~pool.pid.isin(used | set(chosen))]
            need = a.ratio - len(chosen)
            if len(m):
                take = m.pid.values[rng.permutation(len(m))[:need]]
                for p in take:
                    pk = pick_series(pool_ser[pool_ser.PatientID == p], [0, 1, 2])
                    if pk is not None and len(chosen) < a.ratio:
                        chosen.append(p)
                        ctrl_ser.append(pk.assign(pid=p, group="control"))
            if len(chosen) >= a.ratio:
                break
        used |= set(chosen)
        for p in chosen:
            ctrl_rows.append(dict(pid=p, matched_case=r.pid, index_days=r.candx_days))
    ctrl = pd.DataFrame(ctrl_rows).merge(
        pool[["pid", "age", "gender", "cigsmok", "race", "canc_free_days", "vendor_first",
              "scr_days0", "scr_days1", "scr_days2"]], on="pid")

    keep = ["pid", "group", "round", "SeriesInstanceUID", "SeriesDescription", "kernel", "soft", "thickness",
            "vendor", "instanceCount", "series_size_MB", "kernel_consistent"]
    allser = pd.concat(sel_rows + ctrl_ser)[keep].rename(columns={"instanceCount": "n_images"})
    allser["case_id"] = allser.pid + "_T" + allser["round"].astype(str)

    cases["n_rounds"] = cases.rounds.map(len)
    cases["rounds"] = cases.rounds.map(lambda t: "".join(map(str, t)))
    cases.to_csv(os.path.join(a.out, f"{a.prefix}_cases.csv"), index=False)
    ctrl.to_csv(os.path.join(a.out, f"{a.prefix}_controls.csv"), index=False)
    allser.to_csv(os.path.join(a.out, f"{a.prefix}_series.csv"), index=False)
    n_match = ctrl.groupby("matched_case").size().reindex(cases.pid, fill_value=0)
    print(f"[{a.scope}] imaged first primaries: {n_all}; with >=2 pre-dx rounds and series: {len(cases)} "
          f"(3 rounds: {(cases.n_rounds == 3).sum()}; aged >=65 at dx: {(cases.age_at_dx >= 65).sum()})\n"
          f"controls: {len(ctrl)} (cases with {a.ratio} controls: {(n_match == a.ratio).sum()}, "
          f"<{a.ratio}: {(n_match < a.ratio).sum()}, 0: {(n_match == 0).sum()})\n"
          f"series: {len(allser)} ({allser.series_size_MB.sum() / 1e3:.0f} GB DICOM); kernel-consistent "
          f"participants: cases {allser[allser.group == 'case'].groupby('pid').kernel_consistent.first().mean():.0%}, "
          f"controls {allser[allser.group == 'control'].groupby('pid').kernel_consistent.first().mean():.0%}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
