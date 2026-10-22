"""Paper 2: participant-flow counts for Figure 1 / the STROBE diagram, recomputed from the public IDC tables with the
same definitions as p2_select_serial.py (nothing here changes the selection).

Usage: python p2_flow_counts.py --out results/p2/p2_flow_counts.csv
"""
import argparse, os, sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nlst_select_download import table  # noqa: E402
from p2_select_serial import series_index  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from idc_index import IDCClient
    c = IDCClient()
    nl_all = c.index[(c.index.collection_id == "nlst") & (c.index.Modality == "CT")]
    ser = series_index(c)
    prsn, canc = table("nlst_prsn"), table("nlst_canc")
    for t in (prsn, canc):
        t["pid"] = t["pid"].astype(str)
        for col in ["candx_days", "scr_days0", "scr_days1", "scr_days2"]:
            if col in t:
                t[col] = pd.to_numeric(t[col], errors="coerce")
    img = ser.groupby("PatientID")["round"].apply(frozenset)
    rows = [("NLST participants in the public person table", prsn.pid.nunique()),
            ("participants with CT in IDC", nl_all.PatientID.astype(str).nunique()),
            ("participants with an original axial screening series (>=40 images)", len(img))]
    first = canc[(canc.first_lc == 1)].merge(prsn[["pid", "scr_days0", "scr_days1", "scr_days2"]], on="pid")
    first = first[first.pid.isin(img.index)]
    rows.append(("first primary lung cancer, imaged", len(first)))
    n_prior = first.apply(lambda r: sum(1 for k in (0, 1, 2) if pd.notna(r[f"scr_days{k}"])
                                        and r[f"scr_days{k}"] <= r.candx_days and k in img[r.pid]), axis=1)
    rows += [("  with 0-1 imaged round on/before diagnosis (excluded)", int((n_prior < 2).sum())),
             ("  with >=2 imaged rounds on/before diagnosis", int((n_prior >= 2).sum()))]
    pool = prsn[~prsn.pid.isin(set(canc.pid)) & (prsn.cancyr.astype(str) == "")]
    rows += [("never diagnosed with lung cancer", len(pool)),
             ("  of whom all three rounds imaged (control pool)",
              int(pool.pid.map(lambda p: img.get(p, frozenset()) >= {0, 1, 2}).sum()))]
    for s in ("p2", "p2b"):
        for kind in ("cases", "controls", "series"):
            f = os.path.join(ROOT, "data", "nlst_serial", f"{s}_{kind}.csv")
            rows.append((f"selected {s} {kind}", len(pd.read_csv(f))))
    out = pd.DataFrame(rows, columns=["step", "n"])
    out.to_csv(a.out, index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
