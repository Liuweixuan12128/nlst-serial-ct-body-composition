"""NLST (CT arm) surgical candidates: pick the last screening LDCT before lung-cancer diagnosis and download it.

Candidates (public IDC clinical tables, no application needed): first primary lung cancer, participant has CT in
IDC, and a pathological T category (T1-T4) is recorded -- pathological staging requires a resection specimen, so
this is a proxy for surgical resection that is CONFIRMED later against the CDAS treatment dataset.
Screen rounds map to IDC StudyDate years: T0 = 1999, T1 = 2000, T2 = 2001.
Series: original primary axial ("OPA"), soft kernel preferred, thinnest <= 3 mm.
Series description format: "<screen>,<OPA|OPL..>,<vendor>,<model>,<kernel>,<fov>,<thickness>,<kVp>,<mA>,...".

Usage: python nlst_select_download.py --out data/nlst [--download] [--rounds last|all]
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd

CLIN = os.path.expanduser("~/.local/share/IDC/v24/clinical_data")
SOFT = {"STANDARD", "SOFT", "B30F", "B31F", "B20F", "B30S", "B35F", "FC01", "FC02", "FC10", "C", "B", "STND"}
YEAR_OF_ROUND = {0: 1999, 1: 2000, 2: 2001}


def table(name):
    return pd.read_parquet(glob.glob(os.path.join(CLIN, name, "*"))[0])


def parse_desc(s):
    f = str(s).split(",")
    if len(f) < 7:
        return None
    try:
        thick = float(f[6])
    except ValueError:
        thick = np.nan
    return dict(screen=f[0], img_type=f[1], vendor=f[2], kernel=f[4].upper(), thickness=thick)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--rounds", choices=["last", "all"], default="last")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from idc_index import IDCClient
    c = IDCClient()
    prsn, canc = table("nlst_prsn"), table("nlst_canc")
    for t in (prsn, canc):
        t["pid"] = t["pid"].astype(str)
        for col in ["age", "candx_days", "scr_days0", "scr_days1", "scr_days2"]:
            if col in t:
                t[col] = pd.to_numeric(t[col], errors="coerce")
    nl = c.index[c.index.collection_id == "nlst"].copy()
    nl["PatientID"] = nl["PatientID"].astype(str)
    first = canc[canc.first_lc == 1]
    cand = first[first.path_t.astype(str).isin(["100", "200", "300", "400"])]
    cand = cand.merge(prsn[["pid", "age", "gender", "scr_days0", "scr_days1", "scr_days2"]], on="pid", how="left")
    cand = cand[cand.pid.isin(set(nl.PatientID))]
    cand["age_at_dx"] = cand["age"] + cand["candx_days"] / 365.25

    ser = nl[nl.Modality == "CT"].copy()
    parsed = ser.SeriesDescription.map(parse_desc)
    ser = ser[parsed.notna()].copy()
    ser = pd.concat([ser, pd.DataFrame(list(parsed.dropna()), index=ser.index)], axis=1)
    ser = ser[(ser.img_type == "OPA") & (ser.instanceCount >= 40)]
    ser["year"] = pd.to_datetime(ser.StudyDate).dt.year
    ser["soft"] = ser.kernel.isin(SOFT)
    ser["rank"] = ser.soft.astype(int) * 10 + (ser.thickness <= 3).astype(int) * 5 - ser.thickness.fillna(9) / 10

    rows = []
    for _, r in cand.iterrows():
        s = ser[ser.PatientID == r.pid]
        rounds = []
        for k in (0, 1, 2):
            d = r.get(f"scr_days{k}")
            if pd.notna(d) and d <= r.candx_days:  # screen performed on or before the diagnosis date
                rounds.append(k)
        if not rounds:
            continue
        use = rounds if a.rounds == "all" else [max(rounds)]
        for k in use:
            sk = s[s.year == YEAR_OF_ROUND[k]].sort_values("rank", ascending=False)
            if sk.empty:
                continue
            b = sk.iloc[0]
            rows.append(dict(pid=r.pid, round=k, days_screen_to_dx=r.candx_days - r[f"scr_days{k}"],
                             age_at_dx=round(r.age_at_dx, 1), gender=r.gender, path_t=r.path_t, path_n=r.path_n,
                             path_stag=r.path_stag, SeriesInstanceUID=b.SeriesInstanceUID,
                             SeriesDescription=b.SeriesDescription, kernel=b.kernel, soft=bool(b.soft),
                             thickness=b.thickness, vendor=b.vendor, n_images=int(b.instanceCount),
                             size_MB=round(b.series_size_MB, 1)))
    sel = pd.DataFrame(rows)
    sel.to_csv(os.path.join(a.out, f"nlst_selected_series_{a.rounds}.csv"), index=False)
    print(f"candidates with pT: {len(cand)}; selected series: {len(sel)} from {sel.pid.nunique()} patients; "
          f"aged >=65 at dx: {sel[sel.age_at_dx >= 65].pid.nunique()}; soft kernel {sel.soft.mean():.0%}; "
          f"{sel.size_MB.sum() / 1e3:.1f} GB", file=sys.stderr)
    if a.download:
        dcm = os.path.join(a.out, "dicom")
        os.makedirs(dcm, exist_ok=True)
        c.download_from_selection(seriesInstanceUID=sel.SeriesInstanceUID.tolist(), downloadDir=dcm,
                                  dirTemplate="%PatientID/%SeriesInstanceUID", quiet=True)
        nii = os.path.join(a.out, "nifti")
        os.makedirs(nii, exist_ok=True)
        exe = os.path.join(os.path.dirname(sys.executable), "dcm2niix")
        for _, r in sel.iterrows():
            name = f"{r.pid}_T{r['round']}"
            if os.path.exists(os.path.join(nii, f"{name}.nii.gz")):
                continue
            src = os.path.join(dcm, str(r.pid), r.SeriesInstanceUID)
            subprocess.run([exe, "-z", "y", "-f", name, "-o", nii, src], capture_output=True)
            extra = sorted(glob.glob(os.path.join(nii, f"{name}*.nii.gz")), key=os.path.getsize)
            for e in extra[:-1]:
                os.remove(e)  # keep the largest volume if dcm2niix split the series
            if extra and extra[-1] != os.path.join(nii, f"{name}.nii.gz"):
                os.replace(extra[-1], os.path.join(nii, f"{name}.nii.gz"))
        print("download + conversion finished", file=sys.stderr)


if __name__ == "__main__":
    main()
