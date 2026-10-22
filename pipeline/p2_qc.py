"""Paper 2 QC flags: scan-level and change-level, then montages for visual review.

Scan level: measurement error, paraspinal slab T6-T11 incomplete, compartment touching the image edge, values outside
the 0.25-99.75th percentile of all scans, and disagreement between the two segmentation models (paraspinal volume vs
TotalSegmentator autochthon volume ratio outside its 0.25-99.75th percentile).
Change level: consecutive-round change in volume (100 x log), radiodensity or fat % larger than 4 robust SDs
(1.4826 x MAD of consecutive changes in controls), i.e. far beyond what normal ageing and noise produce.
Level check (automatic exclusion, objective): the T6-T11 slab length (slices x spacing) is anatomical and should not
change between a participant's rounds; a scan whose slab length differs from the participant's median by >8% has a
vertebral level misidentified and is excluded (qc_pass = 0 in p2_qc_review.csv, reason recorded).
Position check (automatic exclusion): a NIfTI affine with the left-right and anterior-posterior axes both reversed
(dcm2niix output of an image orientation of -1,0,0,0,-1,0) marks a prone acquisition; posture changes paraspinal shape
and attenuation, so these scans are excluded.
Other flags select participants for visual review; they do not exclude anyone by themselves.

Usage: python p2_qc.py --metrics results/p2/p2_paraspinal.csv --out results/p2 [--montages]
"""
import argparse, glob, os, subprocess, sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--montages", action="store_true")
    a = ap.parse_args()
    m = pd.read_csv(a.metrics)
    if "error" not in m:
        m["error"] = np.nan
    ser = pd.concat([pd.read_csv(f, dtype={"pid": str})
                     for f in sorted(glob.glob(os.path.join(ROOT, "data", "nlst_serial", "p2*_series.csv")))])
    d = ser[["case_id", "pid", "group", "round"]].merge(m, on="case_id", how="left")
    d["psm_lvol"] = 100 * np.log(d.psm_vol_cm3)
    d["ratio_aut"] = d.aut_vol_cm3 / d.psm_vol_cm3
    flags = []
    for _, r in d.iterrows():
        f = []
        if pd.isna(r.get("psm_available")):
            f.append("not measured")
        elif isinstance(r.error, str):
            f.append("measurement error")
        else:
            if r.psm_available != 1:
                f.append("T6-T11 slab incomplete")
            elif r.psm_edge == 1:
                f.append("paraspinal at image edge")
        flags.append(f)
    for v in ["psm_hu", "psm_lvol", "psm_fat_pct", "ratio_aut"]:
        lo, hi = d[v].quantile([0.0025, 0.9975])
        for i in np.where((d[v] < lo) | (d[v] > hi))[0]:
            flags[i].append(f"{v} extreme ({d[v].iloc[i]:.1f})")
    d["scan_flags"] = ["; ".join(f) for f in flags]
    d["slab_mm"] = d.psm_nslices * d.slice_mm
    d["slab_dev"] = (d.slab_mm / d.pid.map(d[d.scan_flags == ""].groupby("pid").slab_mm.median()) - 1).abs()
    d["level_flag"] = (d.slab_dev > 0.08) & (d.scan_flags == "")

    # change level
    ok = d[(d.scan_flags.str.len() == 0) & ~d.level_flag].sort_values(["pid", "round"])
    ch = []
    for v in ["psm_lvol", "psm_hu", "psm_fat_pct"]:
        x = ok.groupby("pid")[v].diff()
        ok[f"d_{v}"] = x
    ctrl = ok[ok.group == "control"]
    thr = {v: 4 * 1.4826 * (ctrl[f"d_{v}"] - ctrl[f"d_{v}"].median()).abs().median()
           for v in ["psm_lvol", "psm_hu", "psm_fat_pct"]}
    for _, r in ok.iterrows():
        f = [f"change {v} {r[f'd_{v}']:+.1f} (>|{thr[v]:.1f}|)" for v in thr
             if pd.notna(r[f"d_{v}"]) and abs(r[f"d_{v}"]) > thr[v]]
        if f:
            ch.append(dict(case_id=r.case_id, pid=r.pid, change_flags="; ".join(f)))
    ch = pd.DataFrame(ch, columns=["case_id", "pid", "change_flags"])
    out = d[["case_id", "pid", "group", "round", "scan_flags", "slab_mm", "slab_dev", "level_flag"]].merge(
        ch[["case_id", "change_flags"]], on="case_id", how="left")
    out["change_flags"] = out.change_flags.fillna("")
    os.makedirs(a.out, exist_ok=True)
    out.to_csv(os.path.join(a.out, "p2_qc_flags.csv"), index=False)
    import nibabel as nib
    nii = os.path.join(ROOT, "data", "nlst_serial", "nifti")
    prone = [c for c in out.case_id
             if (lambda a: a[0, 0] > 0 and a[1, 1] < 0)(nib.load(os.path.join(nii, f"{c}.nii.gz")).affine)]
    auto = out[out.level_flag].assign(qc_pass=0, reviewer="automatic rule",
                                      reason=lambda x: "T6-T11 slab length differs from participant median by "
                                      + (100 * x.slab_dev).round(1).astype(str) + "% (vertebral level misidentified)")
    auto = pd.concat([auto, pd.DataFrame({"case_id": prone, "qc_pass": 0, "reviewer": "automatic rule",
                                          "reason": "prone acquisition (affine left-right and anterior-posterior axes "
                                                    "reversed)"})]).drop_duplicates("case_id")
    print(f"prone acquisitions excluded: {len(prone)}", file=sys.stderr)
    auto[["case_id", "qc_pass", "reviewer", "reason"]].to_csv(os.path.join(a.out, "p2_qc_auto_exclusions.csv"),
                                                              index=False)
    pd.Series(thr).rename("threshold_4_robust_sd").to_csv(os.path.join(a.out, "p2_qc_change_thresholds.csv"))
    rev = out[((out.scan_flags != "") & ~out.scan_flags.str.contains("not measured|slab incomplete|image edge"))
              | (out.change_flags != "")]
    print(f"scans {len(out)}; scan-level flags {(out.scan_flags != '').sum()}; level-check exclusions "
          f"{int(out.level_flag.sum())}; change flags {(out.change_flags != '').sum()}; "
          f"participants to review {rev.pid.nunique()}\nthresholds {thr}\n"
          + out.scan_flags.str.split('; ').explode().replace('', np.nan).dropna()
          .str.replace(r' \(.*', '', regex=True).value_counts().to_string(), file=sys.stderr)
    if a.montages:
        pids = sorted(rev.pid.unique())
        for i in range(0, len(pids), 40):
            subprocess.run([sys.executable, os.path.join(ROOT, "pipeline", "p2_qc_montage.py"), "--pids",
                            ",".join(pids[i:i + 40]), "--nii-dir", os.path.join(ROOT, "data", "nlst_serial", "nifti"),
                            "--seg-dir", os.path.join(ROOT, "data", "seg", "nlst_serial"),
                            "--out", os.path.join(a.out, "qc_montages")], check=False)


if __name__ == "__main__":
    main()
