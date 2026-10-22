"""Paper 2: truncation-robust paraspinal muscle measures for within-person change on serial screening CT.

Why a separate measure: NLST screening scans are reconstructed with a lung field of view, so the body outline is cut
at T8 in 89% and at T12 in 65% of scans (paper 1), and the field of view changes between rounds. Whole-trunk muscle and
fat volumes would therefore change for technical reasons. The paraspinal compartment (erector spinae +
transversospinalis) lies next to the spine and is rarely cut; we measure it over whole vertebral-body levels and flag
any slice where it touches the image edge.

Per scan, over the slab from the T11 to the T6 vertebral-body centroid (primary; present in >99% of NLST scans) and
from T12 to T5 (secondary), using the harmonised image of bodycomp_metrics.Case (5 mm slab, 0.75 mm in-plane
Gaussian):
  psm_vol_cm3        paraspinal compartment volume (anatomical mask, no HU window)
  psm_lean_cm3       voxels at -29..150 HU inside the compartment
  psm_fat_pct        100 x voxels at -190..-30 HU / compartment voxels (intramuscular fat fraction)
  psm_hu             mean HU of lean-window voxels (harmonised); psm_hu_raw on the native image
  psm_hu_all         mean HU of every compartment voxel (no window)
  psm_edge           1 if the compartment touches the image edge on any slab slice
The same measures from TotalSegmentator `total` autochthon (an independent model) are written with prefix `aut_`.
Also: per-level paraspinal area/HU at T8 and T10, trabecular HU at T8/T10/T12, liver HU SD (noise proxy), and lung
volume, %LAA-950 and 15th-percentile lung HU (emphysema proxy).

Usage: python p2_metrics.py --nii-dir data/nlst_serial/nifti --seg-dir data/seg/nlst_serial
                            --out results/p2/p2_paraspinal.csv [--workers 24]
"""
import argparse, os, sys, traceback
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bodycomp_metrics import HU_FAT, HU_MUSCLE, Case, in_range, labels, vertebra_info  # noqa: E402

PSM = ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"]
SLABS = {"": ("T6", "T11"), "_t5t12": ("T5", "T12")}


def slab_measures(c, mask, top, bottom, prefix):
    zt, ok_t = vertebra_info(c.segs["total"], top)
    zb, ok_b = vertebra_info(c.segs["total"], bottom)
    if zt is None or zb is None or not (ok_t and ok_b):
        return {f"{prefix}_available": 0}
    z0, z1 = sorted((int(round(zb)), int(round(zt))))
    m = mask[:, :, z0:z1 + 1]
    img, raw = c.ct[:, :, z0:z1 + 1], c.raw[:, :, z0:z1 + 1]
    lean = m & in_range(img, HU_MUSCLE)
    fat = m & in_range(img, HU_FAT)
    n = int(m.sum())
    edge = bool(m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any())
    return {f"{prefix}_available": 1, f"{prefix}_nslices": z1 - z0 + 1,
            f"{prefix}_vol_cm3": n * c.vox_cm3, f"{prefix}_lean_cm3": float(lean.sum()) * c.vox_cm3,
            f"{prefix}_fat_pct": 100.0 * fat.sum() / max(n, 1),
            f"{prefix}_hu": float(img[lean].mean()) if lean.sum() > 200 else np.nan,
            f"{prefix}_hu_raw": float(raw[lean].mean()) if lean.sum() > 200 else np.nan,
            f"{prefix}_hu_all": float(img[m].mean()) if n > 200 else np.nan,
            f"{prefix}_edge": int(edge)}


def one(args):
    cid, nii, seg = args
    try:
        c = Case(nii, seg, harmonize=True)
        aut = labels(c.segs["total"], "total", ["autochthon_left", "autochthon_right"])
        row = {"case_id": cid, "px_mm": c.zooms[0], "slice_mm": c.zooms[2]}
        for suffix, (top, bottom) in SLABS.items():
            row.update(slab_measures(c, c.esm, top, bottom, "psm" + suffix))
            row.update(slab_measures(c, aut, top, bottom, "aut" + suffix))
        for lv in ("T8", "T10", "T12"):
            zc, ok = vertebra_info(c.segs["total"], lv)
            if zc is None or not ok:
                continue
            z = int(round(zc))
            if lv != "T12":
                m = c.esm[:, :, z]
                lean = m & in_range(c.ct[:, :, z], HU_MUSCLE)
                row[f"{lv}_psm_area"] = float(m.sum()) * c.px_cm2
                row[f"{lv}_psm_hu"] = float(c.ct[:, :, z][lean].mean()) if lean.sum() > 50 else np.nan
            row[f"{lv}_trab_hu"] = c.trabecular_hu(lv, z)
        row.update({k: v for k, v in c.organs().items() if k in ("liver_hu", "liver_sd_raw", "liver_sd_harmonised")})
        # lungs (emphysema proxy; public tables have no pack-years): volume, % voxels < -950 HU and 15th percentile
        # HU on the harmonised image
        if c.lungs.sum() > 10000:
            v = c.ct[c.lungs]
            row.update(lung_vol_l=float(c.lungs.sum()) * c.vox_cm3 / 1000, lung_laa950_pct=float(100 * np.mean(v < -950)),
                       lung_perc15_hu=float(np.percentile(v, 15)))
        return row
    except Exception:
        return {"case_id": cid, "error": traceback.format_exc()[-500:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nii-dir", required=True)
    ap.add_argument("--seg-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--only", default="", help="comma-separated case ids (testing)")
    ap.add_argument("--shard", default="0/1", help="i/n: measure every n-th pending scan starting at i")
    ap.add_argument("--skip", default="", help="CSV of scans already measured elsewhere (case_id column)")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    done, prev = set(), pd.DataFrame(columns=["case_id"])
    if os.path.exists(a.out):
        prev = pd.read_csv(a.out)
        done = set(prev.loc[prev["error"].isna() if "error" in prev else slice(None), "case_id"])
    if a.skip and os.path.exists(a.skip):
        sk = pd.read_csv(a.skip)
        done |= set(sk.loc[sk["error"].isna() if "error" in sk else slice(None), "case_id"])
    jobs = []
    for f in sorted(os.listdir(a.nii_dir)):
        if not f.endswith(".nii.gz"):
            continue
        cid = f[:-7]
        if (a.only and cid not in a.only.split(",")) or cid in done:
            continue
        seg = os.path.join(a.seg_dir, cid)
        if all(os.path.exists(os.path.join(seg, f"{t}.nii.gz")) for t in ("total", "abdominal_muscles", "body")):
            jobs.append((cid, os.path.join(a.nii_dir, f), seg))
    si, sn = map(int, a.shard.split("/"))
    jobs = jobs[si::sn]
    with Pool(a.workers) as p:
        rows = p.map(one, jobs, chunksize=1)
    new = pd.DataFrame(rows)
    keep = prev[prev.case_id.isin(done)]
    if len(keep):
        new = pd.concat([keep, new], ignore_index=True)
    new.sort_values("case_id").to_csv(a.out, index=False)
    err = new["error"].notna().sum() if "error" in new else 0
    print(f"measured {len(rows)} new scans; total {len(new)}; errors {err}", file=sys.stderr)


if __name__ == "__main__":
    main()
