"""Paper 2: shape-normalised intramuscular-fat maps of the paraspinal compartment, one per QC-passing scan.

For every axial slice between the T11 and T6 vertebral-body centroids, each side's paraspinal compartment (erector
spinae + transversospinalis) is mapped to a normalised frame: u = medial (0) -> lateral (1), from the midline (spinal
cord centroid) to the compartment's lateral edge on that slice; v = anterior (0) -> posterior (1) across the
compartment's extent on that slice. Left and right are mirrored onto one frame. Voxels are binned on a 24 (v) x 32 (u)
grid and accumulated over the slab:
  cnt   compartment voxels per cell          fat   voxels at -190..-30 HU (harmonised image, as in p2_metrics)
  hu    sum of harmonised HU                 lev_cnt / lev_fat   the same totals in 6 craniocaudal bins (T11 -> T6)
Harmonisation is the one used for every measure (5-mm slab mean along z, 0.75-mm in-plane Gaussian), applied to a crop.

Usage: python p2_maps.py --out results/p2/maps [--shard i/n] [--workers 44]
"""
import argparse, os, sys
from multiprocessing import Pool

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from bodycomp_metrics import HU_FAT, SLAB_MM, SMOOTH_MM, labels, vertebra_info  # noqa: E402

NV, NU, NL = 24, 32, 6
SIDES = {"left": ["erector_spinae_left", "transversospinalis_left"],
         "right": ["erector_spinae_right", "transversospinalis_right"]}


def one(args):
    cid, out = args
    dst = os.path.join(out, f"{cid}.npz")
    if os.path.exists(dst):
        return cid, "exists"
    try:
        img = nib.load(os.path.join(ROOT, "data", "nlst_serial", "nifti", f"{cid}.nii.gz"))
        zooms = [float(z) for z in img.header.get_zooms()[:3]]
        seg = os.path.join(ROOT, "data", "seg", "nlst_serial", cid)
        tot = np.asanyarray(nib.load(os.path.join(seg, "total.nii.gz")).dataobj)
        am = np.asanyarray(nib.load(os.path.join(seg, "abdominal_muscles.nii.gz")).dataobj)
        zt, _ = vertebra_info(tot, "T6")
        zb, _ = vertebra_info(tot, "T11")
        z0, z1 = sorted((int(round(zb)), int(round(zt))))
        side_m = {s: labels(am, "abdominal_muscles", n) for s, n in SIDES.items()}
        anym = side_m["left"] | side_m["right"]
        sub = anym[:, :, z0:z1 + 1]
        xs, ys = np.nonzero(sub.any(axis=2))
        k = max(1, int(round(SLAB_MM / zooms[2])))
        pad = 6
        xa, xb = max(xs.min() - pad, 0), min(xs.max() + pad + 1, anym.shape[0])
        ya, yb = max(ys.min() - pad, 0), min(ys.max() + pad + 1, anym.shape[1])
        za, zbb = max(z0 - k, 0), min(z1 + k + 1, anym.shape[2])
        ct = np.asanyarray(img.dataobj)
        ct = (ct[..., 0] if ct.ndim == 4 else ct)[xa:xb, ya:yb, za:zbb].astype(np.float32)
        h = ndi.uniform_filter1d(ct, k, axis=2, mode="nearest") if k > 1 else ct
        h = ndi.gaussian_filter(h, sigma=(SMOOTH_MM / zooms[0], SMOOTH_MM / zooms[1], 0))
        cord = labels(tot, "total", ["spinal_cord"])[xa:xb, ya:yb, za:zbb]
        cnt = np.zeros((NV, NU)); fat = np.zeros((NV, NU)); hu = np.zeros((NV, NU))
        lev_cnt = np.zeros(NL); lev_fat = np.zeros(NL)
        for z in range(z0, z1 + 1):
            zl = z - za
            img_s = h[:, :, zl]
            c = cord[:, :, zl]
            xm = np.nonzero(c.any(axis=1))[0].mean() if c.any() else None
            lev = min(int((z - z0) / max(z1 - z0, 1) * NL), NL - 1)
            for s, mfull in side_m.items():
                m = mfull[xa:xb, ya:yb, z]
                if m.sum() < 30:
                    continue
                px, py = np.nonzero(m)
                mid = xm if xm is not None else (np.nonzero(anym[xa:xb, ya:yb, z])[0].mean())
                d = np.abs(px - mid)
                if d.max() - d.min() < 3 or py.max() - py.min() < 3:
                    continue
                u = (d - d.min()) / (d.max() - d.min())
                v = (py.max() - py) / (py.max() - py.min())          # +axis1 = anterior -> v = 0 anterior
                iu = np.minimum((u * NU).astype(int), NU - 1)
                iv = np.minimum((v * NV).astype(int), NV - 1)
                vals = img_s[px, py]
                f = (vals >= HU_FAT[0]) & (vals <= HU_FAT[1])
                np.add.at(cnt, (iv, iu), 1)
                np.add.at(fat, (iv, iu), f)
                np.add.at(hu, (iv, iu), vals)
                lev_cnt[lev] += len(vals)
                lev_fat[lev] += f.sum()
        np.savez_compressed(dst, cnt=cnt.astype(np.float32), fat=fat.astype(np.float32), hu=hu.astype(np.float32),
                            lev_cnt=lev_cnt, lev_fat=lev_fat)
        return cid, "ok"
    except Exception as e:
        return cid, f"error {type(e).__name__}: {str(e)[:120]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--workers", type=int, default=44)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ids = pd.read_csv(os.path.join(ROOT, "results", "p2", "partA_long.csv"), usecols=["case_id"]).case_id.tolist()
    i, n = map(int, a.shard.split("/"))
    ids = sorted(ids)[i::n]
    with Pool(a.workers) as p:
        res = p.map(one, [(c, a.out) for c in ids], chunksize=4)
    bad = [r for r in res if not r[1] in ("ok", "exists")]
    print(f"shard {a.shard}: {len(res)} scans, {len(bad)} errors", bad[:5], file=sys.stderr)


if __name__ == "__main__":
    main()
