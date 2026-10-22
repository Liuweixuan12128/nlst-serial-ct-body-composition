"""Paper 2 Figure 2a: 3D rendering of the three reserve measures on one screening CT, plus small 3D icons of each tissue.

Hero view: lungs as translucent glass with the darkest 15% of lung voxels (the lung-density measure) rendered solid,
thoracic vertebrae with T10 highlighted (the trabecular-bone measure), and the T6-T11 paraspinal compartment (the
muscle measures). Icons: each structure alone, same camera, for labels across the figure set.
Writes figures/assets/p2_reserve3d_<view>.png (RGBA) and p2_icon_<tissue>.png.

  env -u DISPLAY PYTHONPATH=env/vtk97 VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow \
      env/bin/python figures/p2_render_reserve3d.py 100658_T0
"""
import os
import sys

import numpy as np
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
sys.path.insert(0, HERE)
from bodycomp_metrics import CAGE, HU_FAT, INV, LUNGS, VERTEBRAE, Case, in_range, labels, vertebra_info  # noqa: E402
from render3d import actor, render, surface  # noqa: E402

BONE_COL, MUSCLE_COL, FAT_COL, LUNG_COL = "#d99000", "#d63c3c", "#5b45b0", "#11968f"
IVORY, LUNG_GLASS, LUNG_DARK = "#efe9da", "#9fd3cf", "#0c6f6a"
OUT = os.path.join(ROOT, "figures", "assets")
VIEWS = [(-75, "left"), (-75, "right"), (75, "left"), (75, "right")]
PSM = ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"]


def clean(mask, min_vox):
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    sz = ndi.sum(mask, lab, index=np.arange(1, n + 1))
    return np.isin(lab, 1 + np.nonzero(sz >= min_vox)[0])


def main(cid):
    os.makedirs(OUT, exist_ok=True)
    c = Case(os.path.join(ROOT, "data", "nlst_serial", "nifti", f"{cid}.nii.gz"),
             os.path.join(ROOT, "data", "seg", "nlst_serial", cid))
    tot, am = c.segs["total"], c.segs["abdominal_muscles"]
    zz = c.zooms[2]
    info = {lv: vertebra_info(tot, lv)[0] for lv in ["T1", "T3", "T6", "T10", "T11", "T12", "L1"]}
    lungs_full = labels(tot, "total", LUNGS)
    lz = np.nonzero(lungs_full.any(axis=(0, 1)))[0]
    zlo = max(int(min(lz.min(), (info["L1"] or info["T12"]))) - 3, 0)
    zhi = min(int(lz.max()) + 3, tot.shape[2] - 1)
    sl = (slice(None), slice(None), slice(zlo, zhi + 1))
    lungs = lungs_full[sl]
    p15 = float(np.percentile(c.ct[lungs_full], 15))
    low = lungs & (c.ct[sl] < p15)
    # the darkest voxels are scattered; render where they cluster (smoothed density > 0.5 within the lungs)
    dens = ndi.gaussian_filter(low.astype(np.float32), sigma=[1.2 / c.zooms[0], 1.2 / c.zooms[1], max(1.2 / zz, 0.6)])
    low_r = clean((dens > 0.55) & lungs, 30)
    vert = labels(tot, "total", VERTEBRAE)[sl]
    ribs = labels(tot, "total", [x for x in CAGE if x not in VERTEBRAE])[sl]
    t10 = (tot == INV["total"]["vertebrae_T10"])[sl]
    z11, z6 = int(round(info["T11"])), int(round(info["T6"]))
    span = np.zeros(lungs.shape[2], bool)
    span[z11 - zlo:z6 - zlo + 1] = True
    psm = labels(am, "abdominal_muscles", PSM)[sl] & span[None, None, :]
    fat = clean(ndi.binary_opening(psm & in_range(c.ct[sl], HU_FAT), structure=np.ones((2, 2, 1), bool)), 8)
    zooms = c.zooms
    print("surfaces ...", flush=True)
    S = {"lung": surface(lungs, zooms, 1.6, decimate=0.6), "low": surface(low_r, zooms, 0.6, decimate=0.5),
         "vert": surface(vert & ~t10, zooms, 0.8), "t10": surface(t10, zooms, 0.8),
         "ribs": surface(ribs, zooms, 0.8, decimate=0.5), "psm": surface(psm, zooms, 0.9),
         "fat": surface(fat, zooms, 0.5, decimate=0.2)}
    xs, ys = np.nonzero(lungs.any(axis=2))
    centre = [(xs.min() + xs.max()) / 2 * zooms[0], (ys.min() + ys.max()) / 2 * zooms[1], (lungs.shape[2] - 1) / 2 * zz]
    xs2, ys2 = np.nonzero(psm.any(axis=2))
    zp = np.nonzero(psm.any(axis=(0, 1)))[0]
    c_psm = [(xs2.min() + xs2.max()) / 2 * zooms[0], (ys2.min() + ys2.max()) / 2 * zooms[1], zp.mean() * zz]
    zt = np.nonzero(t10.any(axis=(0, 1)))[0]
    xt, yt = np.nonzero(t10.any(axis=2))
    c_t10 = [xt.mean() * zooms[0], yt.mean() * zooms[1], zt.mean() * zz]
    # cut-away hero views: the lung (and ribs) nearer the viewer are removed so the spine and T10 stay visible
    side = {k: labels(tot, "total", v)[sl] for k, v in
            {"left": ["lung_upper_lobe_left", "lung_lower_lobe_left"],
             "right": ["lung_upper_lobe_right", "lung_middle_lobe_right", "lung_lower_lobe_right"]}.items()}
    rib_side = {k: labels(tot, "total", [f"rib_{k}_{i}" for i in range(1, 13)])[sl] for k in ("left", "right")}
    for az, keep in VIEWS:
        acts = [actor(surface(side[keep], zooms, 1.6, decimate=0.6), LUNG_GLASS, opacity=0.16, spec=0.4, amb=0.25),
                actor(surface(low_r & side[keep], zooms, 0.6, decimate=0.5), LUNG_DARK, opacity=0.92, spec=0.2),
                actor(surface(rib_side[keep], zooms, 0.8, decimate=0.5), IVORY, opacity=0.35, spec=0.2),
                actor(S["vert"], IVORY, spec=0.25), actor(S["t10"], BONE_COL, spec=0.3),
                actor(S["psm"], MUSCLE_COL, spec=0.15)]
        fn = os.path.join(OUT, f"p2_reserve3d_az{az}_{keep}.png")
        render(acts, centre, az, 8, (1500, 1800), 1.15, fn)
        print("wrote", fn, flush=True)
    icons = {"bone": ([actor(S["t10"], BONE_COL, spec=0.3)], c_t10, -35, 20),
             "muscle": ([actor(S["psm"], MUSCLE_COL, spec=0.15)], c_psm, -35, 12),
             "fat": ([actor(S["psm"], MUSCLE_COL, opacity=0.3, spec=0.05), actor(S["fat"], FAT_COL, spec=0.2)],
                     c_psm, -35, 12),
             "lung": ([actor(S["lung"], LUNG_GLASS, opacity=0.22, spec=0.4, amb=0.25),
                       actor(S["low"], LUNG_DARK, opacity=0.9, spec=0.2)], centre, 0, 0)}
    for k, (acts, cc, az, el) in icons.items():
        fn = os.path.join(OUT, f"p2_icon_{k}.png")
        render(acts, cc, az, el, (700, 700) if k != "muscle" and k != "fat" else (500, 900), 1.0, fn, ssao=k != "lung")
        print("wrote", fn, flush=True)
    print("p15", p15, "low voxels", int(low.sum()), "rendered", int(low_r.sum()))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "100658_T0")
