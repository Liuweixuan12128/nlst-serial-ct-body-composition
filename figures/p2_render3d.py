"""Paper 2 Figure 1 3D renders of one screening CT (VTK, GPU offscreen via EGL), reusing figures/render3d.py.

Views of the same scan, identical scale:
  spine   posterior-oblique: thoracic skeleton; paraspinal compartment (erector spinae + transversospinalis) inside the
          measured T6-T11 slab in tissue red, outside the slab translucent neutral
  fat     the T6-T11 compartment alone, translucent red, with the voxels in the fat window (-190 to -30 HU on the
          harmonised image) rendered solid violet: what the primary measure (intramuscular fat %) counts
Writes figures/assets/p2_render3d_<case>_<view>.png (RGBA) and a JSON with T6/T8/T11 screen positions.

  env -u DISPLAY PYTHONPATH=env/vtk97 VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow \
      env/bin/python figures/p2_render3d.py 100658_T0
"""
import json
import os
import sys

import numpy as np
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
sys.path.insert(0, HERE)
from bodycomp_metrics import CAGE, HU_FAT, Case, in_range, labels, vertebra_info  # noqa: E402
from render3d import actor, render, surface  # noqa: E402

RED, BONE, NEUTRAL, VIOLET = "#e34948", "#f2ecdc", "#b9a7a3", "#4a3aa7"
OUT = os.path.join(ROOT, "figures", "assets")
PSM = ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"]


def main(cid):
    os.makedirs(OUT, exist_ok=True)
    c = Case(os.path.join(ROOT, "data", "nlst_serial", "nifti", f"{cid}.nii.gz"),
             os.path.join(ROOT, "data", "seg", "nlst_serial", cid))
    tot, am = c.segs["total"], c.segs["abdominal_muscles"]
    info = {lv: vertebra_info(tot, lv)[0] for lv in ["T3", "T5", "T6", "T8", "T11", "T12"]}
    zz = c.zooms[2]
    z11, z6 = int(round(info["T11"])), int(round(info["T6"]))
    zlo = max(int(min(info["T12"], info["T11"]) - 25 / zz), 0)
    zhi = min(int((info["T3"] or info["T5"]) + 15 / zz), tot.shape[2] - 1)
    sl = (slice(None), slice(None), slice(zlo, zhi + 1))
    bone = labels(tot, "total", CAGE)[sl]
    psm = labels(am, "abdominal_muscles", PSM)[sl]
    span = np.zeros(psm.shape[2], bool)
    span[z11 - zlo:z6 - zlo + 1] = True
    psm_in, psm_out = psm & span[None, None, :], psm & ~span[None, None, :]
    fat = psm_in & in_range(c.ct[sl], HU_FAT)
    fat = ndi.binary_opening(fat, structure=np.ones((2, 2, 1), bool))
    lab, n = ndi.label(fat)
    sz = ndi.sum(fat, lab, index=np.arange(1, n + 1))
    fat = np.isin(lab, 1 + np.nonzero(sz >= 8)[0])
    zooms = c.zooms
    print("surfaces ...", flush=True)
    S = {"bone": surface(bone, zooms, 0.8), "in": surface(psm_in, zooms, 0.9), "out": surface(psm_out, zooms, 0.9),
         "fat": surface(fat, zooms, 0.5, decimate=0.2),
         "bone_in": surface(bone & span[None, None, :], zooms, 0.8)}
    xs, ys = np.nonzero(psm.any(axis=2))
    centre = [(xs.min() + xs.max()) / 2 * zooms[0], (ys.min() + ys.max()) / 2 * zooms[1], (psm.shape[2] - 1) / 2 * zz]
    cord = labels(tot, "total", ["spinal_cord"])[sl]
    sx = np.median(np.nonzero(cord.any(axis=(1, 2)))[0]) * zooms[0]
    sy = np.median(np.nonzero(cord.any(axis=(0, 2)))[0]) * zooms[1]
    levels = {lv: [sx, sy, (info[lv] - zlo) * zz] for lv in ["T6", "T8", "T11"]}
    views = {
        "spine": ([actor(S["bone"], BONE, spec=0.25), actor(S["in"], RED), actor(S["out"], NEUTRAL, opacity=0.4)],
                  -30, 12, (1300, 1900), 1.25),
        "fat": ([actor(S["in"], RED, opacity=0.18, spec=0.05), actor(S["fat"], VIOLET, spec=0.2)],
                -30, 12, (1300, 1900), 1.25),
    }
    meta = {"case": cid, "views": {}}
    for name, (acts, az, el, size, zoom) in views.items():
        fn = os.path.join(OUT, f"p2_render3d_{cid}_{name}.png")
        to_px, scale = render(acts, centre, az, el, size, zoom, fn)
        meta["views"][name] = {"file": os.path.relpath(fn, ROOT), "size_px": size, "mm_per_px": 2 * scale / size[1],
                               "levels_px": {lv: to_px(p) for lv, p in levels.items()}}
        print("wrote", fn, flush=True)
    meta["fat_voxels"] = int(fat.sum())
    meta["psm_voxels_in_slab"] = int(psm_in.sum())
    json.dump(meta, open(os.path.join(OUT, f"p2_render3d_{cid}.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "100658_T0")
