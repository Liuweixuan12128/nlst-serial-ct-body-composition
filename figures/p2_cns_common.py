"""Shared helpers for the paper-2 figure set (CNS layout): data access, image crops, small plotting primitives."""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
import style  # noqa: E402
from style import INK, INK2, MUTED, HAIR, WASH, SERIES, TISSUE, MM, fnum  # noqa: E402

RES = os.environ.get("P2_RES", os.path.join(ROOT, "results", "p2"))
style.OUT = os.environ.get("P2_FIG_OUT", os.path.join(ROOT, "figures", "p2", "out"))
NII = os.path.join(ROOT, "data", "nlst_serial", "nifti")
SEG = os.path.join(ROOT, "data", "seg", "nlst_serial")
CASE, CTRL = SERIES[1], SERIES[0]
MEAS = ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu"]
ALL = MEAS + ["lung_perc15_hu"]
LABEL = {"psm_hu": "Muscle radiodensity", "psm_lvol": "Muscle volume", "psm_fat_pct": "Intramuscular fat",
         "T10_trab_hu": "Trabecular bone (T10)", "lung_perc15_hu": "Lung density (15th pct)"}
UNIT = {"psm_hu": "HU", "psm_lvol": "%", "psm_fat_pct": "pp", "T10_trab_hu": "HU", "lung_perc15_hu": "HU"}
PSM = ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"]


def rd(name, res=RES):
    return pd.read_csv(os.path.join(res, name))


def sel(stem):
    return pd.concat([pd.read_csv(os.path.join(ROOT, "data", "nlst_serial", f"{s}_{stem}.csv"),
                                  dtype={"pid": str, "matched_case": str}) for s in ("p2", "p2b")])


def primary():
    return open(os.path.join(RES, "partA_primary_exposure.txt")).read().split()[1]


def axial(cid, level="T8", half_mm=80, shift_mm=20):
    """Native axial slice at a vertebral-body centroid, cropped around the spine, with the paraspinal mask."""
    import nibabel as nib
    from bodycomp_metrics import labels, vertebra_info
    img = nib.load(os.path.join(NII, f"{cid}.nii.gz"))
    ct = np.asanyarray(img.dataobj)
    ct = ct[..., 0] if ct.ndim == 4 else ct
    tot = np.asanyarray(nib.load(os.path.join(SEG, cid, "total.nii.gz")).dataobj)
    am = np.asanyarray(nib.load(os.path.join(SEG, cid, "abdominal_muscles.nii.gz")).dataobj)
    zc = vertebra_info(tot, level)[0]
    if zc is None:
        zc = ct.shape[2] / 2
    z = int(round(zc))
    m = labels(am[:, :, z], "abdominal_muscles", PSM)
    v = tot[:, :, z] == style_inv("total", f"vertebrae_{level}")
    ref = m | v
    px = float(img.header.get_zooms()[0])
    if ref.any():
        ys, xs = np.nonzero(ref)
        cy, cx = int(ys.mean()), int(xs.mean())
    else:
        cy, cx = ct.shape[0] // 2, ct.shape[1] // 2
    h = int(round(half_mm / px))
    sh = int(round(shift_mm / px))
    sl = (slice(max(cy - h, 0), cy + h), slice(max(cx - int(0.72 * h) + sh, 0), cx + int(0.72 * h) + sh))
    return np.rot90(ct[:, :, z][sl]), np.rot90(m[sl]), px


def style_inv(task, name):
    from bodycomp_metrics import INV
    return INV[task][name]


def sagittal(cid, half_w_mm=70):
    """Midline sagittal slice through the spinal canal with vertebra labels T5-T12 and the T6-T11 slab limits."""
    import nibabel as nib
    from bodycomp_metrics import INV, labels, vertebra_info
    img = nib.load(os.path.join(NII, f"{cid}.nii.gz"))
    ct = np.asanyarray(img.dataobj)
    ct = ct[..., 0] if ct.ndim == 4 else ct
    tot = np.asanyarray(nib.load(os.path.join(SEG, cid, "total.nii.gz")).dataobj)
    cord = labels(tot, "total", ["spinal_cord"])
    x = int(np.median(np.nonzero(cord.any(axis=(1, 2)))[0]))
    zooms = img.header.get_zooms()
    yc = int(np.median(np.nonzero(cord.any(axis=(0, 2)))[0]))
    hw = int(round(half_w_mm / zooms[1]))
    ys = slice(max(yc - int(0.6 * hw), 0), yc + int(1.4 * hw))   # more room anterior (+y) than posterior
    sag = ct[x, ys, :]
    lab = tot[x, ys, :]
    lv = {k: vertebra_info(tot, k)[0] for k in ["T6", "T11"]}
    vert = {k: (lab == INV["total"][f"vertebrae_{k}"]) for k in [f"T{i}" for i in range(5, 13)]}
    return sag, vert, lv, zooms


def imshow_ct(ax, img, mask=None, color=None, lw=0.6, window=(-160, 240)):
    ax.imshow(img, cmap="gray", vmin=window[0], vmax=window[1], interpolation="lanczos")
    if mask is not None and mask.any():
        ax.contour(mask, levels=[0.5], colors=color or TISSUE["muscle"], linewidths=lw)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def scalebar(ax, px_mm, length_mm=20, color="white"):
    h, w = ax.images[0].get_array().shape[:2]
    L = length_mm / px_mm
    ax.plot([w * 0.06, w * 0.06 + L], [h * 0.93, h * 0.93], color=color, lw=1.2, solid_capstyle="butt")
    ax.text(w * 0.06 + L / 2, h * 0.9, f"{length_mm} mm", color=color, fontsize=5, ha="center", va="bottom")


def dot_ci(ax, x, y, lo, hi, color, horizontal=True, ms=3.2, mfc=None, marker="o", lw=0.8, zorder=3, label=None):
    if horizontal:
        ax.errorbar(x, y, xerr=np.vstack([np.asarray(x) - lo, np.asarray(hi) - x]), fmt=marker, ms=ms, color=color,
                    mfc=mfc or color, elinewidth=lw, capsize=0, zorder=zorder, label=label)
    else:
        ax.errorbar(x, y, yerr=np.vstack([np.asarray(y) - lo, np.asarray(hi) - y]), fmt=marker, ms=ms, color=color,
                    mfc=mfc or color, elinewidth=lw, capsize=0, zorder=zorder, label=label)


def crop_alpha(path, pad=10):
    from PIL import Image
    im = Image.open(path)
    a = np.asarray(im)[..., 3]
    ys, xs = np.nonzero(a > 8)
    return im.crop((max(xs.min() - pad, 0), max(ys.min() - pad, 0), xs.max() + pad, ys.max() + pad))
