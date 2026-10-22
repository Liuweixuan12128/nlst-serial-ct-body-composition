"""Paper 2 QC: one PNG per participant, rounds side by side, paraspinal mask at T8 and T10 (cropped to the spine).

Used to check extreme within-person changes (|MAA z| > 3 on any candidate measure) and a random sample, so that a
segmentation failure in one round cannot pass as "accelerated ageing".

Usage: python p2_qc_montage.py --pids 100012,100147 --nii-dir data/nlst_serial/nifti --seg-dir data/seg/nlst_serial
                               --out results/p2/qc
"""
import argparse, glob, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bodycomp_metrics import INV, labels, vertebra_info  # noqa: E402

PSM = ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"]


def panel(nii, seg_dir, level):
    ct = np.asanyarray(nib.load(nii).dataobj)
    if ct.ndim == 4:
        ct = ct[..., 0]
    tot = np.asanyarray(nib.load(os.path.join(seg_dir, "total.nii.gz")).dataobj)
    am = np.asanyarray(nib.load(os.path.join(seg_dir, "abdominal_muscles.nii.gz")).dataobj)
    zc, ok = vertebra_info(tot, level)
    if zc is None:
        return None
    z = int(round(zc))
    m = labels(am[:, :, z], "abdominal_muscles", PSM)
    v = tot[:, :, z] == INV["total"][f"vertebrae_{level}"]
    ys, xs = np.nonzero(m | v) if (m | v).any() else (np.array([ct.shape[0] // 2]), np.array([ct.shape[1] // 2]))
    cy, cx, h = int(ys.mean()), int(xs.mean()), 140
    sl = (slice(max(cy - h, 0), cy + h), slice(max(cx - h, 0), cx + h))
    return np.rot90(ct[:, :, z][sl]), np.rot90(m[sl]), ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pids", required=True)
    ap.add_argument("--nii-dir", required=True)
    ap.add_argument("--seg-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="", help="text shown in the title, e.g. the measure and z")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for pid in a.pids.split(","):
        scans = sorted(glob.glob(os.path.join(a.nii_dir, f"{pid}_T*.nii.gz")))
        if not scans:
            continue
        fig, axes = plt.subplots(2, len(scans), figsize=(2.6 * len(scans), 5.4), squeeze=False)
        for j, f in enumerate(scans):
            cid = os.path.basename(f)[:-7]
            for i, lv in enumerate(("T8", "T10")):
                ax = axes[i, j]
                ax.axis("off")
                try:
                    p = panel(f, os.path.join(a.seg_dir, cid), lv)
                except Exception as e:  # missing mask etc.
                    ax.set_title(f"{cid} {lv}: {type(e).__name__}", fontsize=7)
                    continue
                if p is None:
                    ax.set_title(f"{cid} {lv}: level missing", fontsize=7)
                    continue
                img, m, ok = p
                ax.imshow(img, cmap="gray", vmin=-160, vmax=240)
                ax.contour(m, levels=[0.5], colors="#e8590c", linewidths=0.6)
                ax.set_title(f"{cid} {lv}" + ("" if ok else " (incomplete)"), fontsize=7)
        fig.suptitle(f"{pid} {a.tag}", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(a.out, f"{pid}.png"), dpi=110)
        plt.close(fig)


if __name__ == "__main__":
    main()
