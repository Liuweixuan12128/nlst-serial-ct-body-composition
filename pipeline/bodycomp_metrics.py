"""Extract body-composition metrics from a CT + TotalSegmentator masks.

Arrays are reoriented to RAS+ (axis0 -> patient Right, axis1 -> Anterior, axis2 -> Superior).

Noise harmonisation. Staging CTs come with thin sharp-kernel (lung/bone) and 5-mm
soft-kernel reconstructions; HU-threshold areas are noise dependent. Every CT is
therefore converted to a 5-mm moving-average volume and smoothed in-plane with a
Gaussian (sigma 0.75 mm) before thresholding. On a 1.25-mm BONEPLUS series this
brings liver noise from SD 98 HU to ~12 HU (the level of a 5-mm soft-kernel image).

Per vertebral level (T4, T8, T10, T12, L1, L3; only when the vertebra lies wholly
inside the scan) on the harmonised slice through the vertebral centroid:
  sma       skeletal-muscle area, cm2 (DL muscle masks, HU -29..150)
  smd       skeletal-muscle radiodensity, mean HU (lower = myosteatosis) on the harmonised image
  smd_raw   same voxels, native-image HU (no slab partial volume); choice between the two is made on
            scan-rescan reproducibility in the public cohort
  lama/nama low (-29..29 HU) / normal (30..150 HU) attenuation muscle area
  imat      intramuscular fat area (fat HU inside the muscle masks)
  sat       subcutaneous fat: outside the convex hull of muscles + thoracic cage
  mwf       muscle-wall fat: between that outer hull and the inner cage (inter-muscular,
            paravertebral, axillary fat)
  ift       internal fat: inside the convex hull of cage + internal organs, lungs
            (dilated 3 mm) and GI-lumen masks excluded -> mediastinal/epicardial fat on thoracic levels,
            visceral fat on T12-L3
  esm_*     erector spinae + transversospinalis;  pm_* pectoralis major
  trab_hu   vertebral-body trabecular HU (opportunistic BMD), T12 and L1
  fov_trunc 1 when the trunk touches the image edge (SAT/MWF not valid on that level)
Volumetric (T5 centroid -> T12 centroid): muscle volume and HU, SAT/MWF/IFT volumes; T5-T12 vertebral
volume as a height-free body-size reference (muscle-to-vertebra ratio).
Organs: liver, spleen mean HU (5-mm eroded); liver-spleen difference for steatosis.
With a tissue_4_types mask (licensed model) its areas are added as ts4_<level>_*.
"""
import os

import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from skimage.morphology import convex_hull_image
from totalsegmentator.map_to_binary import class_map

HU_MUSCLE = (-29, 150)
HU_FAT = (-190, -30)
HU_LAMA = (-29, 29)
HU_NAMA = (30, 150)
LEVELS = ["T4", "T8", "T10", "T12", "L1", "L3"]
BMD_LEVELS = ["T12", "L1"]
SLAB_MM = 5.0
SMOOTH_MM = 0.75
LUNG_DILATE_MM = 3.0

TASKS = ["total", "abdominal_muscles", "body", "tissue_4_types"]
INV = {t: {v: k for k, v in class_map[t].items()} for t in TASKS}
VERTEBRAE = [f"vertebrae_T{i}" for i in range(1, 13)] + [f"vertebrae_L{i}" for i in range(1, 6)]
CAGE = VERTEBRAE + [f"rib_{s}_{i}" for s in ("left", "right") for i in range(1, 13)] + ["sternum", "costal_cartilages"]
LUNGS = ["lung_upper_lobe_left", "lung_lower_lobe_left", "lung_upper_lobe_right", "lung_middle_lobe_right",
         "lung_lower_lobe_right", "trachea"]
ORGANS = LUNGS + ["heart", "aorta", "esophagus", "liver", "spleen", "stomach", "pancreas", "kidney_left",
                  "kidney_right", "gallbladder", "duodenum", "small_bowel", "colon", "inferior_vena_cava",
                  "superior_vena_cava", "pulmonary_vein", "atrial_appendage_left", "portal_vein_and_splenic_vein",
                  "adrenal_gland_left", "adrenal_gland_right"]
TOTAL_MUSCLES = ["autochthon_left", "autochthon_right", "iliopsoas_left", "iliopsoas_right"]
GI_LUMEN = ["stomach", "colon", "small_bowel", "duodenum", "esophagus"]


def load(path):
    return nib.as_closest_canonical(nib.load(path))


def labels(arr, task, names):
    return np.isin(arr, [INV[task][n] for n in names if n in INV[task]])


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def in_range(img, lo_hi):
    return (img >= lo_hi[0]) & (img <= lo_hi[1])


def hull(m):
    return convex_hull_image(m) if m.sum() >= 3 else np.zeros_like(m)


def vertebra_info(total, name):
    """(centroid z index, fully inside the scan?) of one vertebral BODY; (None, False) if absent.

    Only voxels anterior to the spinal canal are used: the spinous process and laminae slope caudally
    to the level of the next vertebral body, so whole-vertebra extents flag fully imaged bodies as cut
    off and pull the centroid below the mid-body level conventionally used for measurement.
    """
    m = total == INV["total"][f"vertebrae_{name}"]
    if m.sum() < 500:  # stray voxels are not a vertebra
        return None, False
    zs = np.nonzero(m.any(axis=(0, 1)))[0]
    cord = total[:, :, zs.min():zs.max() + 1] == INV["total"]["spinal_cord"]
    body = m.copy()
    if cord.any():
        front = int(np.percentile(np.nonzero(cord)[1], 95))  # anterior edge of the canal (+axis1 = anterior)
    else:
        front = int(np.median(np.nonzero(m)[1]))
    body[:, : front + 1, :] = False
    if body.sum() < 200:
        body = m
    counts = body.sum(axis=(0, 1))
    zb = np.nonzero(counts)[0]
    complete = zb.min() > 0 and zb.max() < total.shape[2] - 1
    return float((np.arange(len(counts)) * counts).sum() / counts.sum()), bool(complete)


class Case:
    def __init__(self, ct_path, seg_dir, harmonize=True):
        img = load(ct_path)
        ct = np.asanyarray(img.dataobj).astype(np.float32)
        if ct.ndim == 4:  # dcm2niix stacked repeated acquisitions; TotalSegmentator segments the first volume
            ct = ct[..., 0]
        self.zooms = tuple(float(z) for z in img.header.get_zooms()[:3])
        self.raw = ct
        k = max(1, int(round(SLAB_MM / self.zooms[2])))
        if harmonize:
            h = ndi.uniform_filter1d(ct, k, axis=2, mode="nearest") if k > 1 else ct
            self.ct = ndi.gaussian_filter(h, sigma=(SMOOTH_MM / self.zooms[0], SMOOTH_MM / self.zooms[1], 0))
        else:  # native slice, no smoothing (only for the harmonisation ablation)
            self.ct, k = ct, 1
        self.slab_k = k
        self.px_cm2 = self.zooms[0] * self.zooms[1] / 100.0
        self.vox_cm3 = self.px_cm2 * self.zooms[2] / 10.0

        self.segs = {}
        for t in TASKS:
            p = os.path.join(seg_dir, f"{t}.nii.gz")
            if os.path.exists(p):
                s = load(p)
                assert s.shape == ct.shape, f"{t} shape {s.shape} != CT {ct.shape}"
                self.segs[t] = np.asanyarray(s.dataobj).astype(np.int16)
        tot, am = self.segs["total"], self.segs["abdominal_muscles"]
        self.muscle = (am > 0) | labels(tot, "total", TOTAL_MUSCLES)
        self.esm = labels(am, "abdominal_muscles", ["erector_spinae_left", "erector_spinae_right",
                                                    "transversospinalis_left", "transversospinalis_right"])
        self.pm = labels(am, "abdominal_muscles", ["pectoralis_major_left", "pectoralis_major_right"])
        self.cage = labels(tot, "total", CAGE)
        self.organs_m = labels(tot, "total", ORGANS)
        self.lungs = labels(tot, "total", LUNGS)
        self.gi = labels(tot, "total", GI_LUMEN)  # food/faeces can fall in the fat HU window
        self.cord = labels(tot, "total", ["spinal_cord"])
        self.body = self.segs["body"] == INV["body"]["body_trunc"]
        self.lung_r = max(1, int(round(LUNG_DILATE_MM / self.zooms[0])))

    # ---------- tissue compartments on one axial slice ----------
    def slice_masks(self, z):
        img = self.ct[:, :, z]
        mus = self.muscle[:, :, z]
        body = ndi.binary_fill_holes(self.body[:, :, z])
        cage = self.cage[:, :, z]
        outer = hull(mus | cage) & body
        inner = hull(cage | self.organs_m[:, :, z]) & outer
        lung = ndi.binary_dilation(self.lungs[:, :, z], structure=disk(self.lung_r))
        fat = in_range(img, HU_FAT) & body & ~mus
        hu_mus = in_range(img, HU_MUSCLE)
        return dict(
            muscle=mus & hu_mus,
            lama=mus & in_range(img, HU_LAMA),
            nama=mus & in_range(img, HU_NAMA),
            imat=mus & in_range(img, HU_FAT),
            sat=fat & ~outer,
            mwf=fat & outer & ~inner,
            ift=fat & inner & ~lung & ~self.gi[:, :, z],
            esm=self.esm[:, :, z] & hu_mus,
            pm=self.pm[:, :, z] & hu_mus,
            body=body,
        )

    def level_metrics(self, name):
        zc, complete = vertebra_info(self.segs["total"], name)
        out = {f"{name}_present": int(zc is not None), f"{name}_complete": int(complete)}
        if zc is None or not complete:
            return out
        z = int(round(zc))
        m = self.slice_masks(z)
        img = self.ct[:, :, z]
        for k in ["muscle", "lama", "nama", "imat", "sat", "mwf", "ift", "esm", "pm"]:
            key = {"muscle": "sma", "esm": "esm_area", "pm": "pm_area"}.get(k, k)
            out[f"{name}_{key}"] = float(m[k].sum() * self.px_cm2)
        for k, key in [("muscle", "smd"), ("esm", "esm_hu"), ("pm", "pm_hu")]:
            out[f"{name}_{key}"] = float(img[m[k]].mean()) if m[k].sum() > 50 else np.nan
        # same voxels, native (unsmoothed) HU: no slab partial volume; noise averages out over the mask
        raw = self.raw[:, :, z]
        out[f"{name}_smd_raw"] = float(raw[m["muscle"]].mean()) if m["muscle"].sum() > 50 else np.nan
        b = m["body"]
        out[f"{name}_fov_trunc"] = int(b[0, :].any() or b[-1, :].any() or b[:, 0].any() or b[:, -1].any())
        out[f"{name}_z"] = z
        if "tissue_4_types" in self.segs:
            t4 = self.segs["tissue_4_types"][:, :, z]
            for lab, key in [("subcutaneous_fat", "sat"), ("torso_fat", "ift"),
                             ("skeletal_muscle", "sma"), ("intermuscular_fat", "imat")]:
                msk = t4 == INV["tissue_4_types"][lab]
                if key == "sma":
                    msk &= in_range(img, HU_MUSCLE)
                    out[f"ts4_{name}_smd"] = float(img[msk].mean()) if msk.sum() > 50 else np.nan
                out[f"ts4_{name}_{key}"] = float(msk.sum() * self.px_cm2)
        if name in BMD_LEVELS:
            out[f"{name}_trab_hu"] = self.trabecular_hu(name, z)
        return out

    def trabecular_hu(self, name, zc):
        """Mean raw HU of the eroded vertebral body anterior to the spinal canal (3 slices around centroid)."""
        lab = INV["total"][f"vertebrae_{name}"]
        er = max(1, int(round(3.0 / self.zooms[0])))
        vals = []
        for z in range(max(0, zc - 1), min(self.raw.shape[2], zc + 2)):
            v = self.segs["total"][:, :, z] == lab
            if v.sum() < 50:
                continue
            c = self.cord[:, :, z]
            front = np.nonzero(c.any(axis=0))[0].max() if c.any() else int(np.nonzero(v.any(axis=0))[0].mean())
            body = v.copy()
            body[:, : front + er + 1] = False  # keep voxels anterior of the canal (+axis1 = anterior)
            lab_cc, n = ndi.label(body)
            if n == 0:
                continue
            body = lab_cc == (np.argmax(np.bincount(lab_cc.ravel())[1:]) + 1)
            core = ndi.binary_erosion(body, structure=disk(er))
            if core.sum() >= 20:
                vals.append(self.raw[:, :, z][core])
        return float(np.concatenate(vals).mean()) if vals else np.nan

    def volumetric(self, top="T5", bottom="T12"):
        zt, ok_t = vertebra_info(self.segs["total"], top)
        zb, ok_b = vertebra_info(self.segs["total"], bottom)
        if zt is None or zb is None or not (ok_t and ok_b):
            return {"vol_available": 0}
        n = dict(muscle=0, sat=0, mwf=0, ift=0)
        hu_sum = hu_raw_sum = 0.0
        for z in range(int(round(zb)), int(round(zt)) + 1):
            m = self.slice_masks(z)
            for k in n:
                n[k] += int(m[k].sum())
            hu_sum += float(self.ct[:, :, z][m["muscle"]].sum())
            hu_raw_sum += float(self.raw[:, :, z][m["muscle"]].sum())
        vert = labels(self.segs["total"], "total", [f"vertebrae_T{i}" for i in range(5, 13)])
        out = {"vol_available": 1, "vol_span_cm": (zt - zb) * self.zooms[2] / 10.0,
               "vol_muscle_hu": hu_sum / max(n["muscle"], 1),
               "vol_muscle_hu_raw": hu_raw_sum / max(n["muscle"], 1),
               "vol_vert_T5_T12_cm3": float(vert.sum() * self.vox_cm3)}  # body-size proxy when height is missing
        out.update({f"vol_{k}_cm3": v * self.vox_cm3 for k, v in n.items()})
        return out

    def organs(self):
        out = {}
        for org in ["liver", "spleen"]:
            m = self.segs["total"] == INV["total"][org]
            core = ndi.distance_transform_edt(m, sampling=self.zooms) > 5.0 if m.sum() > 1000 else m & False
            out[f"{org}_hu"] = float(self.raw[core].mean()) if core.sum() > 100 else np.nan
            if org == "liver":  # image-noise proxy before/after harmonisation (homogeneous parenchyma)
                out["liver_sd_raw"] = float(self.raw[core].std()) if core.sum() > 100 else np.nan
                out["liver_sd_harmonised"] = float(self.ct[core].std()) if core.sum() > 100 else np.nan
        out["liver_minus_spleen_hu"] = out["liver_hu"] - out["spleen_hu"]
        return out

    def run(self):
        row = {"px_mm": self.zooms[0], "slice_mm": self.zooms[2], "nz": int(self.raw.shape[2]),
               "slab_k": self.slab_k}
        for lv in LEVELS:
            row.update(self.level_metrics(lv))
        row.update(self.volumetric())
        row.update(self.organs())
        return row
