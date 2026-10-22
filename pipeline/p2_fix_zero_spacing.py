"""Paper 2: repair NIfTIs whose slice axis has zero length (duplicated slice positions in the DICOM series).

Some NLST series end with two images at the same ImagePositionPatient; dcm2niix then writes a zero slice spacing and
TotalSegmentator fails ("SVD did not converge"). For each such NIfTI: download the series again, keep one image per
slice position (lowest InstanceNumber), convert, check the slice spacing, and remove the old segmentation and lock so
that a segmentation worker processes the case again. Every repaired case is logged.

Usage: python p2_fix_zero_spacing.py --out data/nlst_serial
"""
import argparse, glob, os, shutil, subprocess, sys

import nibabel as nib
import numpy as np
import pandas as pd
import pydicom

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def zero_spacing(f):
    try:
        return np.linalg.norm(nib.load(f).affine[:3, 2]) < 0.3
    except Exception:
        return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from idc_index import IDCClient
    c = IDCClient()
    exe = os.path.join(os.path.dirname(sys.executable), "dcm2niix")
    nii = os.path.join(a.out, "nifti")
    seg = os.path.join(ROOT, "data", "seg", "nlst_serial")
    man = pd.concat([pd.read_csv(f, dtype={"pid": str}) for f in glob.glob(os.path.join(a.out, "p2*_series.csv"))])
    uid = dict(zip(man.case_id, man.SeriesInstanceUID))
    bad = [os.path.basename(f)[:-7] for f in sorted(glob.glob(os.path.join(nii, "*.nii.gz"))) if zero_spacing(f)]
    print(f"zero-spacing NIfTIs: {len(bad)}", flush=True)
    log = []
    for cid in bad:
        if os.path.islink(os.path.join(nii, f"{cid}.nii.gz")):  # paper-1 file: unlink here, never modify the source
            os.remove(os.path.join(nii, f"{cid}.nii.gz"))
            open(os.path.join(nii, f"{cid}.nii.gz"), "w").close()
            if os.path.islink(os.path.join(seg, cid)):
                os.remove(os.path.join(seg, cid))
        tmp = os.path.join(a.out, "dicom_fix", cid)
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        c.download_from_selection(seriesInstanceUID=[uid[cid]], downloadDir=tmp, dirTemplate="%SeriesInstanceUID",
                                  quiet=True)
        files = glob.glob(os.path.join(tmp, "*", "*.dcm"))
        meta = []
        for f in files:
            d = pydicom.dcmread(f, stop_before_pixels=True)
            meta.append((tuple(round(float(x), 3) for x in d.ImagePositionPatient), int(d.InstanceNumber), f))
        meta.sort()
        seen, dropped = set(), 0
        for pos, _, f in meta:
            if pos in seen:
                os.remove(f)
                dropped += 1
            else:
                seen.add(pos)
        out = os.path.join(nii, f"{cid}.nii.gz")
        os.remove(out)
        for j in glob.glob(os.path.join(nii, f"{cid}*.json")):
            os.remove(j)
        src = glob.glob(os.path.join(tmp, "*"))[0]
        subprocess.run([exe, "-z", "y", "-f", cid, "-o", nii, src], capture_output=True)
        made = sorted(glob.glob(os.path.join(nii, f"{cid}*.nii.gz")), key=os.path.getsize)
        for e in made[:-1]:
            os.remove(e)
        if made and made[-1] != out:
            os.replace(made[-1], out)
        ok = os.path.exists(out) and not zero_spacing(out)
        shutil.rmtree(os.path.join(seg, cid), ignore_errors=True)
        lock = os.path.join(seg, "_locks", f"{cid}.lock")
        if os.path.exists(lock):
            os.remove(lock)
        shutil.rmtree(tmp, ignore_errors=True)
        dz = float(np.linalg.norm(nib.load(out).affine[:3, 2])) if os.path.exists(out) else np.nan
        log.append(dict(case_id=cid, status="repaired" if ok else "FAILED", n_images=len(files),
                        duplicates_dropped=dropped, slice_spacing_mm=round(dz, 3)))
        print(log[-1], flush=True)
    shutil.rmtree(os.path.join(a.out, "dicom_fix"), ignore_errors=True)
    path = os.path.join(a.out, "p2_zero_spacing_repairs.csv")
    old = pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()
    pd.concat([old, pd.DataFrame(log)]).drop_duplicates("case_id", keep="last").to_csv(path, index=False)


if __name__ == "__main__":
    main()
