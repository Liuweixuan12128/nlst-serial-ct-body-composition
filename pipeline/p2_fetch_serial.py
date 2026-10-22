"""Paper 2: download the series listed in p2_series.csv chunk by chunk, convert to NIfTI, delete the DICOM.

Series already processed for paper 1 (same SeriesInstanceUID in data/nlst/nlst_selected_series_last.csv) are linked,
NIfTI and segmentation both, instead of downloaded again. DICOM is removed after each chunk converts, so the peak
footprint is one chunk; every series can be fetched again from IDC by its SeriesInstanceUID.

Usage: python p2_fetch_serial.py --manifest data/nlst_serial/p2_series.csv --out data/nlst_serial [--chunk 150]
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def convert(exe, src, nii, name):
    subprocess.run([exe, "-z", "y", "-f", name, "-o", nii, src], capture_output=True)
    out = sorted(glob.glob(os.path.join(nii, f"{name}*.nii.gz")), key=os.path.getsize)
    for e in out[:-1]:
        os.remove(e)  # keep the largest volume if dcm2niix split the series
    for e in glob.glob(os.path.join(nii, f"{name}?*.json")):
        os.remove(e)
    if out and out[-1] != os.path.join(nii, f"{name}.nii.gz"):
        os.replace(out[-1], os.path.join(nii, f"{name}.nii.gz"))
    return os.path.exists(os.path.join(nii, f"{name}.nii.gz"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk", type=int, default=150)
    ap.add_argument("--tmp", default="dicom_tmp", help="scratch folder name under --out (one per concurrent run)")
    a = ap.parse_args()
    from idc_index import IDCClient
    c = IDCClient()
    exe = os.path.join(os.path.dirname(sys.executable), "dcm2niix")
    man = pd.read_csv(a.manifest, dtype={"pid": str})
    nii = os.path.join(a.out, "nifti")
    seg = os.path.join(ROOT, "data", "seg", "nlst_serial")
    os.makedirs(nii, exist_ok=True)
    os.makedirs(seg, exist_ok=True)

    old = pd.read_csv(os.path.join(ROOT, "data", "nlst", "nlst_selected_series_last.csv"), dtype={"pid": str})
    old["old_id"] = old.pid + "_T" + old["round"].astype(str)
    reuse = man.merge(old[["SeriesInstanceUID", "old_id"]], on="SeriesInstanceUID")
    for _, r in reuse.iterrows():
        src_nii = os.path.join(ROOT, "data", "nlst", "nifti", f"{r.old_id}.nii.gz")
        src_seg = os.path.join(ROOT, "data", "seg", "nlst", r.old_id)
        for src, dst in ((src_nii, os.path.join(nii, f"{r.case_id}.nii.gz")), (src_seg, os.path.join(seg, r.case_id))):
            if os.path.exists(src) and not os.path.lexists(dst):
                os.symlink(os.path.relpath(src, os.path.dirname(dst)), dst)
    todo = man[~man.SeriesInstanceUID.isin(set(reuse.SeriesInstanceUID))]
    todo = todo[[not os.path.exists(os.path.join(nii, f"{cid}.nii.gz")) for cid in todo.case_id]]
    print(f"linked from paper 1: {len(reuse)}; to download: {len(todo)}", flush=True)

    tmp = os.path.join(a.out, a.tmp)
    failed = []
    for i in range(0, len(todo), a.chunk):
        part = todo.iloc[i:i + a.chunk]
        shutil.rmtree(tmp, ignore_errors=True)
        os.makedirs(tmp)
        c.download_from_selection(seriesInstanceUID=part.SeriesInstanceUID.tolist(), downloadDir=tmp,
                                  dirTemplate="%SeriesInstanceUID", quiet=True)
        for _, r in part.iterrows():
            if not convert(exe, os.path.join(tmp, r.SeriesInstanceUID), nii, r.case_id):
                failed.append(r.case_id)
        shutil.rmtree(tmp, ignore_errors=True)
        print(f"chunk {i // a.chunk + 1}/{-(-len(todo) // a.chunk)} done; failed so far {len(failed)}", flush=True)
    stem = os.path.basename(a.manifest).replace("_series.csv", "")
    pd.Series(failed, name="case_id").to_csv(os.path.join(a.out, f"{stem}_fetch_failed.csv"), index=False)
    print("FETCH DONE", flush=True)


if __name__ == "__main__":
    main()
