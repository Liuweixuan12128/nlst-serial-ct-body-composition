"""Paper 2 Extended Data: how often the field of view cuts the body outline vs the paraspinal compartment, per
vertebral level, on the first screening CT of every control (the justification for the paraspinal measure).

For levels T4, T6, T8, T10, T11, T12 (vertebral-body centroid slice): body_trunc = body mask (TotalSegmentator `body`)
touches any image edge; psm_trunc = paraspinal compartment touches any image edge. Levels not imaged are recorded.

Usage: python p2_fov_truncation.py --out results/p2/p2_fov_truncation.csv [--workers 40]
"""
import argparse, glob, os, sys
from multiprocessing import Pool

import nibabel as nib
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
from bodycomp_metrics import INV, labels, vertebra_info  # noqa: E402

LEVELS = ["T4", "T6", "T8", "T10", "T11", "T12"]
PSM = ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"]


def edge(m):
    return bool(m[0, :].any() or m[-1, :].any() or m[:, 0].any() or m[:, -1].any())


def one(cid):
    seg = os.path.join(ROOT, "data", "seg", "nlst_serial", cid)
    try:
        tot = np.asanyarray(nib.load(os.path.join(seg, "total.nii.gz")).dataobj)
        am = nib.load(os.path.join(seg, "abdominal_muscles.nii.gz"))
        body = nib.load(os.path.join(seg, "body.nii.gz"))
        row = {"case_id": cid}
        for lv in LEVELS:
            zc, ok = vertebra_info(tot, lv)
            if zc is None:
                row[f"{lv}_imaged"] = 0
                continue
            z = int(round(zc))
            b = np.asanyarray(body.dataobj[:, :, z]) == INV["body"]["body_trunc"]
            p = labels(np.asanyarray(am.dataobj[:, :, z]), "abdominal_muscles", PSM)
            row.update({f"{lv}_imaged": 1, f"{lv}_body_trunc": int(edge(b)), f"{lv}_psm_trunc": int(edge(p))})
        return row
    except Exception as e:
        return {"case_id": cid, "error": str(e)[:200]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=40)
    a = ap.parse_args()
    ser = pd.concat([pd.read_csv(f, dtype={"pid": str}) for f in glob.glob(os.path.join(ROOT, "data", "nlst_serial", "p2*_series.csv"))])
    first = ser[ser.group == "control"].sort_values("round").groupby("pid").head(1)
    with Pool(a.workers) as p:
        rows = p.map(one, list(first.case_id), chunksize=4)
    d = pd.DataFrame(rows)
    d.to_csv(a.out, index=False)
    s = {lv: (d[f"{lv}_body_trunc"].mean(), d[f"{lv}_psm_trunc"].mean(), d[f"{lv}_imaged"].mean()) for lv in LEVELS}
    print(len(d), {k: tuple(round(x, 3) for x in v) for k, v in s.items()}, file=sys.stderr)


if __name__ == "__main__":
    main()
