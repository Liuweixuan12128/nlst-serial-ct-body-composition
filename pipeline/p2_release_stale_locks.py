"""Release segmentation locks left by killed workers: lock older than --minutes and the case's masks incomplete.
Each released case id is appended to logs/p2_seg_released_locks.txt.

Usage: python p2_release_stale_locks.py [--minutes 10]
"""
import argparse, glob, os, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEG = os.path.join(ROOT, "data", "seg", "nlst_serial")
TASKS = ["total", "abdominal_muscles", "body"]

ap = argparse.ArgumentParser()
ap.add_argument("--minutes", type=float, default=10)
a = ap.parse_args()
released = []
for lk in glob.glob(os.path.join(SEG, "_locks", "*.lock")):
    cid = os.path.basename(lk)[:-5]
    complete = all(os.path.exists(os.path.join(SEG, cid, f"{t}.nii.gz")) for t in TASKS)
    if not complete and time.time() - os.path.getmtime(lk) > a.minutes * 60:
        os.remove(lk)
        released.append(cid)
with open(os.path.join(ROOT, "logs", "p2_seg_released_locks.txt"), "a") as f:
    f.write("".join(c + "\n" for c in released))
print(f"released {len(released)} stale locks")
