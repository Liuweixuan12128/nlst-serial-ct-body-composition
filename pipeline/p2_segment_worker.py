"""Paper 2: segmentation worker that keeps polling a NIfTI folder while it is still being filled.

Same TotalSegmentator tasks and settings as segment.py. Each case is claimed with an exclusive lock file, so any
number of workers on any number of nodes can share one folder without processing a case twice. The worker exits when
the fetch log reports FETCH DONE and no unclaimed case is left.

Usage: python p2_segment_worker.py --nii-dir data/nlst_serial/nifti --seg-dir data/seg/nlst_serial
                                   --done-flag logs/p2_fetch.log
"""
import argparse, glob, os, socket, sys, time, traceback

os.environ.setdefault("TOTALSEG_HOME_DIR", os.path.join(os.path.dirname(sys.executable), "..", "totalseg_weights"))
from totalsegmentator.python_api import totalsegmentator  # noqa: E402

TASKS = ["total", "abdominal_muscles", "body"]


def complete(d):
    return all(os.path.exists(os.path.join(d, f"{t}.nii.gz")) for t in TASKS)


def claim(lock_dir, cid):
    try:
        fd = os.open(os.path.join(lock_dir, f"{cid}.lock"), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, f"{socket.gethostname()} {os.getpid()} {time.ctime()}\n".encode())
        os.close(fd)
        return True
    except FileExistsError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nii-dir", required=True)
    ap.add_argument("--seg-dir", required=True)
    ap.add_argument("--done-flag", required=True)
    a = ap.parse_args()
    lock_dir = os.path.join(a.seg_dir, "_locks")
    os.makedirs(lock_dir, exist_ok=True)
    while True:
        fetch_done = os.path.exists(a.done_flag) and "FETCH DONE" in open(a.done_flag).read()
        did = 0
        for f in sorted(glob.glob(os.path.join(a.nii_dir, "*.nii.gz"))):
            cid = os.path.basename(f)[:-7]
            out_dir = os.path.join(a.seg_dir, cid)
            if complete(out_dir) or not claim(lock_dir, cid):
                continue
            os.makedirs(out_dir, exist_ok=True)
            for t in TASKS:
                out = os.path.join(out_dir, f"{t}.nii.gz")
                if os.path.exists(out):
                    continue
                t0 = time.time()
                for attempt in range(3):  # concurrent workers can catch TotalSegmentator's config.json mid-write
                    try:
                        totalsegmentator(f, out, task=t, ml=True, device="gpu", quiet=True,
                                         nr_thr_resamp=4, nr_thr_saving=4)
                        print(f"{cid}\t{t}\tok\t{time.time() - t0:.0f}s", flush=True)
                        break
                    except Exception:
                        print(f"{cid}\t{t}\tFAILED (attempt {attempt + 1})\n{traceback.format_exc()}", flush=True)
                        time.sleep(5 + 10 * attempt)
            did += 1
        if fetch_done and did == 0:
            break
        if did == 0:
            time.sleep(60)
    print("WORKER DONE", flush=True)


if __name__ == "__main__":
    main()
