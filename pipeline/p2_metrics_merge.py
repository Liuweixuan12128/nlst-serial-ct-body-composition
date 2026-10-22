"""Merge sharded p2_metrics outputs into results/p2/p2_paraspinal.csv (latest error-free row per scan wins)."""
import glob, os, sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
files = sorted(glob.glob(os.path.join(ROOT, "results", "p2", "shards", "pass*_*.csv")))
d = pd.concat([pd.read_csv(f).assign(_src=os.path.basename(f)) for f in files], ignore_index=True)
if "error" not in d:
    d["error"] = pd.NA
d["_ok"] = d["error"].isna()
d = d.sort_values(["case_id", "_ok", "_src"]).drop_duplicates("case_id", keep="last").drop(columns=["_ok", "_src"])
out = os.path.join(ROOT, "results", "p2", "p2_paraspinal.csv")
d.to_csv(out, index=False)
print(f"{len(files)} shard files -> {len(d)} scans ({d['error'].notna().sum()} with error) -> {out}", file=sys.stderr)
