"""Paper 2 QC: tile per-participant montages (results/p2/qc_montages) into review sheets of 4, each captioned with
its flags from p2_qc_flags.csv, so that every flagged participant is looked at.

Usage: python p2_qc_sheets.py --qc results/p2
"""
import argparse, glob, os, textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--qc", required=True)
a = ap.parse_args()
fl = pd.read_csv(os.path.join(a.qc, "p2_qc_flags.csv"), dtype={"pid": str}).fillna("")
pngs = sorted(glob.glob(os.path.join(a.qc, "qc_montages", "*.png")))
out = os.path.join(a.qc, "qc_sheets")
os.makedirs(out, exist_ok=True)
for i in range(0, len(pngs), 4):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax in axes.ravel():
        ax.axis("off")
    for ax, f in zip(axes.ravel(), pngs[i:i + 4]):
        pid = os.path.basename(f)[:-4]
        ax.imshow(plt.imread(f))
        r = fl[fl.pid == pid]
        txt = "; ".join(f"T{int(x['round'])}: {x.scan_flags} {x.change_flags}".strip()
                        for _, x in r.iterrows() if (x.scan_flags or x.change_flags))
        grp = r.group.iloc[0] if len(r) else ""
        ax.set_title(f"{pid} ({grp})\n" + "\n".join(textwrap.wrap(txt, 90)[:3]), fontsize=8, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(out, f"sheet_{i // 4 + 1:03d}.png"), dpi=80)
    plt.close(fig)
print(f"{len(pngs)} montages -> {-(-len(pngs) // 4)} sheets in {out}")
