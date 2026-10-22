"""Paper 2 Extended Data figures (image-based and descriptive):
ED1 field-of-view truncation by level; ED2 quality control; ED3 segmentation across scanners and kernels;
ED6 matching balance and case characteristics.  (ED4/5/7/8 are in p2_figs_main.py with the analysis figures.)
Usage: python figures/p2_figs_ed.py [ed1 ed2 ed3 ed6]
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from p2_cns_common import (ALL, CASE, INK, INK2, LABEL, MM, MUTED, WASH, axial, dot_ci, imshow_ct, rd, sagittal,
                           scalebar, sel, style)
from p2_premium import BAND
from style import W2, panel_label


def ed1():
    """FOV truncation of the body outline vs the paraspinal compartment, controls' first scan, by vertebral level."""
    style.setup()
    d = rd("p2_fov_truncation.csv")
    d = d[d.get("error").isna()] if "error" in d else d
    lv = ["T4", "T6", "T8", "T10", "T11", "T12"]
    rows = []
    for L in lv:
        im = d[d[f"{L}_imaged"] == 1]
        rows.append(dict(level=L, n=len(im), body=100 * im[f"{L}_body_trunc"].mean(), psm=100 * im[f"{L}_psm_trunc"].mean(),
                         imaged=100 * d[f"{L}_imaged"].mean()))
    s = pd.DataFrame(rows)
    fig = plt.figure(figsize=(W2, 70 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.4, 1.0], wspace=0.25, left=0.07, right=0.97, top=0.88, bottom=0.17)
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(lv))
    ax.bar(x - 0.19, s.body, 0.36, color=WASH, edgecolor=INK2, lw=0.5, label="Body outline cut by the field of view")
    ax.bar(x + 0.19, s.psm, 0.36, color=style.TISSUE["muscle"], lw=0, label="Paraspinal compartment cut")
    for xi, p in enumerate(s.psm):   # paraspinal bars are too short to read off the axis
        ax.text(xi + 0.19, p + 1.5, f"{p:.1f}", ha="center", fontsize=5.3, color=INK2)
    ax.axvspan(0.5, 4.5, color=BAND, lw=0, zorder=0)
    ax.text(2.5, 104, "T6–T11", ha="center", fontsize=6, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(lv)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Scans truncated (%)")
    ax.set_xlabel("Vertebral level (first screening CT of each control)")
    ax.legend(loc="upper left", bbox_to_anchor=(0, 1.16), ncol=2)
    panel_label(ax, "a", x=-0.09, y=1.1)
    style.source_data(s, "P2_EDFig1", "a")
    # b, c: one truncated example at T8 (body cut, paraspinal intact) and its full-FOV slice
    ex = d[(d.T8_body_trunc == 1) & (d.T8_psm_trunc == 0)].case_id.sort_values().iloc[0]
    for k, (cid, lab, title) in enumerate([(ex, "b", "Body cut, paraspinal intact (T8)")]):
        import nibabel as nib
        from p2_cns_common import NII, SEG, style_inv
        from bodycomp_metrics import labels, vertebra_info
        img = nib.load(f"{NII}/{cid}.nii.gz")
        ct = np.asanyarray(img.dataobj)
        tot = np.asanyarray(nib.load(f"{SEG}/{cid}/total.nii.gz").dataobj)
        am = np.asanyarray(nib.load(f"{SEG}/{cid}/abdominal_muscles.nii.gz").dataobj)
        body = np.asanyarray(nib.load(f"{SEG}/{cid}/body.nii.gz").dataobj)
        z = int(round(vertebra_info(tot, "T8")[0]))
        ax = fig.add_subplot(gs[0, 1 + k])
        imshow_ct(ax, np.rot90(ct[:, :, z]), np.rot90(labels(am[:, :, z], "abdominal_muscles",
                  ["erector_spinae_left", "erector_spinae_right", "transversospinalis_left", "transversospinalis_right"])))
        ax.contour(np.rot90(body[:, :, z] == style_inv("body", "body_trunc")), levels=[0.5], colors="#eda100",
                   linewidths=0.6)
        ax.set_title(title, fontsize=6.5, pad=3)
        scalebar(ax, float(img.header.get_zooms()[0]), 50)
        panel_label(ax, lab, x=-0.05, y=1.08)
    return style.save(fig, "P2_EDFig1")


def ed2():
    """Quality control: exclusions by reason, and one example of each failure type."""
    style.setup()
    q = rd("p2_qc_review.csv")
    fl = rd("p2_qc_flags.csv").fillna("")
    reason = q.reason.str.replace(r" by [0-9.]+%.*", "", regex=True)
    cat = np.select([reason.str.contains("slab length"), reason.str.contains("prone"), reason.str.contains("field of view"),
                     reason.str.contains("scoliosis"), reason.str.contains("streak"), reason.str.contains("segmentation failure")],
                    ["Vertebral level misidentified", "Prone acquisition", "Posterior field-of-view cut",
                     "Severe spinal deformity", "Streak artefact", "Segmentation failure"], "Other")
    auto_inc = (fl.scan_flags.str.contains("slab incomplete")).sum()
    auto_edge = (fl.scan_flags.str.contains("image edge") & ~fl.scan_flags.str.contains("slab incomplete")).sum()
    counts = pd.Series(cat).value_counts()
    counts["T6–T11 slab not fully imaged"] = auto_inc
    counts["Paraspinal at image edge"] = auto_edge
    counts = counts.sort_values()
    FH = 100
    fig = plt.figure(figsize=(W2, FH * MM))
    ax = fig.add_axes([0.24, 0.47, 0.34, 0.47])
    y = np.arange(len(counts))
    ax.barh(y, counts.values, color=[style.TISSUE["muscle"] if "image" in k or "slab not" in k else INK2
                                     for k in counts.index], height=0.62)
    for yi, v in zip(y, counts.values):
        ax.text(v + 2, yi, f"{v}", va="center", fontsize=5.5)
    ax.set_yticks(y)
    ax.set_yticklabels(counts.index, fontsize=6)
    ax.set_xlabel("Scans excluded")
    ax.set_xlim(0, counts.max() * 1.22)
    panel_label(fig, "a", x=0.01, y=0.985, fig=fig)
    style.source_data(counts.rename("n").reset_index().rename(columns={"index": "reason"}), "P2_EDFig2", "a")
    ex = [("129379_T0", "T8", "Severe spinal deformity"), ("116433_T0", "T8", "Posterior field-of-view cut"),
          ("105514_T0", "T8", "Streak artefact"), ("132231_T2", "T8", "Segmentation failure"),
          ("212584_T0", "T8", "Prone acquisition"), ("101347_T0", "T8", "Noisy but retained")]
    panel_label(fig, "b", x=0.01, y=0.36, fig=fig)
    w = 0.158
    for k, (cid, lv, title) in enumerate(ex):
        img, m, px = axial(cid, lv, half_mm=95, shift_mm=0)
        hfrac = w * W2 / MM * 0.72 / FH          # nominal crop aspect; images keep true aspect inside the box
        ax = fig.add_axes([0.015 + k * (w + 0.006), 0.02, w, hfrac])
        imshow_ct(ax, img, m)
        ax.set_facecolor("black")
        ax.set_title(title, fontsize=6, pad=2.5)
        scalebar(ax, px)
    # c: level misidentification, sagittal of two rounds of the same participant
    panel_label(fig, "c", x=0.66, y=0.985, fig=fig)
    for k, cid in enumerate(["111243_T0", "111243_T1"]):
        ax = fig.add_axes([0.68 + k * 0.16, 0.44, 0.15, 0.5])
        sag, vert, lvz, zooms = sagittal(cid)
        img = np.rot90(sag)
        ax.imshow(img, cmap="gray", vmin=-160, vmax=240, aspect=zooms[2] / zooms[1], interpolation="lanczos")
        H = img.shape[0]
        for name in ["T6", "T11"]:
            if lvz[name] is not None:
                yy = H - 1 - lvz[name]
                ax.axhline(yy, color=style.TISSUE["muscle"], lw=0.7)
                ax.text(img.shape[1] * 0.97, yy - 4, name, color="white", fontsize=6, ha="right", va="bottom",
                        fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        sl = fl[fl.case_id == cid].slab_mm.iloc[0]
        ax.set_title(f"{cid.split('_')[1]}: slab {sl:.0f} mm", fontsize=6, pad=2)

    return style.save(fig, "P2_EDFig2")


def ed3():
    """Segmentation of the paraspinal compartment across scanner vendors and reconstruction kernels (controls)."""
    style.setup()
    s = sel("series")
    m = rd("p2_paraspinal.csv").set_index("case_id")
    combos = [("GE", "STANDARD", "GE · STANDARD"), ("SE", "B30F", "Siemens · B30f"), ("PH", "C", "Philips · C"),
              ("TO", "FC10", "Toshiba · FC10"), ("SE", "B50F", "Siemens · B50f (sharp)"), ("GE", "BONE", "GE · BONE (sharp)")]
    rng = np.random.default_rng(20261009)
    fig, axes = plt.subplots(2, 6, figsize=(W2, 62 * MM), gridspec_kw=dict(hspace=0.04, wspace=0.04,
                             left=0.01, right=0.99, top=0.9, bottom=0.08))
    src = []
    for j, (v, k, title) in enumerate(combos):
        pool = s[(s.group == "control") & (s.vendor == v) & (s.kernel == k) & (s["round"] == 0)].case_id.tolist()
        pool = [c for c in pool if c in m.index and m.loc[c, "psm_edge"] == 0 and m.loc[c, "psm_available"] == 1]
        cid = pool[rng.integers(len(pool))]
        for i, lv in enumerate(["T8", "T10"]):
            ax = axes[i, j]
            img, msk, px = axial(cid, lv, half_mm=75, shift_mm=0)
            imshow_ct(ax, img, msk)
            if i == 0:
                ax.set_title(title, fontsize=6, pad=3)
            ax.text(0.03, 0.97, lv, transform=ax.transAxes, color="white", fontsize=5.5, va="top")
            if i == 1:
                ax.text(0.5, -0.08, f"fat {m.loc[cid, 'psm_fat_pct']:.1f}% · {m.loc[cid, 'psm_hu']:.0f} HU",
                        transform=ax.transAxes, ha="center", va="top", fontsize=5.5, color=INK2)
        src.append(dict(case_id=cid, vendor=v, kernel=k, fat_pct=m.loc[cid, "psm_fat_pct"], hu=m.loc[cid, "psm_hu"]))
    style.source_data(pd.DataFrame(src), "P2_EDFig3", "a")
    return style.save(fig, "P2_EDFig3")


def ed6():
    """Matching balance and case characteristics."""
    style.setup()
    cas, ctl = sel("cases"), sel("controls")
    maa = rd("p2_maa.csv").astype({"pid": str})
    ids = set(maa.pid)
    cas, ctl = cas[cas.pid.isin(ids)], ctl[ctl.pid.isin(ids)]
    ser = sel("series")
    v0 = ser.sort_values("round").groupby("pid").vendor.first()
    def smd(a, b):
        a, b = np.asarray(a, float), np.asarray(b, float)
        return (a.mean() - b.mean()) / np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    rows = [("Age at entry", smd(cas.age, ctl.age)), ("Men", smd(cas.gender == 1, ctl.gender == 1)),
            ("Current smoker", smd(cas.cigsmok == 1, ctl.cigsmok == 1))]
    for v, lab in [("GE", "Vendor GE"), ("SE", "Vendor Siemens"), ("PH", "Vendor Philips"), ("TO", "Vendor Toshiba")]:
        rows.append((lab, smd(cas.pid.map(v0) == v, ctl.pid.map(v0) == v)))
    rows.append(("White race (not matched)", smd(cas.race == 1, ctl.race == 1)))
    b = pd.DataFrame(rows, columns=["variable", "smd"])
    fig = plt.figure(figsize=(W2, 68 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.5, left=0.165, right=0.99, top=0.86, bottom=0.2)
    ax = fig.add_subplot(gs[0, 0])
    y = -np.arange(len(b))
    ax.axvspan(-0.1, 0.1, color=BAND, lw=0)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.plot(b.smd, y, "o", ms=3.5, color=INK)
    ax.set_yticks(y); ax.set_yticklabels(b.variable, fontsize=6)
    ax.set_xlim(-0.25, 0.25)
    ax.set_xlabel("Standardised mean difference\n(cases − controls)")
    panel_label(fig, "a", x=0.01, y=0.97, fig=fig)
    style.source_data(b, "P2_EDFig6", "a")
    st = cas.de_stag_7thed.map({110: "IA", 120: "IB", 210: "IIA", 220: "IIB", 310: "IIIA", 320: "IIIB", 400: "IV"}).fillna("Unk")
    ax = fig.add_subplot(gs[0, 1])
    order = ["IA", "IB", "IIA", "IIB", "IIIA", "IIIB", "IV", "Unk"]
    vc = st.value_counts().reindex(order, fill_value=0)
    p2ids = set(pd.read_csv(f"{style.ROOT}/data/nlst_serial/p2_cases.csv", dtype={"pid": str}).pid)
    rs = cas.assign(st=st, r=cas.pid.isin(p2ids)).groupby("st").r.sum().reindex(order, fill_value=0)
    ax.bar(order, rs.values, color=CASE, lw=0, label="Resected")
    ax.bar(order, (vc - rs).values, bottom=rs.values, color=WASH, edgecolor=CASE, lw=0.5, label="Not resected")
    ax.set_ylabel("Cases")
    ax.set_xlabel("Stage (7th edition)")
    ax.tick_params(axis="x", labelsize=5.5)
    ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.16))
    panel_label(ax, "b", x=-0.3, y=1.08)
    style.source_data(pd.DataFrame({"stage": order, "n": vc.values, "resected": rs.values}), "P2_EDFig6", "b")
    ax = fig.add_subplot(gs[0, 2])
    last = cas.apply(lambda r: r.candx_days - max(r[f"scr_days{k}"] for k in (0, 1, 2)
                                                   if pd.notna(r[f"scr_days{k}"]) and r[f"scr_days{k}"] <= r.candx_days), axis=1) / 365.25
    ax.hist(last, bins=np.arange(0, 5.75, 0.25), color=CASE, alpha=0.85, lw=0)
    ax.axvline(np.median(last), color=INK, lw=0.7, ls=(0, (3, 2)))
    ax.set_xlabel("Last screening CT to diagnosis (years)")
    ax.set_ylabel("Cases")
    panel_label(ax, "c", x=-0.3, y=1.08)
    style.source_data(pd.DataFrame({"years": last}), "P2_EDFig6", "c")
    return style.save(fig, "P2_EDFig6")


if __name__ == "__main__":
    todo = sys.argv[1:] or ["ed1", "ed2", "ed3", "ed6"]
    for t in todo:
        print(globals()[t]())
