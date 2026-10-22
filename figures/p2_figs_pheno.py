"""Paper 2 Figure 2: the CT reserve score (bone + muscle + lung) and later lung cancer diagnosis.
Reads results/p2/pheno_*.csv from stats/p2_phenotype.R (plan: docs/paper2/P2_03_phenotype_prespec.md) and the 3D
assets of figures/p2_render_reserve3d.py.
Usage: python figures/p2_figs_pheno.py
"""
import colorsys
import os

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from p2_cns_common import CASE, CTRL, INK, INK2, MM, MUTED, NII, ROOT, SEG, axial, crop_alpha, fnum, imshow_ct, rd, \
    scalebar, style
from p2_premium import TIS, forest, leader, raincloud
from style import W2

EX = "100658_T0"   # the example scan of Figure 1
HERO = "p2_reserve3d_az-75_right"


def lab(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=8, fontweight="bold", va="top", ha="left", color=INK)


def t10_crop(cid, half_mm=30):
    """T10 axial slice in bone window with the trabecular core used for measurement (bodycomp_metrics.trabecular_hu)"""
    from bodycomp_metrics import INV, disk, labels, vertebra_info
    img = nib.load(os.path.join(NII, f"{cid}.nii.gz"))
    ct = np.asanyarray(img.dataobj)
    tot = np.asanyarray(nib.load(os.path.join(SEG, cid, "total.nii.gz")).dataobj)
    px = float(img.header.get_zooms()[0])
    z = int(round(vertebra_info(tot, "T10")[0]))
    v = tot[:, :, z] == INV["total"]["vertebrae_T10"]
    c = labels(tot, "total", ["spinal_cord"])[:, :, z]
    er = max(1, int(round(3.0 / px)))
    front = np.nonzero(c.any(axis=0))[0].max() if c.any() else int(np.nonzero(v.any(axis=0))[0].mean())
    body = v.copy()
    body[:, : front + er + 1] = False
    lab_cc, _ = ndi.label(body)
    body = lab_cc == (np.argmax(np.bincount(lab_cc.ravel())[1:]) + 1)
    core = ndi.binary_erosion(body, structure=disk(er))
    xs, ys = np.nonzero(v)
    cx, cy = int(xs.mean()), int(ys.mean())
    h = int(round(half_mm / px))
    sl = (slice(cx - h, cx + h), slice(cy - h, cy + h))
    return np.rot90(ct[:, :, z][sl]), np.rot90(core[sl]), px


def lung_slice(cid, level="T8"):
    """full axial slice through both lungs; mask of lung voxels below -950 HU (their share varies between people,
    unlike the share below the scan's own 15th percentile)"""
    from bodycomp_metrics import LUNGS, labels, vertebra_info
    img = nib.load(os.path.join(NII, f"{cid}.nii.gz"))
    ct = np.asanyarray(img.dataobj)
    tot = np.asanyarray(nib.load(os.path.join(SEG, cid, "total.nii.gz")).dataobj)
    px = float(img.header.get_zooms()[0])
    z = int(round(vertebra_info(tot, level)[0]))
    lungs = labels(tot, "total", LUNGS)
    p15 = float(np.percentile(ct[lungs], 15))
    lz = lungs[:, :, z]
    xs, ys = np.nonzero(lz)
    pad = int(round(8 / px))
    sl = (slice(max(xs.min() - pad, 0), xs.max() + pad), slice(max(ys.min() - pad, 0), ys.max() + pad))
    low = lz & (ct[:, :, z] < -950)
    return np.rot90(ct[:, :, z][sl]), np.rot90(low[sl]), np.rot90(lz[sl]), px, p15


def hue_centroid(img, hue_deg, tol=14, sat_min=0.45):
    """centroid (fractions of width/height from top-left) of saturated pixels of one hue in an RGBA image"""
    a = np.asarray(img).astype(float) / 255
    r, g, b, al = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    mx, mn = np.maximum(np.maximum(r, g), b), np.minimum(np.minimum(r, g), b)
    sat = np.where(mx > 0, (mx - mn) / np.maximum(mx, 1e-6), 0)
    ys, xs = np.nonzero((al > 0.8) & (sat > sat_min))
    hh = np.array([colorsys.rgb_to_hsv(r[y, x], g[y, x], b[y, x])[0] * 360 for y, x in zip(ys, xs)])
    d = np.minimum(np.abs(hh - hue_deg), 360 - np.abs(hh - hue_deg))
    k = d < tol
    return xs[k].mean() / a.shape[1], ys[k].mean() / a.shape[0], int(k.sum())


def pick(r):
    return {"est": float(r["or"]), "lo": float(r["lo"]), "hi": float(r["hi"])}


def fig2():
    style.setup()
    sc, est, qt = rd("pheno_scores.csv"), rd("pheno_estimates.csv"), rd("pheno_quartiles.csv")
    sc = sc[sc.S.notna()]
    sc = sc[sc.set.isin(sc[sc.group == "case"].set)]
    met = rd("p2_paraspinal.csv").set_index("case_id").loc[EX]
    er = lambda name: est[est.analysis == name].iloc[0]
    FH, FX = 186, 40      # existing panels a-f, extra row g-i (mm); 226 mm in all, to fit a page with its title
    root = plt.figure(figsize=(W2, (FH + FX) * MM))
    fig, sf2 = root.subfigures(2, 1, height_ratios=[FH, FX])
    asp = W2 / (FH * MM)

    # ================= a: the three measures in 3D and on CT
    lab(fig, 0.004, 0.992, "a")
    hero = crop_alpha(os.path.join(ROOT, "figures", "assets", f"{HERO}.png"), pad=4)
    hh = 0.3
    hw = hh * hero.width / hero.height / asp
    hx, hy = 0.05, 0.655
    ax = fig.add_axes([hx, hy, hw, hh]); ax.imshow(hero, interpolation="lanczos"); ax.axis("off")
    small = hero.resize((hero.width // 4, hero.height // 4))
    for hue, tis, text, dy in [(38, "bone", "T10 trabecular\nbone", -0.035),
                               (2, "muscle", "Paraspinal\nmuscle, T6–T11", 0.05),
                               (176, "lung", "Lung voxels below\nthe 15th percentile", 0.0)]:
        fx, fy, n = hue_centroid(small, hue, tol=12 if tis != "lung" else 16)
        px_, py_ = hx + fx * hw, hy + (1 - fy) * hh
        if tis == "lung":
            leader(fig, (px_, py_), (0.03, 0.972), text.replace("\n", " "), TIS[tis], ha="left", fs=5.4)
        else:
            leader(fig, (px_, py_), (hx + hw + 0.012, py_ + dy), text, TIS[tis], ha="left", fs=5.4)
    # CT strip: one column of three native slices
    t10, core, px_b = t10_crop(EX)
    psm, msk, px_m = axial(EX, "T8", half_mm=58, shift_mm=0)
    lung, low, lz, px_l, p15 = lung_slice(EX)
    col_x = hx + hw + 0.105
    strip = [(t10, core, px_b, "bone", "Trabecular bone", f"{fnum(met.T10_trab_hu, '{:.0f}')} HU", (-450, 1050), 20),
             (psm, msk, px_m, "muscle", "Paraspinal muscle", f"{met.psm_vol_cm3:.0f} cm³, T6–T11", (-160, 240), 20),
             (lung, None, px_l, "lung", "Lung density", f"p15 {fnum(met.lung_perc15_hu, '{:.0f}')} HU\n{met.lung_laa950_pct:.1f}% < −950 HU", (-1050, -350),
              50)]
    h = 0.088
    yy = hy + hh - h - 0.004
    for im, m, px, tis, title, val, win, bar in strip:
        wa = min(h * im.shape[1] / im.shape[0] / asp, 0.15)
        hu = wa * asp * im.shape[0] / im.shape[1]
        a2 = fig.add_axes([col_x, yy + (h - hu), wa, hu])
        imshow_ct(a2, im, m, color=TIS[tis], lw=0.7, window=win)
        if tis == "lung":
            ov = np.zeros(low.shape + (4,))
            ov[low] = [0.03, 0.47, 0.45, 0.95]
            a2.imshow(ov, interpolation="nearest")
        scalebar(a2, px, bar, color="white" if tis != "lung" else INK)
        a2.text(1.04, 0.5, title, transform=a2.transAxes, fontsize=5.6, color=TIS[tis], fontweight="bold", va="top")
        a2.text(1.04, 0.3, val, transform=a2.transAxes, fontsize=5.6, color=INK, va="top")
        yy -= h + 0.016
    # score definition band
    bx = fig.add_axes([0.035, 0.603, 0.575, 0.042]); bx.axis("off")
    bx.text(0.02, 0.5, r"Reserve score = mean(z$_{\rm bone}$, z$_{\rm muscle}$, z$_{\rm lung}$)", fontsize=6.3,
            va="center", transform=bx.transAxes)

    # ================= b: raincloud of the score
    lab(fig, 0.62, 0.992, "b")
    ax = fig.add_axes([0.675, 0.665, 0.305, 0.285])
    g = {k: sc[sc.group == k].S.values for k in ("control", "case")}
    raincloud(ax, [g["control"], g["case"]], [CTRL, CASE], ["Controls", "Cases"], [1.0, 0.0], bw=0.22, cloud_h=0.55,
              xlim=(-4, 3.2), point_alpha=0.22)
    ax.set_yticks([0.1, 1.1])
    ax.set_yticklabels([f"Cases\n(n = {len(g['case'])})", f"Controls\n(n = {len(g['control']):,})"], fontsize=5.8)
    for t, c in zip(ax.get_yticklabels(), [CASE, CTRL]):
        t.set_color(c)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-4, 3.2); ax.set_ylim(-0.42, 1.9)
    ax.axvline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=0)
    ax.set_xlabel("Reserve score on the first scan (control SD)")
    m0 = er("Primary: composite score (bone + muscle + lung), per 1 SD lower")
    style.source_data(sc[["pid", "set", "group", "S"]], "P2_Fig2", "b")

    # ================= c: quartiles, two tiers
    lab(fig, 0.004, 0.565, "c")
    q = qt.set_index("quartile").loc[["Q4", "Q3", "Q2", "Q1"]]
    xx = np.arange(4)
    ramp = ["#f6d3c3", "#efa786", "#e7794f", CASE]
    ax = fig.add_axes([0.075, 0.43, 0.29, 0.12])
    ax.plot(xx, q["or"].values, color=CASE, lw=0.8, zorder=2)
    for i, (qq, r) in enumerate(q.iterrows()):
        if qq != "Q4":
            ax.plot([i, i], [r.lo, r.hi], color=CASE, lw=1.0, solid_capstyle="butt")
        ax.plot(i, r["or"], "o", ms=4.2, color=ramp[i], mec=CASE, mew=0.8, zorder=3)
        if qq == "Q4":
            ax.text(i, 0.95, "ref.", fontsize=5.2, ha="center", va="top", color=INK2)
    ax.axhline(1, color=INK2, lw=0.6)
    ax.set_yscale("log"); ax.set_ylim(0.85, 2.4)
    ax.set_yticks([1, 1.5, 2]); ax.set_yticklabels(["1", "1.5", "2"]); ax.minorticks_off()
    ax.set_xticks(xx); ax.set_xticklabels([]); ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel("Odds ratio")
    ax2 = fig.add_axes([0.075, 0.335, 0.29, 0.085])
    share = 100 * q.n_cases / q.n_cases.sum()
    ax2.bar(xx, share.values, width=0.62, color=ramp, edgecolor=CASE, lw=0.5)
    ax2.axhline(25, color=INK2, lw=0.6, ls=(0, (3, 2)))
    ax2.set_ylim(0, 38); ax2.set_yticks([0, 25]); ax2.set_ylabel("Cases (%)")
    ax2.set_xlim(-0.5, 3.5); ax2.set_xticks(xx)
    ax2.set_xticklabels([f"{qq}\n{int(r.n_cases)}/{int(r.n_controls)}" for qq, r in q.iterrows()], fontsize=5.5)
    ax2.set_xlabel("Control quartile of the score, highest (Q4) to lowest (Q1); cases/controls", fontsize=5.8)
    style.source_data(qt, "P2_Fig2", "c")

    # ================= d: consistency across case series and stability
    lab(fig, 0.44, 0.565, "d")
    ax = fig.add_axes([0.665, 0.335, 0.165, 0.2])
    rows = [dict(kind="primary", label="All matched sets", n=int(m0.n_cases), **pick(m0))]
    rows.append(dict(kind="header", label="Case series",
                     pint=est[est.analysis.str.startswith("Series:")].p_interaction.iloc[0]))
    for s in ["Resected", "Not resected"]:
        r = er(f"Series: {s}"); rows.append(dict(kind="sub", label=s, n=int(r.n_cases), **pick(r)))
    rows.append(dict(kind="header", label="First scan to diagnosis",
                     pint=est[est.analysis.str.startswith("Lead time:")].p_interaction.iloc[0]))
    for s in ["<2 y", "2–<4 y", "≥4 y"]:
        r = er(f"Lead time: {s}")
        rows.append(dict(kind="sub", label=s.replace("2–<4", "2–4"), n=int(r.n_cases), **pick(r)))
    r = er("Composite score at latest pre-diagnostic scan, per 1 SD lower")
    rows.append(dict(kind="row", label="Latest pre-diagnostic scan", n=int(r.n_cases), **pick(r)))
    forest(ax, rows, xlim=(0.85, 1.75), ticks=[1, 1.25, 1.5], label_w=1.3,
           right_cols=(("OR (95% CI)", 1.05), ("Cases", 1.85)), xlabel="OR per 1 SD lower score")
    style.source_data(pd.DataFrame(rows), "P2_Fig2", "d")

    # ================= e: components and score variants
    lab(fig, 0.004, 0.27, "e")
    ax = fig.add_axes([0.25, 0.035, 0.15, 0.215])
    LAA = "laa950"
    comp = [("Reserve score (3 measures)", "Primary: composite score (bone + muscle + lung), per 1 SD lower",
             "primary", None),
            ("Single measures", None, "header", None),
            ("Lung density", "Component: lung_perc15_hu, per 1 SD lower reserve", "row", "lung"),
            ("Trabecular bone", "Component: T10_trab_hu, per 1 SD lower reserve", "row", "bone"),
            ("Muscle volume", "Component: psm_lvol, per 1 SD lower reserve", "row", "muscle"),
            ("Muscle radiodensity", "Component: psm_hu, per 1 SD lower reserve", "row", "muscle_h"),
            ("Intramuscular fat (higher)", "Component: psm_fat_pct, per 1 SD lower reserve", "row", "fat"),
            ("Bone + muscle score", None, "header", None),
            ("Alone", "Musculoskeletal score (bone + muscle), per 1 SD lower", "sub", None),
            ("Adjusted for lung p15", "Musculoskeletal score, adjusted for lung density z", "sub", None),
            ("Adjusted for % lung < −950 HU", LAA, "sub", None),
            ("Other", None, "header", None),
            ("Lung density, adj. bone + muscle", "Lung density z, adjusted for musculoskeletal score", "sub", "lung"),
            ("Five-measure score", "Five-measure score (adds intramuscular fat, lean radiodensity), per 1 SD lower",
             "sub", None)]
    # %LAA-adjusted bone + muscle: sensitivity analysis 6 of docs/paper2/P2_04 (stats/p2_review_sens.R)
    rv = rd("rev_sens_estimates.csv")
    rv = rv[rv.analysis.str.startswith("6.") & (rv.model == "Musculoskeletal score adjusted for lung z")]
    assert len(rv) == 1
    rows = []
    for labl, name, kind, tis in comp:
        if kind == "header":
            rows.append(dict(kind="header", label=labl)); continue
        d = dict(kind=kind, label=labl, **pick(rv.iloc[0] if name == LAA else er(name)))
        if tis:
            t = tis.replace("_h", "")
            d.update(color=TIS[t], hollow=tis.endswith("_h"))
        rows.append(d)
    forest(ax, rows, xlim=(0.85, 1.75), ticks=[1, 1.25, 1.5], label_w=1.55, right_cols=(("OR (95% CI)", 1.05),),
           xlabel="OR per 1 SD lower reserve")
    style.source_data(pd.DataFrame(rows), "P2_Fig2", "e")

    # ================= f: subgroups
    lab(fig, 0.5, 0.27, "f")
    ax = fig.add_axes([0.735, 0.035, 0.11, 0.215])
    rows = []
    for fac, levels in [("Histology", ["Adenocarcinoma", "Squamous", "Small cell", "Other"]),
                        ("Stage", ["Stage I–II", "Stage III–IV"]),
                        ("Detection", ["Screen-detected", "Interval/post-screening"]),
                        ("Sex", ["Men", "Women"]), ("Smoking", ["Current smoker", "Former smoker"])]:
        sub = est[est.analysis.str.startswith(fac + ":")]
        rows.append(dict(kind="header", label=fac, pint=sub.p_interaction.iloc[0]))
        for lv in levels:
            r = sub[sub.analysis == f"{fac}: {lv}"].iloc[0]
            rows.append(dict(kind="sub", label=lv.replace("Stage ", "").replace("Interval/post-screening",
                                                                                 "Interval/post-screen"),
                             n=int(r.n_cases), **pick(r)))
    forest(ax, rows, xlim=(0.75, 1.8), ticks=[0.8, 1, 1.25, 1.5], label_w=1.8,
           right_cols=(("OR (95% CI)", 1.05), ("Cases", 2.22)), xlabel="OR per 1 SD lower score", fs=5.4)
    style.source_data(pd.DataFrame(rows), "P2_Fig2", "f")
    from p2_figs_extra_rows import row_fig2
    row_fig2(sf2)
    return style.save(root, "P2_Fig2")


if __name__ == "__main__":
    print(fig2())
