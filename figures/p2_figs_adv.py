"""Paper 2 advanced figures: Fig. 7 spatial maps of intramuscular fat; Fig. 8 aligned raster and detectability;
Extended Data Fig. 9 specification curve.  Inputs: stats/p2_map_stats.py, stats/p2_power.py outputs.
Usage: python figures/p2_figs_adv.py [fig7 fig8 ed9]
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

from p2_cns_common import (CASE, CTRL, HAIR, INK, INK2, MM, MUTED, RES, ROOT, WASH, axial, dot_ci, fnum, imshow_ct,
                           rd, style)
from p2_premium import BAND, TIS
from style import BLUE_RAMP, W2

DIV = LinearSegmentedColormap.from_list("div", ["#184f95", "#3987e5", "#9ec5f4", "#efeeea", "#f6b896", "#eb6834",
                                                "#a8401a"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#f4f3f8", "#cfc9e6", "#9a8fcf", "#6a5bb3", "#4a3aa7", "#2c2170"])


def lab(fig, ax, letter, dx=0.0, dy=0.0):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=8, fontweight="bold", va="bottom", ha="left", color=INK)


def bilateral(m, support):
    """mirror the medial->lateral map onto both sides of a midline gap (anterior up)"""
    m = np.where(support, m, np.nan)
    gap = np.full((m.shape[0], 2), np.nan)
    return np.hstack([np.fliplr(m), gap, m])


def show_map(ax, m, support, cmap, norm=None, vmin=None, vmax=None, sig=None, orient=True):
    """orient: write the orientation words (only on the first map; the others share its frame)"""
    b = bilateral(m, support)
    im = ax.imshow(b, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, interpolation="nearest", aspect="equal")
    if sig is not None:
        s = bilateral(np.where(np.isnan(sig), 0, sig).astype(float), support)
        s = np.nan_to_num(s)
        if s.any():
            ax.contour(s, levels=[0.5], colors=INK, linewidths=0.7)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    w = b.shape[1]
    if orient:
        ax.text(w / 2 - 0.5, -1.2, "midline", fontsize=4.8, color=INK2, ha="center", va="bottom")
        ax.text(-0.5, b.shape[0] + 0.6, "lateral", fontsize=4.8, color=INK2, ha="left", va="top")
        ax.text(w - 0.5, b.shape[0] + 0.6, "lateral", fontsize=4.8, color=INK2, ha="right", va="top")
        ax.text(-1.2, 0, "anterior", fontsize=4.8, color=INK2, rotation=90, ha="right", va="top")
        ax.text(-1.2, b.shape[0] - 1, "posterior", fontsize=4.8, color=INK2, rotation=90, ha="right", va="bottom")
    return im


def cbar(fig, im, ax, label, ticks=None):
    bb = ax.get_position()
    cax = fig.add_axes([bb.x0 + bb.width * 0.15, bb.y0 - 0.035, bb.width * 0.7, 0.012])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=ticks)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=5, length=2, width=0.4)
    cb.set_label(label, fontsize=5.5, labelpad=1.5)
    return cb


# ====================================================================== Figure 7: spatial maps
def fig7():
    style.setup()
    from p2_cns_common import crop_alpha
    z = np.load(os.path.join(RES, "map_stats.npz"))
    sup = z["support"]
    prof = rd("map_level_profiles.csv")
    FH = 150
    fig = plt.figure(figsize=(W2, FH * MM))
    vage = float(np.nanpercentile(z["age_b"][sup], 98))
    vcc = float(np.nanpercentile(np.abs(z["cc_d"][sup]), 99))
    vch = float(np.nanpercentile(np.abs(z["ch_d"][sup]), 99))
    # a: normalisation schematic on the example scan
    ax = fig.add_axes([0.01, 0.69, 0.19, 0.25])
    img, msk, px = axial("100658_T0", "T8", half_mm=62, shift_mm=0)
    imshow_ct(ax, img, msk)
    ax.set_facecolor("black")
    ys, xs = np.nonzero(msk)
    mid = img.shape[1] / 2
    ax.annotate("", xy=(xs.max(), ys.max() + 6), xytext=(mid + 3, ys.max() + 6),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=0.7, mutation_scale=5))
    ax.text((mid + xs.max()) / 2, ys.max() + 9, "u: medial to lateral", color="white", fontsize=5, ha="center", va="top")
    ax.annotate("", xy=(xs.max() + 6, ys.max()), xytext=(xs.max() + 6, ys.min()),
                arrowprops=dict(arrowstyle="-|>", color="white", lw=0.7, mutation_scale=5))
    ax.text(xs.max() + 8, (ys.min() + ys.max()) / 2, "v", color="white", fontsize=5, va="center")
    ax.set_title("Shape-normalised frame", fontsize=6, pad=3)
    lab(fig, ax, "a", dy=0.045)
    # b: control mean fat map
    ax = fig.add_axes([0.24, 0.73, 0.22, 0.16])
    im = show_map(ax, z["mean_fat"], sup, SEQ, vmin=0, vmax=np.nanpercentile(z["mean_fat"][sup], 98))
    ax.set_title("Controls: mean fat", fontsize=6, pad=8)
    cbar(fig, im, ax, "fat (% of voxels)")
    lab(fig, ax, "b", dx=-0.025, dy=0.035)
    # c: age map
    ax = fig.add_axes([0.5, 0.73, 0.22, 0.16])
    im = show_map(ax, z["age_b"], sup, SEQ, vmin=0, vmax=vage, sig=z["age_p"] < 0.05, orient=False)
    ax.set_title("Age effect per 10 years", fontsize=6, pad=8)
    cbar(fig, im, ax, "pp per 10 years")
    lab(fig, ax, "c", dx=-0.025, dy=0.035)
    # d: age map projected on 3D anatomy
    ax = fig.add_axes([0.75, 0.64, 0.24, 0.33])
    ax.imshow(crop_alpha(os.path.join(ROOT, "figures", "assets", "p2_maps3d_age.png"), pad=4)); ax.axis("off")
    ax.set_title("Age effect, 3D", fontsize=6, pad=2)
    lab(fig, ax, "d", dy=-0.005)
    # e: case-control level difference
    ax = fig.add_axes([0.04, 0.42, 0.24, 0.16])
    im = show_map(ax, z["cc_d"], sup, DIV, norm=TwoSlopeNorm(0, -vcc, vcc), sig=z["cc_p"] < 0.05, orient=False)
    ax.set_title("Case \u2212 control: level", fontsize=6, pad=8)
    cbar(fig, im, ax, "pp (case \u2212 controls)")
    lab(fig, ax, "e", dx=-0.03, dy=0.035)
    # f: case-control level difference in 3D
    ax = fig.add_axes([0.3, 0.34, 0.22, 0.31])
    ax.imshow(crop_alpha(os.path.join(ROOT, "figures", "assets", "p2_maps3d_cc.png"), pad=4)); ax.axis("off")
    ax.set_title("Level difference, 3D", fontsize=6, pad=2)
    lab(fig, ax, "f", dy=-0.005)
    # g: case-control change difference
    ax = fig.add_axes([0.56, 0.42, 0.24, 0.16])
    im = show_map(ax, z["ch_d"], sup, DIV, norm=TwoSlopeNorm(0, -vch, vch), sig=z["ch_p"] < 0.05, orient=False)
    ax.set_title("Case \u2212 control: annual change", fontsize=6, pad=8)
    cbar(fig, im, ax, "pp per year (case \u2212 controls)")
    lab(fig, ax, "g", dx=-0.03, dy=0.035)
    # h: controls' observed annual change
    ax = fig.add_axes([0.83, 0.42, 0.16, 0.16])
    vo = float(np.nanpercentile(np.abs(z["ctrl_change"][sup]), 98))
    im = show_map(ax, z["ctrl_change"], sup, DIV, norm=TwoSlopeNorm(0, -vo, vo), orient=False)
    ax.set_title("Controls: annual change", fontsize=6, pad=8)
    cbar(fig, im, ax, "pp per year")
    lab(fig, ax, "h", dx=-0.02, dy=0.035)
    # i: craniocaudal profile
    ax = fig.add_axes([0.08, 0.07, 0.36, 0.2])
    x = np.arange(len(prof))
    ax.fill_between(x, prof.ctrl_lo, prof.ctrl_hi, color=CTRL, alpha=0.12, lw=0)
    ax.plot(x, prof.ctrl_mean, color=CTRL, lw=1.1, label="controls (mean, IQR)")
    ax.plot(x, prof.case_mean, color=CASE, lw=1.1, label="cases (mean)")
    ax.set_xticks(x); ax.set_xticklabels(["T11", "", "", "", "", "T6"])
    ax.set_xlabel("Craniocaudal position (6 bins)"); ax.set_ylabel("Intramuscular fat (%)")
    ax.legend(loc="upper left", fontsize=5.2)
    lab(fig, ax, "i", dx=-0.07, dy=0.02)
    # j: case-control difference by level
    ax = fig.add_axes([0.58, 0.07, 0.38, 0.2])
    dot_ci(ax, x, prof["diff"].values, prof.diff_lo.values, prof.diff_hi.values, CASE, horizontal=False)
    ax.axhline(0, color=HAIR, lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(["T11", "", "", "", "", "T6"])
    ax.set_xlabel("Craniocaudal position (6 bins)"); ax.set_ylabel("Case \u2212 controls (pp)")
    lab(fig, ax, "j", dx=-0.07, dy=0.02)
    for k, arr in [("mean_fat", z["mean_fat"]), ("age_b", z["age_b"]), ("cc_d", z["cc_d"]), ("ch_d", z["ch_d"]),
                   ("ctrl_change", z["ctrl_change"])]:
        style.source_data(pd.DataFrame(arr), "P2_Fig7", k)
    style.source_data(prof, "P2_Fig7", "ij")
    return style.save(fig, "P2_Fig7")


# ====================================================================== Figure 8: aligned raster and detectability
def fig8():
    style.setup()
    iv = rd("raster_intervals.csv").astype({"pid": str})
    pg = rd("power_grid.csv")
    fig = plt.figure(figsize=(W2, 135 * MM))
    # a: raster of every case interval aligned to diagnosis
    cas = iv[iv.group == "case"].copy()
    order = cas.groupby("pid").end_tau.max().sort_values().index
    yi = {p: i for i, p in enumerate(order)}
    ax = fig.add_axes([0.07, 0.07, 0.38, 0.66])
    from matplotlib.collections import LineCollection
    segs = [[(r.tau0, yi[r.pid]), (r.end_tau, yi[r.pid])] for r in cas.itertuples()]
    lc = LineCollection(segs, cmap=DIV, norm=TwoSlopeNorm(0, -2, 2), linewidths=0.55)
    lc.set_array(np.clip(cas.z.values, -2, 2))
    ax.add_collection(lc)
    ax.axvline(0, color=INK, lw=0.6)
    ax.axvspan(-1, 0, color=BAND, lw=0, zorder=0)
    ax.set_xlim(-6.2, 0.15); ax.set_ylim(-5, len(order) + 5)
    ax.set_xlabel("Years relative to lung cancer diagnosis")
    ax.set_ylabel(f"Cases (n = {len(order)}), ordered by time from last screen to diagnosis")
    ax.set_yticks([])
    ax.text(-0.5, len(order) + 8, "<1 y", fontsize=5.5, color=INK2, ha="center")
    cax = ax.inset_axes([0.47, 0.1, 0.4, 0.018])
    cb = fig.colorbar(lc, cax=cax, orientation="horizontal", ticks=[-2, -1, 0, 1, 2])
    cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5, length=2)
    cb.set_label("Interval change (z)", fontsize=5.2, labelpad=2)
    cb.ax.xaxis.set_label_position("top")
    fig.text(0.005, 0.985, "a", fontsize=8, fontweight="bold", va="top")
    # marginal: binned mean z by interval end, cases and aligned controls
    axm = fig.add_axes([0.07, 0.76, 0.38, 0.17], sharex=ax)
    bins = np.arange(-6, 0.01, 0.5)
    for g, col in [("control", CTRL), ("case", CASE)]:
        q = iv[iv.group == g]
        mids, mu, lo, hi = [], [], [], []
        for a, b in zip(bins[:-1], bins[1:]):
            s = q[(q.end_tau > a) & (q.end_tau <= b)].z
            if len(s) >= 20:
                mids.append((a + b) / 2); mu.append(s.mean()); se = s.std() / np.sqrt(len(s))
                lo.append(s.mean() - 1.96 * se); hi.append(s.mean() + 1.96 * se)
        axm.fill_between(mids, lo, hi, color=col, alpha=0.18, lw=0)
        axm.plot(mids, mu, color=col, lw=1, label="cases" if g == "case" else "matched controls (aligned)")
    axm.axhline(0, color=HAIR, lw=0.5); axm.axvline(0, color=INK, lw=0.6)
    axm.set_ylim(-0.45, 0.45); axm.set_ylabel("mean z")
    axm.tick_params(labelbottom=False)
    axm.legend(loc="upper left", fontsize=5.2, ncol=2, bbox_to_anchor=(0, 1.25))
    style.source_data(cas[["pid", "tau0", "end_tau", "z"]], "P2_Fig8", "a")
    # b: detectability heatmap (intramuscular fat)
    ax = fig.add_axes([0.56, 0.56, 0.42, 0.37])
    q = pg[pg.measure == "psm_fat_pct"].pivot(index="onset_years", columns="multiple_of_normal", values="power")
    ks = q.columns.values; Ls = q.index.values
    im = ax.imshow(q.values, cmap=LinearSegmentedColormap.from_list("pw", ["#f4f3f8"] + BLUE_RAMP[1:]), vmin=0, vmax=1,
                   aspect="auto", origin="lower")
    for i in range(len(Ls)):
        for j in range(len(ks)):
            v = q.values[i, j]
            ax.text(j, i, "1" if v >= 0.9995 else (">.99" if v >= 0.995 else f"{v:.2f}".lstrip("0")), ha="center",
                    va="center", fontsize=5,
                    color="white" if v > 0.6 else INK)
    ax.contour(q.values, levels=[0.8], colors=CASE, linewidths=1.0)
    ax.set_xticks(range(len(ks))); ax.set_xticklabels([f"{k:g}×" for k in ks], fontsize=5.5)
    ax.set_yticks(range(len(Ls))); ax.set_yticklabels([f"{L:g}" for L in Ls])
    ax.set_xlabel("Injected extra fat gain (× normal annual rate)")
    ax.set_ylabel("Onset before diagnosis (years)")
    lab(fig, ax, "b", dx=-0.06, dy=0.01)
    style.source_data(pg, "P2_Fig8", "b")
    # c: power curves at onset 1 year, fat and volume
    ax = fig.add_axes([0.56, 0.08, 0.18, 0.33])
    for m, col, nm in [("psm_fat_pct", TIS["fat"], "Intramuscular fat"), ("psm_lvol", TIS["muscle"], "Muscle volume")]:
        q = pg[(pg.measure == m) & (pg.onset_years == 1.0)]
        ax.plot(q.multiple_of_normal, q.power, "o-", color=col, ms=2.5, lw=1, label=nm)
    ax.axhline(0.8, color=HAIR, lw=0.6, ls=(0, (3, 2)))
    ax.axhline(0.025, color=HAIR, lw=0.5)
    ax.set_xlabel("Multiple of normal rate"); ax.set_ylabel("Power (onset 1 y)")
    ax.set_xlim(-0.3, 6.3); ax.set_ylim(0, 1.03)
    ax.legend(loc="lower right", fontsize=5.2)
    lab(fig, ax, "c", dx=-0.06, dy=0.02)
    # d: observed estimate against injected shifts
    ax = fig.add_axes([0.82, 0.08, 0.16, 0.33])
    div = rd("partA_divergence.csv")
    obs = div[(div.measure == "psm_fat_pct") & (div.term == "windowcase: ends <1 y before dx")].iloc[0]
    q = pg[(pg.measure == "psm_fat_pct") & (pg.onset_years == 1.0)]
    ax.plot(q.multiple_of_normal, q.injected_shift, color=INK2, lw=1, marker="o", ms=2.2, label="expected shift")
    ax.axhspan(obs.lo, obs.hi, color=CASE, alpha=0.15, lw=0)
    ax.axhline(obs.estimate, color=CASE, lw=1, label="observed (95% CI)")
    ax.set_xlabel("Multiple of normal rate"); ax.set_ylabel("<1-y window z")
    ax.set_xlim(-0.3, 6.3)
    ax.legend(loc="upper left", fontsize=5)
    lab(fig, ax, "d", dx=-0.06, dy=0.02)
    return style.save(fig, "P2_Fig8")


# ====================================================================== Extended Data Fig. 9: specification curve
def ed9():
    style.setup()
    sc = rd("spec_curve.csv").sort_values("estimate").reset_index(drop=True)
    fig = plt.figure(figsize=(W2, 100 * MM))
    ax = fig.add_axes([0.2, 0.52, 0.78, 0.42])
    x = np.arange(len(sc))
    ax.axhspan(-0.2, 0.2, color=BAND, lw=0, zorder=0)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.vlines(x, sc.lo, sc.hi, color=TIS["fat"], lw=0.6, alpha=0.7)
    ax.plot(x, sc.estimate, "o", ms=1.8, color=TIS["fat"])
    ax.set_xlim(-1, len(sc)); ax.set_xticks([])
    ax.set_ylabel("Late-window excess z\n(intramuscular fat)")
    fig.text(0.005, 0.985, "a", fontsize=8, fontweight="bold", va="top")
    rows = [("segmentation", "aut", "Second segmentation model"), ("noise_adjusted", False, "No noise adjustment"),
            ("covariate_reference", False, "Unadjusted control reference"), ("kernel_consistent_only", True, "Same kernel only"),
            ("exclude_scans_90d", True, "Exclude scans <90 d before dx"), ("window_cut_years", 0.5, "Late window <0.5 y"),
            ("window_cut_years", 1.5, "Late window <1.5 y")]
    axb = fig.add_axes([0.2, 0.06, 0.78, 0.4], sharex=ax)
    for i, (col, val, name) in enumerate(rows):
        on = (sc[col] == val).values
        axb.scatter(x[on], np.full(on.sum(), -i), s=3, color=INK, lw=0)
        axb.scatter(x[~on], np.full((~on).sum(), -i), s=3, color=HAIR, lw=0)
    axb.set_yticks(-np.arange(len(rows))); axb.set_yticklabels([r[2] for r in rows], fontsize=5.8)
    axb.set_ylim(-len(rows) + 0.5, 0.5)
    axb.set_xticks([]); axb.spines["bottom"].set_visible(False)
    axb.set_xlabel("Specification")
    fig.text(0.005, 0.47, "b", fontsize=8, fontweight="bold", va="top")
    style.source_data(sc, "P2_EDFig9", "a")
    return style.save(fig, "P2_EDFig9")


if __name__ == "__main__":
    for t in (sys.argv[1:] or ["fig7", "fig8", "ed9"]):
        print(globals()[t]())
