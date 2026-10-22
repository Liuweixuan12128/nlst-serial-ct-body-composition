"""Shared design layer for the paper-2 result figures (2-6): tissue colours, 3D tissue icons, journal-style forest
plots, rainclouds and ribbons. Built on style.py; every function draws real numbers passed in by the caller.

Tissue palette (dataviz validator, white surface, --pairs all: worst CVD dE 12.0, normal-vision dE 19.0; amber is
below 3:1 contrast, so amber marks are always labelled): bone amber, muscle red, intramuscular fat violet, lung teal.
Lean-muscle radiodensity shares the muscle hue and is told apart by hollow markers and dashed lines.
Icons are small renders of the example scan (figures/p2_render_reserve3d.py).
"""
import os

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib.patches import Circle, Rectangle
from matplotlib.transforms import blended_transform_factory
from scipy.stats import gaussian_kde

from p2_cns_common import CASE, CTRL, HAIR, INK, INK2, MUTED, ROOT, crop_alpha

TIS = {"bone": "#d99000", "muscle": "#d63c3c", "fat": "#5b45b0", "lung": "#11968f"}
MEAS_TIS = {"T10_trab_hu": "bone", "psm_lvol": "muscle", "psm_hu": "muscle", "psm_fat_pct": "fat",
            "lung_perc15_hu": "lung"}
HOLLOW = {"psm_hu"}
ZEBRA = "#f3f2ee"
BAND = "#f1f0eb"   # equivalence / reference bands
GRID = "#e6e5df"
_ICONS = {}


def mcol(m):
    return TIS[MEAS_TIS[m]]


def mstyle(m):
    """marker face colour and line style for a measure (radiodensity hollow/dashed)"""
    return ("white", (0, (3, 1.6))) if m in HOLLOW else (mcol(m), "-")


def icon_img(tissue, h_px=220):
    if tissue not in _ICONS:
        from PIL import Image
        im = crop_alpha(os.path.join(ROOT, "figures", "assets", f"p2_icon_{tissue}.png"), pad=6)
        w = int(im.width * h_px / im.height)
        _ICONS[tissue] = np.asarray(im.resize((w, h_px), Image.LANCZOS))
    return _ICONS[tissue]


def add_icon(ax, tissue, xy, height_pt=11, xycoords="data", box_alignment=(0.5, 0.5), zorder=5):
    img = icon_img(tissue)
    oi = OffsetImage(img, zoom=height_pt / img.shape[0], interpolation="lanczos")
    ab = AnnotationBbox(oi, xy, xycoords=xycoords, frameon=False, box_alignment=box_alignment, pad=0,
                        annotation_clip=False, zorder=zorder)
    ax.add_artist(ab)
    return ab


def forest(ax, rows, xlim, ticks, label_w=1.15, right_cols=(("OR (95% CI)", 1.04),), ref=1.0, log=True,
           xlabel=None, header=True, fs=5.6, zebra=False, ms_range=(2.0, 4.6)):
    """Journal-style forest plot drawn from the top.
    rows: dicts with kind ('header' | 'primary' | 'row' | 'sub'), label, and for data rows est, lo, hi;
          optional color, hollow, icon (tissue), n, pint (header rows), txt (right-column text override).
    label_w: width of the label column in axes-width units (labels are right-aligned at x = -0.03).
    Marker area scales with inverse variance (log scale) so precise estimates read heavier."""
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    data = [r for r in rows if r["kind"] != "header"]
    se = [(np.log(r["hi"]) - np.log(r["lo"])) / 3.92 if log else (r["hi"] - r["lo"]) / 3.92 for r in data]
    wmax = max(1 / s ** 2 for s in se)
    y, k = 0, 0
    for r in rows:
        if zebra and r["kind"] != "header":
            if k % 2 == 0:
                ax.add_patch(Rectangle((-label_w, y - 0.5), label_w + 1 + 0.42, 1, transform=tr, color=ZEBRA, lw=0,
                                       zorder=0, clip_on=False))
            k += 1
        if r["kind"] == "header":
            ax.text(-label_w + 0.02, y, r["label"], transform=tr, ha="left", va="center", fontsize=fs,
                    fontweight="bold", color=INK)
            if r.get("pint") is not None:
                ax.text(right_cols[0][1], y, rf"$P_{{\rm interaction}}$ = {r['pint']:.2f}", transform=tr, ha="left",
                        va="center", fontsize=fs - 0.6, color=INK2)
            k = 0
        else:
            prim = r["kind"] == "primary"
            col = r.get("color", CASE if prim else INK)
            s = (np.log(r["hi"]) - np.log(r["lo"])) / 3.92 if log else (r["hi"] - r["lo"]) / 3.92
            ms = ms_range[0] + (ms_range[1] - ms_range[0]) * np.sqrt((1 / s ** 2) / wmax)
            ax.plot([r["lo"], r["hi"]], [y, y], color=col, lw=1.1 if prim else 0.85, solid_capstyle="butt", zorder=3)
            if prim:
                ax.plot(r["est"], y, marker="D", ms=ms + 1.2, color=col, mec="white", mew=0.5, zorder=4)
            else:
                ax.plot(r["est"], y, marker="s", ms=ms, color=col, mfc="white" if r.get("hollow") else col,
                        mec=col, mew=0.8, zorder=4)
            x_lab = -0.03
            if r.get("icon"):   # icons form their own column at the left edge of the label area
                add_icon(ax, r["icon"], (-label_w + 0.03, y), height_pt=8.5, xycoords=tr, box_alignment=(0, 0.5))
            indent = "   " if r["kind"] == "sub" else ""
            ax.text(x_lab, y, indent + r["label"], transform=tr, ha="right", va="center", fontsize=fs,
                    color=INK, fontweight="bold" if prim else "normal")
            txt = r.get("txt") or f"{r['est']:.2f} ({r['lo']:.2f}–{r['hi']:.2f})"
            ax.text(right_cols[0][1], y, txt, transform=tr, ha="left", va="center", fontsize=fs - 0.4,
                    color=INK if prim else INK2, fontweight="bold" if prim else "normal",
                    family="sans-serif")
            if len(right_cols) > 1 and r.get("n") is not None:
                ax.text(right_cols[1][1], y, f"{r['n']:,}", transform=tr, ha="right", va="center",
                        fontsize=fs - 0.4, color=INK2)
        y -= 1
    if header:
        for name, xpos in right_cols:
            ha = "right" if name in ("n", "Cases") else "left"
            ax.text(xpos, 1.05, name, transform=tr, ha=ha, va="center", fontsize=fs - 0.4, color=INK2,
                    fontweight="bold")
    ax.axvline(ref, color=INK2, lw=0.6, zorder=2)
    for t in ticks:
        if t != ref:
            ax.axvline(t, color=GRID, lw=0.5, zorder=1)
    if log:
        ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.set_ylim(y + 0.45, 0.55 + (0.6 if header else 0))
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)


def raincloud(ax, groups, colors, labels, ypos, bw=0.25, cloud_h=0.42, rng_seed=1, point_alpha=0.18, xlim=None):
    """horizontal raincloud: half-violin above, box (IQR + median) on the line, jittered points below"""
    rng = np.random.default_rng(rng_seed)
    for v, c, lab, y0 in zip(groups, colors, labels, ypos):
        v = np.asarray(v)
        xs = np.linspace(*(xlim or (v.min(), v.max())), 400)
        k = gaussian_kde(v, bw_method=bw)(xs)
        k = k / k.max() * cloud_h
        ax.fill_between(xs, y0 + 0.06, y0 + 0.06 + k, color=c, alpha=0.28, lw=0)
        ax.plot(xs, y0 + 0.06 + k, color=c, lw=0.9)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.add_patch(Rectangle((q1, y0 - 0.035), q3 - q1, 0.07, fc="white", ec=c, lw=0.8, zorder=4))
        ax.plot([med, med], [y0 - 0.035, y0 + 0.035], color=c, lw=1.4, zorder=5, solid_capstyle="butt")
        lo, hi = np.percentile(v, [2.5, 97.5])
        ax.plot([lo, q1], [y0, y0], color=c, lw=0.7, zorder=3)
        ax.plot([q3, hi], [y0, y0], color=c, lw=0.7, zorder=3)
        jy = y0 - 0.1 - rng.uniform(0, 0.22, len(v))
        ax.scatter(v, jy, s=0.6, color=c, alpha=point_alpha, lw=0, rasterized=True, zorder=2)
        ax.plot(v.mean(), y0, marker="o", ms=3.2, mfc="white", mec=c, mew=0.9, zorder=6)


def ribbon(ax, x, y, lo, hi, color, ls="-", mfc=None, label=None, ms=2.8, alpha=0.14):
    ax.fill_between(x, lo, hi, color=color, alpha=alpha, lw=0)
    ax.plot(x, y, color=color, lw=1.0, ls=ls)
    ax.plot(x, y, "o", ms=ms, color=color, mfc=mfc or color, mec=color, mew=0.8, label=label, zorder=4)


def panel_title(ax, text, y=1.1, fs=6.2):
    """plain title, left-aligned above the axes"""
    ax.text(0.0, y, text, transform=ax.transAxes, fontsize=fs, va="center", color=INK)


def leader(fig, xy_from, xy_to, text, color, ha="left", fs=5.8, weight="bold"):
    """callout from a point on an image (figure fraction) to a label (figure fraction)"""
    tr = getattr(fig, "transSubfigure", fig.transFigure)   # works for a Figure and a SubFigure
    fig.add_artist(Line2D([xy_from[0], xy_to[0]], [xy_from[1], xy_to[1]], transform=tr, color=color, lw=0.6))
    fig.add_artist(Circle(xy_from, 0.0028, transform=tr, color=color, zorder=10))
    fig.text(xy_to[0] + (0.004 if ha == "left" else -0.004), xy_to[1], text, ha=ha, va="center", fontsize=fs,
             color=color, fontweight=weight)
