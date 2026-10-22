"""Paper 2 main figures and analysis Extended Data figures, CNS layout.
Display numbers (after the 2026-10-11 reordering; function names keep the original order):
  fig1 -> Fig. 1, fig6 -> Fig. 3 (stable level offset), fig2 -> Fig. 4 (normative), fig3 -> Fig. 5 (reliability),
  fig4 -> Fig. 6 (no acceleration), fig5 -> Extended Data Fig. 10 (subgroups); ED4, ED5, ED7, ED8 unchanged.
Fig. 2 (reserve score) is in p2_figs_pheno.py.
All panels read results/p2 outputs of stats/p2_partA.R and stats/p2_extra.R; nothing is typed by hand.
Usage: python figures/p2_figs_main.py [fig1 ... fig6 ed4 ed5 ed7 ed8]
"""
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from p2_cns_common import (ALL, CASE, CTRL, HAIR, INK, INK2, LABEL, MEAS, MM, MUTED, RES, ROOT, UNIT, WASH, axial,
                           crop_alpha, dot_ci, fnum, imshow_ct, primary, rd, scalebar, sel, style)
from style import BLUE_RAMP, W2, panel_label
from p2_premium import BAND, HOLLOW, TIS, forest, mcol, mstyle, panel_title, raincloud, ribbon

SHORT = {"psm_fat_pct": "Intramuscular\nfat", "psm_hu": "Muscle\nradiodensity", "psm_lvol": "Muscle\nvolume",
         "T10_trab_hu": "Trabecular\nbone (T10)", "lung_perc15_hu": "Lung\ndensity"}
WIN = ["case: ends >=2 y before dx", "case: ends 1-2 y before dx", "case: ends <1 y before dx"]
WIN_LAB = ["≥2", "1–2", "<1"]


def long_data():
    return pd.read_csv(os.path.join(RES, "partA_long.csv"), dtype={"pid": str, "set": str})


def control_sd():
    d = long_data()
    first = d[d.group == "control"].sort_values("t").groupby("pid").head(1)
    return {v: first[v].std() for v in ALL}


def lab(fig, ax, letter, dx=-0.0, dy=0.0):
    """panel letter at the top-left corner of an axes, in figure coordinates (keeps letters aligned)"""
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, letter, fontsize=8, fontweight="bold", va="bottom", ha="left", color=INK)


# ====================================================================== Figure 1
def fig1():
    style.setup()
    flow = rd("p2_flow_counts.csv").set_index("step").n
    cnt = rd("partA_counts.csv").set_index("item").n
    p = primary()
    FH = 130
    st = 138 / FH          # positions below were laid out for 138 mm; this keeps sizes and closes the gap between rows
    top = lambda y: 1 - (1 - y) * st
    fig = plt.figure(figsize=(W2, FH * MM))
    # a: flow, CONSORT style (white boxes, thin outline, straight arrows, exclusions in side boxes)
    ax = fig.add_axes([0.0, top(0.47), 0.46, 0.52 * st]); ax.axis("off")
    from matplotlib.patches import Rectangle

    def box(x, y, w, h, text, fs=5.8):
        ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc="white", ec=INK, lw=0.6, transform=ax.transAxes))
        ax.text(x, y, text, ha="center", va="center", fontsize=fs, transform=ax.transAxes, linespacing=1.25)

    def line(xs, ys):
        ax.plot(xs, ys, color=INK, lw=0.6, transform=ax.transAxes, solid_capstyle="butt")

    def arr(x, y0, y1):
        ax.annotate("", xy=(x, y1), xytext=(x, y0), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK, mutation_scale=6, shrinkA=0, shrinkB=0))

    def excl(x, y, text, w=0.23):   # side branch from the vertical arrow at x to an exclusion box on its right
        line([x, x + 0.03], [y, y])
        box(x + 0.03 + w / 2, y, w, 0.075, text, fs=5.1)

    lost_c = int(flow["selected p2 cases"] + flow["selected p2b cases"] - cnt["cases_analysed"])
    lost_k = int(flow["selected p2 controls"] + flow["selected p2b controls"] - cnt["controls_analysed"])
    box(0.5, 0.93, 0.8, 0.09, f"NLST CT arm, screening CT in IDC\nn = {int(flow['participants with CT in IDC']):,}")
    line([0.5, 0.5], [0.885, 0.85]); line([0.25, 0.75], [0.85, 0.85])
    arr(0.25, 0.85, 0.8125); arr(0.75, 0.85, 0.8125)
    box(0.25, 0.765, 0.47, 0.095, f"First primary lung cancer\nn = {int(flow['first primary lung cancer, imaged']):,}")
    box(0.75, 0.765, 0.47, 0.095, f"Never lung cancer, 3 rounds imaged\nn = {int(flow['  of whom all three rounds imaged (control pool)']):,}")
    arr(0.25, 0.7175, 0.605); arr(0.75, 0.7175, 0.605)
    excl(0.25, 0.661, f"<2 pre-diagnosis CTs\nn = {int(flow['  with 0-1 imaged round on/before diagnosis (excluded)'])}")
    box(0.25, 0.55, 0.47, 0.11, f"Cases, 2–3 pre-diagnosis CTs\nn = {int(flow['selected p2 cases'] + flow['selected p2b cases'])}\n"
        f"(resected {int(flow['selected p2 cases'])}; not {int(flow['selected p2b cases'])})")
    box(0.75, 0.55, 0.47, 0.11, f"Matched 1:3 (sex, smoking,\nage, scanner vendor)\nn = {int(flow['selected p2 controls'] + flow['selected p2b controls']):,}")
    arr(0.25, 0.495, 0.395); arr(0.75, 0.495, 0.395)
    # participants with <2 valid scans (STROBE 13); scan-level exclusions are in the legend (manuscript/p2/build.py)
    excl(0.25, 0.447, f"<2 valid scans\nn = {lost_c}", w=0.19)
    excl(0.75, 0.447, f"<2 valid scans\nn = {lost_k}", w=0.19)
    box(0.25, 0.345, 0.47, 0.1, f"Analysed: {int(cnt['cases_analysed'])} cases\n{int(cnt['case_scans']):,} scans")
    box(0.75, 0.345, 0.47, 0.1, f"Analysed: {int(cnt['controls_analysed']):,} controls\n{int(cnt['control_scans']):,} scans")
    fig.text(0.005, top(0.985), "a", fontsize=8, fontweight="bold", va="top")
    style.source_data(pd.DataFrame({"step": list(flow.index) + list(cnt.index), "n": list(flow.values) + list(cnt.values)}),
                      "P2_Fig1", "a")
    # b, c: 3D renders of the example scan (case, T0)
    meta = json.load(open(os.path.join(ROOT, "figures", "assets", "p2_render3d_100658_T0.json")))
    for k, (view, letter, title) in enumerate([("spine", "b", "Paraspinal compartment, T6–T11"),
                                                ("fat", "c", "Intramuscular fat")]):
        img = crop_alpha(os.path.join(ROOT, meta["views"][view]["file"]), pad=4)
        ax = fig.add_axes([0.47 + k * 0.265, top(0.53), 0.25, 0.42 * st])
        ax.imshow(img); ax.axis("off")
        ax.set_title(title, fontsize=6.3, pad=2)
        fig.text(0.47 + k * 0.265, top(0.985), letter, fontsize=8, fontweight="bold", va="top")
    # d: timeline
    bx = fig.add_axes([0.05, 0.09 * st, 0.38, 0.26 * st])
    fig.text(0.005, 0.45 * st, "d", fontsize=8, fontweight="bold", va="top")
    dx = 2.6
    for k in range(3):
        bx.plot([k, k], [0.22, 0.8], color=INK, lw=1)
        bx.text(k, 0.86, f"T{k}", ha="center", fontsize=6.3)
    bx.plot([dx, dx], [0.12, 0.9], color=CASE, lw=1.2)
    bx.text(dx, 0.94, "diagnosis", color=CASE, ha="center", fontsize=6.3)
    for a0, a1, t in [(dx - 1, dx, "<1 y"), (dx - 2, dx - 1, "1–2 y"), (-0.4, dx - 2, "≥2 y")]:
        bx.axvspan(a0, a1, ymin=0, ymax=0.13, color=WASH, lw=0)
        bx.text((a0 + a1) / 2, 0.05, t, ha="center", va="center", fontsize=5.5, color=INK2)
    for k0, k1 in [(0, 1), (1, 2)]:
        bx.annotate("", xy=(k1 - 0.04, 0.5), xytext=(k0 + 0.04, 0.5),
                    arrowprops=dict(arrowstyle="-|>", lw=0.8, color=INK2, mutation_scale=6))
        bx.text((k0 + k1) / 2, 0.56, "annual change", ha="center", fontsize=5.3, color=INK2)
    bx.set_xlim(-0.4, 3.2); bx.set_ylim(0, 1); bx.set_xticks([0, 1, 2, 3])
    bx.set_xlabel("Years since randomisation (annual low-dose CT)")
    bx.spines["left"].set_visible(False); bx.tick_params(left=False, labelleft=False)
    # e: serial axial, case and matched control
    met = rd("p2_paraspinal.csv").set_index("case_id")
    fig.text(0.47, 0.45 * st, "e", fontsize=8, fontweight="bold", va="top")
    for i, (who, pid, col) in enumerate([("Case", "100658", CASE), ("Matched control", "113832", CTRL)]):
        for k in range(3):
            cid = f"{pid}_T{k}"
            img, msk, px = axial(cid, "T8", half_mm=70, shift_mm=0)
            w = 0.168
            ax = fig.add_axes([0.47 + k * (w + 0.004), (0.255 - i * 0.205) * st, w, w * 183 / FH * 0.72])
            imshow_ct(ax, img, msk)
            ax.set_facecolor("black")
            ax.text(0.03, 0.95, f"T{k}", transform=ax.transAxes, color="white", fontsize=5.5, va="top")
            ax.text(0.97, 0.05, f"{met.loc[cid, p]:.1f}%", transform=ax.transAxes, color="white", fontsize=5.5,
                    ha="right", va="bottom")
            if k == 0 and i == 1:
                scalebar(ax, px)
        fig.text(0.5, (0.445 - i * 0.205) * st, who, fontsize=6, color=col, va="top")
    return style.save(fig, "P2_Fig1")


# ====================================================================== Figure 4 (function fig2): normative reference
def fig2():
    style.setup()
    cur = rd("extra_centile_curves.csv")
    cur = cur[(cur.age >= 56) & (cur.age <= 73)].copy()
    cur["value"] = cur.groupby(["measure", "sex", "age"]).value.transform(lambda v: np.sort(v.values))  # rearrangement
    d = long_data()
    first = d[d.group == "control"].sort_values("t").groupby("pid").head(1).copy()
    first["vol_cm3"] = np.exp(first.psm_lvol / 100)
    first["age_scan"] = first.age0
    norm = rd("partA_normative_model.csv")
    CM = [("psm_fat_pct", "Intramuscular fat (%)", "fat"), ("psm_hu", "Lean-muscle radiodensity (HU)", "muscle"),
          ("vol_cm3", "Paraspinal volume (cm³)", "muscle"), ("T10_trab_hu", "T10 trabecular bone (HU)", "bone")]
    root = plt.figure(figsize=(W2, 208 * MM))
    fig, sf2 = root.subfigures(2, 1, height_ratios=[156, 52])
    gs = fig.add_gridspec(2, 4, hspace=0.14, wspace=0.42, left=0.075, right=0.975, top=0.925, bottom=0.485)
    for i, sx in enumerate(["male", "female"]):
        for j, (v, ylab, tis) in enumerate(CM):
            col = TIS[tis]
            ax = fig.add_subplot(gs[i, j])
            pts = first[first.sex == sx]
            ax.scatter(pts.age_scan + np.random.default_rng(1).uniform(-0.4, 0.4, len(pts)), pts[v], s=0.7,
                       color="#cfcdc6", lw=0, rasterized=True, zorder=1)
            qd = {tq: cur[(cur.measure == v) & (cur.sex == sx) & (cur.tau == tq)].sort_values("age")
                  for tq in (0.05, 0.25, 0.5, 0.75, 0.95)}
            age = qd[0.5].age.values
            ax.fill_between(age, qd[0.05].value, qd[0.95].value, color=col, alpha=0.13, lw=0, zorder=2)
            ax.fill_between(age, qd[0.25].value, qd[0.75].value, color=col, alpha=0.3, lw=0, zorder=2)
            for tq in (0.05, 0.95):
                ax.plot(age, qd[tq].value, color=col, lw=0.45, alpha=0.7, zorder=3)
            ax.plot(age, qd[0.5].value, color=col, lw=1.4, zorder=4, ls="--" if v == "psm_hu" else "-")
            lo, hi = np.nanpercentile(first[v], [0.5, 99.5])
            ax.set_ylim(lo, hi)
            if j == 0:
                for tq, q in qd.items():
                    if lo < q.value.iloc[-1] < hi:
                        ax.text(73.6, q.value.iloc[-1], f"{int(tq * 100)}", fontsize=4.8, color=INK2, va="center")
            ax.set_xlim(54, 77)
            if i == 0:
                panel_title(ax, ylab, y=1.13)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Age (years)")
            if j == 0:
                ax.set_ylabel("Men" if sx == "male" else "Women", fontsize=7, fontweight="bold")
    lab(fig, fig.axes[0], "a", dx=-0.065, dy=0.05)
    style.source_data(cur, "P2_Fig4", "a")
    # b: longitudinal change terms
    gb = fig.add_gridspec(1, 4, wspace=0.5, left=0.15, right=0.975, top=0.31, bottom=0.075)
    rows = [("t", "Annual change\n(reference)"), ("t:sexfemale", "Women vs men"),
            ("t:smokecurrent", "Current vs former\nsmoker"), ("t:age0c", "Per year older")]
    src = []
    for j, m in enumerate(MEAS):
        ax = fig.add_subplot(gb[0, j])
        sub = norm[norm.measure == m].set_index("term").reindex([r[0] for r in rows]).reset_index()
        y = np.arange(len(rows))[::-1]
        ax.axvline(0, color=INK2, lw=0.5, zorder=1)
        mf, _ = mstyle(m)
        for yk, (_, r) in zip(y, sub.iterrows()):
            ax.plot([r.lo, r.hi], [yk, yk], color=mcol(m), lw=0.9, zorder=2)
            ax.plot(r.estimate, yk, "o", ms=3.8, color=mcol(m), mfc=mf, mec=mcol(m), mew=0.9, zorder=3)
        ax.set_yticks(y); ax.set_yticklabels([r[1] for r in rows] if j == 0 else [], fontsize=5.8)
        ax.tick_params(axis="y", length=0)
        panel_title(ax, LABEL[m])
        ax.set_xlabel(f"{UNIT[m]} per year", color=INK2)
        ax.set_ylim(-0.5, len(rows) - 0.5)
        ax.spines["left"].set_visible(False)
        sub["measure"] = m
        src.append(sub)
    lab(fig, fig.axes[-4], "b", dx=-0.14, dy=0.035)
    style.source_data(pd.concat(src), "P2_Fig4", "b")
    from p2_figs_extra_rows import row_fig4
    row_fig4(sf2)
    return style.save(root, "P2_Fig4")


# ====================================================================== Figure 5 (function fig3): reliability of change
def fig3():
    style.setup()
    p = primary()
    d = long_data()
    rel, norm, xm = rd("partA_reliability.csv"), rd("partA_normative_model.csv"), rd("partA_crossmodel_slopes.csv")
    root = plt.figure(figsize=(W2, 176 * MM))
    fig, sf2 = root.subfigures(2, 1, height_ratios=[124, 52])
    FAT = TIS["fat"]
    # a: change from first scan, controls with 3 scans, with the marginal distribution of the final change
    ax = fig.add_axes([0.07, 0.56, 0.34, 0.37])
    c = d[d.group == "control"].copy()
    gam = norm[(norm.measure == p) & (norm.term == "noise_c")].estimate.iloc[0]
    c["y"] = c[p] - gam * c.noise_c
    c = c.sort_values(["pid", "t"])
    c["dy"] = c.y - c.groupby("pid").y.transform("first")
    three = c.groupby("pid").filter(lambda g: len(g) == 3)
    ids = np.random.default_rng(20261009).choice(three.pid.unique(), 300, replace=False)
    for pid, g in three[three.pid.isin(ids)].groupby("pid"):
        ax.plot(g.t, g.dy, color=FAT, lw=0.35, alpha=0.22, zorder=1)
    slope = norm[(norm.measure == p) & (norm.term == "t")].estimate.iloc[0]
    mdc = rel[(rel.measure == p) & (rel.subset == "all controls")].mdc95_2y.iloc[0]
    tt = np.linspace(0, 2.2, 10)
    ax.fill_between(tt, -mdc, mdc, color=FAT, alpha=0.07, zorder=0, lw=0)
    for s_ in (-mdc, mdc):
        ax.plot(tt, [s_] * len(tt), color=FAT, lw=0.5, ls=(0, (3, 2)), alpha=0.8)
    ax.plot(tt, slope * tt, color=INK, lw=1.4, zorder=3)
    ax.text(2.25, slope * 2.2, "normal\nageing", fontsize=5.5, va="center")
    ax.text(0.05, mdc * 1.05, "2-year MDC₉₅", fontsize=5.4, color=FAT, va="bottom")
    ax.set_xlim(-0.05, 2.75); ax.set_ylim(-mdc * 1.6, mdc * 1.6)
    ax.set_xlabel("Years since first screening CT"); ax.set_ylabel("Change in intramuscular fat (pp)")
    lab(fig, ax, "a", dx=-0.06, dy=0.01)
    axm = fig.add_axes([0.415, 0.56, 0.045, 0.37], sharey=ax)
    last = three.groupby("pid").tail(1).dy.values
    axm.hist(last[np.abs(last) < mdc * 1.6], bins=np.linspace(-mdc * 1.6, mdc * 1.6, 31), orientation="horizontal",
             color=FAT, alpha=0.55, lw=0)
    axm.axhline(slope * np.median(three.groupby("pid").t.max()), color=INK, lw=0.8)
    for s_ in (-mdc, mdc):
        axm.axhline(s_, color=FAT, lw=0.5, ls=(0, (3, 2)))
    axm.axis("off")
    axm.text(0.5, 1.0, "final\nchange", transform=axm.transAxes, fontsize=4.8, ha="center", va="bottom", color=INK2)
    style.source_data(three[three.pid.isin(ids)][["pid", "t", "dy"]], "P2_Fig5", "a")
    # b: reliability of a single scan versus of change (dumbbell)
    ax = fig.add_axes([0.57, 0.6, 0.41, 0.33])
    r = rel[rel.subset == "all controls"].set_index("measure").reindex(ALL).reset_index()
    x = np.arange(len(ALL))
    ax.axhspan(0, 0.4, color="#f7f6f2", lw=0, zorder=0)
    ax.axhline(0.40, color=INK2, lw=0.6, ls=(0, (3, 2)))
    for xi, (_, rr) in zip(x, r.iterrows()):
        col = mcol(rr.measure)
        mf, _ = mstyle(rr.measure)
        ax.plot([xi, xi], [max(rr.rel_hi, rr.slope_reliability) + 0.02, rr.level_icc - 0.025], color="#c9c7bf",
                lw=0.6, ls=(0, (1, 1.5)), zorder=1)
        ax.plot(xi, rr.level_icc, "o", ms=4.2, mfc="white", mec=INK2, mew=0.8, zorder=3)
        ax.plot([xi, xi], [rr.rel_lo, rr.rel_hi], color=col, lw=1.0, zorder=3)
        ax.plot(xi, rr.slope_reliability, "o", ms=4.6, color=col, mfc=mf, mec=col, mew=0.9, zorder=4)
        ax.text(xi + 0.12, rr.slope_reliability, f"{rr.slope_reliability:.2f}", fontsize=5.2, va="center", color=INK2)
    ax.set_xticks(x); ax.set_xticklabels([SHORT[m] for m in ALL], fontsize=5.5)
    ax.tick_params(axis="x", length=0, pad=1)
    ax.set_ylim(0, 1.05); ax.set_xlim(-0.5, len(ALL) - 0.3); ax.set_ylabel("Reliability")
    ax.plot([], [], "o", ms=4, mfc="white", mec=INK2, label="single scan (ICC)")
    ax.plot([], [], "o", ms=4, color=INK2, label="individual change (slope reliability)")
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=5.3)
    lab(fig, ax, "b", dx=-0.06, dy=0.035)
    style.source_data(r, "P2_Fig5", "b")
    # c: years of ageing to exceed MDC
    ax = fig.add_axes([0.13, 0.08, 0.19, 0.33])
    ref = norm[norm.term == "t"].set_index("measure").reindex(MEAS)
    yrs = (r.set_index("measure").mdc95_2y.reindex(MEAS) / ref.estimate.abs())
    cap = 40
    yy = np.arange(len(MEAS))[::-1]
    SH = ["Intramuscular fat", "Radiodensity", "Volume", "Trabecular bone"]
    for yi, v, m in zip(yy, yrs.values, MEAS):
        ax.barh(yi, min(v, cap), height=0.56, color=mcol(m), alpha=0.22, edgecolor=mcol(m), lw=0.8,
                hatch="////" if m in HOLLOW else None)
        ax.text(min(v, cap) + 0.8, yi, f">{cap}" if v > cap else f"{v:.0f} y", va="center", fontsize=5.4, color=INK2)
    ax.set_yticks(yy); ax.set_yticklabels(SH, fontsize=5.8)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.set_xlim(0, cap * 1.25); ax.set_xlabel("Years of normal ageing to\nexceed the 2-year MDC₉₅")
    ax.spines["left"].set_visible(False)
    lab(fig, ax, "c", dx=-0.12, dy=0.02)
    style.source_data(pd.DataFrame({"measure": MEAS, "years": yrs.values}), "P2_Fig5", "c")
    # d: kernel switch (dumbbell)
    ax = fig.add_axes([0.47, 0.08, 0.2, 0.33])
    kc = rel[rel.subset == "kernel consistent"].set_index("measure").slope_reliability.reindex(MEAS)
    ks = rel[rel.subset == "kernel switched"].set_index("measure").slope_reliability.reindex(MEAS)
    nks = int(rel[rel.subset == "kernel switched"].n.iloc[0])
    for yi, a, b, m in zip(yy, kc.values, ks.values, MEAS):
        ax.plot([b, a], [yi, yi], color=mcol(m), lw=1.6, alpha=0.35, zorder=1, solid_capstyle="round")
        ax.plot(a, yi, "o", ms=4.4, color=mcol(m), zorder=3)
        ax.plot(b, yi, "o", ms=4.4, mfc="white", mec=mcol(m), mew=1.0, zorder=3)
    ax.plot([], [], "o", ms=4, color=INK2, label="same kernel in all rounds")
    ax.plot([], [], "o", ms=4, mfc="white", mec=INK2, label=f"kernel switched (n = {nks})")
    ax.set_yticks(yy); ax.set_yticklabels(SH, fontsize=5.8)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, len(MEAS) - 0.5)
    ax.set_xlim(-0.02, 0.6); ax.set_xlabel("Slope reliability")
    ax.spines["left"].set_visible(False)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.05, 1.0), fontsize=5.2)
    lab(fig, ax, "d", dx=-0.12, dy=0.08)
    style.source_data(pd.DataFrame({"measure": MEAS, "kernel_consistent": kc.values, "kernel_switched": ks.values}),
                      "P2_Fig5", "d")
    # e: two segmentation models, density
    ax = fig.add_axes([0.765, 0.08, 0.21, 0.33])
    a_, b_ = xm[p], xm[p.replace("psm", "aut")]
    ok = a_.notna() & b_.notna()
    lim = np.nanpercentile(np.r_[a_[ok], b_[ok]], [1, 99])
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("fat", ["#efecf7", "#b9addf", FAT, "#2c2170"])
    hb = ax.hexbin(a_[ok], b_[ok], gridsize=38, extent=(*lim, *lim), cmap=cmap, mincnt=1, bins="log", lw=0.1,
                   edgecolors="white")
    ax.plot(lim, lim, color=INK2, lw=0.6, ls=(0, (3, 2)))
    rho = pd.Series(a_[ok]).corr(pd.Series(b_[ok]), method="spearman")
    ax.text(0.04, 0.96, f"Spearman ρ = {rho:.2f}\nn = {ok.sum():,}", transform=ax.transAxes, fontsize=5.5, va="top")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Slope, model 1 (pp/y)"); ax.set_ylabel("Slope, model 2 (pp/y)")
    cax = ax.inset_axes([0.86, 0.06, 0.035, 0.32])
    cb = fig.colorbar(hb, cax=cax); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=4.6, length=1.5)
    cb.set_ticks([1, 10, 100]); cb.set_ticklabels(["1", "10", "100"]); cb.ax.minorticks_off()
    cax.set_title("n", fontsize=4.8, pad=1.5)
    cax.yaxis.set_ticks_position("left")
    lab(fig, ax, "e", dx=-0.06, dy=0.02)
    style.source_data(pd.DataFrame({"model1": a_, "model2": b_}), "P2_Fig5", "e")
    from p2_figs_extra_rows import row_fig5
    row_fig5(sf2)
    return style.save(root, "P2_Fig5")


# ====================================================================== Figure 6 (function fig4): primary hypothesis
def fig4():
    style.setup()
    p = primary()
    div, ivs, als, maa = rd("partA_divergence.csv"), rd("partA_interval_windows.csv"), rd("partA_aligned_bins.csv"), \
        rd("p2_maa.csv").astype({"pid": str})
    long = long_data()
    grp = long[["pid", "group"]].drop_duplicates()
    root = plt.figure(figsize=(W2, 184 * MM))
    fig, sf2 = root.subfigures(2, 1, height_ratios=[130, 54])
    gs = fig.add_gridspec(1, 4, wspace=0.12, left=0.08, right=0.985, top=0.88, bottom=0.6)
    ax0 = None
    src = []
    for j, m in enumerate(MEAS):
        ax = fig.add_subplot(gs[0, j], sharey=ax0); ax0 = ax0 or ax
        w = div[(div.measure == m) & (div.analysis == "interval z: case window vs controls")].copy()
        w["window"] = w.term.str.replace("^window", "", regex=True)
        w = w.set_index("window").reindex(WIN).reset_index()
        n = ivs[ivs.measure == m].set_index("window").reindex(WIN).n.values
        ax.axhspan(-0.2, 0.2, color="#f1f0eb", lw=0, zorder=0)
        ax.axhline(0, color=INK2, lw=0.6)
        mf, ls = mstyle(m)
        xs = np.arange(3)
        ax.plot(xs, w.estimate.values, color=mcol(m), lw=0.7, ls=ls, alpha=0.6, zorder=2)
        for xi, (_, r) in zip(xs, w.iterrows()):
            ax.plot([xi, xi], [r.lo, r.hi], color=mcol(m), lw=1.1, zorder=3, solid_capstyle="butt")
            ax.plot(xi, r.estimate, "o", ms=4.6, color=mcol(m), mfc=mf, mec=mcol(m), mew=1.0, zorder=4)
        ax.set_xticks(xs); ax.set_xticklabels([f"{a}\n({b})" for a, b in zip(WIN_LAB, n)], fontsize=5.5)
        ax.set_xlim(-0.5, 2.5); ax.set_ylim(-0.3, 0.3)
        panel_title(ax, LABEL[m] + (" (primary)" if m == p else ""), y=1.12)
        if j:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("Excess ageing in cases (z)")
        w["measure"] = m
        src.append(w)
    fig.text(0.53, 0.52, "Years from the end of the interval to diagnosis (number of case intervals)", ha="center",
             fontsize=6, color=INK2)
    lab(fig, fig.axes[0], "a", dx=-0.065, dy=0.075)
    style.source_data(pd.concat(src), "P2_Fig6", "a")
    # b: raincloud of acceleration, cases vs controls
    ax = fig.add_axes([0.105, 0.075, 0.24, 0.33])
    mm = maa[maa.measure == p].merge(grp, on="pid")
    zc0, zk0 = mm[mm.group == "control"].maa_z.dropna().values, mm[mm.group == "case"].maa_z.dropna().values
    zc, zk = zc0.clip(-4, 4), zk0.clip(-4, 4)
    raincloud(ax, [zc, zk], [CTRL, CASE], ["Controls", "Cases"], [1.0, 0.0], bw=0.2, cloud_h=0.6, xlim=(-4, 4),
              point_alpha=0.2)
    ax.set_yticks([0.12, 1.12])
    ax.set_yticklabels([f"Cases\n(n = {len(zk)})", f"Controls\n(n = {len(zc):,})"], fontsize=5.6)
    for t, c in zip(ax.get_yticklabels(), [CASE, CTRL]):
        t.set_color(c)
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.42, 1.95); ax.set_xlim(-4, 4)
    ax.axvline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=0)
    ax.set_xlabel("Muscle-ageing acceleration, intramuscular fat (z)")
    lab(fig, ax, "b", dx=-0.09, dy=0.03)
    style.source_data(mm[["pid", "group", "maa_z"]], "P2_Fig6", "b")
    # c: same-round difference over time, primary measure
    ax = fig.add_axes([0.45, 0.075, 0.22, 0.33])
    b = als[als.measure == p]
    ax.axhline(0, color=INK2, lw=0.6)
    ribbon(ax, b.mean_tau.values, b.mean_diff.values, b.lo.values, b.hi.values, mcol(p), ms=3.4)
    ax.set_xlabel("Years relative to diagnosis"); ax.set_ylabel("Case minus matched controls (pp)")
    lab(fig, ax, "c", dx=-0.07, dy=0.03)
    style.source_data(b, "P2_Fig6", "c")
    # d: conditional logistic OR per SD of acceleration
    ax = fig.add_axes([0.83, 0.075, 0.07, 0.3])
    cl = div[div.analysis == "conditional logistic"].set_index("measure").reindex(ALL)
    rows = [dict(kind="row", label=LABEL[m].replace(" (15th pct)", "").replace(" (T10)", ""), est=cl.loc[m, "estimate"],
                 lo=cl.loc[m, "lo"], hi=cl.loc[m, "hi"], color=mcol(m), hollow=m in HOLLOW)
            for m in ALL]
    forest(ax, rows, xlim=(0.85, 1.18), ticks=[0.9, 1, 1.1], label_w=2.15, right_cols=(("OR (95% CI)", 1.08),),
           xlabel="OR per SD of\nacceleration", fs=5.4)
    lab(fig, ax, "d", dx=-0.16, dy=0.06)
    style.source_data(cl.reset_index(), "P2_Fig6", "d")
    from p2_figs_extra_rows import row_fig6
    row_fig6(sf2)
    return style.save(root, "P2_Fig6")




# ====================================================================== Figure 5: subgroups and tumour burden
def fig5():
    style.setup()
    p = primary()
    sg, st, maa = rd("extra_subgroups.csv"), rd("partA_case_strata.csv"), rd("p2_maa.csv").astype({"pid": str})
    fig = plt.figure(figsize=(W2, 125 * MM))
    # a: subgroup forest
    ax = fig.add_axes([0.25, 0.07, 0.28, 0.86])
    order = ["Age", "Sex", "Smoking", "Stage", "Histology", "Detection", "Resection", "Scanner vendor"]
    rows, y, yt, ylab = [], 0, [], []
    for f in order:
        q = sg[(sg.factor == f) & (sg.window == "all")]
        ax.text(-0.74, -y, f"{f}  (P interaction {q.p_interaction.iloc[0]:.2f})", fontsize=5.8, fontweight="bold",
                va="center", transform=ax.get_yaxis_transform())
        y += 1
        for _, r in q.iterrows():
            q1 = sg[(sg.factor == f) & (sg.window == "<1 y") & (sg.subgroup == r.subgroup)].iloc[0]
            rows.append((y, r, q1)); yt.append(-y); ylab.append(f"{r.subgroup} ({r.n_cases})")
            y += 1
        y += 0.4
    ax.axvspan(-0.2, 0.2, color=BAND, lw=0, zorder=0)
    ax.axvline(0, color=INK2, lw=0.6)
    dot_ci(ax, [r[1].estimate for r in rows], [-r[0] + 0.15 for r in rows], [r[1].lo for r in rows],
           [r[1].hi for r in rows], mcol(p), ms=3, label="all intervals")
    dot_ci(ax, [r[2].estimate for r in rows], [-r[0] - 0.15 for r in rows], [r[2].lo for r in rows],
           [r[2].hi for r in rows], INK2, ms=2.6, mfc="white", marker="s", label="interval ending <1 y before diagnosis")
    ax.set_yticks(yt); ax.set_yticklabels(ylab, fontsize=5.5)
    ax.set_ylim(-y + 0.5, 0.8); ax.set_xlim(-0.75, 0.75)
    ax.set_xlabel("Excess ageing in cases (z), intramuscular fat")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=5.3, ncol=1)
    fig.text(0.005, 0.985, "a", fontsize=8, fontweight="bold", va="top")
    style.source_data(sg, "P2_EDFig10", "a")
    # b: MAA by stage, histology, detection
    ax = fig.add_axes([0.71, 0.5, 0.27, 0.43])
    s = st[st.measure == p]
    spec = [("stage4", ["I", "II", "III", "IV"]), ("histo", ["adenocarcinoma", "squamous", "small cell", "other"]),
            ("detection", ["screen-detected", "interval/post-screening"])]
    rr, yl, y = [], [], 0
    for fac, levels in spec:
        for lv in levels:
            q = s[(s.factor == fac) & (s.stratum == lv)]
            if len(q):
                rr.append((y, q.mean_maa.iloc[0], q.se.iloc[0]))
                yl.append(f"{'Stage ' + lv if fac == 'stage4' else lv[0].upper() + lv[1:]} ({int(q.n.iloc[0])})")
                y += 1
        y += 0.6
    yy = -np.array([a[0] for a in rr])
    ax.axvline(0, color=INK2, lw=0.6)
    dot_ci(ax, np.array([a[1] for a in rr]), yy, np.array([a[1] - 1.96 * a[2] for a in rr]),
           np.array([a[1] + 1.96 * a[2] for a in rr]), mcol(p))
    ax.set_yticks(yy); ax.set_yticklabels(yl, fontsize=5.5)
    ax.set_xlabel("Acceleration (z, mean ± 95% CI)")
    fig.text(0.56, 0.985, "b", fontsize=8, fontweight="bold", va="top")
    style.source_data(s, "P2_EDFig10", "b")
    # c: MAA vs tumour size
    ax = fig.add_axes([0.71, 0.07, 0.27, 0.32])
    cas = sel("cases")
    x = maa[maa.measure == p].merge(cas[["pid", "lesionsize"]], on="pid")
    x["size"] = pd.to_numeric(x.lesionsize, errors="coerce")
    x = x[(x["size"] > 0) & (x["size"] <= 120) & x.maa_z.notna()]
    ax.scatter(x["size"], x.maa_z.clip(-4, 4), s=3, color=mcol(p), alpha=0.4, lw=0, rasterized=True)
    b1, b0 = np.polyfit(x["size"], x.maa_z, 1)
    xs = np.linspace(x["size"].min(), x["size"].max(), 50)
    boot = np.array([np.polyval(np.polyfit(*x.sample(len(x), replace=True, random_state=i)[["size", "maa_z"]].values.T, 1), xs)
                     for i in range(300)])
    ax.fill_between(xs, *np.percentile(boot, [2.5, 97.5], axis=0), color=mcol(p), alpha=0.15, lw=0)
    ax.plot(xs, b0 + b1 * xs, color=mcol(p), lw=1)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlabel("Tumour size (mm)"); ax.set_ylabel("Acceleration (z)")
    ax.set_ylim(-4.2, 4.2)
    fig.text(0.56, 0.43, "c", fontsize=8, fontweight="bold", va="top")
    style.source_data(x[["pid", "size", "maa_z"]], "P2_EDFig10", "c")
    return style.save(fig, "P2_EDFig10")


# ====================================================================== Figure 3 (function fig6): stable level offsets
def fig6():
    style.setup()
    sd = control_sd()
    div, als, cen, off = rd("partA_divergence.csv"), rd("partA_aligned_bins.csv"), rd("extra_case_centiles.csv"), \
        rd("extra_offsets_by_subgroup.csv")
    root = plt.figure(figsize=(W2, 186 * MM))
    fig, sf2 = root.subfigures(2, 1, height_ratios=[128, 58])
    # a: level offsets in control SD, with and without adjustment for lung density
    ax = fig.add_axes([0.235, 0.57, 0.17, 0.36])
    order = ["lung_perc15_hu", "T10_trab_hu", "psm_lvol", "psm_hu", "psm_fat_pct"]
    rows = []
    for m in order:
        for an, kind, labl in [("mixed model", "row", LABEL[m].replace(" (15th pct)", "").replace(" (T10)", "")),
                               ("mixed model, adjusted for lung density", "sub", "+ lung density")]:
            q = div[(div.term == "groupcase") & (div.analysis == an) & (div.measure == m)]
            if not len(q):
                continue
            q = q.iloc[0]
            e, lo, hi = q.estimate / sd[m], q.lo / sd[m], q.hi / sd[m]
            rows.append(dict(kind=kind, label=labl, est=e, lo=lo, hi=hi, color=mcol(m),
                             hollow=(kind == "sub") or m in HOLLOW,
                             txt=f"{fnum(e)} ({fnum(lo)} to {fnum(hi)})"))
    forest(ax, rows, xlim=(-0.32, 0.17), ticks=[-0.3, -0.2, -0.1, 0, 0.1], label_w=1.3, ref=0, log=False,
           right_cols=(("SD (95% CI)", 1.04),), xlabel="Case minus matched controls (control SD)", fs=5.5)
    ax.set_xticklabels([fnum(t, "{:g}") for t in [-0.3, -0.2, -0.1, 0, 0.1]])
    lab(fig, ax, "a", dx=-0.23, dy=0.04)
    srcA = div[(div.term == "groupcase") & div.analysis.isin(["mixed model", "mixed model, adjusted for lung density"])].copy()
    srcA["control_sd"] = srcA.measure.map(sd)
    for c_ in ["estimate", "lo", "hi"]:
        srcA[c_ + "_sd"] = srcA[c_] / srcA.control_sd
    style.source_data(srcA, "P2_Fig3", "a")
    # b: same-round differences over time, ribbons
    ax = fig.add_axes([0.66, 0.57, 0.24, 0.36])
    offx = {"psm_lvol": -0.09, "T10_trab_hu": -0.03, "lung_perc15_hu": 0.03, "psm_fat_pct": 0.09}
    for m in ["psm_fat_pct", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]:
        b = als[als.measure == m]
        ribbon(ax, b.mean_tau.values + offx[m], (b.mean_diff / sd[m]).values, (b.lo / sd[m]).values,
               (b.hi / sd[m]).values, mcol(m), ms=2.8, alpha=0.06)
        yl = (b.mean_diff / sd[m]).iloc[-1] + {"psm_lvol": -0.006, "T10_trab_hu": -0.02, "lung_perc15_hu": -0.03,
                                                "psm_fat_pct": 0.0}[m]
        ax.text(b.mean_tau.iloc[-1] + 0.25, yl, LABEL[m].replace(" (15th pct)", "").replace(" (T10)", ""), fontsize=5.4,
                color=mcol(m), va="center")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlim(-5.6, 0.2)
    ax.set_xlabel("Years relative to diagnosis"); ax.set_ylabel("Case minus controls (control SD)")
    lab(fig, ax, "b", dx=-0.07, dy=0.04)
    style.source_data(als, "P2_Fig3", "b")
    # c: case centiles on the control reference
    gc = fig.add_gridspec(1, 3, wspace=0.38, left=0.075, right=0.5, top=0.4, bottom=0.085)
    for j, (m, t, tis) in enumerate([("vol_cm3", "Muscle volume", "muscle"), ("T10_trab_hu", "Trabecular bone", "bone"),
                                     ("psm_fat_pct", "Intramuscular fat", "fat")]):
        ax = fig.add_subplot(gc[0, j])
        q = cen[(cen.measure == m) & (cen.group == "case")].centile
        h, e = np.histogram(q, bins=np.linspace(0, 1, 11))
        pct = 100 * h / h.sum()
        cols = [TIS[tis] if k < 2 else "#d9d7d0" for k in range(10)]
        ax.bar(e[:-1] * 100 + 5, pct, width=8.6, color=cols, lw=0, alpha=0.9)
        ax.axhline(10, color=INK2, lw=0.6, ls=(0, (3, 2)))
        ax.set_ylim(0, 16); ax.set_xticks([0, 50, 100])
        panel_title(ax, t, y=1.12, fs=5.9)
        ax.set_xlabel("Centile vs controls")
        if j == 0:
            ax.set_ylabel("Cases (%)")
    lab(fig, fig.axes[-3], "c", dx=-0.065, dy=0.04)
    style.source_data(cen[cen.group == "case"], "P2_Fig3", "c")
    # d: offsets by stage
    ax = fig.add_axes([0.62, 0.085, 0.36, 0.3])
    stg = ["I", "II", "III", "IV"]
    for k, m in enumerate(["psm_lvol", "T10_trab_hu", "lung_perc15_hu"]):
        q = off[(off.factor == "stage") & (off.measure == m)].set_index("subgroup").reindex(stg)
        xx = np.arange(4) + (k - 1) * 0.2
        for xi, (_, r) in zip(xx, q.iterrows()):
            ax.plot([xi, xi], [r.lo / sd[m], r.hi / sd[m]], color=mcol(m), lw=1.0, zorder=2)
        ax.plot(xx, (q.estimate / sd[m]).values, "o", ms=4.0, color=mcol(m), zorder=3,
                label=LABEL[m].replace(" (15th pct)", "").replace(" (T10)", ""))
    ax.axhline(0, color=INK2, lw=0.6)
    ns = off[(off.factor == "stage") & (off.measure == "psm_lvol")].set_index("subgroup").reindex(stg).n_cases
    ax.set_xticks(np.arange(4)); ax.set_xticklabels([f"{s}\n({int(n)})" for s, n in zip(stg, ns)])
    ax.set_xlim(-0.5, 3.5)
    ax.set_xlabel("Stage (cases)"); ax.set_ylabel("Level offset (control SD)")
    ax.legend(loc="lower left", fontsize=5.3, ncol=3, bbox_to_anchor=(0, 1.0), handletextpad=0.3, columnspacing=1.0)
    lab(fig, ax, "d", dx=-0.07, dy=0.09)
    style.source_data(off[off.factor == "stage"], "P2_Fig3", "d")
    from p2_figs_extra_rows import row_fig3
    row_fig3(sf2)
    return style.save(root, "P2_Fig3")


# ====================================================================== ED4: observed slopes vs true-slope spread
def ed4():
    style.setup()
    rel = rd("partA_reliability.csv")
    xm = rd("partA_crossmodel_slopes.csv").astype({"pid": str})
    long = long_data()
    grp = long[["pid", "group"]].drop_duplicates()
    fig, axes = plt.subplots(1, 3, figsize=(W2, 58 * MM), gridspec_kw=dict(wspace=0.35, left=0.07, right=0.98,
                                                                          top=0.85, bottom=0.22))
    for ax, m in zip(axes, ["psm_fat_pct", "psm_hu", "psm_lvol"]):
        s = xm.merge(grp, on="pid")
        z = s[s.group == "control"][m].dropna()
        r = rel[(rel.measure == m) & (rel.subset == "all controls")].iloc[0]
        lo, hi = np.percentile(z, [0.5, 99.5])
        ax.hist(z.clip(lo, hi), bins=60, density=True, color=WASH, edgecolor=INK2, lw=0.3)
        xx = np.linspace(lo, hi, 300)
        mu = z.mean()
        ax.plot(xx, np.exp(-0.5 * ((xx - mu) / r.slope_true_sd) ** 2) / (r.slope_true_sd * np.sqrt(2 * np.pi)), color=CTRL,
                lw=1, label="true-slope spread")
        ax.set_title(LABEL[m], fontsize=6.5, pad=3)
        ax.set_xlabel(f"Individual slope ({UNIT[m]}/year)")
    axes[0].legend(loc="upper right", fontsize=5.2, handlelength=1.0)
    axes[0].set_ylabel("Density (controls)")
    lab(fig, axes[0], "a", dx=-0.06, dy=0.03)
    style.source_data(xm, "P2_EDFig4", "a")
    return style.save(fig, "P2_EDFig4")


# ====================================================================== ED5: sensitivity, all measures
def ed5():
    style.setup()
    runs = [("Main analysis", ""), ("Same kernel all rounds", "sens_kernel_consistent"), ("Three rounds only", "sens_three_rounds"),
            ("Excluding scans <90 d before dx", "sens_gap90"), ("Second segmentation model", "sens_autochthon"),
            ("Resected cancers only", "sens_resected_only"), ("Non-resected cancers only", "sens_unresected_only")]
    terms = [("windowcase: ends >=2 y before dx", "≥2 y"), ("windowcase: ends 1-2 y before dx", "1–2 y"),
             ("windowcase: ends <1 y before dx", "<1 y")]
    rows = []
    for name, sub in runs:
        d = rd("partA_divergence.csv", os.path.join(RES, sub) if sub else RES)
        for m in MEAS:
            for t, tl in terms:
                q = d[(d.measure == m) & (d.term == t)]
                if len(q):
                    rows.append(dict(run=name, measure=m, window=tl, estimate=q.estimate.iloc[0], lo=q.lo.iloc[0], hi=q.hi.iloc[0]))
    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(len(MEAS), 3, figsize=(W2, 170 * MM), sharex=True, sharey=True,
                             gridspec_kw=dict(hspace=0.38, wspace=0.08, left=0.22, right=0.985, top=0.93, bottom=0.06))
    names = [r[0] for r in runs]
    for i, m in enumerate(MEAS):
        for j, (_, tl) in enumerate(terms):
            ax = axes[i, j]
            q = df[(df.measure == m) & (df.window == tl)].set_index("run").reindex(names)
            yy = -np.arange(len(names))
            ax.axvspan(-0.2, 0.2, color=BAND, lw=0, zorder=0)
            ax.axvline(0, color=INK2, lw=0.6)
            dot_ci(ax, q.estimate.values, yy, q.lo.values, q.hi.values, mcol(m), ms=2.8,
                   mfc="white" if m in HOLLOW else None)
            ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=5.5)
            if i == 0:
                ax.set_title(f"Interval ends {tl} before diagnosis", fontsize=6.3, pad=12)
            if j == 0:
                ax.text(-0.78, 1.06, LABEL[m], transform=ax.transAxes, fontsize=7, fontweight="bold", va="bottom",
                        ha="left")
            if i == len(MEAS) - 1:
                ax.set_xlabel("Excess ageing in cases (z)")
    ax.set_xlim(-0.5, 0.5)
    style.source_data(df, "P2_EDFig5", "a")
    return style.save(fig, "P2_EDFig5")


# ====================================================================== ED7: level offsets by histology and detection
def ed7():
    style.setup()
    sd = control_sd()
    off = rd("extra_offsets_by_subgroup.csv")
    fig, axes = plt.subplots(1, len(ALL), figsize=(W2, 72 * MM), sharey=True,
                             gridspec_kw=dict(wspace=0.12, left=0.2, right=0.985, top=0.86, bottom=0.18))
    spec = [("histo", ["Adenocarcinoma", "Squamous", "Small cell", "Other"]),
            ("detection", ["Screen-detected", "Interval/post-screening"]), ("stage", ["I", "II", "III", "IV"])]
    labels, ys, y = [], [], 0
    for f, lv in spec:
        for l in lv:
            labels.append((f, l)); ys.append(-y); y += 1
        y += 0.6
    for ax, m in zip(axes, ALL):
        for (f, l), yy in zip(labels, ys):
            q = off[(off.factor == f) & (off.subgroup == l) & (off.measure == m)]
            if len(q):
                dot_ci(ax, [q.estimate.iloc[0] / sd[m]], [yy], [q.lo.iloc[0] / sd[m]], [q.hi.iloc[0] / sd[m]], mcol(m),
                       ms=2.8, mfc="white" if m in HOLLOW else None)
        ax.axvline(0, color=INK2, lw=0.6)
        ax.set_title(LABEL[m].replace(" (15th pct)", "\n(15th pct)").replace(" (T10)", "\n(T10)"), fontsize=6, pad=3)
        ax.set_xlabel("Offset (SD)")
    nmap = off[off.measure == "psm_lvol"].set_index(["factor", "subgroup"]).n_cases
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels([f"{'Stage ' + l if f == 'stage' else l} ({int(nmap.get((f, l), 0))})" for f, l in labels],
                            fontsize=5.5)
    style.source_data(off, "P2_EDFig7", "a")
    return style.save(fig, "P2_EDFig7")


# ====================================================================== ED8: lung density positive control
def ed8():
    style.setup()
    div, norm, ivs = rd("partA_divergence.csv"), rd("partA_normative_model.csv"), rd("partA_interval_windows.csv")
    long = long_data()
    fig = plt.figure(figsize=(W2, 60 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 0.8, 0.9], wspace=0.75, left=0.08, right=0.985, top=0.85, bottom=0.25)
    ax = fig.add_subplot(gs[0, 0])
    c = long[(long.group == "control")].sort_values("t").groupby("pid").head(1)
    for g, col in [("former", CTRL), ("current", CASE)]:
        v = c[c.smoke == g].lung_perc15_hu
        ax.hist(v, bins=np.arange(-1000, -700, 8), density=True, histtype="step", color=col, lw=1, label=f"{g} smokers")
    ax.set_xlabel("Lung density, 15th percentile (HU)"); ax.set_ylabel("Density (controls)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.3)
    lab(fig, ax, "a", dx=-0.06, dy=0.03)
    ax = fig.add_subplot(gs[0, 1])
    w = div[(div.measure == "lung_perc15_hu") & (div.analysis == "interval z: case window vs controls")].copy()
    w["window"] = w.term.str.replace("^window", "", regex=True)
    w = w.set_index("window").reindex(WIN).reset_index()
    ax.axhline(0, color=INK2, lw=0.6)
    dot_ci(ax, np.arange(3), w.estimate.values, w.lo.values, w.hi.values, mcol("lung_perc15_hu"), horizontal=False)
    ax.set_xticks(np.arange(3)); ax.set_xticklabels(WIN_LAB)
    ax.set_xlabel("Years before diagnosis"); ax.set_ylabel("Excess change in cases (z)")
    lab(fig, ax, "b", dx=-0.06, dy=0.03)
    ax = fig.add_subplot(gs[0, 2])
    q = div[(div.measure == "lung_perc15_hu") & (div.term == "groupcase") & (div.analysis == "mixed model")].iloc[0]
    t = norm[(norm.measure == "lung_perc15_hu") & (norm.term == "t:smokecurrent")].iloc[0]
    items = [("Cases vs controls\n(level)", q), ("Current vs former\nsmokers (change/y)", t)]
    for k, (name, r) in enumerate(items):
        dot_ci(ax, [r.estimate], [-k], [r.lo], [r.hi], mcol("lung_perc15_hu") if k == 0 else INK2)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_yticks([0, -1]); ax.set_yticklabels([i[0] for i in items], fontsize=5.5)
    ax.set_ylim(-1.6, 0.6); ax.set_xlabel("HU (95% CI)")
    lab(fig, ax, "c", dx=-0.12, dy=0.03)
    style.source_data(w, "P2_EDFig8", "b")
    return style.save(fig, "P2_EDFig8")


if __name__ == "__main__":
    todo = sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "ed4", "ed5", "ed7", "ed8"]
    for t in todo:
        print(globals()[t]())
