"""Additional panel rows for paper-2 Figures 2-6, drawn into a matplotlib SubFigure under the existing panels.
Data: results/p2/extra_*.csv (stats/p2_extra_panels.py) and existing analysis outputs; descriptive summaries only.
"""
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle

from p2_cns_common import CASE, CTRL, INK, INK2, MUTED, LABEL, axial, fnum, imshow_ct, rd, scalebar, style
from p2_premium import HOLLOW, TIS, forest, mcol, mstyle, raincloud, ribbon

SHORTL = {"lung_perc15_hu": "Lung density", "T10_trab_hu": "Trabecular bone", "psm_lvol": "Muscle volume",
          "psm_hu": "Muscle radiodensity", "psm_fat_pct": "Intramuscular fat"}


def letter(sf, x, y, s):
    sf.text(x, y, s, fontsize=8, fontweight="bold", va="top", ha="left", color=INK)


# ============================================================ Figure 2: g, h, i
def row_fig2(sf):
    from p2_figs_pheno import lung_slice, t10_crop
    # g: reserve score at every scan, by time to diagnosis
    letter(sf, 0.004, 0.99, "g")
    t = rd("extra_score_by_time.csv")
    ax = sf.add_axes([0.075, 0.2, 0.25, 0.66])
    x = t.mean_tau.values
    for grp, col in [("control", CTRL), ("case", CASE)]:
        ribbon(ax, x, t[f"{grp}_mean"].values, t[f"{grp}_lo"].values, t[f"{grp}_hi"].values, col, ms=3.0, alpha=0.15)
    ax.axhline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlabel("Years relative to the case’s diagnosis")
    ax.set_ylabel("Reserve score (control SD)")
    ax.set_ylim(-0.6, 0.25)
    ax.text(0.02, 0.97, "Controls", color=CTRL, transform=ax.transAxes, fontsize=5.6, va="top")
    ax.text(0.02, 0.1, "Cases", color=CASE, transform=ax.transAxes, fontsize=5.6, va="bottom")
    style.source_data(t, "P2_Fig2", "g")
    # h: correlations between the components (controls)
    letter(sf, 0.39, 0.99, "h")
    cm = rd("extra_component_corr.csv").set_index("Unnamed: 0")
    order = ["z_lung_perc15_hu", "z_T10_trab_hu", "z_psm_lvol", "z_psm_hu", "z_psm_fat_pct"]
    M = cm.loc[order, order].values
    ax = sf.add_axes([0.49, 0.2, 0.14, 0.66])
    cmap = LinearSegmentedColormap.from_list("rb", ["#2a78d6", "#a9c8ee", "#ffffff", "#f5b99d", "#eb6834"])
    ax.imshow(M, cmap=cmap, vmin=-0.8, vmax=0.8)
    for i in range(5):
        for j in range(5):
            if i != j:
                ax.text(j, i, f"{M[i, j]:.2f}".replace("-", "−").replace("0.", ".").replace("−.", "−."),
                        ha="center", va="center", fontsize=5.0, color=INK if abs(M[i, j]) < 0.55 else "white")
            else:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, color="#f3f2ee", lw=0))
    ax.set_xticks([]); ax.set_yticks([])
    names = ["Lung", "Bone", "Volume", "Radiod.", "Fat"]
    for k, nm in enumerate(names):
        ax.text(-0.65, k, nm, ha="right", va="center", fontsize=5.0, color=INK)
        ax.text(k, 4.75, nm, ha="right", va="top", fontsize=5.0, color=INK, rotation=40, rotation_mode="anchor")
    for s_ in ax.spines.values():
        s_.set_visible(False)
    ax.text(2, -0.8, "Spearman ρ", ha="center", fontsize=5.4, color=INK2)
    style.source_data(cm.loc[order, order].reset_index(), "P2_Fig2", "h")
    # i: what the score range looks like (rule-chosen controls at the 10th and 90th centiles)
    letter(sf, 0.665, 0.99, "i")
    ex = rd("extra_examples.csv").sort_values("centile")
    for r_i, (_, r) in enumerate(ex.iterrows()):
        y0 = 0.53 - r_i * 0.43
        t10, core, pxb = t10_crop(r.case_id)
        psm, msk, pxm = axial(r.case_id, "T8", half_mm=58, shift_mm=0)
        lung, low, lz, pxl, _ = lung_slice(r.case_id)
        x0 = 0.725
        for k, (im, m, px, ti, win, z) in enumerate([(t10, core, pxb, "bone", (-450, 1050), r.z_bone),
                                                     (psm, msk, pxm, "muscle", (-160, 240), r.z_muscle),
                                                     (lung, None, pxl, "lung", (-1050, -350), r.z_lung)]):
            hh = 0.34
            ww = min(hh * (sf.bbox.height / sf.bbox.width) * im.shape[1] / im.shape[0], 0.078)
            hu = ww * (sf.bbox.width / sf.bbox.height) * im.shape[0] / im.shape[1]
            a2 = sf.add_axes([x0, y0 + (hh - hu) / 2, ww, hu])
            imshow_ct(a2, im, m, color=TIS[ti], lw=0.5, window=win)
            if ti == "lung":
                ov = np.zeros(low.shape + (4,))
                ov[low] = [0.03, 0.47, 0.45, 0.95]
                a2.imshow(ov, interpolation="nearest")
            a2.text(0.5, -0.04, f"z {fnum(z)}", transform=a2.transAxes, ha="center", va="top", fontsize=5.0,
                    color=INK2)
            x0 += ww + 0.008
        sf.text(0.668, y0 + 0.24, f"{int(r.centile)}th\ncentile", fontsize=5.6, va="center", color=INK)
        sf.text(0.668, y0 + 0.07, f"score {fnum(r.S)}", fontsize=5.0, va="center", color=INK2)
    style.source_data(ex, "P2_Fig2", "i")


# ============================================================ Figure 3: e, f
def row_fig3(sf):
    from p2_figs_main import control_sd
    sd = control_sd()
    # e: shift function
    letter(sf, 0.004, 0.99, "e")
    s = rd("extra_shift_function.csv")
    ax = sf.add_axes([0.08, 0.2, 0.36, 0.66])
    for v, col, lab_ in [("S", INK, "Reserve score"), ("z_lung_perc15_hu", TIS["lung"], "Lung density"),
                         ("z_T10_trab_hu", TIS["bone"], "Trabecular bone"), ("z_psm_lvol", TIS["muscle"], "Muscle volume")]:
        q = s[s.measure == v]
        x = q.decile.values * 100
        if v == "S":
            ax.fill_between(x, q.lo.values, q.hi.values, color=INK, alpha=0.1, lw=0)
        ax.plot(x, q["diff"].values, color=col, lw=1.4 if v == "S" else 0.9, zorder=3)
        ax.plot(x, q["diff"].values, "o", ms=3.0 if v == "S" else 2.4, color=col, zorder=4, label=lab_)
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_xlim(5, 95); ax.set_xticks(range(10, 100, 10))
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, fontsize=5.2, handletextpad=0.3, columnspacing=1.0)
    ax.set_xlabel("Decile of the distribution (%)"); ax.set_ylabel("Case minus control (control SD)")
    style.source_data(s, "P2_Fig3", "e")
    # f: offsets by histology
    letter(sf, 0.53, 0.99, "f")
    off = rd("extra_offsets_by_subgroup.csv")
    hs = ["Adenocarcinoma", "Squamous", "Small cell", "Other"]
    ax = sf.add_axes([0.6, 0.2, 0.38, 0.66])
    for k, m in enumerate(["psm_lvol", "T10_trab_hu", "lung_perc15_hu"]):
        q = off[(off.factor == "histo") & (off.measure == m)].set_index("subgroup").reindex(hs)
        xx = np.arange(4) + (k - 1) * 0.2
        for xi, (_, r) in zip(xx, q.iterrows()):
            ax.plot([xi, xi], [r.lo / sd[m], r.hi / sd[m]], color=mcol(m), lw=1.0)
        ax.plot(xx, (q.estimate / sd[m]).values, "o", ms=3.8, color=mcol(m), label=SHORTL[m], zorder=3)
    ax.axhline(0, color=INK2, lw=0.6)
    ns = off[(off.factor == "histo") & (off.measure == "psm_lvol")].set_index("subgroup").reindex(hs).n_cases
    ax.set_xticks(range(4)); ax.set_xticklabels([f"{h}\n({int(n)})" for h, n in zip(hs, ns)], fontsize=5.6)
    ax.set_xlim(-0.5, 3.5); ax.set_ylabel("Level offset (control SD)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=5.2, handletextpad=0.3, columnspacing=1)
    style.source_data(off[off.factor == "histo"], "P2_Fig3", "f")


# ============================================================ Figure 4: c, d
def row_fig4(sf):
    c = rd("extra_cross_vs_long.csv")
    order = ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]
    letter(sf, 0.004, 0.99, "c")
    ax = sf.add_axes([0.2, 0.22, 0.27, 0.62])
    yy = -np.arange(len(order))
    for y, m in zip(yy, order):
        r = c.set_index("measure").loc[m]
        f = 10 / r.control_sd   # control SD per decade
        ax.plot([r.cross_lo * f, r.cross_hi * f], [y + 0.15] * 2, color=mcol(m), lw=0.9)
        ax.plot(r.cross * f, y + 0.15, "o", ms=3.8, mfc="white", mec=mcol(m), mew=0.9)
        ax.plot([r.long_lo * f, r.long_hi * f], [y - 0.15] * 2, color=mcol(m), lw=0.9)
        ax.plot(r.long * f, y - 0.15, "o", ms=3.8, color=mcol(m))
        ax.text(-0.04, y, SHORTL[m], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.6)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_ylim(yy[-1] - 0.6, 0.6)
    ax.set_xlabel("Change per decade (control SD)")
    ax.plot([], [], "o", ms=3.5, mfc="white", mec=INK2, label="cross-sectional (age at entry)")
    ax.plot([], [], "o", ms=3.5, color=INK2, label="longitudinal (within person)")
    ax.legend(loc="lower left", bbox_to_anchor=(-0.6, 1.0), ncol=2, fontsize=5.2)
    style.source_data(c, "P2_Fig4", "c")
    # d: variance explained
    letter(sf, 0.53, 0.99, "d")
    v = rd("extra_variance_explained.csv").set_index("measure").reindex(order)
    ax = sf.add_axes([0.7, 0.22, 0.27, 0.62])
    parts = [("Sex", "#3d3c38"), ("Age", "#7a7871"), ("Smoking", "#b3b1a8"), ("Scanner", "#dcdad2")]
    for y, m in zip(yy, order):
        left = 0
        for name, col in parts:
            w = max(v.loc[m, name], 0) * 100
            ax.barh(y, w, left=left, height=0.6, color=col, lw=0.4, edgecolor="white")
            left += w
        ax.text(v.loc[m, "full"] * 100 + 1.5, y, f"R² {v.loc[m, 'full'] * 100:.0f}%", va="center", fontsize=5.3)
        ax.text(-0.04, y, SHORTL[m], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.6)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 70); ax.set_ylim(yy[-1] - 0.6, 0.6)
    ax.set_xlabel("Between-person variance explained (%, unique share)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=col, label=name) for name, col in parts], loc="lower left",
              bbox_to_anchor=(-0.55, 1.0), ncol=4, fontsize=5.2, handlelength=1.0)
    style.source_data(v.reset_index(), "P2_Fig4", "d")


# ============================================================ Figure 5: f, g
def row_fig5(sf):
    pr = rd("extra_projected_reliability.csv")
    letter(sf, 0.004, 0.99, "f")
    ax = sf.add_axes([0.08, 0.2, 0.34, 0.66])
    ax.axhspan(0, 0.4, color="#f7f6f2", lw=0, zorder=0)
    for yv, txt in [(0.4, "0.40"), (0.7, "0.70")]:
        ax.axhline(yv, color=INK2, lw=0.5, ls=(0, (3, 2)))
        ax.text(0.8, yv + 0.012, txt, fontsize=5.0, color=INK2, va="bottom")
    for m in ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]:
        q = pr[pr.measure == m]
        mf, ls = mstyle(m)
        ax.plot(q.years, q.reliability, color=mcol(m), lw=1.1, ls=ls)
        ax.plot(q.years, q.reliability, "o", ms=2.4, color=mcol(m), mfc=mf, mec=mcol(m), mew=0.7, label=SHORTL[m])
    ax.axvline(2, color=CASE, lw=0.7)
    ax.text(2.15, 0.03, "NLST", fontsize=5.0, color=CASE, va="bottom")
    ax.set_xlim(0.7, 10.5); ax.set_ylim(0, 1)
    ax.legend(loc="lower right", fontsize=5.0, handlelength=1.6, labelspacing=0.3, borderaxespad=0.4)
    ax.set_xticks([1, 2, 4, 6, 8, 10])
    ax.set_xlabel("Years of annual low-dose CT"); ax.set_ylabel("Projected slope reliability")
    style.source_data(pr, "P2_Fig5", "f")
    # g: kernel switch, individual radiodensity slopes
    letter(sf, 0.53, 0.99, "g")
    k = rd("extra_kernel_change.csv")
    ax = sf.add_axes([0.66, 0.2, 0.32, 0.66])
    a = k[k.kernel_consistent == 1].slope.clip(-12, 8).values
    b = k[k.kernel_consistent == 0].slope.clip(-12, 8).values
    raincloud(ax, [a, b], [TIS["muscle"], INK2], ["kept", "switched"], [1.0, 0.0], bw=0.25, cloud_h=0.55,
              xlim=(-12, 8), point_alpha=0.25)
    ax.set_yticks([0.12, 1.12])
    ax.set_yticklabels([f"Kernel switched\n(n = {len(b)})", f"Same kernel\n(n = {len(a):,})"], fontsize=5.5)
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.42, 1.75); ax.set_xlim(-12, 8)
    ax.axvline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlabel("Individual radiodensity slope, controls (HU/year)")
    style.source_data(k, "P2_Fig5", "g")


# ============================================================ Figure 6: e, f
def row_fig6(sf):
    from p2_figs_main import control_sd
    sd = control_sd()
    div = rd("partA_divergence.csv")
    order = ["psm_fat_pct", "psm_hu", "psm_lvol", "T10_trab_hu", "lung_perc15_hu"]
    letter(sf, 0.004, 0.99, "e")
    ax = sf.add_axes([0.26, 0.2, 0.14, 0.62])
    rows = []
    for m in order:
        q = div[(div.measure == m) & (div.term == "groupcase:tau") & (div.analysis == "mixed model")].iloc[0]
        e, lo, hi = q.estimate / sd[m], q.lo / sd[m], q.hi / sd[m]
        rows.append(dict(kind="row", label=SHORTL[m], est=e, lo=lo, hi=hi, color=mcol(m),
                         hollow=m in HOLLOW, txt=f"{fnum(e, '{:.3f}')} ({fnum(lo, '{:.3f}')} to {fnum(hi, '{:.3f}')})"))
    forest(ax, rows, xlim=(-0.035, 0.035), ticks=[-0.03, 0, 0.03], label_w=1.75, ref=0, log=False,
           right_cols=(("Per year (95% CI)", 1.06),), xlabel="Case − control slope (control SD/year)", fs=5.4)
    ax.set_xticklabels([fnum(t, "{:g}") for t in [-0.03, 0, 0.03]])
    style.source_data(pd.DataFrame(rows), "P2_Fig6", "e")
    # f: largest true acceleration compatible with the data
    letter(sf, 0.6, 0.99, "f")
    key = rd("p2_key_numbers.csv").set_index("measure").reindex(order)
    ax = sf.add_axes([0.75, 0.2, 0.22, 0.62])
    yy = -np.arange(len(order))
    for y, m in zip(yy, order):
        v = key.loc[m, "max_upper_bound_true_sd"]
        ax.barh(y, v, height=0.56, color=mcol(m), alpha=0.25, edgecolor=mcol(m), lw=0.8,
                hatch="////" if m in HOLLOW else None)
        ax.text(v + 0.015, y, f"{v:.2f}", va="center", fontsize=5.3, color=INK2)
        ax.text(-0.04, y, SHORTL[m], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=5.6)
    for xv, txt in [(0.2, "small"), (0.5, "moderate")]:
        ax.axvline(xv, color=INK2, lw=0.5, ls=(0, (3, 2)))
        ax.text(xv, 0.75, txt, fontsize=4.9, color=INK2, ha="center", va="bottom")
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 0.75); ax.set_ylim(yy[-1] - 0.6, 0.9)
    ax.set_xlabel("Largest compatible acceleration\n(true-change SD, upper 95% limit)")
    style.source_data(key.reset_index()[["measure", "max_upper_bound_z", "max_upper_bound_true_sd"]], "P2_Fig6", "f")
