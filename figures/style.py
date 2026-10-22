"""Shared publication style (Nature / Cell / Science artwork specifications).

Sizes: single column 89 mm, 1.5 column 120 mm, double column 183 mm, max height 247 mm.
Type: Helvetica; 7 pt labels, 6 pt ticks/annotations, 8 pt bold lowercase panel letters.
Lines: axes 0.5 pt, data lines 1 pt; no gridlines; outward ticks. Export: vector PDF with embedded
TrueType text (editable in Illustrator) + 600-dpi PNG, and one Source Data CSV per panel.

Colours (validated with the dataviz validator on a white surface, --pairs all):
  SERIES  blue / orange / aqua          worst CVD dE 9.2, normal dE 24.0  (aqua < 3:1 -> always labelled)
  TISSUE  muscle red / SAT yellow / internal-fat blue / IMAT violet
                                        worst CVD dE 13.0, normal dE 16.3 (yellow < 3:1 -> always labelled)
Series colours and tissue colours never share a figure panel with different meanings.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

MM = 1 / 25.4
W1, W15, W2 = 89 * MM, 120 * MM, 183 * MM

INK = "#0b0b0b"         # primary text / key marks
INK2 = "#52514e"        # secondary text
MUTED = "#898781"       # axis labels of lesser importance, reference marks
HAIR = "#c3c2b7"        # axes, reference lines
WASH = "#e1e0d9"        # light fills / backgrounds
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
TISSUE = {"muscle": "#e34948", "sat": "#eda100", "ift": "#2a78d6", "imat": "#4a3aa7"}
TISSUE_LABEL = {"muscle": "Skeletal muscle", "sat": "Subcutaneous fat", "ift": "Visceral/mediastinal fat",
                "imat": "Intramuscular fat"}
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "figures", "out")


def setup():
    for f in os.listdir(os.path.expanduser("~/.fonts")) if os.path.isdir(os.path.expanduser("~/.fonts")) else []:
        if f.lower().startswith("helvetica"):
            font_manager.fontManager.addfont(os.path.join(os.path.expanduser("~/.fonts"), f))
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Noto Sans CJK SC", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "legend.fontsize": 6, "legend.frameon": False, "legend.handlelength": 1.2, "legend.borderaxespad": 0.2,
        "axes.linewidth": 0.5, "axes.edgecolor": INK, "axes.labelcolor": INK, "axes.titleweight": "normal",
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "axes.labelpad": 2,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.major.pad": 1.5, "ytick.major.pad": 1.5, "xtick.color": INK, "ytick.color": INK,
        "xtick.direction": "out", "ytick.direction": "out",
        "lines.linewidth": 1.0, "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round",
        "patch.linewidth": 0.5, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 600, "figure.dpi": 150, "savefig.bbox": "standard",
        # italic P and superscripts in Helvetica, not DejaVu
        "mathtext.fontset": "custom", "mathtext.rm": "Helvetica", "mathtext.it": "Helvetica:italic",
        "mathtext.bf": "Helvetica:bold", "mathtext.sf": "Helvetica", "mathtext.default": "regular",
        "axes.unicode_minus": True,
    })
    for f in os.listdir(os.path.join(ROOT, "figures", "fonts")) if os.path.isdir(os.path.join(ROOT, "figures", "fonts")) else []:
        font_manager.fontManager.addfont(os.path.join(ROOT, "figures", "fonts", f))  # e.g. Noto Sans CJK SC Bold
    from i18n import LANG
    if LANG == "zh":  # Chinese review version: the Helvetica file maps CJK to empty boxes, so lead with Noto CJK
        mpl.rcParams["font.sans-serif"] = ["Noto Sans CJK SC", "Droid Sans Fallback", "DejaVu Sans"]


def fnum(v, fmt="{:.2f}"):
    """Number formatted with a true minus sign (U+2212), as journals typeset it."""
    return fmt.format(v).replace("-", "\u2212")


def panel_label(ax_or_fig, text, x=0.0, y=1.0, fig=None, **kw):
    """Bold lowercase panel letter at the top-left of an axes (axes-fraction) or figure (figure-fraction)."""
    if fig is not None:
        fig.text(x, y, text, fontsize=8, fontweight="bold", va="top", ha="left", color=INK, **kw)
    else:
        ax_or_fig.text(x, y, text, transform=ax_or_fig.transAxes, fontsize=8, fontweight="bold",
                       va="bottom", ha="right", color=INK, **kw)


def save(fig, name, synthetic=False):
    """Write <name>.pdf and <name>.png. Synthetic layout tests go to the scratch folder, never figures/out.
    Chinese review versions (FIG_LANG=zh) get the suffix _zh."""
    from i18n import LANG
    if LANG == "zh":
        name += "_zh"
    out = OUT
    if synthetic:
        watermark(fig)
        out = os.environ.get("FIG_SCRATCH", "/tmp")
    os.makedirs(out, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"))
    plt.close(fig)
    return os.path.join(out, f"{name}.pdf")


def watermark(fig):
    from i18n import T
    fig.text(0.5, 0.5, T("SYNTHETIC DATA — LAYOUT TEST ONLY"), fontsize=22, color="#d03b3b", alpha=0.35,
             ha="center", va="center", rotation=25, fontweight="bold", zorder=1000)


def source_data(df, fig_name, panel):
    """Nature requires Source Data for every graph: one CSV per panel next to the figure."""
    from i18n import LANG
    if LANG == "zh":
        return  # the English run writes the source data
    d = os.path.join(OUT, "source_data")
    os.makedirs(d, exist_ok=True)
    df.to_csv(os.path.join(d, f"{fig_name}_{panel}.csv"), index=False)


def despine(ax, left=True, bottom=True):
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    if not left:
        ax.tick_params(left=False, labelleft=False)
    if not bottom:
        ax.tick_params(bottom=False, labelbottom=False)
