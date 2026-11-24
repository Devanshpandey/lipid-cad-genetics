"""One figure theme for every figure of the lipid-axes paper (revision r4).

Geometry and type follow the Nature Communications and Nature figure guides: 180 mm wide,
at most 170 mm tall for a main figure, Arial, panel letters 8 pt bold lowercase, labels 7 pt,
ticks and keys 6 pt, nothing below 5 pt; open L-shaped axes, no gridlines, CIs without caps,
ratios on log axes, a solid grey null line behind the data.

One encoding everywhere:
  hue      the exposure, fixed across figures (Okabe-Ito; validated with the dataviz skill's
           validate_palette.js --pairs all: worst CVD dE 11.0, normal-vision 15.6)
  shape    what the estimate describes within a figure (stage, cohort or MR method), named in its key;
           a diamond is always a randomized trial
  fill     filled = the 95% CI excludes the null (or the significance rule named in the legend);
           open = it does not
  text     always ink or muted grey, never the data colour
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator

MM = 72 / 25.4
W = 180 * MM                    # 510.24 pt
H_MAX = 170 * MM                # 481.9 pt

FS = dict(letter=8, title=7, label=7, tick=6, key=6, note=6, small=5.5)
LW = dict(axis=0.6, ci=1.0, ref=0.6, link=0.5)
MS = 3.8                        # marker diameter, pt

INK = "#1A1A1A"
MUTED = "#595959"
NULL = "#8C8C8C"
BAND = "#F2F2F2"
NS = "#BDBDBD"
TINT = "#EDEDED"

TRAIT = {
    "LDL-C": "#0072B2", "ApoB": "#56B4E9", "Lp(a)": "#D55E00",
    "TG": "#E69F00", "HDL-C": "#009E73", "non-lipid": "#6E6E6E",
}
TRIAL = INK


def tint(hex_, f):
    """f of the colour over white (f = 0.18 gives a fill for boxes)."""
    h = hex_.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return "#%02X%02X%02X" % tuple(round(255 - f * (255 - c)) for c in (r, g, b))


matplotlib.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS["tick"], "axes.labelsize": FS["label"], "axes.titlesize": FS["title"],
    "xtick.labelsize": FS["tick"], "ytick.labelsize": FS["tick"], "legend.fontsize": FS["key"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
    "axes.linewidth": LW["axis"], "xtick.major.width": LW["axis"], "ytick.major.width": LW["axis"],
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5, "ytick.minor.size": 1.5, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "lines.linewidth": LW["ci"], "lines.markersize": MS, "errorbar.capsize": 0,
    "legend.frameon": False, "legend.handlelength": 1.2, "legend.handletextpad": 0.4,
    "legend.labelspacing": 0.3, "legend.borderaxespad": 0.3,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.unicode_minus": True, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 300, "savefig.bbox": "standard", "figure.facecolor": "white", "savefig.facecolor": "white",
})


# ------------------------------------------------------------------ page and layout, in points from the top left
class Page:
    def __init__(self, h, w=W):
        assert h <= H_MAX + 0.5 or w < W, "a main figure is at most 170 mm tall"
        self.w, self.h = w, h
        self.fig = plt.figure(figsize=(w / 72, h / 72))

    def ax(self, x, ytop, w, h):
        return self.fig.add_axes([x / self.w, 1 - (ytop + h) / self.h, w / self.w, h / self.h])

    def text(self, x, y, s, **kw):
        kw.setdefault("color", INK)
        kw.setdefault("va", "baseline")
        kw.setdefault("ha", "left")
        kw.setdefault("fontsize", FS["note"])
        return self.fig.text(x / self.w, 1 - y / self.h, s, **kw)

    def panel(self, x, y, letter, title=None):
        """Panel letter (8 pt bold) and a short title (7 pt), on one baseline."""
        self.text(x, y, letter, fontsize=FS["letter"], fontweight="bold")
        if title:
            self.text(x + 9, y, title, fontsize=FS["title"], fontweight="bold")

    def key(self, x, y, items, gap=10.0):
        """One row of key entries, left to right: (marker, facecolour, edgecolour, label)."""
        rend = self.fig.canvas.get_renderer()
        for m, fc, ec, lab in items:
            if m == "line":
                self.fig.add_artist(Line2D([(x) / self.w, (x + 9) / self.w], [1 - (y - 2) / self.h] * 2,
                                           transform=self.fig.transFigure, color=ec, lw=LW["ci"]))
                x += 11
            else:
                self.fig.add_artist(Line2D([(x + 3) / self.w], [1 - (y - 2) / self.h], transform=self.fig.transFigure,
                                           ls="none", marker=m, ms=MS, mfc=fc, mec=ec, mew=0.9 if fc == "white" else 0.5))
                x += 8
            t = self.text(x, y, lab, fontsize=FS["key"])
            x += t.get_window_extent(rend).width * 72 / self.fig.dpi + gap
        return x

    def save(self, stem):
        os.makedirs(os.path.dirname(stem), exist_ok=True)
        self.fig.savefig(stem + ".pdf")
        self.fig.savefig(stem + ".png", dpi=300)
        plt.close(self.fig)


# ------------------------------------------------------------------ axes helpers
def clean(ax, left=False, bottom=True):
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    if not left:
        ax.tick_params(axis="y", length=0)


def logx(ax, lo, hi, ticks):
    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FixedFormatter([("%g" % t) for t in ticks]))
    ax.xaxis.set_minor_locator(NullLocator())


def null(ax, at=1.0, axis="x"):
    (ax.axvline if axis == "x" else ax.axhline)(at, color=NULL, lw=LW["ref"], zorder=1)


def band(ax, y0, y1):
    ax.axhspan(y0, y1, color=BAND, lw=0, zorder=0)


def point(ax, y, est, lo, hi, colour, marker="s", filled=True, ms=MS, ci_colour=None, z=3, lw=None):
    """A 95% CI without caps, then the marker: filled with a white ring, or open."""
    ax.plot([lo, hi], [y, y], color=ci_colour or colour, lw=lw or LW["ci"], solid_capstyle="butt", zorder=z)
    if filled:
        ax.plot([est], [y], marker=marker, ms=ms, mfc=colour, mec="white", mew=0.5, ls="none", zorder=z + 1)
    else:
        ax.plot([est], [y], marker=marker, ms=ms - 0.3, mfc="white", mec=colour, mew=0.9, ls="none", zorder=z + 1)


def excludes(lo, hi, at=1.0):
    return hi < at or lo > at


def ci_text(est, lo, hi, d=2):
    f = "%." + str(d) + "f"
    return (f % est) + " (" + (f % lo) + "–" + (f % hi) + ")"


def p_text(p, sig=1):
    """P in mathtext: 0.03, 9.6 x 10^-9."""
    if p >= 0.001:
        return "%.2g" % p if p < 0.1 else "%.2f" % p
    m, e = ("%.*e" % (sig - 1 if sig > 1 else 0, p)).split("e")
    return r"%s$\times$10$^{%d}$" % (m, int(e))


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def num(v):
    return float(v) if v not in ("", "NA", None) else float("nan")


def lower(b, se):
    """Hazard ratio per unit LOWER exposure, with its 95% interval."""
    return math.exp(-b), math.exp(-b - 1.96 * se), math.exp(-b + 1.96 * se)


SRC = []


def src(fig, panel, row, series, est, lo, hi, **extra):
    SRC.append(dict(figure=fig, panel=panel, row=row, series=series, estimate="%.4f" % est,
                    lower="%.4f" % lo, upper="%.4f" % hi, **extra))


def write_src(path):
    keys = ["figure", "panel", "row", "series", "estimate", "lower", "upper"]
    extra = sorted({k for r in SRC for k in r} - set(keys))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys + extra)
        w.writeheader()
        w.writerows(SRC)
