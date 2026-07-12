"""figstyle.py — shared publication house style for all manuscript figures.

Importing this module applies a clean, striking theme (Helvetica Neue,
refined typography, thin elegant axes) and exports a colorblind-safe palette
plus helpers. Import at the top of every figure script:

    from figstyle import CB, INK, GREY, panel, title_block, save, despine
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---- typography ----
mpl.rcParams.update({
    "font.family": "Helvetica Neue",
    "font.size": 9.5,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlepad": 9,
    "axes.labelsize": 9.5, "axes.labelweight": "regular",
    "axes.linewidth": 0.9,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3.2, "ytick.major.size": 3.2,
    "xtick.major.width": 0.9, "ytick.major.width": 0.9,
    "xtick.major.pad": 3, "ytick.major.pad": 3,
    "legend.frameon": False, "legend.fontsize": 8,
    "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.06,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "axes.grid": False,
    "axes.unicode_minus": False,  # ASCII hyphen for PDF text-layer extraction
    "pdf.fonttype": 42, "ps.fonttype": 42,  # embed TrueType (editable/searchable text)
})

# ---- palette (Wong colorblind-safe, refined) ----
INK  = "#20303C"   # near-black slate for text / spines
GREY = "#8A94A6"   # muted grey for secondary elements
FAINT = "#E7EBF0"  # very light fill
CB = {
    "blue":   "#1F77B4",   # protective / cool
    "red":    "#D1495B",   # risk / warm
    "green":  "#2A9D8F",   # evidence
    "amber":  "#E9A13B",   # drug / highlight
    "purple": "#8E6BAF",
    "sky":    "#4FA6D6",
    "rose":   "#C77CA0",
    "slate":  INK,
    "grey":   GREY,
    "ink":    INK,
    "yellow": "#E9C46A",
    "orange": "#E9A13B",   # alias of amber
}
mpl.rcParams.update({
    "text.color": INK, "axes.labelcolor": INK,
    "axes.edgecolor": INK, "xtick.color": INK, "ytick.color": INK,
})

# ---- GLOBAL FIGURE GRAMMAR (one visual language across all figures) ----
# Fixed biological colors: same trait axis -> same color everywhere.
SUBTLE = "#5B6675"   # darker grey for subtitles (stays legible in print/PDF)
LIPID = {
    "LDL-C":  "#D1495B",  # LDL-C / ApoB axis  -> red
    "ApoB":   "#D1495B",
    "non-HDL-C": "#D1495B",
    "TG":     "#E08A2B",  # triglyceride / remnant axis -> orange
    "Triglycerides": "#E08A2B",
    "remnant": "#E08A2B",
    "HDL-C":  "#1F77B4",  # HDL-C / ApoA1 axis -> blue
    "ApoA1":  "#1F77B4",
    "Lp(a)":  "#8E6BAF",  # Lp(a) -> purple
    "rare":   "#2A9D8F",  # rare coding evidence -> teal
    "outcome": "#7B2D42",  # CAD / outcome evidence -> burgundy
    "ns":     "#C3C9D2",  # nonsignificant / background -> light grey
}
# Fixed evidence shapes: source encoded by marker, not only colour.
SHAPE = {
    "discovery": "o",   # UK Biobank discovery
    "eu_rep":    "D",   # independent European replication
    "noneu":     "s",   # non-European analysis
    "rare":      "^",   # rare-variant evidence
}


def lipid_color(name):
    """Return the fixed biological colour for a trait label (falls back to grey)."""
    return LIPID.get(name, LIPID["ns"])


def despine(ax, left=True, bottom=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)


def panel(ax, letter, dx=-0.02, dy=1.06):
    """Bold panel label (a, b, c) in figure-corner style."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top", ha="right", color=INK)


def title_block(ax, title, subtitle=None, y=1.10, sy=1.045):
    """Bold title with an optional grey subtitle for clear hierarchy."""
    ax.set_title("")
    ax.text(0.5, y, title, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=11, fontweight="bold", color=INK)
    if subtitle:
        ax.text(0.5, sy, subtitle, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=8.2, color=GREY)


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}")
    print(f"wrote {stem}.png / .pdf")
