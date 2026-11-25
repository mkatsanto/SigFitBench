"""Overview figure of the SigFitBench workflow: a simple version of plot_workflow_schematic.py.

Real or simulated mutation catalogues go to the fitting tools, which all use the same COSMIC
signatures. With real data the tools are compared with each other; with simulated data
their results are scored against the known signature activities, which the tools never see.

The small spectra and bars are drawn from the data in the repo, as in the detailed figure.

Usage (from the SigFitBench root):
    python scripts/plot_workflow_overview.py --out sigfitbench_workflow_overview.pdf
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

from plot_workflow_schematic import (  # shared style, icons and example data
    BLUE, INK, MM, MUTED, ORANGE, ORANGE_TXT, PURPLE, RED, TEAL, activities, arrow,
    icon_heatmap, icon_performance, icon_reference, load_examples, path, rbox, spectrum, txt)

W, H = 183.0, 78.0  # mm

BLUE_FILL = "#E9F1F9"
PURPLE_FILL = "#F2EBF6"
ORANGE_FILL = "#FDF1E4"
FRAME_FILL = "#F7F7F7"

TOOLS = ["SigProfilerAssignment", "MuSiCal", "MutationalPatterns", "deconstructSigs",
         "sigfit", "mmsig", "sigminer", "SigsPack", "YAPSA", "sigLASSO", "signature.tools.lib",
         "SigTrim", "SigID"]
NEW_METHODS = ["SigTrim", "SigID"]

T_TITLE, T_TEXT, T_TAG = 9.5, 7.3, 6.5


def card(ax, x, y, w, h, fill):
    rbox(ax, x, y, w, h, fill, "none", lw=0, r=2.4, z=1)


def tag(ax, x, y, text, color):
    """Small rounded label, left edge at x."""
    wd = 2.6 + 1.22 * len(text)
    rbox(ax, x, y - 2.2, wd, 4.4, "white", color, lw=0.6, r=2.2, z=4)
    txt(ax, x + wd / 2, y + 0.05, text, size=T_TAG, color=color, weight="bold")


def tool_dots(ax, cx, cy, n_known, new_names):
    """One dot per tool; the new methods are larger, red and named."""
    per_row = 7
    names = ["known"] * n_known + new_names
    rows = (len(names) + per_row - 1) // per_row
    gap = 3.4
    for i, name in enumerate(names):
        r, c = divmod(i, per_row)
        in_row = min(per_row, len(names) - r * per_row)
        x = cx + (c - (in_row - 1) / 2) * gap
        y = cy + (r - (rows - 1) / 2) * gap
        new = name != "known"
        ax.add_patch(Circle((x, y), 1.25 if new else 1.0, fc=RED if new else ORANGE,
                            ec="white", lw=0.6, zorder=4))


def draw(ax, args):
    real_spec, sig_names, sig_vals, unk, expected, counts = load_examples(args)
    new = [t for t in args.tools if t in args.new_methods]
    n_known = len(args.tools) - len(new)

    # Workflow frame
    rbox(ax, 2.5, 2.5, W - 5, H - 5, FRAME_FILL, "none", lw=0, r=3.0, z=0)
    txt(ax, 7, 7.2, "SigFitBench", size=T_TITLE, weight="bold", ha="left")
    txt(ax, 30.5, 7.2, "one config file · one Snakemake run", size=T_TEXT - 0.5,
        color=MUTED, ha="left")

    lx, lw_ = 7.0, 56.0          # inputs
    fx, fw = 75.0, 36.0          # fitting
    rx, rw = 123.0, 53.0         # results
    top, h1 = 13.0, 23.0         # real-data row
    bot, h2 = 40.0, 27.0         # simulated-data row
    y_real, y_sim = top + h1 / 2, bot + h2 / 2

    # --- Real data ---
    card(ax, lx, top, lw_, h1, PURPLE_FILL)
    spectrum(ax, lx + 3.5, top + 5.5, 19, 9, real_spec)
    txt(ax, lx + 26, top + 6.5, "Real data", size=T_TITLE, weight="bold", color=PURPLE,
        ha="left")
    txt(ax, lx + 26, top + 11.5, "PCAWG, your own\ncatalogue or VCF", size=T_TEXT,
        color=MUTED, ha="left", va="top")
    tag(ax, lx + 3.5, top + h1 - 4.0, "no ground truth", PURPLE)

    # --- Simulated data ---
    card(ax, lx, bot, lw_, h2, BLUE_FILL)
    activities(ax, lx + 3.5, bot + 4.5, 9.5, 9, sig_names, list(sig_vals), unk,
               label_size=6.0, show_names=False)
    arrow(ax, (lx + 14.0, bot + 9.5), (lx + 17.2, bot + 9.5), color=BLUE, ms=5)
    spectrum(ax, lx + 18.0, bot + 4.5, 12, 9, counts)
    txt(ax, lx + 32, bot + 6.0, "Simulated data", size=T_TITLE, weight="bold", color=BLUE,
        ha="left")
    txt(ax, lx + 32, bot + 11.0, "real tumour\nmixtures, then\nN mutations drawn",
        size=T_TEXT - 0.4, color=MUTED, ha="left", va="top")
    tag(ax, lx + 3.5, bot + h2 - 4.0, "ground truth", BLUE)

    # --- Fitting ---
    fy, fh = top + 4.0, bot + h2 - top - 8.0
    card(ax, fx, fy, fw, fh, ORANGE_FILL)
    icon_reference(ax, fx + fw / 2 - 1.5, fy + 9.5)
    txt(ax, fx + fw / 2, fy + 20.5, "Fit signatures", size=T_TITLE, weight="bold",
        color=ORANGE_TXT)
    txt(ax, fx + fw / 2, fy + 25.5, f"{len(args.tools)} tools, same\nCOSMIC reference",
        size=T_TEXT, color=MUTED, va="top")
    tool_dots(ax, fx + fw / 2, fy + fh - 7.0, n_known, new)
    if new:
        txt(ax, fx + fw / 2, fy + fh - 2.2, "new: " + ", ".join(new), size=T_TAG,
            color=RED)

    # Inputs -> fitting -> results
    ym = fy + fh / 2
    path(ax, [(lx + lw_, y_real), (lx + lw_ + 5.5, y_real), (lx + lw_ + 5.5, ym - 3),
              (fx, ym - 3)], color=PURPLE, lw=1.3)
    path(ax, [(lx + lw_, y_sim), (lx + lw_ + 5.5, y_sim), (lx + lw_ + 5.5, ym + 3),
              (fx, ym + 3)], color=BLUE, lw=1.3)
    path(ax, [(fx + fw, ym - 3), (fx + fw + 5.5, ym - 3), (fx + fw + 5.5, y_real),
              (rx, y_real)], color=PURPLE, lw=1.3)
    path(ax, [(fx + fw, ym + 3), (fx + fw + 5.5, ym + 3), (fx + fw + 5.5, y_sim),
              (rx, y_sim)], color=BLUE, lw=1.3)

    # --- Results ---
    card(ax, rx, top, rw, h1, PURPLE_FILL)
    icon_heatmap(ax, rx + 10.5, y_real, s=2.4)
    txt(ax, rx + 20, top + 7.0, "Compare tools", size=T_TITLE, weight="bold", color=PURPLE,
        ha="left")
    txt(ax, rx + 20, top + 12.0, "how well they\nagree with each other", size=T_TEXT,
        color=MUTED, ha="left", va="top")

    card(ax, rx, bot, rw, h2, BLUE_FILL)
    icon_performance(ax, rx + 10.5, y_sim, s=1.0)
    txt(ax, rx + 20, bot + 7.0, "Score accuracy", size=T_TITLE, weight="bold", color=BLUE,
        ha="left")
    txt(ax, rx + 20, bot + 12.0, "error, precision\nand recall against\nthe ground truth",
        size=T_TEXT, color=MUTED, ha="left", va="top")

    # Ground truth goes straight to scoring, past the tools
    gy = bot + h2 + 3.6
    path(ax, [(lx + lw_ / 2, bot + h2), (lx + lw_ / 2, gy), (rx + rw / 2, gy),
              (rx + rw / 2, bot + h2)], color=TEAL, lw=1.1, ls=(0, (3, 1.8)))
    txt(ax, fx + fw / 2, gy - 1.9, "ground truth, never shown to the tools", size=T_TAG,
        color=TEAL, style="italic")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tissue-dir", default="cosmic_tissue_data")
    ap.add_argument("--cosmic", default="input/COSMIC_v3.3.1_SBS_GRCh38.txt")
    ap.add_argument("--real-catalog",
                    default="real_mutational_catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat")
    ap.add_argument("--real-type", default="Skin-Melanoma")
    ap.add_argument("--sim-type", default="Breast")
    ap.add_argument("--unknown-sig", default="SBS90")
    ap.add_argument("--sim-n", type=int, default=500)
    ap.add_argument("--tools", nargs="+", default=TOOLS, help="Tools that are run")
    ap.add_argument("--new-methods", nargs="*", default=NEW_METHODS,
                    help="Tools to highlight as the new methods")
    ap.add_argument("--out", default="sigfitbench_workflow_overview.pdf")
    ap.add_argument("--png", default=None, help="Optional PNG preview path")
    args = ap.parse_args()

    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    draw(ax, args)

    fig.savefig(args.out, dpi=600)
    if args.png:
        fig.savefig(args.png, dpi=220)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
