"""Overview figure of SigFitBench: the workflow on a central band, with numbered panels around it.

Center: the Snakemake workflow as labeled steps. Simulated catalogs (upper path) or real
catalogs (lower path) are fitted by every tool with the same COSMIC signatures; simulated
results are scored against the known signature activities, real results are compared
between tools, everything ends in one HTML report, and GitHub Actions rebuilds a live
version of it from simulated data.

Around it, each numbered like the step it explains: the simulation, the real-data inputs,
the fitting tools, the metrics and the report, plus how the workflow meets the FAIR
principles. The spectra and the two report thumbnails are drawn from data and results
already in the repo; the thumbnails show only the tools listed in the figure.

Usage (from the SigFitBench root):
    python scripts/plot_workflow_figure.py --out sigfitbench_workflow_figure.pdf
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

from plot_workflow_schematic import (  # shared helpers
    INK, MM, MUTED, SBS_COL, rbox, spectrum, txt)

W, H = 183.0, 112.0  # mm

# palette
SIM = "#2A9D8F"       # simulated data and scoring against the truth
REAL = "#E76F51"      # real data and comparing tools
REF = "#C9A227"       # COSMIC reference
FIT = "#3D5A80"       # fitting tools
METRIC = "#6C757D"    # metrics: same grey as the FAIR panel, so the green and orange inside stand out
REPORT = "#8E5572"    # report
LIVE = "#2B2D42"      # live report built by GitHub Actions
NEUTRAL = "#6C757D"   # FAIR panel
NEW = "#D62828"       # new method in the tool list
BAND = "#EEF2F6"
BAR_COLS = [FIT, SIM, REF, REPORT, "#98C1D9"]
THUMB_COLS = [SIM, FIT, REF, NEW, REPORT, REAL]
FAIR_COLS = {"F": "#1F4E79", "A": "#2E7D32", "I": "#A66A00", "R": "#9E2A2B"}

TOOLS = ["SigProfilerAssignment", "MuSiCal", "MutationalPatterns", "deconstructSigs",
         "sigfit", "mmsig", "sigminer", "SigsPack", "YAPSA", "sigLASSO", "signature.tools.lib",
         "mSigAct", "StarSignDNA", "SPA-tree"]
NEW_METHODS = ["SPA-tree"]
# names the result tables use for tools that are listed under another name
RESULT_NAMES = {"SigProfilerAssignment": "SPA", "signature.tools.lib": "signature_tools",
                "SPA-tree": "SPA_tree"}

# type sizes (pt) at final print size: 8 for titles, 7 for headings and key labels,
# and nothing below 6, the smallest size most journals accept
T_TITLE, T_HEAD, T_TEXT, T_SMALL, T_TINY = 8.0, 7.0, 7.0, 6.0, 6.0


def tint(color, a=0.12):
    """Color mixed with white, as a solid fill."""
    r, g, b, _ = mpl.colors.to_rgba(color)
    return (1 - a + a * r, 1 - a + a * g, 1 - a + a * b)


# ---------------------------------------------------------------------------
# Building blocks (mm, origin top-left)
# ---------------------------------------------------------------------------
def badge(ax, x, y, n, color, r=1.9, size=T_SMALL):
    ax.add_patch(Circle((x, y), r, fc=color, ec="white", lw=0.6, zorder=9))
    txt(ax, x, y + 0.08, str(n), size=size, weight="bold", color="white", z=10)


def panel(ax, x, y, w, h, color, n, title):
    """Tinted panel with a numbered badge (or none) and its title."""
    rbox(ax, x, y, w, h, tint(color, 0.10), "none", lw=0, r=2.0, z=1)
    ax.add_patch(Rectangle((x, y + 1.6), 1.1, h - 3.2, fc=color, ec="none", zorder=2))
    if n is None:
        txt(ax, x + 3.6, y + 4.0, title, size=T_HEAD, weight="bold", color=color, ha="left")
    else:
        badge(ax, x + 5.0, y + 4.0, n, color, r=2.0, size=T_SMALL + 0.3)
        txt(ax, x + 8.2, y + 4.0, title, size=T_HEAD, weight="bold", color=color, ha="left")


def chip(ax, x, y, text, color, size=T_TINY):
    """Rounded label with its left edge at x and center at y; returns its right edge."""
    wd = 1.8 + 0.188 * len(text) * size
    rbox(ax, x, y - 1.6, wd, 3.2, "white", "none", lw=0, r=1.6, z=3)
    txt(ax, x + wd / 2, y + 0.05, text, size=size, color=color, z=4)
    return x + wd


def pill(ax, cx, cy, w, text, color, n=None, sub=None, outline=False, h=6.0):
    """One workflow step: a rounded box with its name, an optional number and a note below."""
    if outline:
        rbox(ax, cx - w / 2, cy - h / 2, w, h, "white", INK, lw=0.8, r=h / 2, z=6)
        txt(ax, cx, cy + 0.08, text, size=T_SMALL, color=INK, z=7)
    else:
        rbox(ax, cx - w / 2, cy - h / 2, w, h, color, "none", lw=0, r=h / 2, z=6)
        txt(ax, cx, cy + 0.08, text, size=T_TEXT, weight="bold", color="white", z=7)
    if n is not None:
        badge(ax, cx - w / 2 + 0.5, cy - h / 2 + 0.3, n, color, r=1.7, size=T_TINY)
    if sub:
        txt(ax, cx, cy + h / 2 + 2.0, sub, size=T_TINY, color=MUTED, z=7)
    return (cx - w / 2, cx + w / 2)


def link(ax, p0, p1, color=MUTED, ls="-", rad=0.0, lw=0.8):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=5.5, color=color,
                                 lw=lw, ls=ls, shrinkA=0, shrinkB=0,
                                 connectionstyle=f"arc3,rad={rad}", zorder=5))


def bars(ax, x, y, w, h, vals):
    """Signature activities as colored bars."""
    bw = w / len(vals)
    for k, v in enumerate(vals):
        bh = h * v / max(vals)
        ax.add_patch(Rectangle((x + k * bw + bw * 0.15, y + h - bh), bw * 0.7, bh,
                               fc=BAR_COLS[k % len(BAR_COLS)], ec="none", zorder=3))
    ax.plot([x, x + w], [y + h, y + h], color=INK, lw=0.4, zorder=3)


def fair_letter(ax, x, y, letter, s=2.9):
    rbox(ax, x, y - s / 2, s, s, FAIR_COLS[letter], "none", lw=0, r=0.5, z=4)
    txt(ax, x + s / 2, y + 0.05, letter, size=T_TINY, weight="bold", color="white", z=5)


# ---------------------------------------------------------------------------
# Data shown in the small pictures
# ---------------------------------------------------------------------------
def standard_order(index):
    """COSMIC files list contexts by 5' base; plots group them by mutation type."""
    return sorted(range(len(index)), key=lambda i: (index[i][2:5], index[i][0], index[i][6]))


def load_pictures(args):
    cosmic = pd.read_csv(args.cosmic, sep="\t", index_col=0)
    order = standard_order(list(cosmic.index))

    real = pd.read_csv(args.real_catalog, sep="\t", index_col=0, comment="#")
    col = next(c for c in real.columns if c.startswith(args.real_type))
    real_spec = real[col].reindex(cosmic.index).fillna(0).values[order]

    tum = pd.read_csv(os.path.join(args.tissue_dir, f"signature_contributions_WGS_{args.sim_type}.dat"),
                      sep="\t", index_col=0)
    row = tum.loc[tum.sum(axis=1).sort_values().index[len(tum) // 2]]
    row = row[[s for s in row.index if s in cosmic.columns and row[s] > 0]]
    w = (row / row.sum()).sort_values(ascending=False)
    p = (cosmic[w.index] * w).sum(axis=1)
    p = (p / p.sum()).values
    counts = np.bincount(np.random.default_rng(2).choice(96, p=p, size=args.sim_n), minlength=96)

    # thumbnails use only the tools shown in the figure
    shown = {RESULT_NAMES.get(t, t) for t in args.tools}
    metrics = pd.read_csv(os.path.join(args.accuracy_run, "all_evaluation_metrics.tsv"), sep="\t")
    metrics = metrics[metrics.tool.isin(shown)]
    tae = metrics.groupby(["num_muts", "tool"]).TAE.mean().unstack()

    agree = pd.read_csv(os.path.join(args.agreement_run, "all_tool_agreement.tsv"), sep="\t")
    agree = agree[agree.tool_1.isin(shown) & agree.tool_2.isin(shown)]
    pair = agree.groupby(["tool_1", "tool_2"]).mean_cosine_similarity.mean()
    names = sorted(set(agree.tool_1) | set(agree.tool_2))
    mat = pd.DataFrame(1.0, index=names, columns=names)
    for (a, b), v in pair.items():
        mat.loc[a, b] = mat.loc[b, a] = v
    # a synthetic sample: the first mix in the user's table of signature weights
    mix = pd.read_csv(args.custom, sep="\t", index_col=0).iloc[0]
    mix = list(mix[mix > 0].values)
    return real_spec, list(w.values[:5]), mix, p[order], counts[order], tae, mat


# ---------------------------------------------------------------------------
def draw(fig, ax, args):
    real_spec, acts, mix, expected, counts, tae, agree = load_pictures(args)
    new = [t for t in args.tools if t in args.new_methods]

    # ---- central band: the workflow ---------------------------------------
    by0, by1 = 38.0, 74.0
    rbox(ax, 2.0, by0, W - 4.0, by1 - by0, BAND, "none", lw=0, r=3.0, z=0)
    txt(ax, 5.0, by0 + 4.0, "SigFitBench", size=T_TITLE, weight="bold", ha="left", z=7)
    txt(ax, 24.0, by0 + 4.0, "one Snakemake workflow", size=T_SMALL, color=MUTED, ha="left",
        z=7)

    yc, dy = by0 + 20.0, 8.5
    cfg = pill(ax, 13.0, yc, 15.0, "config.yaml", INK, outline=True)
    sim = pill(ax, 37.0, yc - dy, 17.0, "Simulate", SIM, n=1, sub="known truth")
    real = pill(ax, 37.0, yc + dy, 17.0, "Real data", REAL, n=2, sub="no truth")
    ref = pill(ax, 61.0, yc, 15.0, "COSMIC", REF, sub="reference")
    fit = pill(ax, 84.0, yc, 12.0, "Fit", FIT, n=3, sub="each sample alone")
    score = pill(ax, 109.0, yc - dy, 16.0, "Score", SIM, n=4, sub="vs ground truth")
    comp = pill(ax, 109.0, yc + dy, 16.0, "Compare", REAL, n=4, sub="between tools")
    rep = pill(ax, 136.0, yc, 15.0, "Report", REPORT, n=5, sub="HTML + TSV")
    live = pill(ax, 163.0, yc, 18.0, "Live HTML", LIVE, n=5, sub="fixed dataset")

    link(ax, (cfg[1], yc), (sim[0], yc - dy), SIM)
    link(ax, (cfg[1], yc), (real[0], yc + dy), REAL)
    link(ax, (sim[1], yc - dy), (ref[0], yc - 1.2), SIM)
    link(ax, (real[1], yc + dy), (ref[0], yc + 1.2), REAL)
    link(ax, (ref[1], yc), (fit[0], yc))
    link(ax, (fit[1], yc - 1.2), (score[0], yc - dy), SIM)
    link(ax, (fit[1], yc + 1.2), (comp[0], yc + dy), REAL)
    link(ax, (score[1], yc - dy), (rep[0], yc - 1.2), SIM)
    link(ax, (comp[1], yc + dy), (rep[0], yc + 1.2), REAL)
    # the live report is built separately, from one fixed simulated dataset
    # dashed line with a solid head, so the head is drawn whole
    ax.plot([rep[1], live[0] - 1.5], [yc, yc], color=LIVE, lw=0.8, ls=(0, (3, 2)), zorder=5)
    link(ax, (live[0] - 1.5, yc), (live[0], yc), LIVE)
    # the ground truth goes straight to scoring
    # drawn with straight lines: up from Simulate, across, down into Score
    gy = yc - dy - 5.5
    ax.plot([37.0, 37.0, 109.0, 109.0], [yc - dy - 3.0, gy, gy, gy + 1.0], color=SIM, lw=0.8,
            ls=(0, (3, 2)), zorder=5)
    link(ax, (109.0, gy + 1.0), (109.0, yc - dy - 3.0), SIM)
    txt(ax, 73.0, by0 + 3.0, "ground truth, never shown to the tools", size=T_SMALL,
        color=SIM, style="italic", z=7)

    # ---- top panels -------------------------------------------------------
    ty0, th = 2.0, 33.0

    # 1 simulation
    x0, w0 = 2.0, 73.0
    panel(ax, x0, ty0, w0, th, SIM, 1, "Simulations")
    # two ways to choose the activities, both sampled the same way
    bw, bh = 8.5, 6.0
    rows = [("Realistic", "real tumor combinations", acts), ("Synthetic", "user-defined signature combinations", mix)]
    row_y = [ty0 + 11.5, ty0 + 20.0]
    # activities × COSMIC signatures, then N mutations are sampled
    sx, sw, sh = x0 + 52.0, 15.0, 8.0
    sy = (row_y[0] + row_y[1]) / 2 - sh / 2
    for (name, note, vals), y in zip(rows, row_y):
        bars(ax, x0 + 4.5, y - bh / 2, bw, bh, vals)
        txt(ax, x0 + 15.0, y - 1.5, name, size=T_TEXT, weight="bold", color=SIM, ha="left")
        txt(ax, x0 + 15.0, y + 1.6, note, size=T_TINY, color=MUTED, ha="left")
        link(ax, (x0 + 28.0, y - 1.5), (sx - 0.8, sy + sh / 2 + (y - sy - sh / 2) * 0.25), SIM,
             lw=0.7)
    txt(ax, sx + sw / 2, sy - 3.0, "× COSMIC, sample", size=T_TINY, color=SIM)
    spectrum(ax, sx, sy, sw, sh, counts)
    txt(ax, sx + sw / 2, sy + sh + 1.4, "N mutations", size=T_SMALL, color=MUTED, va="top")
    txt(ax, x0 + 4.5, ty0 + 26.6, "Ground truth: the known activities", size=T_TEXT,
        color=SIM, weight="bold", ha="left")
    xx = x0 + 4.5
    for s in ["cancer type", "N mutations", "N samples", "replicates per sample"]:
        xx = chip(ax, xx, ty0 + 30.4, s, SIM) + 0.8

    # FAIR
    x1, w1 = 77.5, 46.0
    panel(ax, x1, ty0, w1, th, NEUTRAL, None, "FAIR by design")
    checks = [("Open code, MIT license", "FA"), ("One conda env per tool", "AR"),
              ("One config, fixed seeds", "R"), ("VCF, TSV, COSMIC formats", "I"),
              ("Laptop or SLURM cluster", "A"), ("Easy to add a new tool", "R"),
              ("Live benchmark, fixed dataset", "FR")]
    for i, (s, letters) in enumerate(checks):
        y = ty0 + 8.6 + i * 3.6
        txt(ax, x1 + 3.6, y, s, size=T_SMALL, ha="left")
        for k, L in enumerate(reversed(letters)):
            fair_letter(ax, x1 + w1 - 4.4 - k * 3.4, y, L)

    # 3 tools
    x2 = 126.0
    w2 = W - 2.0 - x2
    panel(ax, x2, ty0, w2, th, FIT, 3, f"{len(args.tools)} fitting tools")
    txt(ax, x2 + 3.6, ty0 + 8.4, "each sample fitted on its own, same signatures", size=T_TINY,
        color=MUTED,
        ha="left")
    rows = (len(args.tools) + 1) // 2
    for i, t in enumerate(args.tools):
        c, r = divmod(i, rows)
        is_new = t in new
        txt(ax, x2 + 3.6 + c * 26.0, ty0 + 12.2 + r * 3.1, t, size=T_SMALL, ha="left",
            color=NEW if is_new else INK, weight="bold" if is_new else "normal")

    # ---- bottom panels ----------------------------------------------------
    py0, ph = 77.0, 33.0

    # 2 real data
    x3, w3 = 2.0, 57.0
    panel(ax, x3, py0, w3, ph, REAL, 2, "Real data")
    spectrum(ax, x3 + 4.5, py0 + 8.0, 21, 10, real_spec)
    for k in range(6):
        ax.add_patch(Rectangle((x3 + 4.5 + k * 21 / 6 + 0.1, py0 + 18.6), 21 / 6 - 0.2, 0.7,
                               fc=SBS_COL[k], ec="none", zorder=3))
    txt(ax, x3 + 15.0, py0 + 20.6, f"e.g. PCAWG\n{args.real_type}", size=T_TINY, color=MUTED,
        style="italic", va="top")
    txt(ax, x3 + 29.0, py0 + 8.0, "Real WGS or\nWES catalog", size=T_TEXT,
        weight="bold", ha="left", va="top")
    txt(ax, x3 + 29.0, py0 + 16.4, "formats:", size=T_TINY, color=MUTED, ha="left")
    for i, s in enumerate(["VCF + genome", "SBS count matrix"]):
        y = py0 + 19.8 + i * 3.4
        ax.add_patch(Circle((x3 + 29.9, y), 0.5, fc=REAL, ec="none", zorder=4))
        txt(ax, x3 + 31.4, y, s, size=T_SMALL, ha="left")
    # worded like the matching line in panel 1
    txt(ax, x3 + 4.5, py0 + ph - 3.8, "Ground truth: unknown", size=T_TEXT, color=REAL,
        weight="bold", ha="left")

    # 4 metrics
    x4, w4 = 61.5, 57.0
    panel(ax, x4, py0, w4, ph, METRIC, 4, "Metrics")
    groups = [("Simulated, vs ground truth", SIM, ["Estimation Error", "precision", "recall", "F1", "MCC"]),
              ("Real, between tools", REAL, ["cosine similarity", "shared signatures"]),
              ("Every run", NEUTRAL, ["runtime", "memory"])]
    y = py0 + 8.6
    for title, colr, items in groups:
        txt(ax, x4 + 3.6, y, title, size=T_SMALL, weight="bold", color=colr, ha="left")
        xx = x4 + 3.6
        for it in items:
            xx = chip(ax, xx, y + 3.6, it, colr) + 0.8
        y += 8.2

    # 5 report, with two thumbnails from results on disk
    x5 = 121.0
    w5 = W - 2.0 - x5
    panel(ax, x5, py0, w5, ph, REPORT, 5, "HTML report")

    def mm_axes(x_, y_, w_, h_):
        return fig.add_axes([x_ / W, 1 - (y_ + h_) / H, w_ / W, h_ / H])

    a1 = mm_axes(x5 + 8.0, py0 + 9.5, 20.0, 12.5)
    for k, tool in enumerate(tae.columns):
        a1.plot(tae.index, tae[tool], lw=0.7, marker="o", ms=1.1,
                color=THUMB_COLS[k % len(THUMB_COLS)])
    a1.set_xscale("log")
    a1.set_xticks([])
    a1.set_yticks([])
    a1.minorticks_off()
    a1.patch.set_alpha(0)
    for sp in ("top", "right"):
        a1.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        a1.spines[sp].set_linewidth(0.5)
    a1.set_xlabel("mutations", fontsize=T_TINY, labelpad=0.8)
    a1.set_ylabel("error", fontsize=T_TINY, labelpad=0.8)
    a1.set_title("accuracy", fontsize=T_SMALL, pad=1.5, color=SIM)

    a2 = mm_axes(x5 + 37.0, py0 + 9.5, 12.5, 12.5)
    a2.imshow(agree.values, cmap=mpl.colors.LinearSegmentedColormap.from_list(
        "r", ["#FFFFFF", REAL]), vmin=0, vmax=1)
    a2.set_xticks([])
    a2.set_yticks([])
    for sp in a2.spines.values():
        sp.set_visible(False)
    a2.set_title("agreement", fontsize=T_SMALL, pad=1.5, color=REAL)
    txt(ax, x5 + w5 / 2, py0 + ph - 5.6, "TSV tables for each sample", size=T_SMALL,
        color=MUTED)
    txt(ax, x5 + w5 / 2, py0 + ph - 2.6, "Live benchmark on a fixed dataset",
        size=T_SMALL, color=LIVE, weight="bold")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tissue-dir", default="cosmic_tissue_data")
    ap.add_argument("--cosmic", default="input/COSMIC_v3.3.1_SBS_GRCh38.txt")
    ap.add_argument("--real-catalog",
                    default="real_mutational_catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat")
    ap.add_argument("--real-type", default="Skin-Melanoma")
    ap.add_argument("--sim-type", default="Breast")
    ap.add_argument("--sim-n", type=int, default=500)
    ap.add_argument("--custom", default="config/custom.tsv",
                    help="Table of signature mixes; its first row is the synthetic example")
    ap.add_argument("--accuracy-run", default="results_unknown_0_all_tools",
                    help="Results folder of a simulated run, for the accuracy thumbnail")
    ap.add_argument("--agreement-run", default="results_wgs_pcawg",
                    help="Results folder of a real-data run, for the agreement thumbnail")
    ap.add_argument("--tools", nargs="+", default=TOOLS, help="Tools listed, in order")
    ap.add_argument("--new-methods", nargs="*", default=NEW_METHODS,
                    help="Tools shown in red as new methods")
    ap.add_argument("--out", default="sigfitbench_workflow_figure.pdf")
    ap.add_argument("--png", default=None, help="Optional PNG preview path")
    args = ap.parse_args()

    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    draw(fig, ax, args)

    fig.savefig(args.out, dpi=600)
    if args.png:
        fig.savefig(args.png, dpi=220)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
