#!/usr/bin/env python3
"""
Figure 3 combined:
  Panel A – best-module tissue-relevance precision (existing fig3)
  Panel B – runtime boxplot across tissues: CorrAdjust vs COBRA MLE
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.image import imread

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# works from src/ or from project root (flat copy on AWS)
BASE = (os.path.dirname(_SCRIPT_DIR)
        if os.path.basename(_SCRIPT_DIR) == "src"
        else _SCRIPT_DIR)
FIG_DIR = os.path.join(BASE, "figures")
SUMMARY = os.path.join(BASE, "results", "benchmark_medcount_summary.tsv")

def main():
    # ── load runtime data ──
    df = pd.read_csv(SUMMARY, sep="\t")
    ca = df[df["method"] == "CorrAdjust"]["elapsed_s"].values / 3600  # hours
    cobra = df[df["method"] == "COBRA_deprecated"]["elapsed_s"].values / 3600

    # ── load panel A image ──
    panel_a = imread(os.path.join(FIG_DIR, "fig_3bar_precision.png"))

    # ── build combined figure (side by side) ──
    fig = plt.figure(figsize=(22, 7))
    # Panel A on left (image)
    ax_a = fig.add_axes([0.01, 0.02, 0.62, 0.95])
    ax_a.imshow(panel_a)
    ax_a.axis("off")
    ax_a.text(0.01, 0.97, "A", transform=ax_a.transAxes,
              fontsize=18, fontweight="bold", va="top")

    # Panel B on right (boxplot)
    ax_b = fig.add_axes([0.70, 0.12, 0.27, 0.78])

    bp = ax_b.boxplot(
        [cobra, ca],
        tick_labels=["COBRA MLE", "CorrAdjust"],
        widths=0.4,
        patch_artist=True,
        showmeans=True,
        meanprops=dict(marker="D", markerfacecolor="black", markersize=6),
    )
    colors = ["#d62728", "#1f77b4"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    for median in bp["medians"]:
        median.set_color("black")
        median.set_linewidth(1.5)

    # jitter individual tissue points
    for i, (data, color) in enumerate(zip([cobra, ca], colors)):
        jitter = np.random.default_rng(42).uniform(-0.08, 0.08, len(data))
        ax_b.scatter(np.full_like(data, i + 1) + jitter, data,
                     color=color, edgecolor="white", s=50, zorder=5, alpha=0.85)

    ax_b.set_ylabel("Runtime (hours)", fontsize=16)
    ax_b.tick_params(axis="x", labelsize=14)
    ax_b.tick_params(axis="y", labelsize=14)
    ax_b.grid(axis="y", alpha=0.3)
    ax_b.text(-0.05, 1.05, "B", transform=ax_b.transAxes,
              fontsize=18, fontweight="bold", va="top")

    out = os.path.join(FIG_DIR, "fig3_combined.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
