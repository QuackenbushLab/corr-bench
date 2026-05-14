#!/usr/bin/env python3
"""
Figure 2 (variant 3): Side-by-side stacked bars showing WGCNA module size
distributions for CorrAdjust vs COBRA MLE across 11 GTEx tissues.

Each tissue gets two columns: left=CorrAdjust, right=COBRA MLE.
Each segment = one module (colored), grey = unassigned (module 0).

Output:
    data/figures/fig2_module_stacked.png / .pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CA_DIR = os.path.join(BASE, "results", "tpmfilter_data", "wgcna_ora")
COBRA_DIR = os.path.join(BASE, "results", "wgcna_ora_cov_all")
HARDY_DIR = os.path.join(BASE, "results", "wgcna_ora_hardy")
FIG_DIR = os.path.join(BASE, "figures")

HARDY_TISSUES = {"GTEx-Muscle-Skeletal", "GTEx-Adipose-Subcutaneous", "GTEx-Skin-Sun-Exposed-Lower-leg"}

TISSUES = [
    "GTEx-Adipose-Subcutaneous",
    "GTEx-Artery-Tibial",
    "GTEx-Heart-Left-Ventricle",
    "GTEx-Lung",
    "GTEx-Muscle-Skeletal",
    "GTEx-Nerve-Tibial",
    "GTEx-Skin-Sun-Exposed-Lower-leg",
    "GTEx-Small-Intestine-Terminal-Ileum",
    "GTEx-Spleen",
    "GTEx-Thyroid",
    "GTEx-Whole-Blood",
]

SHORT = {
    "GTEx-Adipose-Subcutaneous": "Adipose",
    "GTEx-Artery-Tibial": "Artery",
    "GTEx-Heart-Left-Ventricle": "Heart",
    "GTEx-Lung": "Lung",
    "GTEx-Muscle-Skeletal": "Muscle",
    "GTEx-Nerve-Tibial": "Nerve",
    "GTEx-Skin-Sun-Exposed-Lower-leg": "Skin",
    "GTEx-Small-Intestine-Terminal-Ileum": "SmallInt.",
    "GTEx-Spleen": "Spleen",
    "GTEx-Thyroid": "Thyroid",
    "GTEx-Whole-Blood": "WholeBlood",
}


def get_module_sizes(path):
    """Return dict {module_id: gene_count} sorted by size descending (0=unassigned last)."""
    df = pd.read_csv(path, sep="\t")
    counts = df["module"].value_counts().to_dict()
    unassigned = counts.pop(0, 0)
    # Sort assigned modules by size descending
    assigned = sorted(counts.items(), key=lambda x: -x[1])
    return assigned, unassigned


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    # Gather data
    data = []
    for tissue in TISSUES:
        ca_path = os.path.join(CA_DIR, tissue, "CorrAdjust", "module_assignments.tsv")
        ca_assigned, ca_unassigned = get_module_sizes(ca_path)

        # For Hardy tissues, use Hardy MLE instead of COBRA MLE
        if tissue in HARDY_TISSUES:
            right_path = os.path.join(HARDY_DIR, tissue, "MLE", "module_assignments.tsv")
            right_label = "Hardy"
        else:
            right_path = os.path.join(COBRA_DIR, tissue, "MLE", "module_assignments.tsv")
            right_label = "COBRA"
        rt_assigned, rt_unassigned = get_module_sizes(right_path)

        data.append({
            "tissue": SHORT[tissue],
            "ca_assigned": ca_assigned, "ca_unassigned": ca_unassigned,
            "rt_assigned": rt_assigned, "rt_unassigned": rt_unassigned,
            "right_label": right_label,
        })
        print(f"  {SHORT[tissue]:12s}  CA: {len(ca_assigned):2d} modules  "
              f"{right_label}: {len(rt_assigned):2d} modules")

    # Color palette: use tab20 for up to 20 modules
    max_modules = max(
        max(len(d["ca_assigned"]) for d in data),
        max(len(d["rt_assigned"]) for d in data),
    )
    cmap = plt.cm.tab20
    colors = [cmap(i / max(max_modules - 1, 1)) for i in range(max_modules)]

    n = len(data)
    bar_width = 0.35
    gap = 0.15
    x_positions = np.arange(n) * (2 * bar_width + gap + 0.15)

    fig, ax = plt.subplots(figsize=(14, 5.5))

    for i, d in enumerate(data):
        x_ca = x_positions[i] - bar_width / 2 - 0.02
        x_rt = x_positions[i] + bar_width / 2 + 0.02
        is_hardy = d["right_label"] == "Hardy"

        # CorrAdjust
        bottom = 0
        for j, (mod_id, count) in enumerate(d["ca_assigned"]):
            ax.bar(x_ca, count, bar_width, bottom=bottom, color=colors[j],
                   edgecolor="white", linewidth=0.3)
            bottom += count
        if d["ca_unassigned"] > 0:
            ax.bar(x_ca, d["ca_unassigned"], bar_width, bottom=bottom,
                   color="#d9d9d9", edgecolor="white", linewidth=0.3)

        # Right bar: COBRA MLE or Hardy MLE
        bottom = 0
        for j, (mod_id, count) in enumerate(d["rt_assigned"]):
            ax.bar(x_rt, count, bar_width, bottom=bottom, color=colors[j],
                   edgecolor="white", linewidth=0.3)
            bottom += count
        if d["rt_unassigned"] > 0:
            ax.bar(x_rt, d["rt_unassigned"], bar_width, bottom=bottom,
                   color="#d9d9d9", edgecolor="white", linewidth=0.3)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([d["tissue"] for d in data], rotation=45, ha="right",
                       fontsize=16)
    ax.set_ylabel("Number of genes", fontsize=16)
    ax.tick_params(axis="y", labelsize=14)

    ax.set_xlim(x_positions[0] - 0.5, x_positions[-1] + 0.5)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()

    for ext in ["png", "pdf"]:
        out = os.path.join(FIG_DIR, f"fig2_module_stacked.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
