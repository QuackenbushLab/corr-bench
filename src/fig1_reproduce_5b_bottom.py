#!/usr/bin/env python3
"""fig1_reproduce_5b_bottom.py — Figure 1: Reproduction of CorrAdjust Fig 5B bottom panel.

Horizontal dot plot showing Global Enrichment Scores (GES) for:
  - Raw (uncorrected)
  - CorrAdjust (corrected)
  - COBRA MLE (deprecated mode, corr)
per GTEx tissue.  Right-side annotations: #PCs selected and n_train.
Uses Thyroid override (4 PCs) instead of auto (0 PCs).
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

# ── Data ──────────────────────────────────────────────────────────
# CorrAdjust results (from benchmark_gtex_full.tsv)
# Thyroid uses override row (4 PCs)
corradjust_data = {
    "GTEx-Lung":                          {"raw": 1.331897, "ca": 2.458416, "pcs": 7, "n_train": 289, "n_covs": 8},
    "GTEx-Adipose-Subcutaneous":          {"raw": 1.866807, "ca": 2.160726, "pcs": 2, "n_train": 332, "n_covs": 8},
    "GTEx-Thyroid":                       {"raw": 1.715600, "ca": 2.084186, "pcs": 4, "n_train": 327, "n_covs": 8},
    "GTEx-Small-Intestine-Terminal-Ileum":{"raw": 1.270549, "ca": 2.019464, "pcs": 2, "n_train": 94, "n_covs": 8},
    "GTEx-Spleen":                        {"raw": 1.944790, "ca": 2.010413, "pcs": 1, "n_train": 121, "n_covs": 8},
    "GTEx-Muscle-Skeletal":               {"raw": 1.444600, "ca": 1.997370, "pcs": 2, "n_train": 402, "n_covs": 8},
    "GTEx-Heart-Left-Ventricle":          {"raw": 0.917951, "ca": 1.768638, "pcs": 6, "n_train": 216, "n_covs": 8},
    "GTEx-Artery-Tibial":                 {"raw": 1.667506, "ca": 1.667506, "pcs": 0, "n_train": 332, "n_covs": 8},
    "GTEx-Skin-Sun-Exposed-Lower-leg":    {"raw": 1.345542, "ca": 1.345542, "pcs": 0, "n_train": 351, "n_covs": 8},
    "GTEx-Nerve-Tibial":                  {"raw": 1.253402, "ca": 1.253402, "pcs": 0, "n_train": 310, "n_covs": 8},
    "GTEx-Whole-Blood":                   {"raw": 0.402540, "ca": 1.177090, "pcs": 3, "n_train": 378, "n_covs": 8},
}

# COBRA MLE (deprecated) scores (from cobra_3modes_summary.tsv)
cobra_mle = {
    "GTEx-Adipose-Subcutaneous":          1.419211,
    "GTEx-Artery-Tibial":                 0.904588,
    "GTEx-Heart-Left-Ventricle":          0.588835,
    "GTEx-Lung":                          0.046956,
    "GTEx-Muscle-Skeletal":               0.006333,
    "GTEx-Nerve-Tibial":                  0.864968,
    "GTEx-Skin-Sun-Exposed-Lower-leg":    0.917739,
    "GTEx-Small-Intestine-Terminal-Ileum": 1.175685,
    "GTEx-Spleen":                        0.009413,
    "GTEx-Thyroid":                       0.968359,
    "GTEx-Whole-Blood":                   0.318948,
}

# Build dataframe
rows = []
for tissue, d in corradjust_data.items():
    rows.append({
        "tissue": tissue,
        "short": tissue.replace("GTEx-", ""),
        "raw": d["raw"],
        "corradjust": d["ca"],
        "cobra": cobra_mle.get(tissue, 0.0),
        "n_pcs": d["pcs"],
        "n_covs": d["n_covs"],
        "n_train": d["n_train"],
    })

df = pd.DataFrame(rows)
# Sort by CorrAdjust score descending (matching paper order)
df = df.sort_values("corradjust", ascending=False).reset_index(drop=True)

# ── Plot ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5.5))

n = len(df)
y = np.arange(n)

# Raw — grey open circles
ax.scatter(df["raw"], y, marker="o", facecolors="none", edgecolors="#888888",
           s=60, linewidths=1.2, zorder=3, label="Raw")

# CorrAdjust — green stars
ax.scatter(df["corradjust"], y, marker="*", color="#2ca02c",
           s=140, zorder=4, label="CorrAdjust")

# COBRA MLE — blue diamonds
ax.scatter(df["cobra"], y, marker="D", color="#1f77b4",
           s=50, zorder=4, label="COBRA")

# Connecting lines (raw → corradjust)
for i in range(n):
    xmin = min(df.loc[i, "raw"], df.loc[i, "corradjust"], df.loc[i, "cobra"])
    xmax = max(df.loc[i, "raw"], df.loc[i, "corradjust"], df.loc[i, "cobra"])
    ax.plot([xmin, xmax], [i, i], color="#dddddd", lw=0.7, zorder=1)

ax.set_yticks(y)
ax.set_yticklabels(df["short"], fontsize=16)
ax.invert_yaxis()
ax.set_xlabel("Global Enrichment Score", fontsize=16)
ax.tick_params(axis="x", labelsize=14)
ax.grid(axis="x", alpha=0.2, lw=0.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(-0.1, None)

# Right-side annotations: PCs, n_covariates, n_train
xmax_lim = ax.get_xlim()[1]
for i in range(n):
    ax.annotate(
        f"  PCs={int(df.loc[i, 'n_pcs']):>2d}   covs={int(df.loc[i, 'n_covs']):>1d}   n={int(df.loc[i, 'n_train']):>4d}",
        xy=(xmax_lim, i), fontsize=7, va="center", family="monospace",
        annotation_clip=False)

# Legend
ax.legend(loc="lower right", fontsize=9, framealpha=0.9, edgecolor="lightgrey")

fig.tight_layout()
os.makedirs("figures", exist_ok=True)
out = "figures/fig1_ges_gtex.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"Saved: {out}")

out_pdf = "figures/fig1_ges_gtex.pdf"
fig.savefig(out_pdf, bbox_inches="tight")
print(f"Saved: {out_pdf}")
plt.close()
