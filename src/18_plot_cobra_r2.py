#!/usr/bin/env python3
"""18_plot_cobra_r2.py — Plot COBRA train/test eigenvalue prediction R² per tissue."""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

df = pd.read_csv(os.path.join(RESULTS_DIR, "cobra_train_test_r2.tsv"), sep="\t")

# Short tissue names
df["tissue"] = df["dataset"].str.replace("GTEx-", "")

tissues = df["tissue"].unique()
modes = ["MLE", "nnls", "nnlasso"]
mode_labels = {"MLE": "COBRA MLE", "nnls": "COBRA NNLS", "nnlasso": "COBRA NNLASSO"}
colors = {"MLE": "#1f77b4", "nnls": "#ff7f0e", "nnlasso": "#2ca02c"}

x = np.arange(len(tissues))
width = 0.25

# Plot 1: R² diagonal (actual 1 - SS_res/SS_tot)
fig, ax = plt.subplots(figsize=(14, 5))

for i, mode in enumerate(modes):
    mask = df["cobra_mode"] == mode
    vals = df.loc[mask, "corr_diag"].astype(float).values
    ax.bar(x + i * width, vals, width, label=mode_labels[mode], color=colors[mode])

ax.set_ylabel("Pearson's r", fontsize=16)
ax.set_xlabel("Tissue", fontsize=16)
ax.set_xticks(x + width)
ax.set_xticklabels(tissues, rotation=45, ha="right", fontsize=14)
ax.tick_params(axis="y", labelsize=14)
ax.legend(fontsize=14, loc="upper left", bbox_to_anchor=(1.01, 1.0),
          borderaxespad=0., frameon=False)
ax.axhline(y=0, color="gray", linewidth=0.5)

plt.tight_layout()
out = os.path.join(FIGURES_DIR, "fig4_cobra_train_test_r2.png")
plt.savefig(out, dpi=150)
print(f"Saved: {out}")
plt.close()

# Also plot corr_mat (full matrix Pearson)
fig, ax = plt.subplots(figsize=(14, 5))

for i, mode in enumerate(modes):
    mask = df["cobra_mode"] == mode
    vals = df.loc[mask, "corr_matrix"].astype(float).values
    ax.bar(x + i * width, vals, width, label=mode_labels[mode], color=colors[mode])

ax.set_ylabel("Pearson's r", fontsize=16)
ax.set_xlabel("Tissue", fontsize=16)
ax.set_xticks(x + width)
ax.set_xticklabels(tissues, rotation=45, ha="right", fontsize=14)
ax.tick_params(axis="y", labelsize=14)
ax.legend(fontsize=14, loc="upper left", bbox_to_anchor=(1.01, 1.0),
          borderaxespad=0., frameon=False)
ax.axhline(y=0, color="gray", linewidth=0.5)

plt.tight_layout()
out2 = os.path.join(FIGURES_DIR, "fig4_cobra_train_test_r2.png")
plt.savefig(out2, dpi=150)
print(f"Saved: {out2}")
plt.close()
