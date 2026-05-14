#!/usr/bin/env python3
"""
Combined summary dotplot: CorrAdjust (left) and COBRA MLE (right) side by side.
Top 6 GO BP terms per tissue-module (FDR < 5%).
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COBRA_DIR = os.path.join(BASE, "results", "wgcna_ora_cov_all")
CA_DIR = os.path.join(BASE, "results", "tpmfilter_data", "wgcna_ora")
FIG_DIR = os.path.join(BASE, "figures")

PADJ_CUTOFF = 0.05
TOP_K = 6

COBRA_PANELS = [
    ("GTEx-Adipose-Subcutaneous", 3),
    ("GTEx-Artery-Tibial", 7),
    ("GTEx-Heart-Left-Ventricle", 4),
    ("GTEx-Lung", 3),
    ("GTEx-Muscle-Skeletal", 5),
    ("GTEx-Nerve-Tibial", 2),
    ("GTEx-Small-Intestine-Terminal-Ileum", 1),
    ("GTEx-Thyroid", 1),
    ("GTEx-Whole-Blood", 2),
]

CA_PANELS = [
    ("GTEx-Adipose-Subcutaneous", 1),
    ("GTEx-Heart-Left-Ventricle", 5),
    ("GTEx-Muscle-Skeletal", 2),
    ("GTEx-Nerve-Tibial", 7),
    ("GTEx-Small-Intestine-Terminal-Ileum", 7),
    ("GTEx-Thyroid", 4),
    ("GTEx-Whole-Blood", 2),
]


def clean_term(t):
    for prefix in ("GOBP_", "GOCC_", "GOMF_", "KEGG_", "REACTOME_",
                    "WP_", "BIOCARTA_", "PID_", "NABA_",
                    "KEGG_MEDICUS_REFERENCE_", "SA_", "SIG_"):
        if t.startswith(prefix):
            t = t[len(prefix):]
    t = t.replace("_", " ").title()
    if len(t) > 50:
        t = t[:47] + "..."
    return t


def load_module_terms(base_dir, method_sub, tissue, module):
    path = os.path.join(base_dir, tissue, method_sub, "ora_gobp.tsv")
    if not os.path.exists(path):
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    sig = df[(df["module"] == module) & (df["FDR"] < PADJ_CUTOFF)].copy()
    if sig.empty:
        return pd.DataFrame()
    sig = sig.sort_values(["FDR", "overlap_count"], ascending=[True, False]).head(TOP_K)
    sig["neg_log_fdr"] = sig["FDR"].clip(lower=1e-300).apply(lambda x: -np.log10(x))
    sig["clean_term"] = sig["Term"].apply(clean_term)
    return sig


def build_panel_data(base_dir, method_sub, panels):
    x_labels = []
    tissue_data = []
    for tissue, module in panels:
        short = tissue.replace("GTEx-", "")
        x_labels.append(f"{short}\n(M{module})")
        df = load_module_terms(base_dir, method_sub, tissue, module)
        tissue_data.append(df)
    return x_labels, tissue_data


def plot_dotplot(ax, x_labels, tissue_data, title, norm, max_ov):
    # Build ordered term list
    term_order = []
    for df in tissue_data:
        if df is not None and not df.empty:
            for t in df["clean_term"]:
                if t not in term_order:
                    term_order.append(t)

    term_idx = {t: i for i, t in enumerate(term_order)}
    n_tissues = len(x_labels)
    n_terms = len(term_order)

    for xi, df in enumerate(tissue_data):
        if df is None or df.empty:
            continue
        for _, r in df.iterrows():
            y = term_idx[r["clean_term"]]
            sz = 30 + 200 * (r["overlap_count"] / max_ov)
            c = plt.cm.YlOrRd(norm(r["neg_log_fdr"]))
            ax.scatter(xi, y, s=sz, c=[c], edgecolors="grey", linewidths=0.5, zorder=3)

    ax.set_xticks(range(n_tissues))
    ax.set_xticklabels(x_labels, fontsize=7, ha="center")
    ax.set_yticks(range(n_terms))
    ax.set_yticklabels(term_order, fontsize=6.5)
    ax.set_xlim(-0.5, n_tissues - 0.5)
    ax.set_ylim(-0.5, n_terms - 0.5)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.15, zorder=0)
    ax.set_title(title, fontsize=10, fontweight="bold")

    return n_terms


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    # Build data for both panels
    cobra_labels, cobra_data = build_panel_data(COBRA_DIR, "MLE", COBRA_PANELS)
    ca_labels, ca_data = build_panel_data(CA_DIR, "CorrAdjust", CA_PANELS)

    # Global normalization across both panels
    all_nlf = []
    all_ov = []
    for dlist in [cobra_data, ca_data]:
        for df in dlist:
            if df is not None and not df.empty:
                all_nlf.extend(df["neg_log_fdr"].tolist())
                all_ov.extend(df["overlap_count"].tolist())
    max_nlf = max(all_nlf) if all_nlf else 1
    max_ov = max(all_ov) if all_ov else 1
    norm = Normalize(vmin=0, vmax=max(max_nlf, 1))

    # Count terms for height estimation
    def count_terms(tissue_data):
        terms = []
        for df in tissue_data:
            if df is not None and not df.empty:
                for t in df["clean_term"]:
                    if t not in terms:
                        terms.append(t)
        return len(terms)

    n_cobra = count_terms(cobra_data)
    n_ca = count_terms(ca_data)
    max_terms = max(n_cobra, n_ca)
    fig_height = max(10, max_terms * 0.27 + 3)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22, fig_height),
                                    gridspec_kw={"wspace": 0.45})

    plot_dotplot(ax1, ca_labels, ca_data,
                 f"CorrAdjust — Top {TOP_K} GO BP terms per module (FDR < 5%)",
                 norm, max_ov)
    plot_dotplot(ax2, cobra_labels, cobra_data,
                 f"COBRA MLE — Top {TOP_K} GO BP terms per module (FDR < 5%)",
                 norm, max_ov)

    # Colorbar attached to right panel only (avoids stealing space between panels)
    sm = ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax2, shrink=0.25, pad=0.08)
    cbar.set_label("-log\u2081\u2080(FDR)", fontsize=9)

    # Size legend on right panel
    for frac, lab in [(0.25, str(int(max_ov * 0.25))),
                       (0.5, str(int(max_ov * 0.5))),
                       (1.0, str(int(max_ov)))]:
        ax2.scatter([], [], s=30 + 200 * frac, c="grey", alpha=0.5,
                    edgecolors="grey", label=f"{lab} genes")
    ax2.legend(title="Overlap", loc="lower right", fontsize=7, title_fontsize=8,
               framealpha=0.9)

    # Panel labels
    ax1.text(-0.02, 1.02, "A", transform=ax1.transAxes, fontsize=14,
             fontweight="bold", va="bottom", ha="right")
    ax2.text(-0.02, 1.02, "B", transform=ax2.transAxes, fontsize=14,
             fontweight="bold", va="bottom", ha="right")

    fig.subplots_adjust(left=0.12, right=0.92)

    for ext in ["png", "pdf"]:
        out = os.path.join(FIG_DIR, f"fig3_summary_gobp_combined.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
