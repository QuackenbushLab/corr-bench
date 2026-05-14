#!/usr/bin/env python3
"""
Generate supplementary figures S1 (ORA top-100) and S2 (GSEA prerank).

Each figure is a 2-panel dotplot (CorrAdjust | COBRA) for Canonical Pathways
across 11 GTEx tissues.

Thyroid CorrAdjust uses the 4-PC override model instead of auto (0 PCs).
COBRA uses the MLE mode (formerly "deprecated").

Outputs:
    data/figures/figs1_ora_cp.pdf / .png
    data/figures/figs2_gsea_cp.pdf / .png
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from collections import OrderedDict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# ── Paths ─────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ORA_DIR = os.path.join(BASE, "figures", "ora_top100_all")
GSEA_DIR = os.path.join(BASE, "results", "tpmfilter_data", "gsea_all")
CORRADJUST_DIR = os.path.join(
    BASE, "results", "tpmfilter_data", "corradjust_gtex_outputs_models"
)
GMT_CP = os.path.join(BASE, "data", "input_data", "GMT_files",
                      "c2.cp.v2023.2.Hs.symbols.gmt")
FIG_DIR = os.path.join(BASE, "data", "figures")

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

TOP_N = 100
PADJ_CUTOFF = 0.05
TOP_K_PLOT = 5
FDR_CUTOFF = 0.05


# ── Helpers ───────────────────────────────────────────────────────
def load_gmt(path):
    gene_sets = OrderedDict()
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            gene_sets[parts[0]] = set(parts[2:])
    return gene_sets


def run_ora(query_genes, gene_sets, background_genes):
    query = set(query_genes) & background_genes
    N = len(background_genes)
    n = len(query)
    rows = []
    for gs_name, gs_genes in gene_sets.items():
        K_full = gs_genes & background_genes
        K = len(K_full)
        if K < 5 or K > 2000:
            continue
        overlap = query & K_full
        k = len(overlap)
        if k == 0:
            continue
        a, b = k, n - k
        c, d = K - k, N - K - n + k
        if d < 0:
            continue
        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "term": gs_name, "overlap": k, "query_size": n,
            "term_size": K, "pvalue": pval,
            "genes": ",".join(sorted(overlap)),
        })
    if not rows:
        return pd.DataFrame(columns=["term", "overlap", "query_size",
                                      "term_size", "pvalue", "genes", "padj"])
    df = pd.DataFrame(rows).sort_values("pvalue").reset_index(drop=True)
    m = len(df)
    df["padj"] = df["pvalue"] * m / (np.arange(1, m + 1))
    df["padj"] = df["padj"].clip(upper=1.0)
    df["padj"] = df["padj"][::-1].cummin()[::-1].values
    return df


def clean_term(name):
    for prefix in ["REACTOME_", "KEGG_", "WP_", "BIOCARTA_", "PID_", "SA_",
                    "SIG_", "KEGG_MEDICUS_REFERENCE_"]:
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    name = name.replace("_", " ").title()
    if len(name) > 55:
        name = name[:52] + "..."
    return name


def short_tissue(name):
    return name.replace("GTEx-", "")


# ── ORA: load existing + re-run Thyroid override ─────────────────
def load_ora_results():
    """Return {method_label: [(tissue, df_enr), ...]}"""
    gene_sets = load_gmt(GMT_CP)
    results = {"CorrAdjust": [], "COBRA": []}

    for tissue in TISSUES:
        # -- CorrAdjust --
        if tissue == "GTEx-Thyroid":
            # Re-run with override feature scores
            override_path = os.path.join(
                CORRADJUST_DIR, "GTEx-Thyroid-override",
                "feature_scores_Clean_Canonical_Pathways.tsv",
            )
            df_feat = pd.read_csv(override_path, sep="\t")
            df_feat = df_feat.sort_values("enrichment", ascending=False)
            background = set(df_feat["feature_name"].tolist())
            top_genes = df_feat["feature_name"].head(TOP_N).tolist()
            df_enr = run_ora(top_genes, gene_sets, background)
            print(f"  ORA Thyroid CorrAdjust (override 4 PCs): "
                  f"{(df_enr['padj'] < PADJ_CUTOFF).sum()} sig terms")
        else:
            tsv = os.path.join(ORA_DIR, tissue, "CorrAdjust_cp.tsv")
            df_enr = pd.read_csv(tsv, sep="\t")
        results["CorrAdjust"].append((tissue, df_enr))

        # -- COBRA (was "deprecated", now just "COBRA") --
        tsv = os.path.join(ORA_DIR, tissue, "COBRA_deprecated_cp.tsv")
        df_enr = pd.read_csv(tsv, sep="\t")
        results["COBRA"].append((tissue, df_enr))

    return results


# ── GSEA: load existing + re-run Thyroid override ────────────────
def load_gsea_results():
    """Return {method_label: [(tissue, df_res), ...]}"""
    import gseapy as gp

    gene_sets = load_gmt(GMT_CP)
    # Convert sets to lists for gseapy
    gene_sets_list = {k: list(v) for k, v in gene_sets.items()}

    results = {"CorrAdjust": [], "COBRA": []}

    for tissue in TISSUES:
        # -- CorrAdjust --
        if tissue == "GTEx-Thyroid":
            # Re-run GSEA with override feature scores
            override_path = os.path.join(
                CORRADJUST_DIR, "GTEx-Thyroid-override",
                "feature_scores_Clean_Canonical_Pathways.tsv",
            )
            df_feat = pd.read_csv(override_path, sep="\t")
            rnk = (df_feat.set_index("feature_name")["enrichment"]
                   .dropna().sort_values(ascending=False))
            rnk = rnk[~rnk.index.duplicated(keep="first")]
            print(f"  GSEA Thyroid CorrAdjust (override 4 PCs): "
                  f"{len(rnk)} genes, running prerank...")
            res = gp.prerank(
                rnk=rnk, gene_sets=gene_sets_list,
                min_size=15, max_size=500,
                permutation_num=1000, outdir=None,
                no_plot=True, seed=42, threads=4, verbose=False,
            )
            df_res = res.res2d.copy()
            sig = df_res[df_res["FDR q-val"].astype(float) < FDR_CUTOFF]
            print(f"    {len(sig)} sig terms (FDR<{FDR_CUTOFF})")
        else:
            tsv = os.path.join(GSEA_DIR, tissue, "CorrAdjust_gsea_cp.tsv")
            df_res = pd.read_csv(tsv, sep="\t")
        results["CorrAdjust"].append((tissue, df_res))

        # -- COBRA --
        tsv = os.path.join(GSEA_DIR, tissue, "COBRA_deprecated_gsea_cp.tsv")
        df_res = pd.read_csv(tsv, sep="\t")
        results["COBRA"].append((tissue, df_res))

    return results


# ── Dotplot: ORA ──────────────────────────────────────────────────
def build_ora_plot_data(tissue_results):
    records = []
    for tissue, df_enr in tissue_results:
        if df_enr.empty:
            continue
        sig = df_enr[df_enr["padj"] < PADJ_CUTOFF].head(TOP_K_PLOT)
        for _, row in sig.iterrows():
            records.append({
                "tissue": short_tissue(tissue),
                "term": clean_term(row["term"]),
                "padj": row["padj"],
                "overlap": row["overlap"],
                "neg_log_padj": min(-np.log10(max(row["padj"], 1e-300)), 150),
            })
    return pd.DataFrame(records) if records else pd.DataFrame()


# ── Dotplot: GSEA ─────────────────────────────────────────────────
def build_gsea_plot_data(tissue_results):
    records = []
    for tissue, df_res in tissue_results:
        if df_res.empty:
            continue
        df_res = df_res.copy()
        df_res["FDR q-val"] = df_res["FDR q-val"].astype(float)
        df_res["NES"] = df_res["NES"].astype(float)
        sig = df_res[df_res["FDR q-val"] < FDR_CUTOFF].copy()
        sig["abs_NES"] = sig["NES"].abs()
        sig = sig.sort_values("abs_NES", ascending=False).head(TOP_K_PLOT)
        for _, row in sig.iterrows():
            records.append({
                "tissue": short_tissue(tissue),
                "term": clean_term(row["Term"]),
                "fdr": row["FDR q-val"],
                "nes": row["NES"],
                "abs_nes": abs(row["NES"]),
                "neg_log_fdr": min(-np.log10(max(row["FDR q-val"], 1e-300)), 150),
            })
    return pd.DataFrame(records) if records else pd.DataFrame()


# ── Generic 2-panel dotplot ───────────────────────────────────────
def plot_two_panel_dotplot(data_left, data_right, label_left, label_right,
                           size_col, color_col, color_label, size_label,
                           suptitle, cmap="YlOrRd", panel_labels=None):
    """
    data_left / data_right: DataFrames with 'tissue', 'term' columns
                            plus size_col and color_col.
    panel_labels: optional tuple like ("A", "B") for panel titles.
    """
    fig, axes = plt.subplots(1, 2, figsize=(18, 10))

    for idx, (ax, df_plot, method_label) in enumerate([
        (axes[0], data_left, label_left),
        (axes[1], data_right, label_right),
    ]):
        if panel_labels:
            ax.text(-0.02, 1.05, panel_labels[idx], transform=ax.transAxes,
                    fontsize=18, fontweight="bold", va="bottom", ha="right")
        if df_plot.empty:
            ax.text(0.5, 0.5, "No significant terms",
                    ha="center", va="center", fontsize=12, transform=ax.transAxes)
            ax.set_title(method_label, fontsize=11, fontweight="bold")
            continue

        tissues = sorted(df_plot["tissue"].unique(),
                         key=lambda t: TISSUES.index("GTEx-" + t)
                         if "GTEx-" + t in TISSUES else 999)
        term_order = (
            df_plot.groupby("term")
            .agg(count=(color_col, "size"), best=(color_col, "max"))
            .sort_values(["count", "best"], ascending=[False, False])
            .index.tolist()
        )

        ds_idx = {d: i for i, d in enumerate(tissues)}
        term_idx = {t: i for i, t in enumerate(term_order)}

        max_color = df_plot[color_col].max()
        norm = Normalize(vmin=0, vmax=max_color)
        max_size_val = max(df_plot[size_col].max(), 1)

        for _, row in df_plot.iterrows():
            x = ds_idx[row["tissue"]]
            y = term_idx[row["term"]]
            s = 30 + 200 * (row[size_col] / max_size_val)
            c = plt.cm.get_cmap(cmap)(norm(row[color_col]))
            ax.scatter(x, y, s=s, c=[c], edgecolors="grey",
                       linewidths=0.5, zorder=3)

        ax.set_xticks(range(len(tissues)))
        ax.set_xticklabels(tissues, rotation=60, ha="right", fontsize=8)
        ax.set_yticks(range(len(term_order)))
        ax.set_yticklabels(term_order, fontsize=7)
        ax.set_xlim(-0.5, len(tissues) - 0.5)
        ax.set_ylim(-0.5, len(term_order) - 0.5)
        ax.invert_yaxis()
        ax.grid(True, alpha=0.15, zorder=0)
        ax.set_title(method_label, fontsize=11, fontweight="bold")

        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
        cbar.set_label(color_label, fontsize=9)

        for frac, label in [(0.2, str(int(max_size_val * 0.2))),
                             (0.6, str(int(max_size_val * 0.6))),
                             (1.0, str(int(max_size_val)))]:
            ax.scatter([], [], s=30 + 200 * frac, c="grey", alpha=0.5,
                       edgecolors="grey", label=f"{label} {size_label}")
        ax.legend(title=size_label.capitalize(), loc="lower right", fontsize=7,
                  title_fontsize=8, framealpha=0.9)

    fig.suptitle(suptitle, fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    return fig


# ── Main ──────────────────────────────────────────────────────────
def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    # ── Figure S1: ORA top-100 ────────────────────────────────────
    print("=" * 60)
    print("  Figure S1: ORA top-100 (Canonical Pathways)")
    print("=" * 60)
    ora = load_ora_results()
    df_left = build_ora_plot_data(ora["CorrAdjust"])
    df_right = build_ora_plot_data(ora["COBRA"])
    print(f"  CorrAdjust: {len(df_left)} plot records")
    print(f"  COBRA:      {len(df_right)} plot records")

    fig1 = plot_two_panel_dotplot(
        df_left, df_right,
        label_left="CorrAdjust", label_right="COBRA",
        size_col="overlap", color_col="neg_log_padj",
        color_label="-log\u2081\u2080(padj)", size_label="genes",
        suptitle="Figure S1 — ORA (top 100 genes) — Canonical Pathways",
        panel_labels=("A", "B"),
    )
    for ext in ["png", "pdf"]:
        path = os.path.join(FIG_DIR, f"figs1_ora_cp.{ext}")
        fig1.savefig(path, dpi=200, bbox_inches="tight")
        print(f"  Saved → {path}")
    plt.close(fig1)

    # ── Figure S2: GSEA prerank ───────────────────────────────────
    print("\n" + "=" * 60)
    print("  Figure S2: GSEA prerank (Canonical Pathways)")
    print("=" * 60)
    gsea = load_gsea_results()
    df_left = build_gsea_plot_data(gsea["CorrAdjust"])
    df_right = build_gsea_plot_data(gsea["COBRA"])
    print(f"  CorrAdjust: {len(df_left)} plot records")
    print(f"  COBRA:      {len(df_right)} plot records")

    fig2 = plot_two_panel_dotplot(
        df_left, df_right,
        label_left="CorrAdjust", label_right="COBRA",
        size_col="abs_nes", color_col="neg_log_fdr",
        color_label="-log\u2081\u2080(FDR)", size_label="|NES|",
        suptitle="Figure S2 — GSEA prerank — Canonical Pathways",
    )
    for ext in ["png", "pdf"]:
        path = os.path.join(FIG_DIR, f"figs2_gsea_cp.{ext}")
        fig2.savefig(path, dpi=200, bbox_inches="tight")
        print(f"  Saved → {path}")
    plt.close(fig2)

    print("\nDone!")


if __name__ == "__main__":
    main()
