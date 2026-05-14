#!/usr/bin/env python3
"""plot_3bar_precision.py — 3-bar tissue-relevance precision: CorrAdjust, COBRA, COBRA+Hardy"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.expanduser("~/corrrepro")
CA_DIR     = os.path.join(BASE, "results", "wgcna_ora")
COBRA_DIR  = os.path.join(BASE, "results", "wgcna_ora_cov_all")
HARDY_DIR  = os.path.join(BASE, "results", "wgcna_ora_hardy")
FIG_DIR    = os.path.join(BASE, "figures")
PADJ = 0.05

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

TISSUE_KEYWORDS = {
    "GTEx-Adipose-Subcutaneous": ["adipos", "adipocyte", "adipogenesis", "lipid", "lipolysis",
        "fatty_acid", "fat_cell", "triglyceride", "lipase", "lipid_storage", "lipid_droplet",
        "lipid_metabol", "lipid_biosynth", "lipid_catabol", "lipid_transport", "cholesterol",
        "insulin", "glucose", "thermogen"],
    "GTEx-Artery-Tibial": ["arteri", "vascul", "angiogen", "endotheli", "smooth_muscle",
        "blood_vessel", "blood_pressure", "vasodilat", "vasoconstrict", "vessel_morpho",
        "vessel_develop", "atherosclero", "arterial", "vessel_remodel", "circulatory"],
    "GTEx-Heart-Left-Ventricle": ["cardiac", "heart", "cardiomyocyte", "myocardi", "ventricl",
        "atrial", "atrioventricular", "cardiac_muscle", "heart_develop", "heart_morpho",
        "cardiac_conduct", "sarcomere", "cardiomyopath", "cardiac_contract", "muscle_contract"],
    "GTEx-Lung": ["lung", "pulmonary", "alveol", "bronch", "respirat", "airway", "surfactant",
        "pneumocyte", "trachea", "cilium", "cilia", "mucociliary"],
    "GTEx-Muscle-Skeletal": ["skeletal_muscle", "muscle_contract", "muscle_cell", "muscle_develop",
        "muscle_differ", "muscle_fiber", "muscle_tissue", "myofibril", "myotube", "myoblast",
        "sarcomere", "striated_muscle", "muscle_struct", "myogenesis", "muscle_organ",
        "muscle_system", "muscle_adapt"],
    "GTEx-Nerve-Tibial": ["neuro", "nerve", "axon", "synap", "dendrit", "myelin", "glia",
        "schwann", "neurotransmit", "neural", "neurogenesis", "neuron", "ensheathment",
        "nerve_develop", "gliogenesis", "oligodendrocyte", "peripheral_nervous",
        "nervous_system", "action_potential", "glutamate_receptor",
        "membrane_potential", "neurite", "l1cam", "ncam",
        "locomotor", "walking_behavior",
        "cell_morphogenesis_involved_in_neuron", "kinesin"],
    "GTEx-Skin-Sun-Exposed-Lower-leg": ["skin", "epiderm", "keratinocyte", "melanocyte", "melanin",
        "cornif", "dermis", "hair_follicle", "wound_heal", "keratiniz", "pigment",
        "skin_develop", "skin_barrier", "epithelial_cell_differ"],
    "GTEx-Small-Intestine-Terminal-Ileum": ["intestin", "digest", "gut", "entero", "absorpt",
        "intestinal_epithe", "villus", "ileum", "ileal", "colon", "bowel", "mucosal",
        "gastrointestin", "brush_border", "digestive"],
    "GTEx-Spleen": ["spleen", "splenic", "lymphocyte", "immune", "t_cell", "b_cell", "leukocyte",
        "hematopoie", "immunoglobulin", "antigen", "lymph", "cytokine", "innate_immune",
        "adaptive_immune", "immune_response", "complement", "phagocyt", "myeloid"],
    "GTEx-Thyroid": ["thyroid", "thyrocyte", "thyroxine", "iodide", "thyroid_hormone",
        "thyroglobulin", "endocrine", "hormone_secret", "hormone_metabol", "hormone_biosynth",
        "skeletal", "ossification"],
    "GTEx-Whole-Blood": ["erythrocyte", "hemoglobin", "heme", "blood_cell", "red_blood",
        "platelet", "megakaryocyte", "coagulat", "thrombocyte", "erythropoie", "oxygen_transport",
        "gas_transport", "porphyrin", "iron_homeostasis", "hemostasis", "blood_coagul"],
}


def load_sig(base_dir, tissue, method_sub):
    frames = []
    search_dir = os.path.join(base_dir, tissue, method_sub)
    for fname in ("ora_gobp.tsv", "ora_cp.tsv"):
        p = os.path.join(search_dir, fname)
        if os.path.exists(p):
            df = pd.read_csv(p, sep="\t")
            if "FDR" in df.columns:
                frames.append(df[df["FDR"] < PADJ])
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def is_relevant(term, keywords):
    tl = term.lower()
    return any(k in tl for k in keywords)


def best_module(df, kws):
    if df.empty:
        return None
    best = None
    for mod in df["module"].unique():
        mdf = df[df["module"] == mod]
        total = len(mdf)
        if total == 0:
            continue
        rel = sum(1 for t in mdf["Term"] if is_relevant(t, kws))
        if rel == 0:
            continue
        prec = rel / total * 100
        if best is None or prec > best[1] or (prec == best[1] and rel > best[2]):
            best = (mod, prec, rel, total)
    return best


def main():
    os.makedirs(FIG_DIR, exist_ok=True)

    methods = [
        ("CorrAdjust",    CA_DIR,    "CorrAdjust", "#1f77b4"),
        ("COBRA",         COBRA_DIR, "MLE",        "#d62728"),
        ("COBRA + Hardy", HARDY_DIR, "MLE",        "#ff7f0e"),
    ]

    n_methods = len(methods)
    n_tissues = len(TISSUES)

    prec_data = {}
    for method_name, base_dir, sub, color in methods:
        prec_list = []
        for tissue in TISSUES:
            kws = TISSUE_KEYWORDS[tissue]
            df = load_sig(base_dir, tissue, sub)
            b = best_module(df, kws)
            if b:
                prec_list.append(b[1])
            elif df.empty:
                prec_list.append(None)
            else:
                prec_list.append(0)
        prec_data[method_name] = prec_list

    # Print summary
    header = "{:<40s}".format("Tissue")
    for m, _, _, _ in methods:
        header += "{:>18s}".format(m)
    print("\nPrecision Summary:")
    print(header)
    for i, tissue in enumerate(TISSUES):
        short = tissue.replace("GTEx-", "")
        row = "{:<40s}".format(short)
        for m, _, _, _ in methods:
            v = prec_data[m][i]
            if v is None:
                row += "{:>18s}".format("N/A")
            else:
                row += "{:>17.1f}%".format(v)
        print(row)

    # Plot
    short_labels = [t.replace("GTEx-", "") for t in TISSUES]
    x = np.arange(n_tissues)
    total_width = 0.7
    bar_width = total_width / n_methods

    fig, ax = plt.subplots(figsize=(16, 6))

    for j, (method_name, _, _, color) in enumerate(methods):
        positions = x - total_width/2 + bar_width/2 + j * bar_width
        values = [v if v is not None else 0 for v in prec_data[method_name]]

        bars = ax.bar(positions, values, bar_width * 0.9, label=method_name,
                      color=color, alpha=0.85, edgecolor="white", linewidth=0.5)

        for i, v in enumerate(prec_data[method_name]):
            if v is None:
                ax.text(positions[i], 1, "N/A", ha="center", va="bottom",
                        fontsize=6, color="gray", rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, fontsize=16, rotation=35, ha="right")
    ax.set_ylabel("Precision (%)", fontsize=16)
    ax.tick_params(axis="y", labelsize=14)
    ax.legend(fontsize=14, loc="upper right")
    ax.grid(axis="y", alpha=0.3)

    all_vals = [v for m in prec_data for v in prec_data[m] if v is not None and v > 0]
    ax.set_ylim(0, max(all_vals) * 1.15 if all_vals else 10)

    fig.tight_layout()
    out_png = os.path.join(FIG_DIR, "fig_3bar_precision.png")
    out_pdf = os.path.join(FIG_DIR, "fig_3bar_precision.pdf")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print("\nSaved:", out_png)
    print("Saved:", out_pdf)


if __name__ == "__main__":
    main()
