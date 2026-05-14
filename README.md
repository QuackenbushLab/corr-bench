# corr-bench

Plotting code for the figures of *"Quality not quantity of enriched pathways
should guide benchmarks of batch correction methods in biological networks"*
(Babaie, Saha, Quackenbush, Ben Guebila).

This repository contains only the final plotting scripts and their output
figures. Analysis results (per-tissue WGCNA modules, ORA/GSEA tables, COBRA
train/test R², runtime benchmarks, GMT files) are released on Zenodo
(DOI: TBD) and should be placed under `results/` and `data/` matching the
paths below.

### Gene filtering note

All COBRA and CorrAdjust results in `results/tpmfilter_data/`,
`results/wgcna_ora_cov_all/`, and `results/wgcna_ora_hardy/` use the **same
gene set per tissue**: genes with median raw count ≥ 1 across the training
samples (the legacy folder name `tpmfilter_data` predates the switch from a
TPM-based filter; despite the name, it is a count-based filter). Comparisons
in Figs 2 and S3 are therefore apples-to-apples on genes; observed
differences reflect the methods.

## Data (Zenodo)

The deposit contains 5 tarballs and 2 small TSVs. After download, extract
each tarball at the indicated location:

| Archive | Extract at | Used by |
|---|---|---|
| `tpmfilter_data.tar.gz` (~2.2 GB) | `results/` → `results/tpmfilter_data/` | Figs 1, 2, 3, S1, S2, S3 |
| `wgcna_ora_cov_all.tar.gz` (~19 MB) | `results/` → `results/wgcna_ora_cov_all/` | Figs 2, 3, S3 |
| `wgcna_ora_hardy.tar.gz` (~13 MB) | `results/` → `results/wgcna_ora_hardy/` | Fig 2 (Hardy tissues) |
| `GMT_files.tar.gz` (~2.2 MB) | `data/input_data/` → `data/input_data/GMT_files/` | Figs S1, S2 |
| `ora_top100_all.tar.gz` (~1.9 MB) | `figures/misc/` → `figures/misc/ora_top100_all/` | Figs S1, S2 |
| `cobra_train_test_r2.tsv` (4 KB) | `results/` | Fig 4 |
| `benchmark_medcount_summary.tsv` (4 KB) | `results/` | Fig 3 (panel B, runtime) |

Quick extract (run from repo root):

```bash
tar -xzf tpmfilter_data.tar.gz       -C results/
tar -xzf wgcna_ora_cov_all.tar.gz    -C results/
tar -xzf wgcna_ora_hardy.tar.gz      -C results/
tar -xzf GMT_files.tar.gz            -C data/input_data/
mkdir -p figures/misc && tar -xzf ora_top100_all.tar.gz -C figures/misc/
mv cobra_train_test_r2.tsv benchmark_medcount_summary.tsv results/
```

## Figures

| Figure | Script | Generates |
|---|---|---|
| Fig 1  | [src/fig1_reproduce_5b_bottom.py](src/fig1_reproduce_5b_bottom.py)   | [figures/fig1_ges_gtex.png](figures/fig1_ges_gtex.png) |
| Fig 2  | [src/plot_fig2_module_stacked.py](src/plot_fig2_module_stacked.py)   | [figures/fig2_module_stacked.png](figures/fig2_module_stacked.png) |
| Fig 3A | [src/plot_3bar_precision.py](src/plot_3bar_precision.py)             | `figures/fig_3bar_precision.png` (intermediate) |
| Fig 3  | [src/plot_fig3_combined.py](src/plot_fig3_combined.py)               | [figures/fig3_combined.png](figures/fig3_combined.png) |
| Fig 4  | [src/18_plot_cobra_r2.py](src/18_plot_cobra_r2.py)                   | [figures/fig4_cobra_train_test_r2.png](figures/fig4_cobra_train_test_r2.png) |
| Fig S1 | [src/plot_figs1_figs2.py](src/plot_figs1_figs2.py)                   | [figures/figs1_ora_cp.png](figures/figs1_ora_cp.png) |
| Fig S2 | [src/plot_figs1_figs2.py](src/plot_figs1_figs2.py)                   | [figures/figs2_gsea_cp.png](figures/figs2_gsea_cp.png) |
| Fig S3 | [src/plot_fig3_summary_combined.py](src/plot_fig3_summary_combined.py) | [figures/figs3_summary_gobp_combined.png](figures/figs3_summary_gobp_combined.png) |

## Expected data layout

After extracting the Zenodo tarballs (see [Data (Zenodo)](#data-zenodo)),
files sit at:

```
results/
  benchmark_medcount_summary.tsv             # fig 3 (panel B)
  cobra_train_test_r2.tsv                    # fig 4
  wgcna_ora_cov_all/<tissue>/MLE/{module_assignments,ora_gobp,ora_cp}.tsv
  wgcna_ora_hardy/<tissue>/MLE/{module_assignments,ora_gobp,ora_cp}.tsv
  tpmfilter_data/
    wgcna_ora/<tissue>/CorrAdjust/{module_assignments,ora_gobp,ora_cp}.tsv
    gsea_all/<tissue>/...
    corradjust_gtex_outputs_models/<tissue>/...

data/
  input_data/GMT_files/c2.cp.v2023.2.Hs.symbols.gmt
```

`<tissue>` is one of the 11 GTEx tissues used in the paper
(e.g. `GTEx-Adipose-Subcutaneous`, `GTEx-Nerve-Tibial`).

## Install

```bash
pip install -r requirements.txt
```

## Reproduce

```bash
python src/fig1_reproduce_5b_bottom.py
python src/plot_fig2_module_stacked.py
python src/plot_3bar_precision.py        # produces fig 3 panel A
python src/plot_fig3_combined.py         # combines panel A + runtime panel B
python src/18_plot_cobra_r2.py
python src/plot_figs1_figs2.py
python src/plot_fig3_summary_combined.py
```

## License

MIT — see [LICENSE](LICENSE).
