# corr-bench

Plotting code for the figures of *"Quality not quantity of enriched pathways
should guide benchmarks of batch correction methods in biological networks"*
(Babaie, Saha, Quackenbush, Ben Guebila).

This repository contains only the final plotting scripts and their output
figures. Intermediate analysis results (per-tissue WGCNA modules, ORA/GSEA
tables, COBRA train/test R², runtime benchmarks) will be released on Zenodo
(DOI: TBD) and should be placed under `results/` and `data/` matching the
paths below.

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

After downloading the Zenodo archive, files should sit at:

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
