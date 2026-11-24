# glycoQTL Manuscript Resources

Code and cached data to reproduce the analyses and figures in the glycoQTL manuscript.

## About glycoQTL

glycoQTL maps genetic regulation of site-specific protein glycosylation (glycopeptiforms) in aging human brain. The study uses glycoproteomics and genotype data from ROSMAP, and integrates glycoQTL with fine-mapping, multi-omics colocalization (ColocBoost) and GlycoWAS (TWAS framework) for Alzheimer's disease.

## Repository Structure

This repository contains the code and cached data used to generate the analyses and figures in the manuscript, available at: https://github.com/rfeng2023/glycoQTL

- **1_data_analysis**: Scripts for SuSiE fine-mapping with GlycoWAS weights (`2026_run_mnm_array.slurm`), genome-wide GlycoWAS (`2026_run_twas.slurm`) and multi-omics colocalization with ColocBoost (`2026_run_colocboost.sh`).

- **2_main_figures**: One folder per figure (`Figure1` to `Figure6`). Each folder contains a `*_preparation` notebook that builds the cached RDS files and a `*_plots` notebook that assembles the figure from `data/*.rds`.

- **3_supplementary_figures**: One folder per figure (`Supplementary_Figure1` to `Supplementary_Figure6`), same layout as the main figures.

- **4_supplementary_tables**: Scripts for Supplementary Tables S1 to S5 and the Excel builder.

Cached inputs for each figure are in the `data/` folder of that figure.

## Getting Started

Suggested order: `1_data_analysis` (fine-mapping, GlycoWAS), then the `*_preparation` notebooks, then the `*_plots` notebooks.

## Computational Requirements

The analyses were performed using:

- R version 4.4 (data.table, ggplot2, dplyr, patchwork, cowplot)
- Key packages: pecotmr, susieR, colocboost
- xqtl-protocol workflows (https://github.com/statfungen/xqtl-protocol)

## Citation

If you use these resources, please cite the glycoQTL manuscript (citation to be added).
