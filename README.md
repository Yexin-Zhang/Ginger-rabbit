# Ginger rabbits: analysis code

Code for Alves *et al.*, "Camouflage against eagles favours a domestic allele in wild rabbits" (submitted).

Sequence processing, variant calling, phasing and population-genetic analyses used published software, with versions and parameters given in the Methods and Supplementary Information. This repository contains the R scripts for the predation, selection, substrate-colour, PCA and MC1R tree analyses.

## Contents

| Script | Analysis | Supports |
| --- | --- | --- |
| `scripts/01_selection_model.R` | Eagle predation rate required for the ginger allele to reach its observed frequency | Extended Data Fig. 6, Supplementary Table 6 |
| `scripts/02_foam_model_tests.R` | Fisher's exact tests and relative risk for the foam-model predation experiment | Fig. 5d |
| `scripts/03_eagle_prey_glmm.R` | Binomial GLMM comparing ginger frequency in eagle prey and field surveys | Fig. 6c |
| `scripts/04_substrate_colour_correlation.R` | Weighted correlations between ginger frequency and substrate redness (a*) | Fig. 4a,b |
| `scripts/05_pca.R` | PCA plots of genome-wide and MC1R-region variation | Fig. 3a |
| `scripts/06_upgma_tree.R` | UPGMA tree of the MC1R region | Fig. 2c,d |

Input data are in `data/`. The location tables for script 04 are taken from Supplementary Data 6, and the input tables for script 03 from Supplementary Data 8 and 2. Scripts 01 and 02 need no input files.

## System requirements

- R 4.3.3, tested on Linux and macOS.
- R packages: `ggplot2`, `dplyr`, `glmmTMB` (tested with 1.1.8 and 1.1.11), `cowplot`, `ape` and `phangorn`.
- No non-standard hardware is required.

## Installation

```r
install.packages(c("ggplot2", "dplyr", "glmmTMB", "cowplot", "ape", "phangorn"))
```

Typical install time on a normal desktop computer: under 5 minutes.

## Demo

Run from the repository root:

```bash
Rscript scripts/01_selection_model.R
Rscript scripts/02_foam_model_tests.R
Rscript scripts/03_eagle_prey_glmm.R
Rscript scripts/04_substrate_colour_correlation.R
Rscript scripts/05_pca.R
Rscript scripts/06_upgma_tree.R
```

Expected run time: under 60 seconds per script on a normal desktop computer. Outputs are written to `results/`; the expected outputs are included in this repository.

MIT. See `LICENSE`.
