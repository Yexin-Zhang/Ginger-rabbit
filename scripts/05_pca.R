###############################################################################
## 05_pca.R
## Principal component analysis of 88 rabbits with known phenotype, using
## genome-wide intergenic capture targets (chromosome 18 excluded) and the
## region flanking MC1R (NC_091449.1:68.7-69.4 Mb). PCs were computed in
## PLINK v1.9; this script plots them.
## Supports: Fig. 3a
## Input:  data/pca/  (PLINK .eigenvec/.eigenval files and sample information)
## Output: results/05_pca/
## Run from the repository root: Rscript scripts/05_pca.R
###############################################################################

library(dplyr)
library(ggplot2)
library(cowplot)

out_dir <- "results/05_pca"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.delim("data/pca/pca_samples.txt")

regions <- list(
  genome_wide = list(prefix = "data/pca/genome_wide_excl_chr18",
                     title  = "Genome-wide (excluding chr 18)"),
  mc1r        = list(prefix = "data/pca/mc1r_region",
                     title  = "MC1R region")
)

genotype_colours <- c(
  "Homozygous ginger" = "#E67F52",
  "Heterozygous grey" = "#F6D6BC",
  "Homozygous grey"   = "#BDBDBD"
)
population_shapes <- c("Domestic" = 17, "Australia" = 16)

plot_pca <- function(region) {
  pcs <- read.table(paste0(region$prefix, ".eigenvec"))
  colnames(pcs) <- c("FID", "Sample", paste0("PC", seq_len(ncol(pcs) - 2)))
  eigenvalues <- scan(paste0(region$prefix, ".eigenval"), quiet = TRUE)
  pct <- eigenvalues / sum(eigenvalues) * 100

  d <- pcs %>%
    inner_join(samples, by = "Sample") %>%
    mutate(
      Genotype = case_when(
        Population == "Domestic" | Phenotype == "Ginger" ~ "Homozygous ginger",
        Pheno_Category == "Hetero_Grey" ~ "Heterozygous grey",
        Pheno_Category == "Homo_Grey"   ~ "Homozygous grey"
      ),
      Genotype   = factor(Genotype, levels = names(genotype_colours)),
      Population = factor(Population, levels = names(population_shapes))
    )

  ggplot(d, aes(PC1, PC2, colour = Genotype, shape = Population)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_colour_manual(values = genotype_colours) +
    scale_shape_manual(values = population_shapes) +
    labs(title = region$title,
         x = sprintf("PC1 (%.1f%%)", pct[1]),
         y = sprintf("PC2 (%.1f%%)", pct[2])) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
}

plots <- lapply(regions, plot_pca)
legend <- get_legend(plots[[1]])
plots <- lapply(plots, function(p) p + theme(legend.position = "none"))

figure <- plot_grid(plots$genome_wide, plots$mc1r, legend,
                    nrow = 1, rel_widths = c(1, 1, 0.5))
ggsave(file.path(out_dir, "pca.pdf"), figure, width = 10, height = 3.5)
