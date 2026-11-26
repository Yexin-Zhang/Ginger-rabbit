###############################################################################
## 04_substrate_colour_correlation.R
## Correlation between ginger frequency and substrate redness (a*) across
## Australian sampling locations, weighted by the number of rabbits.
## Supports: Fig. 4a,b (scatter plots)
## Input:  data/fig4_genotype_locations.csv  (Supplementary Data 6)
##         data/fig4_phenotype_locations.csv (Supplementary Data 6)
## Output: results/04_substrate_colour/
## Run from the repository root: Rscript scripts/04_substrate_colour_correlation.R
###############################################################################

library(ggplot2)

out_dir <- "results/04_substrate_colour"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genotype <- read.csv("data/fig4_genotype_locations.csv", check.names = FALSE)
phenotype <- read.csv("data/fig4_phenotype_locations.csv", check.names = FALSE)

datasets <- list(
  genotype = data.frame(
    a_star    = genotype[["a*"]],
    frequency = genotype[["Ginger allele frequency (HW-corrected)"]],
    n         = genotype[["Individuals genotyped"]]
  ),
  phenotype = data.frame(
    a_star    = phenotype[["a*"]],
    frequency = phenotype[["Ginger frequency"]],
    n         = phenotype[["Rabbits recorded"]]
  )
)

y_labels <- c(
  genotype  = "Ginger allele frequency (%)",
  phenotype = "Ginger rabbit frequency (%)"
)

# Pearson correlation weighted by sample size; P from the weighted
# least-squares slope
weighted_correlation <- function(d) {
  r <- cov.wt(cbind(d$a_star, d$frequency), wt = d$n / sum(d$n),
              cor = TRUE)$cor[1, 2]
  p <- summary(lm(frequency ~ a_star, data = d, weights = n))$coefficients[2, 4]
  c(r = r, p = p)
}

results <- data.frame()

for (nm in names(datasets)) {
  d  <- datasets[[nm]]
  wc <- weighted_correlation(d)

  results <- rbind(results, data.frame(
    dataset   = nm,
    locations = nrow(d),
    rabbits   = sum(d$n),
    r         = round(wc[["r"]], 2),
    p         = signif(wc[["p"]], 2)
  ))

  label <- sprintf("R = %.2f, P = %s", wc[["r"]],
                   formatC(wc[["p"]], format = "e", digits = 1))

  plot <- ggplot(d, aes(x = a_star, y = frequency * 100)) +
    geom_smooth(aes(weight = n), method = "lm", formula = y ~ x, se = FALSE,
                colour = "black", linewidth = 0.4, linetype = "dotted") +
    geom_point(aes(size = n, fill = frequency * 100), shape = 21,
               colour = "white", stroke = 0.5) +
    scale_fill_gradient(low = "#BDBDBD", high = "#E67F52", limits = c(0, 50),
                        oob = scales::squish, name = "Frequency (%)") +
    scale_size_area(max_size = 10, name = "Rabbits") +
    annotate("text", x = Inf, y = Inf, label = label,
             hjust = 1.05, vjust = 1.6, size = 3.8) +
    labs(x = "Substrate redness (a*)", y = y_labels[[nm]]) +
    theme_classic()

  ggsave(file.path(out_dir, sprintf("fig4_%s_correlation.pdf", nm)),
         plot, width = 6, height = 5)
}

print(results, row.names = FALSE)
write.csv(results, file.path(out_dir, "fig4_correlations.csv"), row.names = FALSE)
